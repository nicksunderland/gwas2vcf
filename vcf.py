import logging
import pickle
from heapq import heappop

import numpy as np
import pysam


class Vcf:

    # no longer necessary
    # @staticmethod
    # def convert_pval_to_neg_log10(p):
    #     # prevent negative 0 output
    #     if p == 1:
    #         return 0
    #     # prevent Inf output
    #     if p == 0:
    #         return 999
    #     return -np.log10(p)

    @staticmethod
    def is_float32_lossy(input_float):
        if (
            input_float == 0
            or input_float is None
            or input_float == np.inf
            or input_float == -np.inf
        ):
            return False

        # convert val to float32
        out_float = np.float32(input_float)

        return out_float in [0, np.inf, 0, -np.inf]

    @staticmethod
    def remove_illegal_chars(input_string):
        if input_string is None:
            return None
        illegal_chars = (" ", ";", ":", "=", ",")
        out_string = input_string.strip()
        for bad_char in illegal_chars:
            out_string = out_string.replace(bad_char, "_")
        return out_string

    """
    Write GWAS file to VCF
    Expects an open file handle to a Pickle file of GWAS results & file index dict(chromosome[(position, offset)]) 
    """

    @staticmethod
    def write_to_file(
        gwas_file,
        gwas_idx,
        path,
        fasta,
        build,
        trait_id,
        sample_metadata=None,
        file_metadata=None,
        csi=False,
    ):
        logging.info(f"Writing headers to BCF/VCF: {path}")

        header = pysam.VariantHeader()

        # INFO
        header.add_line(
            '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">'
        )

        # FORMAT
        header.add_line(
            '##FORMAT=<ID=ES,Number=A,Type=Float,Description="Effect size estimate relative to the alternative allele">'
        )
        header.add_line(
            '##FORMAT=<ID=SE,Number=A,Type=Float,Description="Standard error of effect size estimate">'
        )
        header.add_line(
            '##FORMAT=<ID=LP,Number=A,Type=Float,Description="-log10 p-value for effect estimate">'
        )
        header.add_line(
            '##FORMAT=<ID=AF,Number=A,Type=Float,Description="Alternate allele frequency in the association study">'
        )
        header.add_line(
            '##FORMAT=<ID=SS,Number=A,Type=Integer,Description="Sample size used to estimate genetic effect">'
        )
        header.add_line(
            '##FORMAT=<ID=EZ,Number=A,Type=Float,Description="Z-score provided if it was used to derive the EFFECT and SE fields">'
        )
        header.add_line(
            '##FORMAT=<ID=SI,Number=A,Type=Float,Description="Accuracy score of summary data imputation">'
        )
        header.add_line(
            '##FORMAT=<ID=NC,Number=A,Type=Integer,Description="Number of cases used to estimate genetic effect">'
        )
        header.add_line(
            '##FORMAT=<ID=ID,Number=1,Type=String,Description="Study variant identifier">'
        )

        # META
        header.add_line(
            '##META=<ID=TotalVariants,Number=1,Type=Integer,Description="Total number of variants in input">'
        )
        header.add_line(
            '##META=<ID=VariantsNotRead,Number=1,Type=Integer,Description="Number of variants that could not be read">'
        )
        header.add_line(
            '##META=<ID=HarmonisedVariants,Number=1,Type=Integer,Description="Total number of harmonised variants">'
        )
        header.add_line(
            '##META=<ID=VariantsNotHarmonised,Number=1,Type=Integer,Description="Total number of variants that could not be harmonised">'
        )
        header.add_line(
            '##META=<ID=SwitchedAlleles,Number=1,Type=Integer,Description="Total number of variants strand switched">'
        )
        header.add_line(
            '##META=<ID=TotalControls,Number=1,Type=Integer,Description="Total number of controls in the association study">'
        )
        header.add_line(
            '##META=<ID=TotalCases,Number=1,Type=Integer,Description="Total number of cases in the association study">'
        )
        header.add_line(
            '##META=<ID=StudyType,Number=1,Type=String,Description="Type of GWAS study [Continuous or CaseControl]">'
        )

        # SAMPLES
        header.samples.add(trait_id)
        if file_metadata is not None:
            meta_string = "".join(f",{k}={sample_metadata[k]}" for k in sample_metadata)
            header.add_line(f"##SAMPLE=<ID={trait_id}{meta_string}>")

        # CONTIG
        assert len(fasta.references) == len(fasta.lengths)
        for n, contig in enumerate(fasta.references):
            header.add_line(
                f"##contig=<ID={contig},length={fasta.lengths[n]}, assembly={build}>"
            )

        # add metadata
        if file_metadata is not None:
            for k in file_metadata:
                header.add_line(f"##{k}={file_metadata[k]}")

        vcf = pysam.VariantFile(path, "w", header=header)

        # recall variant objects in chromosome position order
        logging.info(f"Writing variants to BCF/VCF: {path}")
        for contig in fasta.references:
            if contig not in gwas_idx:
                continue

            while gwas_idx[contig]:
                chr_pos = heappop(gwas_idx[contig])

                # load GWAS result
                gwas_file.seek(chr_pos[1])
                result = pickle.load(gwas_file)

                result.nlog_pval = result.nlog_pval

                # check floats
                if Vcf.is_float32_lossy(result.b):
                    logging.warning(
                        f"Effect field cannot fit into float32. Expect loss of precision for: {result.b}"
                    )

                if Vcf.is_float32_lossy(result.se):
                    result.se = np.float64(np.finfo(np.float32).tiny).item()
                    logging.warning(
                        f"Standard error field cannot fit into float32. Expect loss of precision for: {result.se}"
                    )
                if Vcf.is_float32_lossy(result.nlog_pval):
                    logging.warning(
                        f"-log10(pval) field cannot fit into float32. Expect loss of precision for: {result.nlog_pval}"
                    )
                if Vcf.is_float32_lossy(result.alt_freq):
                    logging.warning(
                        f"Allele frequency field cannot fit into float32. Expect loss of precision for: {result.alt_freq}"
                    )
                if Vcf.is_float32_lossy(result.imp_z):
                    logging.warning(
                        f"Imputation Z score field cannot fit into float32. Expect loss of precision for: {result.imp_z}"
                    )
                if Vcf.is_float32_lossy(result.imp_info):
                    logging.warning(
                        f"Imputation INFO field cannot fit into float32. Expect loss of precision for: {result.imp_info}"
                    )

                record = vcf.new_record()
                record.chrom = result.chrom
                assert " " not in record.chrom
                record.pos = result.pos
                assert record.pos > 0
                record.id = Vcf.remove_illegal_chars(result.dbsnpid)
                record.alleles = (result.ref, result.alt)
                record.filter.add(result.vcf_filter)

                if result.alt_freq is not None:
                    record.info["AF"] = result.alt_freq

                if result.b is not None:
                    record.samples[trait_id]["ES"] = result.b
                if result.se is not None:
                    record.samples[trait_id]["SE"] = result.se
                if result.nlog_pval is not None:
                    record.samples[trait_id]["LP"] = result.nlog_pval
                if result.alt_freq is not None:
                    record.samples[trait_id]["AF"] = result.alt_freq
                if result.n is not None:
                    record.samples[trait_id]["SS"] = round(result.n)
                if result.imp_z is not None:
                    record.samples[trait_id]["EZ"] = result.imp_z
                if result.imp_info is not None:
                    record.samples[trait_id]["SI"] = result.imp_info
                if result.ncase is not None:
                    record.samples[trait_id]["NC"] = round(result.ncase)
                if result.dbsnpid is not None:
                    record.samples[trait_id]["ID"] = record.id

                # write to file
                vcf.write(record)

        vcf.close()

        # index output file
        logging.info("Indexing output file")
        pysam.tabix_index(path, preset="vcf", force=True, csi=csi)


    """
    Write GWAS file to VCF for multisample input
    Expects an open file handle to a Pickle file of GWAS results & file index dict(chromosome[(position, offset)]) 
    """

    @staticmethod
    def write_to_file_multi_sample(
        gwas_file_dict,
        gwas_idx_dict,
        path,
        fasta,
        build,
        sample_metadata_dict=None,
        file_metadata=None,
        csi=False,
    ):
        logging.info(f"Writing headers to BCF/VCF: {path}")

        header = pysam.VariantHeader()

        # INFO
        header.add_line(
            '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">'
        )
        header.add_line(
            '##INFO=<ID=CR,Number=A,Type=Float,Description="Variant call rate">'
        )
        header.add_line(
            '##INFO=<ID=AM,Number=0,Type=Flag,Description="Variant ambiguous">'
        )
        header.add_line(
            '##INFO=<ID=SR,Number=A,Type=Character,Description="The strand">'
        )
        header.add_line(
            '##INFO=<ID=TY,Number=A,Type=Character,Description="Variant type e.g. Biallelic">'
        )

        # FORMAT
        header.add_line(
            '##FORMAT=<ID=ES,Number=A,Type=Float,Description="Effect size estimate relative to the alternative allele">'
        )
        header.add_line(
            '##FORMAT=<ID=SE,Number=A,Type=Float,Description="Standard error of effect size estimate">'
        )
        header.add_line(
            '##FORMAT=<ID=LP,Number=A,Type=Float,Description="-log10 p-value for effect estimate">'
        )
        header.add_line(
            '##FORMAT=<ID=AF,Number=A,Type=Float,Description="Alternate allele frequency in the association study">'
        )
        header.add_line(
            '##FORMAT=<ID=SS,Number=A,Type=Integer,Description="Sample size used to estimate genetic effect">'
        )
        header.add_line(
            '##FORMAT=<ID=EZ,Number=A,Type=Float,Description="Z-score provided if it was used to derive the EFFECT and SE fields">'
        )
        header.add_line(
            '##FORMAT=<ID=SI,Number=A,Type=Float,Description="Accuracy score of summary data imputation">'
        )
        header.add_line(
            '##FORMAT=<ID=NC,Number=A,Type=Integer,Description="Number of cases used to estimate genetic effect">'
        )
        header.add_line(
            '##FORMAT=<ID=ID,Number=1,Type=String,Description="Study variant identifier">'
        )

        # META
        header.add_line(
            '##META=<ID=TotalVariants,Number=1,Type=Integer,Description="Total number of variants in input">'
        )
        header.add_line(
            '##META=<ID=VariantsNotRead,Number=1,Type=Integer,Description="Number of variants that could not be read">'
        )
        header.add_line(
            '##META=<ID=HarmonisedVariants,Number=1,Type=Integer,Description="Total number of harmonised variants">'
        )
        header.add_line(
            '##META=<ID=VariantsNotHarmonised,Number=1,Type=Integer,Description="Total number of variants that could not be harmonised">'
        )
        header.add_line(
            '##META=<ID=SwitchedAlleles,Number=1,Type=Integer,Description="Total number of variants strand switched">'
        )
        header.add_line(
            '##META=<ID=TotalControls,Number=1,Type=Integer,Description="Total number of controls in the association study">'
        )
        header.add_line(
            '##META=<ID=TotalCases,Number=1,Type=Integer,Description="Total number of cases in the association study">'
        )
        header.add_line(
            '##META=<ID=StudyType,Number=1,Type=String,Description="Type of GWAS study [Continuous or CaseControl]">'
        )

        # CONTIG metadata
        assert len(fasta.references) == len(fasta.lengths)
        for n, contig in enumerate(fasta.references):
            header.add_line(
                f"##contig=<ID={contig},length={fasta.lengths[n]}, assembly={build}>"
            )

        # add file metadata
        if file_metadata is not None:
            for k in file_metadata:
                header.add_line(f"##{k}={file_metadata[k]}")

        # add SAMPLES metadata
        for trait_id in gwas_file_dict.keys():
            header.samples.add(trait_id)
            if sample_metadata_dict is not None:
                meta_string = "".join(f",{k}={sample_metadata_dict[trait_id][k]}" for k in sample_metadata_dict[trait_id])
                header.add_line(f"##SAMPLE=<ID={trait_id}{meta_string}>")

        # create the vcf out
        vcf = pysam.VariantFile(path, "w", header=header)

        # recall variant objects in chromosome position order
        logging.info(f"Writing variants to BCF/VCF: {path}")

        for contig in fasta.references:

            if not any([gwas_i.get(contig) for gwas_i in gwas_idx_dict.values()]):
                continue

            while any([gwas_i.get(contig) for gwas_i in gwas_idx_dict.values()]):

                record = vcf.new_record()

                for trait_id in gwas_file_dict.keys():

                    if not gwas_idx_dict[trait_id]:
                        continue

                    chr_pos = heappop(gwas_idx_dict[trait_id][contig])
                    #print("chr_pos: " + str(chr_pos))

                    # load GWAS result
                    #print("seeking to chr_pos[1]: " + str(chr_pos[1]))
                    gwas_file_dict[trait_id].seek(chr_pos[1])
                    result = pickle.load(gwas_file_dict[trait_id])
                    #print("result: " + str(result))

                    result.nlog_pval = result.nlog_pval
                    #print("result.nlog_pval: " + str(result.nlog_pval))

                    ## PUT BACK CHECKS HERE
                    # check floats
                    if Vcf.is_float32_lossy(result.b):
                        logging.warning(
                            f"Effect field cannot fit into float32. Expect loss of precision for: {result.b}"
                        )

                    if Vcf.is_float32_lossy(result.se):
                        result.se = np.float64(np.finfo(np.float32).tiny).item()
                        logging.warning(
                            f"Standard error field cannot fit into float32. Expect loss of precision for: {result.se}"
                        )
                    if Vcf.is_float32_lossy(result.nlog_pval):
                        logging.warning(
                            f"-log10(pval) field cannot fit into float32. Expect loss of precision for: {result.nlog_pval}"
                        )
                    if Vcf.is_float32_lossy(result.alt_freq):
                        logging.warning(
                            f"Allele frequency field cannot fit into float32. Expect loss of precision for: {result.alt_freq}"
                        )
                    if Vcf.is_float32_lossy(result.imp_z):
                        logging.warning(
                            f"Imputation Z score field cannot fit into float32. Expect loss of precision for: {result.imp_z}"
                        )
                    if Vcf.is_float32_lossy(result.imp_info):
                        logging.warning(
                            f"Imputation INFO field cannot fit into float32. Expect loss of precision for: {result.imp_info}"
                        )

                    record.chrom = result.chrom
                    assert " " not in record.chrom
                    #print("record.chrom: " + str(record.chrom))

                    record.pos = result.pos
                    assert record.pos > 0
                    #print("record.pos: " + str(record.pos))

                    record.id = Vcf.remove_illegal_chars(result.dbsnpid)
                    #print("record.id: " + str(record.id))

                    record.alleles = (result.ref, result.alt)
                    #print("record.alleles: " + str(record.alleles))

                    if result.alt_freq is not None:
                        record.info["AF"] = result.alt_freq
                    #print("result.alt_freq: " + str(result.alt_freq))

                    if result.call_rate is not None:
                        record.info["CR"] = result.call_rate
                    #print("result.call_rate: " + str(result.call_rate))

                    if result.ambig is not None:
                        record.info["AM"] = result.ambig
                    #print("result.ambiguous: " + str(result.ambiguous))

                    if result.strand is not None:
                        record.info["SR"] = result.strand
                    # print("result.strand: " + str(result.strand))

                    if result.vtype is not None:
                        record.info["TY"] = result.vtype
                    #print("result.type: " + str(result.type))

                    if result.b is not None:
                        record.samples[trait_id]["ES"] = result.b
                    #print("result.b: " + str(result.b))

                    if result.se is not None:
                        record.samples[trait_id]["SE"] = result.se
                    #print("result.se: " + str(result.se))

                    if result.nlog_pval is not None:
                        record.samples[trait_id]["LP"] = result.nlog_pval
                    #print("result.nlog_pval: " + str(result.nlog_pval))

                    if result.alt_freq is not None:
                        record.samples[trait_id]["AF"] = result.alt_freq
                    #print("result.alt_freq: " + str(result.alt_freq))

                    if result.n is not None:
                        record.samples[trait_id]["SS"] = round(result.n)
                    #print("result.n: " + str(result.n))

                    if result.imp_z is not None:
                        record.samples[trait_id]["EZ"] = result.imp_z
                    #print("result.imp_z: " + str(result.imp_z))

                    if result.imp_info is not None:
                        record.samples[trait_id]["SI"] = result.imp_info
                    #print("result.imp_info: " + str(result.imp_info))

                    if result.ncase is not None:
                        record.samples[trait_id]["NC"] = round(result.ncase)
                    #print("result.ncase: " + str(result.ncase))

                    if result.dbsnpid is not None:
                        record.samples[trait_id]["ID"] = record.id
                    #print("result.dbsnpid: " + str(result.dbsnpid))

                    # print("record: " + str(record))


                # write to file
                vcf.write(record)

        vcf.close()

        # index output file
        logging.info("Indexing output file")
        pysam.tabix_index(path, preset="vcf", force=True, csi=csi)