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
            '##INFO=<ID=RSID,Number=1,Type=String,Description="dbSNP identifier">'
        )
        header.add_line(
            '##INFO=<ID=AF,Number=A,Type=Float,Description="Allele Frequency">'
        )
        header.add_line(
            '##INFO=<ID=CR,Number=A,Type=Float,Description="Variant call rate">'
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

        # add study metadata
        # if study_metadata is not None:
        #     meta_string = ",".join(f"{k}={study_metadata[k]}" for k in study_metadata)
        #     header.add_line(f"##STUDY=<{meta_string}>")

        # add SAMPLES metadata
        for trait_id in gwas_file_dict.keys():
            header.samples.add(trait_id)
            if sample_metadata_dict is not None:
                meta_string = "".join(
                    f",{k}={sample_metadata_dict[trait_id][k]}" for k in sample_metadata_dict[trait_id])
                header.add_line(f"##SAMPLE=<ID={trait_id}{meta_string}>")

        # create the vcf out
        vcf = pysam.VariantFile(path, "w", header=header)

        # recall variant objects in chromosome position order
        logging.info(f"Writing variants to BCF/VCF: {path}")

        # cycles the reference chromosomes
        for contig in fasta.references:

            # while there is chromosome data in at least one trait heap, process (gets popped off as we go)
            while any([gwas_index.get(contig) for gwas_index in gwas_idx_dict.values()]):

                heap_top_results = []
                # find the top of each trait's heap
                for trait_id in gwas_idx_dict.keys():
                    try:
                        heap_top = gwas_idx_dict[trait_id][contig][0]   # get the top of the heap
                        gwas_file_dict[trait_id].seek(heap_top[1])      # seek the byte offset
                        result = pickle.load(gwas_file_dict[trait_id])  # load the result
                        heap_top_results.append(heap_top + (trait_id,) + (result,))  # (chr_pos, byte, trait_id, result)
                    except (KeyError, IndexError):
                        continue

                # pop the minimum (matching) positions off each heap, ensure the alleles match
                min_res = min(heap_top_results)
                heap_pops = [heappop(gwas_idx_dict[trait_id][contig]) + (trait_id,) + (result, )
                             for chr_pos, byte, trait_id, result in heap_top_results
                             if chr_pos == min_res[0] and
                             result.ref == min_res[3].ref and
                             result.alt == min_res[3].alt]

                # create an empty record
                record = vcf.new_record()

                # process each of the traits into a record (chromosome positions are the same)
                for chr_pos, byte_offset, trait_id, result in heap_pops:

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

                    # TODO: overwriting stuff here, maybe check that each trait has the same info (if it should)
                    # and error if not.
                    # optimise by taking out of the loop?
                    record.chrom = result.chrom
                    assert " " not in record.chrom

                    record.pos = result.pos
                    assert record.pos > 0

                    record.id = Vcf.remove_illegal_chars(result.dbsnpid)

                    record.alleles = (result.ref, result.alt)

                    if result.dbsnpid is not None:
                        record.info["RSID"] = record.id

                    if result.alt_freq is not None:
                        record.info["AF"] = result.alt_freq

                    if result.call_rate is not None:
                        record.info["CR"] = result.call_rate

                    if result.strand is not None:
                        record.info["SR"] = result.strand

                    if result.vtype is not None:
                        record.info["TY"] = result.vtype

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

                # write record to file
                vcf.write(record)

        vcf.close()

        # index output file
        logging.info("Indexing output file")
        pysam.tabix_index(path, preset="vcf", force=True, csi=csi)
