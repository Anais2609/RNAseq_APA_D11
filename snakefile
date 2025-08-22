################## Command to launch ##################

# snakemake --config projf="project.json" conff="config.json" -rp -j60 --use-conda

################## Importations ##################

import os
import sys
import errno
import json

################## Functions ##################

def getConditions():
	sampleList = list()
	for cond in config["replicates"].keys():
		for s in config["replicates"][cond]["samples"]:
			sd = dict()
			sd["samplename"] = s
			sd["filename"] = s
			sd["condition"] = cond
			sampleList.append(sd)
	return sampleList

def getFastqFileFW(wildcards):
	for s in config["samples"]:
		if (s["sample"]== wildcards.sampleName):
			return s["files"]["forward"]

def getFastqFileRV(wildcards):
	for s in config["samples"]:
		if (s["sample"]== wildcards.sampleName):
			return s["files"]["reverse"]

################## Configurations variables ##################

configfile: config["projf"]
configfile: config["conff"]

LOGFCTHRESHOLD = 1
FDRTHRESHOLD = 0.05

REFERENCE = "/LAB-DATA/GLiCID/projects/CRCI2NA_DATA/CHILD/genomes/hg38/hg38.analysisSet.fa"
ANNOTATION = "/LAB-DATA/GLiCID/projects/CRCI2NA_DATA/CHILD/genomes/hg38/GCF_000001405.40_GRCh38.p14_genomic_with_chr.gtf"
OUTPUTDIR = "/scratch/nautilus/projects/CRCI2NA_DATA/CHILD/RNAseq_APA_D11"
STAR_INDEX = "/LAB-DATA/GLiCID/projects/CRCI2NA_DATA/CHILD/genomes/genome_index/STAR/genomeDir_hg38_p14"
KALLISTO_INDEX = "/LAB-DATA/GLiCID/projects/CRCI2NA_DATA/CHILD/genomes/genome_index/kallisto/human_refseq_GCF_000001405.40_with_PAX3FOXO1/transcriptome.idx"
TRANSCRIPT_GENE_TAB = "/LAB-DATA/GLiCID/projects/CRCI2NA_DATA/CHILD/genomes/hg38/table_symbol_refseq.tsv"
GSEA_PATH = "/home/postec-a-1@univ-nantes.fr/GSEA_4.3.2/gsea-cli.sh"

sample_all = list()
for s in config["samples"]:
	sample_all.append(s["sample"])

################## Rules ##################

rule all:
	input:
		expand(OUTPUTDIR + "/Samples/{sampleName}/STAR/{sampleName}.Aligned.sortedByCoord.out.bam.bai", sampleName=sample_all),
		OUTPUTDIR + "/multiqc_report.html",
		OUTPUTDIR + "/GSEA_outs/check_GSEA.txt",
		OUTPUTDIR + "/deseq2_kallisto/results_all/PCAplot_with_names.pdf"

rule filter_fastp:
	input: 
		R1 = getFastqFileFW,
		R2 = getFastqFileRV
	output: 
		R1 = temp(OUTPUTDIR + "/Samples/{sampleName}/FASTP/{sampleName}_R1.fastq.gz"),
		R2 = temp(OUTPUTDIR + "/Samples/{sampleName}/FASTP/{sampleName}_R2.fastq.gz"),
		html = OUTPUTDIR + "/Samples/{sampleName}/FASTP/{sampleName}.html",
		json = OUTPUTDIR + "/Samples/{sampleName}/FASTP/{sampleName}.json"
	threads:5
	resources:
		tmpdir=OUTPUTDIR_TMP,
		mem_mb=10000
	shell:
		"""
		guix time-machine -C {OUTPUTDIR}/.guix/fastp_0.23.3/channels.scm -- shell -m {OUTPUTDIR}/.guix/fastp_0.23.3/manifest.scm -- \
		fastp -i {input.R1} -o {output.R1} -I {input.R2} -O {output.R2} --detect_adapter_for_pe -q 30  -j {output.json} --html {output.html} --thread {threads}
	   	"""

rule align_reads:
	input: 
		R1 = OUTPUTDIR + "/Samples/{sampleName}/FASTP/{sampleName}_R1.fastq.gz",
		R2 = OUTPUTDIR + "/Samples/{sampleName}/FASTP/{sampleName}_R2.fastq.gz"
	output:
		OUTPUTDIR + "/Samples/{sampleName}/STAR/{sampleName}.Aligned.sortedByCoord.out.bam"
	threads:30
	resources:
		tmpdir=OUTPUTDIR_TMP,
		parallel_star=1,
		mem_mb=40000
	shell:
		"""
		guix time-machine -C {OUTPUTDIR}/.guix/star_2.7.8a/channels.scm -- shell -m {OUTPUTDIR}/.guix/star_2.7.8a/manifest.scm -- \
		STAR \
		--runThreadN {threads} \
		--genomeDir {STAR_INDEX} \
		--readFilesIn {input.R1} {input.R2} \
		--outFileNamePrefix {OUTPUTDIR}/Samples/{wildcards.sampleName}/STAR/{wildcards.sampleName}. \
		--outSAMtype BAM Unsorted \
		--outBAMsortingThreadN 1 \
		--readFilesCommand zcat \
		--twopassMode Basic

		guix time-machine -C {OUTPUTDIR}/.guix/samtools_1.19/channels.scm -- shell -m {OUTPUTDIR}/.guix/samtools_1.19/manifest.scm -- \
		samtools sort -@ {threads} {OUTPUTDIR}/Samples/{wildcards.sampleName}/STAR/{wildcards.sampleName}.Aligned.out.bam -o {output}
		rm {OUTPUTDIR}/Samples/{wildcards.sampleName}/STAR/{wildcards.sampleName}.Aligned.out.bam
		"""

rule index_bam:
	input:
		OUTPUTDIR + "/Samples/{sampleName}/STAR/{sampleName}.Aligned.sortedByCoord.out.bam"
	output:
		OUTPUTDIR + "/Samples/{sampleName}/STAR/{sampleName}.Aligned.sortedByCoord.out.bam.bai"
	threads:5
	resources:
		tmpdir=OUTPUTDIR_TMP,
		mem_mb=5000
	shell:
		"""
		guix time-machine -C {OUTPUTDIR}/.guix/samtools_1.19/channels.scm -- shell -m {OUTPUTDIR}/.guix/samtools_1.19/manifest.scm -- \
		samtools index -@ {threads} {input}
		"""

rule htseq_count:
	input:
		reads = OUTPUTDIR + "/Samples/{sampleName}/STAR/{sampleName}.Aligned.sortedByCoord.out.bam",
		annotation = ANNOTATION
	output:
		OUTPUTDIR + "/Samples/{sampleName}/{sampleName}_count.txt"
	threads:1
	resources:
		tmpdir=OUTPUTDIR_TMP,
		mem_mb=15000
	shell:
		"""
		guix time-machine -C {OUTPUTDIR}/.guix/htseq_2.0.2/channels.scm -- shell -m {OUTPUTDIR}/.guix/htseq_2.0.2/manifest.scm -- \
		htseq-count \
		--format bam \
		--order pos \
		--quiet \
		-s reverse \
		{input.reads} \
		{input.annotation} \
	   	> {output}
		"""

rule multiqc:
	input:
		expand(OUTPUTDIR + "/Samples/{sampleName}/{sampleName}_count.txt", sampleName=sample_all)
	output:
		OUTPUTDIR + "/multiqc_report.html"
	threads:1
	resources:
		tmpdir=OUTPUTDIR_TMP,
		mem_mb=15000
	shell:
		"""
		guix time-machine -C {OUTPUTDIR}/.guix/multiqc_1.14/channels.scm -- shell -m {OUTPUTDIR}/.guix/multiqc_1.14/manifest.scm -- \
		multiqc -f .
		"""

rule import_count:
	input:
		OUTPUTDIR + "/Samples/{sampleName}/{sampleName}_count.txt"
	output:
		OUTPUTDIR + "/count_dir/{sampleName}.txt"
	threads:1
	resources:
		tmpdir=OUTPUTDIR_TMP,
		mem_mb=5000
	shell:
		"""
		mkdir -p {OUTPUTDIR}/count_dir
		cp {input} {output}
		mkdir -p {OUTPUTDIR}/DESEQ2/results_all
		"""

rule deseq2_conditions:
	input:
		"config.json"
	output: 
		tab = OUTPUTDIR + "/DESEQ2/DESEQ2_conditions.tab"
	run:
		condArray = getConditions()
		with open(output.tab, 'w') as condfile:
			condfile.write("samplename" + "\t" + "filename" + "\t" + "condition" + "\n")
			for s in condArray:
				condfile.write(s["samplename"] + "\t" + s["filename"] + ".txt" + "\t" + s["condition"] + "\n")

rule kallisto:
	input:
		R1 = OUTPUTDIR + "/Samples/{sampleName}/FASTP/{sampleName}_R1.fastq.gz",
		R2 = OUTPUTDIR + "/Samples/{sampleName}/FASTP/{sampleName}_R2.fastq.gz"
	output:
		OUTPUTDIR + "/Samples/{sampleName}/kallisto/abundance.tsv"
	threads:5
	resources:
		tmpdir=OUTPUTDIR_TMP,
		mem_mb=15000
	shell:
		"""
		this_dirname=$(dirname {output})
		guix time-machine -C {OUTPUTDIR}/.guix/kallisto_0.50.1/channels.scm -- shell -m {OUTPUTDIR}/.guix/kallisto_0.50.1/manifest.scm -- \
		kallisto quant -i {KALLISTO_INDEX} -o $this_dirname -b 100 <(zcat {input.R1}) <(zcat {input.R2})
		"""

rule kallisto_conditions:
	input:
		"config.json"
	output: 
		tab = OUTPUTDIR + "/deseq2_kallisto/kallisto_conditions.tab"
	threads:1
	resources:
		tmpdir=OUTPUTDIR_TMP,
		mem_mb=5000
	run:
		condArray = getConditions()
		with open(output.tab, 'w') as condfile:
			condfile.write("samplename" + "\t" + "filename" + "\t" + "condition" + "\n")
			for s in condArray:
				condfile.write(s["samplename"] + "\t" + OUTPUTDIR + "/Samples/" + s["filename"] + "/kallisto/abundance.tsv" + "\t" + s["condition"] + "\n")

rule deseq2:
	input:
		conditions = OUTPUTDIR + "/DESEQ2/DESEQ2_conditions.tab",
		counts = expand(OUTPUTDIR + "/count_dir/{sampleName}.txt", sampleName=sample_all),
		script = OUTPUTDIR + "/r_scripts/runDESeq2.R"
	output: 
		OUTPUTDIR + "/DESEQ2/results_all/PCAplot_with_names.pdf",
		temp(OUTPUTDIR + "/RNAseq_APA_D11_expression_v0.gct")
	threads:5
	resources:
		tmpdir=OUTPUTDIR_TMP,
		mem_mb=15000
	shell:	
		"""
		guix time-machine -C {OUTPUTDIR}/.guix/r_4.4.2_renv_1.0.11/channels.scm -- shell -m {OUTPUTDIR}/.guix/r_4.4.2_renv_1.0.11/manifest.scm -- \
		Rscript {input.script} {input.conditions} {OUTPUTDIR}/count_dir {OUTPUTDIR}/DESEQ2/results_all
		"""

rule deseq2_kallisto:
	input:
		script = OUTPUTDIR + "/r_scripts/runDESeq2_kallisto.R",
		conditions = OUTPUTDIR + "/deseq2_kallisto/kallisto_conditions.tab",
		count_kallisto = expand(OUTPUTDIR + "/Samples/{sampleName}/kallisto/abundance.tsv", sampleName=sample_all)
	output:
		OUTPUTDIR + "/deseq2_kallisto/results_all/PCAplot_with_names.pdf"
	threads:5
	resources:
		tmpdir=OUTPUTDIR_TMP,
		mem_mb=15000
	shell:
		"""
		guix time-machine -C {OUTPUTDIR}/.guix/r_4.4.2_renv_1.0.11/channels.scm -- shell -m {OUTPUTDIR}/.guix/r_4.4.2_renv_1.0.11/manifest.scm -- \
		Rscript {input.script} {input.conditions} {OUTPUTDIR}/deseq2_kallisto/results_all {TRANSCRIPT_GENE_TAB}
		"""

rule convert_gct_for_gsea:
	input:
		OUTPUTDIR + "/RNAseq_APA_D11_expression_v0.gct"
	output:
		OUTPUTDIR + "/RNAseq_APA_D11_expression_v1.gct"
	params:
		nb_samples=len(sample_all)
	threads:1
	resources:
		tmpdir=OUTPUTDIR_TMP,
		mem_mb=10000
	shell:
		"""
		echo -n "#1.2\n43236\t{params.nb_samples}\nNAME\tDESCRIPTION" > {output}
		cat {input} >> {output}
		"""

rule make_cls_for_gsea:
	input:
		OUTPUTDIR + "/RNAseq_APA_D11_expression_v1.gct"
	output:
		OUTPUTDIR + "/samples.cls"
	threads:1
	resources:
		tmpdir=OUTPUTDIR_TMP,
		mem_mb=5000
	shell:
		"""
		read -r -a array_samples <<< $(awk 'NR==3 {{for (i=3; i<=NF; i++) print $i}}' {input} | sed 's/_rep[0-9]*$//' | tr "\n" " ")
		read -r -a array_cond <<< $(awk 'NR==3 {{for (i=3; i<=NF; i++) print $i}}' {input} | sed 's/_rep[0-9]*$//' | uniq | tr "\n" " ")
		echo "${{#array_samples[@]}} ${{#array_cond[@]}} 1" > {output}
		echo "#  ${{array_cond[@]}}" >> {output}
		echo "${{array_samples[@]}}" >> {output}
		"""

rule gsea:
	input:
		gct_file = OUTPUTDIR + "/RNAseq_APA_D11_expression_v1.gct",
		cls_file = OUTPUTDIR + "/samples.cls",
		gmt_file = OUTPUTDIR + "/GO_human_2025_genes_only.gmt"
	output:	
		temp(OUTPUTDIR + "/GSEA_outs/check_GSEA.txt")
	threads:5
	resources:
		tmpdir=OUTPUTDIR_TMP,
		mem_mb=15000
	shell:
		"""
		this_dirname=$(dirname {output})
		guix time-machine -C {OUTPUTDIR}/.guix/openjdk_18.0.2.1/channels.scm -- shell -m {OUTPUTDIR}/.guix/openjdk_18.0.2.1/manifest.scm -- \
		bash {GSEA_PATH} GSEA -res {input.gct_file} -cls {input.cls_file}#RH30_APA_14h_versus_RH30_CTL_14h -gmx {input.gmt_file} -collapse false -permute gene_set -rpt_label RH30_APA_14h_vs_RH30_CTL_14h -out $this_dirname
		guix time-machine -C {OUTPUTDIR}/.guix/openjdk_18.0.2.1/channels.scm -- shell -m {OUTPUTDIR}/.guix/openjdk_18.0.2.1/manifest.scm -- \
		bash {GSEA_PATH} GSEA -res {input.gct_file} -cls {input.cls_file}#RH30_D11_14h_versus_RH30_CTL_14h -gmx {input.gmt_file} -collapse false -permute gene_set -rpt_label RH30_D11_14h_vs_RH30_CTL_14h -out $this_dirname
		guix time-machine -C {OUTPUTDIR}/.guix/openjdk_18.0.2.1/channels.scm -- shell -m {OUTPUTDIR}/.guix/openjdk_18.0.2.1/manifest.scm -- \
		bash {GSEA_PATH} GSEA -res {input.gct_file} -cls {input.cls_file}#RH30_APA_48h_versus_RH30_CTL_48h -gmx {input.gmt_file} -collapse false -permute gene_set -rpt_label RH30_APA_48h_vs_RH30_CTL_48h -out $this_dirname
		guix time-machine -C {OUTPUTDIR}/.guix/openjdk_18.0.2.1/channels.scm -- shell -m {OUTPUTDIR}/.guix/openjdk_18.0.2.1/manifest.scm -- \
		bash {GSEA_PATH} GSEA -res {input.gct_file} -cls {input.cls_file}#RH30_D11_48h_versus_RH30_CTL_48h -gmx {input.gmt_file} -collapse false -permute gene_set -rpt_label RH30_D11_48h_vs_RH30_CTL_48h -out $this_dirname
		echo "GSEA done !" > $this_dirname/check_GSEA.txt
		"""
