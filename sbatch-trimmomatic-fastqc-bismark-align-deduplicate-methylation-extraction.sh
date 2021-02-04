#!/bin/bash 

###############################################################################
#
# important settings for sbatch, bismark alignment to execute
#
###############################################################################
# --cpus-per-task=50
# --time=96:00:00
# --mem=156g
# --gres=lscratch:800
###############################################################################

function do_fastqc_trimmomatic_bismark_align_deduplicate_methylation_extraction () {
	basedir="/data/NHLBI_BCB/Sean_MethylSeq/17_MKJ5433"
	datadir="$basedir/01-fastqs"
	outdir="$basedir/02_methylseq_analysis_pipeline"
	bisulfite_converted_reference="/data/NHLBI_BCB/bin/hg38_Bisulfite_Genome"
	bisulfite_converted_phage="/data/NHLBI_BCB/Sean_MethylSeq/11_set_MKJ5201/03-PhageDNA"
	numcpus=50

	while read SAMPLE READ1 READ2; do
	    echo $SAMPLE,$READ1,$READ2
#		mkdir -p "$outdir/$SAMPLE"
#		sbatch  --job-name="${SAMPLE}" \
#			--partition=norm \
#			--time=96:00:00 \
#			--mem=156g \
#			--cpus-per-task=$numcpus \
#			--gres=lscratch:800 \
#			--error="$outdir/$SAMPLE.slurm.COVID19.methylseq.err.txt" \
#			--output="$outdir/$SAMPLE.slurm.COVID19.methylseq.out.txt" \
sh			$outdir/trimmomatic-fastqc-bismark-align-deduplicate-methylation-extraction.sh $SAMPLE $datadir $READ1 $READ2 $bisulfite_converted_reference $bisulfite_converted_phage $outdir $numcpus
	done < "/data/NHLBI_BCB/Sean_MethylSeq/17_MKJ5433/covid19_batch3_sampleIDs_read1_read2.txt"
}

function do_delete_sam_files () {
	basedir="/data/NHLBI_BCB/Sean_MethylSeq/17_MKJ5433"
	datadir="$basedir/01-fastqs"
	outdir="$basedir/02_methylseq_analysis_pipeline"

	while read SAMPLE READ1 READ2; do
		echo $SAMPLE,$READ1,$READ2,${SAMPLE}.bismark_pe.sam
		rm -f $outdir/${SAMPLE}/bismark_alignment/${SAMPLE}.bismark_pe.sam
	done < "/data/NHLBI_BCB/Sean_MethylSeq/17_MKJ5433/covid19_batch3_sampleIDs_read1_read2.txt"
}

# do_fastqc_trimmomatic_bismark_align_deduplicate_methylation_extraction
do_delete_sam_files
