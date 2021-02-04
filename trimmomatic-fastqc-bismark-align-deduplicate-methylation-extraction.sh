#!/bin/bash

set -e        # stop the script if a command fails

function fail {
    echo "FAIL: $@" >&2
    exit 1  # signal failure
}

SAMPLE=$1
DATAPATH=$2
READ1=$3
READ2=$4
bisulfite_converted_refgenome=$5
bisulfite_converted_phagegenome=$6
outdir=$7
numcpus=$8

trimthreads=4
alignthreads=2

# module load trimmomatic
module load trimgalore
module load fastqc
module load bowtie
module load bismark

function do_fastqc () {
	date
	########################################################################################################################
	trimmed=$1

	if [[ "$trimmed" == "no" ]]; 
	then
		DATAPATH=$DATAPATH
		READ1=$READ1
		READ2=$READ2
		out_dir="fastqc"
		#statements
	else
		# if trimming, change DATAPATH (scratch OR saving trimmed fastq files in local directory)
#		DATAPATH="/lscratch/${SLURM_JOBID}"
		DATAPATH="$outdir/$SAMPLE/TrimGalore"

		# if trimming, change $READ1 and $READ2, using Trimmomatic
		#	READ1="${SAMPLE}_1P.fastq.gz"
		#	READ2="${SAMPLE}_2P.fastq.gz"

		# if trimming, change $READ1 and $READ2, using TrimGalore
		READ1="${SAMPLE}_val_1.fq.gz"
		READ2="${SAMPLE}_val_2.fq.gz"

		# if trimming, change $out_dir to something like "fastqc_post_trimming" if you prefer
		out_dir="fastqc_posttrimming"
	fi

#	fastqc -o output_dir [-f fastq|bam|sam] -c contaminant_file seqfile1 .. seqfileN

	mkdir -p $outdir/$SAMPLE/$out_dir

	fastqc -o "$outdir/$SAMPLE/$out_dir"  \
	--nogroup \
	"$DATAPATH/$READ1"  \
	"$DATAPATH/$READ2"  \
	|| fail "fastqc failed"

	echo "fastqc done"
	########################################################################################################################
	date
}

function do_get_fastqc_stats () {
	date
	########################################################################################################################
	trimmed=$1

	if [[ "$trimmed" == "no" ]]; 
	then
		fastqc_READ1=`echo $READ1 | sed 's/.fastq.gz/_fastqc.zip/g'`
		fastqc_READ2=`echo $READ2 | sed 's/.fastq.gz/_fastqc.zip/g'`
		fastqc_READ1_dir=`echo $READ1 | sed 's/.fastq.gz/_fastqc/g'`
		out_dir="fastqc"
		unzip -n $outdir/$SAMPLE/$out_dir/$fastqc_READ1 -d $outdir/$SAMPLE/$out_dir/
		unzip -n $outdir/$SAMPLE/$out_dir/$fastqc_READ2 -d $outdir/$SAMPLE/$out_dir/
		total_sequences=`grep "Total Sequences" $outdir/$SAMPLE/$out_dir/$fastqc_READ1_dir/fastqc_data.txt`
		total_sequences_num=`echo $total_sequences | cut -d" " -f3`
		echo -e $SAMPLE'\t'$total_sequences_num >> "$outdir/total_sequences_num_notrim.txt"
	else
		# if trimming, change $READ1 and $READ2, using TrimGalore
		READ1="${SAMPLE}_val_1.fq.gz"
		READ2="${SAMPLE}_val_2.fq.gz"
		fastqc_READ1=`echo $READ1 | sed 's/.fq.gz/_fastqc.zip/g'`
		fastqc_READ2=`echo $READ2 | sed 's/.fq.gz/_fastqc.zip/g'`
		fastqc_READ1_dir=`echo $READ1 | sed 's/.fq.gz/_fastqc/g'`
		# if trimming, change $out_dir to something like "fastqc_post_trimming" if you prefer
		out_dir="fastqc_posttrimming"
		unzip -n $outdir/$SAMPLE/$out_dir/$fastqc_READ1 -d $outdir/$SAMPLE/$out_dir/
		unzip -n $outdir/$SAMPLE/$out_dir/$fastqc_READ2 -d $outdir/$SAMPLE/$out_dir/
		total_sequences=`grep "Total Sequences" $outdir/$SAMPLE/$out_dir/$fastqc_READ1_dir/fastqc_data.txt`
		total_sequences_num=`echo $total_sequences | cut -d" " -f3`
		echo -e $SAMPLE'\t'$total_sequences_num >> "$outdir/total_sequences_num.txt"
	fi
	########################################################################################################################
	date
}

function do_trimmomatic () {
	date
	########################################################################################################################
	java -jar $TRIMMOJAR PE \
	            -threads $numcpus \
	            "$DATAPATH/$READ1" \
	            "$DATAPATH/$READ2" \
	 	    -baseout "/lscratch/${SLURM_JOBID}/${SAMPLE}.fastq.gz" \
	           ILLUMINACLIP:"/usr/local/apps/trimmomatic/Trimmomatic-0.36/adapters/TruSeq3-PE-2.fa":2:30:10 \
	           MINLEN:50 \
		   HEADCROP:10
	
	# if trimming bases from the start of reads, change HEADCROP
	# HEADCROP:10	\
	echo "trimmomatic done"
	########################################################################################################################
	date
}

function do_trimgalore () {
	date
	########################################################################################################################
	out_dir="TrimGalore"
	mkdir -p $outdir/$SAMPLE/$out_dir

	trim_galore 	--paired \
			--cores $trimthreads \
			--basename ${SAMPLE} \
			--output_dir "$outdir/$SAMPLE/$out_dir" \
			--clip_R1 10 \
			--clip_R2 10 \
			--three_prime_clip_R1 10 \
			--three_prime_clip_R2 10 \
			--length 50 \
			--gzip \
			"$DATAPATH/$READ1" "$DATAPATH/$READ2"

	echo "TrimGalore done"
	########################################################################################################################
	date
#	--output_dir "/lscratch/${SLURM_JOBID}/" \
}

function do_bismark_genome_preparation () {
	numcpus=8
	bismark_genome_preparation	--verbose \
					--parallel $numcpus \
					--single_fasta \
					/data/NHLBI_BCB/Sean_MethylSeq/11_set_MKJ5201/03-PhageDNA

	echo -e "

	DESCRIPTION

	This script is supposed to convert a specified reference genome into two different bisulfite
	converted versions and index them for alignments with Bowtie 2 (default), or HISAT2. The first
	bisulfite genome will have all Cs converted to Ts (C->T), and the other one will have all Gs
	converted to As (G->A). Both bisulfite genomes will be stored in subfolders within the reference
	genome folder. Once the bisulfite conversion has been completed, the program will fork and launch
	two simultaneous instances of the Bowtie 2 or HISAT2 indexer (bowtie2-build or hisat2-build,
	resepctively). Be aware that the indexing process can take up to several hours; this will mainly
	depend on genome size and system resources.

	The new, and still experimental, --slam mode will produce T->C and A->G converted genomes instead.
	The structure of the genome folder will remain the same as for BS-Seq data. This might, or might not,
	change in a future release.

	The following is a brief description of command line options and arguments to control the
	Bismark Genome Preparation:


	USAGE: bismark_genome_preparation [options] <arguments>


	OPTIONS:

	--help                   Displays this help file and exits.

	--version                Displays version information and exits.

	--verbose                Print verbose output for more details or debugging.

	--path_to_aligner </../> The full path to the Bowtie 2 or HISAT2 installation folder on your system
		                 (depending on which aligner/indexer you intend to use; please note that this
		                 is the folder and not any executable). Unless this path is specified, it is
		                 assumed that the aligner in question (Bowtie 2/HISAT2) is in the PATH.

	--bowtie2                This will create bisulfite indexes for use with Bowtie 2. Recommended for most bisulfite
		                 sequencing applications (Default: ON).

	--hisat2                 This will create bisulfite indexes for use with HISAT2. At the time of writing, this is
		                 still unchartered territory, and only recommended for specialist applications such as
		                 RNA-methylation analyses or SLAM-seq type applications (see also: --slam). (Default: OFF).

	--parallel INT           Use several threads for each indexing process to speed up the genome preparation step.
		                 Remember that the indexing is run twice in parallel already (for the top and bottom strand
		                 separately), so e.g. '--parallel 4' will use 8 threads in total. Please also see --large-index
		                 for parallel processing of VERY LARGE genomes (e.g. the axolotl)

	--single_fasta           Instruct the Bismark Indexer to write the converted genomes into
		                 single-entry FastA files instead of making one multi-FastA file (MFA)
		                 per chromosome. This might be useful if individual bisulfite converted
		                 chromosomes are needed (e.g. for debugging), however it can cause a
		                 problem with indexing if the number of chromosomes is vast (this is likely
		                 to be in the range of several thousand files; the operating system can
		                 only handle lists up to a certain length, and some newly assembled
		                 genomes may contain 20000-500000 contigs of scaffold files which do exceed
		                 this list length limit).

	--genomic_composition    Calculate and extract the genomic sequence composition for mono and di-nucleotides
		                 and write the genomic composition table 'genomic_nucleotide_frequencies.txt' to the
		                 genome folder. This may be useful later on when using bam2nuc or the Bismark option
		                 --nucleotide_coverage.
							 
							 
	--slam                   Instead of performing an in-silico bisulfite conversion, this mode transforms T to C (forward strand),
		                 or A to G (reverse strand). The folder structure and rest of the indexing process is currently
		                 exactly the same as for bisulfite sequences, but this might change at some point. This means
		                 that a genome prepared in --slam mode is currently indistinguishable from a true Bisulfite Genome,
		                 so please make sure you name the genome folder appropriately to avoid confusion.

	--large-index            Force generated index to be 'large', even if reference has fewer than 4 billion nucleotides. At the time
		                 of writing this is required for parallel processing of VERY LARGE genomes (e.g. the axolotl)

	ARGUMENTS:

	<path_to_genome_folder>  The path to the folder containing the genome to be bisulfite converted. The Bismark Genome
		                 Preparation expects one or more fastA files in the folder (with the file extension: .fa or
		                 .fasta (also ending in .gz)). Specifying this path is mandatory.


	This script was last modified on 14 April 2019.

	" > /dev/null
}


function do_bismark_align () {
	date
	########################################################################################################################
	out_dir="bismark_alignment"

	mkdir -p $outdir/$SAMPLE/$out_dir

	tempdir="/lscratch/${SLURM_JOB_ID}"

# 		if trimming, change DATAPATH (scratch OR saving trimmed fastq files in local directory)
#		DATAPATH="/lscratch/${SLURM_JOBID}"
		DATAPATH="$outdir/$SAMPLE/TrimGalore"

	# if trimming, change $READ1 and $READ2
#	READ1="${SAMPLE}_1P.fastq.gz"
#	READ2="${SAMPLE}_2P.fastq.gz"

	# if trimming, change $READ1 and $READ2, using Trim Galore
	READ1="${SAMPLE}_val_1.fq.gz"
	READ2="${SAMPLE}_val_2.fq.gz"

	bismark	\
		--bowtie2	\
		-N 1	\
		--output_dir $outdir/$SAMPLE/$out_dir	\
		--bam	\
		-p $alignthreads	\
		-L 22	\
		--score_min L,-0.6,-0.6	\
		--X 1000	\
		--un	\
		--ambiguous	\
		--multicore 8 \
		--temp_dir $tempdir	\
		--genome_folder $bisulfite_converted_refgenome	\
		-1 "$DATAPATH/$READ1"  -2 "$DATAPATH/$READ2"

	echo -e "
	     DESCRIPTION

		The following is a brief description of command line options and arguments to control the Bismark
		bisulfite mapper and methylation caller. Bismark takes in FastA or FastQ files and aligns the
		reads to a specified bisulfite genome. Sequence reads are transformed into a bisulfite converted forward strand
		version (C->T conversion) or into a bisulfite treated reverse strand (G->A conversion of the forward strand).
		Each of these reads are then aligned to bisulfite treated forward strand index of a reference genome
		(C->T converted) and a bisulfite treated reverse strand index of the genome (G->A conversion of the
		forward strand, by doing this alignments will produce the same positions). These 4 instances of Bowtie 2 or HISAT2
		are run in parallel. The sequence file(s) are then read in again sequence by sequence to pull out the original
		sequence from the genome and determine if there were any protected C's present or not.

		The final output of Bismark is in BAM/SAM format by default, described in more detail below.


		USAGE: bismark [options] <genome_folder> {-1 <mates1> -2 <mates2> | <singles>}


		ARGUMENTS:

		<genome_folder>          The path to the folder containing the unmodified reference genome
				         as well as the subfolders created by the Bismark_Genome_Preparation
				         script (/Bisulfite_Genome/CT_conversion/ and /Bisulfite_Genome/GA_conversion/).
				         Bismark expects one or more fastA files in this folder (file extension: .fa, .fa.gz
				         or .fasta or .fasta.gz). The path can be relative or absolute. The path may also be set
				         as '--genome_folder /path/to/genome/folder/'.

		-1 <mates1>              Comma-separated list of files containing the #1 mates (filename usually includes
				         _1), e.g. flyA_1.fq,flyB_1.fq). Sequences specified with this option must
				         correspond file-for-file and read-for-read with those specified in <mates2>.
				         Reads may be a mix of different lengths. Bismark will produce one mapping result
				         and one report file per paired-end input file pair.

		-2 <mates2>              Comma-separated list of files containing the #2 mates (filename usually includes
				         _2), e.g. flyA_1.fq,flyB_1.fq). Sequences specified with this option must
				         correspond file-for-file and read-for-read with those specified in <mates1>.
				         Reads may be a mix of different lengths.

		<singles>                A comma- or space-separated list of files containing the reads to be aligned (e.g.
				         lane1.fq,lane2.fq lane3.fq). Reads may be a mix of different lengths. Bismark will
				         produce one mapping result and one report file per input file. Please note that
				         one should supply a list of files in conjunction with --basename as the output files
				         will constantly overwrite each other...



		OPTIONS:


		Input:

		--se/--single_end <list> Sets single-end mapping mode explicitly giving a list of file names as <list>.
				         The filenames may be provided as a comma [,] or colon [:] separated list.

		-q/--fastq               The query input files (specified as <mate1>,<mate2> or <singles> are FASTQ
				         files (usually having extension .fg or .fastq). This is the default. See also
				         --solexa-quals.

		-f/--fasta               The query input files (specified as <mate1>,<mate2> or <singles> are FASTA
				         files (usually having extensions .fa, .mfa, .fna or similar). All quality values
				         are assumed to be 40 on the Phred scale. FASTA files are expected to contain both
				         the read name and the sequence on a single line (and not spread over several lines).

		-s/--skip <int>          Skip (i.e. do not align) the first <int> reads or read pairs from the input.

		-u/--upto <int>          Only aligns the first <int> reads or read pairs from the input. Default: no limit.

		--phred33-quals          FASTQ qualities are ASCII chars equal to the Phred quality plus 33. Default: ON.

		--phred64-quals          FASTQ qualities are ASCII chars equal to the Phred quality plus 64. Default: off.

		--path_to_bowtie2        The full path </../../> to the Bowtie 2 installation on your system. If not
				         specified it is assumed that Bowtie 2 is in the PATH.

		--path_to_hisat2         The full path </../../> to the HISAT2 installation on your system. If not
				         specified it is assumed that HISAT2 is in the PATH.

		Alignment:


		-N <int>                 Sets the number of mismatches to allowed in a seed alignment during multiseed alignment.
				         Can be set to 0 or 1. Setting this higher makes alignment slower (often much slower)
				         but increases sensitivity. Default: 0. This option is only available for Bowtie 2 (for
				         Bowtie 1 see -n).

		-L <int>                 Sets the length of the seed substrings to align during multiseed alignment. Smaller values
				         make alignment slower but more senstive. Default: the --sensitive preset of Bowtie 2 is
				         used by default, which sets -L to 20. maximum of L can be set to 32. The length of the seed
				         would effect the alignment speed dramatically while the larger L, the faster the aligment.
				         This option is only available for Bowtie 2 (for Bowtie 1 see -l).

		--ignore-quals           When calculating a mismatch penalty, always consider the quality value at the mismatched
				         position to be the highest possible, regardless of the actual value. I.e. input is treated
				         as though all quality values are high. This is also the default behavior when the input
				         doesn't specify quality values (e.g. in -f mode). This option is invariable and on by default.

		-I/--minins <int>        The minimum insert size for valid paired-end alignments. E.g. if -I 60 is specified and
				         a paired-end alignment consists of two 20-bp alignments in the appropriate orientation
				         with a 20-bp gap between them, that alignment is considered valid (as long as -X is also
				         satisfied). A 19-bp gap would not be valid in that case. Default: 0.

		-X/--maxins <int>        The maximum insert size for valid paired-end alignments. E.g. if -X 100 is specified and
				         a paired-end alignment consists of two 20-bp alignments in the proper orientation with a
				         60-bp gap between them, that alignment is considered valid (as long as -I is also satisfied).
				         A 61-bp gap would not be valid in that case. Default: 500.

		--parallel <int>         (May also be --multicore <int>) Sets the number of parallel instances of Bismark to be run concurrently.
				         This forks the Bismark alignment step very early on so that each individual Spawn of Bismark processes
				         only every n-th sequence (n being set by --parallel). Once all processes have completed,
				         the individual BAM files, mapping reports, unmapped or ambiguous FastQ files are merged
				         into single files in very much the same way as they would have been generated running Bismark
				         conventionally with only a single instance.

				         If system resources are plentiful this is a viable option to speed up the alignment process
				         (we observed a near linear speed increase for up to --parallel 8 tested). However, please note
				         that a typical Bismark run will use several cores already (Bismark itself, 2 or 4 threads of
				         Bowtie2/HISAT2, Samtools, gzip etc...) and ~10-16GB of memory depending on the choice of aligner
				         and genome. WARNING: Bismark Parallel (BP?) is resource hungry! Each value of --parallel specified
				         will effectively lead to a linear increase in compute and memory requirements, so --parallel 4 for
				         e.g. the GRCm38 mouse genome will probably use ~20 cores and eat ~40GB or RAM, but at the same time
				         reduce the alignment time to ~25-30%. You have been warned.

		--local                  In this mode, it is not required that the entire read aligns from one end to the other. Rather, some
				         characters may be omitted (“soft-clipped”) from the ends in order to achieve the greatest possible
				         alignment score. For Bowtie 2, the match bonus --ma (default: 2) is used in this mode, and the best possible
				         alignment score is equal to the match bonus (--ma) times the length of the read. This is mutually exclusive with
				         end-to-end alignments. For HISAT2, it is currently not exactly known how the best alignment is calculated.
				         DEFAULT: OFF.


		Output:

		--non_directional        The sequencing library was constructed in a non strand-specific manner, alignments to all four
				         bisulfite strands will be reported. Default: OFF.

				         (The current Illumina protocol for BS-Seq is directional, in which case the strands complementary
				         to the original strands are merely theoretical and should not exist in reality. Specifying directional
				         alignments (which is the default) will only run 2 alignment threads to the original top (OT)
				         or bottom (OB) strands in parallel and report these alignments. This is the recommended option
				         for sprand-specific libraries).

		--pbat                   This options may be used for PBAT-Seq libraries (Post-Bisulfite Adapter Tagging; Kobayashi et al.,
				         PLoS Genetics, 2012). This is essentially the exact opposite of alignments in 'directional' mode,
				         as it will only launch two alignment threads to the CTOT and CTOB strands instead of the normal OT
				         and OB ones. Use this option only if you are certain that your libraries were constructed following
				         a PBAT protocol (if you don't know what PBAT-Seq is you should not specify this option). The option
				         --pbat works only for FastQ files (in both Bowtie and Bowtie 2 mode) and using uncompressed
				         temporary files only).

		--sam-no-hd              Suppress SAM header lines (starting with @). This might be useful when very large input files are
				         split up into several smaller files to run concurrently and the output files are to be merged.

		--rg_tag                 Write out a Read Group tag to the resulting SAM/BAM file. This will write the following line to the
				         SAM header: @RG PL: ILLUMINA ID:SAMPLE SM:SAMPLE ; to set ID and SM see --rg_id and --rg_sample.
				         In addition each read receives an RG:Z:RG-ID tag. Default: OFF.

		--rg_id <string>         Sets the ID field in the @RG header line. The default is 'SAMPLE'.

		--rg_sample <string>     Sets the SM field in the @RG header line; can't be set without setting --rg_id as well. The default is
				         'SAMPLE'.

		-un/--unmapped           Write all reads that could not be aligned to a file in the output directory. Written reads will
				         appear as they did in the input, without any translation of quality values that may have
				         taken place within Bowtie or Bismark. Paired-end reads will be written to two parallel files with _1
				         and _2 inserted in their filenames, i.e. _unmapped_reads_1.txt and unmapped_reads_2.txt. Reads
				         with more than one valid alignment with the same number of lowest mismatches (ambiguous mapping)
				         are also written to _unmapped_reads.txt unless the option --ambiguous is specified as well.

		--ambiguous              Write all reads which produce more than one valid alignment with the same number of lowest
				         mismatches or other reads that fail to align uniquely to a file in the output directory.
				         Written reads will appear as. they did in the input, without any of the translation of quality
				         values that may have taken place within Bowtie or Bismark. Paired-end reads will be written to two
				         parallel files with _1 and _2 inserted in their filenames, i.e. _ambiguous_reads_1.txt and
				         _ambiguous_reads_2.txt. These reads are not written to the file specified with --un.

		-o/--output_dir <dir>    Write all output files into this directory. By default the output files will be written into
				         the same folder as the input file(s). If the specified folder does not exist, Bismark will attempt
				         to create it first. The path to the output folder can be either relative or absolute.

		--temp_dir <dir>         Write temporary files to this directory instead of into the same directory as the input files. If
				         the specified folder does not exist, Bismark will attempt to create it first. The path to the
				         temporary folder can be either relative or absolute.

		--non_bs_mm              Optionally outputs an extra column specifying the number of non-bisulfite mismatches a read during the
				         alignment step. This option is only available for SAM format. In Bowtie 2 context, this value is
				         just the number of actual non-bisulfite mismatches and ignores potential insertions or deletions.
				         The format for single-end reads and read 1 of paired-end reads is 'XA:Z:number of mismatches'
				         and 'XB:Z:number of mismatches' for read 2 of paired-end reads.

		--gzip                   Temporary bisulfite conversion files will be written out in a GZIP compressed form to save disk
				         space. This option is available for most alignment modes but is not available for paired-end FastA
				         files. This option might be somewhat slower than writing out uncompressed files, but this awaits
				         further testing.

		--sam                    The output will be written out in SAM format instead of the default BAM format. Bismark will
				         attempt to use the path to Samtools that was specified with '--samtools_path', or, if it hasn't
				         been specified, attempt to find Samtools in the PATH. If no installation of Samtools can be found,
				         the SAM output will be compressed with GZIP instead (yielding a .sam.gz output file).

		--cram                   Writes the output to a CRAM file instead of BAM. This requires the use of Samtools 1.2 or higher.

		--cram_ref <ref_file>    CRAM output requires you to specify a reference genome as a single FastA file. If this single-FastA
				         reference file is not supplied explicitly it will be regenerated from the genome .fa sequence(s)
				         used for the Bismark run and written to a file called 'Bismark_genome_CRAM_reference.mfa' into the
				         oputput directory.

		--samtools_path          The path to your Samtools installation, e.g. /home/user/samtools/. Does not need to be specified
				         explicitly if Samtools is in the PATH already.

		--prefix <prefix>        Prefixes <prefix> to the output filenames. Trailing dots will be replaced by a single one. For
				         example, '--prefix test' with 'file.fq' would result in the output file 'test.file.fq_bismark.sam' etc.

		-B/--basename <basename> Write all output to files starting with this base file name. For example, '--basename foo'
				         would result in the files 'foo.bam' and 'foo_SE_report.txt' (or its paired-end equivalent). Takes
				         precedence over --prefix. Be advised that you should not use this option in conjunction with supplying
				         lists of files to be processed consecutively, as all output files will constantly overwrite each other.

		--old_flag               Only in paired-end SAM mode, uses the FLAG values used by Bismark v0.8.2 and before. In addition,
				         this options appends /1 and /2 to the read IDs for reads 1 and 2 relative to the input file. Since
				         both the appended read IDs and custom FLAG values may cause problems with some downstream tools
				         such as Picard, new defaults were implemented as of version 0.8.3.


				                             default                         old_flag
				                       ===================              ===================
				                       Read 1       Read 2              Read 1       Read 2

				              OT:         99          147                  67          131

				              OB:         83          163                 115          179

				              CTOT:      147           99                  67          131

				              CTOB:      163           83                 115          179

		--ambig_bam              For reads that have multiple alignments a random alignment is written out to a special file ending in
				         '.ambiguous.bam'. The alignments are in Bowtie2 format and do not any contain Bismark specific
				         entries such as the methylation call etc. These ambiguous BAM files are intended to be used as
				         coverage estimators for variant callers.

		--nucleotide_coverage    Calculates the mono- and di-nucleotide sequence composition of covered positions in the analysed BAM
				         file and compares it to the genomic average composition once alignments are complete by calling 'bam2nuc'.
				         Since this calculation may take a while, bam2nuc attempts to write the genomic sequence composition
				         into a file called 'genomic_nucleotide_frequencies.txt' indside the reference genome folder so it can
				         be re-used the next time round instead of calculating it once again. If a file 'nucleotide_stats.txt' is
				         found with the Bismark reports it will be automatically detected and used for the Bismark HTML report.
				         This option works only for BAM or CRAM files.
								 
		--icpc                   This option will truncate read IDs at the first space or tab it encounters, which are sometimes used to add
				         comments to a FastQ entry (instead of replacing them with underscores (_) as is the default behaviour). The
				         opion is deliberately somewhat cryptic (I couldnt possibly comment), as it only becomes relevant when R1 and R2
				         of read pairs are mapped separately in single-end mode, and then re-paired afterwards (the SAM format dictates
				         that R1 and R2 have the same read ID). Paired-end mapping already creates BAM files with identical read IDs.
				         For more information please see here: https://github.com/FelixKrueger/Bismark/issues/236. Default: OFF.


		OTHER:

		-h/--help                Displays this help file.

		-v/--version             Displays version information.




		BOWTIE 2 SPECIFIC OPTIONS:

		--bowtie2                Default: ON. Uses Bowtie 2 as default aligner. Bismark limits Bowtie 2 to only perform end-to-end
				         alignments, i.e. searches for alignments involving all read characters (also called
				         untrimmed or unclipped alignments). Bismark assumes that raw sequence data is adapter
				         and/or quality trimmed where appropriate. Both small (.bt2) and large (.bt2l) Bowtie 2
				         indexes are supported. To use HISAT2 instead of Bowtie 2 please see option --hisat2.
								 
								 
		--no_dovetail            It is possible, though unusual, for the mates to dovetail, with the mates seemingly extending
				         past each other as in this example:

				         Mate 1:                 GTCAGCTACGATATTGTTTGGGGTGACACATTACGC
				         Mate 2:            TATGAGTCAGCTACGATATTGTTTGGGGTGACACAT
				         Reference: GCAGATTATATGAGTCAGCTACGATATTGTTTGGGGTGACACATTACGCGTCTTTGAC

				         By default, dovetailing is considered inconsistent with concordant alignment, but by default
				         Bismark calls Bowtie 2 with --dovetail, causing it to consider dovetailing alignments as
				         concordant. This becomes relevant whenever reads are clipped from their 5 end prior to mapping,
				         e.g. because of quality or bias issues.

				         Specify --no_dovetail to turn off this behaviour for paired-end libraries. Default: OFF.


		HISAT2 SPECIFIC OPTIONS:


		--hisat2                 Uses HISAT2 instead of Bowtie 2. Bismark uses HISAT2 in end-to-end mode by default,
				         i.e. searches for alignments involving all read characters (also called untrimmed or unclipped alignments)
				         using the option --no-softclipping. Bismark assumes that raw sequence data is adapter and/or quality
				         trimmed where appropriate. From v0.22.0 onwards, Bismark also supports the local alignment mode of
				         HISAT2 (please see --local). Both small (.ht2) and large (.ht2l) HISAT2 indexes are supported. Default: OFF. 

		--no-spliced-alignment   Disable spliced alignment. Default: spliced-alignments are performed.

		--known-splicesite-infile <path>   Provide a list of known splice sites.


		Paired-end options:

		--no-mixed               This option disables the behavior to try to find alignments for the individual mates if
				         it cannot find a concordant or discordant alignment for a pair. This option is invariably on by default.

		--no-discordant          Normally, Bowtie 2 or HISAT2 look for discordant alignments if it cannot find any concordant alignments.
				         A discordant alignment is an alignment where both mates align uniquely, but that does not
				         satisfy the paired-end constraints (--fr/--rf/--ff, -I, -X). This option disables that behavior
				         and it is on by default.



		Bowtie 2 effort options:

		-D <int>                 Up to <int> consecutive seed extension attempts can fail before Bowtie 2 moves on, using
				         the alignments found so far. A seed extension fails if it does not yield a new best or a
				         new second-best alignment. Default: 15.

		-R <int>                 <int> is the maximum number of times Bowtie 2 will re-seed reads with repetitive seeds.
				         When re-seeding, Bowtie 2 simply chooses a new set of reads (same length, same number of
				         mismatches allowed) at different offsets and searches for more alignments. A read is considered
				         to have repetitive seeds if the total number of seed hits divided by the number of seeds
				         that aligned at least once is greater than 300. Default: 2.

		Bowtie 2/ HISAT2 parallelization options:


		-p NTHREADS              Launch NTHREADS parallel search threads (default: 1). Threads will run on separate processors/cores
				         and synchronize when parsing reads and outputting alignments. Searching for alignments is highly
				         parallel, and speedup is close to linear. Increasing -p increases Bowtie 2's memory footprint.
				         E.g. when aligning to a human genome index, increasing -p from 1 to 8 increases the memory footprint
				         by a few hundred megabytes. This option is only available if bowtie is linked with the pthreads
				         library (i.e. if BOWTIE_PTHREADS=0 is not specified at build time). In addition, this option will
				         automatically use the option '--reorder', which guarantees that output SAM records are printed in
				         an order corresponding to the order of the reads in the original input file, even when -p is set
				         greater than 1 (Bismark requires the Bowtie 2 output to be this way). Specifying --reorder and
				         setting -p greater than 1 causes Bowtie 2 to run somewhat slower and use somewhat more memory then
				         if --reorder were not specified. Has no effect if -p is set to 1, since output order will naturally
				         correspond to input order in that case.

		Scoring options:

		--score_min <func>       Sets a function governing the minimum alignment score needed for an alignment to be considered
				         valid (i.e. good enough to report). This is a function of read length. 

				         In end-to-end mode (default), and --local mode for HISAT2 only, --score_min is set as a linear function
				         and is set as <L,value,value>.
				         For instance, specifying L,0,-0.2 sets the minimum-score function f to f(x) = 0 + (-0.2) * x, where x
				         is the read length. The default for end-to-end (global) alignments is: L,0,-0.2.
				         
				         In --local mode for Bowtie 2, the function is logarithmic and is set as <G,value,value>. For instance, specifying
				         G,20,8 sets the minimum-score function f to f(x) = 20 + 8 * ln(x), where x is the read length.
				         The default is for local alignments in Bowtie 2 mode is: G,20,8.
								 
				         See also: setting function options at http://bowtie-bio.sourceforge.net/bowtie2.

		--rdg <int1>,<int2>      Sets the read gap open (<int1>) and extend (<int2>) penalties. A read gap of length N gets a penalty
				         of <int1> + N * <int2>. Default: 5, 3.

		--rfg <int1>,<int2>      Sets the reference gap open (<int1>) and extend (<int2>) penalties. A reference gap of length N gets
				         a penalty of <int1> + N * <int2>. Default: 5, 3.



		Bismark BAM/SAM OUTPUT (default):

		 (1) QNAME  (seq-ID)
		 (2) FLAG   (this flag tries to take the strand a bisulfite read originated from into account (this is different from ordinary DNA alignment flags!))
		 (3) RNAME  (chromosome)
		 (4) POS    (start position)
		 (5) MAPQ   (always 255 for use with Bowtie)
		 (6) CIGAR
		 (7) RNEXT
		 (8) PNEXT
		 (9) TLEN
		(10) SEQ
		(11) QUAL   (Phred33 scale)
		(12) NM-tag (edit distance to the reference)
		(13) MD-tag (base-by-base mismatches to the reference (handles indels)
		(14) XM-tag (methylation call string)
		(15) XR-tag (read conversion state for the alignment)
		(16) XG-tag (genome conversion state for the alignment)
		(17) XA/XB-tag (non-bisulfite mismatches) (optional!)

		Each read of paired-end alignments is written out in a separate line in the above format.


		Last modified on 21 April 2019

	" > /dev/null
	echo "bismark_alignment done"
	########################################################################################################################
	date
}

function do_bismark_align_to_phageDNA () {
	date
	########################################################################################################################
	out_dir="bismark_alignment_to_phageDNA"

	mkdir -p $outdir/$SAMPLE/$out_dir

	tempdir="/lscratch/${SLURM_JOB_ID}"

# 	if trimming, change DATAPATH (scratch OR saving trimmed fastq files in local directory)
#	DATAPATH="/lscratch/${SLURM_JOBID}"
	DATAPATH="$outdir/$SAMPLE/TrimGalore"

	# if trimming, change $READ1 and $READ2
#	READ1="${SAMPLE}_1P.fastq.gz"
#	READ2="${SAMPLE}_2P.fastq.gz"

	# if trimming, change $READ1 and $READ2, using Trim Galore
	READ1="${SAMPLE}_val_1.fq.gz"
	READ2="${SAMPLE}_val_2.fq.gz"

	bismark	\
		--bowtie2	\
		-N 1	\
		--output_dir $outdir/$SAMPLE/$out_dir	\
		--bam	\
		-p $alignthreads	\
		-L 22	\
		--score_min L,-0.6,-0.6	\
		--X 1000	\
		--un	\
		--ambiguous	\
		--multicore 8 \
		--temp_dir $tempdir	\
		--genome_folder $bisulfite_converted_phagegenome	\
		-1 "$DATAPATH/$READ1"  -2 "$DATAPATH/$READ2"

	echo "bismark_alignment to PhageDNA done"
	########################################################################################################################
	date
}

function do_bismark_deduplicate () {
	date
	########################################################################################################################
	out_dir="bismark_deduplicated"
	mkdir -p $outdir/$SAMPLE/$out_dir
		
	bismark_align_out_dir="bismark_alignment"

	samtools view -h $outdir/$SAMPLE/$bismark_align_out_dir/${SAMPLE}_*.bam > $outdir/$SAMPLE/$bismark_align_out_dir/$SAMPLE.bismark_pe.sam

	cd $outdir/$SAMPLE/$out_dir

	deduplicate_bismark -p --bam "$outdir/$SAMPLE/$bismark_align_out_dir/$SAMPLE.bismark_pe.sam"
	echo "bismark_deduplication done"
	########################################################################################################################
	date
}

function do_bismark_methylation_extract () {
	date
	########################################################################################################################
	out_dir="bismark_methylation_deduplicated_extracted"
	mkdir -p $outdir/$SAMPLE/$out_dir

	bismark_align_deduplicated_out_dir="bismark_deduplicated"

	bismark_methylation_extractor	\
		--paired-end 	\
		--no_overlap	\
		--report	\
		--bedGraph	\
		--counts	\
		--buffer_size 100G	\
		--no_header	\
		--cytosine_report	\
		--output "$outdir/$SAMPLE/$out_dir"	\
		--multicore $numcpus	\
		--gzip	\
		--scaffolds	\
		--genome_folder	$bisulfite_converted_refgenome	\
		"$outdir/$SAMPLE/$bismark_align_deduplicated_out_dir/$SAMPLE.bismark_pe.deduplicated.bam"

	echo -e "
	DESCRIPTION

	The following is a brief description of all options to control the Bismark
	methylation extractor. The script reads in a bisulfite read alignment results file
	produced by the Bismark bisulfite mapper and extracts the methylation information
	for individual cytosines. This information is found in the methylation call field
	which can contain the following characters:

	       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
	       ~~~   X   for methylated C in CHG context                      ~~~
	       ~~~   x   for not methylated C CHG                             ~~~
	       ~~~   H   for methylated C in CHH context                      ~~~
	       ~~~   h   for not methylated C in CHH context                  ~~~
	       ~~~   Z   for methylated C in CpG context                      ~~~
	       ~~~   z   for not methylated C in CpG context                  ~~~
	       ~~~   U   for methylated C in Unknown context (CN or CHN       ~~~
	       ~~~   u   for not methylated C in Unknown context (CN or CHN)  ~~~
	       ~~~   .   for any bases not involving cytosines                ~~~
	       ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

	The methylation extractor outputs result files for cytosines in CpG, CHG and CHH
	context (this distinction is actually already made in Bismark itself). As the methylation
	information for every C analysed can produce files which easily have tens or even hundreds of
	millions of lines, file sizes can become very large and more difficult to handle. The C
	methylation info additionally splits cytosine methylation calls up into one of the four possible
	strands a given bisulfite read aligned against:

		     OT      original top strand
		     CTOT    complementary to original top strand

		     OB      original bottom strand
		     CTOB    complementary to original bottom strand

	Thus, by default twelve individual output files are being generated per input file (unless
	--comprehensive is specified, see below). The output files can be imported into a genome
	viewer, such as SeqMonk, and re-combined into a single data group if desired (in fact
	unless the bisulfite reads were generated preserving directionality it doesn't make any
	sense to look at the data in a strand-specific manner). Strand-specific output files can
	optionally be skipped, in which case only three output files for CpG, CHG or CHH context
	will be generated. For both the strand-specific and comprehensive outputs there is also
	the option to merge both non-CpG contexts (CHG and CHH) into one single non-CpG context.


	The output files are in the following format (tab delimited):

	<sequence_id>     <strand>      <chromosome>     <position>     <methylation call>


	USAGE: methylation_extractor [options] <filenames>


	ARGUMENTS:
	==========

	<filenames>              A space-separated list of Bismark result files in SAM format from
		                 which methylation information is extracted for every cytosine in
		                 the reads. For alignment files in the older custom Bismark output
		                 see option '--vanilla'.

	OPTIONS:

	-s/--single-end          Input file(s) are Bismark result file(s) generated from single-end
		                 read data. Specifying either --single-end or --paired-end is
		                 mandatory.

	-p/--paired-end          Input file(s) are Bismark result file(s) generated from paired-end
		                 read data. Specifying either --paired-end or --single-end is
		                 mandatory.

	--vanilla                The Bismark result input file(s) are in the old custom Bismark format
		                 (up to version 0.5.x) and not in SAM format which is the default as
		                 of Bismark version 0.6.x or higher. Default: OFF.

	--no_overlap             For paired-end reads it is theoretically possible that read_1 and
		                 read_2 overlap. This option avoids scoring overlapping methylation
		                 calls twice (only methylation calls of read 1 are used for in the process
		                 since read 1 has historically higher quality basecalls than read 2).
		                 Whilst this option removes a bias towards more methylation calls
		                 in the center of sequenced fragments it may de facto remove a sizable
		                 proportion of the data. This option is highly recommended for paired-end
		                 data.

	--ignore <int>           Ignore the first <int> bp from the 5' end of Read 1 when processing the
		                 methylation call string. This can remove e.g. a restriction enzyme site
		                 at the start of each read or any other source of bias (e.g. PBAT-Seq data).

	--ignore_r2 <int>        Ignore the first <int> bp from the 5' end of Read 2 of paired-end sequencing
		                 results only. Since the first couple of bases in Read 2 of BS-Seq experiments
		                 show a severe bias towards non-methylation as a result of end-repairing
		                 sonicated fragments with unmethylated cytosines (see M-bias plot), it is
		                 recommended that the first couple of bp of Read 2 are removed before
		                 starting downstream analysis. Please see the section on M-bias plots in the
		                 Bismark User Guide for more details.

	--comprehensive          Specifying this option will merge all four possible strand-specific
		                 methylation info into context-dependent output files. The default

		                 contexts are:
		                  - CpG context
		                  - CHG context
		                  - CHH context

	--merge_non_CpG          This will produce two output files (in --comprehensive mode) or eight
		                 strand-specific output files (default) for Cs in
		                  - CpG context
		                  - non-CpG context

	--report                 Prints out a short methylation summary as well as the paramaters used to run
		                 this script.

	--no_header              Suppresses the Bismark version header line in all output files for more convenient
		                 batch processing.

	-o/--output DIR          Allows specification of a different output directory (absolute or relative
		                 path). If not specified explicitely, the output will be written to the current directory.

	--samtools_path          The path to your Samtools installation, e.g. /home/user/samtools/. Does not need to be specified
		                 explicitly if Samtools is in the PATH already.

	--gzip                   The methylation extractor files (CpG_OT_..., CpG_OB_... etc) will be written out in
		                 a GZIP compressed form to save disk space. This option does not work on bedGraph and
		                 genome-wide cytosine reports as they are 'tiny' anyway.

	--version                Displays version information.

	-h/--help                Displays this help file and exits.

	--mbias_only             The methylation extractor will read the entire file but only output the M-bias table and plots as
		                 well as a report (optional) and then quit. Default: OFF.



	bedGraph specific options:
	==========================

	--bedGraph               After finishing the methylation extraction, the methylation output is written into a
		                 sorted bedGraph file that reports the position of a given cytosine and its methylation
		                 state (in %, see details below). The methylation extractor output is temporarily split up into
		                 temporary files, one per chromosome (written into the current directory or folder
		                 specified with -o/--output); these temp files are then used for sorting and deleted
		                 afterwards. By default, only cytosines in CpG context will be sorted. The option
		                 '--CX_context' may be used to report all cytosines irrespective of sequence context
		                 (this will take MUCH longer!). The default folder for temporary files during the sorting
		                 process is the output directory. The bedGraph conversion step is performed by the external
		                 module 'bismark2bedGraph'; this script needs to reside in the same folder as the
		                 bismark_methylation_extractor itself.


	--cutoff [threshold]     The minimum number of times a methylation state has to be seen for that nucleotide
		                 before its methylation percentage is reported. Default: 1.

	--remove_spaces          Replaces whitespaces in the sequence ID field with underscores to allow sorting.


	--CX/--CX_context        The sorted bedGraph output file contains information on every single cytosine that was covered
		                 in the experiment irrespective of its sequence context. This applies to both forward and
		                 reverse strands. Please be aware that this option may generate large temporary and output files
		                 and may take a long time to sort (up to many hours). Default: OFF.
		                 (i.e. Default = CpG context only).

	--buffer_size <string>   This allows you to specify the main memory sort buffer when sorting the methylation information.
		                 Either specify a percentage of physical memory by appending % (e.g. --buffer_size 50%) or
		                 a multiple of 1024 bytes, e.g. 'K' multiplies by 1024, 'M' by 1048576 and so on for 'T' etc.
		                 (e.g. --buffer_size 20G). For more information on sort type 'info sort' on a command line.
		                 Defaults to 2G.

	--scaffolds/--gazillion  Users working with unfinished genomes sporting tens or even hundreds of thousands of
		                 scaffolds/contigs/chromosomes frequently encountered errors with pre-sorting reads to
		                 individual chromosome files. These errors were caused by the operating system's limit
		                 of the number of filehandle that can be written to at any one time (typically 1024; to
		                 find out this limit on Linux, type: ulimit -a).
		                 To bypass the limitation of open filehandles, the option --scaffolds does not pre-sort
		                 methylation calls into individual chromosome files. Instead, all input files are
		                 temporarily merged into a single file (unless there is only a single file), and this
		                 file will then be sorted by both chromosome AND position using the Unix sort command.
		                 Please be aware that this option might take a looooong time to complete, depending on
		                 the size of the input files, and the memory you allocate to this process (see --buffer_size).
		                 Nevertheless, it seems to be working.

	--ample_memory           Using this option will not sort chromosomal positions using the UNIX 'sort' command, but will
		                 instead use two arrays to sort methylated and unmethylated calls. This may result in a faster
		                 sorting process of very large files, but this comes at the cost of a larger memory footprint
		                 (two arrays of the length of the largest human chromosome 1 (~250M bp) consume around 16GB
		                 of RAM). Due to overheads in creating and looping through these arrays it seems that it will
		                 actually be *slower* for small files (few million alignments), and we are currently testing at
		                 which point it is advisable to use this option. Note that --ample_memory is not compatible
		                 with options '--scaffolds/--gazillion' (as it requires pre-sorted files to begin with).



	Genome-wide cytosine methylation report specific options:
	=========================================================

	--cytosine_report        After the conversion to bedGraph has completed, the option '--cytosine_report' produces a
		                 genome-wide methylation report for all cytosines in the genome. By default, the output uses 1-based
		                 chromosome coordinates (zero-based cords are optional) and reports CpG context only (all
		                 cytosine context is optional). The output considers all Cs on both forward and reverse strands and
		                 reports their position, strand, trinucleotide content and methylation state (counts are 0 if not
		                 covered). The cytsoine report conversion step is performed by the external module
		                 'bedGraph2cytosine'; this script needs to reside in the same folder as the bismark_methylation_extractor
		                 itself.

	--CX/--CX_context        The output file contains information on every single cytosine in the genome irrespective of
		                 its context. This applies to both forward and reverse strands. Please be aware that this will
		                 generate output files with > 1.1 billion lines for a mammalian genome such as human or mouse.
		                 Default: OFF (i.e. Default = CpG context only).

	--zero_based             Uses zero-based coordinates like used in e.g. bed files instead of 1-based coordinates. Default: OFF.

	--genome_folder <path>   Enter the genome folder you wish to use to extract sequences from (full path only). Accepted
		                 formats are FastA files ending with '.fa' or '.fasta'. Specifying a genome folder path is mandatory.

	--split_by_chromosome    Writes the output into individual files for each chromosome instead of a single output file. Files
		                 will be named to include the input filename and the chromosome number.



	OUTPUT:

	The bismark_methylation_extractor output is in the form:
	========================================================
	<seq-ID>  <methylation state*>  <chromosome>  <start position (= end position)>  <methylation call>

	* Methylated cytosines receive a '+' orientation,
	* Unmethylated cytosines receive a '-' orientation.



	The bedGraph output (optional) looks like this (tab-delimited; 0-based start coords, 1-based end coords):
	=========================================================================================================

	track type=bedGraph (header line)

	<chromosome>  <start position>  <end position>  <methylation percentage>



	The coverage output looks like this (tab-delimited, 1-based genomic coords):
	============================================================================

	<chromosome>  <start position>  <end position>  <methylation percentage>  <count methylated>  <count non-methylated>



	The genome-wide cytosine methylation output file is tab-delimited in the following format:
	==========================================================================================
	<chromosome>  <position>  <strand>  <count methylated>  <count non-methylated>  <C-context>  <trinucleotide context>

	" > /dev/null
	echo "bismark_methylation_extraction done"
	########################################################################################################################
	date
}

# do_fastqc no
# do_trimmomatic
# do_trimgalore
# do_fastqc yes
# do_get_fastqc_stats no
do_get_fastqc_stats yes
# do_bismark_genome_preparation
# do_bismark_align
# do_bismark_align_to_phageDNA
# do_bismark_deduplicate
# do_bismark_methylation_extract
