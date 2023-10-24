#!/usr/local/bin/perl
use strict;

my $trimmomatic = "/usr/local/devel/ANNOTATION/hlorenzi/Trimmomatic-0.33/trimmomatic-0.33.jar";
my $bowtie = "/usr/local/devel/ANNOTATION/hlorenzi/bin/bowtie2";
my $sortsam = "/usr/local/devel/ANNOTATION/hlorenzi/picard-master/dist/picard.jar SortSam";
my $markduplicates = "/usr/local/devel/ANNOTATION/hlorenzi/picard-master/dist/picard.jar MarkDuplicates";
my $samtools = "/usr/local/bin/samtools";
my $GATK = "/usr/local/devel/ANNOTATION/hlorenzi/bin/GenomeAnalysisTK.jar";
my $java2 = "/usr/local/devel/ANNOTATION/hlorenzi/jre1.7.0_75//bin/java";
my $java = "/usr/local/bin/java";
my $bcftools = "/usr/local/bin/bcftools";
my $filter = "-sLowQual -g3 -G10 -e'%QUAL<10 || (%MAX(DV)<=3 && FMT/GT=\"./1\")|| (%MAX(DV)/%MAX(DP)<=0.3 && FMT/GT=\"./1\") || %MAX(FMT/DP)<5 || (RPB>=0 && RPB<0.1 && %QUAL<15) || (MQB>=0 && MQB<0.05) || (BQB>=0 && BQB<0.05) || (MQSB>=0 && MQSB<0.0)'";
my $snpEff = "/usr/local/devel/ANNOTATION/hlorenzi/snpEff/snpEff.jar";
my $summarize_coding_snps = "/home/hlorenzi/PERL/SCRIPTS/summarize_coding_snps_v2.pl";
my $adapters = "/usr/local/devel/ANNOTATION/hlorenzi/Trimmomatic-0.33/adapters/adapters";

my $usage = "$0 -r <reference_genome> -o <oligos_file> -R <parent/reference fastq prefix> -f <other fastq file prefixes [ONE,TWO,THREE> -F <bcftool filter> -A <reference annotation [GT1/ME49/other path to file]> \n\n";

my %arg = @ARGV;
die $usage unless ($arg{-r} && $arg{-o} && $arg{-R} && $arg{-A});

$filter = $arg{-F} if $arg{-F};
my %annotation = (	"GT1" => "/usr/local/projdata/700030/projects/GCID_PARASITES/KATIE_CLONES/tggt1_names.txt",
			"ME49" => "/usr/local/devel/ANNOTATION/hlorenzi/snpEff/data/ToxoDB-13.0_TgondiiME49/tga4_com_name");
my $annotation_file = $annotation{ $arg{-A} }? $annotation{ $arg{-A} } : $arg{-A};
die "ERROR, I cannot find file $arg{-A}\n" unless (-e $annotation_file);

my @prefix = ($arg{-R});
if ($arg{-f}){
	foreach my $f (split(",",$arg{-f})){
		push  @prefix, $f;
	}
}

## Trim fastq files
my $gz;
my $bowtie_ref_prefix = $arg{-r}; $bowtie_ref_prefix =~ s/\.(fasta|fa|fsa)$//;
foreach my $prefix (@prefix){
	$gz = (-e $prefix.".R2.fastq.gz") ? "fastq.gz" : "fastq";
	next if (-e $prefix.".paired.R2.$gz"); ## skip file if trimmed file already exists
	my $CMD = "$java -jar $trimmomatic PE -threads 10 -phred33 $prefix.R1.$gz $prefix.R2.$gz $prefix.paired.R1.$gz $prefix.unpaired.R1.$gz $prefix.paired.R2.$gz $prefix.unpaired.R2.$gz ILLUMINACLIP:$adapters:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:50"; 
	&print_stderr ("## Trimming reads from $prefix.$gz file",$CMD);

} 
## Map reads to reference
foreach my $prefix (@prefix){
	my $CMD = "$bowtie  --rg-id $prefix --rg PL:Illumina --rg SM:$prefix --rg LB:$prefix -p 10 --phred33 -x $bowtie_ref_prefix -1 $prefix.paired.R1.$gz -2 $prefix.paired.R2.$gz -U $prefix.unpaired.R1.$gz".",$prefix.unpaired.R2.$gz | $samtools view -bSu - > $prefix.bam";
	next if (-e $prefix.".sam"); ## skip file if trimmed file already exists
	#print STDERR "## Mapping reads from $prefix to reference\n$CMD\n\n"; 
	&print_stderr ("## Mapping reads from $prefix to reference",$CMD);
}

## eliminate duplicates from sam
my @files_to_map;
foreach my $prefix (@prefix){
	push @files_to_map, "$prefix.sorted.dedup.realigned.reads.bam";
	
	next if (-e "$prefix.sorted.sam"); ## skip file if trimmed file already exists
	## Sort sam file
	my $CMD = "$java -jar $sortsam INPUT=$prefix.bam OUTPUT=$prefix.sorted.bam SO=coordinate TMP_DIR=/usr/local/scratch/EUK/hlorenzi/TMP/";
	&print_stderr ("## Sorting $prefix.sam file",$CMD) unless (-e "$prefix.sorted.sam");

	next if (-e "$prefix.sorted.dedup.bam"); ## skip file if trimmed file already exists
	## Mark duplicated reads
	my $CMD = "$java -jar $markduplicates INPUT=$prefix.sorted.bam OUTPUT=$prefix.sorted.dedup.bam METRICS_FILE=$prefix.sorted.metrics";
	&print_stderr ("## Markind duplicates in $prefix.sorted.sam file",$CMD) unless (-e "$prefix.sorted.dedup.bam");

	next if (-e "$prefix.sorted.dedup.bam.bai"); ## skip file if trimmed file already exists
	## indexing 
	my $CMD = "$samtools index $prefix.sorted.dedup.bam";
	&print_stderr ("## Markind duplicates in $prefix.sorted.dedup.bam file",$CMD) unless (-e "$prefix.sorted.dedup.bam.bai");

	next if (-e "$prefix.target.intervals.list"); ## skip file if trimmed file already exists
	## Realign reads around indels
	my $CMD = "$java2 -jar $GATK -T RealignerTargetCreator -R $arg{-r} -I $prefix.sorted.dedup.bam -o $prefix.target.intervals.list";
	&print_stderr ("## Realign reads in $prefix.sorted.dedup.bam file",$CMD) unless (-e "$prefix.target.intervals.list");

	next if (-e "$prefix.sorted.dedup.realigned.reads.bam"); ## skip file if trimmed file already exists
	my $CMD = "$java2 -jar $GATK -T IndelRealigner -R $arg{-r} -I $prefix.sorted.dedup.bam -targetIntervals $prefix.target.intervals.list -o $prefix.sorted.dedup.realigned.reads.bam";
	&print_stderr ("## Realign reads in $prefix.sorted.dedup.bam file",$CMD) unless (-e "$prefix.sorted.dedup.realigned.reads.bam");
}



## SNP calls with mpileup
my $p = join("_",@prefix);

## chech if file name will be too long
$p = length($p) > 50 ? "ALL_SAMPLES" : $p;
my $CMD = "$samtools mpileup -o $p.mpileup.bcf -t DP,DPR,DV,DP4,INFO/DPR,SP -g -s -u -f $arg{-r} ".join(" ",@files_to_map);

&print_stderr ("## Running mpileup on @files_to_map",$CMD) unless (-e "$p.mpileup.bcf");

## Bcf call
my $CMD = "$bcftools call -cv -p 0.05 $p.mpileup.bcf > $p.mpileup.vcf";
&print_stderr ("## Running bcftools call on $p.mpileup.bcf",$CMD) unless (-e "$p.mpileup.vcf");

## Bcf filter
my $CMD = "$bcftools filter $filter $p.mpileup.vcf > $p.mpileup.filter.vcf";
&print_stderr ("## Running bcftools filter on $p.mpileup.vcf",$CMD) unless (-e "$p.mpileup.filter.vcf");

## snpEff
my $ref = $arg{-A} eq "GT1" ? "gt1_genome_w_RH_api" : $arg{-A} eq "ME49" ? "ToxoDB-13.0_TgondiiME49" : "OTHER";
if ($ref eq "OTHER" ){
	print STDERR "No available library for $arg{-A}. SnpEff will not run\n\n";
	exit(0);
}
my $CMD = "$java -Xmx4g -jar $snpEff -ud 1000 $ref $p.mpileup.filter.vcf > $p.mpileup.filter.annotation.vcf";
&print_stderr ("## Running snpEff on $p.mpileup.filter.vcf",$CMD) unless (-e "$p.mpileup.filter.annotation.vcf");


my $CMD = "perl $summarize_coding_snps -i $p.mpileup.filter.annotation.vcf -p F -m F -n $annotation_file > $p.mpileup.filter.annotation.SNP.summary.txt";
&print_stderr ("## Running summarize_coding_snps.pl on $p.mpileup.filter.annotation.vcf",$CMD) unless (-e "$p.mpileup.filter.annotation.SNP.summary.txt");

exit(0);

####################################

sub print_stderr {
	my ($message,$cmd) = @_;
	print STDERR "$message\n$cmd\n\n";
	my $err = system("$cmd");
	#`$cmd`;
	if ($err){
		die "ERROR, the commands below failed with error $err\n$cmd\n\n";
	}
}


