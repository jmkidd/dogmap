# dogmap
Kidd lab draft pipeline for dog Illumina WGS alignment

**Last updated 23 March 2021**

This is a candidate pipeline and associated file set for consideration and discussion by the Dog 10K Project.
This approach is not meant to be definitive and is a work in progress.

## Associated files
Associated files for alignment can be found at https://kiddlabshare.med.umich.edu/public-data/UU_Cfam_GSD_1.0-Y/

```
total 2.4G
  87K Feb 26 14:19 canFam4.chromAlias.txt
  144 Feb 26 14:20 chrY.chromAlias.txt
 898M Mar  3 13:57 SRZ189891_722g.simp.GSD1.0.vcf.gz
 1.6M Mar  3 14:01 SRZ189891_722g.simp.GSD1.0.vcf.gz.tbi
  11M Mar  5 16:00 SRZ189891_722g.simp.header.Axiom_K9_HD.names.GSD_1.0.filter.vcf.gz
 1.2M Mar  5 16:00 SRZ189891_722g.simp.header.Axiom_K9_HD.names.GSD_1.0.filter.vcf.gz.tbi
  11M Mar  5 16:00 SRZ189891_722g.simp.header.CanineHDandAxiom_K9_HD.GSD_1.0.vcf.gz
 1.2M Mar  5 16:00 SRZ189891_722g.simp.header.CanineHDandAxiom_K9_HD.GSD_1.0.vcf.gz.tbi
 2.5M Mar  5 16:00 SRZ189891_722g.simp.header.CanineHD.names.GSD_1.0.filter.vcf.gz
 798K Mar  5 16:00 SRZ189891_722g.simp.header.CanineHD.names.GSD_1.0.filter.vcf.gz.tbi
 680M Mar 23 09:56 UU_Cfam_GSD_1.0.BQSR.DB.bed.gz
 2.0M Mar 23 09:55 UU_Cfam_GSD_1.0.BQSR.DB.bed.gz.tbi
 343K Feb 26 14:20 UU_Cfam_GSD_1.0_ROSY.dict
 104K Feb 26 14:20 UU_Cfam_GSD_1.0_ROSY.fa.fai
 750M Feb 26 14:20 UU_Cfam_GSD_1.0_ROSY.fa.gz
```


The genome reference is based on the German Shepherd assembly from [Wang et. al.](https://www.nature.com/articles/s42003-021-01698-x)
supplemented with Y chromosome sequence from a Labrador retriever [GCF_014441545.1](https://www.ncbi.nlm.nih.gov/assembly/GCF_014441545.1)

Specifically, the UU_Cfam_GSD_1.0 genome assembly was downloaded from [UCSC](http://hgdownload.soe.ucsc.edu/goldenPath/canFam4/bigZips/).
UCSC naming conventions were utilized.  The following Y chromosome sequences were added:

```
sequenceName accession
chrY_NC_051844.1 NW_024010443.1
chrY_unplaced_NW_024010443.1 NW_024010443.1
chrY_unplaced_NW_024010444.1 NW_024010444.1
```
Chromosome order is as indicated in the genome [.fai file](https://kiddlabshare.med.umich.edu/public-data/UU_Cfam_GSD_1.0-Y/UU_Cfam_GSD_1.0_ROSY.fa.fai)

A set of known variants is required for Base Quality Score Recalibration (BQSR). Mismatches
at the known variant positions are ignored when building the recalibration model. We utilize a 
set of SNP and indel positions provided by Reuben Buckley in the Ostander lab. This file is derived
from a set of canine variants that have been hard-filtered using the GATK best practices guidlines.
The variant positions have been liftedOver to UU_Cfam_GSD_1.0 coordinates and are provided in [UU_Cfam_GSD_1.0.BQSR.DB.bed.gz](https://kiddlabshare.med.umich.edu/public-data/UU_Cfam_GSD_1.0-Y/UU_Cfam_GSD_1.0.BQSR.DB.bed.gz)
This file contains a total of 58,325,305 unique SNV and 7,251,525 unique indel positions. Note that there is some redundancy in the file.

We also converted the variant file of 91 million SNV and small indels derived from 772 canines reported in [Plassais et al.](https://www.nature.com/articles/s41467-019-09373-w) The 
variant file was converted to UU_Cfam_GSD_1.0 coordinates using the LiftoverVcf tool from Picard/GATK, resulting in
a total of 71,541,892 SNVs and 16,939,218 indels found in [SRZ189891_722g.simp.GSD1.0.vcf.gz](https://kiddlabshare.med.umich.edu/public-data/UU_Cfam_GSD_1.0-Y/SRZ189891_722g.simp.GSD1.0.vcf.gz).

VQSR, calculation of effective read depth, and calculation of IBS versus existing sample databases
require a set known variant sites.  A collection of known variants has been created for consideration. These files are based
on SNVs genotyped by the Illumina and Axiom genotyping arrays. Successfully unifying data from genotyping
arrays with genome sequence presents many challenges including strand/orientation issues, presence of multiple alleles,
and empirical array performance. The approach taken to deal with these issues is to intersect the variant positions
found on the array with other data, including the variants identified from WGS in [Plassais et al.](https://www.nature.com/articles/s41467-019-09373-w) (note, this
include both PASS and filtered variants from the VQSR employed in that study).

**SRZ189891_722g.simp.header.CanineHD.names.GSD_1.0.filter.vcf.gz** This file contains a subset of sites from the Illumina
CanineHD genotyping array.  The file CanineHD_B.csv was downloaded from the Illumina web page and CanFam3.1 positions were exctracted.
Sites were further filtered for those positions with genotypes reported in [Shannon et al.](https://www.pnas.org/content/112/44/13639)
The resultant set of positions was intersected with the variants from the Plassais et al. 772 canines study, without regard to PASS annotation.  Positions
were then lifted over to UU_Cfam_GSD_1.0 coordinates and further filtered to remove SNVs with more than
two alleles and to remove sites placed on chromosomes other than chr1-chr38 and chrX.  The resultant file 
contains 150,299 SNV positions.

**SRZ189891_722g.simp.header.Axiom_K9_HD.names.GSD_1.0.filter.vcf.gz**  This file contains a subset of 
sites from the Axiom K9 HD array.  Positions were downloaded from the ThermoFisher web page (file: Axiom_K9_HD.na35.r5.a7.annot.csv),
intersected with SNVs from Plassais et al., filtered to remove multi-allelic sites and converted to UU_Cfam_GSD_1.0 coordinates. 
Positions  placed on chromosomes other than chr1-chr38 and chrX were removed.  The resultant file contains 
684,414 SNV postions.

**SRZ189891_722g.simp.header.CanineHDandAxiom_K9_HD.GSD_1.0.vcf**  This file is the union of the sites
on the Illumina CanineHD and the Axiom K9 HD array as described above.  It contains 684,503 sites.

### Recommendations for processing
For BQSR, known positions of variation should be filtered out.  The file `UU_Cfam_GSD_1.0.BQSR.DB.bed.gz` is appropriate for this.

For coverage determination we utilize the effective coverage actually available for SNV
discovery and genotyping. This can be determined by tabulating called coverage at the sites in `SRZ189891_722g.simp.header.CanineHD.names.GSD_1.0.filter.vcf.gz`.
This has the added bonus of producing genotypes at the sites on the Illumina HD array which can be compared with existing
data collections for sample/breed assignment and analysis.

For VQSR, a set of postions highly likely to be truly variable is required. For SNVs, the union set in 
`SRZ189891_722g.simp.header.CanineHDandAxiom_K9_HD.GSD_1.0.vcf.gz` may be appropriate for VQSR training.  A reliable set 
for indels is not yet known.

## Software versions

The following software and versions are used
```
bwa-mem2 version 2.1
gatk version 4.2.0.0
GNU parallel
samtools version >= 1.9
tabix and bgzip from htslib version >= 1.9
```

## Conceptual overview of pipeline

Pipeline steps are implemented in process-illumina.py and process-illumina-file.py.  This script is designed to run on our
[cluster environment](https://arc-ts.umich.edu/greatlakes/configuration/) which features compute cores with local solid state drives.

All paths are given as examples and should be modified for your own use.

**Note: This implementation is designed for use with fast storage.**  If you are running
on a cluster with slower storage, such as network mounted file systems, then steps like MarkDuplicatesSpark will be very slow.
Traditional markduplicates and sorting may work better.

### Step 1: read alignment using bwa-mem2

```
bwa-mem2 mem -K 100000000  -t NUM_THREADS -Y \
-R READGROUPINFO \
PATH_TO_UU_Cfam_GSD_1.0_ROSY.fa_with_index \
fq1.gz fq2.gz | samtools view -bS - >  /tmp/SAMPLE.bam 
```

This command is inspired by [Reiger et al. Functional equivalence of genome sequencing analysis](https://www.nature.com/articles/s41467-018-06159-4) and by 
discussions with colleagues. The -K option removes non-determinism in the fragment length distributions identified using different numbers of threads. -Y uses
soft clipping for supplementary alignments, which may aid down stream structural variation analyses. 

### Step 2: marking duplicates and sorting

Duplicate marking and sorting is performed in one step using MarkDuplicatesSpark.

```
gatk MarkDuplicatesSpark -I /tmp/SAMPLE.bam  \
-O /tmp/SAMPLE.sort.md.bam \
-M /final/SAMPLE.sort.md.metricts.txt \
--tmp-dir /tmp \
--conf 'spark.executor.cores=NUM_THREADS' \
--conf 'spark.local.dir=/tmp'
```

### Step 3: BQSR step 1

BQSR is run in parallel using GNU parallel with 40 seperate jobs for chr1-chr38, chrX, and all other sequences together. The resulting
recalibration data tables are then combined using GatherBQSRReports.

Sample cmd:
```
gatk --java-options "-Xmx4G" BaseRecalibrator \
--tmp-dir /tmp \
-I /tmp/SAMPLE.sort.md.bam \
-R PATH_TO_UU_Cfam_GSD_1.0_ROSY.fa \
--intervals INTERVAL_FILE_TO_PROCESS \
--known-sites PATH_TO_UU_Cfam_GSD_1.0.BQSR.DB.bed.gz \
-O /tmp/bqsr/INTERVAL.recal_data.table
```
Then when completed the recalibration data is combined:

```
gatk --java-options "-Xmx6G" GatherBQSRReports \
--input /tmp/bqsr/INTERVAL-1.recal_data.table \
--input /tmp/bqsr/INTERVAL-2.recal_data.table \
--input /tmp/bqsr/INTERVAL-3.recal_data.table \
...
--output /tmp/bqsrgathered.reports.list
```

### Step 4: BQSR step 2
The base calibration is then performed with ApplyBQSR based on the combined reports. This is done in parallel with
41 separate jobs for chr1-chr38, chrX, all other sequences together, and unmapped reads.
The recalibrated qualities scores are discretized into 3 bins with low quality values unchanged.

Sample cmd:
```
gatk --java-options "-Xmx4G" ApplyBQSR \
--tmp-dir /tmp \
-I /tmp/SAMPLE.sort.md.bam \
-R PATH_TO_UU_Cfam_GSD_1.0_ROSY.fa \
-O /tmp/bqsr/INTERVAL.bqsr.bam \
--intervals INTERVAL_FILE_TO_PROCESS \
--bqsr-recal-file /tmp/bqsrgathered.reports.list \
--preserve-qscores-less-than 6 \
--static-quantized-quals 10 \
--static-quantized-quals 20 \
--static-quantized-quals 30
```
The 41 recalibrated BAMs are then combined using GatherBamFiles

```
gatk --java-options "-Xmx6G" GatherBamFiles \
--CREATE_INDEX true \
-I /tmp/bqsr/INTERVAL-1.bqsr.bam \
-I /tmp/bqsr/INTERVAL-2.bqsr.bam \
-I /tmp/bqsr/INTERVAL-3.bqsr.bam \
...
-O /tmp/SAMPLE.recal.bam
```
### Step 5: Create GVCF file

HaplotypeCaller is then run to create a per-sample GVCF file for subsequent cohort
short variant identification.  This is run in parallel across 39 separate jobs 
for chr1-chr38, chrX, and all other sequences together.  The resulting GVCF files
are then combined using GatherVcfs.

Sample cmd:
```
gatk --java-options "-Xmx4G" HaplotypeCaller \
--tmp-dir /tmp \
-R PATH_TO_UU_Cfam_GSD_1.0_ROSY.fa \
-I /tmp/SAMPLE.recal.bam \
--intervals INTERVAL_FILE_TO_PROCESS \
-O /tmp/GVCF/INTERVAL.g.vcf.gz \
-ERC GVCF 
```

The GVCFs are combined using:
```
gatk --java-options "-Xmx6G" GatherVcfs \
-I /tmp/GVCF/INTERVAL-1.g.vcf.gz \
-I /tmp/GVCF/INTERVAL-2.g.vcf.gz \
-I /tmp/GVCF/INTERVAL-3.g.vcf.gz \
...
-O /tmp/SAMPLE.g.vcf.gz
 ```
The final GCF is then copied to a permanent storage location. 


### Step 6: Convert recalibrated BAM to CRAM

A cram file is created and copied to a permanent storage location

```
gatk --java-options "-Xmx6G" PrintReads \
-R PATH_TO_UU_Cfam_GSD_1.0_ROSY.fa \
-I /tmp/SAMPLE.recal.bam \
-O /permanent/storage/dir/SAMPLE.cram
```

## Use of pipeline script

On our cluster this can be executed automatically as a job that uses multiple processing 
cores on a single compute node. Local fast storage on the node is utilized as tmp space.

Sample cmd:
```
python dogmap/process-illumina.py \
-t 24 \
--sn CH027 \
--lib CH027lib1 \
--fq1 dl/ERR2750983_1.fastq.gz \
--fq2 dl/ERR2750983_2.fastq.gz \
--ref UU_Cfam_GSD_1.0_ROSY.fa \
--refBWA bwa-mem2index/UU_Cfam_GSD_1.0_ROSY.fa \
--tmpdir /tmpssd/CH027tmp \
--finaldir genome-processing/aligned \
--knownsites UU_Cfam_GSD_1.0.BQSR.DB.bed.gz
```

The above has been superseded by process-illumina-file.py which reads in fastq and run information
from a file.  This version can handle multiple fastq/runs/libraries for a single sample.  Following
bwa-mem2 alignment BAMs are merged prior duplicate marking and sorting. 

A sample cmd and input file are shown:

```
python dogmap/process-illumina-file.py \
-t 24 \
--table VILLPT49.runs.table.txt
--ref UU_Cfam_GSD_1.0_ROSY.fa \
--refBWA bwa-mem2index/UU_Cfam_GSD_1.0_ROSY.fa \
--tmpdir /tmpssd/CH027tmp \
--finaldir genome-processing/aligned \
--knownsites UU_Cfam_GSD_1.0.BQSR.DB.bed.gz
```

Where `VILLPT49.runs.table.txt` gives the information for the runs to be processed for the sample. 
The format is
```
sampleName  libraryName readGroupID fastq1  fastq2
```
an example is:
```
cat VILLPT49.runs.table.txt
VILLPT49	VILLPT49_lib	SRR3384031	dl/vd/SRR3384031_1.fastq.gz	dl/vd/SRR3384031_2.fastq.gz
VILLPT49	VILLPT49_lib	SRR3384032	dl/vd/SRR3384032_1.fastq.gz	dl/vd/SRR3384032_2.fastq.gz
```

## Reducing GVCF file size
The default options for generating GVCF.gz files in GATK are optimized for performance. The file size
can be reduced by ~37% by recompressing the .g.vcf.gz files using bgzip with the highest compression
setting.  This can be performed using `recompress-gvcf.py`

## Pipeline completion indicators
An empty file named SAMPLENAME.map.complete is created to indicate 
completion of the alignment pipeline.

An empty file named SAMPLENAME.g.vcf.gz.recompressed is created to indicate 
completion of the GVCF recompression.

These indicator files may aid management of sample processing on your cluster.

## Gathering useful QC metrics
Useful QC metrics, including effective read depth and genotypes as sites on the Illumina HD array, 
can be gathered using the script `run-stats.py`.  This will run several tools
to calculate statistics and compile them in a summary file, SAMPLENAME.cram.stats.txt.  The resulting
stats files can be merged into a single table for analysis and visualization using `combine-stats.py`.


## Known issues and next steps
* Calculate IBS metric from collection of known samples, estimate sample sample identity/breed assignment
* Define PAR regions on X chromosome and genotype males properly
* Identify training set for indels
* Define known copy-number 2 regions for normalization in QuicK-mer2 and fastCN depth-of-coverage based pipelines


## Changes
**17 May 2021** Bug fix in X vs Autosome depth stats report

**23 March 2021** Proposed frozen version for use
This is a proposed frozen version for consideration as part of Dog10K comparison.  There are
several changes:

- Python code refactored
- New driver script process-illumina-file.py handles samples with multiple lanes/libraries.  BAMs are merged using samtools cat with proper header edits
- Recommend using UU_Cfam_GSD_1.0.BQSR.DB.bed.gz for BQSR calibration
- Reduce number of qual score bins kept, in line with https://github.com/CCDG/Pipeline-Standardization/blob/master/PipelineStandard.md.  This results in ~25% reduction in final CRAM file size
- Implement recompression of GVCF files.  This results in ~37% reduction in final g.VCF.gz file size.
- Update run-stats.py to record when there are multiple orientations reported in CollectInsertSizeMetrics.  Reports fracton of read-pairs with the first orientation 
- Switch to bwa-mem2 version 2.1

**8 March 2021** Initial version for testing
