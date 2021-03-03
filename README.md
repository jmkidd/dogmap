# dogmap
Kidd lab draft pipeline for dog Illumina WGS alignment

**Last updated 3 March 2021**

This is a candidate pipeline and associated file set for consideration by the Dog 10K Project
This approach is not meant to be definitive and is a work in progress.

## Associated files
Associated files for alignment can be found at https://kiddlabshare.med.umich.edu/public-data/UU_Cfam_GSD_1.0-Y/

total 1.7G
```
 87K Feb 26 14:19 canFam4.chromAlias.txt
 144 Feb 26 14:20 chrY.chromAlias.txt
898M Mar  3 13:57 SRZ189891_722g.simp.GSD1.0.vcf.gz
1.6M Mar  3 14:01 SRZ189891_722g.simp.GSD1.0.vcf.gz.tbi
343K Feb 26 14:20 UU_Cfam_GSD_1.0_ROSY.dict
104K Feb 26 14:20 UU_Cfam_GSD_1.0_ROSY.fa.fai
750M Feb 26 14:20 UU_Cfam_GSD_1.0_ROSY.fa.gz
```


The genome reference is based on the German Shepherd assembly from [Wang et. al.](https://www.nature.com/articles/s42003-021-01698-x)
supplemented with Y chromosome sequence from a Labrador retriever [GCF_014441545.1](https://www.ncbi.nlm.nih.gov/assembly/GCF_014441545.1)

Specifically, the UU_Cfam_GSD_1.0 was downloaded from [UCSC](http://hgdownload.soe.ucsc.edu/goldenPath/canFam4/bigZips/).
UCSC naming conventions were utilized.  The following Y chromosome sequences were added:

```
sequenceName accession
chrY_NC_051844.1 NW_024010443.1
chrY_unplaced_NW_024010443.1 NW_024010443.1
chrY_unplaced_NW_024010444.1 NW_024010444.1
```
Chromosome order is as indicated in the genome [.fai file](https://kiddlabshare.med.umich.edu/public-data/UU_Cfam_GSD_1.0-Y/UU_Cfam_GSD_1.0_ROSY.fa.fai)

A set of known variants is required for Base Quality Score Recalibration (BQSR). Mismatches
at the known variant positions are ignored when building the quality recalibration model. For BQSR we utilized
the variant file of 91 million SNV and small indels derived from 772 canines reported in [Plassais et al.](https://www.nature.com/articles/s41467-019-09373-w).

The variant file was converted to UU_Cfam_GSD_1.0 coordinates using the LiftoverVcf from Picard/GATK, resulting in
a total of 71,541,892 SNVs and 16,939,218 indels found in [SRZ189891_722g.simp.GSD1.0.vcf.gz](https://kiddlabshare.med.umich.edu/public-data/UU_Cfam_GSD_1.0-Y/SRZ189891_722g.simp.GSD1.0.vcf.gz).

## Software versions

The following software and versions are used
```
bwa-mem2 version 2.1
gatk version 4.2.0.0
GNU parallel
samtools version >= 1.9
```

## Conceptual overview of pipeline

Pipeline steps are implemented in (process-illumina.py).  This script is designed to run on our
[cluster environment](https://arc-ts.umich.edu/greatlakes/configuration/) which features compute cores with local solid state drives.

All paths are given as examples, and should be modified for your own use.

### Step 1: read alignment using bwa-mem2

```
bwa-mem2 mem -K 100000000  -t NUM_THREADS -Y \
-R READGROUPINFO \
PATH_TO_UU_Cfam_GSD_1.0_ROSY.fa_with_index \
fq1.gz fq2.gz | samtools view -bS - >  /tmp/SAMPLE.bam 
```

This command is inspired by [Reiger et al. Functional equivalence of genome sequencing analysis](https://www.nature.com/articles/s41467-018-06159-4) and by 
discussions with colleagues. The -K option removes non-determinism in framgent length distributions identified using different numbers of threads. -Y uses
soft clipping for supplementary alignments, which may aid down stream structural variation analyses. 

### Step 2: marking duplicates and sorting.

Duplicate marking and sorting is performed in one step using MarkDuplicatesSpark.

```
gatk MarkDuplicatesSpark -I /tmp/SAMPLE.bam  -O /tmp/SAMPLE.sort.md.bam -M /final/SAMPLE.sort.md.metricts.txt \
--tmp-dir /tmp \
--conf 'spark.executor.cores=NUM_THREADS
--conf ''spark.local.dir=/tmp
```
