library(QuasR); library(Rsubread); library(DESeq2)

## upload data
adapter <- scan('adapter.txt', character())
## adapter removal and trimming
# obtain a list of fastq file paths
fastqFiles <- list.files(path = ".", pattern = 'fastq.gz')

# defined processed fastq file names
outfiles <- paste("processed_", fastqFiles,".fastq",sep="")

# process fastq files
# trim 4 bases from the end of the reads (truncateEndBases)
# Remove TGGAATTCTCGGGTGCCAAGG patern if it occurs at the start (Lpattern)
# remove reads shorter than 23 base-pairs (minLength)
preprocessReads(fastqFiles, outfiles, 
                truncateEndBases=4,
                Lpattern=adapter,
                minLength=23)

## align against reference HG38 canonical
processed <- list.files(path = ".", pattern = '.fastq$')

# buildindex(basename="reference_index",reference='GCA_000001405.15_GRCh38_genomic.fna.gz')
# buildindex(basename="hg38_ref",reference='references/GRCh38.primary_assembly.genome.fa.gz')
gc()
# align.stat <- align(index = "./reference_index", annot.inbuilt = 'GCA_000001405.15_GRCh38_genomic.fna.gz', readfile1 = processed, 
#                     output_file =  paste(processed,"subread",output_format='BAM',sep="."), phredOffset = 33)
align.stat <- align(index = "./hg38_ref", annot.inbuilt = 'hg38', readfile1 = processed,
                    output_file =  paste(processed,"subread",output_format='BAM',sep="."), phredOffset = 33, nthreads = -1)
gc()
## extract counts and prepare count matrix
files <- list.files(path = ".", pattern = '.BAM$')
props <- propmapped(files=files) # proportion mapped
sink('mapped_proportions.txt')
props
sink()

fc <- featureCounts(files, annot.inbuilt="hg38", nthreads = -1)
# fc <- featureCounts(files, annot.ext="annotation_files/hg38_gencode.GTF", isGTFAnnotationFile=T, 
#               nthreads=16, allowMultiOverlap=TRUE, fraction=TRUE, isPairedEnd=FALSE)
## DESeq2 Comparison conditions
counts <- as.data.frame(fc$counts)
counts <- counts[]
write.table(counts, 'rawcounts.txt', sep = '\t')

mtdt <- read.table('mtdt.txt', header = T, row.names = 1)
## clean-up samples
dim(counts)
counts <- counts[rowSums(counts[,-1])>0,] ## quitamos todos los reads que sean 0 para todas las muestras
dim(counts)
write.table(counts, 'cleanedcounts.txt', sep = '\t')

dds <- DESeqDataSetFromMatrix(countData = counts, colData = mtdt, design = ~ condition) ## deseq2 normalización y pipeline
dds <- DESeq(dds)
results <- results(dds) ## resultados (log2FC, pvalue y adjusted pvalue)
write.table(results, 'fullresults_deseq2.txt', sep = '\t')

results.significativos <- results[which(results$padj<.05),]
write.table(results.significativos, 'resultados_significativos.txt', sep = '\t')
norm.counts <- counts(dds, normalized = T) ## dataframe con los counts normalizados (ILR)
write.table(norm.counts, 'normalized_counts.txt', sep = '\t')
