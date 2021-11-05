# DADA2 pipeline

### Load libraries

```r
library("tidyverse")
library("gridExtra")
library("dada2")
library("phyloseq")
library("DECIPHER")
library("phangorn")
library("ape")
```

### Set up the environment

```r
cat("Getting ready ...", format(Sys.time(), "%c"), "\n")
n_cores <- 16
indir <- "0_data"
filter_dir <- 'fastq_filtered'
outdir <- 'ASVs'
tax <- 'SILVA'
tax_key <- 'SILVA/silva_nr99_v138.1_wSpecies_train_set.fa.gz'
dir.create(filter_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)
dir.create(tax, showWarnings = FALSE, recursive = TRUE)
download.file("https://zenodo.org/record/4587955/files/silva_nr99_v138.1_wSpecies_train_set.fa.gz", file.path(tax, 'silva_nr99_v138.1_wSpecies_train_set.fa.gz'))
fastqs_raw <- sort(list.files(indir, pattern = '.fastq', full.names = TRUE))
sampleIDs <- sub("*.fastq", "", basename(fastqs_raw))
fastqs_filt <- file.path(filter_dir, paste0(sampleIDs, '.filt.fastq'))
```

### Filter and trim

```r
cat("Filterning and trimming ...", format(Sys.time(), "%c"), "\n")
filter_results <- filterAndTrim(fastqs_raw, fastqs_filt, maxN = 0, truncQ = 2, rm.phix = FALSE, multithread = n_cores, compress = FALSE, verbose = TRUE) 
exists <- file.exists(fastqs_raw) & file.exists(fastqs_filt)
fastqs_filt <- fastqs_filt[exists]
```

### Learn errors

```r
cat("Learning errors ...", format(Sys.time(), "%c"), "\n")
filts <- list.files(filter_dir, pattern=".filt.fastq", full.names=TRUE)
sample.names <- basename(filts)
names(filts) <- sample.names
set.seed(100)
err <- learnErrors(filts, nbases = 1e8, multithread=n_cores, randomize=TRUE)
```

### Dereplicate and infer samples

```r
cat("Dereplicating ...", format(Sys.time(), "%c"), "\n")
dds <- vector("list", length(sample.names))
names(dds) <- sample.names

for(sam in sample.names) {
  cat("Processing of file", sam, "begins on:", format(Sys.time(), "%c"), "\n")
  derep <- derepFastq(filts[[sam]])
  dds[[sam]] <- dada(derep, err=err, multithread=n_cores)
}
```

### Create ASV table

```r
cat("Creating ASV table ...", format(Sys.time(), "%c"), "\n")
seqtab_all <- makeSequenceTable(dds)
```

### Remove chimeras

```r
cat("Removing chimeras ...", format(Sys.time(), "%c"), "\n")
seqtab <- removeBimeraDenovo(seqtab_all, method = 'consensus', multithread = n_cores, verbose = TRUE)
otutable <- seqtab
colnames(otutable) <- paste('ASV', 1:ncol(seqtab), sep = '_')
```

### Assign taxonomy

```r
cat("Assigning taxonomy ...", format(Sys.time(), "%c"), "\n")
df <- as.data.frame(t(seqtab))
n <- 1000
nr <- nrow(df)
a <- split(df, rep(1:ceiling(nr/n), each=n, length.out=nr))
dds.taxa <- vector("list", length(a))
for(i in 1:length(a)) {
  cat("Processing part", i, "of", length(a), "begins on:", format(Sys.time(), "%c"), "\n")
  tmp <- t(a[[i]])
  dds.taxa[[i]] <- assignTaxonomy(tmp, tax_key, multithread = n_cores)
}
taxa <- do.call("rbind", dds.taxa)
colnames(taxa) <- c('Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Genus', "Species")
row.names(taxa) <- paste('ASV', 1:ncol(seqtab), sep = '_')
```

### Save everything

```r
cat("Saving ASVs sequences ...", format(Sys.time(), "%c"), "\n")
asv_seqs <- colnames(seqtab)
asv_headers <- paste('ASV', 1:ncol(seqtab), sep = '_')
asv_fasta <- c(rbind(asv_headers, asv_seqs))
write(asv_fasta, file = file.path(outdir, 'ASVs.fa'))
seqs <- getSequences(seqtab)
names(seqs) <- asv_headers   

cat("Creating phyloseq object ...", format(Sys.time(), "%c"), "\n")

ps <- phyloseq(otu_table(otutable, taxa_are_rows = FALSE), tax_table(taxa))

rm(list=setdiff(ls(), c("ps")))

cat("Saving everything ...", format(Sys.time(), "%c"), "\n")

save.image(file = 'ASV_table.rds')å
```
