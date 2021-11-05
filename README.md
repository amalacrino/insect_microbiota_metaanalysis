# Host phylogeny drives the structure and diversity of insect microbiota  

### Antonino Malacrinò

#### *Institute for Evolution and Biodiversity, Westfälische Wilhelms-Universität Münster, Münster, Germany*

## Abstract 

As for most of the life that inhabits our planet, microorganisms play an essential role in the fitness of insects, including nutrition, reproduction, defence and many other functions. More recently, we assisted to an exponential growth of studies describing the taxonomical composition of bacterial communities across insects’ phylogeny. However, there is still an outstanding question that needs to be answered: which factors contribute most in shaping insects’ microbiomes? This study tries to find an answer to this question by taking advantage of publicly available sequencing data and reanalysing over 4,000 samples of insect-associated bacterial communities under a common framework. Results suggest that insect taxonomy has a wider impact on the structure and diversity of their associated microbial communities than the other factors considered (diet, sex, life stage, sample origin and treatment). However, when specifically testing for signatures of co-diversification of insect species and their microbiota, analyses found a weak support for this, suggesting that while insect species strongly drives the structure and diversity of insect microbiota, the diversification of those microbial communities did not follow their host’s phylogeny. Furthermore, a parallel survey of the literature highlights several methodological limitations that needs to be considered in future research endeavours.


## 1. Data

SRA accession number are provided in the file SRA_Acc_List.txt

## 2. Merge PE reads

PE reads were merged using FLASH 1.2.11 (Magoč and Salzberg, 2011)

```bash
for i in *_1.fastq; do FILENAME=`basename $i _1.fastq`; echo $FILENAME; flash ${FILENAME}_1.fastq ${FILENAME}_2.fastq --output-prefix=${FILENAME}.merged; done

find . ! -name '*.merged.extendedFrags.fastq' -delete
for file in *; do mv ${file} ${file/.merged.extendedFrags/}; done
```

## 3a. Process data with VSEARCH

Data handling was carried out using VSEARCH (Rognes et al. 2016), assigning taxonomy using SILVA v132 database (Quast et al. 2012). 

Run the pipeline according to https://github.com/torognes/vsearch/wiki/VSEARCH-pipeline and assign taxonomy using the command

```bash
vsearch --sintax $PRJ_DIR/otus.fasta --db $PRJ_DIR/silva_16s_v123.fa --tabbedout $PRJ_DIR/ASV_tax_raw.txt --sintax_cutoff 0.5 --threads $NTHREADS
```

## 3b. Process data with DADA2

The pipeline used to process data with DADA2 (Callahan et al., 2016) is reported here.

## 4. Data analysis

Load packages

```R
library("dplyr")
library("phyloseq") 
library("vegan")
library("data.table") 
library("car")
library("lme4")
library("parallel")
library("tidyr")
library("stringr")
library("stringi")
library("ape")
library("ggplot2")
library("MuMIn")
library("emmeans")
library("ggpubr")
```

Load data

```R
data <- read.table("data.txt", header = T, sep = "\t", row.names = 1)
map <- read.table("metadata.txt", header = T, sep = "\t")
taxa <- read.table("taxonomy.txt", header = F, sep = ";")
MicroDF <- merge_phyloseq(data, map, taxa)
```

Clean database

```R
MicroDF <- subset_taxa(MicroDF, Class !="Chloroplast") 
MicroDF <- subset_taxa(MicroDF, Order !="Chloroplast")
MicroDF <- filter_taxa(MicroDF, function (x) {sum(x > 0) > 1}, prune=TRUE)
MicroDF <- prune_samples(sample_sums(MicroDF) >= 5000, MicroDF)

GM.whole <- subset_samples(GM, Sample_tissue == "whole" | Sample_tissue == "surface_sterilized")
GM.gut <- subset_samples(GM, Sample_tissue == "gut")
```

Test the influence of factors on the diversity of insect-associated microbial communities

```R
calculate.shannon.div <- function(x){
  diversity <- estimate_richness(x, split = TRUE, measures = "Shannon")
  div <- cbind(sample_data(x), diversity)
  return(div)
}

model.shannon.div <- function(x){
  shannon.model <- lmer(Shannon ~ Insect_ID + Sex + Life_stage + Diet_study + Sample_treatment + Sample_origin + (1|Study_ID) + (Insect_Order/Insect_ID), data = x)
  shannon.anova <- Anova(shannon.model)
  shannon.anova <- data.table::setDT(shannon.anova, keep.rownames = TRUE)[]
  r2lmer <- data.frame(Factor = c("Insect_ID", "Sex", "Life_stage", "Diet_study", "Sample_treatment", "Sample_origin"),
                       R2 = c(r.squaredGLMM(lmer(Shannon ~ Insect_ID + (1|Study_ID), data = x))[1],
                             r.squaredGLMM(lmer(Shannon ~ Sex + (1|Study_ID), data = x))[1],
                             r.squaredGLMM(lmer(Shannon ~ Life_stage + (1|Study_ID), data = x))[1],
                             r.squaredGLMM(lmer(Shannon ~ Diet_study + (1|Study_ID), data = x))[1],
                             r.squaredGLMM(lmer(Shannon ~ Sample_treatment + (1|Study_ID), data = x))[1],
                             r.squaredGLMM(lmer(Shannon ~ Sample_origin + (1|Study_ID), data = x))[1]))
  model.sum <- merge(shannon.anova, r2lmer, by.x = "rn", by.y = "Factor")
  return(model.sum)
}

div.w <- calculate.shannon.div(GM.whole)
mod.shannon.w <- model.shannon.div(div.w)
div.g <- calculate.shannon.div(GM.gut)
mod.shannon.g <- model.shannon.div(div.g)
```

Test the influence of factors on the structure of insect-associated microbial communities

```R
calculate.permanova <- function(x, y){
  sampledf <- data.frame(sample_data(x))
  perm <- how(nperm = 999)
  set.seed(100)
  setBlocks(perm) <- with(sampledf, Study_ID)
  dist.mat <- phyloseq::distance(x, method = y)
  pmv <- adonis2(dist.mat ~ Insect_ID + Sex + Life_stage + Diet_study + Sample_treatment + Sample_origin, data = sampledf, permutations = perm, parallel = NTHREADS)
  return(pmv)
}

permanova.bray.w <- calculate.permanova(GM.whole, "bray")
permanova.jacc.w <- calculate.permanova(GM.whole, "jaccard")
permanova.bray.g <- calculate.permanova(GM.gut, "bray")
permanova.jacc.g <- calculate.permanova(GM.gut, "jaccard")
```

Infer phylosymbiosis (only insect guts)

```R
ref_tree <- ape::read.tree('Insecta_species.nwk')
list.tree <- as.vector(ref_tree$tip.label)
list.gm.g <- sample_data(GM.gut)$Insect_ID
common.list.g <- intersect(list.gm.g, list.tree)
merged.GM.g <- subset_samples(GM.gut, Insect_ID %in% common.list.g)

generate.filt.tree <- function(x){
  mergedGM.tree <- merge_samples(x, "Insect_ID")
  dist.mat.bray.tree <- phyloseq::distance(mergedGM.tree, method = "bray")
  dist.mat.bray.tree <- as.matrix(dist.mat.bray.tree)
  tree.filt <- picante::prune.sample(phylo = ref_tree, samp = dist.mat.bray.tree)
  return(tree.filt)
}

filt.tree.g <- generate.filt.tree(merged.GM.g)

tip.age <- function(x){
  data.tree <- setNames(x$edge.length[sapply(1:length(x$tip.label),function(x,y) which (y==x),y=x$edge[,2])],x$tip.label) 
  data.tree <- as.data.frame(data.tree)
  data.tree <- data.table::setDT(data.tree, keep.rownames = TRUE)[]
}

tip.age.g <- tip.age(filt.tree.g)
dv.g <- merge(div.g, tip.age.g, by.x = "Insect_ID", by.y = "rn")
tree.g <- ape::cophenetic.phylo(filt.tree.g)

### Blomberg’s K
lngest = dv.g$Shannon
names(lngest) = dv.g$Insect.ID
set.seed(100)
picante::phylosignal(lngest[filt.tree.g$tip.label], filt.tree.g, reps = 9999)

### Pagel's Lambda
p.lambda <- comparative.data(phy = filt.tree.g, data = dv.g, names.col = Insect.ID, vcv = TRUE, na.omit = FALSE, warn.dropped = TRUE)
model.lambda <- pgls(Shannon ~ 1, data = p.lambda, lambda='ML')
summary(model.lambda)

### Mantel test
dist.mat.g.a <- generate.dist.mat(merged.GM.g)
mantel.g <- mantel(tree.g, dist.mat.g.a, method = "pearson", permutations = 9999, na.rm = TRUE)
mantel.g
```

## References

Callahan, B. J., McMurdie, P. J., Rosen, M. J., Han, A. W., Johnson, A. J. A., & Holmes, S. P. (2016). DADA2: high-resolution sample inference from Illumina amplicon data. Nature methods, 13(7), 581-583.

Magoč, T., & Salzberg, S. L. (2011). FLASH: fast length adjustment of short reads to improve genome assemblies. *Bioinformatics*, *27*(21), 2957-2963.

Rognes, T., Flouri, T., Nichols, B., Quince, C., Mahé, F. (2016). VSEARCH: a versatile open source tool for metagenomics. *PeerJ* *4*:e2584.
