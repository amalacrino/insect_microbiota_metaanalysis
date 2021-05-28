# Host phylogeny drives the structure and diversity of insect microbiota  

### Antonino Malacrinò

#### *Institute for Evolution and Biodiversity, Westfälische Wilhelms-Universität Münster, Münster, Germany*

## Abstract 

Microorganisms have an enormous impact on most of the life that inhabits our planet. Insects are an excellent example, as research showed that several microbial species are essential for insect nutrition, reproduction, fitness, defence and many other functions. More recently, we assisted to an exponential growth of studies describing the taxonomical composition of bacterial communities across insects’ phylogeny. However, there is still an outstanding question that needs to be answered: which factors contribute most in shaping insects’ microbiomes? This study tries to find an answer to this question by taking advantage of publicly available sequencing data and reanalysing over 2,500 samples of insect-associated bacterial communities under a common framework. Results suggest that insect taxonomy has a wider impact on the structure and diversity of their associated microbial communities than the other factors considered (diet, sex, life stage, sample origin and treatment). Also, a survey of the literature highlights several methodological limitations that needs to be considered in future research endeavours. This study proofs the amount of collective effort that lead to the current understanding of insect-microbiota interactions and their influence on insect biology, ecology and evolution with potential impact on insect conservation and management practices.


## 1. Data

SRA accession number are provided in the file SRA_Acc_List.txt

## 2. Merge PE reads

PE reads were merged using FLASH 1.2.11 (Magoč and Salzberg, 2011)

```bash
for i in *_1.fastq; do FILENAME=`basename $i _1.fastq`; echo $FILENAME; flash ${FILENAME}_1.fastq ${FILENAME}_2.fastq --output-prefix=${FILENAME}.merged; done

find . ! -name '*.merged.extendedFrags.fastq' -delete
for file in *; do mv ${file} ${file/.merged.extendedFrags/}; done
```

## 3. Process data with VSEARCH

Data handling was carried out using VSEARCH (Rognes et al. 2016), assigning taxonomy using SILVA v132 database (Quast et al. 2012). 

Run the pipeline according to https://github.com/torognes/vsearch/wiki/VSEARCH-pipeline and assign taxonomy using the command

```bash
vsearch --sintax $PRJ_DIR/otus.fasta --db $PRJ_DIR/silva_16s_v123.fa --tabbedout $PRJ_DIR/ASV_tax_raw.txt --sintax_cutoff 0.5 --threads $NTHREADS
```

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

Infer phylosymbiosis

```R
ref_tree <- ape::read.tree('Insecta_species.nwk')
list.tree <- as.vector(ref_tree$tip.label)
list.gm.w <- sample_data(GM.whole)$Insect_ID
common.list.w <- intersect(list.gm.w, list.tree)
merged.GM.w <- subset_samples(GM.whole, Insect_ID %in% common.list.w)
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

filt.tree.w <- generate.filt.tree(merged.GM.w)
filt.tree.g <- generate.filt.tree(merged.GM.g)

tip.age <- function(x){
  data.tree <- setNames(x$edge.length[sapply(1:length(x$tip.label),function(x,y) which (y==x),y=x$edge[,2])],x$tip.label) 
  data.tree <- as.data.frame(data.tree)
  data.tree <- data.table::setDT(data.tree, keep.rownames = TRUE)[]
}

tip.age.w <- tip.age(filt.tree.w)
tip.age.g <- tip.age(filt.tree.g)
dv.w <- merge(div.w, tip.age.w, by.x = "Insect_ID", by.y = "rn")
dv.g <- merge(div.g, tip.age.g, by.x = "Insect_ID", by.y = "rn")
tree.w <- ape::cophenetic.phylo(filt.tree.w) 
tree.g <- ape::cophenetic.phylo(filt.tree.g)

cor.test( ~ Shannon + data.tree, data = dv.w, method = "pearson", continuity = FALSE)
cor.test( ~ Shannon + data.tree, data = dv.g, method = "pearson", continuity = FALSE)

generate.dist.mat <- function(x){
  mergedGM.tree <- merge_samples(x, "Insect_ID")
  dist.mat <- phyloseq::distance(mergedGM.tree, method = "bray")
  dist.mat <- as.matrix(dist.mat)
  return(dist.mat)
}

dist.mat.w.a <- generate.dist.mat(merged.GM.w)
mantel.w <- mantel(tree.w, dist.mat.w.a, method = "pearson", permutations = 9999, na.rm = TRUE)
dist.mat.g.a <- generate.dist.mat(merged.GM.g)
mantel.g <- mantel(tree.g, dist.mat.g.a, method = "pearson", permutations = 9999, na.rm = TRUE)
                                             
plot.corr.shannon <- function(x, y, z){
  plot <- ggplot(x, aes(x=data.tree, y=Shannon)) + 
          theme_bw(base_size = 15) +
          geom_point(size = 3, alpha = 0.5) +
          ggtitle(y) +
          theme(legend.position="none",
                legend.text = element_text(size = 12),
                panel.background = element_rect(size = 0.5, linetype = "solid"),
                panel.grid.major = element_blank(), 
                panel.grid.minor = element_blank(),
                plot.title = element_text(size = 14, hjust = 0.5),
                axis.title=element_text(size=12)) +
            xlab("Divergence time (Myr)") +
            ylab(z) +
            geom_smooth(method=lm, se = FALSE) +
            xlim(0, 400) +
            ylim(0, 5)
  return(plot)
}

plot.corr.mantel <- function(x, y, z){
  aa <- as.vector(x)
  tt <- as.vector(y)
  mat <- data.frame(aa,tt)
  plot <- ggplot(mat, aes(x=tt, y=aa)) + 
                    theme_bw(base_size = 15) +
                    geom_point(size = 3, alpha = 0.5) +
                    theme(legend.position="none",
                          legend.text = element_text(size = 12),
                          panel.background = element_rect(size = 0.5, linetype = "solid"),
                          panel.grid.major = element_blank(), 
                          panel.grid.minor = element_blank(),
                          plot.title = element_text(size = 14, hjust = 0.5),
                          axis.title=element_text(size=12)) +
                      xlab("Host phylogenetic distance") +
                      ylab(z) +
                      geom_smooth(method=lm, se = FALSE) +
                      xlim(0, 850) +
                      ylim(0, 1.05) +
                      scale_y_continuous(labels = scales::number_format(accuracy = 0.1))
  return(plot)
}

plot.sh.w <- plot.corr.shannon(dv.w, "Whole insects", "Shannon index")
plot.sh.g <- plot.corr.shannon(dv.g, "Insect guts", " ")
plot.mantel.w <- plot.corr.mantel(dist.mat.w.a, tree.w, "Bray-Curtis distance")
plot.mantel.g <- plot.corr.mantel(dist.mat.g.a, tree.g, " ")

p <- ggpubr::ggarrange(plot.sh.w, plot.sh.g,
                       plot.mantel.w, plot.mantel.g,
               ncol = 2, nrow = 2,  align = "hv", 
               widths = 1, heights = 1,
               labels = c("A)", "C)", "B)", "D)"),
               common.legend = F)
p                                             
```

## References

Magoč, T., & Salzberg, S. L. (2011). FLASH: fast length adjustment of short reads to improve genome assemblies. *Bioinformatics*, *27*(21), 2957-2963.

Rognes, T., Flouri, T., Nichols, B., Quince, C., Mahé, F. (2016). VSEARCH: a versatile open source tool for metagenomics. *PeerJ* *4*:e2584.
