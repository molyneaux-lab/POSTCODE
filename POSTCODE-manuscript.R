# setwd()
path <- getwd()

library("phyloseq")
library("dplyr")
library("tidyverse")
library("vegan")
library("reshape2")
library("qiime2R")
library("decontam")
#BiocManager::install("microbiome")
library("microbiome")
library(pacman)
#if(!require(pacman))install.packages("pacman")
pacman::p_load('dplyr', 'tidyr', 'gapminder',
               'ggplot2',  'ggalt',
               'forcats', 'R.utils', 'png', 
               'grid', 'ggpubr', 'scales',
               'bbplot')
library(Gmisc)
library("FSA")
library(ggbiplot)
##### Load qiime2 artefacts into R ####
# Sequencing data quality control 
# Load the metadata file into the environment
metadata <- read_csv("Original-files/metadata.csv")
metadata <- metadata %>%
  filter(!is.na(ddPCR)) %>% # drop any samples with no ddPCR data, this would've resulted in no microbiota reads anyway
  filter(!SampleType=="OR") # only select BAL samples
metadata <- column_to_rownames(metadata, var = "sample-id")

list.files(path)

# Load the qiime2 artefacts into the environment
physeq<-qza_to_phyloseq(
  features="Original-files/filtered-metadata-table.qza",
  tree="Original-files/merged-midrooted-tree.qza",
  taxonomy= "Original-files/taxonomy-merged-final.qza")
physeq 

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 945 taxa and 353 samples ]
#tax_table()   Taxonomy Table:    [ 945 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 945 tips and 934 internal nodes ]

mapping=metadata 
# 269 samples, 13 metadata columns
sample_data(physeq)<-mapping # metadata dataframe becomes the sample_data slot of physeq object
view(sample_data(physeq)) # check to see if it's correct, no OR samples. 
gplots::venn(list(mapping=rownames(mapping), physeq=sample_names(physeq))) # All samples but one match
setdiff(rownames(mapping), sample_names(physeq)) #mock samples and OR samples missing which is correct. 
  # [1] "POST.01.007a.BAL" this is a duplicate BAL sample, but this has never been sent for sequencing.

# Contamination check
ps<-physeq
summarize_phyloseq(ps)
# "1] Min. number of reads = 93"
# "2] Max. number of reads = 1339062"
# "3] Total number of reads = 9383077"
# "4] Average number of reads = 35011.4813432836"
# "5] Median number of reads = 31268"
# "6] Any OTU sum to 1 or less? YES"
# "7] Sparsity = 0.962042959804154"
# "8] Number of singletons = 161"
# "9] Percent of OTUs that are singletons \n 
# "10] Number of sample variables are: 13"

## Plot to look at library sizes per sample 
# Separated by Diagnosis and sample type i.e., negatives vs. true samples. 
df <- as.data.frame(sample_data(ps)) 
df$LibrarySize <- sample_sums(ps) # similar to rowSums/colSums but automated
df <- df[order(df$LibrarySize),]
df$Index <- seq(nrow(df))
ggplot(data=df, aes(x=Index, y=LibrarySize, color=`Diagnosis`, shape=`SampleType`)) + geom_point() +
  ylim(lower = 0, upper = 80000) # one sample is removed, BRU.1000. Is clearly an outlier due to the ridiculous librarysize (1339062).
  # Most negative samples have low library sizes, which is reassuring. 

df$Diagnosis <- as.factor(df$Diagnosis)

## Boxplot of library size by diagnosis 
ggplot(df, aes(x=Diagnosis, y = LibrarySize, 
               color = Diagnosis)) + 
  geom_point(aes(fill=Diagnosis, colour=Diagnosis),position = "identity") +
  geom_boxplot(aes(color = Diagnosis), alpha=0.8, outlier.colour = "grey") + 
  scale_color_manual(values=c("#FAAB18", "#1380A1", "#990000", "#588300", "#778899")) +
  ylim(lower=0, upper=80000) + # remove that one outlier. 
  facet_grid(.~SampleType,
             scales = "free") + theme(legend.position="right") + 
  labs(x = "", 
       y = "Library size of samples (bp)",
       title = "Library size by Diagnosis",
       subtitle = "Unfiltered data") + 
  theme_classic()+
  theme(axis.ticks.x = element_blank(),
        axis.text.x = element_text(angle = 30, hjust=1, vjust=1, size=8))
ggsave("Figures/Library sizes by Diagnosis groups.pdf", width=11, height=8, unit="in")

# Rarefaction curve; show how much % of taxa is actually covered. 
  # Made with qiime2 on HPC, or on R
genfac = factor(tax_table(ps)[, "Genus"]) # extract the genera from the ps object
# Tabulate the counts for each genera in each sample
gentab = apply(otu_table(ps), MARGIN = 2, function(x) {
  tapply(x, INDEX = genfac, FUN = sum, na.rm = TRUE, simplify = TRUE)
})
head(gentab)[, 1:10]
gentab <- data.frame(gentab)

S <- specnumber(gentab)
(raremax <- min(colSums(gentab)))
Srare <- rarefy(gentab, raremax)
plot(S, Srare, xlab = "Observed No. of genera", ylab = "Rarefied No. of genera")
abline(0, 1)
rarecurve(t(gentab), step = 20, sample = 18000, col = "darkblue", cex = 0.3, xlim= c(0,100000), label=TRUE,
          ylab = "Genera")
# manually export if needed

#### Normalisation WITHOUT removing contaminants ####
sample_data(ps)$is.neg <- sample_data(ps)$Diagnosis=="Negative Control" #Only considered reagent control as negative
contamdf.freq <- isContaminant(ps, method="combined", conc="ddPCR", neg="is.neg", threshold=0.1) 
# Gets data frame with a list of potential contaminants

table(contamdf.freq$contaminant)#FALSE=925, TRUE=20
head(which(contamdf.freq$contaminant)) # not the highest abundant taxa
# [1] 395 445 491 493 504 552

# To make a presence/absence plot
ps.pa <- transform_sample_counts(ps, function(abund) 1*(abund>0))
# Transforms data to absence or presence 
ps.pa.neg <- prune_samples(sample_data(ps.pa)$Diagnosis =="Negative Control", ps.pa)
ps.pa.pos <- prune_samples(!sample_data(ps.pa)$Diagnosis == "Negative Control", ps.pa)
# Make data.frame of prevalence in positive and negative samples
df.pa <- data.frame(pa.pos=taxa_sums(ps.pa.pos), pa.neg=taxa_sums(ps.pa.neg),
                    contaminant=contamdf.freq$contaminant)
ggplot(data=df.pa, aes(x=pa.neg, y=pa.pos, color=contaminant)) + geom_point() +
  xlab("Prevalence (Negative Controls)") + ylab("Prevalence (True Samples)") 
# Up to four identified contaminants may not be true contaminants. 

tax <- as(tax_table(ps), "matrix")
contaminants<-tax[which(contamdf.freq$contaminant),]
#write.table(contaminants,file="Original-files/Contaminants-list.txt", col.names=NA, row.names=T,sep="\t")

ps.contam <- prune_taxa(contamdf.freq$contaminant, ps) # only subset contaminants from the dataset

# Identify ASVs above 1000 reads in the list of contaminants. These are big influencers. 
filter <- phyloseq::genefilter_sample(ps.contam, filterfun_sample(function(x) x >= 1000))
ps.contam.1k <- prune_taxa(filter, ps.contam)
otu_table<-as.data.frame(ps.contam@otu_table)
tax_table<-as.data.frame(ps.contam@tax_table) 

library(tibble)
library(dplyr)

## Back to the main UNFILTERED dataset i.e., ps. 
sum(taxa_sums(ps) == 0) # how many taxa aren't present in ANY samples
summarize_phyloseq(ps)
ps

#phyloseq-class experiment-level object
#otu_table()   OTU Table:         [ 945 taxa and 268 samples ]
#sample_data() Sample Data:       [ 268 samples by 13 sample variables ]
#tax_table()   Taxonomy Table:    [ 945 taxa by 7 taxonomic ranks ]
#phy_tree()    Phylogenetic Tree: [ 945 tips and 934 internal nodes ]

# Extract taxa information to get the number of unique Phyla
tax <- as(tax_table(ps), "matrix")
tax_df <- as.data.frame(tax)
filterPhyla = unique(tax_df$Phylum)
filterPhyla <- na.omit(filterPhyla)

ps1 = subset_taxa(ps, !(!Phylum %in% filterPhyla))
ps1 # Only keep the Phyla in filterPhyla in the filtered reads dataset
summarize_phyloseq(ps1)

# Check if there are unique Phyla names
unique(as.data.frame(as(tax_table(ps1), "matrix"))$Phylum)

#Filter at 0.005% 
minTotRelAbun = 0.00005
x = taxa_sums(ps1)
keepTaxa = which((x / sum(x)) > minTotRelAbun)
prunedSet = prune_taxa(names(keepTaxa), ps1)
prunedSet # Only keeping phyla that have a minimum relative abundance of 0.005%. 

thrownTaxa = which((x / sum(x)) < minTotRelAbun)
trashSet = prune_taxa(names(thrownTaxa), ps1)
view(tax_table(trashSet)) # just cross check that those that are removed are indeed in low abundance/not normal

sum(taxa_sums(prunedSet)==0) # should be 0

summarize_phyloseq(prunedSet)
# 1] Min. number of reads = 79"
# "2] Max. number of reads = 1337785"
# "3] Total number of reads = 9325479"
# "4] Average number of reads = 34796.5634328358"
# "5] Median number of reads = 31159"
# "6] Any OTU sum to 1 or less? NO"
# "7] Sparsity = 0.844968120562237"
# "8] Number of singletons = 0"
# "9] Percent of OTUs that are singletons \n
# "10] Number of sample variables are: 14

# Normalising samples i.e., relative abundances
normalizeSample = function(x) {
  x/sum(x)}

Controls_relative = transformSampleCounts(prunedSet, normalizeSample)
otu_table(Controls_relative)
OTU1 = as(otu_table(Controls_relative), "matrix")
OTUdf = as.data.frame(OTU1)
write.csv(OTUdf, "Original-files/OTUdf.csv")

TAXdf = as(tax_table(Controls_relative), "matrix")
TAXdf = as.data.frame(TAXdf)
write.csv(TAXdf, "Original-files/tax_table.csv")

Controls_Phylum <- aggregate_taxa(Controls_relative, 'Phylum') #7 phyla most likely a contaminant
otu_table<-as.data.frame(Controls_Phylum@otu_table)
write.table(otu_table,file="Unfiltered/Phylum/Phylum-relative-abundance.txt", col.names=NA, row.names=T,sep="\t")

Controls_Class <- aggregate_taxa(Controls_relative, 'Class')
otu_table<-as.data.frame(Controls_Class@otu_table)
write.table(otu_table,file="Unfiltered/Class/Class-relative-abundance.txt", col.names=NA, row.names=T,sep="\t")

Controls_Family <- aggregate_taxa(Controls_relative, 'Family')
otu_table<-as.data.frame(Controls_Family@otu_table)
write.table(otu_table,file="Unfiltered/Family/Family-relative-abundance.txt", col.names=NA, row.names=T,sep="\t")

Controls_Genus <- aggregate_taxa(Controls_relative, 'Genus')
otu_table<-as.data.frame(Controls_Genus@otu_table)
write.table(otu_table,file="Unfiltered/Genus/Genus-relative-abundance.txt", col.names=NA, row.names=T,sep="\t")

## Have joined clinical metadata to the taxanomic reads. 
# every row is a sample and every column is a unique taxa member e.g., genus

#### Plots for manuscript ####
##### Boxplot of bacterial burden by diagnosis ####
bacterial_burden <- metadata %>%
  select(ddPCR, Diagnosis, SampleType)

bacterial_burden_stats <- bacterial_burden %>%
  group_by(Diagnosis) %>%
  mutate(median_burden = median(ddPCR),
         mean_burden = mean(ddPCR),
         q1 = quantile(ddPCR, 0.25),  # 1st quartile
         q3 = quantile(ddPCR, 0.75), # 3rd quartile
         norm_test = shapiro.test(ddPCR)$p.value) # Get p-value from shapiro.test

kruskal.test(ddPCR ~ Diagnosis, data = bacterial_burden)
dunn_Test <- dunnTest(ddPCR ~ Diagnosis, data=bacterial_burden, method = "holm")

#Comparison          Z      P.unadj        P.adj
#1                 CHP - COVID -0.3562195 7.216762e-01 1.000000e+00
#2               CHP - Healthy -0.6743500 5.000888e-01 1.000000e+00
#4                   CHP - IPF -3.1490803 1.637852e-03 9.827112e-03 # sig
#3             COVID - Healthy -0.1851490 8.531122e-01 8.531122e-01
#5                 COVID - IPF -2.0550020 3.987881e-02 1.595152e-01
#6               Healthy - IPF -2.1347475 3.278164e-02 1.639082e-01
#7      CHP - Negative Control  8.4067875 4.213936e-17 3.792543e-16 # ignore Negative
#8    COVID - Negative Control  6.6133977 3.755976e-11 2.629183e-10 # ignore Negative
#9  Healthy - Negative Control  7.7604859 8.460461e-15 6.768369e-14 # ignore Negative
#10     IPF - Negative Control  9.7464031 1.911201e-22 1.911201e-21 # ignore Negative

p <- ggplot(bacterial_burden_stats, aes(x=Diagnosis, y=ddPCR, fill=Diagnosis)) +
  geom_dotplot(binaxis="y", stackdir = "center", binwidth = 0.12, position="dodge") +
  geom_errorbar(aes(x=Diagnosis, ymin=q1, ymax=q3), width=0.3, color='black', linewidth=1) +
  scale_y_log10() + theme_classic() + 
  labs(title="Bacterial burden across disease groups",
       y="Bacteridal burden (16S rRNA gene/mL of BAL)",
       x="") +
  annotate("text", x = 3, y=18000000, label= "p<0.001", size = 4) +
  annotate("segment", x=1, xend=5, y=1e07, yend=1e07)+
  scale_fill_manual(values=c("#FAAB18", "#1380A1", "#990000", "#588300", "#778899")) 

p

ggsave("Figures/Bacterial burden by Diagnosis groups.pdf", width=11, height=8, unit="in")

##### Negative Control images #####
#### Negative PCOA Plot ####
df <- read_csv("Unfiltered/Genus/Genus-normalised-metadata.csv")
df <- df %>%
  filter(!is.na(Actinomyces)) %>%
  select(!Unknown) %>% # drop unknown genera
  filter(!SampleType=="OR") # only BAL and negative control samples (BAL flush and reagent)
abund_table <- df %>%
  select(1,15:ncol(df))
abund_table <- column_to_rownames(abund_table, var="sample-id")
metadata <- df %>%
  select(`sample-id`:`FVC %`)

## PcoA plot of BAL samples
Negative_PCA <- prcomp(abund_table, scale.=TRUE)
pc1_variance <- summary(Negative_PCA)$importance[2, 1] * 100  # Proportion of Variance for PC1
pc2_variance <- summary(Negative_PCA)$importance[2, 2] * 100  # Proportion of Variance for PC2

str(Negative_PCA$x) #x is the coordinates of each sample on the plot
Negative_PCA2 <- cbind(df, Negative_PCA$x[,1:2])

# is ddPCR a significant factor in data spread
adonis2(abund_table ~ ddPCR, data=metadata, permutations=999, method="bray")

ggplot(Negative_PCA2, aes(PC1, PC2, colour=Diagnosis))+
  geom_point() + stat_ellipse()+
  labs(title="PCoA of BAL and reagent samples by diagnosis",
       x = paste("PC1 (", pc1_variance, "% explained)", sep=""),
       y = paste("PC2 (", pc2_variance, "% explained)", sep="")) + 
  scale_color_manual(values=c("#FAAB18", "#1380A1", "#990000", "#588300", "#778899")) +
  theme_classic()+
  theme(legend.position = c(0.87,0.26),
        #legend.text = element_text(size=5),
        #legend.title = element_text(size=5.5),
        legend.key.size = unit(0.2, "cm"),
        legend.spacing = unit(0.1,"cm")) +
  annotate("text", x=5, y=4, label="Permanova: p<0.05", size=3, fontface="italic")

ggsave("Figures/PCoA of BAL samples by Diagnosis groups.pdf", width=11, height=8, unit="in")

#### Negative Stacked bar plot #####
## Presence of contaminants across TRUE samples based on negative controls
df <- read_csv("Unfiltered/Genus/Genus-normalised-metadata.csv")
df <- df %>%
  filter(!is.na(Actinomyces)) %>%
  select(!Unknown) %>% # drop unknown genera
  filter(!SampleType=="OR") # only BAL samples

metadata <- df %>%
  select(`sample-id`:`FVC %`)

negative <- df %>%
  filter(Diagnosis == "Negative Control")
negative_genera <- negative %>%
  select(Actinomyces:ncol(negative))

negative_genera <- negative_genera[,order(colSums(negative_genera),decreasing=TRUE)]
negative_genera <- negative_genera/rowSums(negative_genera)*100 # change to percentages
rowSums(negative_genera)
# only select genera where total abundance exceeds 0.01%
negative_genera <- names(negative_genera)[colSums(negative_genera) > (sum(negative_genera)*0.01)]

df_negative <- df[, negative_genera]
df_negative <- df_negative/rowSums(df_negative)*100
df_negative <- cbind(metadata, df_negative)
df_negative <- df_negative %>%
  filter(rowSums(df_negative[,15:ncol(df_negative)])>0)
df_dropped <- df_negative %>%
  filter(!rowSums(df_negative[,15:ncol(df_negative)])>0)

df_negative_long <- melt(df_negative, id.vars = c("sample-id",
                                                  "ddPCR",
                                                  "PatientID",
                                                  "SampleType",
                                                  "Diagnosis",
                                                  "Severe FVC",
                                                  "Severe TLCO",
                                                  "Severe Either PFTs",
                                                  "Severe CT",
                                                  "Severe anything",
                                                  "PFT Improvement",
                                                  "CT improvement",
                                                  "Fibrosis",
                                                  "FVC %"),
                                  variable.name = "Genus")

df_negative_long_summarised <- df_negative_long %>%
  group_by(Diagnosis, SampleType, Genus) %>%
  summarise(mean_value=mean(value),
            sd=sd(value))

df_negative_long_summarised$Genus <- as.factor(df_negative_long_summarised$Genus)
genera_contaminant <- df_negative_long_summarised %>% # only keep taxa at a 5% more abundance, more likely to be true contaminant
  filter(mean_value > 5 & Diagnosis == "Negative Control")
genera_contaminant <- as.character(genera_contaminant$Genus)

df_negative_long_summarised <- df_negative_long_summarised %>%
  filter(Genus %in% genera_contaminant)

unique(genera_contaminant)
Palette <- c(Sphingomonas = "#483d8b",
             Halomonas = "#306030",
             Bradyrhizobium = "#8e4585",
             Streptococcus = "#cb410b",
             Staphylococcus = "#f4c430")

ggplot(df_negative_long_summarised,
       aes(x = mean_value,
           y = Genus,
           fill = Genus)) +
  #Set the width of the bars in the plot
  geom_bar(stat = "identity",
           width = 0.7) +
  geom_errorbar(aes(xmax=mean_value+sd, xmin=mean_value-sd),
                width=.2,
                position=position_dodge(.9)) +
  facet_grid(SampleType ~ Diagnosis,
             scales = "free", 
             drop = TRUE)+ 
  scale_x_continuous(expand = c(0, 0),
                     limits = c(-20, 100))+
  scale_fill_manual(values = Palette) +
  labs(title="Relative abundance comparison ranked by most abundant taxa in negative control samples.",
       subtitle = "Negative control specimens were BAL Flush controls and reagent controls.",
       x = "Mean relative abundance (%)",
       y = "Genus") +
  theme(#Set the title font size
    plot.title = element_text(size=10),
    plot.subtitle = element_text(size = 8),
    legend.position = "right",
    legend.title = element_text(size=8),
    legend.text = element_text(size=8,
                               face = "italic"),
    legend.background = element_blank(),
    legend.key = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle=90,
                               hjust=1,
                               vjust=0.5,
                               size=8),
    axis.title.x = element_text(size=8),
    axis.title = element_text(size=8),
    axis.text.y = element_text(size=8),
    axis.line = element_line(size = 0.5,
                             linetype = "solid",
                             colour = "black"),
    aspect.ratio = 1.7
  )

#Save as a pdf for size to go into Inkspace figure
ggsave("Figures/Negative genera across BAL samples.pdf", width = 300, height = 150, units = c("mm"), dpi = 300)

#### Alpha-plot ASV ####
# Is the lung microbiota of COVID different than IPF, CHP and HC?
physeq_filtered <- subset_samples(physeq, !Diagnosis=="Negative Control") 
physeq_filtered <- subset_samples(physeq_filtered, !PatientID=="POST.02.005") # remove outlier sample, see "top ten taxa" image
p = plot_richness(physeq_filtered, x="Diagnosis", color="Diagnosis", measures=c("Observed","Shannon", "Chao1"))

alpha_diversity <- p$data
alpha_diversity_shannon <- alpha_diversity %>% filter(variable=="Shannon")
kruskal.test(value ~ Diagnosis, alpha_diversity_shannon)

alpha_diversity_chao <- alpha_diversity %>% filter(variable=="Chao1")
kruskal.test(value ~ Diagnosis, alpha_diversity_chao)
dunnTest(value ~ Diagnosis, alpha_diversity_chao)
# Only to COVID is it significant p = 0.02

print(p)

p + geom_boxplot(data = p$data, aes(x = Diagnosis, y = value, color = NULL), alpha = 0.1) +
  scale_colour_manual(values=c("#FAAB18", "#1380A1", "#990000", "#588300", "#778899")) + 
  theme_classic()+ theme(axis.text.x=element_blank())

ggsave("Figures/Alpha Diversity by Diagnosis.pdf", width = 11, height = 8, units = "in")

p = plot_richness(physeq_filtered, x="Diagnosis", color="Diagnosis", measures=c("Chao1"))
p + geom_boxplot(data = p$data, aes(x = Diagnosis, y = value, color = NULL), alpha = 0.1) +
  scale_colour_manual(values=c("#FAAB18", "#1380A1", "#990000", "#588300", "#778899")) + 
  theme_classic()+ theme(axis.text.x=element_blank()) + 
  annotate("text", x=1.5, y=85, label="p=0.032", size=3) + 
  annotate("segment", x=1,xend=2, y=80, yend=80) + 
  annotate("text", x=3, y=90, label="p=0.027", size=3) + 
  annotate("segment", x=2,xend=4, y=85, yend=85)
ggsave("Figures/Chao1 Diversity (ASV) by Diagnosis.pdf", width = 11, height = 8, units = "in")

#### Beta-plot ASV ####
ordu = ordinate(physeq_filtered, "PCoA", "unifrac", weighted=TRUE)
#Temporary fix to colour bug is to add a "dummy variable" in sample_data:
#sample_data(physeq)[ , 2] <- sample_data(physeq)[ ,1]
p = plot_ordination(physeq, ordu, color="Diagnosis") + geom_point(size = 3) +
  ggtitle("PCoA on weighted-UniFrac distance") + (scale_colour_brewer(type="qual", palette="Set1"))
print(p)
myplotdiagnosis <- p
myplotdiagnosis + theme_bw() + stat_ellipse() + theme_classic() + 
  theme(panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(), 
        axis.line = element_line(colour = "black")) + 
  scale_colour_manual(values=c("#FAAB18", "#1380A1", "#990000", "#588300", "#778899"))
ggsave("Figures/Beta Diversity by Diagnosis.pdf", width = 11, height = 8, units = "in")

##### Top ten taxa: Stacked bar plot, COVID ####
df <- read_csv("Unfiltered/Genus/Genus-normalised-metadata.csv")
df <- df %>% filter(Diagnosis=="COVID" | Diagnosis=="Healthy") %>% filter(SampleType=="BAL") 
df_filtered <- df %>% drop_na(Streptococcus) # 8 missing metadata entries, 1 missing microbiota data.

abund_table <- df_filtered %>%
  select(1,15:ncol(df_filtered))
abund_table <- column_to_rownames(abund_table, var="sample-id")
rowSums(abund_table)
abund_table <- abund_table/rowSums(abund_table)*100 # change to percentages

meta_table <- df_filtered %>%
  select(`sample-id`, PatientID, Diagnosis, `Severe FVC`:Fibrosis)

# Make a dataframe "top" that contains the top 10 most abundant genera. You can alter these as you want.
top <- abund_table[,order(colSums(abund_table),decreasing=TRUE)]
N <- 10
taxa_list <- colnames(top)[1:N]
N <- length(taxa_list)
top <- data.frame(top[,colnames(top) %in% taxa_list])
top

# Make a dataframe "other" that contains all of the other genera. You can alter these as you want.
other <- abund_table[,order(colSums(abund_table),decreasing=TRUE)]
# Extract list of top N Taxa
N <- dim(df_filtered)[2]
taxa_list2 <- colnames(other)[11:N]
N <- length(taxa_list2)
Others <- data.frame(other[,colnames(other) %in% taxa_list2])
Others
# Sum all the other columns into one.
Others <- rowSums(Others)
Others <- as.data.frame(Others)
Others

# Combine the top 10 genus and the other column and check the rows all still add up to 100.
top_other <- cbind(top, Others)
top_other
top_other_rowsums <- rowSums(top_other)
top_other_rowsums <- as.data.frame(top_other_rowsums)
top_other_rowsums

sample_data <- cbind(meta_table, top_other)
sample_data <- sample_data %>%
  arrange(Streptococcus, #desc() # Using desc() If you want to arrange in descending order.
  )

# OPTIONAL: Fix the order of sample IDs in the order of genus proportion.
sample_data$`sample-id` <- factor(sample_data$`sample-id`,
                                  levels=unique(sample_data$`sample-id`))

sample_data_long <- melt(sample_data, id.vars = c("sample-id",
                                                  "PatientID",
                                                  "Diagnosis",
                                                  "Severe FVC",
                                                  "Severe TLCO",
                                                  "Severe Either PFTs",
                                                  "Severe CT",
                                                  "Severe anything",
                                                  "PFT Improvement",
                                                  "CT improvement",
                                                  "Fibrosis"), variable.name = "Genus")
sample_data_long
taxa_list

Palette <- c(Streptococcus = "#990000",
             Prevotella = "#e32636",
             Veillonella = "#f94d00",
             Actinomyces = "#e25822",
             Sphingomonas = "#f4c430",
             Haemophilus = "#588300", 
             Rothia = "#306030",
             Neisseria = "#0a7e8c",
             Gemella = "#2a52be",
             Granulicatella = "#c54b8c",
             Others = "grey")

ggplot(sample_data_long,
       aes(x = `sample-id`,
           y = value,
           fill = Genus)) +
  geom_bar(stat = "identity",
           width = 0.7) +
  scale_fill_manual(values = Palette) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 100.1)) +
  facet_grid(.~Diagnosis,
             scales = "free_x",
             drop=TRUE)+
  labs(title="Ten most abundant genera in COVID patients",
       x = "Sample ID",
       y = "Relative abundance (%)")+
  theme(#Set the title font size
    plot.title = element_text(size=10),
    plot.subtitle = element_text(size = 8),
    legend.position = "right",
    legend.title = element_text(size=8),
    legend.text = element_text(size=8,
                               face = "italic"),
    legend.background = element_blank(),
    legend.key = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle=90,
                               hjust=1,
                               vjust=0.5,
                               size=8),
    axis.title.x = element_text(size=8),
    axis.title = element_text(size=8),
    axis.text.y = element_text(size=8),
    axis.line = element_line(size = 0.5,
                             linetype = "solid",
                             colour = "black"),
    aspect.ratio = 0.5
  )

ggsave("Figures/Ten most abundant Genus unfiltered_contaminant.pdf", width = 300, height = 150, units = c("mm"), dpi = 300)

# Clearly one COVID sample is an outlier. POST.02.005.BAL will be excluded.
# BRU.1032, BRU0.3853 and BRU.03815 to be excluded too for Healthy.


#### Top ten genera: outlier removed ####
df_filtered <- df_filtered %>%
  filter(Diagnosis=="COVID") %>%
  filter(!`sample-id`=="POST.02.005.BAL")
  
abund_table <- df_filtered %>%
  select(1,15:ncol(df_filtered))
abund_table <- column_to_rownames(abund_table, var="sample-id")
rowSums(abund_table)
abund_table <- abund_table/rowSums(abund_table)*100 # change to percentages

meta_table <- df_filtered %>%
  select(`sample-id`, PatientID, Diagnosis, `Severe FVC`:Fibrosis)

# Make a dataframe "top" that contains the top 10 most abundant genera. You can alter these as you want.
top <- abund_table[,order(colSums(abund_table),decreasing=TRUE)]
N <- 10
taxa_list <- colnames(top)[1:N]
N <- length(taxa_list)
top <- data.frame(top[,colnames(top) %in% taxa_list])
top

# Make a dataframe "other" that contains all of the other genera. You can alter these as you want.
other <- abund_table[,order(colSums(abund_table),decreasing=TRUE)]
# Extract list of top N Taxa
N <- dim(df_filtered)[2]
taxa_list2 <- colnames(other)[11:N]
N <- length(taxa_list2)
Others <- data.frame(other[,colnames(other) %in% taxa_list2])
Others
# Sum all the other columns into one.
Others <- rowSums(Others)
Others <- as.data.frame(Others)
Others

# Combine the top 10 genus and the other column and check the rows all still add up to 100.
top_other <- cbind(top, Others)
top_other
top_other_rowsums <- rowSums(top_other)
top_other_rowsums <- as.data.frame(top_other_rowsums)
top_other_rowsums

sample_data <- cbind(meta_table, top_other)
sample_data <- sample_data %>%
  arrange(Streptococcus, desc(Corynebacterium) # Using desc() If you want to arrange in descending order.
  )

# OPTIONAL: Fix the order of sample IDs in the order of genus proportion.
sample_data$`sample-id` <- factor(sample_data$`sample-id`,
                                  levels=unique(sample_data$`sample-id`))

sample_data_long <- melt(sample_data, id.vars = c("sample-id",
                                                  "PatientID",
                                                  "Diagnosis",
                                                  "Severe FVC",
                                                  "Severe TLCO",
                                                  "Severe Either PFTs",
                                                  "Severe CT",
                                                  "Severe anything",
                                                  "PFT Improvement",
                                                  "CT improvement",
                                                  "Fibrosis"), variable.name = "Genus")
sample_data_long
taxa_list

Palette <- c(Streptococcus = "#990000",
             Prevotella = "#e32636",
             Veillonella = "#f94d00",
             Actinomyces = "#e25822",
             Haemophilus = "#f4c430",
             Rothia = "#588300", 
             Neisseria = "#306030",
             Gemella = "#0a7e8c",
             Corynebacterium = "#2a52be",
             Granulicatella = "#c54b8c",
             Others = "grey")

ggplot(sample_data_long,
       aes(x = `sample-id`,
           y = value,
           fill = Genus)) +
  geom_bar(stat = "identity",
           width = 0.7) +
  scale_fill_manual(values = Palette) +
  scale_y_continuous(expand = c(0, 0),
                     limits = c(0, 100.1)) +
  labs(title="Ten most abundant genera in COVID patients",
       x = "Sample ID",
       y = "Relative abundance (%)")+
  theme(#Set the title font size
    plot.title = element_text(size=10),
    plot.subtitle = element_text(size = 8),
    legend.position = "right",
    legend.title = element_text(size=8),
    legend.text = element_text(size=8,
                               face = "italic"),
    legend.background = element_blank(),
    legend.key = element_blank(),
    panel.background = element_blank(),
    panel.border = element_blank(),
    panel.grid.major = element_blank(),
    panel.grid.minor = element_blank(),
    axis.text.x = element_text(angle=90,
                               hjust=1,
                               vjust=0.5,
                               size=8),
    axis.title.x = element_text(size=8),
    axis.title = element_text(size=8),
    axis.text.y = element_text(size=8),
    axis.line = element_line(size = 0.5,
                             linetype = "solid",
                             colour = "black"),
    aspect.ratio = 0.5
  )

ggsave("Figures/Ten most abundant Genus-outlier removed.pdf", width = 300, height = 150, units = c("mm"), dpi = 300)

#### Metadata analysis ####
# Severe COVID vs. non-severe COVID (are any of the variables significant, then focus on this)
df_filtered <- df_filtered %>%
  filter(Diagnosis=="COVID") %>%
  filter(!`sample-id`=="POST.02.005.BAL") %>% drop_na()

abund_table <- df_filtered %>%
  select(1,15:ncol(df_filtered))
abund_table <- column_to_rownames(abund_table, var="sample-id")
rowSums(abund_table)
abund_table <- abund_table/rowSums(abund_table)*100 # change to percentages

meta_table <- df_filtered %>%
  select(`sample-id`, PatientID, Diagnosis, `Severe FVC`:Fibrosis)
adonis2(abund_table ~ 
          `Severe CT`,
          #`Severe FVC`,
          #`Severe TLCO`,
          #`Severe Either PFTs`,
          #`Severe anything`,
          #`PFT Improvement`,
          #`CT improvement`,
          #`Fibrosis`,
          data = meta_table, permutations = 9999, method = "bray")
# Only Severe CT is p<0.1 (p=0.07)

#### Heatmap ####
# Heatmap related to all clinical parameters that are relevant
# dataframe is top_other from the stacked bar plot figure
library(ComplexHeatmap)
df <- read_csv("Unfiltered/Genus/Genus-normalised-metadata.csv")
df <- df %>% filter(Diagnosis=="COVID") %>% filter(SampleType=="BAL") 
df_filtered <- df %>% drop_na() # 8 missing metadata entries, 1 missing microbiota data.

abund_table <- df_filtered %>%
  select(1,15:ncol(df_filtered))
abund_table <- column_to_rownames(abund_table, var="sample-id")
rowSums(abund_table)
abund_table <- abund_table/rowSums(abund_table)*100 # change to percentages

meta_table <- df_filtered %>%
  select(`sample-id`, PatientID, Diagnosis, `Severe FVC`:Fibrosis)
# Make a dataframe "top" that contains the top 10 most abundant genera. You can alter these as you want.
top <- abund_table[,order(colSums(abund_table),decreasing=TRUE)]
N <- 10
taxa_list <- colnames(top)[1:N]
N <- length(taxa_list)
top <- data.frame(top[,colnames(top) %in% taxa_list])
top

# Make a dataframe "other" that contains all of the other genera. You can alter these as you want.
other <- abund_table[,order(colSums(abund_table),decreasing=TRUE)]
# Extract list of top N Taxa
N <- dim(df_filtered)[2]
taxa_list2 <- colnames(other)[11:N]
N <- length(taxa_list2)
Others <- data.frame(other[,colnames(other) %in% taxa_list2])
Others
# Sum all the other columns into one.
Others <- rowSums(Others)
Others <- as.data.frame(Others)
Others

# Combine the top 10 genus and the other column and check the rows all still add up to 100.
top_other <- cbind(top, Others)
top_other
# Heatmap will be made using a distance matrix
data_dist_rows <- vegdist(top_other, method = "bray")
row_clustering <- hclust(data_dist_rows, "average")
data_dist_columns <- vegdist(t(top_other), method = "bray")
col_clustering <- hclust(data_dist_columns, "average")

colour_palette <- colorRampPalette(colors=c("white",
                                            "orange",
                                            "darkorange",
                                            "red",
                                            "darkred",
                                            "black"))(100)
# Add annotation to heatmap AFTER distance matrix has been calculated
# the metadata does not influence the distance matrix. 
# only interested in: FVC, TLCO, CT, Fibrosis
meta_table_COVID <- df_filtered %>%
  select(`sample-id`:`FVC %`)
meta_table_COVID <- meta_table_COVID %>%
  mutate(`Severe FVC` = ifelse(`Severe CT`=="1", "FVC<80%", "FVC>80%")) %>%
  mutate(`Severe TLCO` = ifelse(`Severe TLCO`=="1", "TLCO<50%", "TLCO>50%")) %>%
  mutate(`Severe CT` = ifelse(`Severe CT`=="1", "Disease extent>20%", "Disease extent<20%")) %>%
  mutate(Fibrosis = ifelse(Fibrosis=="1", "Fibrosis on CT", "No Fibrosis on CT"))

meta_table_COVID$`Severe FVC` <- as.factor(meta_table_COVID$`Severe FVC`)
meta_table_COVID$`Severe TLCO` <- as.factor(meta_table_COVID$`Severe TLCO`)
meta_table_COVID$`Severe CT` <- as.factor(meta_table_COVID$`Severe CT`)
meta_table_COVID$Fibrosis <- as.factor(meta_table_COVID$Fibrosis)

annot_df1 <- data.frame(FVC = meta_table_COVID$`Severe FVC`)
annot_df2 <- data.frame(TLCO = meta_table_COVID$`Severe TLCO`)
annot_df3 <- data.frame(CT = meta_table_COVID$`Severe CT`)
annot_df4 <- data.frame(Fibrosis = meta_table_COVID$Fibrosis)

col1 = list(FVC = c("FVC<80%" = "#FAAB18",
                    "FVC>80%" = "#1380A1"))
col2 = list(TLCO = c("TLCO<50%" = "#FAAB18",
                     "TLCO>50%" = "#1380A1"))
col3 = list(CT = c("Disease extent<20%" = "#FAAB18",
                   "Disease extent>20%" = "#1380A1"))
col4 = list(Fibrosis = c("Fibrosis on CT" = "#FAAB18",
                         "No Fibrosis on CT" = "#1380A1"))

sidebar_annotation1 <- rowAnnotation(df = annot_df1, # Dataframe containing treatment groups
                                     col = col1, # The list of treatment groups and their assigned colours
                                     show_annotation_name = TRUE,
                                     annotation_width = unit(c(.2), "cm"), # Set the width of the side bar
                                     annotation_legend_param = list(title = "Severe FVC", # Sidebar legend title
                                                                    title_gp = gpar(fontsize = 7), # Sidebar legend title font size
                                                                    labels_gp = gpar(fontsize = 7))) # Sidebar legend label font size
sidebar_annotation2 <- rowAnnotation(df = annot_df2, # Dataframe containing treatment groups
                                     col = col2, # The list of treatment groups and their assigned colours
                                     show_annotation_name = TRUE,
                                     annotation_width = unit(c(.2), "cm"), # Set the width of the side bar
                                     annotation_legend_param = list(title = "Severe TLCO", # Sidebar legend title
                                                                    title_gp = gpar(fontsize = 7), # Sidebar legend title font size
                                                                    labels_gp = gpar(fontsize = 7))) # Sidebar legend label font size
sidebar_annotation3 <- rowAnnotation(df = annot_df3, # Dataframe containing treatment groups
                                     col = col3, # The list of treatment groups and their assigned colours
                                     show_annotation_name = TRUE, 
                                     annotation_width = unit(c(.2), "cm"), # Set the width of the side bar
                                     annotation_legend_param = list(title = "CT Severity", # Sidebar legend title
                                                                    title_gp = gpar(fontsize = 7), # Sidebar legend title font size
                                                                    labels_gp = gpar(fontsize = 7))) # Sidebar legend label font size
sidebar_annotation4 <- rowAnnotation(df = annot_df4, # Dataframe containing treatment groups
                                     col = col4, # The list of treatment groups and their assigned colours
                                     show_annotation_name = TRUE,
                                     annotation_width = unit(c(.2), "cm"), # Set the width of the side bar
                                     annotation_legend_param = list(title = "Fibrosis", # Sidebar legend title
                                                                    title_gp = gpar(fontsize = 7), # Sidebar legend title font size
                                                                    labels_gp = gpar(fontsize = 7))) # Sidebar legend label font size

heatmap <- Heatmap(as.matrix(top_other), # The dataframe containing the heatmap data
                   name = "Proportion", # Name is used as the title of the heatmap legend if shown
                   col = colour_palette, # The predefined colour palette
                   cluster_rows = row_clustering, # Cluster the rows using the predefined clustering
                   #cluster_columns = col_clustering, # Cluster the columns using the predefined clustering
                   show_row_names = TRUE, # Show or hide the row names, TRUE to show rownames
                   row_names_gp = gpar(fontsize = 6), # Row name font size
                   column_names_gp = gpar(fontsize = 8, # Column name font size
                                          fontface = "italic"), # Column names in italics
                   column_title_gp = gpar(fontsize = 8), # Column title font size
                   row_title = "Samples", # Set row title
                   row_labels = rownames(top_other),
                   row_names_side = c("left"),
                   row_title_gp = gpar(fontsize = 8), # Set row title font size
                   column_title = "", # Set column title
                   column_title_side = "bottom", # Set column title font size
                   heatmap_legend_param = list(title = "Relative\nabundance\n(%)", # Set legend title
                                               at = c(0,20,40,60,80,100), # Set legend scale breaks
                                               labels = c("0","20","40","60","80","100"), # Set legend scale labels
                                               title_gp = gpar(fontsize = 8), # Set legend title font size
                                               labels_gp = gpar(fontsize = 8))) # Set legend label font size

p <- heatmap + sidebar_annotation1 + sidebar_annotation2 + sidebar_annotation3 + sidebar_annotation4
p

pdf(file = "Figures/Heatmap-Genus_Diagnosis.pdf", width = 11, height = 8)

# print(p) saves the figure into a file
print(p)
dev.off()

#### Taxa specific analysis ####
library("dplyr")
library("vegan")
library("ggplot2")
library("tidyverse")

data <- read_csv("Unfiltered/Genus/Genus-normalised-metadata.csv")
data <- data %>%
  filter(rowSums(data[,15:ncol(data)])>0) %>%
  filter(Diagnosis=="COVID" & SampleType=="BAL") 
  #filter(Diagnosis == "COVID" & SampleType=="BAL" |Diagnosis=="Healthy" & SampleType=="BAL")


abund_table <- data %>%
  select(15:ncol(data)) # Only select the reads
abund_table$Prevotella <- abund_table$Prevotella + abund_table$`[Prevotella]`
abund_table <- abund_table %>%
  select(-`[Prevotella]`)

rownames(abund_table) <- data$`sample-id`

meta_table <- data %>%
  select(`sample-id`:`FVC %`)

meta_table$`Severe FVC` <- as.factor(meta_table$`Severe FVC`)
meta_table$`Severe anything` <- as.factor(meta_table$`Severe anything`)
meta_table$`Severe CT` <- as.factor(meta_table$`Severe CT`)
meta_table$`Severe Either PFTs` <- as.factor(meta_table$`Severe Either PFTs`)
meta_table$`Severe TLCO` <- as.factor(meta_table$`Severe TLCO`)
meta_table$`PFT Improvement` <- as.factor(meta_table$`PFT Improvement`)
meta_table$Fibrosis <- as.factor(meta_table$Fibrosis)

# bv.step ========================================================================================================
# Create bv.step function to use later
# Code sourced from: userweb.eng.gla.ac.uk/umer.ijaz/bioinformatics/ecological.html

# This R script is an extension of vegan library's bioenv()
# function and uses the bio.env() and bio.step() of
#	http://menugget.blogspot.co.uk/2011/06/clarke-and-ainsworths-bioenv-and-bvstep.html
#	The original author suggested these functions to overcome
#	the inflexibility of the bioenv() function which uses
#	a similarity matrix based on normalized "euclidean" distance.
# The new functions are given below and implement the following algorithms:
# Clarke, K. R & Ainsworth, M. 1993. A method of linking multivariate community structure to environmental variables. Marine Ecology Progress Series, 92, 205-219.
# Clarke, K. R., Gorley, R. N., 2001. PRIMER v5: User Manual/Tutorial. PRIMER-E, Plymouth, UK.
# Clarke, K. R., Warwick, R. M., 2001. Changes in Marine Communities: An Approach to Statistical Analysis and Interpretation, 2nd edition. PRIMER-E Ltd, Plymouth, UK.
# Clarke, K. R., Warwick, R. M., 1998. Quantifying structural redundancy in ecological communities. Oecologia, 113:278-289.

# Create bv.step function.
bv.step <- function(fix.mat, var.mat,
                    fix.dist.method="gower", var.dist.method="euclidean", correlation.method="spearman",
                    scale.fix=FALSE, scale.var=TRUE,
                    max.rho=0.95,
                    min.delta.rho=0.001,
                    random.selection=TRUE,
                    prop.selected.var=0.2,
                    num.restarts=10,
                    var.always.include=NULL,
                    var.exclude=NULL,
                    output.best=10
){
  
  if(dim(fix.mat)[1] != dim(var.mat)[1]){stop("fixed and variable matrices must have the same number of rows")}
  if(sum(var.always.include %in% var.exclude) > 0){stop("var.always.include and var.exclude share a variable")}
  require(vegan)
  
  if(scale.fix){fix.mat<-scale(fix.mat)}else{fix.mat<-fix.mat}
  if(scale.var){var.mat<-scale(var.mat)}else{var.mat<-var.mat}
  
  fix.dist <- vegdist(as.matrix(fix.mat), method=fix.dist.method)
  
  #an initial removal phase
  var.dist.full <- vegdist(as.matrix(var.mat), method=var.dist.method)
  full.cor <- suppressWarnings(cor.test(fix.dist, var.dist.full, method=correlation.method))$estimate
  var.comb <- combn(1:ncol(var.mat), ncol(var.mat)-1)
  RES <- data.frame(var.excl=rep(NA,ncol(var.comb)), n.var=ncol(var.mat)-1, rho=NA)
  for(i in 1:dim(var.comb)[2]){
    var.dist <- vegdist(as.matrix(var.mat[,var.comb[,i]]), method=var.dist.method)
    temp <- suppressWarnings(cor.test(fix.dist, var.dist, method=correlation.method))
    RES$var.excl[i] <- c(1:ncol(var.mat))[-var.comb[,i]]
    RES$rho[i] <- temp$estimate
  }
  delta.rho <- RES$rho - full.cor
  exclude <- sort(unique(c(RES$var.excl[which(abs(delta.rho) < min.delta.rho)], var.exclude)))
  
  if(random.selection){
    num.restarts=num.restarts
    prop.selected.var=prop.selected.var
    prob<-rep(1,ncol(var.mat))
    if(prop.selected.var< 1){
      prob[exclude]<-0
    }
    n.selected.var <- min(sum(prob),prop.selected.var*dim(var.mat)[2])
  } else {
    num.restarts=1
    prop.selected.var=1
    prob<-rep(1,ncol(var.mat))
    n.selected.var <- min(sum(prob),prop.selected.var*dim(var.mat)[2])
  }
  
  RES_TOT <- c()
  for(i in 1:num.restarts){
    step=1
    RES <- data.frame(step=step, step.dir="F", var.incl=NA, n.var=0, rho=0)
    attr(RES$step.dir, "levels") <- c("F","B")
    best.comb <- which.max(RES$rho)
    best.rho <- RES$rho[best.comb]
    delta.rho <- Inf
    selected.var <- sort(unique(c(sample(1:dim(var.mat)[2], n.selected.var, prob=prob), var.always.include)))
    while(best.rho < max.rho & delta.rho > min.delta.rho & RES$n.var[best.comb] < length(selected.var)){
      #forward step
      step.dir="F"
      step=step+1
      var.comb <- combn(selected.var, RES$n.var[best.comb]+1, simplify=FALSE)
      if(RES$n.var[best.comb] == 0){
        var.comb.incl<-1:length(var.comb)
      } else {
        var.keep <- as.numeric(unlist(strsplit(RES$var.incl[best.comb], ",")))
        temp <- NA*1:length(var.comb)
        for(j in 1:length(temp)){
          temp[j] <- all(var.keep %in% var.comb[[j]])
        }
        var.comb.incl <- which(temp==1)
      }
      
      RES.f <- data.frame(step=rep(step, length(var.comb.incl)), step.dir=step.dir, var.incl=NA, n.var=RES$n.var[best.comb]+1, rho=NA)
      for(f in 1:length(var.comb.incl)){
        var.incl <- var.comb[[var.comb.incl[f]]]
        var.incl <- var.incl[order(var.incl)]
        var.dist <- vegdist(as.matrix(var.mat[,var.incl]), method=var.dist.method)
        temp <- suppressWarnings(cor.test(fix.dist, var.dist, method=correlation.method))
        RES.f$var.incl[f] <- paste(var.incl, collapse=",")
        RES.f$rho[f] <- temp$estimate
      }
      
      last.F <- max(which(RES$step.dir=="F"))
      RES <- rbind(RES, RES.f[which.max(RES.f$rho),])
      best.comb <- which.max(RES$rho)
      delta.rho <- RES$rho[best.comb] - best.rho
      best.rho <- RES$rho[best.comb]
      
      if(best.comb == step){
        while(best.comb == step & RES$n.var[best.comb] > 1){
          #backward step
          step.dir="B"
          step <- step+1
          var.keep <- as.numeric(unlist(strsplit(RES$var.incl[best.comb], ",")))
          var.comb <- combn(var.keep, RES$n.var[best.comb]-1, simplify=FALSE)
          RES.b <- data.frame(step=rep(step, length(var.comb)), step.dir=step.dir, var.incl=NA, n.var=RES$n.var[best.comb]-1, rho=NA)
          for(b in 1:length(var.comb)){
            var.incl <- var.comb[[b]]
            var.incl <- var.incl[order(var.incl)]
            var.dist <- vegdist(as.matrix(var.mat[,var.incl]), method=var.dist.method)
            temp <- suppressWarnings(cor.test(fix.dist, var.dist, method=correlation.method))
            RES.b$var.incl[b] <- paste(var.incl, collapse=",")
            RES.b$rho[b] <- temp$estimate
          }
          RES <- rbind(RES, RES.b[which.max(RES.b$rho),])
          best.comb <- which.max(RES$rho)
          best.rho<- RES$rho[best.comb]
        }
      } else {
        break()
      }
      
    }
    
    RES_TOT <- rbind(RES_TOT, RES[2:dim(RES)[1],])
    print(paste(round((i/num.restarts)*100,3), "% finished"))
  }
  
  RES_TOT <- unique(RES_TOT[,3:5])
  
  
  if(dim(RES_TOT)[1] > output.best){
    order.by.best <- RES_TOT[order(RES_TOT$rho, decreasing=TRUE)[1:output.best],]
  } else {
    order.by.best <-  RES_TOT[order(RES_TOT$rho, decreasing=TRUE), ]
  }
  rownames(order.by.best)<-NULL
  
  order.by.i.comb <- c()
  for(i in 1:length(selected.var)){
    f1 <- which(RES_TOT$n.var==i)
    f2 <- which.max(RES_TOT$rho[f1])
    order.by.i.comb <- rbind(order.by.i.comb, RES_TOT[f1[f2],])
  }
  rownames(order.by.i.comb)<-NULL
  
  if(length(exclude)<1){var.exclude=NULL} else {var.exclude=exclude}
  out <- list(
    order.by.best=order.by.best,
    order.by.i.comb=order.by.i.comb,
    best.model.vars=paste(colnames(var.mat)[as.numeric(unlist(strsplit(order.by.best$var.incl[1], ",")))], collapse=","),
    best.model.rho=order.by.best$rho[1],
    var.always.include=var.always.include,
    var.exclude=var.exclude
  )
  out
  
}
#res.bv.step.biobio=====================================================================================

# Define parameter commands (other options are available for each one).
# Bray (Bray-Curtis) is the most commonly used one and in usually best to use.
# You can change these without changing the command in the functions.
cmethod<-"pearson" # Correlation method to use: pearson, spearman, kendall
fmethod<-"bray" # Fixed distance method: euclidean, manhattan, gower, altGower, canberra, bray, kulczynski, morisita,horn, binomial, and cao
vmethod<-"bray" # Variable distance method: euclidean, manhattan, gower, altGower, canberra, bray, kulczynski, morisita,horn, binomial, and cao
nmethod<-"bray" # NMDS distance method:  euclidean, manhattan, gower, altGower, canberra, bray, kulczynski, morisita,horn, binomial, and cao

# Create the command to find the top 10 genus driving the separation of samples.
res.bv.step.biobio <- bv.step(wisconsin(abund_table),
                              wisconsin(abund_table),
                              fix.dist.method=fmethod,
                              var.dist.method=vmethod,
                              correlation.method=cmethod,
                              scale.fix=FALSE,
                              scale.var=FALSE,
                              max.rho=0.95,
                              min.delta.rho=0.001,
                              random.selection=TRUE,
                              prop.selected.var=0.3,
                              num.restarts=10,
                              output.best=10,
                              var.always.include=NULL)

res.bv.step.biobio
# WARNINGS: These refer to empty rows and can be ignored.
# Sub-samplING the community matrix in bv.step to only include a subset of genus
# means there will be cases when the "abund_table" will have empty rows,
# especially when we are selecting only 10 genus.


# Find the 10 best subset of genus driving the separation of samples
taxaNames<-colnames(abund_table)
bestTaxaFit <-""
for(i in (1:length(res.bv.step.biobio$order.by.best$var.incl)))
{
  bestTaxaFit[i] <- paste(paste(taxaNames[as.numeric(unlist(strsplit(res.bv.step.biobio$order.by.best$var.incl[i],
                                                                     split = ",")))],
                                collapse = ' + '),
                          " = ",
                          res.bv.step.biobio$order.by.best$rho[i],
                          sep = "")
}

# Select the column (genus) names of the top ten genus
bestTaxaFit <- data.frame(bestTaxaFit)
colnames(bestTaxaFit) <- "Best combination of taxa with similarity score"
bestTaxaFit
# Generate the NMDS plot
MDS_res = metaMDS(abund_table,
                  distance = nmethod,
                  k = 3, # Set the number of dimensions
                  smin = 0.0001, # 0.0001 is the default
                  sfgrmin = 0.0000001, # 0.0000001 is the default
                  sratmax = 0.99999, # 0.99999 is the default
                  try = 200, # Minimum number of random starts
                  maxit = 300, # Maximum number of random starts
                  trymax = 300, # Number of attempts at calculating the NMDS plot to find the best one
                  plot = TRUE, # TRUE shows the plots developing as the code is running.
                  engine = "monoMDS")
MDS_res


# Test the significance of the taxa in the NMDS plot
bio.keep <- as.numeric(unlist(strsplit(res.bv.step.biobio$order.by.best$var.incl[1], ",")))
bio.fit <- envfit(MDS_res, abund_table[,bio.keep,drop=F], perm = 999)
bio.fit

# Get the vectors for bio.fit
# Limit the taxa shown to only those with significant p-values.
# This can be useful to avoid overcrowding the NMDS plot with overlapping names.
# Extract the data for the arrows.
vectors <- as.list(bio.fit$vectors)
vectors
# r2 is the degree to which the taxa explain differences in the NMDS plot.
# Scale the length of the arrows by the r2.
arrows <- as.data.frame(vectors$arrows*sqrt(vectors$r))
arrows
# Extract the p-values for the arrows.
pvals <- as.data.frame(vectors$pvals)
pvals
# Create a dataframe with the arrows and p-values
significant_arrows <- cbind(arrows, pvals)
significant_arrows
# Subset the dataframe to only those with p<0.05
significant_taxa <- subset(significant_arrows, pvals < 0.05)
significant_taxa


# Test the significance of the metadata variables in the NMDS plot
en = envfit(MDS_res, meta_table, permutations = 9999, na.rm = TRUE)
en

# Test the nMDS plot and factors using base R.
plot(MDS_res)
plot(en)

# Extract the categorical factor metadata variables that explain NMDS plot variance.
en_factor <- as.list(en$factors)
en_factor
# Create a dataframe with the centroids (the center point of each factor level) and p-values.
en_factor_pvals <- as.data.frame(en_factor$pvals)
en_factor_pvals <- rownames_to_column(en_factor_pvals, "group")
en_factor_var_id <- as.data.frame(en_factor$var.id)
en_factor_centroids <- as.data.frame(en_factor$centroids)
en_factor_centroids <- rownames_to_column(en_factor_centroids, "centroid")
bind <- cbind(en_factor_var_id, en_factor_centroids)
names(bind)[names(bind) == 'en_factor$var.id'] <- 'group'
merged <- merge(bind, en_factor_pvals, by = "group")
names(merged)[names(merged) == 'en_factor$pvals'] <- 'pvals'
merged

# Limit metadata factors shown to only those with significant p-values
# This is to avoid overcrowding the NMDS plot with overlapping names
# Subset the dataframe to only those with p < 0.05
en_coord_cat_sig <- subset(merged, merged$pvals < 0.05) 
en_coord_cat_sig <- en_coord_cat_sig[,2:4]
en_coord_cat_sig <- rownames_to_column(en_coord_cat_sig, "none")
en_coord_cat_sig <- column_to_rownames(en_coord_cat_sig, "centroid")
en_coord_cat_sig <- en_coord_cat_sig[,-1]
# Use en_coord_cat_sig to annotate the ggplot.
en_coord_cat_sig


# Prepare the dataframe to make the ggplot.
# Extract the coordinates to draw the NMDS plot.
dataframe <- scores(MDS_res, display=c("sites"))
dataframe

# Convert into a dataframe.
dataframe <- data.frame(dataframe)

# Add Treatment group information to the dataframe.
dataframe <- cbind(dataframe, meta_table)
# Set the levels of the metadata factor that is used to color the plot.
# This will set the order groups will appear in the legend and when assigning colors and symbols.
dataframe$`Severe CT` <- ifelse(dataframe$`Severe CT`=="1", "Disease extent>20%", "Disease extent<20%")
dataframe$`Severe CT` <- factor(dataframe$`Severe CT`,
                                levels = c("Disease extent>20%",
                                           "Disease extent<20%"))
head(dataframe)
dataframe <- dataframe %>%
  select(-PatientID)

#Create the ggplot
p <- ggplot(data = dataframe,
            aes(x = NMDS1,
                y = NMDS2)) +
  # Set the metadata to be used to change the colour and shape of the plot symbols.
  geom_point(aes(colour = `Severe CT`,
                 #shape = `Severe CT`
  )) +
  # Set the color of the symbols used for each group.
  scale_colour_manual(name = "Severe CT",
                      labels = c("Disease extent>20%",
                                 "Disease extent<20%"),
                      values = c("red1",
                                 "steelblue3")) +
  # Set the shape of the symbols used for each group.
  scale_shape_manual(name = "CT severity",
                     labels = c("Disease extent>20%",
                                "Disease extent<20%"),
                     # Define the shapes of the symbols used in the plot (each number is a different R symbol).
                     values = c(22, 22,22,22)) +
  # Defining the x axis limits helps stop taxa names being cut off.
  scale_x_continuous(limits = c(-1.4, 1.1)) +
  # Set the plot title.
  labs(title = "Taxa analysis (Unfiltered 16S data)")+
  #subtitle = "At Genus level, p=0.0514") +
  theme(# Define the plot title.
    plot.title = element_text(size=10),
    plot.subtitle = element_text(size=8),
    # Define the plot margin size.
    plot.margin = margin(t = 4, r = 4, b = 4, l = 4),
    # Define the legend position.
    legend.position = "right",
    # Define the x axis legend spacing.
    legend.spacing.x = unit(.1, 'cm'),
    # Define the y axis legend spacing.
    legend.spacing.y = unit(.1, 'cm'),
    # Format the legend title.
    legend.title = element_text(size=8),
    # Format the size of the legend text.
    legend.text = element_text(size=8),
    # Remove the grey background.
    legend.background = element_blank(),
    # Remove rectangle around the legend
    legend.box.background = element_blank(),
    # Remove the box from the legend.
    legend.key = element_blank(),
    # Remove the grey plot background.
    panel.background = element_blank(),
    # Remove the plot border.
    panel.border = element_blank(),
    # Remove the major plot grid lines.
    panel.grid.major = element_blank(),
    # Remove the minor plot grid lines.
    panel.grid.minor = element_blank(),
    # Format x axis title.
    axis.title.x = element_text(size=8, colour = "black"),
    # Format x axis labels.
    axis.text.x = element_text(size=8, colour = "black"),
    # Format the y axis title text size.
    axis.title.y = element_text(size=8, colour = "black"),
    # Format the axis label text size.
    axis.text.y = element_text(size=8, colour = "black"),
    # Format x and y axis ticks.
    axis.ticks = element_line(size = 0.35, colour = "black"),
    # Format x and y axis lines.
    axis.line = element_line(size = 0.35, colour = "black"),
    # Define the plot aspect ratio.
    aspect.ratio = 1) +
  
  # Add the taxa arrows to the ggplot.
  geom_segment(data = significant_taxa,
               aes(x = 0,
                   y = 0,
                   xend = NMDS1,
                   yend = NMDS2),
               # Format the size of the arrow head.
               arrow = arrow(length = unit(0.2, "cm")),
               # Format the color of the arrows.
               color = "black",
               # Format the thickness of the arrows.
               alpha = .8) +
  
  # Add the taxa names to the ggplot.
  geom_text(data = as.data.frame(significant_taxa * 1), # "significant_taxa*1.5" moves the names 1.5 x the length of the arrow.
            aes(NMDS1,
                NMDS2,
                label = rownames(significant_taxa)),
            # Format genus font.
            fontface = "italic",
            # Format genus names.
            color = "black",
            # Format the size of the genus names.
            size = 3) +

# Add the categorical metadata variable points to the ggplot.
geom_point(data = en_coord_cat_sig,
           aes(x = NMDS1, y = NMDS2),
           # Format the shape of the symbol.
           shape = "diamond",
           # Format the color of the symbol.
           colour = "darkorchid4",
           # Format the size of the symbol.
           size = 5) +
  
  # Add the categorical metadata variable names to the ggplot.
  geom_text(data = en_coord_cat_sig,
            aes(x = NMDS1, y = NMDS2 + 0.2), # Add +0.08 to move the names above the diamonds.
            label = row.names(en_coord_cat_sig),
            # Format the color of the variable names.
            colour = "darkorchid4",
            # Format the font of the variable names.
            fontface = "bold",
            # Format the size of the variable names.
            size = 3,)

# Print the plot within R.
p

# This will print into the current working directory file.
# Set the saved pdf file name and define the size of the plot (height and width are in inches).
pdf("Figures/16S Unfiltered-Genus analysis.pdf", width = 12, height = 12)

# print(p) saves the figure as a pdf file.
print(p)
dev.off()

#### Bar chart COVID vs. HC ####
## At genus level
df_genus <- read_csv("Unfiltered/Genus/Genus-normalised-metadata.csv")
df_genus <- df_genus %>%
  filter(Diagnosis == "COVID" | Diagnosis=="Healthy"|Diagnosis=="IPF") %>%
  #filter(!Diagnosis=="COVID") %>% filter(!Diagnosis=="Negative Control") %>%
  filter(!SampleType=="OR") %>%
  select(-Unknown) %>%
  drop_na(ddPCR | Actinomyces) %>%
  filter(!PatientID=="POST.02.005") %>% filter(!PatientID == "BRU.1032") %>% 
  filter(!PatientID == "BRU.03853") %>% filter(!PatientID == "BRU.03853") %>%
  filter(!PatientID=="BRU.03815") 

abund_table <- df_genus %>%
  select(1,15:ncol(df_genus))
abund_table <- column_to_rownames(abund_table, var="sample-id")
rowSums(abund_table)
abund_table <- abund_table/rowSums(abund_table)*100 # change to percentages

meta_table <- df_genus %>%
  select(1,5)
# top ten family taxa
top <- abund_table[,order(colSums(abund_table),decreasing=TRUE)]
N <- 10
taxa_list <- colnames(top)[1:N]
N <- length(taxa_list)
top <- data.frame(top[,colnames(top) %in% taxa_list])
top
df_genus_top <- cbind(meta_table, top)
taxa_list

kruskal.test(Streptococcus ~ Diagnosis, df_genus_top) #p=0.033
dunnTest(Streptococcus ~ Diagnosis, df_genus_top, method="holm") # COVID-IPF

kruskal.test(Prevotella ~ Diagnosis, df_genus_top)
kruskal.test(Veillonella ~ Diagnosis, df_genus_top)
kruskal.test(Actinomyces ~ Diagnosis, df_genus_top)
kruskal.test(Neisseria ~ Diagnosis, df_genus_top)
kruskal.test(Haemophilus ~ Diagnosis, df_genus_top)
kruskal.test(Sphingomonas ~ Diagnosis, df_genus_top) #p=0.007
dunnTest(Sphingomonas ~ Diagnosis, df_genus_top, method="holm") # COVID-Healthy, COVID-IPF
kruskal.test(Rothia ~ Diagnosis, df_genus_top)
kruskal.test(Gemella ~ Diagnosis, df_genus_top)
kruskal.test(Granulicatella ~ Diagnosis, df_genus_top)

df_long <- melt(df_genus_top, id.vars = "Diagnosis", 
                measure.vars = c("Streptococcus", "Prevotella", 
                                 "Veillonella", "Actinomyces", 
                                 "Neisseria", "Haemophilus", 
                                 "Sphingomonas", "Rothia", 
                                 "Gemella", "Granulicatella"), 
                variable.name = "Genus",
                factorsAsStrings = TRUE, na.rm = TRUE)
df_long_summarised <- df_long %>%
  group_by(Diagnosis, Genus) %>%
  summarise(mean_value=mean(value),
            sd=sd(value),
            q1 = quantile(value, 0.25),  # 1st quartile
            q3 = quantile(value, 0.75)) # 3rd quartile

dat_text <- data.frame(
  label=c(" ", "*","*"),
  Diagnosis=c("Healthy","COVID","IPF"),
  x = c(1,1,1), 
  y = c(40,50,40)
)

dat_text2 <- data.frame(
  label=c("**", "**",""),
  Diagnosis=c("Healthy","COVID","IPF"),
  x = c(7,7,7), 
  y = c(20,20,20)
)

p <- ggplot(df_long_summarised, aes(x=Genus, y=mean_value)) + 
  geom_bar(aes(y = mean_value, x = Genus, fill = Genus),
           stat="identity") +
  geom_errorbar(aes(x=Genus, ymin=(mean_value-sd), ymax=(mean_value+sd)), width=0.3, color='black', linewidth=0.5)
df_long_summarised$Diagnosis <- ordered(df_long_summarised$Diagnosis, levels = c("Healthy", "COVID","IPF", "CHP"))
levels(df_long_summarised$Diagnosis)
p + scale_fill_manual(values = c("darkseagreen", "forestgreen", "darkgreen", "gold", "sandybrown", "coral", "lightskyblue", "royalblue", "darkblue", "firebrick")) + 
  theme_classic() + 
  geom_text(data=dat_text, mapping= aes(x=x, y=y, label=label), size = 6) +
  geom_text(data=dat_text2, mapping= aes(x=x, y=y, label=label), size = 6) +
  labs(x = "Genus", y = "Relative abundance (%)") + 
  theme(axis.text.x = element_text(angle=90)) +
  guides(fill = guide_legend(title = "Genus")) + 
  facet_grid(Diagnosis~.) 

ggsave("Figures/Ordered barplot at Genus by Diagnosis.pdf", width=11, height = 8, unit="cm")

#### Stacked bar plot by Diagnosis ####
## At phylum level
df_phylum <- read_csv("Unfiltered/phylum/phylum-normalised-metadata.csv")
df_phylum <- df_phylum %>%
  filter(Diagnosis=="COVID"|Diagnosis=="Healthy"|Diagnosis=="IPF") %>%
  #filter(!Diagnosis=="COVID") %>% filter(!Diagnosis=="Negative Control") %>%
  filter(!SampleType=="OR") %>%
  drop_na(ddPCR | Bacteroidetes) %>%
  filter(!PatientID=="POST.02.005") %>% filter(!PatientID == "BRU.1032") %>% 
  filter(!PatientID == "BRU.03853") %>% filter(!PatientID == "BRU.03853") %>%
  filter(!PatientID=="BRU.03815")

meta_table <- df_phylum %>%
  select(1,5)
abund_table <- df_phylum %>%
  select(1,15:ncol(df_phylum))
abund_table <- column_to_rownames(abund_table, var="sample-id")
rowSums(abund_table)
abund_table <- abund_table/rowSums(abund_table)*100 # change to percentages

# top ten family taxa
top <- abund_table[,order(colSums(abund_table),decreasing=TRUE)]
N <- 10
taxa_list <- colnames(top)[1:N]
N <- length(taxa_list)
top <- data.frame(top[,colnames(top) %in% taxa_list])
top
taxa_list

df_phylum_top <- cbind(meta_table, top)
df_long <- melt(df_phylum_top, id.vars = "Diagnosis", 
                measure.vars = c("Firmicutes", "Bacteroidetes", 
                                 "Proteobacteria", "Actinobacteria", 
                                 "Fusobacteria", "Verrucomicrobia"), 
                variable.name = "Phylum",
                factorsAsStrings = TRUE, na.rm = TRUE)
df_long_summarised <- df_long %>%
  group_by(Diagnosis, Phylum) %>%
  summarise(mean_value=mean(value),
            sd=sd(value))
p <- ggplot(df_long_summarised, aes(x=Diagnosis, y=mean_value)) + 
  geom_bar(aes(y = mean_value, x = Diagnosis, fill = Phylum),
           stat="identity", position=position_stack())
df_long_summarised$Diagnosis <- ordered(df_long_summarised$Diagnosis, levels = c("Healthy", "COVID", "IPF"))
levels(df_long_summarised$Diagnosis)
p + scale_fill_manual(values = c("darkblue", "royalblue", "lightskyblue", "gold", "sandybrown", "coral", "darkgreen", "forestgreen", "palegreen", "firebrick")) + 
  theme_classic() + labs(x = "Diagnosis", y = "Relative abundance (%)") + 
  theme(axis.text.x = element_text(angle=0)) +
  guides(fill = guide_legend(title = "Phylum")) 
ggsave("Figures/Stacked barplot at Phylum by Diagnosis.pdf", width=11, height = 8, unit="cm")

