#load packages
library(Biostrings)
library(ggplot2)
library(Peptides)
library(motifRG)
library(rtracklayer)
library(stringr)
library(peplib)
library(extrafont)
library(ape)
library(reshape2)

font_import()

loadfonts(device="win")

#load palettes
cbPalette <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")
cbscale <- c('#8e0152','#c51b7d','#de77ae','#f1b6da','#fde0ef','#e6f5d0','#b8e186','#7fbc41','#4d9221','#276419')

#define ggplot2 theme
theme_science <- function (base_size = 12, base_family = "Arial Rounded MT Bold") 
{
  theme_bw(base_size = base_size, base_family = base_family) %+replace% 
    theme(panel.border = element_blank(), axis.line = element_line(colour = "black", size=2), 
          panel.grid.major = element_line(), panel.grid.major.x = element_blank(),
          axis.line.x = element_line(colour= "black", size=1),  axis.line.y = element_line(colour= "black", size=1),
          panel.grid.major.y = element_blank(), panel.grid.minor = element_line(), 
          panel.grid.minor.x = element_blank(), panel.grid.minor.y = element_blank(), 
          strip.background = element_rect(colour = "black", 
                                          size = 0.5), legend.key = element_blank())
}

#load plasmid sequencing results
Plasmid_L1 = read.table("/Users/Kyle/Dropbox/TI-seq/NextSeq/Plasmid_L1_01_counts")
Plasmid_L2 = read.table("/Users/Kyle/Dropbox/TI-seq/NextSeq/Plasmid_L2_01_counts")
Plasmid_L3 = read.table("/Users/Kyle/Dropbox/TI-seq/NextSeq/Plasmid_L3_01_counts")
Plasmid_L4 = read.table("/Users/Kyle/Dropbox/TI-seq/NextSeq/Plasmid_L4_01_counts")

#bind results from all four lanes
Plasmid = rbind(Plasmid_L1, Plasmid_L2, Plasmid_L3, Plasmid_L4)
#Group by sequence
Plasmid <- dplyr::group_by(Plasmid, V2)
#get total counts per sequence
Plasmid <- dplyr::summarise(Plasmid, Counts=sum(V1))

#summary statistics for plasmid counts
sum(Plasmid$Counts)
median(Plasmid$Counts)
mean(Plasmid$Counts)

#plot histogram of plasmid counts
hist(log10(Plasmid$Counts))

#load counts for each bin in experiment 1
Bin1_e1 <- read.table("/Users/Kyle/Dropbox/TI-seq/NextSeq/Exp_1_01_counts")
Bin2_e1 <- read.table("/Users/Kyle/Dropbox/TI-seq/NextSeq/Exp_1_02_counts")
Bin3_e1 <- read.table("/Users/Kyle/Dropbox/TI-seq/NextSeq/Exp_1_03_counts")
Bin4_e1 <- read.table("/Users/Kyle/Dropbox/TI-seq/NextSeq/Exp_1_04_counts")
Bin5_e1 <- read.table("/Users/Kyle/Dropbox/TI-seq/NextSeq/Exp_1_05_counts")

#rename columns in data.frames
colnames(Plasmid) <- c("Sequence", "Reads_Plasmid")
colnames(Bin1_e1) <- c("Reads_B1", "Sequence")
colnames(Bin2_e1) <- c("Reads_B2", "Sequence")
colnames(Bin3_e1) <- c("Reads_B3", "Sequence")
colnames(Bin4_e1) <- c("Reads_B4", "Sequence")
colnames(Bin5_e1) <- c("Reads_B5", "Sequence")

#merge each dataset sequentially by sequence
Plasmid_bin1 <- merge(Plasmid, Bin1_e1, by="Sequence", all = TRUE)
Plasmid_bin1_2  <- merge(Plasmid_bin1, Bin2_e1, by="Sequence", all = TRUE)
Plasmid_bin1_2_3  <- merge(Plasmid_bin1_2, Bin3_e1, by="Sequence", all = TRUE)
Plasmid_bin1_2_3_4 <- merge(Plasmid_bin1_2_3, Bin4_e1, by="Sequence", all = TRUE)
Plasmid_bin1_2_3_4_5 <- merge(Plasmid_bin1_2_3_4, Bin5_e1, by="Sequence", all = TRUE)

#rename merged data.frame
All_reads <- Plasmid_bin1_2_3_4_5

#make new column containing sequence as a DNAStringSet
All_reads$DNA <- DNAStringSet(All_reads$Sequence)

#make NAs (no reads in one sample but reads in at least one other) == 0
All_reads$Reads_Plasmid[is.na(All_reads$Reads_Plasmid)] <- 0
All_reads$Reads_B1[is.na(All_reads$Reads_B1)] <- 0
All_reads$Reads_B2[is.na(All_reads$Reads_B2)] <- 0
All_reads$Reads_B3[is.na(All_reads$Reads_B3)] <- 0
All_reads$Reads_B4[is.na(All_reads$Reads_B4)] <- 0
All_reads$Reads_B5[is.na(All_reads$Reads_B5)] <- 0

#get sum of all bins
All_reads$total_bins <- All_reads$Reads_B1 + All_reads$Reads_B2 + All_reads$Reads_B3 + All_reads$Reads_B4 + All_reads$Reads_B5

#write counts for exp1
write.table(All_reads, "all_reads_e1.txt")

#sequences with less than 10 total counts
All_reads <- subset(All_reads, !total_bins < 10)

#get summary statistics
sum(All_reads$total_bins)
mean(All_reads$total_bins)

#plot histogram of total counts
ggplot(All_reads, aes(log10(total_bins))) + 
  xlab("Total Counts per Sequence (log10)") + 
  ylab("Frequency") + geom_histogram(fill="white", colour="black") + 
  scale_y_continuous(limits=c(0,20000), expand = c(0,0)) +
  annotate("text", x=3.5, y=15000, label=paste("Mean =", round(mean(All_reads$total_bins),2))) +
  theme_science()
ggsave("hist_total_counts.tiff", height =4, width=4, units="in")

#get ratios of each bin to total
All_reads$B1_ratio <- All_reads$Reads_B1/All_reads$total_bins
All_reads$B2_ratio <- All_reads$Reads_B2/All_reads$total_bins
All_reads$B3_ratio <- All_reads$Reads_B3/All_reads$total_bins
All_reads$B4_ratio <- All_reads$Reads_B4/All_reads$total_bins
All_reads$B5_ratio <- All_reads$Reads_B5/All_reads$total_bins

#adjust ratio to score or RFU based on FACS results
All_reads$Score <- (All_reads$B1_ratio*1) + (All_reads$B2_ratio*2) + (All_reads$B3_ratio*3) + (All_reads$B4_ratio*4) + (All_reads$B5_ratio*5)
All_reads$RFU <- (All_reads$B1_ratio*20) + (All_reads$B2_ratio*120) + (All_reads$B3_ratio*600) + (All_reads$B4_ratio*3600) + (All_reads$B5_ratio*12000)

#load tAI of sequence library
tai <- read.table("/Users/Kyle/Dropbox/TI-seq/NextSeq/cds_lib_tai.txt")

#make new tAI data.frame with different labels and merge with All_reads
tai <- data.frame(Sequence = tai$sequence, tAI = tai$sequences.tai, AA = tai$AA)
All_reads <- merge(All_reads, tai, by="Sequence")


#Make amino acid sequence a character string
All_reads$AA_character <- as.character(All_reads$AA)

#identify sequences with stop codons
All_reads$Stop <- grepl("......TAA|......TGA|......TAG|...TAA...|...TAG...|...TGA...|TAA......|TGA......|TAG......", All_reads$Sequence)
All_reads$Amber <- grepl("......TAG|...TAG...|TAG......", All_reads$Sequence)
All_reads$Ochre <- grepl("......TAA|...TAA...|TAA......", All_reads$Sequence)
All_reads$Opal <- grepl("......TGA|...TGA...|TGA......", All_reads$Sequence)

#define function for finding stop codons
stop <- function(x) { 
  if(!grepl("......TAA|......TGA|......TAG|...TAA...|...TAG...|...TGA...|TAA......|TGA......|TAG......", x)) y <- "Absent"
  if(grepl("......TAG|...TAG...|TAG......", x)) y <- "Amber (UAG)"
  if(grepl("......TAA|...TAA...|TAA......", x)) y <- "Ochre (UAA)"
  if(grepl("......TGA|...TGA...|TGA......", x)) y <- "Opal (UGA"
  return(y)
}

#apply stop codon function
All_reads$Stop_codon <- sapply(All_reads$Sequence, stop)

#find sequences with rare codons
All_reads$Rare <- grepl("......AGG|...AGG...|AGG......|......AGA|...AGA...|AGA......", All_reads$Sequence)
All_reads$Rare_2 <- grepl("AGGAGA...|AGG...AGG|AGG...AGA|AGGAGG...|AGAAGG...|AGAAGA...|AGA...AGG|AGA...AGA|...AGGAGG|...AGGAGA|...AGAAGG|...AGAAGA", All_reads$Sequence)
All_reads$Rare_3 <- grepl("AGGAGGAGG|AGGAGGAGA|AGGAGAAGG|AGGAGAAGA|AGAAGGAGG|AGAAGGAGA|AGAAGAAGG|AGAAGAAGA", All_reads$Sequence)

rare_codons <- read.table("Rare_codon_combinations.txt")

All_reads$RLI <- grepl(paste(rare_codons$V1, collapse = "|"), All_reads$Sequence)

#find sequences with NcoI and XhoI cutsites
All_reads$xho <- grepl("CTCGAG...|.CTCGAG..|..CTCGAG.|...CTCGAG", All_reads$Sequence)
All_reads$nco <- grepl("CCATGG...|.CCATGG..|..CCATGG.|...CCATGG", All_reads$Sequence)

#remove sequences with NcoI or XhoI cutsites, likely arrise from sequencing errors as those cuts sites are used during library production
All_reads <- subset(All_reads, xho == "FALSE")
All_reads <- subset(All_reads, nco == "FALSE")

#find hydrophobicity and charge of peptides
All_reads$hydrophobicity <- hydrophobicity(All_reads$AA, scale = "KyteDoolittle")
All_reads$charge <- pI(All_reads$AA, pKscale = "EMBOSS")

write.table(All_reads, "all_reads_e1_annotated.txt")

#make violin plot of GFP Score for sequences with and without a stop codon
ggplot(All_reads, aes(Stop, Score)) + 
  geom_violin(fill="grey") + 
  theme_science() + xlab("Contains Stop Codon") + 
  ylab("GFP Score")
ggsave("Stop_exp1.tiff", height=4, width=4, units="in", dpi=300)

#make violin plot of GFP RFU for sequences with and without a stop codon
ggplot(All_reads, aes(Stop, RFU)) + 
  geom_violin(fill="grey") + 
  theme_science() + xlab("Contains Stop Codon") + 
  ylab("GFP RFU")
ggsave("Stop_exp1_RFU.tiff", height=4, width=4, units="in", dpi=300)

#make violin plot of GFP Score for sequences by stop codon type
ggplot(All_reads, aes(Stop_codon, Score)) + 
  geom_violin(fill="grey") + 
  theme_science() + xlab("Stop Codon Type") + 
  ylab("GFP Score")
ggsave("Stop_codon_exp1.tiff", height=4, width=4, units="in", dpi=300)

#make violin plot of GFP RFU for sequences by stop codon type
ggplot(All_reads, aes(Stop_codon, RFU)) + 
  geom_violin(fill="grey") + 
  theme_science() + xlab("Stop Codon Type") + 
  ylab("GFP RFU")
ggsave("Stop_codon_exp1_RFU.tiff", height=4, width=4, units="in", dpi=300)


#define function for finding amber vs CAG codons
amber_gln <- function(x) { 
  if(!grepl("......TAG|...TAG...|TAG......|......CAG|...CAG...|CAG......", x)) y <- "Absent"
  if(grepl("TAG......", x)) y <- "Amber Codon 3"
  if(grepl("...TAG...", x)) y <- "Amber Codon 4"
  if(grepl("......TAG", x)) y <- "Amber Codon 5"
  if(grepl("CAG......", x)) y <- "CAG Codon 3"
  if(grepl("...CAG...", x)) y <- "CAG Codon 4"
  if(grepl("......CAG", x)) y <- "CAG Codon 5"
  return(y)
}

#applies Amber_gln function
All_reads$Amber <- sapply(All_reads$Sequence, amber_gln)

#make violin plot of GFP Score for sequences with Amber or CAG codon
ggplot(All_reads, aes(Amber, Score)) + 
  geom_violin(fill="grey") + 
  theme_science() + xlab("") +
  ylab("GFP Score") + theme(axis.text.x=element_text(angle = 45, hjust = 1))
ggsave("amber_gln_exp1.tiff", height=4, width=4, units="in", dpi=300)

#make boxplot of GFP Score for sequences with Amber or CAG codon
ggplot(All_reads, aes(Amber, Score)) + 
  geom_boxplot(fill="grey", notch = TRUE) + 
  theme_science() + xlab("") +
  ylab("GFP Score") + theme(axis.text.x=element_text(angle = 45, hjust = 1))
ggsave("amber_gln_exp1_box.tiff", height=4, width=4, units="in", dpi=300)

#make violin plot of GFP RFU for sequences with Amber or CAG codon
ggplot(All_reads, aes(Amber, RFU)) + 
  geom_violin(fill="grey") + 
  theme_science() + xlab("") +
  ylab("GFP RFU") + theme(axis.text.x=element_text(angle = 45, hjust = 1))
ggsave("amber_gln_exp1_RFU.tiff", height=4, width=4, units="in", dpi=300)

#make jitter plot of GFP Score for sequences with Amber or CAG codon
ggplot(All_reads, aes(Amber, Score)) + 
  geom_jitter(alpha=0.2) + 
  theme_science() + xlab("") + 
  ylab("GFP Score") + theme(axis.text.x=element_text(angle = 45, hjust = 1))
ggsave("amber_gln_exp1_jitter.tiff", height=4, width=4, units="in", dpi=300)

#make jitter plot of GFP RFU for sequences with Amber or CAG codon
ggplot(All_reads, aes(Amber, RFU)) + 
  geom_jitter(alpha=0.2) + 
  theme_science() + xlab("") + 
  ylab("GFP RFU") + theme(axis.text.x=element_text(angle = 45, hjust = 1))
ggsave("amber_gln_exp1_jitter_RFU.tiff", height=4, width=4, units="in", dpi=300)

#make new data.frame for studying Amber and CAG codons
Amber_gln_test <- All_reads

#split sequence into individual nucleotides
split_amber_gln <- colsplit(Amber_gln_test$Sequence, "", c("N1", "N2", "N3", "N4", "N5", "N6", "N7", "N8", "N9"))

#bind split to data.frame
Amber_gln_test <- cbind(Amber_gln_test, split_amber_gln)

#paste N1-3, N4-6 and N7-9 to make codon sequences
Amber_gln_test$C1 <- do.call(paste, c(Amber_gln_test[c("N1", "N2", "N3")], sep = "")) 
Amber_gln_test$C2 <- do.call(paste, c(Amber_gln_test[c("N4", "N5", "N6")], sep = "")) 
Amber_gln_test$C3 <- do.call(paste, c(Amber_gln_test[c("N7", "N8", "N9")], sep = "")) 

#Copy above columns into new columns named C#_amber or C#_cag
Amber_gln_test$C1_amber <- Amber_gln_test$C1
Amber_gln_test$C2_amber <- Amber_gln_test$C2
Amber_gln_test$C3_amber <- Amber_gln_test$C3

Amber_gln_test$C1_cag <- Amber_gln_test$C1
Amber_gln_test$C2_cag <- Amber_gln_test$C2
Amber_gln_test$C3_cag <- Amber_gln_test$C3

#substitute TAG with X
Amber_gln_test$C1_amber <- gsub("TAG", "X", Amber_gln_test$C1_amber)
Amber_gln_test$C2_amber <- gsub("TAG", "X", Amber_gln_test$C2_amber)
Amber_gln_test$C3_amber <- gsub("TAG", "X", Amber_gln_test$C3_amber)

#substitute CAG with X
Amber_gln_test$C1_cag <- gsub("CAG", "X", Amber_gln_test$C1_cag)
Amber_gln_test$C2_cag <- gsub("CAG", "X", Amber_gln_test$C2_cag)
Amber_gln_test$C3_cag <- gsub("CAG", "X", Amber_gln_test$C3_cag)

#Combine codon sequences into a single string
Amber_gln_test$Amber_X <- do.call(paste, c(Amber_gln_test[c("C1_amber", "C2_amber", "C3_amber")], sep = "")) 
Amber_gln_test$CAG_X <- do.call(paste, c(Amber_gln_test[c("C1_cag", "C2_cag", "C3_cag")], sep = "")) 

#make data.frame for Amber or CAG containing sequences
Amber_test <- data.frame(Identity = Amber_gln_test$Amber_X, Score_Amber = Amber_gln_test$Score)
CAG_test <- data.frame(Identity = Amber_gln_test$CAG_X, Score_CAG = Amber_gln_test$Score)

#merge by sequence, X now refers to CAG or TAG
Amber_CAG <- merge(Amber_test, CAG_test, by="Identity")

#find all sequences with X
Amber_CAG$X <- grepl("X", Amber_CAG$Identity)

#plot scatter of Score with GFP score for sequence with Amber on y and CAG on X
ggplot(subset(Amber_CAG, X ==TRUE), aes(Score_CAG, Score_Amber)) + geom_point(alpha=0.2) +geom_smooth(method="loess") + 
  xlab("GFP Score for CAG Containing Sequence") + ylab("GFP Score for Amber Containing Sequence") + theme_science()
ggsave("amber_cag_scatter.tiff")

#remove sequences with stop codons
All_reads_nostop <- subset(All_reads, Stop == "FALSE")

#plot density and histograms of GFP score or RFU
ggplot(All_reads_nostop, aes(Score)) + 
  geom_density(fill="#56B4E9") + 
  xlab("GFP Score") + 
  ylab("Density") + theme_science() +
  scale_y_continuous(limits=c(0,0.7), expand = c(0,0))
ggsave("density_nostop_exp1.tiff", height =4, width=4, units="in")

ggplot(All_reads_nostop, aes(RFU)) + 
  geom_density(fill="#56B4E9") + 
  xlab("GFP RFU") + 
  ylab("Density") + theme_science()  +
  scale_y_continuous(limits=c(0,0.00032), expand = c(0,0))
ggsave("density_nostop_exp1_RFU.tiff", height =4, width=4, units="in")

ggplot(All_reads_nostop, aes(Score)) + geom_histogram() + theme_science()
ggsave("hist_nostop_exp1.tiff")

ggplot(All_reads_nostop, aes(RFU)) + geom_histogram() + theme_science()
ggsave("hist_nostop_exp1_RFU.tiff")

ggplot(All_reads, aes(Score)) + 
  geom_density(fill="#E69F00") + 
  theme_science() + 
  xlab("GFP Score") + 
  ylab("Density") + 
  theme_science() + 
  scale_y_continuous(limits=c(0,0.7), expand = c(0,0))
ggsave("density_exp1.tiff", height =4, width=4, units="in")

ggplot(All_reads, aes(RFU)) + 
  geom_density(fill="#E69F00") + 
  theme_science() + 
  xlab("GFP Score") + 
  ylab("Density") + 
  theme_science() + 
  scale_y_continuous(limits=c(0,0.00032), expand = c(0,0))
ggsave("density_exp1_RFU.tiff", height =4, width=4, units="in")

ggplot(All_reads, aes(Score)) + geom_histogram() + theme_science()
ggsave("hist_exp1.tiff")

ggplot(All_reads, aes(RFU)) + geom_histogram() + theme_science()
ggsave("hist_exp1_RFU.tiff")

#plot density for GFP Score or RFU for all sequence types
ggplot(All_reads, aes(Score, group=Stop_codon, colour=Stop_codon)) + 
  geom_density(size=2) + 
  theme_science() +
  xlab("GFP Score") + 
  ylab("Density") + 
  scale_colour_manual(values=cbPalette) +
  scale_y_continuous(limits=c(0,1.1), expand = c(0,0))
ggsave("hist_all_exp1.tiff")

ggplot(All_reads, aes(RFU, group=Stop_codon, colour=Stop_codon)) + 
  geom_density(size=2) + 
  theme_science() +
  xlab("GFP RFU") + 
  ylab("Density") + 
  scale_colour_manual(values=cbPalette) +
  scale_y_continuous(limits=c(0,0.0008), expand = c(0,0))
ggsave("hist_all_exp1_RFU.tiff")

ggplot() + geom_density(data=All_reads, aes(Score), fill="#E69F00", colour="#E69F00", alpha=0.8) + 
  geom_density(data=All_reads_nostop, aes(Score), fill="#56B4E9", colour="#56B4E9", alpha=0.8) + 
  xlab("GFP Score") + ylab("Density") + theme_science() + scale_y_continuous(limits=c(0,0.7), expand = c(0,0))

#write All_reads_nostop as .txt
write.table(All_reads_nostop, "all_reads_nostop.txt")

#make boxplots for sequences with and without rare codons for GFP Score or RFU
ggplot(All_reads_nostop, aes(Rare, Score)) + geom_boxplot(notch=TRUE) + xlab("Contains 1 Rare Codon") + ylab("GFP Score") + theme_science() 
ggsave("rare_exp1.tiff")
ggplot(All_reads_nostop, aes(Rare_2, Score)) + geom_boxplot(notch=TRUE)+ xlab("Contains 2 Rare Codons") + ylab("GFP Score") + theme_science()
ggsave("rare_2_exp1.tiff")
ggplot(All_reads_nostop, aes(Rare_3, Score)) + geom_boxplot(notch=TRUE)+scale_x_discrete(labels = c("Absent", "Present")) + xlab("3 Rare Arg Codons")+ ylab("GFP Score") + theme_science()
ggsave("rare_3_exp1.tiff",height = 3, width = 3, units= "in")

ggplot(All_reads_nostop, aes(Rare, RFU)) + geom_boxplot(notch=TRUE)+ xlab("Contains 1 Rare Codon") + ylab("GFP Score") + theme_science()
ggsave("rare_exp1_RFU.tiff")
ggplot(All_reads_nostop, aes(Rare_2, RFU)) + geom_boxplot(notch=TRUE)+ xlab("Contains 2 Rare Codons") + ylab("GFP Score") + theme_science()
ggsave("rare_2_exp1_RFU.tiff")
ggplot(All_reads_nostop, aes(Rare_3, RFU)) + geom_boxplot(notch=TRUE)+ scale_x_discrete(labels = c("Absent", "Present")) +xlab("3 Rare Arg Codons") + ylab("GFP Score") +theme_science()
ggsave("rare_3_exp1_RFU.tiff")

t.test(All_reads_nostop$Score[All_reads_nostop$Rare_3 == TRUE], All_reads_nostop$Score[All_reads_nostop$Rare_3 == FALSE])

ggplot(All_reads_nostop, aes(RLI, Score)) + geom_boxplot(notch=TRUE) + scale_x_discrete(labels = c("Absent", "Present")) + xlab("3 Rare RLI Codons") + ylab("GFP Score") + theme_science() 
ggsave("RLI_exp1.tiff", height = 3, width = 3, units= "in")
t.test(All_reads_nostop$Score[All_reads_nostop$RLI == TRUE], All_reads_nostop$Score[All_reads_nostop$RLI == FALSE])


# determine correlation between GFP Score and hydrophobicity of the tripeptide and plot
cor.test(All_reads_nostop$Score, All_reads_nostop$hydrophobicity)
ggplot(All_reads_nostop, aes(hydrophobicity, Score)) + xlab("Hydrophobicity (Kyte Doolittle scale") + ylab("GFP Score") + geom_point(alpha = 0.1) + theme_science()
ggsave("hydrophobicity_exp1.tiff")
# determine correlation between GFP Score and isolectric point of the tripeptide and plot
cor.test(All_reads_nostop$Score, All_reads_nostop$charge)
ggplot(All_reads_nostop, aes(charge, Score)) + xlab("pI") + ylab("GFP Score") + geom_point(alpha = 0.1) + theme_science()
ggsave("charge_exp1.tiff")
# determine correlation between GFP Score and tAI of the tripeptide and plot
cor.test(All_reads_nostop$Score, All_reads_nostop$tAI)
ggplot(All_reads_nostop, aes(tAI, Score)) + geom_point(alpha=0.1)+ xlab("tAI") + ylab("GFP Score")+ theme_science() + ylab("GFP Score") 
ggsave("tai_exp1.tiff")

#same as above but for RFU instead of GFP score
cor.test(All_reads_nostop$RFU, All_reads_nostop$hydrophobicity)
ggplot(All_reads_nostop, aes(hydrophobicity, RFU)) + geom_point(alpha=0.1) + xlab("Hydrophobicity (Kyte Doolittle scale") + ylab("GFP Score") + theme_science()
ggsave("hydrophobicity_exp1_RFU.tiff")
cor.test(All_reads_nostop$RFU, All_reads_nostop$charge)
ggplot(All_reads_nostop, aes(charge, RFU)) + geom_point(alpha=0.1) + xlab("pI") + ylab("GFP Score") + theme_science()
ggsave("charge_exp1_RFU.tiff")
cor.test(All_reads_nostop$RFU, All_reads_nostop$tAI)
ggplot(All_reads_nostop, aes(tAI, RFU)) + geom_point(alpha=0.1) + xlab("tAI") + ylab("GFP Score")+ theme_science()
ggsave("tai_exp1_RFU.tiff")

#determine GC frequency
All_reads_nostop$G <- letterFrequency(All_reads_nostop$DNA, "G")
All_reads_nostop$C <- letterFrequency(All_reads_nostop$DNA, "C")
All_reads_nostop$GC <- (All_reads_nostop$G + All_reads_nostop$C)

#determine AT frequency
All_reads_nostop$A <- letterFrequency(All_reads_nostop$DNA, "A")
All_reads_nostop$U <- letterFrequency(All_reads_nostop$DNA, "T")
All_reads_nostop$AU <- (All_reads_nostop$A + All_reads_nostop$U)

ggplot(All_reads_nostop, aes(as.factor(GC), Score)) + geom_boxplot(notch=TRUE) + theme_science()
ggsave("GC_exp1.tiff")
ggplot(All_reads_nostop, aes(Score, colour=as.factor(GC))) + geom_step(stat="ecdf") + theme_science()
ggsave("GC_ecdf_exp1.tiff")

ggplot(All_reads_nostop, aes(as.factor(GC), RFU)) + geom_boxplot(notch=TRUE) + theme_science()
ggsave("GC_exp1_RFU.tiff")
ggplot(All_reads_nostop, aes(RFU, colour=as.factor(GC))) + geom_step(stat="ecdf") + theme_science()
ggsave("GC_ecdf_exp1_RFU.tiff")

ggplot(All_reads_nostop, aes(as.factor(AU), Score)) + 
  geom_boxplot(notch=TRUE, fill="grey") + 
  theme_science() + 
  xlab("Number of A or U nt") + ylab("GFP Score")
ggsave("AU_exp1.tiff", height =3, width=5, units="in")

ggplot(All_reads_nostop, aes(as.factor(AU), RFU)) + 
  geom_boxplot(notch=TRUE, fill="grey") +  
  theme_science() + 
  xlab("Number of A or U nt") + ylab("GFP RFU")
ggsave("AU_exp1_RFU.tiff", height =3, width=5, units="in")

ggplot(All_reads_nostop, aes(Score, colour=as.factor(AU))) + 
  geom_step(stat="ecdf", size=1) + 
  scale_colour_manual(values =cbscale, name="Number of \nA or U nt") + 
  theme_science() + 
  xlab("GFP Score") + ylab("Cummulative Fraction")
ggsave("AU_ecdf_exp1.tiff", height =4, width=4, units="in")

ggplot(All_reads_nostop, aes(RFU, colour=as.factor(AU))) + 
  geom_step(stat="ecdf", size=1) + 
  scale_colour_manual(values =cbscale, name="Number of \nA or U nt") + 
  theme_science() + 
  xlab("GFP RFU") + ylab("Cummulative Fraction")
ggsave("AU_ecdf_exp1_RFU.tiff", height =4, width=4, units="in")

cor.test(All_reads_nostop$Score, All_reads_nostop$GC)
cor.test(All_reads_nostop$Score, All_reads_nostop$GC, method="spearman")

ggplot(All_reads_nostop, aes(Reads_Plasmid, Score)) + geom_point(alpha=0.5) + theme_science()
cor.test(All_reads_nostop$Score, All_reads_nostop$Reads_Plasmid)
ggsave("Plasmid_score_exp1.tiff")

ggplot(All_reads_nostop, aes(Reads_Plasmid, RFU)) + geom_point(alpha=0.5) + theme_science()
cor.test(All_reads_nostop$RFU, All_reads_nostop$Reads_Plasmid)
ggsave("Plasmid_RFU_exp1.tiff")


high_sequences <- subset(All_reads_nostop, Score > 4)
low_sequences <- subset(All_reads_nostop, Score < 3)
high <- high_sequences$DNA
#write.fasta(as.list(high_sequences$DNA), as.list(row.names(high_sequences)),  "high.fasta")
low <- low_sequences$DNA
#write.fasta(as.list(low_sequences$DNA), as.list(row.names(low_sequences)),  "low.fasta")
#write.fasta(low, "low.fasta")
category <- c(rep(1, length(high)), rep(0, length(low))) 
MD.motifs <- findMotif(append(high, low),category,  max.motif=6,enriched=TRUE, start.width = 6, max.width=6)
MD.motifs2 <- findMotif(append(high, low),category,  max.motif=6,enriched=TRUE, start.width = 4, max.width=9)
MD.motifs3 <- findMotif(append(high, low),category,  max.motif=6,enriched=TRUE, start.width = 5, max.width=9)
### Get summary of motifs 
summaryMotif(MD.motifs$motifs, MD.motifs$category)
summaryMotif(MD.motifs2$motifs, MD.motifs2$category)
summaryMotif(MD.motifs3$motifs, MD.motifs3$category)
### Create table of motifs in Latex 
motifLatexTable(MD.motifs, main="MD motifs")
motifLatexTable(MD.motifs2, main="MD motifs")
motifLatexTable(MD.motifs3, main="MD motifs")
### Create table of motifs in Html 
motifHtmlTable(MD.motifs, dir="/Users/Kyle/Dropbox/TI-seq/NextSeq/motif_1_exp1.html")
motifHtmlTable(MD.motifs2, dir="/Users/Kyle/Dropbox/TI-seq/NextSeq/motif_2_exp1.html")
motifHtmlTable(MD.motifs3, dir="/Users/Kyle/Dropbox/TI-seq/NextSeq/motif_3_exp1.html")

category <- c(rep(0, length(high)), rep(1, length(low))) 
MD.motifs_low <- findMotif(append(high, low),category,  max.motif=6,enriched=TRUE, start.width = 6, max.width=6)
MD.motifs2_low <- findMotif(append(high, low),category,  max.motif=6,enriched=TRUE, start.width = 4, max.width=9)
MD.motifs3_low <- findMotif(append(high, low),category,  max.motif=6,enriched=TRUE, start.width = 5, max.width=9)
### Get summary of motifs 
summaryMotif(MD.motifs_low$motifs, MD.motifs_low$category)
summaryMotif(MD.motifs2_low$motifs, MD.motifs2_low$category)
summaryMotif(MD.motifs3_low$motifs, MD.motifs3_low$category)


All_reads_nostop$motif1 <- grepl("AA[G|C|A]ATT...|.AA[G|C|A]ATT..|..AA[G|C|A]ATT.|...AA[G|C|A]ATT", All_reads_nostop$Sequence)

#function to assign motif positon
m1 <- function(x) { 
  if(!grepl("AA[G|C|A]ATT...|.AA[G|C|A]ATT..|..AA[G|C|A]ATT.|...AA[G|C|A]ATT", x)) y <- "Absent"
  if(grepl("AA[G|C|A]ATT...", x)) y <- "Position 1"
  if(grepl(".AA[G|C|A]ATT..", x)) y <- "Position 2"
  if(grepl("..AA[G|C|A]ATT.", x)) y <- "Position 3"
  if(grepl("...AA[G|C|A]ATT", x)) y <- "Position 4"
  return(y)
}
#applies above function to new column
All_reads_nostop$motif1_all <- sapply(All_reads_nostop$Sequence, m1)

All_reads_nostop$motif2 <- grepl("AA[T|G|A]TAT...|.AA[T|G|A]TAT..|..AA[T|G|A]TAT.|...AA[T|G|A]TAT", All_reads_nostop$Sequence)


m2 <- function(x) { 
  if(!grepl("AA[T|G|A]TAT...|.AA[T|G|A]TAT..|..AA[T|G|A]TAT.|...AA[T|G|A]TAT", x)) y <- "Absent"
  if(grepl("AA[T|G|A]TAT...", x)) y <- "Position 1"
  if(grepl(".AA[T|G|A]TAT..", x)) y <- "Position 2"
  if(grepl("..AA[T|G|A]TAT.", x)) y <- "Position 3"
  if(grepl("...AA[T|G|A]TAT", x)) y <- "Position 4"
  return(y)
}

All_reads_nostop$motif2_all <- sapply(All_reads_nostop$Sequence, m2)

All_reads_nostop$motif3 <- grepl("AAATT", All_reads_nostop$Sequence)


m3 <- function(x) { 
  if(!grepl("AAATT", x)) y <- "Absent"
  if(grepl("AAATT....", x)) y <- "Position 1"
  if(grepl(".AAATT...", x)) y <- "Position 2"
  if(grepl("..AAATT..", x)) y <- "Position 3"
  if(grepl("...AAATT.", x)) y <- "Position 4"
  if(grepl("....AAATT", x)) y <- "Position 5"
  return(y)
}

All_reads_nostop$motif3_all <- sapply(All_reads_nostop$Sequence, m3)

All_reads_nostop$motif4 <- grepl("CAGGC", All_reads_nostop$Sequence)

m4 <- function(x) { 
  if(!grepl("CAGGC", x)) y <- "Absent"
  if(grepl("CAGGC....", x)) y <- "Position 1"
  if(grepl(".CAGGC...", x)) y <- "Position 2"
  if(grepl("..CAGGC..", x)) y <- "Position 3"
  if(grepl("...CAGGC.", x)) y <- "Position 4"
  if(grepl("....CAGGC", x)) y <- "Position 5"
  return(y)
}

All_reads_nostop$motif4_all <- sapply(All_reads_nostop$Sequence, m4)


gray = c("White", "Dark Grey", "Dark Grey", "Dark Grey", "Dark Grey", "Dark Grey", "Dark Grey", "Dark Grey")

ggplot(All_reads_nostop, aes(motif1, Score, fill = motif1)) + 
  geom_boxplot(notch=TRUE) + scale_fill_manual(values=gray, guide = FALSE)+
  theme_science() + 
  xlab("Contains AAVATT") + ylab("GFP Score")
ggsave("motif1_exp1.tiff", height=4, width=4, unit="in")
ggplot(All_reads_nostop, aes(motif1_all, Score, fill = motif1_all)) + 
  geom_boxplot(notch=TRUE) + scale_fill_manual(values=gray, guide = FALSE) +
  theme_science() + 
  xlab("Position of AAVATT") + ylab("GFP Score")
ggsave("motif1_all_exp1.tiff", height=4, width=5, unit="in")

ggplot(All_reads_nostop, aes(motif1, RFU, fill = motif1)) + 
  geom_boxplot(notch=TRUE) + scale_fill_manual(values=gray, guide = FALSE) +
  theme_science() + 
  xlab("Contains AAVATT") + ylab("GFP RFU")
ggsave("motif1_exp1_RFU.tiff", height=4, width=4, unit="in")
ggplot(All_reads_nostop, aes(motif1_all, RFU, fill = motif1_all)) + 
  geom_boxplot(notch=TRUE) + scale_fill_manual(values=gray, guide = FALSE) +
  theme_science() + 
  xlab("Position of AAVATT") + ylab("GFP RFU")
ggsave("motif1_all_exp1_RFU.tiff", height=4, width=5, unit="in")

All_reads_nostop_motif1_p2 <- subset(All_reads_nostop, motif1_all == "Position 2")

m1_p2_n1 <- function(x) { 
  if(grepl("AAA[G|C|A]ATT..", x)) y <- "A"
  if(grepl("TAA[G|C|A]ATT..", x)) y <- "T"
  if(grepl("CAA[G|C|A]ATT..", x)) y <- "C"
  if(grepl("GAA[G|C|A]ATT..", x)) y <- "G"
  return(y)
}

All_reads_nostop_motif1_p2$motif1_p2_n1 <- sapply(All_reads_nostop_motif1_p2$Sequence, m1_p2_n1)

ggplot(All_reads_nostop_motif1_p2, aes(motif1_p2_n1, Score)) + 
  geom_boxplot(notch=TRUE, fill="grey") + 
  theme_science() + 
  xlab("Nucleotide at Position 1") + ylab("GFP Score")
ggsave("motif1_p2_n1_exp1.tiff", height=4, width=5, unit="in")

ggplot(All_reads_nostop_motif1_p2, aes(motif1_p2_n1, RFU)) + 
  geom_boxplot(notch=TRUE, fill="grey") + 
  theme_science() + 
  xlab("Nucleotide at Position 1") + ylab("GFP RFU")
ggsave("motif1_p2_n1_exp1_RFU.tiff", height=4, width=5, unit="in")

All_reads_nostop_motif1_p4 <- subset(All_reads_nostop, motif1_all == "Position 4")

AA_split_m1 <- colsplit(All_reads_nostop_motif1_p4$AA_character, "", c("AA_1", "AA_2", "AA_3"))

All_reads_nostop_motif1_p4 <- cbind(All_reads_nostop_motif1_p4, AA_split_m1)

ggplot(All_reads_nostop_motif1_p4, aes(reorder(AA_1, -Score), Score)) + 
  geom_boxplot(notch=TRUE, fill="grey") + 
  theme_science() + 
  xlab("AA at Position 1") + ylab("GFP Score") 
ggsave("motif1_p4_aa1_exp1.tiff", height=4, width=5, unit="in")

ggplot(All_reads_nostop_motif1_p4, aes(reorder(AA_1, -RFU), RFU)) + 
  geom_boxplot(notch=TRUE, fill="grey") + 
  theme_science() + 
  xlab("AA at Position 1") + ylab("GFP RFU") 
ggsave("motif1_p4_aa1_exp1_RFU.tiff", height=4, width=5, unit="in")

m1_ile <- function(x) { 
  if(!grepl("AA[G|C|A]ATT...|.AA[G|C|A]ATT..|..AA[G|C|A]ATT.|...AA[G|C|A]ATT", x)) y <- "Absent"
  if(grepl("AA[G|C|A]ATT", x)) y <- "AAVATT"
  if(grepl("AA[G|C|A]ATC", x)) y <- "AAVATC"
  if(grepl("AA[G|C|A]ATA", x)) y <- "AAVATA"
  return(y)
}

All_reads_nostop$motif1_Ile <- sapply(All_reads_nostop$Sequence, m1_ile)

ggplot(All_reads_nostop, aes(x=reorder(motif1_Ile,motif1_Ile, function(x)-length(x)), Score, fill = reorder(motif1_Ile,motif1_Ile, function(x)-length(x)))) + 
  geom_boxplot(notch=TRUE) + scale_fill_manual(values=gray, guide = FALSE) + 
  theme_science() +
  xlab("Motif 1 Variants") + ylab("GFP Score")
ggsave("motif1_Ile_exp1.tiff", height=4, width=5, unit="in")

ggplot(All_reads_nostop, aes(x=reorder(motif1_Ile,motif1_Ile, function(x)-length(x)), RFU, fill = reorder(motif1_Ile,motif1_Ile, function(x)-length(x)))) + 
  geom_boxplot(notch=TRUE) + scale_fill_manual(values=gray, guide = FALSE) + 
  theme_science() + 
  xlab("Motif 1 Variants") + ylab("GFP RFU")
ggsave("motif1_Ile_exp1_RFU.tiff", height=4, width=5, unit="in")


m1_ile_p <- function(x) { 
  if(!grepl("AA[G|C|A]AT[A|T|C]...|.AA[G|C|A]AT[A|T|C]..|..AA[G|C|A]AT[A|T|C].|...AA[G|C|A]AT[A|T|C]", x)) y <- "Absent"
  if(grepl("AA[G|C|A]ATT...", x)) y <- "AAVATT Position 1"
  if(grepl("AA[G|C|A]ATC...", x)) y <- "AAVATC Position 1"
  if(grepl("AA[G|C|A]ATA...", x)) y <- "AAVATA Position 1"
  if(grepl("...AA[G|C|A]ATT", x)) y <- "AAVATT Position 4"
  if(grepl("...AA[G|C|A]ATC", x)) y <- "AAVATC Position 4"
  if(grepl("...AA[G|C|A]ATA", x)) y <- "AAVATA Position 4"
  if(grepl(".AA[G|C|A]ATT..", x)) y <- "NIF"
  if(grepl(".AA[G|C|A]ATC..", x)) y <- "NIF"
  if(grepl(".AA[G|C|A]ATA..", x)) y <- "NIF"
  if(grepl("..AA[G|C|A]ATT.", x)) y <- "NIF"
  if(grepl("..AA[G|C|A]ATC.", x)) y <- "NIF"
  if(grepl("..AA[G|C|A]ATA.", x)) y <- "NIF"
  return(y)
}

All_reads_nostop$motif1_ile_p <- sapply(All_reads_nostop$Sequence, m1_ile_p)

ggplot(subset(All_reads_nostop, !motif1_ile_p == "NIF"), aes(motif1_ile_p, Score)) + 
  geom_boxplot(notch=TRUE, fill="#009E73") + 
  theme_science() + scale_x_discrete(labels = function(x) str_wrap(x, width = 10), limits=c("Absent", "AAVATT Position 1", "AAVATC Position 1", "AAVATA Position 1", "AAVATT Position 4", "AAVATC Position 4", "AAVATA Position 4")) +
  xlab("Motif 1 Variants") + ylab("GFP Score")
ggsave("motif1_Ile_p_exp1.tiff", height=4, width=5, unit="in")

ggplot(subset(All_reads_nostop, !motif1_ile_p == "NIF"), aes(motif1_ile_p, RFU)) + 
  geom_boxplot(notch=TRUE, fill="#009E73") + 
  theme_science() + scale_x_discrete(labels = function(x) str_wrap(x, width = 10), limits=c("Absent", "AAVATT Position 1", "AAVATC Position 1", "AAVATA Position 1", "AAVATT Position 4", "AAVATC Position 4", "AAVATA Position 4")) +
  xlab("Motif 1 Variants") + ylab("GFP RFU")
ggsave("motif1_Ile_p_exp1_RFU.tiff", height=4, width=5, unit="in")

All_reads_nostop_simple <- data.frame(Sequence =  All_reads_nostop$Sequence, Score = All_reads_nostop$Score, RFU = All_reads_nostop$RFU)

Motif_1_p1_TTA <-subset(All_reads_nostop_simple, grepl("AA[G|C|A]TTA...", All_reads_nostop_simple$Sequence))
Motif_1_p1_TTG <-subset(All_reads_nostop_simple, grepl("AA[G|C|A]TTG...", All_reads_nostop_simple$Sequence))
Motif_1_p1_CTT <-subset(All_reads_nostop_simple, grepl("AA[G|C|A]CTT...", All_reads_nostop_simple$Sequence))
Motif_1_p1_CTG <-subset(All_reads_nostop_simple, grepl("AA[G|C|A]CTG...", All_reads_nostop_simple$Sequence))
Motif_1_p1_CTA <-subset(All_reads_nostop_simple, grepl("AA[G|C|A]CTA...", All_reads_nostop_simple$Sequence))
Motif_1_p1_CTC <-subset(All_reads_nostop_simple, grepl("AA[G|C|A]CTC...", All_reads_nostop_simple$Sequence))
Motif_1_p1_ATT <-subset(All_reads_nostop_simple, grepl("AA[G|C|A]ATT...", All_reads_nostop_simple$Sequence))
Motif_1_p4_TTA <-subset(All_reads_nostop_simple, grepl("...AA[G|C|A]TTA", All_reads_nostop_simple$Sequence))
Motif_1_p4_TTG <-subset(All_reads_nostop_simple, grepl("...AA[G|C|A]TTG", All_reads_nostop_simple$Sequence))
Motif_1_p4_CTT <-subset(All_reads_nostop_simple, grepl("...AA[G|C|A]CTT", All_reads_nostop_simple$Sequence))
Motif_1_p4_CTG <-subset(All_reads_nostop_simple, grepl("...AA[G|C|A]CTG", All_reads_nostop_simple$Sequence))
Motif_1_p4_CTA <-subset(All_reads_nostop_simple, grepl("...AA[G|C|A]CTA", All_reads_nostop_simple$Sequence))
Motif_1_p4_CTC <-subset(All_reads_nostop_simple, grepl("...AA[G|C|A]CTC", All_reads_nostop_simple$Sequence))
Motif_1_p4_ATT <-subset(All_reads_nostop_simple, grepl("...AA[G|C|A]ATT", All_reads_nostop_simple$Sequence))

Motif_1_p1_TTA$Mot_leu <-"AAVTTA Position 1"
Motif_1_p1_TTG$Mot_leu <-"AAVTTG Position 1"
Motif_1_p1_CTT$Mot_leu <-"AAVCTT Position 1"
Motif_1_p1_CTG$Mot_leu <-"AAVCTG Position 1"
Motif_1_p1_CTA$Mot_leu <-"AAVCTA Position 1"
Motif_1_p1_CTC$Mot_leu <-"AAVCTC Position 1"
Motif_1_p1_ATT$Mot_leu <-"AAVATT Position 1"
Motif_1_p4_TTA$Mot_leu <-"AAVTTA Position 4"
Motif_1_p4_TTG$Mot_leu <-"AAVTTG Position 4"
Motif_1_p4_CTT$Mot_leu <-"AAVCTT Position 4"
Motif_1_p4_CTG$Mot_leu <-"AAVCTG Position 4"
Motif_1_p4_CTA$Mot_leu <-"AAVCTA Position 4"
Motif_1_p4_CTC$Mot_leu <-"AAVCTC Position 4"
Motif_1_p4_ATT$Mot_leu <- "AAVATT Position 4"

Motif_Leu <- rbind(Motif_1_p1_TTA, Motif_1_p1_TTG, Motif_1_p1_CTT, Motif_1_p1_CTG, Motif_1_p1_CTA, Motif_1_p1_CTC, Motif_1_p1_ATT, Motif_1_p4_TTA, Motif_1_p4_TTG,
                                       Motif_1_p4_CTT, Motif_1_p4_CTG,
                                       Motif_1_p4_CTA,
                                       Motif_1_p4_CTC,
                                       Motif_1_p4_ATT)

ggplot(Motif_Leu, aes(Mot_leu, Score)) + 
  geom_boxplot(notch=TRUE, fill="#009E73") + 
  theme_science() + scale_x_discrete(labels = function(x) str_wrap(x, width = 10), limits=c("AAVATT Position 1", "AAVTTA Position 1", "AAVTTG Position 1", "AAVCTT Position 1", "AAVCTA Position 1", "AAVCTG Position 1", "AAVCTC Position 1", "AAVATT Position 4", "AAVTTA Position 4", "AAVTTG Position 4", "AAVCTT Position 4", "AAVCTA Position 4", "AAVCTG Position 4", "AAVCTC Position 4")) +
  xlab("Motif 1 Variants") + ylab("GFP Score")
ggsave("motif1_Leu_p_exp1.tiff", height=6, width=12, unit="in")

ggplot(Motif_Leu, aes(Mot_leu, RFU)) + 
  geom_boxplot(notch=TRUE, fill="#009E73") + 
  theme_science() + scale_x_discrete(labels = function(x) str_wrap(x, width = 10), limits=c("AAVATT Position 1", "AAVTTA Position 1", "AAVTTG Position 1", "AAVCTT Position 1", "AAVCTA Position 1", "AAVCTG Position 1", "AAVCTC Position 1", "AAVATT Position 4", "AAVTTA Position 4", "AAVTTG Position 4", "AAVCTT Position 4", "AAVCTA Position 4", "AAVCTG Position 4", "AAVCTC Position 4")) +
  xlab("Motif 1 Variants") + ylab("GFP RFU")
ggsave("motif1_Leu_p_exp1_RFU.tiff", height=6, width=12, unit="in")


m1_V_p <- function(x) { 
  if(!grepl("AA[G|C|A]ATT...|.AA[G|C|A]ATT..|..AA[G|C|A]ATT.|...AA[G|C|A]ATT", x)) y <- "Absent"
  if(grepl("AAGATT...", x)) y <- "AAGATT Position 1"
  if(grepl("AACATT...", x)) y <- "AACATT Position 1"
  if(grepl("AAAATT...", x)) y <- "AAAATT Position 1"
  if(grepl("...AAGATT", x)) y <- "AAGATT Position 4"
  if(grepl("...AACATT", x)) y <- "AACATT Position 4"
  if(grepl("...AAAATT", x)) y <- "AAAATT Position 4"
  if(grepl(".AAGATT..", x)) y <- "NIF"
  if(grepl(".AACATT..", x)) y <- "NIF"
  if(grepl(".AAAATT..", x)) y <- "NIF"
  if(grepl("..AAGATT.", x)) y <- "NIF"
  if(grepl("..AACATT.", x)) y <- "NIF"
  if(grepl("..AAAATT.", x)) y <- "NIF"
  return(y)
}

All_reads_nostop$motif1_v_p <- sapply(All_reads_nostop$Sequence,m1_V_p)

ggplot(subset(All_reads_nostop, !motif1_v_p == "NIF"), aes(x = reorder(motif1_v_p,motif1_v_p, function(x)-length(x)), Score, fill = reorder(motif1_v_p,motif1_v_p, function(x)-length(x)))) + 
  geom_boxplot(notch=TRUE) + scale_fill_manual(values=gray, guide = FALSE) +
  theme_science() + scale_x_discrete(labels = function(x) str_wrap(x, width = 10), limits=c("Absent", "AAGATT Position 1", "AACATT Position 1", "AAAATT Position 1", "AAGATT Position 4", "AACATT Position 4", "AAAATT Position 4")) +
  xlab("Motif 1 Variants") + ylab("GFP Score")
ggsave("motif1_v_p_exp1.tiff", height=4, width=6, unit="in")

ggplot(subset(All_reads_nostop, !motif1_v_p == "NIF"), aes(x = reorder(motif1_v_p,motif1_v_p, function(x)-length(x)), RFU, fill = reorder(motif1_v_p,motif1_v_p, function(x)-length(x)))) + 
  geom_boxplot(notch=TRUE) + scale_fill_manual(values=gray, guide = FALSE) +
  theme_science() + scale_x_discrete(labels = function(x) str_wrap(x, width = 10), limits=c("Absent", "AAGATT Position 1", "AACATT Position 1", "AAAATT Position 1", "AAGATT Position 4", "AACATT Position 4", "AAAATT Position 4")) +
  xlab("Motif 1 Variants") + ylab("GFP RFU")
ggsave("motif1_v_p_exp1_RFU.tiff", height=4, width=5, unit="in")


ggplot(All_reads_nostop, aes(motif2, Score, fill=motif2)) + 
  geom_boxplot(notch=TRUE) + scale_fill_manual(values=gray, guide = FALSE) +
  theme_science() + 
  xlab("Contains AADTAT") + ylab("GFP Score")
ggsave("motif2_exp1.tiff", height=4, width=4, unit="in")
ggplot(All_reads_nostop, aes(motif2_all, Score, fill = motif2_all)) + 
  geom_boxplot(notch=TRUE) + scale_fill_manual(values=gray, guide = FALSE) +
  theme_science() + 
  xlab("Position of AADTAT") + ylab("GFP Score")
ggsave("motif2_all_exp1.tiff", height=4, width=5, unit="in")

ggplot(All_reads_nostop, aes(motif2, RFU, fill= motif2)) + 
  geom_boxplot(notch=TRUE) + scale_fill_manual(values=gray, guide = FALSE) +
  theme_science() + 
  xlab("Contains AADTAT") + ylab("GFP RFU")
ggsave("motif2_exp1_RFU.tiff", height=4, width=4, unit="in")
ggplot(All_reads_nostop, aes(motif2_all, RFU, fill = motif2_all)) + 
  geom_boxplot(notch=TRUE) + scale_fill_manual(values=gray, guide = FALSE) +
  theme_science() + 
  xlab("Position of AADTAT") + ylab("GFP RFU")
ggsave("motif2_all_exp1_RFU.tiff", height=4, width=5, unit="in")

All_reads_nostop_motif2_p2 <- subset(All_reads_nostop, motif2_all == "Position 2")

m2_p2_n1 <- function(x) { 
  if(grepl("AAA[T|G|A]TAT..", x)) y <- "A"
  if(grepl("TAA[T|G|A]TAT..", x)) y <- "T"
  if(grepl("CAA[T|G|A]TAT..", x)) y <- "C"
  if(grepl("GAA[T|G|A]TAT..", x)) y <- "G"
  return(y)
}

All_reads_nostop_motif2_p2$motif2_p2_n1 <- sapply(All_reads_nostop_motif2_p2$Sequence, m2_p2_n1)

ggplot(All_reads_nostop_motif2_p2, aes(motif2_p2_n1, Score)) + 
  geom_boxplot(notch=TRUE, fill="grey") + 
  theme_science() + 
  xlab("Nucleotide at Position 1") + ylab("GFP Score")
ggsave("motif2_p2_n1_exp1.tiff", height=4, width=5, unit="in")

ggplot(All_reads_nostop_motif2_p2, aes(motif2_p2_n1, RFU)) + 
  geom_boxplot(notch=TRUE, fill="grey") + 
  theme_science() + 
  xlab("Nucleotide at Position 1") + ylab("GFP RFU")
ggsave("motif2_p2_n1_exp1_RFU.tiff", height=4, width=5, unit="in")

All_reads_nostop_motif2_p4 <- subset(All_reads_nostop, motif2_all == "Position 4")

AA_split <- colsplit(All_reads_nostop_motif2_p4$AA_character, "", c("AA_1", "AA_2", "AA_3"))

All_reads_nostop_motif2_p4 <- cbind(All_reads_nostop_motif2_p4, AA_split)

ggplot(All_reads_nostop_motif2_p4, aes(reorder(AA_1, -Score), Score)) + 
  geom_boxplot(notch=TRUE, fill="grey") + 
  theme_science() + 
  xlab("AA at Position 1") + ylab("GFP Score") 
ggsave("motif2_p4_aa1_exp1.tiff", height=4, width=5, unit="in")

ggplot(All_reads_nostop_motif2_p4, aes(reorder(AA_1, -RFU), RFU)) + 
  geom_boxplot(notch=TRUE, fill="grey") + 
  theme_science() + 
  xlab("AA at Position 1") + ylab("GFP RFU") 
ggsave("motif2_p4_aa1_exp1_RFU.tiff", height=4, width=5, unit="in")

m2_TAC <- function(x) { 
  if(!grepl("AA[T|G|A]TA[T|C]...|.AA[T|G|A]TA[T|C]..|..AA[T|G|A]TA[T|C].|...AA[T|G|A]TA[T|C]", x)) y <- "Absent"
  if(grepl("AA[T|G|A]TAT", x)) y <- "AADTAT"
  if(grepl("AA[T|G|A]TAC", x)) y <- "AADTAC"
  return(y)
}


All_reads_nostop$motif2_TAC <- sapply(All_reads_nostop$Sequence, m2_TAC)

ggplot(All_reads_nostop, aes(x = reorder(motif2_TAC, motif2_TAC, function(x)-length(x)), Score, fill = reorder(motif2_TAC, motif2_TAC, function(x)-length(x)))) + 
  geom_boxplot(notch=TRUE) + scale_fill_manual(values = gray, guide = FALSE) +
  theme_science() + 
  xlab("Motif 2 Variants") + ylab("GFP Score")
ggsave("motif2_TAC_exp1.tiff", height=4, width=5, unit="in")

ggplot(All_reads_nostop, aes(x = reorder(motif2_TAC, motif2_TAC, function(x)-length(x)), RFU, fill = reorder(motif2_TAC, motif2_TAC, function(x)-length(x)))) + 
  geom_boxplot(notch=TRUE) + scale_fill_manual(values = gray, guide = FALSE) + 
  theme_science() + 
  xlab("Motif 2 Variants") + ylab("GFP RFU")
ggsave("motif2_TAC_exp1_RFU.tiff", height=4, width=5, unit="in")

m2_TAC_p <- function(x) { 
  if(!grepl("AA[T|G|A]TA[T|C]...|.AA[T|G|A]TA[T|C]..|..AA[T|G|A]TA[T|C].|...AA[T|G|A]TA[T|C]", x)) y <- "Absent"
  if(grepl("AA[T|G|A]TAT...", x)) y <- "AADTAT Position 1"
  if(grepl("AA[T|G|A]TAC...", x)) y <- "AADTAC Position 1"
  if(grepl("...AA[T|G|A]TAT", x)) y <- "AADTAT Position 4"
  if(grepl("...AA[T|G|A]TAC", x)) y <- "AADTAC Position 4"
  if(grepl(".AA[T|G|A]TAT..", x)) y <- "NIF"
  if(grepl(".AA[T|G|A]TAC..", x)) y <- "NIF"
  if(grepl("..AA[T|G|A]TAT.", x)) y <- "NIF"
  if(grepl("..AA[T|G|A]TAC.", x)) y <- "NIF"
  if(grepl("..AA[T|G|A]TAT.", x)) y <- "NIF"
  if(grepl("..AA[T|G|A]TAC.", x)) y <- "NIF"
  if(grepl(".AA[T|G|A]TAT..", x)) y <- "NIF"
  if(grepl(".AA[T|G|A]TAC..", x)) y <- "NIF"
  return(y)
}

All_reads_nostop$motif2_TAC_p <- sapply(All_reads_nostop$Sequence, m2_TAC_p)

ggplot(subset(All_reads_nostop, !motif2_TAC_p == "NIF"), aes(motif2_TAC_p, Score)) + 
  geom_boxplot(notch=TRUE, fill="#F0E442") + 
  theme_science() + scale_x_discrete(labels = function(x) str_wrap(x, width = 10), limits=c("Absent", "AADTAT Position 1", "AADTAC Position 1", "AADTAT Position 4", "AADTAC Position 4")) +
  xlab("Motif 2 Variants") + ylab("GFP Score")
ggsave("motif2_TAC_p_exp1.tiff", height=4, width=5, unit="in")

ggplot(subset(All_reads_nostop, !motif2_TAC_p == "NIF"), aes(motif2_TAC_p, RFU)) + 
  geom_boxplot(notch=TRUE, fill="#F0E442") + 
  theme_science() + scale_x_discrete(labels = function(x) str_wrap(x, width = 10), limits=c("Absent", "AADTAT Position 1", "AADTAC Position 1", "AADTAT Position 4", "AADTAC Position 4")) +
  xlab("Motif 2 Variants") + ylab("GFP RFU")
ggsave("motif2_TAC_p_exp1_RFU.tiff", height=4, width=5, unit="in")


Motif_2_p1_TTA <-subset(All_reads_nostop_simple, grepl("AA[T|G|A]TTA...", All_reads_nostop_simple$Sequence))
Motif_2_p1_TTG <-subset(All_reads_nostop_simple, grepl("AA[T|G|A]TTG...", All_reads_nostop_simple$Sequence))
Motif_2_p1_CTT <-subset(All_reads_nostop_simple, grepl("AA[T|G|A]CTT...", All_reads_nostop_simple$Sequence))
Motif_2_p1_CTG <-subset(All_reads_nostop_simple, grepl("AA[T|G|A]CTG...", All_reads_nostop_simple$Sequence))
Motif_2_p1_CTA <-subset(All_reads_nostop_simple, grepl("AA[T|G|A]CTA...", All_reads_nostop_simple$Sequence))
Motif_2_p1_CTC <-subset(All_reads_nostop_simple, grepl("AA[T|G|A]CTC...", All_reads_nostop_simple$Sequence))
Motif_2_p1_ATT <-subset(All_reads_nostop_simple, grepl("AA[T|G|A]TAT...", All_reads_nostop_simple$Sequence))
Motif_2_p4_TTA <-subset(All_reads_nostop_simple, grepl("...AA[T|G|A]TTA", All_reads_nostop_simple$Sequence))
Motif_2_p4_TTG <-subset(All_reads_nostop_simple, grepl("...AA[T|G|A]TTG", All_reads_nostop_simple$Sequence))
Motif_2_p4_CTT <-subset(All_reads_nostop_simple, grepl("...AA[T|G|A]CTT", All_reads_nostop_simple$Sequence))
Motif_2_p4_CTG <-subset(All_reads_nostop_simple, grepl("...AA[T|G|A]CTG", All_reads_nostop_simple$Sequence))
Motif_2_p4_CTA <-subset(All_reads_nostop_simple, grepl("...AA[T|G|A]CTA", All_reads_nostop_simple$Sequence))
Motif_2_p4_CTC <-subset(All_reads_nostop_simple, grepl("...AA[T|G|A]CTC", All_reads_nostop_simple$Sequence))
Motif_2_p4_ATT <-subset(All_reads_nostop_simple, grepl("...AA[T|G|A]TAT", All_reads_nostop_simple$Sequence))

Motif_2_p1_TTA$Mot_leu <-"AADTTA Position 1"
Motif_2_p1_TTG$Mot_leu <-"AADTTG Position 1"
Motif_2_p1_CTT$Mot_leu <-"AADCTT Position 1"
Motif_2_p1_CTG$Mot_leu <-"AADCTG Position 1"
Motif_2_p1_CTA$Mot_leu <-"AADCTA Position 1"
Motif_2_p1_CTC$Mot_leu <-"AADCTC Position 1"
Motif_2_p1_ATT$Mot_leu <-"AADTAT Position 1"
Motif_2_p4_TTA$Mot_leu <-"AADTTA Position 4"
Motif_2_p4_TTG$Mot_leu <-"AADTTG Position 4"
Motif_2_p4_CTT$Mot_leu <-"AADCTT Position 4"
Motif_2_p4_CTG$Mot_leu <-"AADCTG Position 4"
Motif_2_p4_CTA$Mot_leu <-"AADCTA Position 4"
Motif_2_p4_CTC$Mot_leu <-"AADCTC Position 4"
Motif_2_p4_ATT$Mot_leu <- "AADTAT Position 4"

Motif_2_Leu <- rbind(Motif_2_p1_TTA, Motif_2_p1_TTG, Motif_2_p1_CTT, Motif_2_p1_CTG, Motif_2_p1_CTA, Motif_2_p1_CTC, Motif_2_p1_ATT, Motif_2_p4_TTA, Motif_2_p4_TTG,
                   Motif_2_p4_CTT, Motif_2_p4_CTG,
                   Motif_2_p4_CTA,
                   Motif_2_p4_CTC,
                   Motif_2_p4_ATT)

ggplot(Motif_2_Leu, aes(Mot_leu, Score)) + 
  geom_boxplot(notch=TRUE, fill="#F0E442") + 
  theme_science() + scale_x_discrete(labels = function(x) str_wrap(x, width = 10), limits=c("AADTAT Position 1", "AADTTA Position 1", "AADTTG Position 1", "AADCTT Position 1", "AADCTA Position 1", "AADCTG Position 1", "AADCTC Position 1", "AADTAT Position 4", "AADTTA Position 4", "AADTTG Position 4", "AADCTT Position 4", "AADCTA Position 4", "AADCTG Position 4", "AADCTC Position 4")) +
  xlab("Motif 2 Variants") + ylab("GFP Score")
ggsave("motif2_Leu_p_exp1.tiff", height=6, width=12, unit="in")

ggplot(Motif_2_Leu, aes(Mot_leu, RFU)) + 
  geom_boxplot(notch=TRUE, fill="#F0E442") + 
  theme_science() + scale_x_discrete(labels = function(x) str_wrap(x, width = 10), limits=c("AADTAT Position 1", "AADTTA Position 1", "AADTTG Position 1", "AADCTT Position 1", "AADCTA Position 1", "AADCTG Position 1", "AADCTC Position 1", "AADTAT Position 4", "AADTTA Position 4", "AADTTG Position 4", "AADCTT Position 4", "AADCTA Position 4", "AADCTG Position 4", "AADCTC Position 4")) +
  xlab("Motif 2 Variants") + ylab("GFP RFU")
ggsave("motif2_Leu_p_exp1_RFU.tiff", height=6, width=12, unit="in")


m2_d_p <- function(x) { 
  if(!grepl("AA[G|T|A]TAT...|.AA[G|T|A]TAT..|..AA[G|T|A]TAT.|...AA[G|T|A]TAT", x)) y <- "Absent"
  if(grepl("AAGTAT...", x)) y <- "AAGTAT Position 1"
  if(grepl("AATTAT...", x)) y <- "AATTAT Position 1"
  if(grepl("AAATAT...", x)) y <- "AAATAT Position 1"
  if(grepl("...AAGTAT", x)) y <- "AAGTAT Position 4"
  if(grepl("...AATTAT", x)) y <- "AATTAT Position 4"
  if(grepl("...AAATAT", x)) y <- "AAATAT Position 4"
  if(grepl(".AAGTAT..", x)) y <- "NIF"
  if(grepl(".AATTAT..", x)) y <- "NIF"
  if(grepl(".AAATAT..", x)) y <- "NIF"
  if(grepl("..AAGTAT.", x)) y <- "NIF"
  if(grepl("..AATTAT.", x)) y <- "NIF"
  if(grepl("..AAATAT.", x)) y <- "NIF"
  return(y)
}

All_reads_nostop$motif2_d_p <- sapply(All_reads_nostop$Sequence,m2_d_p)

ggplot(subset(All_reads_nostop, !motif2_d_p == "NIF"), aes(x = reorder(motif2_d_p,motif1_v_p, function(x)-length(x)), Score, fill = reorder(motif2_d_p,motif1_v_p, function(x)-length(x)))) + 
  geom_boxplot(notch=TRUE) + scale_fill_manual(values = gray, guide = FALSE) +
  theme_science() + scale_x_discrete(labels = function(x) str_wrap(x, width = 10), limits=c("Absent", "AAGTAT Position 1", "AATTAT Position 1", "AAATAT Position 1", "AAGTAT Position 4", "AATTAT Position 4", "AAATAT Position 4")) +
  xlab("Motif 2 Variants") + ylab("GFP Score")
ggsave("motif2_v_p_exp1.tiff", height=4, width=6, unit="in")

ggplot(subset(All_reads_nostop, !motif2_d_p == "NIF"), aes(x = reorder(motif2_d_p,motif1_v_p, function(x)-length(x)), RFU, fill = reorder(motif2_d_p,motif1_v_p, function(x)-length(x)))) + 
  geom_boxplot(notch=TRUE) + scale_fill_manual(values = gray, guide = FALSE) +
  theme_science() + scale_x_discrete(labels = function(x) str_wrap(x, width = 10), limits=c("Absent", "AAGTAT Position 1", "AATTAT Position 1", "AAATAT Position 1", "AAGTAT Position 4", "AATTAT Position 4", "AAATAT Position 4")) +
  xlab("Motif 2 Variants") + ylab("GFP RFU")
ggsave("motif2_v_p_exp1_RFU.tiff", height=4, width=5, unit="in")

ggplot(All_reads_nostop, aes(motif3, Score)) + geom_boxplot(notch=TRUE) + theme_science()
ggsave("motif3_exp1.tiff")
ggplot(All_reads_nostop, aes(motif3_all, Score)) + geom_boxplot(notch=TRUE) + theme_science()
ggsave("motif3_all_exp1.tiff")

ggplot(All_reads_nostop, aes(motif3, RFU)) + geom_boxplot(notch=TRUE) + theme_science()
ggsave("motif3_exp1_RFU.tiff")
ggplot(All_reads_nostop, aes(motif3_all, RFU)) + geom_boxplot(notch=TRUE) + theme_science()
ggsave("motif3_all_exp1_RFU.tiff")

ggplot(All_reads_nostop, aes(motif4, Score)) + geom_boxplot(notch=TRUE) + theme_science()
ggsave("motif4_exp1.tiff")
ggplot(All_reads_nostop, aes(motif4_all, Score)) + geom_boxplot(notch=TRUE) + theme_science()
ggsave("motif4_all_exp1.tiff")

All_reads_nostop$KNIY <- grepl("[K|N][I|Y].|.[K|N][I|Y]", All_reads_nostop$AA)

ggplot(All_reads_nostop, aes(KNIY, Score, fill = KNIY)) + 
  xlab("K|N-I|Y Presence") + ylab("GFP Score") +
  geom_boxplot(notch=TRUE)+ scale_fill_manual(values=gray, guide = FALSE) + theme_science()
ggsave("KNIY_exp1.tiff", height=4, width=4, units="in")

ggplot(All_reads_nostop, aes(KNIY, RFU, fill = KNIY)) + 
  xlab("K|N-I|Y Presence") + ylab("GFP RFU") + 
  geom_boxplot(notch=TRUE) + scale_fill_manual(values=gray, guide = FALSE) + theme_science()
ggsave("KNIY_exp1_RFU.tiff", height=4, width=4, units="in")

KNIY <- function(x) { 
  if(!grepl("[K|N][I|Y].|.[K|N][I|Y]", x)) y <- "Absent"
  if(grepl("[K|N][I|Y].", x)) y <- "Codons 3 and 4"
  if(grepl(".[K|N][I|Y]", x)) y <- "Codons 4 and 5"
  return(y)
}

All_reads_nostop$KNIY_all <- sapply(All_reads_nostop$AA, KNIY)

ggplot(All_reads_nostop, aes(KNIY_all, Score, fill = KNIY_all)) + 
  xlab("K|N-I|Y Position") + ylab("GFP Score") +
  geom_boxplot(notch=TRUE) + scale_fill_manual(values=gray, guide = FALSE) + theme_science()
ggsave("KNIY_all_exp1.tiff", height=4, width=5, units="in")

ggplot(All_reads_nostop, aes(KNIY_all, RFU, fill = KNIY_all)) + 
  xlab("K|N-I|Y Position") + ylab("GFP RFU") +
  geom_boxplot(notch=TRUE) + scale_fill_manual(values=gray, guide = FALSE) + theme_science()
ggsave("KNIY_all_exp1_RFU.tiff", height=4, width=5, units="in")

t.test(All_reads_nostop$Score[All_reads_nostop$KNIY == "TRUE"], All_reads_nostop$Score[All_reads_nostop$KNIY == "FALSE"])

All_reads_nostop$IYKN <- grepl("[I|Y][K|N].|.[I|Y][K|N]", All_reads_nostop$AA)

ggplot(All_reads_nostop, aes(IYKN, Score)) + geom_boxplot(notch=TRUE) + theme_science()
ggsave("IYKN_exp1.tiff")

ggplot(All_reads_nostop, aes(IYKN, RFU)) + geom_boxplot(notch=TRUE) + theme_science()
ggsave("IYKN_exp1_RFU.tiff")

t.test(All_reads_nostop$Score[All_reads_nostop$IYKN == "TRUE"], All_reads_nostop$Score[All_reads_nostop$IYKN == "FALSE"])


high_aa <- high_sequences$AA
low_aa <- low_sequences$AA
write.table(high_aa, "/Users/Kyle/Dropbox/TI-seq/NextSeq/high_aa_exp1.fa", row.names = FALSE, col.names = FALSE)

high_aa_test <- read.sequences("/Users/Kyle/Dropbox/TI-seq/NextSeq/high_aa_exp1.fa")

motifs_aa <- motifModel(high_aa_test)
print(motifs_aa) 
plotPositions(motifs_aa)
png(filename="aa_frequency_high_exp1.png")
plotPositions(motifs_aa)
dev.off()

motifs_aa_table <- as.data.frame(motifs_aa@mmodel)
motifs_aa_table$aa <- row.names(motifs_aa_table)
motifs_aa_table <- motifs_aa_table[!(motifs_aa_table$aa %in% c("X", "-", "B", "X", "Z")),]
ggplot(motifs_aa_table, aes(aa, V1)) + geom_bar(stat='identity') + theme_science() + xlab("Amino Acid at Position 1") + ylab("Frequency")
ggsave("aa_frequency_high_p1_exp1.tiff")
ggplot(motifs_aa_table, aes(aa, V2)) + geom_bar(stat='identity') + theme_science() + xlab("Amino Acid at Position 2") + ylab("Frequency")
ggsave("aa_frequency_high_p2_exp1.tiff")
ggplot(motifs_aa_table, aes(aa, V3)) + geom_bar(stat='identity') + theme_science() + xlab("Amino Acid at Position 3") + ylab("Frequency")
ggsave("aa_frequency_high_p3_exp1.tiff")

write.table(low_aa, "/Users/Kyle/Dropbox/TI-seq/NextSeq/low_aa_exp1.fa", row.names = FALSE, col.names = FALSE)

low_aa_test <- read.sequences("/Users/Kyle/Dropbox/TI-seq/NextSeq/low_aa_exp1.fa")

motifs_aa_low <- motifModel(low_aa_test)
print(motifs_aa_low) 
plotPositions(motifs_aa_low)
png(filename="aa_frequency_low_exp1.png")
plotPositions(motifs_aa_low)
dev.off()

motifs_aa_low_table <- as.data.frame(motifs_aa_low@mmodel)
motifs_aa_low_table$aa <- row.names(motifs_aa_low_table)
motifs_aa_low_table <- motifs_aa_low_table[!(motifs_aa_low_table$aa %in% c("X", "-", "B", "X", "Z")),]
ggplot(motifs_aa_low_table, aes(aa, V1)) + geom_bar(stat='identity') + theme_science() + xlab("Amino Acid at Position 1") + ylab("Frequency")
ggsave("aa_frequency_low_p1_exp1.tiff")
ggplot(motifs_aa_low_table, aes(aa, V2)) + geom_bar(stat='identity') + theme_science() + xlab("Amino Acid at Position 2") + ylab("Frequency")
ggsave("aa_frequency_low_p2_exp1.tiff")
ggplot(motifs_aa_low_table, aes(aa, V3)) + geom_bar(stat='identity') + theme_science() + xlab("Amino Acid at Position 3") + ylab("Frequency")
ggsave("aa_frequency_low_p3_exp1.tiff")



All_reads_nostop$aa_high <- grepl("[K|N][I|L|S][R|I|L]", All_reads_nostop$AA)

ggplot(All_reads_nostop, aes(aa_high, Score)) + geom_boxplot(notch=TRUE) + theme_science()
ggsave("aa_high_exp1.tiff")

ggplot(All_reads_nostop, aes(aa_high, RFU)) + geom_boxplot(notch=TRUE) + theme_science()
ggsave("aa_high_exp1_RFU.tiff")

All_reads_nostop$KLI <- grepl("KLI", All_reads_nostop$AA)

ggplot(All_reads_nostop, aes(KLI, Score)) + geom_boxplot(notch=TRUE) + theme_science()
ggsave("KLI_exp1.tiff")

ggplot(All_reads_nostop, aes(KLI, RFU)) + geom_boxplot(notch=TRUE) + theme_science()
ggsave("KLI_exp1_RFU.tiff")

All_reads_nostop$RRR <- grepl("RRR", All_reads_nostop$AA)

ggplot(All_reads_nostop, aes(RRR, Score)) + geom_boxplot(notch=TRUE) + theme_science()
ggsave("RRR_exp1.tiff")

ggplot(All_reads_nostop, aes(RRR, RFU)) + geom_boxplot(notch=TRUE) + theme_science()
ggsave("RRR_exp1_RFU.tiff")



#####Experiment 2 (AKA Experiment 3)


Bin1_e2 <- read.table("/Users/Kyle/Dropbox/TI-seq/NextSeq/Exp_3_01_counts")
Bin2_e2 <- read.table("/Users/Kyle/Dropbox/TI-seq/NextSeq/Exp_3_02_counts")
Bin3_e2 <- read.table("/Users/Kyle/Dropbox/TI-seq/NextSeq/Exp_3_03_counts")
Bin4_e2 <- read.table("/Users/Kyle/Dropbox/TI-seq/NextSeq/Exp_3_04_counts")
Bin5_e2 <- read.table("/Users/Kyle/Dropbox/TI-seq/NextSeq/Exp_3_05_counts")


colnames(Plasmid) <- c("Sequence", "Reads_Plasmid")
colnames(Bin1_e2) <- c("Reads_B1", "Sequence")
colnames(Bin2_e2) <- c("Reads_B2", "Sequence")
colnames(Bin3_e2) <- c("Reads_B3", "Sequence")
colnames(Bin4_e2) <- c("Reads_B4", "Sequence")
colnames(Bin5_e2) <- c("Reads_B5", "Sequence")

Plasmid_bin1_e2 <- merge(Plasmid, Bin1_e2, by="Sequence", all = TRUE)
Plasmid_bin1_2_e2  <- merge(Plasmid_bin1, Bin2_e2, by="Sequence", all = TRUE)
Plasmid_bin1_2_3_e2  <- merge(Plasmid_bin1_2, Bin3_e2, by="Sequence", all = TRUE)
Plasmid_bin1_2_3_4_e2 <- merge(Plasmid_bin1_2_3, Bin4_e2, by="Sequence", all = TRUE)
Plasmid_bin1_2_3_4_5_e2 <- merge(Plasmid_bin1_2_3_4, Bin5_e2, by="Sequence", all = TRUE)

All_reads_e2 <- Plasmid_bin1_2_3_4_5_e2

All_reads_e2$DNA <- DNAStringSet(All_reads_e2$Sequence)

All_reads_e2$Reads_Plasmid[is.na(All_reads_e2$Reads_Plasmid)] <- 0
All_reads_e2$Reads_B1[is.na(All_reads_e2$Reads_B1)] <- 0
All_reads_e2$Reads_B2[is.na(All_reads_e2$Reads_B2)] <- 0
All_reads_e2$Reads_B3[is.na(All_reads_e2$Reads_B3)] <- 0
All_reads_e2$Reads_B4[is.na(All_reads_e2$Reads_B4)] <- 0
All_reads_e2$Reads_B5[is.na(All_reads_e2$Reads_B5)] <- 0

All_reads_e2$total_bins <- All_reads_e2$Reads_B1 + All_reads_e2$Reads_B2 + All_reads_e2$Reads_B3 + All_reads_e2$Reads_B4 + All_reads_e2$Reads_B5

write.table(All_reads_e2, "all_reads_e2.txt")

All_reads_e2 <- subset(All_reads_e2, !total_bins < 10)

sum(All_reads_e2$total_bins)
mean(All_reads_e2$total_bins)

ggplot(All_reads_e2, aes(log10(total_bins))) + geom_histogram() 
ggsave("hist_total_counts_exp2.tiff")

All_reads_e2$B1_ratio <- All_reads_e2$Reads_B1/All_reads_e2$total_bins
All_reads_e2$B2_ratio <- All_reads_e2$Reads_B2/All_reads_e2$total_bins
All_reads_e2$B3_ratio <- All_reads_e2$Reads_B3/All_reads_e2$total_bins
All_reads_e2$B4_ratio <- All_reads_e2$Reads_B4/All_reads_e2$total_bins
All_reads_e2$B5_ratio <- All_reads_e2$Reads_B5/All_reads_e2$total_bins

All_reads_e2$Score <- (All_reads_e2$B1_ratio*1) + (All_reads_e2$B2_ratio*2) + (All_reads_e2$B3_ratio*3) + (All_reads_e2$B4_ratio*4) + (All_reads_e2$B5_ratio*5)
All_reads_e2$RFU <- (All_reads_e2$B1_ratio*20) + (All_reads_e2$B2_ratio*120) + (All_reads_e2$B3_ratio*600) + (All_reads_e2$B4_ratio*3600) + (All_reads_e2$B5_ratio*12000)


tai <- read.table("/Users/Kyle/Dropbox/TI-seq/NextSeq/cds_lib_tai.txt")

tai <- data.frame(Sequence = tai$sequence, tAI = tai$sequences.tai, AA = tai$AA)
All_reads_e2 <- merge(All_reads_e2, tai, by="Sequence")



All_reads_e2$AA_character <- as.character(All_reads_e2$AA)

All_reads_e2$Stop <- grepl("......TAA|......TGA|......TAG|...TAA...|...TAG...|...TGA...|TAA......|TGA......|TAG......", All_reads_e2$Sequence)
All_reads_e2$Amber <- grepl("......TAG|...TAG...|TAG......", All_reads_e2$Sequence)
All_reads_e2$Ochre <- grepl("......TAA|...TAA...|TAA......", All_reads_e2$Sequence)
All_reads_e2$Opal <- grepl("......TGA|...TGA...|TGA......", All_reads_e2$Sequence)

stop <- function(x) { 
  if(!grepl("......TAA|......TGA|......TAG|...TAA...|...TAG...|...TGA...|TAA......|TGA......|TAG......", x)) y <- "Absent"
  if(grepl("......TAG|...TAG...|TAG......", x)) y <- "Amber (UAG)"
  if(grepl("......TAA|...TAA...|TAA......", x)) y <- "Ochre (UAA)"
  if(grepl("......TGA|...TGA...|TGA......", x)) y <- "Opal (UGA"
  return(y)
}

All_reads_e2$Stop_codon <- sapply(All_reads_e2$Sequence, stop)

All_reads_e2$Rare <- grepl("......AGG|...AGG...|AGG......|......AGA|...AGA...|AGA......", All_reads_e2$Sequence)
All_reads_e2$Rare_2 <- grepl("AGGAGA...|AGG...AGG|AGG...AGA|AGGAGG...|AGAAGG...|AGAAGA...|AGA...AGG|AGA...AGA|...AGGAGG|...AGGAGA|...AGAAGG|...AGAAGA", All_reads_e2$Sequence)
All_reads_e2$Rare_3 <- grepl("AGGAGGAGG|AGGAGGAGA|AGGAGAAGG|AGGAGAAGA|AGAAGGAGG|AGAAGGAGA|AGAAGAAGG|AGAAGAAGA", All_reads_e2$Sequence)

All_reads_e2$RLI <- grepl(paste(rare_codons$V1, collapse = "|"), All_reads_e2$Sequence)

All_reads_e2$xho <- grepl("CTCGAG...|.CTCGAG..|..CTCGAG.|...CTCGAG", All_reads_e2$Sequence)
All_reads_e2$nco <- grepl("CCATGG...|.CCATGG..|..CCATGG.|...CCATGG", All_reads_e2$Sequence)

All_reads_e2 <- subset(All_reads_e2, xho == "FALSE")
All_reads_e2 <- subset(All_reads_e2, nco == "FALSE")

All_reads_e2$hydrophobicity <- hydrophobicity(All_reads_e2$AA, scale = "KyteDoolittle")
All_reads_e2$charge <- pI(All_reads_e2$AA, pKscale = "EMBOSS")

write.table(All_reads_e2, "all_reads_e2_annotated.txt")

ggplot(All_reads_e2, aes(Stop, Score)) + 
  geom_violin(fill="grey") + 
  theme_science() + xlab("Contains Stop Codon") + 
  ylab("GFP Score")
ggsave("Stop_exp2.tiff", height=4, width=4, units="in", dpi=300)

ggplot(All_reads_e2, aes(Stop, RFU)) + 
  geom_violin(fill="grey") + 
  theme_science() + xlab("Contains Stop Codon") + 
  ylab("GFP RFU")
ggsave("Stop_exp2_RFU.tiff", height=4, width=4, units="in", dpi=300)

ggplot(All_reads_e2, aes(Stop_codon, Score)) + 
  geom_violin(fill="grey") + 
  theme_science() + xlab("Stop Codon Type") + 
  ylab("GFP Score")
ggsave("Stop_codon_exp2.tiff", height=4, width=4, units="in", dpi=300)

ggplot(All_reads_e2, aes(Stop_codon, RFU)) + 
  geom_violin(fill="grey") + 
  theme_science() + xlab("Stop Codon Type") + 
  ylab("GFP RFU")
ggsave("Stop_codon_exp2_RFU.tiff", height=4, width=4, units="in", dpi=300)



amber_gln <- function(x) { 
  if(!grepl("......TAG|...TAG...|TAG......|......CAG|...CAG...|CAG......", x)) y <- "Absent"
  if(grepl("TAG......", x)) y <- "Amber Codon 1"
  if(grepl("...TAG...", x)) y <- "Amber Codon 2"
  if(grepl("......TAG", x)) y <- "Amber Codon 3"
  if(grepl("CAG......", x)) y <- "CAG Codon 1"
  if(grepl("...CAG...", x)) y <- "CAG Codon 2"
  if(grepl("......CAG", x)) y <- "CAG Codon 3"
  return(y)
}

All_reads_e2$Amber <- sapply(All_reads_e2$Sequence, amber_gln)

ggplot(All_reads_e2, aes(Amber, Score)) + 
  geom_violin(fill="grey") + 
  theme_science() + xlab("") +
  ylab("GFP Score") + theme(axis.text.x=element_text(angle = 45, hjust = 1))
ggsave("amber_gln_exp2.tiff", height=4, width=4, units="in", dpi=300)

ggplot(All_reads_e2, aes(Amber, RFU)) + 
  geom_violin(fill="grey") + 
  theme_science() + xlab("") +
  ylab("GFP RFU") + theme(axis.text.x=element_text(angle = 45, hjust = 1))
ggsave("amber_gln_exp2_RFU.tiff", height=4, width=4, units="in", dpi=300)


ggplot(All_reads_e2, aes(Amber, Score)) + 
  geom_jitter(alpha=0.2) + 
  theme_science() + xlab("") + 
  ylab("GFP Score") + theme(axis.text.x=element_text(angle = 45, hjust = 1))
ggsave("amber_gln_exp2_jitter.tiff", height=4, width=4, units="in", dpi=300)

ggplot(All_reads_e2, aes(Amber, RFU)) + 
  geom_jitter(alpha=0.2) + 
  theme_science() + xlab("") + 
  ylab("GFP RFU") + theme(axis.text.x=element_text(angle = 45, hjust = 1))
ggsave("amber_gln_exp2_jitter_RFU.tiff", height=4, width=4, units="in", dpi=300)


All_reads_e2_nostop <- subset(All_reads_e2, Stop == "FALSE")

ggplot(All_reads_e2_nostop, aes(Score)) + 
  geom_density(fill="#56B4E9") + 
  xlab("GFP Score") + 
  ylab("Density") + theme_science() + 
  scale_y_continuous(limits=c(0,0.84), expand = c(0,0))
ggsave("density_nostop_exp2.tiff", height =4, width=4, units="in")

ggplot(All_reads_e2_nostop, aes(RFU)) + 
  geom_density(fill="#56B4E9") + 
  xlab("GFP RFU") + 
  ylab("Density") + theme_science()  +
  scale_y_continuous(limits=c(0,0.00032), expand = c(0,0))
ggsave("density_nostop_exp2_RFU.tiff", height =4, width=4, units="in")

ggplot(All_reads_e2_nostop, aes(Score)) + geom_histogram() + theme_science()
ggsave("hist_nostop_exp2.tiff")

ggplot(All_reads_e2_nostop, aes(RFU)) + geom_histogram() + theme_science()
ggsave("hist_nostop_exp2_RFU.tiff")

ggplot(All_reads_e2, aes(Score)) + 
  geom_density(fill="#E69F00") + 
  theme_science() + 
  xlab("GFP Score") + 
  ylab("Density") + 
  theme_science() + 
  scale_y_continuous(limits=c(0,0.84), expand = c(0,0))
ggsave("density_exp2.tiff", height =4, width=4, units="in")

ggplot(All_reads_e2, aes(RFU)) + 
  geom_density(fill="#E69F00") + 
  theme_science() + 
  xlab("GFP Score") + 
  ylab("Density") + 
  theme_science() + 
  scale_y_continuous(limits=c(0,0.00032), expand = c(0,0))
ggsave("density_exp2_RFU.tiff", height =4, width=4, units="in")

ggplot(All_reads_e2, aes(Score)) + geom_histogram() + theme_science()
ggsave("hist_exp2.tiff")

ggplot(All_reads_e2, aes(RFU)) + geom_histogram() + theme_science()
ggsave("hist_exp2_RFU.tiff")

ggplot(All_reads_e2, aes(Score, group=Stop, colour=Stop)) + geom_density() + theme_science()
ggplot(All_reads_e2, aes(Score, group=Stop_codon, colour=Stop_codon)) + 
  geom_density(size=2) + 
  theme_science() +
  xlab("GFP Score") + 
  ylab("Density") + 
  scale_colour_manual(values=cbPalette) +
  scale_y_continuous(limits=c(0,1.2), expand = c(0,0))
ggsave("hist_all_exp2.tiff")

ggplot(All_reads_e2, aes(RFU, group=Stop_codon, colour=Stop_codon)) + 
  geom_density(size=2) + 
  theme_science() +
  xlab("GFP RFU") + 
  ylab("Density") + 
  scale_colour_manual(values=cbPalette) +
  scale_y_continuous(limits=c(0,0.0008), expand = c(0,0))
ggsave("hist_all_exp2_RFU.tiff")



ggplot() + geom_histogram(data=All_reads_e2, aes(Score), fill="#E69F00", alpha =0.5) + geom_histogram(data=All_reads_e2_nostop, aes(Score), fill="#56B4E9", alpha=0.5) + theme_science()
ggplot() + geom_density(data=All_reads_e2, aes(Score), fill="#E69F00", colour="#E69F00", alpha=0.8) + 
  geom_density(data=All_reads_e2_nostop, aes(Score), fill="#56B4E9", colour="#56B4E9", alpha=0.8) + 
  xlab("GFP Score") + ylab("Density") + theme_science() + scale_y_continuous(limits=c(0,0.84), expand = c(0,0))

write.table(All_reads_e2_nostop, "All_reads_e2_nostop.txt")

ggplot(All_reads_e2_nostop, aes(Rare, Score)) + geom_boxplot(notch=TRUE)+ theme_science()
ggsave("rare_exp2.tiff")
ggplot(All_reads_e2_nostop, aes(Rare_2, Score)) + geom_boxplot(notch=TRUE)+ theme_science()
ggsave("rare_2_exp2.tiff")
ggplot(All_reads_e2_nostop, aes(Rare_3, Score)) + geom_boxplot(notch=TRUE)+ theme_science()
ggsave("rare_3_exp2.tiff")

ggplot(All_reads_e2_nostop, aes(Rare, RFU)) + geom_boxplot(notch=TRUE)+ theme_science()
ggsave("rare_exp2_RFU.tiff")
ggplot(All_reads_e2_nostop, aes(Rare_2, RFU)) + geom_boxplot(notch=TRUE)+ theme_science()
ggsave("rare_2_exp2_RFU.tiff")
ggplot(All_reads_e2_nostop, aes(Rare_3, RFU)) + geom_boxplot(notch=TRUE)+ theme_science()
ggsave("rare_3_exp2_RFU.tiff")

ggplot(All_reads_e2_nostop, aes(RLI, Score)) + geom_boxplot(notch=TRUE) + scale_x_discrete(labels = c("Absent", "Present")) + xlab("3 Rare RLI Codons") + ylab("GFP Score") + theme_science() 
ggsave("RLI_exp2.tiff", height = 3, width = 3, units= "in")
t.test(All_reads_e2_nostop$Score[All_reads_e2_nostop$RLI == TRUE], All_reads_e2_nostop$Score[All_reads_e2_nostop$RLI == FALSE])

cor.test(All_reads_e2_nostop$Score, All_reads_e2_nostop$hydrophobicity)
ggplot(All_reads_e2_nostop, aes(hydrophobicity, Score)) + geom_point() + theme_science()
ggsave("hydrophobicity_exp2.tiff")
cor.test(All_reads_e2_nostop$Score, All_reads_e2_nostop$charge)
ggplot(All_reads_e2_nostop, aes(charge, Score)) + geom_point() + theme_science()
ggsave("charge_exp2.tiff")
cor.test(All_reads_e2_nostop$Score, All_reads_e2_nostop$tAI)
ggplot(All_reads_e2_nostop, aes(tAI, Score)) + geom_point() + theme_science()
ggsave("tai_exp2.tiff")

cor.test(All_reads_e2_nostop$RFU, All_reads_e2_nostop$hydrophobicity)
ggplot(All_reads_e2_nostop, aes(hydrophobicity, RFU)) + geom_point() + theme_science()
ggsave("hydrophobicity_exp2_RFU.tiff")
cor.test(All_reads_e2_nostop$RFU, All_reads_e2_nostop$charge)
ggplot(All_reads_e2_nostop, aes(charge, RFU)) + geom_point() + theme_science()
ggsave("charge_exp2_RFU.tiff")
cor.test(All_reads_e2_nostop$RFU, All_reads_e2_nostop$tAI)
ggplot(All_reads_e2_nostop, aes(tAI, RFU)) + geom_point() + theme_science()
ggsave("tai_exp2_RFU.tiff")

All_reads_e2_nostop$G <- letterFrequency(All_reads_e2_nostop$DNA, "G")
All_reads_e2_nostop$C <- letterFrequency(All_reads_e2_nostop$DNA, "C")
All_reads_e2_nostop$GC <- (All_reads_e2_nostop$G + All_reads_e2_nostop$C)

All_reads_e2_nostop$A <- letterFrequency(All_reads_e2_nostop$DNA, "A")
All_reads_e2_nostop$U <- letterFrequency(All_reads_e2_nostop$DNA, "T")
All_reads_e2_nostop$AU <- (All_reads_e2_nostop$A + All_reads_e2_nostop$U)

ggplot(All_reads_e2_nostop, aes(as.factor(GC), Score)) + geom_boxplot(notch=TRUE) + theme_science()
ggsave("GC_exp2.tiff")
ggplot(All_reads_e2_nostop, aes(Score, colour=as.factor(GC))) + geom_step(stat="ecdf") + theme_science()
ggsave("GC_ecdf_exp2.tiff")

ggplot(All_reads_e2_nostop, aes(as.factor(GC), RFU)) + geom_boxplot(notch=TRUE) + theme_science()
ggsave("GC_exp2_RFU.tiff")
ggplot(All_reads_e2_nostop, aes(RFU, colour=as.factor(GC))) + geom_step(stat="ecdf") + theme_science()
ggsave("GC_ecdf_exp2_RFU.tiff")

ggplot(All_reads_e2_nostop, aes(as.factor(AU), Score)) + 
  geom_boxplot(notch=TRUE, fill="grey") + 
  theme_science() + 
  xlab("Number of A or U nt") + ylab("GFP Score")
ggsave("AU_exp2.tiff", height =3, width=5, units="in")

ggplot(All_reads_e2_nostop, aes(as.factor(AU), RFU)) + 
  geom_boxplot(notch=TRUE, fill="grey") +  
  theme_science() + 
  xlab("Number of A or U nt") + ylab("GFP RFU")
ggsave("AU_exp2_RFU.tiff", height =3, width=5, units="in")

ggplot(All_reads_e2_nostop, aes(Score, colour=as.factor(AU))) + 
  geom_step(stat="ecdf", size=1) + 
  scale_colour_manual(values =cbscale, name="Number of \nA or U nt") + 
  theme_science() + 
  xlab("GFP Score") + ylab("Cummulative Fraction")
ggsave("AU_ecdf_exp2.tiff", height =4, width=4, units="in")

ggplot(All_reads_e2_nostop, aes(RFU, colour=as.factor(AU))) + 
  geom_step(stat="ecdf", size=1) + 
  scale_colour_manual(values =cbscale, name="Number of \nA or U nt") + 
  theme_science() + 
  xlab("GFP RFU") + ylab("Cummulative Fraction")
ggsave("AU_ecdf_exp2_RFU.tiff", height =4, width=4, units="in")

cor.test(All_reads_e2_nostop$Score, All_reads_e2_nostop$GC)
cor.test(All_reads_e2_nostop$Score, All_reads_e2_nostop$GC, method="spearman")

ggplot(All_reads_e2_nostop, aes(Reads_Plasmid, Score)) + geom_point(alpha=0.5) + theme_science()
cor.test(All_reads_e2_nostop$Score, All_reads_e2_nostop$Reads_Plasmid)
ggsave("Plasmid_score_exp2.tiff")

ggplot(All_reads_e2_nostop, aes(Reads_Plasmid, RFU)) + geom_point(alpha=0.5) + theme_science()
cor.test(All_reads_e2_nostop$RFU, All_reads_e2_nostop$Reads_Plasmid)
ggsave("Plasmid_RFU_exp2.tiff")


high_sequences_e2 <- subset(All_reads_e2_nostop, Score > 4)
low_sequences_e2 <- subset(All_reads_e2_nostop, Score < 3)
high_e2 <- high_sequences_e2$DNA
#write.fasta(as.list(high_sequences_e2$DNA), as.list(row.names(high_sequences_e2)),  "high.fasta")
low_e2 <- low_sequences_e2$DNA
#write.fasta(as.list(low_sequences_e2$DNA), as.list(row.names(low_sequences_e2)),  "low.fasta")
#write.fasta(low, "low.fasta")
category <- c(rep(1, length(high)), rep(0, length(low))) 
MD.motifs_e2 <- findMotif(append(high_e2, low_e2),category,  max.motif=6,enriched=TRUE, start.width = 6, max.width=6)
MD.motifs2_e2 <- findMotif(append(high_e2, low_e2),category,  max.motif=6,enriched=TRUE, start.width = 4, max.width=9)
MD.motifs3_e2 <- findMotif(append(high_e2, low_e2),category,  max.motif=6,enriched=TRUE, start.width = 5, max.width=9)
### Get summary of motifs 
summaryMotif(MD.motifs_e2$motifs, MD.motifs_e2$category)
summaryMotif(MD.motifs2_e2$motifs, MD.motifs2_e2$category)
summaryMotif(MD.motifs3_e2$motifs, MD.motifs3_e2$category)
### Create table of motifs in Latex 
motifLatexTable(MD.motifs_e2, main="MD motifs")
motifLatexTable(MD.motifs2_e2, main="MD motifs")
motifLatexTable(MD.motifs3_e2, main="MD motifs")
### Create table of motifs in Html 
motifHtmlTable(MD.motifs_e2, dir="/Users/Kyle/Dropbox/TI-seq/NextSeq/motif_1_exp2.html")
motifHtmlTable(MD.motifs2_e2, dir="/Users/Kyle/Dropbox/TI-seq/NextSeq/motif_2_exp2.html")
motifHtmlTable(MD.motifs3_e2, dir="/Users/Kyle/Dropbox/TI-seq/NextSeq/motif_3_exp2.html")

All_reads_e2_nostop$motif1 <- grepl("AA[G|C|A]ATT...|.AA[G|C|A]ATT..|..AA[G|C|A]ATT.|...AA[G|C|A]ATT", All_reads_e2_nostop$Sequence)

#function to assign motif positon
m1 <- function(x) { 
  if(!grepl("AA[G|C|A]ATT...|.AA[G|C|A]ATT..|..AA[G|C|A]ATT.|...AA[G|C|A]ATT", x)) y <- "Absent"
  if(grepl("AA[G|C|A]ATT...", x)) y <- "Position 1"
  if(grepl(".AA[G|C|A]ATT..", x)) y <- "Position 2"
  if(grepl("..AA[G|C|A]ATT.", x)) y <- "Position 3"
  if(grepl("...AA[G|C|A]ATT", x)) y <- "Position 4"
  return(y)
}
#applies above function to new column
All_reads_e2_nostop$motif1_all <- sapply(All_reads_e2_nostop$Sequence, m1)

All_reads_e2_nostop$motif2 <- grepl("AA[T|G|A]TAT...|.AA[T|G|A]TAT..|..AA[T|G|A]TAT.|...AA[T|G|A]TAT", All_reads_e2_nostop$Sequence)


m2 <- function(x) { 
  if(!grepl("AA[T|G|A]TAT...|.AA[T|G|A]TAT..|..AA[T|G|A]TAT.|...AA[T|G|A]TAT", x)) y <- "Absent"
  if(grepl("AA[T|G|A]TAT...", x)) y <- "Position 1"
  if(grepl(".AA[T|G|A]TAT..", x)) y <- "Position 2"
  if(grepl("..AA[T|G|A]TAT.", x)) y <- "Position 3"
  if(grepl("...AA[T|G|A]TAT", x)) y <- "Position 4"
  return(y)
}

All_reads_e2_nostop$motif2_all <- sapply(All_reads_e2_nostop$Sequence, m2)

All_reads_e2_nostop$motif3 <- grepl("AAATT", All_reads_e2_nostop$Sequence)


m3 <- function(x) { 
  if(!grepl("AAATT", x)) y <- "Absent"
  if(grepl("AAATT....", x)) y <- "Position 1"
  if(grepl(".AAATT...", x)) y <- "Position 2"
  if(grepl("..AAATT..", x)) y <- "Position 3"
  if(grepl("...AAATT.", x)) y <- "Position 4"
  if(grepl("....AAATT", x)) y <- "Position 5"
  return(y)
}

All_reads_e2_nostop$motif3_all <- sapply(All_reads_e2_nostop$Sequence, m3)


ggplot(All_reads_e2_nostop, aes(motif1, Score)) + 
  geom_boxplot(notch=TRUE, fill="#009E73") + 
  theme_science() + 
  xlab("Contains AAVATT") + ylab("GFP Score")
ggsave("motif1_exp2.tiff", height=4, width=4, unit="in")
ggplot(All_reads_e2_nostop, aes(motif1_all, Score)) + 
  geom_boxplot(notch=TRUE, fill="#009E73") + 
  theme_science() + 
  xlab("Position of AAVATT") + ylab("GFP Score")
ggsave("motif1_all_exp2.tiff", height=4, width=5, unit="in")

ggplot(All_reads_e2_nostop, aes(motif1, RFU)) + 
  geom_boxplot(notch=TRUE, fill="#009E73") + 
  theme_science() + 
  xlab("Contains AAVATT") + ylab("GFP RFU")
ggsave("motif1_exp2_RFU.tiff", height=4, width=4, unit="in")
ggplot(All_reads_e2_nostop, aes(motif1_all, RFU)) + 
  geom_boxplot(notch=TRUE, fill="#009E73") + 
  theme_science() + 
  xlab("Position of AAVATT") + ylab("GFP RFU")
ggsave("motif1_all_exp2_RFU.tiff", height=4, width=5, unit="in")

All_reads_e2_nostop_motif1_p2 <- subset(All_reads_e2_nostop, motif1_all == "Position 2")

m1_p2_n1 <- function(x) { 
  if(grepl("AAA[G|C|A]ATT..", x)) y <- "A"
  if(grepl("TAA[G|C|A]ATT..", x)) y <- "T"
  if(grepl("CAA[G|C|A]ATT..", x)) y <- "C"
  if(grepl("GAA[G|C|A]ATT..", x)) y <- "G"
  return(y)
}

All_reads_e2_nostop_motif1_p2$motif1_p2_n1 <- sapply(All_reads_e2_nostop_motif1_p2$Sequence, m1_p2_n1)

ggplot(All_reads_e2_nostop_motif1_p2, aes(motif1_p2_n1, Score)) + 
  geom_boxplot(notch=TRUE, fill="grey") + 
  theme_science() + 
  xlab("Nucleotide at Position 1") + ylab("GFP Score")
ggsave("motif1_p2_n1_exp2.tiff", height=4, width=5, unit="in")

ggplot(All_reads_e2_nostop_motif1_p2, aes(motif1_p2_n1, RFU)) + 
  geom_boxplot(notch=TRUE, fill="grey") + 
  theme_science() + 
  xlab("Nucleotide at Position 1") + ylab("GFP RFU")
ggsave("motif1_p2_n1_exp2_RFU.tiff", height=4, width=5, unit="in")

All_reads_e2_nostop_motif1_p4 <- subset(All_reads_e2_nostop, motif1_all == "Position 4")

AA_split_m1 <- colsplit(All_reads_e2_nostop_motif1_p4$AA_character, "", c("AA_1", "AA_2", "AA_3"))

All_reads_e2_nostop_motif1_p4 <- cbind(All_reads_e2_nostop_motif1_p4, AA_split_m1)

ggplot(All_reads_e2_nostop_motif1_p4, aes(reorder(AA_1, -Score), Score)) + 
  geom_boxplot(notch=TRUE, fill="grey") + 
  theme_science() + 
  xlab("AA at Position 1") + ylab("GFP Score") 
ggsave("motif1_p4_aa1_exp2.tiff", height=4, width=5, unit="in")

ggplot(All_reads_e2_nostop_motif1_p4, aes(reorder(AA_1, -RFU), RFU)) + 
  geom_boxplot(notch=TRUE, fill="grey") + 
  theme_science() + 
  xlab("AA at Position 1") + ylab("GFP RFU") 
ggsave("motif1_p4_aa1_exp2_RFU.tiff", height=4, width=5, unit="in")

ggplot(All_reads_e2_nostop, aes(motif2, Score)) + 
  geom_boxplot(notch=TRUE, fill="#F0E442") + 
  theme_science() + 
  xlab("Contains AADTAT") + ylab("GFP Score")
ggsave("motif2_exp2.tiff", height=4, width=4, unit="in")
ggplot(All_reads_e2_nostop, aes(motif2_all, Score)) + 
  geom_boxplot(notch=TRUE, fill="#F0E442") + 
  theme_science() + 
  xlab("Position of AADTAT") + ylab("GFP Score")
ggsave("motif2_all_exp2.tiff", height=4, width=5, unit="in")

ggplot(All_reads_e2_nostop, aes(motif2, RFU)) + 
  geom_boxplot(notch=TRUE, fill="#F0E442") + 
  theme_science() + 
  xlab("Contains AADTAT") + ylab("GFP RFU")
ggsave("motif2_exp2_RFU.tiff", height=4, width=4, unit="in")
ggplot(All_reads_e2_nostop, aes(motif2_all, RFU)) + 
  geom_boxplot(notch=TRUE, fill="#F0E442") + 
  theme_science() + 
  xlab("Position of AADTAT") + ylab("GFP RFU")
ggsave("motif2_all_exp2_RFU.tiff", height=4, width=5, unit="in")

All_reads_e2_nostop_motif2_p2 <- subset(All_reads_e2_nostop, motif2_all == "Position 2")

m2_p2_n1 <- function(x) { 
  if(grepl("AAA[T|G|A]TAT..", x)) y <- "A"
  if(grepl("TAA[T|G|A]TAT..", x)) y <- "T"
  if(grepl("CAA[T|G|A]TAT..", x)) y <- "C"
  if(grepl("GAA[T|G|A]TAT..", x)) y <- "G"
  return(y)
}

All_reads_e2_nostop_motif2_p2$motif2_p2_n1 <- sapply(All_reads_e2_nostop_motif2_p2$Sequence, m2_p2_n1)

ggplot(All_reads_e2_nostop_motif2_p2, aes(motif2_p2_n1, Score)) + 
  geom_boxplot(notch=TRUE, fill="grey") + 
  theme_science() + 
  xlab("Nucleotide at Position 1") + ylab("GFP Score")
ggsave("motif2_p2_n1_exp2.tiff", height=4, width=5, unit="in")

ggplot(All_reads_e2_nostop_motif2_p2, aes(motif2_p2_n1, RFU)) + 
  geom_boxplot(notch=TRUE, fill="grey") + 
  theme_science() + 
  xlab("Nucleotide at Position 1") + ylab("GFP RFU")
ggsave("motif2_p2_n1_exp2_RFU.tiff", height=4, width=5, unit="in")

All_reads_e2_nostop_motif2_p4 <- subset(All_reads_e2_nostop, motif2_all == "Position 4")

AA_split <- colsplit(All_reads_e2_nostop_motif2_p4$AA_character, "", c("AA_1", "AA_2", "AA_3"))

All_reads_e2_nostop_motif2_p4 <- cbind(All_reads_e2_nostop_motif2_p4, AA_split)

ggplot(All_reads_e2_nostop_motif2_p4, aes(reorder(AA_1, -Score), Score)) + 
  geom_boxplot(notch=TRUE, fill="grey") + 
  theme_science() + 
  xlab("AA at Position 1") + ylab("GFP Score") 
ggsave("motif2_p4_aa1_exp2.tiff", height=4, width=5, unit="in")

ggplot(All_reads_e2_nostop_motif2_p4, aes(reorder(AA_1, -RFU), RFU)) + 
  geom_boxplot(notch=TRUE, fill="grey") + 
  theme_science() + 
  xlab("AA at Position 1") + ylab("GFP RFU") 
ggsave("motif2_p4_aa1_exp2_RFU.tiff", height=4, width=5, unit="in")


ggplot(All_reads_e2_nostop, aes(motif3, Score)) + geom_boxplot(notch=TRUE) + theme_science()
ggsave("motif3_exp2.tiff")
ggplot(All_reads_e2_nostop, aes(motif3_all, Score)) + geom_boxplot(notch=TRUE) + theme_science()
ggsave("motif3_all_exp2.tiff")

ggplot(All_reads_e2_nostop, aes(motif3, RFU)) + geom_boxplot(notch=TRUE) + theme_science()
ggsave("motif3_exp2_RFU.tiff")
ggplot(All_reads_e2_nostop, aes(motif3_all, RFU)) + geom_boxplot(notch=TRUE) + theme_science()
ggsave("motif3_all_exp2_RFU.tiff")

All_reads_e2_nostop$KNIY <- grepl("[K|N][I|Y].|.[K|N][I|Y]", All_reads_e2_nostop$AA)

ggplot(All_reads_e2_nostop, aes(KNIY, Score)) + 
  xlab("K|N-I|Y Presence") + ylab("GFP Score") +
  geom_boxplot(notch=TRUE, fill="grey") + theme_science()
ggsave("KNIY_exp2.tiff", height=4, width=4, units="in")

ggplot(All_reads_e2_nostop, aes(KNIY, RFU)) + 
  xlab("K|N-I|Y Presence") + ylab("GFP RFU") +
  geom_boxplot(notch=TRUE, fill="grey") + theme_science()
ggsave("KNIY_exp2_RFU.tiff", height=4, width=4, units="in")

KNIY <- function(x) { 
  if(!grepl("[K|N][I|Y].|.[K|N][I|Y]", x)) y <- "Absent"
  if(grepl("[K|N][I|Y].", x)) y <- "Codons 3 and 4"
  if(grepl(".[K|N][I|Y]", x)) y <- "Codons 4 and 5"
  return(y)
}

All_reads_e2_nostop$KNIY_all <- sapply(All_reads_e2_nostop$AA, KNIY)

ggplot(All_reads_e2_nostop, aes(KNIY_all, Score)) + 
  xlab("K|N-I|Y Position") + ylab("GFP Score") +
  geom_boxplot(notch=TRUE, fill="grey") + theme_science()
ggsave("KNIY_all_exp2.tiff", height=4, width=5, units="in")

ggplot(All_reads_e2_nostop, aes(KNIY_all, RFU)) + 
  xlab("K|N-I|Y Position") + ylab("GFP RFU") +
  geom_boxplot(notch=TRUE, fill="grey") + theme_science()
ggsave("KNIY_all_exp2_RFU.tiff", height=4, width=5, units="in")

t.test(All_reads_e2_nostop$Score[All_reads_e2_nostop$KNIY == "TRUE"], All_reads_e2_nostop$Score[All_reads_e2_nostop$KNIY == "FALSE"])

All_reads_e2_nostop$IYKN <- grepl("[I|Y][K|N].|.[I|Y][K|N]", All_reads_e2_nostop$AA)

ggplot(All_reads_e2_nostop, aes(IYKN, Score)) + geom_boxplot(notch=TRUE) + theme_science()
ggsave("IYKN_exp2.tiff")

ggplot(All_reads_e2_nostop, aes(IYKN, RFU)) + geom_boxplot(notch=TRUE) + theme_science()
ggsave("IYKN_exp2_RFU.tiff")

t.test(All_reads_e2_nostop$Score[All_reads_e2_nostop$IYKN == "TRUE"], All_reads_e2_nostop$Score[All_reads_e2_nostop$IYKN == "FALSE"])


high_aa <- high_sequences_e2$AA
low_aa <- low_sequences_e2$AA
write.table(high_aa, "/Users/Kyle/Dropbox/TI-seq/NextSeq/high_aa_exp2.fa", row.names = FALSE, col.names = FALSE)

high_aa_test <- read.sequences("/Users/Kyle/Dropbox/TI-seq/NextSeq/high_aa_exp2.fa")

motifs_aa <- motifModel(high_aa_test)
print(motifs_aa) 
plotPositions(motifs_aa)
png(filename="aa_frequency_high_exp2.png")
plotPositions(motifs_aa)
dev.off()

motifs_aa_table <- as.data.frame(motifs_aa@mmodel)
motifs_aa_table$aa <- row.names(motifs_aa_table)
motifs_aa_table <- motifs_aa_table[!(motifs_aa_table$aa %in% c("X", "-", "B", "X", "Z")),]
ggplot(motifs_aa_table, aes(aa, V1)) + geom_bar(stat='identity') + theme_science() + xlab("Amino Acid at Position 1") + ylab("Frequency")
ggsave("aa_frequency_high_p1_exp2.tiff")
ggplot(motifs_aa_table, aes(aa, V2)) + geom_bar(stat='identity') + theme_science() + xlab("Amino Acid at Position 2") + ylab("Frequency")
ggsave("aa_frequency_high_p2_exp2.tiff")
ggplot(motifs_aa_table, aes(aa, V3)) + geom_bar(stat='identity') + theme_science() + xlab("Amino Acid at Position 3") + ylab("Frequency")
ggsave("aa_frequency_high_p3_exp2.tiff")

write.table(low_aa, "/Users/Kyle/Dropbox/TI-seq/NextSeq/low_aa_exp2.fa", row.names = FALSE, col.names = FALSE)

low_aa_test <- read.sequences("/Users/Kyle/Dropbox/TI-seq/NextSeq/low_aa_exp2.fa")

motifs_aa_low <- motifModel(low_aa_test)
print(motifs_aa_low) 
plotPositions(motifs_aa_low)
png(filename="aa_frequency_low_exp2.png")
plotPositions(motifs_aa_low)
dev.off()

motifs_aa_low_table <- as.data.frame(motifs_aa_low@mmodel)
motifs_aa_low_table$aa <- row.names(motifs_aa_low_table)
motifs_aa_low_table <- motifs_aa_low_table[!(motifs_aa_low_table$aa %in% c("X", "-", "B", "X", "Z")),]
ggplot(motifs_aa_low_table, aes(aa, V1)) + geom_bar(stat='identity') + theme_science() + xlab("Amino Acid at Position 1") + ylab("Frequency")
ggsave("aa_frequency_low_p1_exp2.tiff")
ggplot(motifs_aa_low_table, aes(aa, V2)) + geom_bar(stat='identity') + theme_science() + xlab("Amino Acid at Position 2") + ylab("Frequency")
ggsave("aa_frequency_low_p2_exp2.tiff")
ggplot(motifs_aa_low_table, aes(aa, V3)) + geom_bar(stat='identity') + theme_science() + xlab("Amino Acid at Position 3") + ylab("Frequency")
ggsave("aa_frequency_low_p3_exp2.tiff")


All_reads_e2_nostop$aa_high <- grepl("[K|N][I|L|S][R|I|L]", All_reads_e2_nostop$AA)

ggplot(All_reads_e2_nostop, aes(aa_high, Score)) + geom_boxplot(notch=TRUE) + theme_science()
ggsave("aa_high_exp2.tiff")

ggplot(All_reads_e2_nostop, aes(aa_high, RFU)) + geom_boxplot(notch=TRUE) + theme_science()
ggsave("aa_high_exp2_RFU.tiff")

All_reads_e2_nostop$KLI <- grepl("KLI", All_reads_e2_nostop$AA)

ggplot(All_reads_e2_nostop, aes(KLI, Score)) + geom_boxplot(notch=TRUE) + theme_science()
ggsave("KLI_exp2.tiff")

ggplot(All_reads_e2_nostop, aes(KLI, RFU)) + geom_boxplot(notch=TRUE) + theme_science()
ggsave("KLI_exp2_RFU.tiff")

All_reads_e2_nostop$RRR <- grepl("RRR", All_reads_e2_nostop$AA)

ggplot(All_reads_e2_nostop, aes(RRR, Score)) + geom_boxplot(notch=TRUE) + theme_science()
ggsave("RRR_exp2.tiff")

ggplot(All_reads_e2_nostop, aes(RRR, RFU)) + geom_boxplot(notch=TRUE) + theme_science()
ggsave("RRR_exp2_RFU.tiff")


All_reads_both <- merge(All_reads, All_reads_e2, by = "Sequence")

All_reads_both$diff <- abs(All_reads_both$Score.x-All_reads_both$Score.y)

All_reads_both_subset <- subset(All_reads_both, diff < 0.3)

All_reads_both_subset$mean_score <- (All_reads_both_subset$Score.x + All_reads_both_subset$Score.y)/2

All_reads_simple <- data.frame(Score = All_reads$Score, AA = All_reads$AA_character, Total_reads = All_reads$total_bins)

All_reads_simple <- dplyr::group_by(All_reads_simple, AA)

All_reads_simple_sum <- dplyr::summarise(All_reads_simple, Range = max(Score)-min(Score), Mean = mean(Score), Total_counts = sum(Total_reads))

write.csv(All_reads_simple_sum, "Range_of_scores_by_AA_exp1.csv")

LLP <- subset(All_reads, AA_character == "LLP")

ggplot(LLP, aes(reorder(Sequence, -Score), Score)) + geom_bar(stat="identity") + theme_science() + theme(axis.text.x=element_text(angle = 45, hjust = 1))
ggsave("LLP_by_score.tiff")
ggplot(LLP, aes(Sequence, Score)) + geom_bar(stat="identity") + theme_science() + theme(axis.text.x=element_text(angle = 45, hjust = 1))
ggsave("LLP_by_sequence.tiff")

split <- strsplit(as.character(LLP$Sequence), "")

LLP <- data.frame(LLP, do.call(rbind, split))

LLP$C1 <- do.call(paste, c(LLP[c("X1", "X2", "X3")], sep = "")) 
LLP$C2 <- do.call(paste, c(LLP[c("X4", "X5", "X6")], sep = ""))
LLP$C3 <- do.call(paste, c(LLP[c("X7", "X8", "X9")], sep = ""))
ggplot(LLP, aes(C1, Score)) + geom_boxplot() + theme_science() + theme(axis.text.x=element_text(angle = 45, hjust = 1))
ggsave("LLP_C1_score.tiff")

ggplot(LLP, aes(C2, Score)) + geom_boxplot() + theme_science() + theme(axis.text.x=element_text(angle = 45, hjust = 1))
ggsave("LLP_C2_score.tiff")

ggplot(LLP, aes(C3, Score)) + geom_boxplot() + theme_science() + theme(axis.text.x=element_text(angle = 45, hjust = 1))
ggsave("LLP_C3_score.tiff")

LLP$C1_C2 <- do.call(paste, c(LLP[c("C1", "C2")], sep=""))
LLP$C2_C3 <- do.call(paste, c(LLP[c("C2", "C3")], sep=""))

ggplot(LLP, aes(C1_C2, Score)) + geom_jitter() + theme_science() + theme(axis.text.x=element_text(angle = 45, hjust = 1))
ggplot(LLP, aes(C2_C3, Score)) + geom_boxplot() + theme_science() + theme(axis.text.x=element_text(angle = 45, hjust = 1))


All_reads$codon_posibilities <- All_reads$AA_character

All_reads$codon_posibilities <- gsub("L", 6, All_reads$codon_posibilities)
All_reads$codon_posibilities <- gsub("G", 4, All_reads$codon_posibilities)
All_reads$codon_posibilities <- gsub("A", 4, All_reads$codon_posibilities)
All_reads$codon_posibilities <- gsub("P", 4, All_reads$codon_posibilities)
All_reads$codon_posibilities <- gsub("I", 3, All_reads$codon_posibilities)
All_reads$codon_posibilities <- gsub("M", 1, All_reads$codon_posibilities)
All_reads$codon_posibilities <- gsub("F", 2, All_reads$codon_posibilities)
All_reads$codon_posibilities <- gsub("Y", 2, All_reads$codon_posibilities)
All_reads$codon_posibilities <- gsub("W", 1, All_reads$codon_posibilities)
All_reads$codon_posibilities <- gsub("S", 6, All_reads$codon_posibilities)
All_reads$codon_posibilities <- gsub("T", 4, All_reads$codon_posibilities)
All_reads$codon_posibilities <- gsub("C", 2, All_reads$codon_posibilities)
All_reads$codon_posibilities <- gsub("N", 2, All_reads$codon_posibilities)
All_reads$codon_posibilities <- gsub("Q", 2, All_reads$codon_posibilities)
All_reads$codon_posibilities <- gsub("K", 2, All_reads$codon_posibilities)
All_reads$codon_posibilities <- gsub("R", 6, All_reads$codon_posibilities)
All_reads$codon_posibilities <- gsub("H", 2, All_reads$codon_posibilities)
All_reads$codon_posibilities <- gsub("D", 2, All_reads$codon_posibilities)
All_reads$codon_posibilities <- gsub("E", 2, All_reads$codon_posibilities)
All_reads$codon_posibilities <- gsub("V", 4, All_reads$codon_posibilities)


split_codon <- strsplit(All_reads$codon_posibilities, "")

All_reads <- data.frame(All_reads, do.call(rbind, split_codon))

All_reads$codon_sum <- as.numeric(as.character(All_reads$X1.1)) + as.numeric(as.character(All_reads$X2.1)) + as.numeric(as.character(All_reads$X1.1))

All_reads_simple_2 <- data.frame(Score = All_reads$Score, AA = All_reads$AA_character, Codons = All_reads$codon_sum)

All_reads_simple_2 <- dplyr::group_by(All_reads_simple_2, Codons)

All_reads_simple_2_sum <- dplyr::summarise(All_reads_simple_2, Range = max(Score)-min(Score), Mean = mean(Score))

ggplot(All_reads_simple_2_sum, aes(Codons, Range)) + geom_point() + geom_smooth(method="lm")
ggsave("codons_v_range.tiff")

All_reads_no_amber <- subset(All_reads, Amber == FALSE)

All_reads_no_amber_simple <- data.frame(Score = All_reads_no_amber$Score, AA = All_reads_no_amber$AA_character, Total_reads = All_reads_no_amber$total_bins)

All_reads_no_amber_simple <- dplyr::group_by(All_reads_no_amber_simple, AA)

All_reads_no_amber_simple_sum <- dplyr::summarise(All_reads_no_amber_simple, Range = max(Score)-min(Score), Mean = mean(Score), Total_counts = sum(Total_reads))

write.csv(All_reads_no_amber_simple_sum, "Range_of_scores_by_AA_exp1_no_amber.csv")

All_reads_no_amber_simple_sum$Stop <- grepl("[*]", All_reads_no_amber_simple_sum$AA)

All_reads_no_amber_simple_sum_stop <- subset(All_reads_no_amber_simple_sum, Stop == TRUE)

write.csv(All_reads_no_amber_simple_sum_stop, "Range_of_scores_by_AA_exp1_no_amber_contains_stop.csv")

KNstop <- subset(All_reads_no_amber, AA_character == "KN*")

ggplot(KNstop, aes(reorder(Sequence, -Score), Score)) + geom_bar(stat="identity") + theme_science() + theme(axis.text.x=element_text(angle = 45, hjust = 1))
ggsave("KNstop_by_score.tiff")
ggplot(KNstop, aes(Sequence, Score)) + geom_bar(stat="identity") + theme_science() + theme(axis.text.x=element_text(angle = 45, hjust = 1))
ggsave("KNstop_by_sequence.tiff")

split <- strsplit(as.character(KNstop$Sequence), "")

KNstop <- data.frame(KNstop, do.call(rbind, split))

KNstop$C1 <- do.call(paste, c(KNstop[c("X1", "X2", "X3")], sep = "")) 
KNstop$C2 <- do.call(paste, c(KNstop[c("X4", "X5", "X6")], sep = ""))
KNstop$C3 <- do.call(paste, c(KNstop[c("X7", "X8", "X9")], sep = ""))
ggplot(KNstop, aes(C1, Score)) + geom_boxplot() + theme_science() + theme(axis.text.x=element_text(angle = 45, hjust = 1))
ggsave("KNstop_C1_score.tiff")

ggplot(KNstop, aes(C2, Score)) + geom_boxplot() + theme_science() + theme(axis.text.x=element_text(angle = 45, hjust = 1))
ggsave("KNstop_C2_score.tiff")

ggplot(KNstop, aes(C3, Score)) + geom_boxplot() + theme_science() + theme(axis.text.x=element_text(angle = 45, hjust = 1))
ggsave("KNstop_C3_score.tiff")

KNstop$C1_C2 <- do.call(paste, c(KNstop[c("C1", "C2")], sep=""))
KNstop$C2_C3 <- do.call(paste, c(KNstop[c("C2", "C3")], sep=""))

ggplot(KNstop, aes(C1_C2, Score)) + geom_boxplot() + theme_science() + theme(axis.text.x=element_text(angle = 45, hjust = 1))
ggsave("KNstop_C1_C2_score.tiff")
ggplot(KNstop, aes(C2_C3, Score)) + geom_boxplot() + theme_science() + theme(axis.text.x=element_text(angle = 45, hjust = 1))
ggsave("KNstop_C2_C3_score.tiff")


KIstop <- subset(All_reads_no_amber, AA_character == "KI*")

ggplot(KIstop, aes(reorder(Sequence, -Score), Score)) + geom_bar(stat="identity") + theme_science() + theme(axis.text.x=element_text(angle = 45, hjust = 1))
ggsave("KIstop_by_score.tiff")
ggplot(KIstop, aes(Sequence, Score)) + geom_bar(stat="identity") + theme_science() + theme(axis.text.x=element_text(angle = 45, hjust = 1))
ggsave("KIstop_by_sequence.tiff")

split <- strsplit(as.character(KIstop$Sequence), "")

KIstop <- data.frame(KIstop, do.call(rbind, split))

KIstop$C1 <- do.call(paste, c(KIstop[c("X1", "X2", "X3")], sep = "")) 
KIstop$C2 <- do.call(paste, c(KIstop[c("X4", "X5", "X6")], sep = ""))
KIstop$C3 <- do.call(paste, c(KIstop[c("X7", "X8", "X9")], sep = ""))
ggplot(KIstop, aes(C1, Score)) + geom_boxplot() + theme_science() + theme(axis.text.x=element_text(angle = 45, hjust = 1))
ggsave("KIstop_C1_score.tiff")

ggplot(KIstop, aes(C2, Score)) + geom_boxplot() + theme_science() + theme(axis.text.x=element_text(angle = 45, hjust = 1))
ggsave("KIstop_C2_score.tiff")

ggplot(KIstop, aes(C3, Score)) + geom_boxplot() + theme_science() + theme(axis.text.x=element_text(angle = 45, hjust = 1))
ggsave("KIstop_C3_score.tiff")

KIstop$C1_C2 <- do.call(paste, c(KIstop[c("C1", "C2")], sep=""))
KIstop$C2_C3 <- do.call(paste, c(KIstop[c("C2", "C3")], sep=""))

ggplot(KIstop, aes(C1_C2, Score)) + geom_boxplot() + theme_science() + theme(axis.text.x=element_text(angle = 45, hjust = 1))
ggsave("KIstop_C1_C2_score.tiff")
ggplot(KIstop, aes(C2_C3, Score)) + geom_boxplot() + theme_science() + theme(axis.text.x=element_text(angle = 45, hjust = 1))
ggsave("KIstop_C2_C3_score.tiff")

All_reads_no_amber_simple_sum_stop$stop3 <- grepl("..[*]", All_reads_no_amber_simple_sum_stop$AA)
All_reads_no_amber_simple_sum_stop$stop2 <- grepl(".[*].", All_reads_no_amber_simple_sum_stop$AA)
All_reads_no_amber_simple_sum_stop$stop1 <- grepl("[*]..", All_reads_no_amber_simple_sum_stop$AA)

split <- strsplit(as.character(All_reads_no_amber_simple_sum_stop$AA), "")

All_reads_no_amber_simple_sum_stop <- data.frame(All_reads_no_amber_simple_sum_stop, do.call(rbind, split))

AAstop <- subset(All_reads_no_amber_simple_sum_stop, stop3 == TRUE)

AAstop <- subset(AAstop, stop1 == FALSE)
AAstop <- subset(AAstop, stop2 == FALSE)

ggplot(AAstop, aes(reorder(X1, -Mean), Mean)) + geom_boxplot() + theme_science() + theme(axis.text.x=element_text(angle = 45, hjust = 1))
ggsave("XXstop_AA1.tiff")
ggplot(AAstop, aes(reorder(X2, -Mean), Mean)) + geom_boxplot() + theme_science() + theme(axis.text.x=element_text(angle = 45, hjust = 1))
ggsave("XXstop_AA2.tiff")



All_reads_no_amber$NC3 <- grepl("......[G|T|C]TG", All_reads_no_amber$Sequence)

NC3_true <- subset(All_reads_no_amber, NC3 == TRUE)

NC3_true_simple  <- data.frame(Score = NC3_true$Score, AA = NC3_true$AA_character, Total_reads = NC3_true$total_bins)

NC3_true_simple <- dplyr::group_by(NC3_true_simple, AA)

NC3_true_simple_sum <- dplyr::summarise(NC3_true_simple, Range = max(Score)-min(Score), Mean = mean(Score), Total_counts = sum(Total_reads))

NC3_true_simple_sum$stop3 <- grepl("..[*]", NC3_true_simple_sum$AA)
NC3_true_simple_sum$stop2 <- grepl(".[*].", NC3_true_simple_sum$AA)
NC3_true_simple_sum$stop1 <- grepl("[*]..", NC3_true_simple_sum$AA)

split <- strsplit(as.character(NC3_true_simple_sum$AA), "")

NC3_true_simple_sum <- data.frame(NC3_true_simple_sum, do.call(rbind, split))

AstopA <- subset(NC3_true_simple_sum, stop3 == FALSE)

AstopA <- subset(AstopA, stop1 == FALSE)
AstopA <- subset(AstopA, stop2 == TRUE)

ggplot(AstopA, aes(reorder(X1, -Mean), Mean)) + geom_boxplot() + theme_science() + theme(axis.text.x=element_text(angle = 45, hjust = 1))
ggsave("XstopX_AA1_near_cognate.tiff")
ggplot(AstopA, aes(reorder(X3, -Mean), Mean)) + geom_boxplot() + theme_science() + theme(axis.text.x=element_text(angle = 45, hjust = 1))
ggsave("XstopX_AA3_near_cognate.tiff")


NC3_false <- subset(All_reads_no_amber, NC3 == FALSE)

NC3_false_simple  <- data.frame(Score = NC3_false$Score, AA = NC3_false$AA_character, Total_reads = NC3_false$total_bins)

NC3_false_simple <- dplyr::group_by(NC3_false_simple, AA)

NC3_false_simple_sum <- dplyr::summarise(NC3_false_simple, Range = max(Score)-min(Score), Mean = mean(Score), Total_counts = sum(Total_reads))

NC3_false_simple_sum$stop3 <- grepl("..[*]", NC3_false_simple_sum$AA)
NC3_false_simple_sum$stop2 <- grepl(".[*].", NC3_false_simple_sum$AA)
NC3_false_simple_sum$stop1 <- grepl("[*]..", NC3_false_simple_sum$AA)

split <- strsplit(as.character(NC3_false_simple_sum$AA), "")

NC3_false_simple_sum <- data.frame(NC3_false_simple_sum, do.call(rbind, split))

AstopA <- subset(NC3_false_simple_sum, stop3 == FALSE)

AstopA <- subset(AstopA, stop1 == FALSE)
AstopA <- subset(AstopA, stop2 == TRUE)

ggplot(AstopA, aes(reorder(X1, -Mean), Mean)) + geom_boxplot() + theme_science() + theme(axis.text.x=element_text(angle = 45, hjust = 1))
ggsave("XstopX_AA1.tiff")
ggplot(AstopA, aes(reorder(X3, -Mean), Mean)) + geom_boxplot() + theme_science() + theme(axis.text.x=element_text(angle = 45, hjust = 1))
ggsave("XstopX_AA3.tiff")


stopAA <- subset(NC3_true_simple_sum, stop3 == FALSE)

stopAA <- subset(stopAA, stop1 == TRUE)
stopAA <- subset(stopAA, stop2 == FALSE)

ggplot(stopAA, aes(reorder(X2, -Mean), Mean)) + geom_boxplot() + theme_science() + theme(axis.text.x=element_text(angle = 45, hjust = 1))
ggsave("stopXX_AA2_near_cognate.tiff")
ggplot(stopAA, aes(reorder(X3, -Mean), Mean)) + geom_boxplot() + theme_science() + theme(axis.text.x=element_text(angle = 45, hjust = 1))
ggsave("stopXX_AA3_near_cognate.tiff")


stopAA <- subset(NC3_false_simple_sum, stop3 == FALSE)

stopAA <- subset(stopAA, stop1 == TRUE)
stopAA <- subset(stopAA, stop2 == FALSE)

ggplot(stopAA, aes(reorder(X2, -Mean), Mean)) + geom_boxplot() + theme_science() + theme(axis.text.x=element_text(angle = 45, hjust = 1))
ggsave("stopXX_AA2.tiff")
ggplot(stopAA, aes(reorder(X3, -Mean), Mean)) + geom_boxplot() + theme_science() + theme(axis.text.x=element_text(angle = 45, hjust = 1))
ggsave("stopXX_AA3.tiff")


All_reads_e1_e2 <- merge(All_reads_e2, All_reads, by="Sequence")
ggplot(All_reads_e1_e2 , aes(Score.x, Score.y)) + xlab("GFP Score Experiment 2")+ ylab("GFP Score Experiment 1") + geom_point(alpha=0.1) + theme_science()
ggsave("e1_v_e2_final.tiff")
cor.test(All_reads_e1_e2$Score.x, All_reads_e1_e2$Score.y)
