library(dada2)
library(Cairo)

#loading the fastq files into R
file_path <- "/home/nelly-wambui/yosef-data"
head(list.files(file_path))

#sorting the data
dataF <- sort(list.files(file_path, pattern="_L001_R1_001.fastq.gz", full.names = TRUE))
head(dataF)
dataR <- sort(list.files(file_path, pattern="_L001_R2_001.fastq.gz", full.names = TRUE))
head(dataR)

#Extracting the file names
list.sample.names_F <- sapply(strsplit(basename(dataF), "_"), `[`, 1)
list.sample.names_F
list.sample.names_R <- sapply(strsplit(basename(dataR), "_"), `[`, 1)
list.sample.names_R

#visualizing the quality of the plots
plotQualityProfile(dataF[1:5])
plotQualityProfile(dataR[1:5])

library(ShortRead)
library(Biostrings)
library(stringr)

#primer trimming
fwd_primer <- "CCTACGGGNGGCWGCAG"
rev_primer <- "GACTACHVGGGTATCTAATCC"

# Function to check the orientation of these primers in the data
allOrients <- function(primer) {
  # Create all orientations of the input sequence
  require(Biostrings)
  dna <- DNAString(primer)  # The Biostrings works w/ DNAString objects rather than character vectors
  orients <- c(Forward = dna, Complement = Biostrings::complement(dna), Reverse = reverse(dna), 
               RevComp = reverseComplement(dna))
  return(sapply(orients, toString))  # Convert back to character vector
}

# get the posisible orientations of the forward and reverse primers
fwd_primer_orients <- allOrients(fwd_primer)
fwd_primer_orients
rev_primer_orients <- allOrients(rev_primer)
rev_primer_orients
fwd_primer_rev <- as.character(reverseComplement(DNAStringSet(fwd_primer))) # reverse complement of the primers
fwd_primer_rev
rev_primer_rev <- as.character(reverseComplement(DNAStringSet(rev_primer)))
rev_primer_rev

#function for counting the reads containing primers
count_primers <- function(primer, filename) {
  num_hits <- vcountPattern(primer, sread(readFastq(filename)), fixed = FALSE)
  return(sum(num_hits > 0))
}

# counting the primers on one set of paired end FASTQ files
rbind(R1_fwd_primer = sapply(fwd_primer_orients, count_primers, filename = dataF[[1]]), 
      R2_fwd_primer = sapply(fwd_primer_orients, count_primers, filename = dataF[[1]]), 
      R1_rev_primer = sapply(rev_primer_orients, count_primers, filename = dataF[[1]]), 
      R2_rev_primer = sapply(rev_primer_orients, count_primers, filename = dataF[[1]]))

rbind(R1_fwd_primer = sapply(fwd_primer_orients, count_primers, filename = dataR[[1]]), 
      R2_fwd_primer = sapply(fwd_primer_orients, count_primers, filename = dataR[[1]]), 
      R1_rev_primer = sapply(rev_primer_orients, count_primers, filename = dataR[[1]]), 
      R2_rev_primer = sapply(rev_primer_orients, count_primers, filename = dataR[[1]]))

# path to cutadapt
cutadapt <- "/opt/apps/cutadapt/3.4/bin/cutadapt" # change to the cutadapt path on your machine
system2(cutadapt, args = "--version") # confirm if the version matches that on the server

# Create an output directory called cutadapt to store the clipped files
cut_dir <- file.path(file_path, "cutadapt")
if (!dir.exists(cut_dir)) dir.create(cut_dir)

fwd_cut <- file.path(cut_dir, basename(dataF))
rev_cut <- file.path(cut_dir, basename(dataR))

names(fwd_cut) <- list.sample.names_F
head(fwd_cut)
names(rev_cut) <- list.sample.names_R
head(rev_cut)

# function for creating cutadapt trimming log files
cut_logs <- path.expand(file.path(cut_dir, paste0(list.sample.names_F, ".log")))
cut_logs <- path.expand(file.path(cut_dir, paste0(list.sample.names_R, ".log")))

# Function specifying the cutadapt functions to be used in this analysis
cutadapt_args <- c("-g", fwd_primer, "-a", rev_primer_rev, 
                   "-G", rev_primer, "-A", fwd_primer_rev,
                   "-n", 2,"-m",1, "-j",32, "--discard-untrimmed")

# creating a loop over the list of files and running cutadapt on each file.
for (i in seq_along(dataF)) {
  system2(cutadapt, 
          args = c(cutadapt_args,
                   "-o", fwd_cut[i], "-p", rev_cut[i], 
                   dataF[i], dataR[i]),
          stdout = NULL)  
}

#checking if the forward and reverse primers have been trimmed
rbind(R1_fwd_primer = sapply(fwd_primer_orients, count_primers, filename = fwd_cut[[1]]), 
      R2_fwd_primer = sapply(fwd_primer_orients, count_primers, filename = fwd_cut[[1]]), 
      R1_rev_primer = sapply(rev_primer_orients, count_primers, filename = fwd_cut[[1]]), 
      R2_rev_primer = sapply(rev_primer_orients, count_primers, filename = fwd_cut[[1]]))

rbind(R1_fwd_primer = sapply(fwd_primer_orients, count_primers, filename = rev_cut[[1]]), 
      R2_fwd_primer = sapply(fwd_primer_orients, count_primers, filename = rev_cut[[1]]), 
      R1_rev_primer = sapply(rev_primer_orients, count_primers, filename = rev_cut[[1]]), 
      R2_rev_primer = sapply(rev_primer_orients, count_primers, filename = rev_cut[[1]]))

#plotting the quality plots to determine the filtering lengths
plotQualityProfile(fwd_cut[1:5])
plotQualityProfile(rev_cut[1:5])

#assigning where the filtered data should be stored "filtered_directory_
filt.dataF <- file.path(file_path, "filtered", paste0(list.sample.names_F, "_1.fastq.gz"))
filt.dataR <- file.path(file_path, "filtered", paste0(list.sample.names_R, "_2.fastq.gz"))


names(filt.dataF) <- list.sample.names_F
head(filt.dataF)
names(filt.dataR) <- list.sample.names_R
head(filt.dataR)

#check for duplication 
duplicated(filt.dataF, filt.dataR)

#file.exists(filt.dataF)
#file.exists(filt.dataR)


#filtering and trimming data
out <- filterAndTrim(fwd_cut, filt.dataF, rev_cut, filt.dataR, truncLen=c(240,200),
                     maxN=0, maxEE=c(2,5), truncQ=2, rm.phix=TRUE,
                     compress=TRUE, multithread=TRUE)
head(out)

#checking the quality of the quality trimmed reads
plotQualityProfile(filt.dataF[1:5])
plotQualityProfile(filt.dataR[1:5])

#Establishing the error rates in the data for both the forwards and reverses
errF <- learnErrors(filt.dataF, multithread=TRUE)
errR <- learnErrors(filt.dataR, multithread=TRUE)

#plotting the error rates
plotErrors(errF, nominalQ=TRUE)
plotErrors(errR, nominalQ=TRUE)

#denoising
dadaF <- dada(filt.dataF, err=errF, multithread=TRUE)
dadaR <- dada(filt.dataR, err=errR, multithread=TRUE)

dadaF[[1]]
dadaR[[1]]

#merging the forward and reverse reads
merge.reads <- mergePairs(dadaF, filt.dataF, dadaR, filt.dataR, verbose=FALSE)
head(merge.reads[[1]])

seqtab <- makeSequenceTable(merge.reads)
dim(seqtab)

#inspecting the number of samples distributed per length
table(nchar(getSequences(seqtab)))

#Removing chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab, method="consensus", multithread=TRUE, verbose=TRUE)
dim(seqtab.nochim)

sum(seqtab.nochim)/sum(seqtab)

getN <- function(x) sum(getUniques(x))
track.nbr.reads <- cbind(out, sapply(dadaF, getN), sapply(dadaR, getN), sapply(merge.reads, getN), rowSums(seqtab.nochim),
                         final_perc_reads_retained=round(rowSums(seqtab.nochim)/out[,1]*100, 1))
track.nbr.reads
colnames(track.nbr.reads) <- c("input", "filtered", "denoisedF", "denoisedR", "merged", "nonchim", "per_retained")
colnames(track.nbr.reads)
rownames(track.nbr.reads) <- list.sample.names_F
list.sample.names
head(track.nbr.reads)

#Taxonomic classification.
taxa <- assignTaxonomy(seqtab.nochim, "/data/asatsa/dada2_preprocess/dada2_amplicon_ex_workflow/refdb/silva_nr99_v138.1_train_set.fa", multithread=TRUE)
taxa <- addSpecies(taxa, "/data/asatsa/dada2_preprocess/dada2_amplicon_ex_workflow/refdb/silva_species_assignment_v138.1.fa.gz") #classification at species level
taxa.print <- taxa # Reassigning taxa to a new name for downstream analysis

#Defining the rownames for the three table
asv_headers <- vector(dim(seqtab.nochim)[2], mode="character")

for (i in 1:dim(seqtab.nochim)[2]) {
  asv_headers[i] <- paste(">ASV", i, sep="_")
}

head(asv_headers)

#generating a sequence table having the defined row names
library(tidyverse)
seqs <- getSequences(seqtab.nochim)
asv_fasta <- c(rbind(asv_headers, seqs))
head(asv_fasta)
write(asv_fasta, "tick_ASV.fasta")

#creating a sequence table dataframe
names(seqs) <- sub(">", "", asv_headers)
seqs <- as.data.frame(seqs)
seqs <- seqs %>% rownames_to_column(var = "OTU")

#generating a feature table with newly defined row names
count_asv_tab <- t(seqtab.nochim)
row.names(count_asv_tab) <- sub(">", "", asv_headers)
write.csv(count_asv_tab, "ASVs_tick_counts.csv")

#generating a taxonomy table with the newly defined row names
rownames(taxa.print) <- gsub(pattern=">", replacement="", x=asv_headers)
head(taxa.print)
write.csv(taxa.print, file="ASVs_tick_taxonomy.csv")

#Creating a matrix of the taxonomy and feature tables for creation of a phyloseq object

library(phyloseq)
#taxonomy table
TAX = tax_table(taxa.print)

#feature table
OTU = otu_table(count_asv_tab, taxa_are_rows = TRUE)

library(dplyr)
#Reading the sample meta data into R
sample_metadata <- read.csv("/home/nelly-wambui/yosef-metadata/yosef-metadata.csv", sep = ",", header = T)
sample_metadata <- sample_metadata[,-1]
sdata <- sample_metadata %>% remove_rownames %>% column_to_rownames(var="SampleID")
samdata = sample_data(sdata)



###Reading the sample meta data into R
#sdata <- read.csv("tick_sample_metadata.csv", sep = ',', header = TRUE)
#colnames(sdata) <- c("Sample_ID", "Sex", "Infection_status", "Location")
#sdata1 <- sdata %>% remove_rownames %>% column_to_rownames(var="Sample_ID")
#samdata = sample_data(sdata)

#Yosef
#sdata <- read.csv("stingless_bee_sample_metadata.csv", sep = ',', header = TRUE)
#sdata1 <-sdata
#samdata = sample_data(sdata1)


#creating a phyloseq object
physeq = phyloseq(OTU, TAX, samdata)
physeq

##Create output files
## Create BIOM file
biomformat::write_biom(biomformat::make_biom(data = t(as.matrix(otu_table(physeq)))), 
                       "/node/cohort4/nelly/16s-rRNA-Project/stinglessbee.biom")


#filtering the unwanted sequences
physeq0 <- subset_taxa(physeq, (Order!="Chloroplast") | is.na(Order))
ntaxa(physeq0)
physeq0 <- subset_taxa(physeq0, (Phylum!="Chloroflexi") | is.na(Phylum))
ntaxa(physeq0)
physeq0 <- subset_taxa(physeq0, (Family!="Mitochondria") | is.na(Family))
ntaxa(physeq0)
physeq0 <- subset_taxa(physeq0, (Kingdom!="Archaea") | is.na(Kingdom))
ntaxa(physeq0)
physeq0 <- subset_taxa(physeq0, (Kingdom!="Eukaryota") | is.na(Kingdom))
ntaxa(physeq0)

#removing the negative control
newPhyloObject = subset_samples(physeq0, sample_names(physeq0) != "NC")

filtered_physeq <- prune_taxa(taxa_sums(newPhyloObject) > 5, newPhyloObject)
filtered_physeq
physeq1 <- subset_taxa(filtered_physeq, (is.na(Genus)))
ntaxa(physeq1)
physeq2 <- subset_taxa(filtered_physeq, (!is.na(Genus)))
ntaxa(physeq2)

 library(metagMisc)
#Extracting the filtered taxonomy and feature tables for barplot plotting
tax_table <- phyloseq_to_df(physeq1, addtax = F, addtot = F, addmaxrank = F)

library(janitor)
cumulation <- tax_table %>% adorn_totals(c("col"))
cumulation <- cumulation[order(cumulation$Total, decreasing = TRUE),]

silva_classified <- phyloseq_to_df(physeq2, addtax = F, addtot = F, addmaxrank = F)


#BLAST for all Microbiota 
#merging the silva classified taxonomy with the blast classified ones
#silva_blast_raw <- as.data.frame(bind_rows(silva_classified, merged_data))
silva_blast_raw <- silva_classified
silva_blast_raw

library(tidyverse)
library(seqRFLP)
Microbiota_seq <- merge(seqs, silva_blast_raw, by = "OTU", all = FALSE)
Microbiota_seq <- Microbiota_seq %>% select(OTU, seqs)
M_blast_abundance <- silva_blast_raw[,c(1:5)]
M_to_blast_dada2_tick_sequences <- dataframe2fas(Microbiota_seq, file = "M_to_blast_dada2_tick_sequences.fasta")

#Running blast
blastn = "/opt/apps/blast/2.10.0+/bin/blastn"
blast_db = "/home/nelly-wambui/16SMicrobial_v4/16SMicrobial"
input = "M_to_blast_dada2_tick_sequences.fasta"
evalue = 1e-6
format = 6
max_target = 1

colnames <- c("qseqid",
              "sseqid",
              "evalue",
              "bitscore",
              "sgi",
              "sacc")

blast_out <- system2("/opt/apps/blast/2.10.0+/bin/blastn", 
                     args = c("-db", blast_db, 
                              "-query", input, 
                              "-outfmt", format, 
                              "-evalue", evalue,
                              "-max_target_seqs", max_target,
                              "-ungapped"),
                     wait = TRUE,
                     stdout = TRUE) %>%
  as_tibble() %>% 
  separate(col = value, 
           into = colnames,
           sep = "\t",
           convert = TRUE)

#Removing .1-9 string to get the rightful accession numbers
blast_out$sacc <- gsub(".[.1-9]$", "", blast_out$sseqid)

#changing the collumn names
#blast_out_1 <- blast_out 
#colnames(blast_out_1) <- c("OTU", "sseqid", "evalue", "bitscore", "sgi", "sacc")


library(taxonomizr)
sacc <- as.vector(blast_out$sacc)
taxaId<-accessionToTaxa(sacc,"/home/nelly-wambui/16s-rRNA-Project/accessionTaxa.sql",version='base')
print(taxaId)
blast_taxa <- getTaxonomy(taxaId,'/home/nelly-wambui/16s-rRNA-Project/accessionTaxa.sql', rownames = FALSE)
print(blast_taxa)

blast_taxa <- as.data.frame(blast_taxa)
blast_taxa <- rownames_to_column(blast_taxa, "taxaId")
blast_out$OTU <- blast_out$qseqid


library(dplyr)

#Removing the first staxids and making the OTU collumn the first
blast_taxa <- subset(blast_taxa, select = -taxaId)
blast_taxa <-blast_taxa %>%
  select(OTU, everything()) #making OTU the first collumn in blast taxa dataframe


#blast_taxa <- read.csv("blast_taxa.csv", sep = ',', header = TRUE)
#blast_taxa <- blast_taxa[,-1]


blast_results <- merge(blast_taxa, blast_out, by = "OTU", all = FALSE)
blast_results <- merge(blast_results, Microbiota_seq, by = "OTU", all = FALSE)

#checking for the presence of duplicates
anyDuplicated(blast_results$OTU)
blast_results <- blast_results[!duplicated(blast_results$OTU),]

#changing the collumn names
colnames(blast_results) <- c("OTU", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species", "sseqid", "evalue", "bitsocre", "sgi", "sacc", "Seqs")

#merging the blast taxonomic clssification to blast abundance table
merged_data <- merge(blast_results, M_blast_abundance, by = "OTU", all = FALSE)


#grouping the data (entire dataset)
Featured_table <- merged_data[,c(7,16:19)]

group <- Featured_table %>%
  group_by(Genus)%>%
  summarise_if(is.numeric, sum)

#Groups the data in the defined order which eases downstream analysis
group <- Featured_table %>%
  group_by(Genus)%>%
  summarise_each(funs(sum), "X10K","X11K","X12K","X13K")
                 
                 #"KL1026","KL1059","KL1082","KL1087","KL1091","KL1100","KL1101","KL1102","KL1103","KL920","KL933","KL940","KL944",
                 #"KL950","KL955","KL956","KW188","KW200","KW211","KW212","KW213","KW214","KW215", "KW216","KW217","KW218","KW220","KW222","KW240",
                 #"KW285","KW323","KW330","KW345","KW395","KW400","KW487","KW489","KW490","KW527","KW566","KW573","KW578","KW610","KW64",
               #"KW650","KW651","KW657","KW700","KW711","KW740","KW748","KW805","KW806","KW821","KW845","KW847","KW850","KW98")




View (group)

#group <- read.csv("group.csv", sep = ",", header = T)
#group <- group[,-1]
#View(group)
#attach(group)

#creating multiple dataframes for the different treatments
D_schmidti <- group[,c(1,2)]
Hypotrigona_spp_1 <- group[,c(1,3:5)]
#male_uninfected <- group[,c(1,31:45)]
#female_uninfected <- group[,c(1,46:59)]

D_schmidti_total <- D_schmidti %>% adorn_totals(c("col"))
D_schmidti_total <- mutate(D_schmidti_total, D_schmidti=rowSums(D_schmidti_total[3])/1)
D_schmidti_total <- D_schmidti_total[,c(1,4)]

Hypotrigona_spp_1_total <- Hypotrigona_spp_1 %>% adorn_totals(c("col"))
Hypotrigona_spp_1_total <- mutate(Hypotrigona_spp_1_total, Hypotrigona_spp_1=rowSums(Hypotrigona_spp_1_total[5])/3)
Hypotrigona_spp_1_total <- Hypotrigona_spp_1_total[,c(1,6)]


#merging the above dataframes
merged <- Reduce(function(x,y) merge(x,y,by="Genus",all=TRUE), list(Hypotrigona_spp_1_total,D_schmidti_total))


#calculating the total abudance per genus and ordering from the most abudant to the lowest
cumulation <- merged %>% adorn_totals(c("col"))
cumulation <- cumulation[order(cumulation$Total, decreasing = TRUE),]
cumulation$perc = cumulation$Total / sum(cumulation$Total) * 100

#Original_for Figure 1_specifying the taxa to be tabulated
to_represent <- c("Bombilactobacillus","Bifidobacterium","Saccharibacter","Bombella","Acinetobacter","Pantoea","Zymobacter","Rosenbergiella","Gibbsiella","Lentilactobacillus", "Others")


#aggregating the rest of the phyla as others
grouped_data <- aggregate(merged[-1], list(Genus = replace(merged$Genus,!(merged$Genus %in% to_represent), "Others")), sum)
View(grouped_data) 

#converting the abudances into percentage
bar <- adorn_percentages(grouped_data, denominator = "col", na.rm = FALSE)
#barp <- adorn_percentages(merged, denominator = "col", na.rm = FALSE)

#2-way tabyls % - Yosef added
bar %>%
  adorn_totals("row") %>%
  adorn_pct_formatting()

#gathering the data
bar <- bar %>%
  gather(value = "abundance", key = "species", -Genus)
bar <- as.data.frame(gsub("\\(", " (", as.matrix(bar)))

# coerce the dataframe columns into respective data type
bar$Genus <- as.factor(bar$Genus)
bar$species <- as.character(bar$species)
bar$abundance <- as.numeric(bar$abundance)

#ordering the data for plotting
bar$Genus <- reorder(bar$Genus, bar$abundance)
bar$Genus <- factor(bar$Genus, levels=rev(levels(bar$Genus)))
bar$Genus <- factor(bar$Genus, levels=c("Bombilactobacillus","Bifidobacterium","Saccharibacter","Bombella","Acinetobacter","Pantoea","Zymobacter","Rosenbergiella","Gibbsiella","Lentilactobacillus", "Others"))

# Choosing the colours to use in the barplot
myPalette <- c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", "#673770", "#D3D93E", "#38333E")


#myPalette <- c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "#599861")

#myPalette <- scale_fill_manual(values = c("Lactobacillus(Firm-5) = "#65BD62", "Lactobacillus(Firm-4)" = "#F5D0FF",  "Bifidobacterium" = "#4D4D4D", "Wolbachia" = "#2171B5", "Neokomagataea" = "#00E3E3", "Saccharibacter" = "#9EFF73", "Bombella(Alpha2.2)" = "#FBAE17", "Nguyenibacter" = "#FF6345", "Zymobacter" = "#CDFFC8", "Acinetobacter" = "#8A7C64", "Alkanindiges" = "#D1A33D", "Chryseobacterium", "Gluconacetobacter", "Snodgrassella", "Gilliamella" = "#CBD588", "Frischella", "Bartonella", "Apibacter", "Commensalibacter(Alpha2.1)""Others" = "#7E7F33"))

# change the relevant label to "italic non-italic"
lbs = brk = levels(bar$Genus)
lbs[match("Bombilactobacillus", brk)] = expression(italic("Bombilactobacillus"))
lbs[match("Bifidobacterium", brk)] = expression(italic("Bifidobacterium"))
lbs[match("Saccharibacter", brk)] = expression(italic("Saccharibacter"))
lbs[match("Bombella", brk)] = expression(italic("Bombella"))
lbs[match("Acinetobacter", brk)] = expression(italic("Acinetobacter"))
lbs[match("Pantoea", brk)] = expression(italic("Pantoea"))
lbs[match("Zymobacter", brk)] = expression(italic("Zymobacter"))
lbs[match("Rosenbergiella", brk)] = expression(italic("Rosenbergiella"))
lbs[match("Gibbsiella", brk)] = expression(italic("Gibbsiella"))
lbs[match("Lentilactobacillus", brk)] = expression(italic("Lentilactobacillus"))
lbs[match("Others", brk)] = expression(plain("Others"))

#bar_all <- read.csv("bar_all.csv", sep = ",", header = T)
#bar_all <- bar_all[,-1]
#attach(bar_all)
#View(bar_all)


#plotting the barplot 
p <- ggplot(bar,aes(x = fct_inorder(species), y = abundance), labs(fill= Genus), group=row.names(bar))+ xlab("Species")+ ylab("abundance") + geom_col(aes(fill = Genus),position = position_stack(reverse = FALSE))+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic", family = "Arial"))+
  scale_fill_manual(values = myPalette,labels = lbs, breaks = brk)+
  guides(fill = guide_legend(reverse = FALSE))+
  theme(legend.text = element_text(size = 10, colour = "black", face = "italic", family = "Arial"), legend.text.align = 0)+
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, size = 10, family = "Arial"))+
  theme(axis.text = element_text(colour = "black", size = 10, family = "Arial"))+
  theme(axis.line = element_line())+
  theme(panel.background = element_rect(fill = "white"),plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), plot.background = element_rect(colour = NULL, size = 1))+
  theme(axis.ticks.length.y = unit(.15, "cm"), axis.ticks.length.x = unit(.25, "cm"), axis.text.x = element_text(margin = margin(t = .3, unit = "cm")))+
  theme(legend.position = "right", legend.justification = "top", legend.direction = "vertical", legend.text = element_text(size = 10))+
  theme(legend.key = element_rect(fill = "white"))+
  theme(legend.title = element_text(face = NULL, size = 10, family = "Arial"))+
  theme(panel.background = element_blank(), axis.text = element_blank())+
  theme(axis.text = element_text(colour = "black", size = 10, family = "Arial")+
          theme(axis.title.y = element_text(size = 10, face = "plain", angle = 90, family = "Arial"))+
          theme(axis.title.x = element_text(size = 10, angle = 0, family = "Arial")))


p + scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + xlab("Species") + theme(strip.background = element_rect(colour = "black", fill = "white")) + theme(strip.text.x = element_text(colour = "white", face = "bold")) + theme(panel.spacing = unit(1, "lines"))



#plotting the barplot #trial #coord_flip
p <- ggplot(bar,aes(x = fct_inorder(species), y = abundance), labs(fill= Genus), group=row.names(bar))+ xlab("Species")+ ylab("Abundance") + geom_col(aes(fill = Genus),position = position_stack(reverse = FALSE))+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.3, family = "Arial"))+
  scale_fill_manual(values = myPalette,labels = lbs, breaks = brk)+
  guides(fill = guide_legend(reverse = FALSE))+
  theme(legend.text = element_text(size = 10, colour = "black", face = "italic", family = "Arial"), legend.text.align = 0)+
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, size = 10, face = "italic", family = "Arial"))+
  theme(axis.text = element_text(colour = "black", size = 10, family = "Arial"))+
  theme(axis.line = element_line())+
  theme(panel.background = element_rect(fill = "white"),plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), plot.background = element_rect(colour = NULL, size = 1))+
  theme(axis.ticks.length.y = unit(.15, "cm"), axis.ticks.length.x = unit(.15, "cm"), axis.text.x = element_text(margin = margin(t = .3, unit = "cm")))+
  theme(legend.position = "right", legend.justification = "top", legend.direction = "vertical", legend.text = element_text(size = 10))+
  theme(legend.key = element_rect(fill = "white"))+
  theme(legend.title = element_text(face = NULL, size = 10, family = "Arial"))+
  theme(panel.background = element_blank(), axis.text = element_blank())+
  theme(axis.text = element_text(colour = "black", size = 10, family = "Arial")+
          theme(axis.title.y = element_text(size = 10, face = "plain", angle = 90, family = "Arial"))+
          theme(axis.title.x = element_text(size = 10, angle = 0, family = "Arial")))


p + coord_flip() + scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + xlab("Species") + theme(strip.background = element_rect(colour = "black", fill = "white")) + theme(strip.text.x = element_text(colour = "white", face = "bold")) + theme(panel.spacing = unit(1, "lines"))



#Beta diversity estimation
#Drawing the venn diagrams
#Removing all genuses with zero hits in the dataframes
D_schmidti<- D_schmidti_total[!(D_schmidti_total$D_schmidti == 0),]
Hypotrigona_spp_1<- Hypotrigona_spp_1_total[!(Hypotrigona_spp_1_total$Hypotrigona_spp_1 == 0),]



D_schmidti <- D_schmidti %>% select(Genus)
Hypotrigona_spp_1<- Hypotrigona_spp_1 %>% select(Genus)



library(gplots)
#vd <- list(M.bocandei, M.togoensis,M.ferruginea, D.schmidti, M.lendliana)
#names(vd) = c("M.bocandei", "M.togoensis","M.ferruginea", "D.schmidti", "M.lendliana")
#venn(vd) 

#Extracting sequences to be included in the study for plotting phylogenetic trees
seq <- merge(seqs, silva_blast_raw, by = "OTU", all = FALSE)
seq <- seq %>% select(OTU,seqs)

#converting the filtered sequences to fasta format ad writing the fasta file to the working directory
phylo_sequences <- dataframe2fas(seq, file = "/home/nelly-wambui/phylo_sequences.fasta")

#Reading the sequences back to R
phylo_sequences <- readDNAStringSet("/home/nelly-wambui/phylo_sequences.fasta")
names(phylo_sequences)

library(DECIPHER)
#Running multiple sequence alignment
alignment <- AlignSeqs(phylo_sequences, anchor = NA)

library(phangorn)
#constructng the phylogenetic tree
phang.align <- phyDat(as(alignment, "matrix"), type="DNA")
dm <- dist.ml(phang.align)
treeNJ <- NJ(dm) #unrooted tree
treeUPGMA <- upgma(dm) #rooted tree
fit = pml(treeNJ, data=phang.align)

#fitting the tree using the GTR model
fitGTR <- update(fit, k=4, inv=0.2)
fitGTR <- optim.pml(fitGTR, model="GTR", optInv=TRUE, optGamma=TRUE,
                    rearrangement = "stochastic", control = pml.control(trace = 0)) 

saveRDS(fitGTR, "stingless_bee_phangorn_tree.RDS")
phangorn <- readRDS("stingless_bee_phangorn_tree.RDS")

#Extracting the tree from the GTR model
phylo_tree <- phangorn$tree
phylogenetic_tree <- phy_tree(phylo_tree)
phylogenetic_tree <- root(phy_tree(phylo_tree), sample(taxa_names(phylo_tree), 1), resolve.root = TRUE)
is.rooted(phylogenetic_tree)

taxonomy_table <- silva_blast_raw[,c(1:5)] #Extracts all the taxonomic ranks from the dataframe
features <-silva_blast_raw[,c(1,2:5)] #Extracts all the abundance collumn for the different samples from the dataframe

#taxonomy table
my_taxonomy <- taxonomy_table %>% remove_rownames %>% column_to_rownames(var="OTU")
my_taxonomy <- as.matrix(my_taxonomy)
TAX2= tax_table(my_taxonomy)

#feature table
my_feature_table <- features %>% remove_rownames %>% column_to_rownames(var="OTU")
my_feature_table <- as.matrix(my_feature_table)
OTU2 = otu_table(my_feature_table, taxa_are_rows = TRUE)

fil <- readRDS("filtered_physeq.RDS")
physeq3 = merge_phyloseq(fil, phylogenetic_tree)
saveRDS(filtered_physeq, "filtered_physeq.RDS")
#creating a phyloseq object
physeq3 = merge_phyloseq(filtered_physeq, phylogenetic_tree)
physeq3

#sample_data(physeq3)[,2] <- sample_data(physeq3)[,1]

#total = median(sample_sums(physeq3))#finds median sample read count
#standf = function(x, t=total) round(t * (x / sum(x)))#function to standardize to median sample read count
#standardized_physeq = transform_sample_counts(physeq3, standf)#apply to phyloseq object
#ntaxa(standardized_physeq)
#sample_sums(standardized_physeq)

# Choosing the colours to use in the barplot
#Palette <- c('#89C5DA', "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "#599861")

#physeq3$samdata <- factor(physeq3$samdata, levels=c("M. bocandei","M. togoensis","M. ferruginea","D. schmidti","M. lendliana","H. araujoi","H. gribodoi","L. sp"))

#scale_fill_manual(values = c("M. bocandei" = "#65BD62", "M.togoensis" = "#F5D0FF",  "M.ferruginea" = "#4D4D4D", "D. schmidti" = "#2171B5", "M. lendliana" = "#00E3E3", "H.araujoi" = "#9EFF73", "L.sp" = "#FBAE17"))+
                          
 #ordinating the phyloseq object

# Choosing the colours to use in the barplot
myPalette = c("#FF0000", "#19BA24", "#0921FF", "#16A3A3", "#FF00FF", "#F5EB02", "#D36318","#120004") 

library(vegan)
#weighted unifrac
ordu = ordinate(physeq3, "PCoA", "unifrac", weighted = TRUE)
#saveRDS(ordu, "ordu.RDS")
p <- plot_ordination(physeq3, ordu, color = "Species") + geom_point(size=2)+
  scale_fill_manual(values = myPalette)+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, face = NULL, colour = "black", family = "Arial"))+
  theme(legend.text = element_text(size = 10, colour = "black", face = "italic"), legend.text.align = 0)+
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, size = 10, colour = "black", family = "Arial"))+
  theme(axis.text = element_text(colour = "black", size = 10))+
  theme(axis.line = element_line())+
  theme(legend.justification = "top")+
  theme(legend.position = "right")+
  theme(legend.key = element_rect(fill = "white"))+
  theme(legend.title = element_text(face = NULL, size = 10))+
  theme(panel.background = element_blank(), axis.text = element_blank())+
  theme(axis.text = element_text(colour = "black", size = 10, family = "Arial")+
          theme(axis.line = element_line())+
          theme(panel.background = element_rect(fill = "white"),plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), plot.background = element_rect(colour = NULL, size = NULL))+
          theme(axis.ticks.length.y = unit(.25, "cm"), axis.ticks.length.x = unit(.25, "cm"), axis.text.x = element_text(margin = margin(t = .3, unit = "cm")))+
          theme(axis.title.y = element_text(size = 10, face = "plain", angle = 90, colour = "black", family = "Arial"))+
          theme(axis.title.x = element_text(size = 10, angle = 0, face = "plain", colour = "black", family = "Arial")))
 

p + stat_ellipse(type = "norm")+labs(tag = "B", plot.tag.position = c(0.2, -0.1))


#weighted unifrac
ordu = ordinate(physeq3, "PCoA", "bray")

#weighted unifrac
ordu = ordinate(physeq3, "PCoA", "unifrac", weighted = FALSE)
p <- plot_ordination(physeq3, ordu, color="Species")+ geom_point(size=2) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
  scale_fill_manual(values = myPalette)+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, face = NULL, colour = "black", family = "Arial"))+
  theme(legend.text = element_text(size = 10, colour = "black", face = "italic"), legend.text.align = 0)+
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, size = 10, colour = "black", family = "Arial"))+
  theme(axis.text = element_text(colour = "black", size = 10))+
  theme(axis.line = element_line())+
  theme(legend.justification = "top")+
  theme(legend.position = "right")+
  theme(legend.key = element_rect(fill = "white"))+
  theme(legend.title = element_text(face = NULL, size = 10))+
  theme(panel.background = element_blank(), axis.text = element_blank())+
  theme(axis.text = element_text(colour = "black", size = 10, family = "Arial")+
          theme(axis.line = element_line())+
          theme(panel.background = element_rect(fill = "white"),plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), plot.background = element_rect(colour = NULL, size = NULL))+
          theme(axis.ticks.length.y = unit(.25, "cm"), axis.ticks.length.x = unit(.25, "cm"), axis.text.x = element_text(margin = margin(t = .3, unit = "cm")))+
          theme(axis.title.y = element_text(size = 10, face = "plain", angle = 90, colour = "black", family = "Arial"))+
          theme(axis.title.x = element_text(size = 10, angle = 0, face = "plain", colour = "black", family = "Arial")))

p + stat_ellipse(type = "norm")+labs(tag = "C", plot.tag.position = c(0.2, -0.1))

p <- plot_ordination(physeq3, ordu, color="species")+ geom_point(size=2) +
  scale_color_manual(values = my_pallet) + 
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, face = NULL, colour = "black", family = "Arial"))+
  theme(legend.text = element_text(size = 10, colour = "black", face = "italic"), legend.text.align = 0)+
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, size = 10, colour = "black", family = "Arial"))+
  theme(axis.text = element_text(colour = "black", size = 10))+
  theme(axis.line = element_line())+
  theme(legend.justification = "top")+
  theme(legend.position = "right")+
  theme(legend.key = element_rect(fill = "white"))+
  theme(legend.title = element_text(face = NULL, size = 10))+
  theme(panel.background = element_blank(), axis.text = element_blank())+
  theme(axis.text = element_text(colour = "black", size = 10, family = "Arial")+
          theme(axis.line = element_line())+
          theme(panel.background = element_rect(fill = "white"),plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), plot.background = element_rect(colour = NULL, size = NULL))+
          theme(axis.ticks.length.y = unit(.25, "cm"), axis.ticks.length.x = unit(.25, "cm"), axis.text.x = element_text(margin = margin(t = .3, unit = "cm")))+
          theme(axis.title.y = element_text(size = 10, face = "plain", angle = 90, colour = "black", family = "Arial"))+
          theme(axis.title.x = element_text(size = 10, angle = 0, face = "plain", colour = "black", family = "Arial")))

p + stat_ellipse() + labs(tag = "A", plot.tag.position = c(0.2, -0.1))
#bray curtis ordination
ordu = ordinate(physeq3, "PCoA", "bray")
p <- plot_ordination(physeq3, ordu, color="Species")+ geom_point(size=2) +
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5))+
  scale_fill_manual(values = myPalette)+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, face = NULL, colour = "black", family = "Arial"))+
  theme(legend.text = element_text(size = 10, colour = "black", face = "italic"), legend.text.align = 0)+
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, size = 10, colour = "black", family = "Arial"))+
  theme(axis.text = element_text(colour = "black", size = 10))+
  theme(axis.line = element_line())+
  theme(legend.justification = "top")+
  theme(legend.position = "right")+
  theme(legend.key = element_rect(fill = "white"))+
  theme(legend.title = element_text(face = NULL, size = 10))+
  theme(panel.background = element_blank(), axis.text = element_blank())+
  theme(axis.text = element_text(colour = "black", size = 10, family = "Arial")+
          theme(axis.line = element_line())+
          theme(panel.background = element_rect(fill = "white"),plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), plot.background = element_rect(colour = NULL, size = NULL))+
          theme(axis.ticks.length.y = unit(.25, "cm"), axis.ticks.length.x = unit(.25, "cm"), axis.text.x = element_text(margin = margin(t = .3, unit = "cm")))+
          theme(axis.title.y = element_text(size = 10, face = "plain", angle = 90, colour = "black", family = "Arial"))+
          theme(axis.title.x = element_text(size = 10, angle = 0, face = "plain", colour = "black", family = "Arial")))

p + stat_ellipse(type = "norm")+labs(tag = "A", plot.tag.position = c(0.2, -0.1))




#alpha diversity estimation
physeq4 <-physeq3 
physeq4

#checking out the total read counts in the samples
reads <- sample_sums(physeq4)
reads

summary(sample_sums(physeq4))


raremax <- 106894

library(microbiome)
#Extracting the otu table from the phyloseq object and plotting the rarefaction curve
otu_tab <- t(abundances(physeq4))

head(otu_tab)

otu_df <- as.data.frame(otu_tab)

# rarefaction curve
r=rarecurve(t(otu_table(physeq3)), step = 100, sample = raremax,xlab = "Number of reads/sample", ylab = "Number of OTUs",
            label = FALSE, xlim = c(0,100000))

Nmax <- sapply(r, function(x) max(attr(x, "Subsample")))
Smax <- sapply(r, max)
plot(c(1, max(Nmax)), c(1, max(Smax)), xlab = "Number of reads/sample",
     ylab = "Number of OTUs", type = "n")
abline(v = raremax)
p <- for (i in seq_along(r)) {
  N <- attr(r[[i]], "Subsample")
  lines(N, r[[i]], col = "black")
}

p + theme(axis.text.x = element_text(angle = 0, hjust = 0.5, size = 10, face = NULL, colour = "black", family = "Arial"))+
  theme(legend.text = element_text(size = 10, colour = "black", face = "italic"), legend.text.align = 0)+
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, size = 10, colour = "black", family = "Arial"))+
  theme(axis.line = element_line())+
  theme(panel.background = element_blank(), axis.text = element_blank())+
  theme(axis.text = element_text(colour = "black", size = 10, family = "Arial")+
          theme(axis.line = element_line())+
          theme(panel.background = element_rect(fill = "white"),plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), plot.background = element_rect(colour = NULL, size = NULL))+
          theme(axis.ticks.length.y = unit(.25, "cm"), axis.ticks.length.x = unit(.25, "cm"), axis.text.x = element_text(margin = margin(t = .3, unit = "cm")))+
          theme(axis.title.y = element_text(size = 10, face = "plain", angle = 90, colour = "black", family = "Arial"))+
          theme(axis.title.x = element_text(size = 10, angle = 0, face = "plain", colour = "black", family = "Arial")))

p + xlab("Number of reads/sample") + ylab("Number of OTUs")

#Rarefying the data
set.seed(9242)  
rarefied <- rarefy_even_depth(physeq4, sample.size = 106894)
rarefied

#calculating the alpha diversity
diversity <- alpha(rarefied, index = "all")
diversity <- rownames_to_column(diversity, "sample_id")

#Extracting the sample metadata from the phyloseq object
sdata1 <- meta(physeq4)
sdata1 <- rownames_to_column(sdata1, "sample_id")

#Extracting the shannon diversity index
shannon <- diversity %>% select(sample_id, diversity_shannon)
shannon_edited <- merge(shannon, sdata1, by = "sample_id", all = TRUE)
shannon_editted <- shannon_edited[c(28:35, 36:45, 46:55, 25:27, 1:10, 11:18, 19:20, 21:24),]


#confirming if the shannon indices are normally distributed
shapiro.test(shannon_editted$diversity_shannon)


#df.pd_editted_275 <- read.csv("df.pd_editted_275.csv", sep = ",", header = T)
#df.pd_editted_275 <- df.pd_editted_275[,-1]
#attach(df.pd_editted_275)
#View(df.pd_editted_275)

library(ggpubr)
#plotting the boxplots for the shannon index data
p <- ggboxplot(shannon_editted_275, "Species","diversity_shannon",
               color = "Species", palette = c("#F5EB02", "#D36318", "#120004", "#FF00FF", "#16A3A3", "#0921FF", "#FF0000", "#19BA24"), add = "jitter", linetype = "solid", Family = "arial", add.params = list(),
               error.plot = "pointrange", legand = NULL, size = NULL, width = 0.7, notch = FALSE, outlier.shape = 20, facet.by = NULL,
               panel.labs = NULL, short.panel.labs = TRUE,bxp.errorbar = FALSE, bxp.errorbar.width = 0.4, ggtheme = theme_pubr())+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "italic", colour = "black", family = "Arial"))+
  theme(legend.text = element_text(size = 10, colour = "black", face = "italic"), legend.text.align = 0)+
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, size = 10, colour = "black", family = "Arial"))+
  theme(axis.text = element_text(colour = "black", size = 10))+
  theme(axis.line = element_line())+
  theme(legend.justification = "top")+
  theme(legend.position = "right")+
  theme(legend.key = element_rect(fill = "white"))+
  theme(legend.title = element_text(face = NULL, size = 10))+
  theme(panel.background = element_blank(), axis.text = element_blank())+
  theme(axis.text = element_text(colour = "black", size = 10, family = "Arial")+
          theme(axis.line = element_line())+
          theme(panel.background = element_rect(fill = "white"),plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), plot.background = element_rect(colour = NULL, size = NULL))+
          theme(axis.ticks.length.y = unit(.25, "cm"), axis.ticks.length.x = unit(.25, "cm"), axis.text.x = element_text(margin = margin(t = .3, unit = "cm")))+
          theme(axis.title.y = element_text(size = 10, face = "plain", angle = 90, colour = "black", family = "Arial"))+
          theme(axis.title.x = element_text(size = 10, angle = 0, face = "plain", colour = "black", family = "Arial")))

p + aes(x = fct_inorder(Species)) + theme(legend.position = "none") + xlab("Stingless bee species") + ylab("Shannon diversity") + labs(tag = "A", plot.tag.position = c(0.2, -0.1))

#+stat_compare_means()

#old order
#c("#19BA24", "#FF0000", "#0921FF", "#16A3A3", "#FF00FF", "#F5EB02", "#D36318","#120004")

#new order
#c("#120004", "#D36318", "#F5EB02", "#FF00FF", "#16A3A3", "#0921FF", "#FF0000", "#19BA24")

#Evenness diversity estimates
evenness_pielou <- diversity %>% select(sample_id, evenness_pielou)
evenness_pielou_edited <- merge(evenness_pielou, sdata1, by = "sample_id", all = TRUE)
evenness_pielou_editted <- evenness_pielou_edited[c(28:35, 36:45, 46:55, 25:27, 1:10, 11:18, 19:20, 21:24),]

#confirming if the Evenness indices are normally distributed
shapiro.test(evenness_pielou_editted$evenness_pielou)

#plotting the boxplots for the Evenness index data
p <- ggboxplot(evenness_pielou_editted_275, "Species","evenness_pielou",
               color = "Species", palette = c("#F5EB02", "#D36318", "#120004", "#FF00FF", "#16A3A3", "#0921FF", "#FF0000", "#19BA24"), add = "jitter", linetype = "solid", Family = "arial", add.params = list(),
               error.plot = "pointrange", legand = NULL, size = NULL, width = 0.7, notch = FALSE, outlier.shape = 20, facet.by = NULL,
               panel.labs = NULL, short.panel.labs = TRUE,bxp.errorbar = FALSE, bxp.errorbar.width = 0.4, ggtheme = theme_pubr())+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "italic", colour = "black", family = "Arial"))+
  theme(legend.text = element_text(size = 10, colour = "black", face = "italic"), legend.text.align = 0)+
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, size = 10, colour = "black", family = "Arial"))+
  theme(axis.text = element_text(colour = "black", size = 10))+
  theme(axis.line = element_line())+
  theme(legend.justification = "top")+
  theme(legend.position = "right")+
  theme(legend.key = element_rect(fill = "white"))+
  theme(legend.title = element_text(face = NULL, size = 10))+
  theme(panel.background = element_blank(), axis.text = element_blank())+
  theme(axis.text = element_text(colour = "black", size = 10, family = "Arial")+
          theme(axis.line = element_line())+
          theme(panel.background = element_rect(fill = "white"),plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), plot.background = element_rect(colour = NULL, size = NULL))+
          theme(axis.ticks.length.y = unit(.25, "cm"), axis.ticks.length.x = unit(.25, "cm"), axis.text.x = element_text(margin = margin(t = .3, unit = "cm")))+
          theme(axis.title.y = element_text(size = 10, face = "plain", angle = 90, colour = "black", family = "Arial"))+
          theme(axis.title.x = element_text(size = 10, angle = 0, face = "plain", colour = "black", family = "Arial")))

p + aes(x = fct_inorder(Species)) + theme(legend.position = "none") + xlab("Stingless bee species") + ylab("Evenness diversity") + labs(tag = "B", plot.tag.position = c(0.2, -0.1))

#+stat_compare_means()

#Chao1 diversity estimates
chao1 <- diversity %>% select(sample_id, chao1)
chao_edited <- merge(chao1, sdata1, by = "sample_id", all = TRUE)
chao_editted <- chao_edited[c(28:35, 36:45, 46:55, 25:27, 1:10, 11:18, 19:20, 21:24),]

#confirming if the chao1 indices are normally distributed
shapiro.test(chao_editted$chao1)


#We dropped chao1 from Stingless bee paper #plotting the boxplots for the chao1 index data
p <- ggboxplot(chao_editted, "species","chao1",
               color = "species", palette = c("#19BA24", "#FF0000", "#0921FF", "#16A3A3", "#FF00FF", "#F5EB02", "#D36318","#120004"), add = "jitter", linetype = "solid", Family = "arial", add.params = list(),
               error.plot = "pointrange", legand = NULL, size = NULL, width = 0.7, notch = FALSE, outlier.shape = 20, facet.by = NULL,
               panel.labs = NULL, short.panel.labs = TRUE,bxp.errorbar = FALSE, bxp.errorbar.width = 0.4, ggtheme = theme_pubr())+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "italic", colour = "black", family = "Arial"))+
  theme(legend.text = element_text(size = 10, colour = "black", face = "italic"), legend.text.align = 0)+
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, size = 10, colour = "black", family = "Arial"))+
  theme(axis.text = element_text(colour = "black", size = 10))+
  theme(axis.line = element_line())+
  theme(legend.justification = "top")+
  theme(legend.position = "right")+
  theme(legend.key = element_rect(fill = "white"))+
  theme(legend.title = element_text(face = NULL, size = 10))+
  theme(panel.background = element_blank(), axis.text = element_blank())+
  theme(axis.text = element_text(colour = "black", size = 10, family = "Arial")+
          theme(axis.line = element_line())+
          theme(panel.background = element_rect(fill = "white"),plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), plot.background = element_rect(colour = NULL, size = NULL))+
          theme(axis.ticks.length.y = unit(.25, "cm"), axis.ticks.length.x = unit(.25, "cm"), axis.text.x = element_text(margin = margin(t = .3, unit = "cm")))+
          theme(axis.title.y = element_text(size = 10, face = "plain", angle = 90, colour = "black", family = "Arial"))+
          theme(axis.title.x = element_text(size = 10, angle = 0, face = "plain", colour = "black", family = "Arial")))

p + aes(x = fct_inorder(species)) + theme(legend.position = "none") + xlab("Stingless bee species") + ylab("chao1 diversity") + labs(tag = "C", plot.tag.position = c(0.2, -0.1))

#+stat_compare_means()


#Faith's phylogeetic diversity
library(picante)
FD_OTU <- as.data.frame(rarefied@otu_table)

df.pd <- pd(t(FD_OTU), phylogenetic_tree,include.root=F)
df.pd <- rownames_to_column(df.pd, "sample_id")
df.pd_edited <- merge(df.pd, sdata1, by = "sample_id", all = TRUE)
df.pd_editted <- df.pd_edited[c(28:35, 36:45, 46:55, 25:27, 1:10, 11:18, 19:20, 21:24),]

#confirming if the FD indices are normally distributed
shapiro.test(df.pd_editted$PD)

# Choosing the colours to use in the barplot
#myPalette <- c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "#599861")


#plotting the boxplots for the PD index data
p <- ggboxplot(df.pd_editted_275, "Species","PD",
               color = "Species", palette = c("#F5EB02", "#D36318", "#120004", "#FF00FF", "#16A3A3", "#0921FF", "#FF0000", "#19BA24"), add = "jitter", linetype = "solid", add.params = list(),
               error.plot = "pointrange", legand = NULL, size = NULL, width = 0.7, notch = FALSE, outlier.shape = 20, facet.by = NULL,
               panel.labs = NULL, short.panel.labs = TRUE,bxp.errorbar = FALSE, bxp.errorbar.width = 0.4, ggtheme = theme_pubr())+
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "italic", colour = "black", family = "Arial"))+
  theme(legend.text = element_text(size = 10, colour = "black", face = "italic"), legend.text.align = 0)+
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, size = 10, colour = "black", family = "Arial"))+
  theme(axis.text = element_text(colour = "black", size = 10))+
  theme(axis.line = element_line())+
  theme(legend.justification = "top")+
  theme(legend.position = "right")+
  theme(legend.key = element_rect(fill = "white"))+
  theme(legend.title = element_text(face = NULL, size = 10))+
  theme(panel.background = element_blank(), axis.text = element_blank())+
  theme(axis.text = element_text(colour = "black", size = 10, family = "Arial")+
  theme(axis.line = element_line())+
  theme(panel.background = element_rect(fill = "white"),plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), plot.background = element_rect(colour = NULL, size = NULL))+
  theme(axis.ticks.length.y = unit(.25, "cm"), axis.ticks.length.x = unit(.25, "cm"), axis.text.x = element_text(margin = margin(t = .3, unit = "cm")))+
  theme(axis.title.y = element_text(size = 10, face = "plain", angle = 90, colour = "black", family = "Arial"))+
  theme(axis.title.x = element_text(size = 10, angle = 0, face = "plain", colour = "black", family = "Arial")))

p + aes(x = fct_inorder(Species)) + theme(legend.position = "none") + xlab("Stingless bee species") + ylab("Faith's PD") + labs(tag = "C", plot.tag.position = c(0.2, -0.1))

#+stat_compare_means()
#Kruskall walis p=0.00073

#Permanova
physeq5 <- physeq3
physeq5

#Calculate bray curtis distance matrix
physeq_bray <- phyloseq::distance(physeq5, method = "bray")

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(physeq5))

# Adonis test
permanova <- adonis(physeq_bray ~ species, data = sampledf)
permanova
# P-value
print(as.data.frame(permanova$aov.tab)["species", "Pr(>F)"])

source("parwise.adonis.r")
# convert OTU ID and sample id to rownames
permanova_feature <- group %>% column_to_rownames(var= "Genus")
#sample_metadata <- sample_metadata %>% column_to_rownames(var= "sample.id")

# transpose the feature table
permanova_feature_transposed <- t(permanova_feature)

# create sample.ids column from rownames
permanova_feature_transposed <- permanova_feature_transposed %>% as.data.frame() %>% rownames_to_column(var = "sample_id")

# merge the sample metadata with the feature table using the sample ids
permanova_feature_transposed <- merge(sdata1, permanova_feature_transposed, by = "sample_id")

# perform permanova analysis 
permanova <- pairwise.adonis(permanova_feature_transposed[,!names(permanova_feature_transposed) %in% c("sample_id", "species","NA.", "NA..1", "NA..2", "NA..3", "NA..4")],
                             factors = permanova_feature_transposed$species,
                             sim.method = 'bray',
                             p.adjust.m ='bonferroni')
permanova

#install.packages("agricolae")
#install.packages("dunn.test")
library(agricolae)
library(dunn.test) 

#Shapiro-Wilk normality test
shapiro.test(chao_editted)

#Bartlett test of homogeneity of variances
bartlett.test(evenness_pielou, Species)
attach(diversity)
attach(chao_editted)
View(chao_editted)



kruskal.test(chao1~NA..3)
dunn.test(chao_editted,NA..3)

#sdata11

sdata11 <- read.csv("sdata11.csv", sep = ",", header = T)
sdata11 <- sdata11[,-1]
attach(sdata11)
View(sdata11)



#### HIVE - PERMANOVA 

####Permanova for different hive - YOSEF draft
physeq6 <- physeq3
physeq6

#Calculate bray curtis distance matrix
physeq_bray <- phyloseq::distance(physeq6, method = "bray")

# make a data frame from the sample_data
sampledf <- data.frame(sample_data(physeq6))

# Adonis test
permanova <- adonis(physeq_bray ~ NA..3, data = sampledf)
permanova

# 

P-value
print(as.data.frame(permanova$aov.tab)["NA..3", "Pr(>F)"])

source("parwise.adonis.r")
# convert OTU ID and sample id to rownames
permanova_feature <- group %>% column_to_rownames(var= "Genus")
#sample_metadata <- sample_metadata %>% column_to_rownames(var= "sample.id")

# transpose the feature table
permanova_feature_transposed <- t(permanova_feature)

# create sample.ids column from rownames
permanova_feature_transposed <- permanova_feature_transposed %>% as.data.frame() %>% rownames_to_column(var = "sample_id")

# merge the sample metadata with the feature table using the sample ids
permanova_feature_transposed <- merge(sdata11, permanova_feature_transposed, by = "sample_id")

# perform permanova analysis 
permanova <- pairwise.adonis(permanova_feature_transposed[,!names(permanova_feature_transposed) %in% c("sample_id", "code","bsize", "hsize", "hive", "location")],
                             factors = permanova_feature_transposed$hive,
                             sim.method = 'bray',
                             p.adjust.m ='bonferroni')
permanova

####HIVE_Relative abundance 

group_hive <- read.csv("group_hive.csv", sep = ",", header = T)
group_hive <- group_hive[,-1]
attach(group_hive)
View(group_hive)


#creating multiple dataframes for the different species
icipe_2H   <- group_hive[,c(1,2:15)]
icipe_4M <- group_hive[,c(1,20:23)]
icipe_5M  <- group_hive[,c(1,19,24:46)]
icipe_clay_pot <- group_hive[,c(1,16:18)]
Natural_nest1 <- group_hive[,c(1,48:52)]
Natural_nest2 <- group_hive[,c(1,47,53:56)]

icipe_2H_total <- icipe_2H %>% adorn_totals(c("col"))
icipe_2H_total <- mutate(icipe_2H_total, icipe_2H=rowSums(icipe_2H_total[16])/14)
icipe_2H_total <- icipe_2H_total[,c(1,17)]

icipe_4M_total <- icipe_4M %>% adorn_totals(c("col"))
icipe_4M_total <- mutate(icipe_4M_total, icipe_4M=rowSums(icipe_4M_total[6])/4)
icipe_4M_total <- icipe_4M_total[,c(1,7)]

icipe_5M_total <- icipe_5M %>% adorn_totals(c("col"))
icipe_5M_total <- mutate(icipe_5M_total, icipe_5M=rowSums(icipe_5M_total[26])/24)
icipe_5M_total <- icipe_5M_total[,c(1,27)]

icipe_clay_pot_total <- icipe_clay_pot %>% adorn_totals(c("col"))
icipe_clay_pot_total <- mutate(icipe_clay_pot_total, icipe_clay_pot=rowSums(icipe_clay_pot_total[5])/3)
icipe_clay_pot_total <- icipe_clay_pot_total[,c(1,6)]

Natural_nest1_total <- Natural_nest1 %>% adorn_totals(c("col"))
Natural_nest1_total <- mutate(Natural_nest1_total, Natural_nest1=rowSums(Natural_nest1_total[7])/5)
Natural_nest1_total <- Natural_nest1_total[,c(1,8)]

Natural_nest2_total <- Natural_nest2 %>% adorn_totals(c("col"))
Natural_nest2_total <- mutate(Natural_nest2_total, Natural_nest2=rowSums(Natural_nest2_total[7])/5)
Natural_nest2_total <- Natural_nest2_total[,c(1,8)]


#merging the above dataframes
merged_hive <- Reduce(function(x,y) merge(x,y,by="Genus",all=TRUE),list(icipe_2H_total,icipe_4M_total,icipe_5M_total,icipe_clay_pot_total,Natural_nest1_total,Natural_nest2_total))

#calculating the total abudance per genus and ordering from the most abudant to the lowest
cumulation_hive <- merged_hive %>% adorn_totals(c("col"))
cumulation_hive <- cumulation_hive[order(cumulation_hive$Total, decreasing = TRUE),]

#specifying the taxa to be tabulated
to_represent <- c("Bifidobacterium", "Saccharibacter", "Neokomagataea", "Wolbachia", "Bombella", "Nguyenibacter", "Zymobacter", "Acinetobacter", "Alkanindiges", "Chryseobacterium", "Gluconacetobacter", "Snodgrassella", "Gilliamella", "Bartonella", "Frischella", "Commensalibacter", "Apibacter", "Others")

#aggregating the rest of the phyla as others
grouped_data_m <- aggregate(merged_m[-1], list(Genus = replace(merged_m$Genus,!(merged_m$Genus %in% to_represent), "Others")), sum)
View(grouped_data_m) 


###barplot_210705

library(dplyr)
library(janitor)
library(ggplot2)
library(tidyverse)


#grouping the data (entire dataset)
#Featured_table_all <- merged_data_all_1[,c(7,15:69)]

#group_11 <- group_2 %>%
#group_by(Genus)%>%
#summarise_if(is.numeric, sum)

#Groups the data in the defined order which eases downstream analysis
group_y <- group_2 %>%
  mutate(Genus=str_remove(Genus," ")) %>% 
 group_by(Genus)%>%
 summarise_each(funs(sum), "Hs1_1","Hs1_2","Hs1_3","Hs1_4","Hs1_5","Hs1_6","Hs1_7","Hs1_8","Hs2_1","Hs2_2",
                "Ls_1","Ls_2","Ls_3","Ls_4","Ds_1","Ds_2","Ds_3","Ds_4","Ds_5","Ds_6","Ds_7","Ds_8","Ds_9","Ds_10","Ml_1","Ml_2","Ml_3","Mf_1","Mf_2",
                "Mf_3","Mf_4","Mf_5","Mf_6","Mf_7","Mf_8","Mf_9","Mf_10","Mt_1" ,"Mt_2",
                "Mt_3","Mt_4","Mt_5","Mt_6","Mt_7","Mt_8","Mt_9","Mt_10","Mb_1","Mb_2",
                "Mb_3","Mb_4","Mb_5","Mb_6","Mb_7","Mb_8") 
                
                
View (group_y)

#group_2 <- read.csv("group_2.csv", sep = ",", header = T)
#group_2 <- group_2[,-1]
#View(group_2)
#attach(group_2)

#creating multiple dataframes for the different species
Hypotrigona_sp1 <- group_y[,c(1:9)]
Hypotrigona_sp2 <- group_y[,c(1,10,11)]
Liotrigona_sp <- group_y[,c(1,12:15)]
D.schmidti <- group_y[,c(1,16:25)]
M.lendliana <- group_y[,c(1,26:28)]
M.ferruginea <- group_y[,c(1,29:38)]
M.togoensis <- group_y[,c(1,39:48)]
M.bocandei <- group_y[,c(1,49:56)]

Hypotrigona_sp1_total <- Hypotrigona_sp1 %>% adorn_totals(c("col"))
Hypotrigona_sp1_total <- mutate(Hypotrigona_sp1_total, Hypotrigona_sp1=rowSums(Hypotrigona_sp1_total[10])/8)
Hypotrigona_sp1_total <- Hypotrigona_sp1_total[,c(1,11)]

Hypotrigona_sp2_total <- Hypotrigona_sp2 %>% adorn_totals(c("col"))
Hypotrigona_sp2_total <- mutate(Hypotrigona_sp2_total, Hypotrigona_sp2=rowSums(Hypotrigona_sp2_total[4])/2)
Hypotrigona_sp2_total <- Hypotrigona_sp2_total[,c(1,5)]

Liotrigona_sp_total <- Liotrigona_sp %>% adorn_totals(c("col"))
Liotrigona_sp_total <- mutate(Liotrigona_sp_total, Liotrigona_sp=rowSums(Liotrigona_sp_total[6])/4)
Liotrigona_sp_total <- Liotrigona_sp_total[,c(1,7)]

D.schmidti_total <- D.schmidti %>% adorn_totals(c("col"))
D.schmidti_total <- mutate(D.schmidti_total, D.schmidti=rowSums(D.schmidti_total[12])/10)
D.schmidti_total <- D.schmidti_total[,c(1,13)]

M.lendliana_total <- M.lendliana %>% adorn_totals(c("col"))
M.lendliana_total <- mutate(M.lendliana_total, M.lendliana=rowSums(M.lendliana_total[5])/3)
M.lendliana_total <- M.lendliana_total[,c(1,6)]

M.ferruginea_total <- M.ferruginea %>% adorn_totals(c("col"))
M.ferruginea_total <- mutate(M.ferruginea_total, M.ferruginea=rowSums(M.ferruginea_total[12])/10)
M.ferruginea_total <- M.ferruginea_total[,c(1,13)]

M.togoensis_total <- M.togoensis %>% adorn_totals(c("col"))
M.togoensis_total <- mutate(M.togoensis_total, M.togoensis=rowSums(M.togoensis_total[12])/10)
M.togoensis_total <- M.togoensis_total[,c(1,13)]

M.bocandei_total <- M.bocandei %>% adorn_totals(c("col"))
M.bocandei_total <- mutate(M.bocandei_total, M.bocandei=rowSums(M.bocandei_total[10])/8)
M.bocandei_total <- M.bocandei_total[,c(1,11)]

#df <- list(Hypotrigona_sp1_total,Hypotrigona_sp2_total,Liotrigona_sp_total,D.schmidti_total,M.lendliana_total,M.ferruginea_total,M.togoensis_total,M.bocandei_total) %>% 
  #reduce(full_join, by = "Genus")

#merging the above dataframes
merged_y <- Reduce(function(x,y) merge(x,y,by="Genus",all=TRUE), list(Hypotrigona_sp1_total,Hypotrigona_sp2_total,Liotrigona_sp_total,D.schmidti_total,M.lendliana_total,M.ferruginea_total,M.togoensis_total,M.bocandei_total))


#calculating the total abudance per genus and ordering from the most abudant to the lowest
cumulation <- merged_y %>% adorn_totals(c("col"))
cumulation <- cumulation[order(cumulation$Total, decreasing = TRUE),]
cumulation$perc = cumulation$Total / sum(cumulation$Total) * 100

#Original_for Figure 1_specifying the taxa to be tabulated
#Original_for Figure 1_specifying the taxa to be tabulated
to_represent <- c("Acetilactobacillus", "Lactobacillus", "Bifidobacterium", "Bombella", "Saccharibacter", "Bombilactobacillus", "Acetobacter", "Nguyenibacter", "Others")

#to_represent <- c("Acetilactobacillus", "Lactobacillus", "Bifidobacterium", "Bombella", "Saccharibacter", "Wolbachia", "Bombilactobacillus", "Acetobacter", "Nguyenibacter", "Others")

#to_represent <- c("Lactobacillus(others)", "Lactobacillus(Firm-4)", "Bifidobacterium", "Lactobacillus(Firm-5)", "Saccharibacter", "Neokomagataea", "Wolbachia", "Bombella(Alpha2.2)", "Pediococcus", "Ameyamaea","Nguyenibacter", "Zymobacter", "Acinetobacter", "Alkanindiges", "Klebsiella", "Acetobacter", "Snodgrassella", "Gilliamella", "Bartonella", "Frischella", "Commensalibacter(Alpha2.1)", "Apibacter", "Others")

#to_represent <- c("Lactobacillus (others)", "Lactobacillus (Firm-4)", "Bifidobacterium", "Lactobacillus (Firm-5)", "Saccharibacter", "Neokomagataea", "Wolbachia", "Bombella (Alpha2.2)", "Pediococcus", "Ameyamaea","Nguyenibacter", "Zymobacter", "Acinetobacter", "Alkanindiges", "Klebsiella", "Acetobacter", "Snodgrassella", "Gilliamella", "Bartonella", "Frischella", "Commensalibacter (Alpha2.1)", "Apibacter", "Others")

#aggregating the rest of the phyla as others
grouped_data <- aggregate(merged_y[-1], list(Genus = replace(merged_y$Genus,!(merged_y$Genus %in% to_represent), "Others")), sum)
View(grouped_data) 

#converting the abudances into percentage
bar <- adorn_percentages(grouped_data, denominator = "col", na.rm = FALSE)
#barp <- adorn_percentages(merged, denominator = "col", na.rm = FALSE)

#2-way tabyls % - Yosef added
bar %>%
  adorn_totals("row") %>%
  adorn_pct_formatting()

#gathering the data
bar <- bar %>%
  gather(value = "abundance", key = "sample_names", -Genus)
bar <- as.data.frame(gsub("\\(", " (", as.matrix(bar)))

# coerce the dataframe columns into respective data type
bar$Genus <- as.factor(bar$Genus)
bar$sample_names <- as.character(bar$sample_names)
bar$abundance <- as.numeric(bar$abundance)

#ordering the data for plotting
bar$Genus <- reorder(bar$Genus, bar$abundance)
bar$Genus <- factor(bar$Genus, levels=rev(levels(bar$Genus)))
bar$Genus <- factor(bar$Genus, levels=c("Acetilactobacillus", "Lactobacillus", "Bifidobacterium", "Bombella", "Saccharibacter", "Bombilactobacillus", "Acetobacter", "Nguyenibacter", "Others"))

# Choosing the colours to use in the barplot
myPalette <- c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#CBD588", "#5F7FC7", "#BFA19C", "#599861")

#myPalette <- c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", "#599861")

#myPalette <- c("#89C5DA", "#DA5724", "#74D944", "#CE50CA", "#3F4921", "#C0717C", "#CBD588", "#5F7FC7", "#673770", "#D3D93E", "#38333E", "#508578", "#D7C1B1", "#689030", "#AD6F3B", "#CD9BCD", "#D14285", "#6DDE88", "#652926", "#7FDCC0", "#C84248", "#8569D5", "#5E738F", "#D1A33D", "#8A7C64", "#599861")

#myPalette <- scale_fill_manual(values = c("Lactobacillus(Firm-5) = "#65BD62", "Lactobacillus(Firm-4)" = "#F5D0FF",  "Bifidobacterium" = "#4D4D4D", "Wolbachia" = "#2171B5", "Neokomagataea" = "#00E3E3", "Saccharibacter" = "#9EFF73", "Bombella(Alpha2.2)" = "#FBAE17", "Nguyenibacter" = "#FF6345", "Zymobacter" = "#CDFFC8", "Acinetobacter" = "#8A7C64", "Alkanindiges" = "#D1A33D", "Chryseobacterium", "Gluconacetobacter", "Snodgrassella", "Gilliamella" = "#CBD588", "Frischella", "Bartonella", "Apibacter", "Commensalibacter(Alpha2.1)""Others" = "#7E7F33"))

# change the relevant label to "italic non-italic"
lbs = brk = levels(bar$Genus)
lbs[match("Acetilactobacillus", brk)] = expression(italic("Acetilactobacillus"))
lbs[match("Lactobacillus", brk)] = expression(italic("Lactobacillus"))
lbs[match("Bifidobacterium", brk)] = expression(italic("Bifidobacterium"))
lbs[match("Bombella", brk)] = expression(italic("Bombella"))
lbs[match("Saccharibacter", brk)] = expression(italic("Saccharibacter"))
lbs[match("Bombilactobacillus", brk)] = expression(italic("Bombilactobacillus"))
lbs[match("Acetobacter", brk)] = expression(italic("Acetobacter"))
lbs[match("Others", brk)] = expression(plain("Others"))

#plotting the barplot #trial
p <- ggplot(bar,aes(x = fct_inorder(sample_names), y = abundance), labs(fill= Genus), group=row.names(bar))+ xlab("Species")+ ylab("Abundance") + geom_col(aes(fill = Genus),position = position_stack(reverse = FALSE))+
  theme(axis.text.x = element_text(angle = 0, hjust = 0.3, family = "Arial"))+
  scale_fill_manual(values = myPalette,labels = lbs, breaks = brk)+
  guides(fill = guide_legend(reverse = FALSE))+
  theme(legend.text = element_text(size = 10, colour = "black", face = "italic", family = "Arial"), legend.text.align = 0)+
  theme(axis.text.y = element_text(angle = 0, vjust = 0.5, size = 10, face = "italic", family = "Arial"))+
  theme(axis.text = element_text(colour = "black", size = 10, family = "Arial"))+
  theme(axis.line = element_line())+
  theme(panel.background = element_rect(fill = "white"),plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm"), plot.background = element_rect(colour = NULL, size = 1))+
  theme(axis.ticks.length.y = unit(.15, "cm"), axis.ticks.length.x = unit(.15, "cm"), axis.text.x = element_text(margin = margin(t = .3, unit = "cm")))+
  theme(legend.position = "right", legend.justification = "top", legend.direction = "vertical", legend.text = element_text(size = 10))+
  theme(legend.key = element_rect(fill = "white"))+
  theme(legend.title = element_text(face = NULL, size = 10, family = "Arial"))+
  theme(panel.background = element_blank(), axis.text = element_blank())+
  theme(axis.text = element_text(colour = "black", size = 10, family = "Arial")+
          theme(axis.title.y = element_text(size = 10, face = "plain", angle = 90, family = "Arial"))+
          theme(axis.title.x = element_text(size = 10, angle = 0, family = "Arial")))


p + coord_flip() + scale_y_continuous(labels = scales::percent_format(accuracy = 1)) + xlab("Stingless bee species") + theme(strip.background = element_rect(colour = "black", fill = "white")) + theme(strip.text.x = element_text(colour = "white", face = "bold")) + theme(panel.spacing = unit(1, "lines"))



