# Author: Marie-Madlen Pust
# Last update: 10 February 2022
# Data analysis (Longitudinal CF airway metagenome)

# clean R environment
rm(list = ls())

# define global variables
# function for installing/importing packages
ipak <- function(pkg){
  new.pkg <- pkg[!(pkg %in% installed.packages()[, 'Package'])]
  if (length(new.pkg)) 
    install.packages(new.pkg, dependencies = TRUE)
  sapply(pkg, require, character.only = TRUE)}

# store required Rpackages
packages_r <- c('readr','vegan','ggplot2','scales','dplyr','ggpubr','rcompanion','ggdendro','pheatmap',
                'factoextra','matrixStats','ggrepel','tidyr','plyr','purrr','stringr','RVAideMemoire', 'umap', 'usedist')

# load required R packages
ipak(packages_r)

# import datasets
# import categorical growth file
growth_file_00 <- read_delim("python_generated_files/growth_file_categorical.csv", ";", escape_double = FALSE, trim_ws = TRUE)
growth_file_01 <- data.frame(growth_file_00)
colnames(growth_file_01)
rownames(growth_file_01) <- growth_file_01$Species
growth_file_01$Species <- NULL
# keep only CF data
MCF_growthrate <- growth_file_01[,53:111]
MCF_growthrate_t <- t(data.frame(MCF_growthrate))
MCF_growthrate_t <- data.frame(MCF_growthrate_t)

# import metadata
metadata_0 <- read_delim("python_generated_files/metadata_0.csv", ";", escape_double = FALSE, trim_ws = TRUE)
metadata_0 <- data.frame(metadata_0)
MCF_metadata <- metadata_0[1:59,]

colnames(metadata_0) <- c("Name", "Patient", "Sample_No", "Sample_taken", "Birth_date", "Height_cm", "Weight_kg", 
                          "Mutation", "Mutation_class", "P_state", "Gender", "Age_in_weeks")
metadata_1 <- metadata_0[order(metadata_0$Name),]

# import phylum table
count_df_phylum_clr_00 <- read_delim("python_generated_files/count_df_phylum_clr_4.csv", ";", escape_double = FALSE, trim_ws = TRUE)
count_df_phylum_clr_00 <- data.frame(count_df_phylum_clr_00)
rownames(count_df_phylum_clr_00) <- count_df_phylum_clr_00$X1
count_df_phylum_clr_00$X1 <- NULL
count_df_phylum_clr_00_t <- data.frame(t(count_df_phylum_clr_00))
count_df_phylum_clr_00_t <- count_df_phylum_clr_00_t[order(rownames(count_df_phylum_clr_00_t)),]
count_df_phylum_clr_00_t$State <- metadata_1$State
count_df_phylum_clr_00_t$patient <- metadata_1$Patient
count_df_phylum_clr_00_t$sample_no <- metadata_1$Sample_No
count_df_phylum_clr_00_t$Identifier <- metadata_1$Name
count_df_phylum_clr_00_t_long <- gather(count_df_phylum_clr_00_t, key="Phylum", "CLR", 1:8)

# import species table
count_clr <- count_df_sp_clr_4 <- read_delim("python_generated_files/count_df_sp_clr_4.csv", ";", escape_double = FALSE, trim_ws = TRUE)
count_clr <- data.frame(count_clr)
rownames(count_clr) <- count_clr$X1
count_clr$X1 <- NULL
colnames(count_clr) <- str_replace_all(colnames(count_clr), "\\.", "-")

# import quantitative growth table
growth_file_quant <- read_delim("python_generated_files/growth_file_numerical.csv", ";", escape_double = FALSE, trim_ws = TRUE)
growth_file_quant <- data.frame(growth_file_quant)

# data clean up
remove_from_count_file <- setdiff(rownames(count_clr), growth_file_quant$Species)
remove_from_growth_file <- setdiff(growth_file_quant$Species, rownames(count_clr))
remove_from_count_file_sample <- setdiff(colnames(count_clr), growth_file_quant$Name)
remove_from_growth_file_sample <- setdiff(growth_file_quant$Name, colnames(count_clr))

count_clr <- count_clr[!rownames(count_clr)%in%remove_from_count_file,]
growth_file_quant <- growth_file_quant[!(growth_file_quant$Species)%in%remove_from_growth_file,]
count_clr <- count_clr[,!colnames(count_clr)%in%remove_from_count_file_sample]
growth_file_quant <- growth_file_quant[!(growth_file_quant$Name)%in%remove_from_growth_file_sample,]
count_clr$aaa_species <- rownames(count_clr)
count_clr <- count_clr[,order(colnames(count_clr))]
count_clr_l <- gather(count_clr, key="patient", value="clr", c(2:ncol(count_clr)))
colnames(count_clr_l) <- c("Species", "Name", "clr")

# generate phylum names from species/genus information
new_dataset <- growth_file_quant %>% inner_join(count_clr_l, by=c("Species","Name"))
new_dataset_2 <- new_dataset %>% inner_join(metadata_0, by=c("Name"))
new_dataset_2$State <- ifelse(grepl("X",new_dataset_2$Name), "Healthy", "CF")
new_dataset_2$Phylum <- unlist(map(str_split(new_dataset_2$Species, " ",simplify = FALSE),1))
new_dataset_2$Phylum <- gsub("Abiotrophia","Bacillota", new_dataset_2$Phylum)
new_dataset_2$Phylum <- gsub("Actinomyces","Actinomycetota", new_dataset_2$Phylum)
new_dataset_2$Phylum <- gsub("Aggregatibacter","Pseudomonadota", new_dataset_2$Phylum)
new_dataset_2$Phylum <- gsub("Atopobium","Actinomycetota", new_dataset_2$Phylum)
new_dataset_2$Phylum <- gsub("Bifidobacterium","Actinomycetota", new_dataset_2$Phylum)
new_dataset_2$Phylum <- gsub("Campylobacter","Pseudomonadota", new_dataset_2$Phylum)
new_dataset_2$Phylum <- gsub("Capnocytophaga","Bacteroidota", new_dataset_2$Phylum)
new_dataset_2$Phylum <- gsub("Corynebacterium","Actinomycetota", new_dataset_2$Phylum)
new_dataset_2$Phylum <- gsub("Cutibacterium","Actinomycetota", new_dataset_2$Phylum)
new_dataset_2$Phylum <- gsub("Enterococcus","Bacillota", new_dataset_2$Phylum)
new_dataset_2$Phylum <- gsub("Escherichia","Pseudomonadota", new_dataset_2$Phylum)
new_dataset_2$Phylum <- gsub("Eubacterium","Bacillota", new_dataset_2$Phylum)
new_dataset_2$Phylum <- gsub("Enterococcus","Bacillota", new_dataset_2$Phylum)
new_dataset_2$Phylum <- gsub("Fusobacterium","Fusobacteria", new_dataset_2$Phylum)
new_dataset_2$Phylum <- gsub("Gemella","Bacillota", new_dataset_2$Phylum)
new_dataset_2$Phylum <- gsub("Haemophilus","Pseudomonadota", new_dataset_2$Phylum)
new_dataset_2$Phylum <- gsub("Human","Peploviricota", new_dataset_2$Phylum)
new_dataset_2$Phylum <- gsub("Kingella","Pseudomonadota", new_dataset_2$Phylum)
new_dataset_2$Phylum <- gsub("Lactobacillus","Bacillota", new_dataset_2$Phylum)
new_dataset_2$Phylum <- gsub("Lactococcus","Bacillota", new_dataset_2$Phylum)
new_dataset_2$Phylum <- gsub("Lautropia","Pseudomonadota", new_dataset_2$Phylum)
new_dataset_2$Phylum <- gsub("Leptotrichia","Fusobacteria", new_dataset_2$Phylum)
new_dataset_2$Phylum <- gsub("Moraxella","Pseudomonadota", new_dataset_2$Phylum)
new_dataset_2$Phylum <- gsub("Mycoplasma","Bacillota", new_dataset_2$Phylum)
new_dataset_2$Phylum <- gsub("Neisseria","Pseudomonadota", new_dataset_2$Phylum)
new_dataset_2$Phylum <- gsub("Parvimonas","Bacillota", new_dataset_2$Phylum)
new_dataset_2$Phylum <- gsub("Porphyromonas","Bacteroidota", new_dataset_2$Phylum)
new_dataset_2$Phylum <- gsub("Prevotella","Bacteroidota", new_dataset_2$Phylum)
new_dataset_2$Phylum <- gsub("Pseudomonas","Pseudomonadota", new_dataset_2$Phylum)
new_dataset_2$Phylum <- gsub("Pseudopropionibacterium","Actinomycetota", new_dataset_2$Phylum)
new_dataset_2$Phylum <- gsub("Rothia","Actinomycetota", new_dataset_2$Phylum)
new_dataset_2$Phylum <- gsub("Schaalia","Actinomycetota", new_dataset_2$Phylum)
new_dataset_2$Phylum <- gsub("Selenomonas","Bacillota", new_dataset_2$Phylum)
new_dataset_2$Phylum <- gsub("Staphylococcus","Bacillota", new_dataset_2$Phylum)
new_dataset_2$Phylum <- gsub("Streptococcus","Bacillota", new_dataset_2$Phylum)
new_dataset_2$Phylum <- gsub("Tannerella","Bacteroidota", new_dataset_2$Phylum)
new_dataset_2$Phylum <- gsub("Veillonella","Bacillota", new_dataset_2$Phylum)
new_dataset_2$Phylum <- gsub("Xanthomonas","Pseudomonadota", new_dataset_2$Phylum)
new_dataset_2$Phylum <- gsub("Serratia","Pseudomonadota", new_dataset_2$Phylum)
new_dataset_2$Phylum <- gsub("Finegoldia","Bacillota", new_dataset_2$Phylum)
new_dataset_2$Phylum <- gsub("Stenotrophomonas","Pseudomonadota", new_dataset_2$Phylum)

new_dataset_3 <- select(new_dataset_2, c("Name", "Phylum", "Growth_Rate"))
new_dataset_3_clr <- select(new_dataset_2, c("Name", "Phylum", "clr"))

# merge growth rate and abundance information tables
new_dataset_4_clr <-
  new_dataset_3_clr %>% 
  group_by(Name, Phylum) %>%
  dplyr::mutate(clr_median = mean(clr, na.rm = FALSE))
new_dataset_4_clr$clr <- NULL  
new_dataset_5_clr <- unique(new_dataset_4_clr)
new_dataset_5_clr$State <- ifelse(grepl("X", new_dataset_5_clr$Name), "Healthy", "CF")
new_dataset_5_clr <- inner_join(new_dataset_5_clr, metadata_0, by=c("Name"))

new_dataset_4 <-
  new_dataset_3 %>% 
  group_by(Name, Phylum) %>%
  dplyr::mutate(growth_rate_median = mean(Growth_Rate, na.rm = FALSE))
new_dataset_4$Growth_Rate <- NULL  
new_dataset_5 <- unique(new_dataset_4)
new_dataset_5$State <- ifelse(grepl("X", new_dataset_5$Name), "Healthy", "CF")
new_dataset_5 <- inner_join(new_dataset_5, metadata_0, by=c("Name"))

merged <- inner_join(new_dataset_5_clr, new_dataset_5, by=c("Name", "Phylum"))
merged_healthy <- filter(merged, State.x == "Healthy")

# Test correlation between abundance and growth rate (in healthy children)
cor.test(merged_healthy$clr_median, merged_healthy$growth_rate_median, method="spearman")
spearman.ci(merged_healthy$clr_median, merged_healthy$growth_rate_median, nrep = 1000, conf.level = 0.95)
# Test correlation between abundance and growth rate (in CF children)
merged_cf <- filter(merged, State.x == "CF")
cor.test(merged_cf$clr_median, merged_cf$growth_rate_median, method="spearman")
spearman.ci(merged_cf$clr_median, merged_cf$growth_rate_median, nrep = 1000, conf.level = 0.95)


# generate circular plot for growth rates
new_dataset_5$Patient <- as.factor(as.character(new_dataset_5$Patient))
empty_bar = 50
to_add <- data.frame(matrix(NA, empty_bar*nlevels(new_dataset_5$Patient), ncol(new_dataset_5)))
colnames(to_add) <- colnames(new_dataset_5)
to_add$Name <- rep(levels(new_dataset_5$Patient), each=empty_bar)
new_dataset_5_2 <- rbind(new_dataset_5, to_add)
label_data_growth <- select(new_dataset_5_2, c(Name, growth_rate_median))
label_data_growth_2 <- ddply(label_data_growth, "Name", numcolwise(sum))
label_data_growth_2$row_num <- rownames(label_data_growth_2)
label_data_growth_2$row_num <- as.numeric(as.character(label_data_growth_2$row_num))
label_data_growth_2$angle <- with(label_data_growth_2, ifelse(row_num < nrow(label_data_growth_2), 90, 180))
number_of_bar <- nrow(label_data_growth_2)
label_data_growth_2$angle <- 90 - 360 * (label_data_growth_2$row_num-0.5) / number_of_bar
label_data_growth_2$angle_2 <- ifelse(label_data_growth_2$angle <  -110, label_data_growth_2$angle+180, label_data_growth_2$angle)
label_data_growth_2$hjust <- ifelse(label_data_growth_2$angle < -120, 0.8, 0.5)
label_data_growth_2$Identifier_2 <- as.character(as.factor(label_data_growth_2$Name))
label_data_growth_2$label_colour <- ifelse(grepl("X",label_data_growth_2$Identifier_2), 'Healthy', 'CF')

plot_growth <-
  ggplot(new_dataset_5_2) +
  geom_col(aes(x=Name, y=growth_rate_median, fill=reorder(Phylum, -growth_rate_median))) + 
  coord_polar(start=0) + ylim(-1.2,25) + theme_minimal() +
  geom_text(data=label_data_growth_2, aes(x=Identifier_2, 
                                          y=10, 
                                          label=Identifier_2, 
                                          angle=angle_2, 
                                          hjust=hjust, 
                                          colour=label_colour), size=2, alpha=0.6) +
  scale_fill_manual(values=c("Bacillota"="forestgreen",
                             "Pseudomonadota"="lightcyan3",
                             "Actinomycetota"="khaki",
                             "Bacteroidota"="skyblue4",
                             "Fusobacteria"="rosybrown2",
                             "Peploviricota"="sienna3",
                             "Ascomycota"="cyan")) +
  scale_colour_manual(values = c("Healthy"="darkblue", "CF"="darkred")) +
  labs(colour="Phylum", fill="State") +
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        panel.grid = element_blank(), plot.margin = unit(c(-12,-12,-12,-12), "cm"), legend.position = "none")


# generate circular plot for phylum abundance
count_df_phylum_clr_00_t_long$patient <- as.factor(as.character(count_df_phylum_clr_00_t_long$patient))
to_add_phylum <- data.frame(matrix(NA, empty_bar*nlevels(count_df_phylum_clr_00_t_long$patient), ncol(count_df_phylum_clr_00_t_long)))
colnames(to_add_phylum) <- colnames(count_df_phylum_clr_00_t_long)
to_add_phylum$Identifier <- rep(levels(count_df_phylum_clr_00_t_long$patient), each=empty_bar)
count_df_phylum_clr_00_t_long_2 <- rbind(count_df_phylum_clr_00_t_long, to_add_phylum)
label_data <- select(count_df_phylum_clr_00_t_long_2, c(Identifier, CLR))
label_data_2 <- ddply(label_data, "Identifier", numcolwise(sum))
label_data_2$row_num <- rownames(label_data_2)
label_data_2$row_num <- as.numeric(as.character(label_data_2$row_num))
label_data_2$angle <- with(label_data_2, ifelse(row_num < nrow(label_data_2), 90, 180))
number_of_bar_clr <- nrow(label_data_2)
label_data_2$angle <- 90 - 360 * (label_data_2$row_num-0.5) / number_of_bar_clr
label_data_2$angle_2 <- ifelse(label_data_2$angle <  -110, label_data_2$angle+180, label_data_2$angle)
label_data_2$hjust <- ifelse(label_data_2$angle < -120, 0.8, 0.5)
label_data_2$Identifier_2 <- as.character(as.factor(label_data_2$Identifier))
label_data_2$label_colour <- ifelse(grepl("X",label_data_2$Identifier_2), 'Healthy', 'CF')
count_df_phylum_clr_00_t_long_2$CLR_2 <- ifelse(count_df_phylum_clr_00_t_long_2$CLR >= 0, count_df_phylum_clr_00_t_long_2$CLR, 0)
count_df_phylum_clr_00_t_long_2$Phylum <- str_replace_all(count_df_phylum_clr_00_t_long_2$Phylum, "Actinobacteria", "Actinomycetota")
count_df_phylum_clr_00_t_long_2$Phylum <- str_replace_all(count_df_phylum_clr_00_t_long_2$Phylum, "Firmicutes", "Bacillota")
count_df_phylum_clr_00_t_long_2$Phylum <- str_replace_all(count_df_phylum_clr_00_t_long_2$Phylum, "Proteobacteria", "Pseudomonadota")
count_df_phylum_clr_00_t_long_2$Phylum <- str_replace_all(count_df_phylum_clr_00_t_long_2$Phylum, "Bacteroidetes", "Bacteroidota")

phylum_plot <- ggplot(count_df_phylum_clr_00_t_long_2) +
  geom_col(aes(x=Identifier, 
               y=CLR_2, 
               fill=reorder(Phylum, -CLR_2))) + 
  coord_polar(start=0) + 
  theme_minimal() + 
  ylim(-5,111) + 
  geom_text(data=label_data_2, aes(x=Identifier, 
                                   y=45, 
                                   label=Identifier, 
                                   angle=angle_2, 
                                   hjust=hjust, 
                                   colour=label_colour), size=2, alpha=0.6) +
  scale_fill_manual(values = c("Bacillota"="forestgreen",
                               "Pseudomonadota"="lightcyan3",
                               "Actinomycetota"="khaki",
                               "Bacteroidota"="skyblue4",
                               "Fusobacteria"="rosybrown2",
                               "Peploviricota"="sienna3",
                               "Ascomycota"="cyan")) +
  scale_colour_manual(values = c("Healthy"="darkblue", "CF"="darkred")) +
  labs(colour="Phylum", fill="State") +
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        panel.grid = element_blank(), plot.margin = unit(c(-12,-12,-12, -12), "cm"), legend.position = "none")

# generate plot for the legend
phylum_plot_for_legend <- ggplot(count_df_phylum_clr_00_t_long_2) +
  geom_col(aes(x=Identifier, 
               y=CLR_2, 
               fill=reorder(Phylum, -CLR_2))) + 
  coord_polar(start=0) + 
  theme_minimal() + 
  ylim(-11,100) + 
  geom_text(data=label_data_2, aes(x=Identifier, 
                                   y=45, 
                                   label=Identifier, 
                                   angle=angle_2, 
                                   hjust=hjust, 
                                   colour=label_colour), size=2, alpha=0.6) +
  geom_hline(yintercept=0, size=0.1, colour="darkblue") +
  scale_fill_manual("Phylum",values = c("Bacillota"="forestgreen",
                                        "Pseudomonadota"="lightcyan3",
                                        "Actinomycetota"="khaki",
                                        "Bacteroidota"="skyblue4",
                                        "Fusobacteria"="rosybrown2",
                                        "Peploviricota"="sienna3",
                                        "Ascomycota"="cyan")) +
  scale_colour_manual("State",values = c("Healthy"="darkblue", "CF"="darkred"),
                      labels = c("Healthy (cross-sectional)", "CF (longitudinal)")) +
  guides(colour=guide_legend(override.aes = list(size=5))) +
  theme(axis.text = element_blank(), axis.title = element_blank(), 
        panel.grid = element_blank(), plot.margin = unit(c(-12.5,-12.5,-12.5,-12.5), "cm"), legend.position = "right")

# Extract the legend. Returns a gtable
leg <- get_legend(phylum_plot_for_legend)

# Convert to a ggplot and print
leg_ggplot <- as_ggplot(leg)



# investigate growth rate profiles
growth_file_quant_s <- growth_file_quant
growth_file_quant_s$Growth_class <- NULL
growth_file_quant_s_1 <- spread(growth_file_quant_s, Name, Growth_Rate, 0)
rownames(growth_file_quant_s_1) <- growth_file_quant_s_1$Species
growth_file_quant_s_1$Species <- NULL
growth_file_quant_s_1_t <- data.frame(t(growth_file_quant_s_1))
growth_file_quant_s_1_t <- growth_file_quant_s_1_t[order(rownames(growth_file_quant_s_1_t)),]
growth_file_quant_s_1_t_CF <- growth_file_quant_s_1_t[!grepl("X", rownames(growth_file_quant_s_1_t)),]
growth_file_quant_s_1_t_H <- growth_file_quant_s_1_t[grepl("X", rownames(growth_file_quant_s_1_t)),]
metadata_0$State <- ifelse(grepl("X", metadata_0$Name), "Healthy", "CF")

growth_file_quant_s_1_t_CF <- data.frame(t(growth_file_quant_s_1_t_CF))
growth_file_quant_s_1_t_CF[is.na(growth_file_quant_s_1_t_CF)] <- 0
growth_file_quant_s_1_t_H <- data.frame(t(growth_file_quant_s_1_t_H))
growth_file_quant_s_1_t_H[is.na(growth_file_quant_s_1_t_H)] <- 0



# perform principal component analysis (CF)
pca_CF <- prcomp(growth_file_quant_s_1_t_CF, scale. = FALSE, center = TRUE)
pca_df <- data.frame(pca_CF$rotation)
pca_df$Age_in_weeks <- MCF_metadata$age_in_weeks
eigs_CF <- pca_CF$sdev^2
dim1_CF <- round((eigs_CF[1] / sum(eigs_CF))*100)
dim2_CF <- round((eigs_CF[2] / sum(eigs_CF))*100)

# perform principal component analysis (healthy)
pca_H <- prcomp(growth_file_quant_s_1_t_H, scale. = FALSE, center = TRUE)
eigs_H <- pca_H$sdev^2
dim1_H <- round((eigs_H[1] / sum(eigs_H))*100)
dim2_H <- round((eigs_H[2] / sum(eigs_H))*100)
pca_H_df <- data.frame(pca_H$rotation)
H_metadata <- metadata_0[grepl("X", metadata_0$Name),]
pca_H_df$Age_in_weeks <- H_metadata$Age_in_weeks


# PCA CF plot
cf_plot <- ggplot(pca_df) +
  geom_point(aes(x=PC1, y=PC2, colour=Age_in_weeks, size=Age_in_weeks), shape=16) + ylim(-0.35,0.3) + xlim(0,0.21) + 
  scale_colour_gradient2(name="Age (in weeks)", 
                         low=muted("red"), 
                         mid= "khaki", 
                         high=muted("blue"), 
                         midpoint=100,
                         guide = guide_legend(title.position = "top", ncol = 1)) +
  scale_size(name="Age (in weeks)", range = c(1,5)) + theme_bw() + xlab(paste("PC1 ","(",dim1_CF, "%)")) + ylab(paste("PC2 ","(",dim2_CF, "%)")) +
  
  theme(panel.grid = element_blank(), legend.position = "bottom", legend.direction = "vertical") +
  guides(colour=guide_legend(nrow=1,byrow=TRUE, title.hjust = 0.5))

# PCA healthy plot
healthy_plot <-
  ggplot(pca_H_df) +
  geom_point(aes(x=PC1*-1, y=PC2*1, colour=Age_in_weeks, size=Age_in_weeks), shape=17) + ylim(-0.35,0.3) + xlim(0,0.21) + 
  scale_colour_gradient2(name="Age (in weeks)", 
                         low=muted("red"), 
                         mid= "khaki", 
                         high=muted("blue"), 
                         midpoint=100,
                         guide = guide_legend(title.position = "top", ncol = 1)) +
  scale_size(name="Age (in weeks)", range = c(1,5)) + theme_bw() + xlab(paste("PC1 ","(",dim1_H, "%)")) + ylab(paste("PC2 ","(",dim2_H, "%)")) +
  theme(panel.grid = element_blank(), legend.position = "bottom", legend.direction = "vertical") +
  guides(colour=guide_legend(nrow=1,byrow=TRUE, title.hjust = 0.5))


# Calculate distance between all last CF samples
remove_zero = 0.0000000
colnames(growth_file_quant_s_1_t_CF) <- str_replace_all(colnames(growth_file_quant_s_1_t_CF), "\\.","-")

list_last_samples <- c("CF-A-09", "CF-B-04", "CF-C-03", "CF-D-02", "CF-E-06", "CF-F-04", "CF-G-06", "CF-H-06", "CF-I-04", "CF-J-04", "CF-K-05", "CF-L-04", "CF-M-02")
growth_diversity_t_last <- select(growth_file_quant_s_1_t_CF, all_of(list_last_samples))
braycurtis_last_CF = vegdist(t(growth_diversity_t_last), "bray")
braycurtis_df_last <- as.data.frame(as.matrix(braycurtis_last_CF))
braycurtis_df_last_list <- c(braycurtis_df_last$`CF-A-09`, 
                             braycurtis_df_last$`CF-B-04`,
                             braycurtis_df_last$`CF-C-03`,
                             braycurtis_df_last$`CF-D-02`,
                             braycurtis_df_last$`CF-E-06`,
                             braycurtis_df_last$`CF-F-04`,
                             braycurtis_df_last$`CF-G-06`,
                             braycurtis_df_last$`CF-H-06`,
                             braycurtis_df_last$`CF-I-04`,
                             braycurtis_df_last$`CF-J-04`,
                             braycurtis_df_last$`CF-K-05`,
                             braycurtis_df_last$`CF-L-04`,
                             braycurtis_df_last$`CF-M-02`)
braycurtis_df_last_list <- braycurtis_df_last_list[!braycurtis_df_last_list %in% remove_zero]
braycurtis_df_last_list <- data.frame(braycurtis_df_last_list)
braycurtis_df_last_list$Type <- "t=last"
colnames(braycurtis_df_last_list) <- c("Distance", "Type")
braycurtis_df_last_list$State <- "CF"

# # Calculate distance between all first CF samples
list_first_samples <- colnames(select(growth_file_quant_s_1_t_CF, ends_with("01")))
growth_diversity_t_first <- select(growth_file_quant_s_1_t_CF, all_of(list_first_samples))
braycurtis_first_CF = vegdist(t(growth_diversity_t_first), "bray")
braycurtis_df_first <- as.data.frame(as.matrix(braycurtis_first_CF))
braycurtis_df_first_list <- c(braycurtis_df_first$`CF-A-01`, 
                              braycurtis_df_first$`CF-B-01`,
                              braycurtis_df_first$`CF-C-01`,
                              braycurtis_df_first$`CF-D-01`,
                              braycurtis_df_first$`CF-E-01`,
                              braycurtis_df_first$`CF-F-01`,
                              braycurtis_df_first$`CF-G-01`,
                              braycurtis_df_first$`CF-H-01`,
                              braycurtis_df_first$`CF-I-01`,
                              braycurtis_df_first$`CF-J-01`,
                              braycurtis_df_first$`CF-K-01`,
                              braycurtis_df_first$`CF-L-01`,
                              braycurtis_df_first$`CF-M-01`)
braycurtis_df_first_list <- braycurtis_df_first_list[!braycurtis_df_first_list %in% remove_zero]
braycurtis_df_first_list <- data.frame(braycurtis_df_first_list)
braycurtis_df_first_list$Type <- "t=first"
colnames(braycurtis_df_first_list) <- c("Distance", "Type")
braycurtis_df_first_list$State <- "CF"

# Calculate distance between healthy infants (<200 weeks) and (>200 weeks of age)
metadata_0_H <- filter(metadata_0, State == "Healthy")
growth_file_quant_s_1_t_H_2 <- data.frame(t(growth_file_quant_s_1_t_H))
growth_file_quant_s_1_t_H_2$Age_in_weeks <- metadata_0_H$Age_in_weeks
growth_file_quant_s_1_t_H_2_small <- filter(growth_file_quant_s_1_t_H_2, Age_in_weeks < 200)
growth_file_quant_s_1_t_H_2_small$Age_in_weeks <- NULL
growth_file_quant_s_1_t_H_3_small <- data.frame(t(growth_file_quant_s_1_t_H_2_small))
braycurtis_first_h = vegdist(t(growth_file_quant_s_1_t_H_3_small), "bray")
braycurtis_df_h_first <- as.data.frame(as.matrix(braycurtis_first_h))
braycurtis_df_h_first_list <- c(braycurtis_df_h_first$X02,
                                braycurtis_df_h_first$X03,
                                braycurtis_df_h_first$X04,
                                braycurtis_df_h_first$X07,
                                braycurtis_df_h_first$X10,
                                braycurtis_df_h_first$X11,
                                braycurtis_df_h_first$X12,
                                braycurtis_df_h_first$X13,
                                braycurtis_df_h_first$X14,
                                braycurtis_df_h_first$X16,
                                braycurtis_df_h_first$X17,
                                braycurtis_df_h_first$X18,
                                braycurtis_df_h_first$X22,
                                braycurtis_df_h_first$X24,
                                braycurtis_df_h_first$X32,
                                braycurtis_df_h_first$X34,
                                braycurtis_df_h_first$X36,
                                braycurtis_df_h_first$X37,
                                braycurtis_df_h_first$X39,
                                braycurtis_df_h_first$X41,
                                braycurtis_df_h_first$X42,
                                braycurtis_df_h_first$X43,
                                braycurtis_df_h_first$X44,
                                braycurtis_df_h_first$X45,
                                braycurtis_df_h_first$X46,
                                braycurtis_df_h_first$X47,
                                braycurtis_df_h_first$X48,
                                braycurtis_df_h_first$X49,
                                braycurtis_df_h_first$X50,
                                braycurtis_df_h_first$X51,
                                braycurtis_df_h_first$X52,
                                braycurtis_df_h_first$X53,
                                braycurtis_df_h_first$X55,
                                braycurtis_df_h_first$X56,
                                braycurtis_df_h_first$X57,
                                braycurtis_df_h_first$X58)

braycurtis_df_h_first_list <- braycurtis_df_h_first_list[!braycurtis_df_h_first_list %in% remove_zero]
braycurtis_df_h_first_list <- data.frame(braycurtis_df_h_first_list)
braycurtis_df_h_first_list$Type <- "<200 weeks"
colnames(braycurtis_df_h_first_list) <- c("Distance", "Type")
braycurtis_df_h_first_list$State <- "Healthy"

growth_file_quant_s_1_t_H_2_large <- filter(growth_file_quant_s_1_t_H_2, Age_in_weeks > 200)
growth_file_quant_s_1_t_H_2_large$Age_in_weeks <- NULL
growth_file_quant_s_1_t_H_2_large <- data.frame(t(growth_file_quant_s_1_t_H_2_large))
braycurtis_last_h = vegdist(t(growth_file_quant_s_1_t_H_2_large), "bray")
braycurtis_df_h_last <- as.data.frame(as.matrix(braycurtis_last_h))
braycurtis_df_h_last_list <- c(braycurtis_df_h_last$X01,
                               braycurtis_df_h_last$X05,
                               braycurtis_df_h_last$X06,
                               braycurtis_df_h_last$X15,
                               braycurtis_df_h_last$X19,
                               braycurtis_df_h_last$X21,
                               braycurtis_df_h_last$X25,
                               braycurtis_df_h_last$X26,
                               braycurtis_df_h_last$X27,
                               braycurtis_df_h_last$X28,
                               braycurtis_df_h_last$X29,
                               braycurtis_df_h_last$X30,
                               braycurtis_df_h_last$X33,
                               braycurtis_df_h_last$X35,
                               braycurtis_df_h_last$X38,
                               braycurtis_df_h_last$X40)
braycurtis_df_h_last_list <- braycurtis_df_h_last_list[!braycurtis_df_h_last_list %in% remove_zero]
braycurtis_df_h_last_list <- data.frame(braycurtis_df_h_last_list)
braycurtis_df_h_last_list$Type <- ">200 weeks"
colnames(braycurtis_df_h_last_list) <- c("Distance", "Type")
braycurtis_df_h_last_list$State <- "Healthy"

# merge CF and healthy distance plot
bothit <- data.frame(rbind(braycurtis_df_first_list, braycurtis_df_last_list, braycurtis_df_h_first_list, braycurtis_df_h_last_list))

# obtain statistics
df_p_val <- bothit %>% 
  rstatix::group_by(State) %>% 
  wilcox_test(Distance ~ Type) %>% 
  rstatix::add_xy_position()

# plot statistics
my_boxplot <-
  ggplot(bothit, aes(x=Type, y=Distance)) +
  geom_violin(aes(fill=State), width=0.8, alpha=0.4) +
  geom_jitter(width = 0.08, size=0.3, alpha=0.3) + 
  geom_boxplot(width = 0.2, colour="orangered4", alpha=0.1) + 
  ylim(0,1.2) + facet_wrap(~State, scales = "free_x") + xlab(" 
                                                             ") +
  ylab("Bray-Curtis distance") + theme_bw() + xlab(" ") +
  theme(panel.background = element_rect(fill=NA), 
        panel.grid = element_blank(),
        strip.background = element_rect(fill = "white"),
        legend.position = "none",
        plot.margin = unit(x=c(2,5,17,0), units= "mm")) +
  scale_fill_manual(values=c("CF"="gray70", "Healthy"="palegreen1")) +
  geom_text(aes(x=1.5, y=1.1), label="* * * *", size=3) +
  geom_segment(aes(x=1, xend=2, y=1.08, yend=1.08), size=0.2)

# calculate effect sizes
bothit_h <- filter(bothit, State == "Healthy")
bothit_cf <- filter(bothit, State == "CF")
rcompanion::wilcoxonR(bothit_h$Distance, g=bothit_h$Type, ci=TRUE) # R: 0.51, CI=0.45-0.54
wilcox.test(bothit_h$Distance, g=bothit_h$Type) # V = 1125750, p-value < 2.2e-16
rcompanion::wilcoxonR(bothit_cf$Distance, g=bothit_cf$Type, ci=TRUE) # R: 0.24, CI=0.13-0.35
wilcox.test(bothit_cf$Distance, g=bothit_cf$Type) # V = 48828, p-value < 2.2e-16

# merge and generate Figure 1
first_plot <- ggarrange(cf_plot, healthy_plot, nrow=1, common.legend = TRUE, legend="bottom", labels = c("C", "D"))
second_plot <- ggarrange(first_plot, my_boxplot, nrow=1, labels = c("C", "E"), widths = c(1,0.5))
phylum_plot_merged <- ggarrange(phylum_plot, plot_growth, leg_ggplot, nrow=1, widths = c(1,1.1,0.4), labels=c("A", "B"))
final_plot_figure1 <- ggarrange(phylum_plot_merged, second_plot, nrow=2, heights = c(1,0.7))




# Rothia mucilaginosa SNP analysis
rothia_snp <- read_delim("python_generated_files/rothia_snp.csv", ";", escape_double = FALSE, trim_ws = TRUE)
rothia_snp <- data.frame(rothia_snp)
rownames(rothia_snp) <- rothia_snp$Variant
rothia_snp$Variant <- NULL
rothia_snp_df <- rothia_snp[,!grepl("Neg",colnames(rothia_snp))]
rothia_snp_df_t <- data.frame(t(rothia_snp_df))
rothia_snp_df_t + 1
rothia_snp_df_t[rothia_snp_df_t == -99] <- 0

# Perform UMAP
custom.config = umap.defaults
custom.config$random_state = 111 
custom.config$n_neighbors = 9 
custom.config$n_components = 2
custom.config$metric = "euclidean"
custom.config$spread = 8 
rothia.umap = umap(rothia_snp_df_t, config=custom.config)
# store UMAP results in dataframe and add metadata
rothia.umap_df <- data.frame(rothia.umap$layout)
rothia.umap_df <- rothia.umap_df[order(rownames(rothia.umap_df)),]
rothia.umap_df$State <- ifelse(grepl("X", rownames(rothia.umap_df)), "Healthy", "CF")
rothia.umap_df$patient <- metadata_0$Patient
rothia.umap_df$Age <- metadata_0$Age_in_weeks
rothia.umap_df$sample_no <- metadata_0$Sample_No
rothia.umap_df$sample_no_2 <- ifelse(rothia.umap_df$patient == "X", NA, rothia.umap_df$sample_no)
rothia.umap_df$Umap_cluster <- with(rothia.umap_df, ifelse(X2 > -4, "1", ifelse(X1 > 5, "3", "2")))
path_info <- rothia.umap_df[rothia.umap_df$patient!="X",]
path_info_h <- rothia.umap_df[rothia.umap_df$patient=="X",]

# plot UMAP results
part1 <-
  ggplot(rothia.umap_df) +
  geom_jitter(aes(x=X1, y=X2, fill=patient, shape=State, size=Age), alpha=0.6) + 
  geom_text(aes(x=X1, y=X2, label=sample_no), size=3) +
  geom_label(aes(x=9,y=20), label="Cluster 1", size=4) +
  geom_label(aes(x=-12,y=-12), label="Cluster 2", size=4) +
  geom_label(aes(x=20,y=-13), label="Cluster 3", size=4) +  theme_bw() +
  
  scale_size("Age [in weeks]", range=c(1,8), breaks = c(50,100,150,200), guide = guide_legend(title.position = "top", nrow = 2)) +
  
  scale_fill_manual("Patient", values=c("X"="grey80", 
                                        "CF-A"="firebrick1", 
                                        "CF-B"="orange", 
                                        "CF-C"="darkolivegreen1", 
                                        "CF-D"="mediumorchid3", 
                                        "CF-E"="rosybrown2",
                                        "CF-F"="seagreen3",
                                        "CF-G"="tan4",
                                        "CF-H"="violet",
                                        "CF-I"="skyblue4",
                                        "CF-J"="magenta4",
                                        "CF-K"="yellow",
                                        "CF-L"="darkslateblue",
                                        "CF-M"="red")) + 
  
  scale_shape_manual("Condition",values=c(22,24), labels=c("CF (longitudinal)", "Healthy (cross-sectional)")) +
  
  guides(shape = guide_legend(override.aes = list(size = 3), title.position = "top", ncol = 1),
         fill = guide_legend(override.aes = list(size = 3, colour=c("X"="grey80", 
                                                                    "CF-A"="firebrick1", 
                                                                    "CF-B"="orange", 
                                                                    "CF-C"="darkolivegreen1", 
                                                                    "CF-D"="mediumorchid3", 
                                                                    "CF-E"="rosybrown2",
                                                                    "CF-F"="seagreen3",
                                                                    "CF-G"="tan4",
                                                                    "CF-H"="violet",
                                                                    "CF-I"="skyblue4",
                                                                    "CF-J"="magenta4",
                                                                    "CF-K"="yellow",
                                                                    "CF-L"="darkslateblue",
                                                                    "CF-M"="red")), title.position = "top", ncol = 2)) +
  
  xlim(-25,25) + 
  ylim(-25,25) + 
  xlab("UMAP1") + 
  ylab("UMAP2") + 
  theme(panel.grid = element_blank(),
        legend.position = "right",
        legend.direction = "vertical",
        legend.box = "vertical",
        legend.box.just = "left")

# obtain clustering information
path_info$sample_no_factor <- factor(path_info$sample_no, levels=c("1", "2", "3", "4", "5", "6", "7", "8", "9"))

# calculate distance to squared group centroids
umap_dist <- dist(rothia.umap_df[1:2], method = "euclidean")
short_table_rothia <- select(rothia.umap_df, c(Umap_cluster))
short_table_rothia$Item <- rownames(short_table_rothia)
rothia.umap_df$Umap_cluster <- factor(rothia.umap_df$Umap_cluster, levels=c("1", "2", "3"))
distances_group_centroid <- data.frame(dist_to_centroids(umap_dist, rothia.umap_df$Umap_cluster, squared = TRUE))
merge_df <- inner_join(distances_group_centroid, short_table_rothia, by="Item")
merge_df$match <- ifelse(merge_df$CentroidGroup == merge_df$Umap_cluster, "match", "no-match")
merge_df_filtered <- merge_df[merge_df$match=="match",]
merge_df_filtered$State <- ifelse(grepl("X",merge_df_filtered$Item), "Healthy", "CF")
merge_df_filtered$Umap_cluster <- factor(merge_df_filtered$Umap_cluster, levels = c("1", "2", "3"))

# plot and compare distances across UMAP clusters
sq_distance_plot <-
  ggline(merge_df_filtered, x = "Umap_cluster", y = "CentroidDistance", 
         add = c("mean_se", "jitter"),
         color = "State", palette = "jco") +
  stat_compare_means(aes(group=State), label = "p.signif", label.y = c(580, 250, 150), size=4.5) +
  theme_bw() + 
  ylim(0,600) + 
  theme(panel.grid = element_blank(),
        legend.position = c(0.8,0.9),
        legend.title = element_blank()) +
  ylab("Squared distance between centroids") + xlab("UMAP clusters")


# generate cluster positions of CF samples
rothia_cf <-
  ggplot(path_info) +
  geom_tile(aes(x=patient,y=sample_no_factor, fill=Umap_cluster), width=0.5, colour="black", height=0.8) + 
  theme_minimal() + coord_polar(start=0) +
  scale_fill_manual("UMAP cluster", values=c("1"="seagreen4", 
                                             "2"="black", 
                                             "3"="tomato")) +
  guides(fill = guide_legend(title.position = "top", ncol = 3)) +
  theme(axis.title = element_blank(), 
        axis.text.y = element_blank(),
        panel.grid = element_blank(), 
        panel.border = element_rect(colour="black", fill=NA),
        legend.position = "bottom")

# gernate cluster positions of healthy children
path_info_h$Age_group <- ifelse(path_info_h$Age > 100, "old", "young")
table(path_info_h$Umap_cluster, path_info_h$Age_group)

variable1 <- c(20/22, 2/22, 0)
variable2 <- c(11/30, 4/30, 15/30)
healthy_rothia <- data.frame(cbind(variable1, variable2))
colnames(healthy_rothia) <- c(">100", "<100")
rownames(healthy_rothia) <- c("UMAP_cluster1", "UMAP_cluster2", "UMAP_cluster3")
healthy_rothia$Umap_cluster <- rownames(healthy_rothia)
healthy_rothia_l <- gather(healthy_rothia, key="Age", value="Percentage", c(1:2))

# plot healthy clustering positions
relAbundance_healthy <-
  ggplot(healthy_rothia_l) +
  geom_col(aes(x=Age, y=Percentage, fill=reorder(Umap_cluster, Percentage)), width=0.9) +
  scale_fill_manual("UMAP cluster",values=c("UMAP_cluster1"="seagreen4", 
                                            "UMAP_cluster2"="black", 
                                            "UMAP_cluster3"="tomato")) +
  xlab("Age [weeks]") + ylab("Relative abundance") + coord_flip() +
  theme_minimal() + theme(legend.position = "none", 
                          panel.grid = element_blank(),
                          panel.border = element_rect(colour="black", fill=NA)) 

# merge figures
first_plot <- ggarrange(rothia_cf, relAbundance_healthy, nrow=2, ncol=1, common.legend=FALSE, heights = c(1, 0.3), widths = c(1,0.5))
final_plot <- ggarrange(part1, first_plot, nrow=1, labels = c("A", "B"), widths = c(1,0.7))

# export figures
pdf("Figure_1.pdf", width = 14, height = 9.5)
final_plot_figure1
dev.off()

pdf("Figure_2.pdf", width = 14, height=8)
final_plot
dev.off()

pdf("Figure_3.pdf", width = 4, height=5)
sq_distance_plot
dev.off()

