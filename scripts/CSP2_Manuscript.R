# Analysis script for CSP2 Manuscript
library(tidyverse)
library(ggplot2)
library(cowplot)
library(ggpmisc)
library(doParallel)
library(ggridges)
library(ggtree)
library(ggnewscale)
library(ape)
library(phangorn)
library(ggnewscale)
library(RColorBrewer)
library(arrow)
library(outliers)
library(emmeans)
library(ggvenn)
library(viridis)
library(broom)

# Functions
getIQROutlier <- function(x, mult = 1.5) {
  q25 <- quantile(x, 0.25)
  q75 <- quantile(x, 0.75)
  iqr <- q75 - q25
  
  low_outlier <- q25 - mult * iqr
  high_outlier <- q75 + mult * iqr
  
  return(c(low_outlier, high_outlier))
}

print_dup_rows <- function(df){
  duplicate_rows <- df[duplicated(df) | duplicated(df, fromLast = TRUE), ]
  print(duplicate_rows)
}

print_na_rows <- function(df){
  df[apply(df, 1, function(x) any(is.na(x))),]
}

read_snpdiffs <- function(file_path, sample_id){
  df <- read_tsv(file_path, comment = "#",show_col_types = FALSE,
                 col_names = c('Ref_Contig','Start_Ref','Ref_Pos',
                               'Query_Contig','Start_Query','Query_Pos',
                               'Ref_Loc','Query_Loc',
                               'Ref_Start','Ref_End',
                               'Query_Start','Query_End',
                               'Ref_Base','Query_Base',
                               'Dist_to_Ref_End','Dist_to_Query_End',
                               'Ref_Aligned','Query_Aligned',
                               'Query_Direction','Perc_Iden','Cat'))
  if(nrow(df) > 0) {
    df <- read_tsv(file_path, 
                   show_col_types = FALSE,
                   comment = "#",
                   col_types = 'ciiciicciiiicciiiiinc',
                   col_names = c('Ref_Contig','Start_Ref','Ref_Pos',
                                 'Query_Contig','Start_Query','Query_Pos',
                                 'Ref_Loc','Query_Loc',
                                 'Ref_Start','Ref_End',
                                 'Query_Start','Query_End',
                                 'Ref_Base','Query_Base',
                                 'Dist_to_Ref_End','Dist_to_Query_End',
                                 'Ref_Aligned','Query_Aligned',
                                 'Query_Direction','Perc_Iden','Cat'))
    return(df %>% mutate(Isolate_ID = sample_id))
  }
}

modified_z_score_test <- function(x) {
  median_x <- median(x)
  mad_x <- mad(x)
  modified_z_score <- sapply(x, function(value) {
    0.6745 * (value - median_x) / mad_x
  })
  
  return(modified_z_score)
}

calculate_ci <- function(est, se) {
  z <- qnorm(0.975)
  lower <- est - z * se
  upper <- est + z * se
  return(paste0(round(lower,3)," - ",round(upper,3)))
}

##### 01: Simulated Data Analysis #####

# Check mutation counts for raw, mutated genomes

# base_mutation_df <- read_tsv(paste0(base_sim_dir,'/Agona_100X_snpListMutated.txt')) %>%
#   mutate(Isolate_ID = paste0("Mut_", replicate)) %>%
#   select(Isolate_ID,Pos = position,Ref_Base = originalBase,Mut_Base = newBase) %>%
#   mutate(Mut_Type = ifelse(Mut_Base == "_deletion","Deletion",
#                            ifelse(str_detect(Mut_Base,"_insertion"),"Insertion","SNP"))) %>%
#   rename(Ref_Pos = "Pos")
#
base_mutation_df <- read_tsv('data/Simulation/Simulated_Mutations.tsv')

# Partition insertions, deletions, and SNPs
base_insertion_df <- base_mutation_df %>%
  filter(Mut_Type == "Insertion") %>%
  mutate(Mut_Base = str_replace(Mut_Base,"_insertion","")) %>%
  mutate(Base_List = str_split(Mut_Base,pattern = "")) %>%
  mutate(Query_Base = Ref_Base,
         Inserted_Base = ifelse(Ref_Base == Base_List[[1]],Base_List[[2]],Base_List[[1]])) %>%
  select(-Query_Base,-Base_List)

base_deletion_df <- base_mutation_df %>%
  filter(Mut_Type == "Deletion") %>%
  mutate(Mut_Base = ".")

base_snp_df <- base_mutation_df %>%
  filter(Mut_Type == "SNP")

# Get snpdiffs files
sim_raw_snpdiffs_dir <- "data/Simulation/snpdiffs/Raw"
sim_skesa_snpdiffs_dir <- "data/Simulation/snpdiffs/SKESA"
sim_spades_snpdiffs_dir <- "data/Simulation/snpdiffs/SPAdes"

raw_snpdiffs_files <- unlist(lapply(sim_raw_snpdiffs_dir, function(x) list.files(x, pattern = "\\.snpdiffs$", full.names = TRUE))) %>%
  set_names(sapply(basename(.), function(x) str_replace(x, "__vs__Agona.snpdiffs", "")))

skesa_screen_snpdiffs_files <- unlist(lapply(sim_skesa_snpdiffs_dir, function(x) list.files(x, pattern = "\\.snpdiffs$", full.names = TRUE))) %>%
  set_names(sapply(basename(.), function(x) str_replace(x, "__vs__Agona.snpdiffs", "")))

spades_screen_snpdiffs_files <- unlist(lapply(sim_spades_snpdiffs_dir, function(x) list.files(x, pattern = "\\.snpdiffs$", full.names = TRUE))) %>%
  set_names(sapply(basename(.), function(x) str_replace(x, "_Spades__vs__Agona.snpdiffs", "")))

# Set directories for screening results
sim_raw_screen_dir <- "data/Simulation/Screening/CSP2_Raw"
sim_skesa_screen_dir <- "data/Simulation/Screening/CSP2_Screen_SKESA"
sim_spades_screen_dir <- "data/Simulation/Screening/CSP2_Screen_SPADES"

raw_screen_log_dir <- paste0(sim_raw_screen_dir,"/logs/Screening_Logs/")
skesa_screen_log_dir <- paste0(sim_skesa_screen_dir,"/logs/Screening_Logs/")
spades_screen_log_dir <- paste0(sim_spades_screen_dir,"/logs/Screening_Logs/")

raw_screen_snp_files <- unlist(lapply(raw_screen_log_dir, function(x) list.files(x, pattern = "\\.tsv$", full.names = TRUE))) %>%
  set_names(sapply(basename(.), function(x) str_replace(x, "__vs__Agona_SNPs.tsv", "")))

skesa_screen_snp_files <- unlist(lapply(skesa_screen_log_dir, function(x) list.files(x, pattern = "\\.tsv$", full.names = TRUE))) %>%
  set_names(sapply(basename(.), function(x) str_replace(x, "__vs__Agona_SNPs.tsv", "")))

spades_screen_snp_files <- unlist(lapply(spades_screen_log_dir, function(x) list.files(x, pattern = "\\.tsv$", full.names = TRUE))) %>%
  set_names(sapply(basename(.), function(x) str_replace(x, "__vs__Agona_SNPs.tsv", "")))

# Check results for raw, mutated genomes versus Agona reference
raw_screen_mummer_df <- read_tsv("data/Simulation/Screening/CSP2_Raw/Raw_MUMmer_Summary.tsv") %>%
  rename(Isolate_ID = "Query_ID") %>%
  mutate(Assembler = "Raw",Dataset = "Raw",Reference_ID = "Agona")

raw_screen_screening_df <- read_tsv("data/Simulation/Screening/CSP2_Raw//Screening_Results.tsv") %>%
  rename(Isolate_ID = "Query_ID") %>%
  mutate(Assembler = "Raw",Dataset = "Raw",Reference_ID = "Agona")

# raw_screen_snp_df <- lapply(names(raw_screen_snp_files), function(sample_id) {
#   read_tsv(raw_screen_snp_files[[sample_id]]) %>% mutate(Isolate_ID = sample_id)}) %>%
#   bind_rows()

raw_screen_snp_df <- read_tsv('data/Simulation/Simulated_Raw_Reference_Screening.tsv')

sim_raw_density_df <- raw_screen_snp_df %>%
  filter(Cat == "Purged_Density") %>%
  select(Isolate_ID,Ref_Contig,Ref_Pos,Ref_Loc,Query_Loc,Ref_Base,Query_Base,Cat)

# Only sites that failed the default density filter were excluded
sim_raw_snp_df <- raw_screen_snp_df %>%
  filter(Cat == "SNP") %>%
  select(Isolate_ID,Ref_Contig,Ref_Pos,Ref_Loc,Query_Loc,Ref_Base,Query_Base,Cat) %>%
  bind_rows(sim_raw_density_df) %>%
  arrange(Isolate_ID,Ref_Loc)

raw_snp_compare_df <- base_snp_df %>%
  rename(Query_Base = "Mut_Base") %>%
  full_join(sim_raw_snp_df, by = c("Isolate_ID","Ref_Pos")) %>%
  arrange(Isolate_ID,Ref_Loc)

# Sanity check
# all(raw_snp_compare_df$Ref_Base.x == raw_snp_compare_df$Ref_Base.y)
# all(raw_snp_compare_df$Query_Base.x == raw_snp_compare_df$Query_Base.y)
# print_na_rows(raw_snp_compare_df)

raw_snp_count_df <- raw_screen_snp_df %>%
  filter(Cat == "SNP")  %>% 
  group_by(Isolate_ID) %>% 
  summarize(Raw_SNP_Count = n())

# Read in assembled isolate data for SKESA and SPAdes
skesa_isolate_data <- lapply(sim_skesa_screen_dir, function(x) read_tsv(file.path(x, "Isolate_Data.tsv")) %>% mutate(Assembler = "SKESA")) %>%
  bind_rows() %>%
  distinct() %>%
  mutate(Dataset = ifelse(!str_detect(Isolate_ID, "_"), "Agona_Reference", sapply(strsplit(as.character(Isolate_ID), "_"), `[`, 1))) %>%
  mutate(ID = ifelse(!str_detect(Isolate_ID, "_"),"Agona_Reference", sapply(strsplit(as.character(Isolate_ID), "_"), function(x) paste(x[-1], collapse = "_")))) %>%
  mutate(Dataset = factor(Dataset, levels = c("Agona_Reference", "Raw",'100X','80X','60X','40X','20X'))) %>%
  select(-ID) %>%
  distinct(Dataset, Isolate_ID, .keep_all = TRUE)

spades_isolate_data <- lapply(sim_spades_screen_dir, function(x) read_tsv(file.path(x, "Isolate_Data.tsv")) %>% mutate(Assembler = "SPAdes")) %>%
  bind_rows() %>%
  distinct() %>%
  mutate(Dataset = ifelse(!str_detect(Isolate_ID, "_"), "Agona_Reference", sapply(strsplit(as.character(Isolate_ID), "_"), `[`, 1))) %>%
  mutate(ID = ifelse(!str_detect(Isolate_ID, "_"),"Agona_Reference", sapply(strsplit(as.character(Isolate_ID), "_"), function(x) paste(x[-1], collapse = "_")))) %>%
  mutate(Dataset = factor(Dataset, levels = c("Agona_Reference", "Raw",'100X','80X','60X','40X','20X'))) %>%
  mutate(Isolate_ID = str_replace(Isolate_ID,"_Spades","")) %>%
  select(-ID) %>%
  distinct(Dataset, Isolate_ID, .keep_all = TRUE)

sim_isolate_data <- bind_rows(skesa_isolate_data,spades_isolate_data) %>%
  mutate(Assembler = ifelse(str_detect(Isolate_ID,"Agona") | str_detect(Isolate_ID,"Raw"),"Raw",Assembler)) %>%
  distinct(Isolate_ID,Dataset,Assembler,.keep_all = TRUE)

sim_isolate_summary_df <- sim_isolate_data %>% 
  filter(Assembler != "Raw") %>%
  group_by(Assembler,Dataset) %>%
  summarize(Min_Contigs = min(Contig_Count),
            Max_Contigs = max(Contig_Count),
            Median_Conitgs = median(Contig_Count),
            Min_Assembly_Bases = min(Assembly_Bases),
            Median_Assembly_Bases = median(Assembly_Bases),
            Max_Assembly_Bases = max(Assembly_Bases))


# Process screening results for SKESA
skesa_screen_mummer_df <- lapply(sim_skesa_screen_dir, function(x) read_tsv(file.path(x, "Raw_MUMmer_Summary.tsv")) %>% mutate(Assembler = "SKESA")) %>%
  bind_rows() %>%
  separate(Query_ID, c("Dataset","Isolate_ID"), sep = "X_") %>%
  mutate(Dataset = paste0(Dataset,"X")) %>%
  mutate(Dataset = factor(Dataset, levels = c('100X','80X','60X','40X','20X'))) %>%
  mutate(Reference_ID = ifelse(str_detect(Reference_ID,"Agona"),"Agona_Ref","Mutated_Genome"))

skesa_screen_screening_df <- lapply(sim_skesa_screen_dir, function(x) read_tsv(file.path(x, "Screening_Results.tsv")) %>% mutate(Assembler = "SKESA")) %>%
  bind_rows() %>%
  separate(Query_ID, c("Dataset","Isolate_ID"), sep = "X_") %>%
  mutate(Dataset = paste0(Dataset,"X")) %>%
  mutate(Dataset = factor(Dataset, levels = c('100X','80X','60X','40X','20X'))) %>%
  mutate(Reference_ID = ifelse(str_detect(Reference_ID,"Agona"),"Agona_Ref","Mutated_Genome"))

# skesa_screen_snpdiffs_df <- lapply(names(skesa_screen_snpdiffs_files), function(sample_id) {
#   read_snpdiffs(skesa_screen_snpdiffs_files[[sample_id]], sample_id)}) %>%
#   bind_rows() %>%
#   mutate(Cat = ifelse(Ref_Base == ".","Insertion",ifelse(Query_Base == ".","Deletion",Cat))) %>%
#   separate(Isolate_ID, c("Dataset","Isolate_ID"), sep = "X_") %>%
#   mutate(Dataset = paste0(Dataset,"X")) %>%
#   mutate(Dataset = factor(Dataset, levels = c('100X','80X','60X','40X','20X'))) %>%
#   mutate(Assembler = "SKESA")

skesa_screen_snpdiffs_df <- read_tsv("data/Simulation/Simulated_SKESA_snpdiffs.tsv")

# skesa_screen_snp_df <- lapply(names(skesa_screen_snp_files), function(sample_id) {
#   read_tsv(skesa_screen_snp_files[[sample_id]]) %>% mutate(CSP2_ID = sample_id)}) %>%
#   bind_rows() %>%
#   mutate(Cat = ifelse(Ref_Base == ".","Insertion",ifelse(Query_Base == ".","Deletion",Cat))) %>%
#   separate(CSP2_ID, c("Dataset","Isolate_ID"), sep = "X_") %>%
#   mutate(Dataset = paste0(Dataset,"X")) %>%
#   mutate(Assembler = "SKESA")

skesa_screen_snp_df <- read_tsv("data/Simulation/Simulated_SKESA_Screening_SNPs.tsv")

# Process screening results for SPAdes
spades_screen_mummer_df <- lapply(sim_spades_screen_dir, function(x) read_tsv(file.path(x, "Raw_MUMmer_Summary.tsv")) %>% mutate(Assembler = "SPAdes")) %>%
  bind_rows() %>%
  separate(Query_ID, c("Dataset","Isolate_ID"), sep = "X_") %>%
  mutate(Dataset = paste0(Dataset,"X")) %>%
  mutate(Dataset = factor(Dataset, levels = c('100X','80X','60X','40X','20X'))) %>%
  mutate(Isolate_ID = str_replace(Isolate_ID,"_Spades","")) %>%
  mutate(Reference_ID = ifelse(str_detect(Reference_ID,"Agona"),"Agona_Ref","Mutated_Genome"))

spades_screen_screening_df <- lapply(sim_spades_screen_dir, function(x) read_tsv(file.path(x, "Screening_Results.tsv")) %>% mutate(Assembler = "SPAdes")) %>%
  bind_rows() %>%
  separate(Query_ID, c("Dataset","Isolate_ID"), sep = "X_") %>%
  mutate(Isolate_ID = str_replace(Isolate_ID,"_Spades","")) %>%
  mutate(Dataset = paste0(Dataset,"X")) %>%
  mutate(Dataset = factor(Dataset, levels = c('100X','80X','60X','40X','20X'))) %>%
  mutate(Reference_ID = ifelse(str_detect(Reference_ID,"Agona"),"Agona_Ref","Mutated_Genome"))

# spades_screen_snpdiffs_df <- lapply(names(spades_screen_snpdiffs_files), function(sample_id) {
#   read_snpdiffs(spades_screen_snpdiffs_files[[sample_id]], sample_id)}) %>%
#   bind_rows() %>%
#   mutate(Cat = ifelse(Ref_Base == ".","Insertion",ifelse(Query_Base == ".","Deletion",Cat))) %>%
#   separate(Isolate_ID, c("Dataset","Isolate_ID"), sep = "X_") %>%
#   mutate(Dataset = paste0(Dataset,"X")) %>%
#   mutate(Dataset = factor(Dataset, levels = c('100X','80X','60X','40X','20X'))) %>%
#   mutate(Isolate_ID = str_replace(Isolate_ID,"_Spades","")) %>%
#   mutate(Assembler = "SPAdes")

spades_screen_snpdiffs_df <- read_tsv("data/Simulation/Simulated_SPAdes_snpdiffs.tsv")

# spades_screen_snp_df <- lapply(names(spades_screen_snp_files), function(sample_id) {
#   read_tsv(spades_screen_snp_files[[sample_id]]) %>% mutate(CSP2_ID = sample_id)}) %>%
#   bind_rows() %>%
#   mutate(Cat = ifelse(Ref_Base == ".","Insertion",ifelse(Query_Base == ".","Deletion",Cat))) %>%
#   separate(CSP2_ID, c("Dataset","Isolate_ID"), sep = "X_") %>%
#   mutate(Isolate_ID = str_replace(Isolate_ID,"_Spades","")) %>%
#   mutate(Dataset = paste0(Dataset,"X")) %>%
#   mutate(Assembler = "SPAdes")

spades_screen_snp_df <- read_tsv("data/Simulation/Simulated_SPAdes_Screening_SNPs.tsv")

# Compare reference coverage
sim_ref_coverage_df <- bind_rows(skesa_screen_screening_df,spades_screen_screening_df) %>%
  select(Isolate_ID,Dataset,Assembler,Reference_Percent_Aligned)

sim_query_coverage_df <- bind_rows(skesa_screen_screening_df,spades_screen_screening_df) %>%
  select(Isolate_ID,Dataset,Assembler,Query_Percent_Aligned)

# Read in SNP Pipeline output
mut_ids <- base_mutation_df %>% pull(Isolate_ID) %>% unique()

sim_skesa_snp_files <- list.files(path ="data/Simulation/SNP/CSP2_SNP_SKESA",full.names = TRUE,recursive = TRUE)
sim_spades_snp_files <- list.files(path ="data/Simulation/SNP/CSP2_SNP_SPADES",full.names = TRUE,recursive = TRUE)

sim_skesa_snp_distance_files <- sim_skesa_snp_files[str_detect(sim_skesa_snp_files,"snp_distance_pairwise.tsv")] 
sim_spades_snp_distance_files <- sim_spades_snp_files[str_detect(sim_spades_snp_files,"snp_distance_pairwise.tsv")]

sim_skesa_query_coverage_files <- sim_skesa_snp_files[str_detect(sim_skesa_snp_files,"Query_Coverage.tsv")]
sim_spades_query_coverage_files <- sim_spades_snp_files[str_detect(sim_spades_snp_files,"Query_Coverage.tsv")]

sim_skesa_locus_coverage_files <- sim_skesa_snp_files[str_detect(sim_skesa_snp_files,"Locus_Categories.tsv")]
sim_spades_locus_coverage_files <- sim_spades_snp_files[str_detect(sim_spades_snp_files,"Locus_Categories.tsv")]

# Read in SNP distance files

# sim_snp_distance_df <- tibble()
# for(mut_id in mut_ids){
#   
#   skesa_distance_file <- sim_skesa_snp_distance_files[str_detect(sim_skesa_snp_distance_files,paste0("/",mut_id,"/"))]
#   spades_distance_file <- sim_spades_snp_distance_files[str_detect(sim_spades_snp_distance_files,paste0("/",mut_id,"/"))]
#   
#   sim_snp_distance_df <- sim_snp_distance_df %>%
#     bind_rows(read_tsv(skesa_distance_file) %>% mutate(Isolate_ID = mut_id,Assembler = "SKESA")) %>%
#     bind_rows(read_tsv(spades_distance_file) %>% mutate(Isolate_ID = mut_id,Assembler = "SPAdes"))
# }
# 

sim_snp_distance_df <- read_tsv("data/Simulation/Simulated_SNP_Distance.tsv")

# Read in query coverage data

# sim_query_coverage_df <- tibble()
# for(mut_id in mut_ids){
# 
#   skesa_query_file <- sim_skesa_query_coverage_files[str_detect(sim_skesa_query_coverage_files,paste0("/",mut_id,"/"))]
#   spades_query_file <- sim_spades_query_coverage_files[str_detect(sim_spades_query_coverage_files,paste0("/",mut_id,"/"))]
# 
#   sim_query_coverage_df <- sim_query_coverage_df %>%
#     bind_rows(read_tsv(skesa_query_file) %>% mutate(Isolate_ID = mut_id,Assembler = "SKESA")) %>%
#     bind_rows(read_tsv(spades_query_file) %>% mutate(Isolate_ID = mut_id,Assembler = "SPAdes"))
# }

sim_query_coverage_df <- read_tsv("data/Simulation/Simulated_Query_Coverage.tsv")

# Read in locus coverage data

# sim_locus_coverage_df <- tibble()
# for(mut_id in mut_ids){
# 
#   skesa_locus_file <- sim_skesa_locus_coverage_files[str_detect(sim_skesa_locus_coverage_files,paste0("/",mut_id,"/"))]
#   spades_locus_file <- sim_spades_locus_coverage_files[str_detect(sim_spades_locus_coverage_files,paste0("/",mut_id,"/"))]
# 
#   sim_locus_coverage_df <- sim_locus_coverage_df %>%
#     bind_rows(read_tsv(skesa_locus_file) %>% mutate(Isolate_ID = mut_id,Assembler = "SKESA")) %>%
#     bind_rows(read_tsv(spades_locus_file) %>% mutate(Isolate_ID = mut_id,Assembler = "SPAdes"))
# }

sim_locus_coverage_df <- read_tsv("data/Simulation/Simulated_Locus_Coverage.tsv")

# Compare SNP distances from Agona based on Screening vs. SNP Pipeline
sim_agona_distance_df <- sim_snp_distance_df %>%
  filter(Query_2 == "Agona") %>%
  rowwise() %>%
  mutate(Dataset = unlist(str_split(Query_1,"_"))[1]) %>% 
  ungroup() %>%
  select(Isolate_ID,Dataset,Assembler,CSP2_SNP_Distance = SNP_Distance) %>%
  left_join(bind_rows(spades_screen_screening_df,skesa_screen_screening_df) %>% select(Isolate_ID,Dataset,Assembler,CSP2_Screen_SNPs))

sim_self_distance_df <- sim_snp_distance_df %>%
  filter(Query_2 != "Agona") %>% 
  rowwise() %>%
  mutate(Dataset_1 = unlist(str_split(Query_1,"_"))[1],
         Dataset_2 = unlist(str_split(Query_2,"_"))[1]) %>%
  ungroup()

# Compare raw SNPs to detected SNPs
valid_snp_df <- sim_raw_snp_df %>%
  select(Isolate_ID,Ref_Pos,Raw_Ref_Base = Ref_Base,Raw_Query_Base = Query_Base,Raw_Cat = Cat)

datasets <- c("20X","40X","60X","80X",'100X')
sim_full_snp_df <- tibble()
for(dataset in datasets){
  for(assembler in c("SKESA","SPAdes")){
    sim_full_snp_df <- bind_rows(sim_full_snp_df,
                                 valid_snp_df %>% mutate(Dataset = dataset, Assembler = assembler)) %>%
      select(Isolate_ID,Dataset,Assembler,Ref_Pos,Raw_Ref_Base,Raw_Query_Base,Raw_Cat)
  }
}

sim_full_snp_df <- sim_full_snp_df %>%
  full_join(skesa_screen_snp_df %>%
              bind_rows(spades_screen_snp_df) %>%
              select(Isolate_ID,Dataset,Assembler,Ref_Pos,Dist_to_Query_End,Perc_Iden,Ref_Aligned,Ref_Base,Query_Base,Cat),by=c("Isolate_ID","Dataset","Assembler","Ref_Pos"))

sim_valid_snp_df <- sim_full_snp_df %>%
  filter(!is.na(Raw_Cat)) %>% 
  mutate(Cat = ifelse(is.na(Cat),"Missing",Cat)) %>%
  mutate(Cat = ifelse(Cat == "Purged_Density" & Raw_Cat != "Purged_Density","Incorrect_Density",Cat)) %>%
  mutate(Cat = ifelse(Raw_Cat == "Purged_Density" & Cat == "Purged_Density","Correct_Density",Cat)) %>%
  mutate(Dataset = factor(Dataset,levels = c("20X","40X","60X",'80X','100X')))


sim_valid_snp_isolate_df <- sim_valid_snp_df %>%
  group_by(Isolate_ID,Assembler,Dataset) %>%
  count(Cat) %>%
  pivot_wider(names_from = Cat, values_from = n) %>%
  ungroup() %>%
  mutate_if(is.numeric, ~ coalesce(., 0)) %>%
  mutate(SNPs_Covered = 300 - Missing,
         Accurate_SNPs = SNP + Correct_Density,
         SNPs_Purged = Incorrect_Density + Purged_Identity + Purged_Length) %>%
  select(Isolate_ID,Assembler,Dataset,SNPs_Covered,Accurate_SNPs,SNPs_Purged,Filtered_Query_Edge,Missing,Incorrect_Density,Purged_Identity,Purged_Length) %>%
  pivot_longer(cols = c('SNPs_Covered','Accurate_SNPs','SNPs_Purged','Filtered_Query_Edge','Missing','Incorrect_Density','Purged_Identity','Purged_Length'), names_to = 'Cat', values_to = 'Count')

# False positives
false_positive_df <- sim_full_snp_df %>%
  filter((is.na(Raw_Cat) | Raw_Cat == "Purged_Density") & (!(Cat %in% c("Insertion",'Deletion')) & !is.na(Cat) & !(Raw_Cat == "Purged_Density" & Cat == "Purged_Density")))

classic_false_positive_df <- sim_full_snp_df %>%
  filter(is.na(Raw_Cat) & Cat == "SNP")

classic_purged_false_positive_df <- sim_full_snp_df %>%
  filter(is.na(Raw_Cat) & Cat != "SNP" & !Cat %in% c("Insertion",'Deletion'))

# Edge Filtering
all_edge_filter_df <- sim_full_snp_df %>% filter(Cat == "Filtered_Query_Edge") %>%
  mutate(Dataset = factor(Dataset,levels = c("20X","40X","60X",'80X','100X')))

rescued_snp_df <- sim_full_snp_df %>% filter(Cat == "SNP") %>% select(Isolate_ID,Assembler,Ref_Pos) %>% distinct() %>%
  left_join(all_edge_filter_df) %>%
  filter(Cat == "Filtered_Query_Edge") %>%
  mutate(Cat = 'SNP_Rescued') %>%
  select(colnames(sim_full_snp_df))

edge_purged_df <- all_edge_filter_df %>%
  select(-Cat) %>%
  left_join(rescued_snp_df %>% select(Isolate_ID,Assembler,Ref_Pos,Cat) %>% distinct()) %>%
  filter(is.na(Cat)) %>%
  mutate(Cat = 'Edge_Purged') %>%
  select(colnames(sim_full_snp_df))

edge_df <- bind_rows(rescued_snp_df,edge_purged_df) %>%
  select(colnames(sim_full_snp_df))

# Bin all SNPs into single dataframe
final_sim_snp_df <- bind_rows(sim_full_snp_df %>% filter(!Cat %in% c("Filtered_Query_Edge")),edge_df) %>%
  mutate(Cat = ifelse(is.na(Cat),"Missing",Cat)) %>%
  mutate(Cat = ifelse(Cat == "Purged_Density" & (!is.na(Raw_Cat) & !Raw_Cat %in% c("Purged_Density")),"Incorrect_Density",Cat)) %>%
  mutate(Cat = ifelse(Raw_Cat %in% c("Purged_Density") & Cat %in% c("Purged_Density"),"Correct_Density",Cat)) %>%
  mutate(Dataset = factor(Dataset,levels = c("20X","40X","60X",'80X','100X'))) %>%
  filter(!Cat %in% c('Deletion','Insertion')) %>%
  mutate(Final_Cat = ifelse(Cat == "Missing","Missing",
                            ifelse(Cat == "Correct_Density","True_Positive",
                                   ifelse(!is.na(Raw_Cat) & Cat %in% c("SNP"),"True_Positive",
                                          ifelse(!is.na(Raw_Cat) & Cat %in% c("SNP_Rescued"),"Rescued_Positive",
                                                 ifelse(!is.na(Raw_Cat) & Cat %in% c("Edge_Purged"),"Edge_Filtered_Negative",
                                                        ifelse(!is.na(Raw_Cat) & (!is.na(Cat) & !Cat %in% c("SNP","SNP_Rescued")),"False_Negative",
                                                               ifelse(is.na(Raw_Cat) & Cat %in% c("SNP","SNP_Rescued"),"Included_False_Positive","Purged_False_Positive"))))))))

final_sim_isolate_snp_df <- final_sim_snp_df %>% 
  group_by(Isolate_ID,Assembler,Dataset) %>% 
  count(Final_Cat) %>% 
  pivot_wider(names_from = Final_Cat, values_from = n) %>% 
  ungroup() %>% 
  mutate_if(is.numeric, ~ coalesce(., 0)) %>%
  mutate(Positive = True_Positive + Rescued_Positive)

final_isolate_distance_df <- final_sim_isolate_snp_df %>%
  mutate(Distance = True_Positive + Rescued_Positive + Included_False_Positive)

##### 02: Simulated Data Plots #####

assembly_plot <- cowplot::plot_grid(ggplot(sim_isolate_data %>% filter(Isolate_Type != "Reference"), aes(x = Dataset, color=Assembler,fill = Assembler,y=Contig_Count)) +
                                      geom_violin(alpha=0.3,position = position_dodge(width = 1)) +
                                      theme_minimal() +
                                      theme(legend.position = "none") +
                                      facet_wrap(~Dataset,scales = "free_x",nrow=1),
                                    ggplot(sim_isolate_data %>% filter(Isolate_Type != "Reference"), aes(x = Dataset, color=Assembler,fill = Assembler,y=Assembly_Bases)) +
                                      geom_violin(alpha=0.3,position = position_dodge(width = 1)) +
                                      theme_minimal() +
                                      facet_wrap(~Dataset,nrow=1,scales="free_x"))


sim_ref_coverage_plot <- sim_ref_coverage_df %>% 
  ggplot(aes(x=Dataset,y=Reference_Percent_Aligned,color=Assembler,fill=Assembler)) +
  geom_violin(alpha = 0.3,position = position_dodge(width = 1)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1,dodge.width = 1)) + 
  theme_minimal() +
  theme(legend.position = "none")

sim_query_coverage_plot <- sim_query_coverage_df %>% 
  ggplot(aes(x=Dataset,y=Query_Percent_Aligned,color=Assembler,fill=Assembler)) +
  geom_violin(alpha = 0.3,position = position_dodge(width = 1)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1,dodge.width = 1)) + 
  theme_minimal() +
  theme(legend.position = "none")

sim_agona_distance_plot_df <- sim_agona_distance_df %>%
  left_join(raw_snp_count_df) %>%
  pivot_longer(cols=c("CSP2_SNP_Distance","CSP2_Screen_SNPs"),names_to = "Method",values_to = "Distance") %>%
  mutate(Dataset = factor(Dataset,levels = c("100X","80X",'60X','40X','20X'))) %>%
  mutate(Corrected_Distance = Distance - Raw_SNP_Count)

sim_agona_distance_plot <- ggplot(sim_agona_distance_plot_df,aes(x=Dataset,y=Distance,color=Method,fill=Method)) +
  geom_violin(alpha = 0.3,position = position_dodge(width = 1)) +
  geom_jitter(alpha=0.2,size=1,color="black",position = position_jitterdodge(jitter.width = 0.1,dodge.width = 1)) +
  theme_minimal() +
  geom_hline(yintercept = 300,linetype = "dashed") +
  facet_wrap(~Assembler,nrow=1)

sim_agona_corrected_distance_plot <- ggplot(sim_agona_distance_plot_df,aes(x=Dataset,y=Corrected_Distance,color=Method,fill=Method)) +
  geom_violin(alpha = 0.3,position = position_dodge(width = 1)) +
  geom_jitter(alpha=0.2,size=1,color="black",position = position_jitterdodge(jitter.width = 0.1,dodge.width = 1)) +
  theme_minimal() +
  geom_hline(yintercept = 0,linetype = "dashed") +
  facet_wrap(~Assembler,nrow=1)

sim_self_distance_plot_df <- sim_self_distance_df %>%
  mutate(Min_Coverage = ifelse(Dataset_1 == "20X" | Dataset_2 == "20X","20X",
                               ifelse(Dataset_1 == "40X" | Dataset_2 == "40X","40X",
                                      ifelse(Dataset_1 == "60X" | Dataset_2 == "60X","60X","80X"))),
         Max_Coverage = ifelse(Dataset_1 == "100X" | Dataset_2 == "100X","100X",
                               ifelse(Dataset_1 == "80X" | Dataset_2 == "80X","80X",
                                      ifelse(Dataset_1 == "60X" | Dataset_2 == "60X","60X","40X")))) %>%
  mutate(Max_Coverage = factor(Max_Coverage,levels = c("40X","60X","80X","100X")),
         Min_Coverage = factor(Min_Coverage,levels = c("20X","40X","60X","80X"))) %>%
  mutate(Min_Max_Coverage = paste0(Min_Coverage,"-",Max_Coverage)) %>%
  mutate(Min_Max_Coverage = factor(Min_Max_Coverage,levels = c("20X-40X","20X-60X","20X-80X","20X-100X","40X-60X","40X-80X","40X-100X","60X-80X","60X-100X","80X-100X"))) %>%
  mutate(Plot_Group = ifelse(Min_Coverage == "20X","20X",ifelse(Min_Coverage == "40X","40X","60X_and_Above"))) %>%
  mutate(Plot_Group = factor(Plot_Group,levels = c("20X",'40X',"60X_and_Above")))

sim_self_distance_plot <- ggplot(sim_self_distance_plot_df,aes(x=Min_Max_Coverage,y=SNP_Distance,fill=Assembler)) +
  geom_jitter(size=1,position = position_jitterdodge(jitter.width = 0.05,dodge.width = 1)) +
  geom_boxplot(alpha = 0.3,position = position_dodge(width = 1),color="black",notch = TRUE,outlier.size = 1) +
  theme_minimal() +
  scale_y_continuous(breaks = seq(0,12,by=1)) +
  facet_wrap(~Plot_Group,scales="free_x",nrow=3)

final_sim_isolate_long_df <- final_sim_isolate_snp_df %>%
  pivot_longer(cols = c('Missing','Rescued_Positive','Edge_Filtered_Negative','True_Positive','False_Negative','Included_False_Positive','Purged_False_Positive'),names_to = "Category",values_to = "Count") %>%
  mutate(Dataset = factor(Dataset,levels = c("20X","40X","60X",'80X','100X'))) %>%
  mutate(Alpha = ifelse(Category == "Purged_False_Positive",0.3,1)) %>%
  mutate(Category = factor(Category,levels = rev(c("Rescued_Positive","True_Positive","Missing","Edge_Filtered_Negative","False_Negative","Included_False_Positive","Purged_False_Positive"))))


##### 03: Real Data Analysis #####

# Data setup
data_dir = 'data/Real_Data/'

plot_colors <- c(viridis(5)[1:4],"red") 
names(plot_colors) <- c("A","B","C","D","F")

# Read in reference genomes
ref_df <- bind_rows(read_tsv(paste0(data_dir,'CFSAN_Refs.tsv')) %>%
                      rename(Reference_ID = "CFSAN_Reference") %>%
                      mutate(Dataset = "CFSAN"),
                    read_tsv(paste0(data_dir,'CSP2_Refs.tsv')) %>%
                      separate_rows(CSP2_Refs, sep = ",") %>%
                      rename(Reference_ID = "CSP2_Refs") %>%
                      mutate(Dataset = "CSP2"))

cfsan_refs <- ref_df %>% filter(Dataset == "CFSAN") %>% pull(Reference_ID)
csp2_refs <- ref_df %>% filter(Dataset == "CSP2") %>% pull(Reference_ID)

# Read in assembly data

# isolate_df <- lapply(list.files(paste0(data_dir,'Isolate_Data'), full.names = TRUE), function(query) {
#   read_tsv(query,show_col_types = FALSE) %>%
#     mutate(Analysis_ID = str_replace(tools::file_path_sans_ext(basename(query)),"_Isolate_Data",""))}) %>%
#   bind_rows() %>%
#   mutate(Species = ifelse(str_detect(Analysis_ID,"Ec_"),"Ecoli",
#                                  ifelse(str_detect(Analysis_ID,"Lm_"),"Listeria",
#                                         ifelse(str_detect(Analysis_ID,"Se_"),"Salmonella",NA))))
# 

isolate_df <- read_tsv(paste0(data_dir,'Assembly_Stats.tsv')) %>%
  select(Species,Analysis_ID,Isolate_ID,Contig_Count,Assembly_Bases,N50,N90,L50,L90) %>%
  mutate(Isolate_Type = ifelse(Isolate_ID %in% cfsan_refs,"CFSAN_Reference",
                               ifelse(Isolate_ID %in% csp2_refs,"CSP2_Reference","Query")))

simple_id_df <- isolate_df %>%
  select(Species,Analysis_ID,Isolate_ID) %>%
  distinct()

reference_isolate_df <- isolate_df %>%
  filter(Isolate_Type %in% c("CFSAN_Reference",'CSP2_Reference')) %>%
  rename(Reference_ID = 'Isolate_ID') %>%
  select(Species,Analysis_ID,Reference_ID,Reference_Type = Isolate_Type)

# Contig stats
contig_df <- isolate_df %>%
  select(Species,Analysis_ID,Isolate_ID,Contig_Count)

contig_stat_df <- contig_df %>%
  group_by(Species) %>%
  reframe(Contig_Count_Q25 = quantile(Contig_Count,0.25),
          Contig_Count_Q50 = quantile(Contig_Count,0.50),
          Contig_Count_Q75 = quantile(Contig_Count,0.75),
          Contig_Count_Outlier_A = getIQROutlier(Contig_Count)[2],
          Contig_Count_Outlier_B = getIQROutlier(Contig_Count,mult=3)[2],
          Contig_Count_Outlier_C = getIQROutlier(Contig_Count,mult=6)[2])

isolate_contig_df <- contig_df %>% 
  left_join(contig_stat_df) %>%
  group_by(Species) %>%
  reframe(Analysis_ID = Analysis_ID,
          Isolate_ID = Isolate_ID,
          Variable = "Contig_Count",
          Value = Contig_Count,
          Category = ifelse(Contig_Count >= Contig_Count_Outlier_C,"F3",
                            ifelse(Contig_Count >= Contig_Count_Outlier_B,"F2",
                                   ifelse(Contig_Count >= Contig_Count_Outlier_A,"F1",
                                          ifelse(Contig_Count <= Contig_Count_Q25,"A",
                                                 ifelse(Contig_Count <= Contig_Count_Q50,"B",
                                                        ifelse(Contig_Count <= Contig_Count_Q75,"C","D"))))))) %>%
  mutate(Category = factor(Category,levels = c("A","B","C","D","F1",'F2','F3'))) %>%
  group_by(Species) %>%
  mutate(ZScore = modified_z_score_test(Value)) %>%
  ungroup()

fail_contig_isolates <- isolate_contig_df %>% filter(Category %in% c("F1","F2","F3")) %>% pull(Isolate_ID)
fail_contig_isolate_count <- isolate_df %>% filter(Isolate_ID %in% fail_contig_isolates) %>% count(Species,name="Count")

# N50 stats
n50_df <- isolate_df %>% 
  filter(!Isolate_ID %in% fail_contig_isolates) %>%
  select(Species,Analysis_ID,Isolate_ID,N50)

n50_stat_df <- n50_df %>%
  group_by(Species) %>%
  reframe(N50_Q25 = quantile(N50,0.25),
          N50_Q50 = quantile(N50,0.50),
          N50_Q75 = quantile(N50,0.75),
          N50_Outlier_A = getIQROutlier(N50)[1],
          N50_Outlier_B = getIQROutlier(N50,mult=3)[1],
          N50_Outlier_C = getIQROutlier(N50,mult=6)[1])

isolate_n50_df <- isolate_df %>% 
  select(Species,Analysis_ID,Isolate_ID,N50) %>%
  left_join(n50_stat_df) %>%
  group_by(Species) %>%
  reframe(Analysis_ID = Analysis_ID,
          Isolate_ID = Isolate_ID,
          Variable = "N50",
          Value = N50,
          Category = ifelse(N50 <= N50_Outlier_C,"F3",
                            ifelse(N50 <= N50_Outlier_B,"F2",
                                   ifelse(N50 <= N50_Outlier_A,"F1",
                                          ifelse(N50 >= N50_Q75,"A",
                                                 ifelse(N50 >= N50_Q50,"B",
                                                        ifelse(N50 >= N50_Q25,"C","D"))))))) %>%
  mutate(Category = factor(Category,levels = c("A","B","C","D","F1",'F2','F3'))) %>%
  group_by(Species) %>%
  mutate(ZScore = modified_z_score_test(Value)) %>%
  ungroup()

# Read in CFSAN metrics
cfsan_metrics_df <- lapply(list.files(paste0(data_dir,'CFSAN_Metrics'), full.names = TRUE), function(query) {
  read_tsv(query,show_col_types = FALSE) %>%
    mutate(Analysis_ID = str_replace(tools::file_path_sans_ext(basename(query)),"_metrics",""))}) %>%
  bind_rows() %>%
  mutate(Species = ifelse(str_detect(Analysis_ID,"Ec_"),"Ecoli",
                          ifelse(str_detect(Analysis_ID,"Lm_"),"Listeria",
                                 ifelse(str_detect(Analysis_ID,"Se_"),"Salmonella",NA)))) %>%
  rename(Isolate_ID = Sample,
         Percent_Mapped = Percent_of_Reads_Mapped,
         Pileup_Depth = Average_Pileup_Depth,
         CFSAN_Missing = Missing_Preserved_SNP_Matrix_Positions) %>%
  select(Species,Analysis_ID,Isolate_ID,Percent_Mapped,Pileup_Depth,CFSAN_Missing)

percent_mapped_stat_df <- cfsan_metrics_df %>%
  group_by(Species) %>%
  reframe(Percent_Mapped_Q25 = quantile(Percent_Mapped,0.25),
          Percent_Mapped_Q50 = quantile(Percent_Mapped,0.50),
          Percent_Mapped_Q75 = quantile(Percent_Mapped,0.75),
          Percent_Mapped_Outlier_A = getIQROutlier(Percent_Mapped)[1],
          Percent_Mapped_Outlier_B = getIQROutlier(Percent_Mapped,mult=3)[1],
          Percent_Mapped_Outlier_C = getIQROutlier(Percent_Mapped,mult=6)[1])

isolate_percent_mapped_df <- cfsan_metrics_df %>% 
  select(Species,Analysis_ID,Isolate_ID,Percent_Mapped) %>%
  left_join(percent_mapped_stat_df) %>%
  group_by(Species) %>%
  reframe(Analysis_ID = Analysis_ID,
          Isolate_ID = Isolate_ID,
          Variable = "Percent_Mapped",
          Value = Percent_Mapped,
          Category = ifelse(Percent_Mapped <= Percent_Mapped_Outlier_C,"F3",
                            ifelse(Percent_Mapped <= Percent_Mapped_Outlier_B,"F2",
                                   ifelse(Percent_Mapped <= Percent_Mapped_Outlier_A,"F1",
                                          ifelse(Percent_Mapped >= Percent_Mapped_Q75,"A",
                                                 ifelse(Percent_Mapped >= Percent_Mapped_Q50,"B",
                                                        ifelse(Percent_Mapped >= Percent_Mapped_Q25,"C","D"))))))) %>%
  mutate(Category = factor(Category,levels = c("A","B","C","D","F1",'F2','F3'))) %>%
  group_by(Species) %>%
  mutate(ZScore = modified_z_score_test(Value)) %>%
  ungroup()

pileup_stat_df <- cfsan_metrics_df %>%
  group_by(Species) %>%
  reframe(Pileup_Depth_Q25 = quantile(Pileup_Depth,0.25),
          Pileup_Depth_Q50 = quantile(Pileup_Depth,0.50),
          Pileup_Depth_Q75 = quantile(Pileup_Depth,0.75),
          Pileup_Depth_Outlier_A = getIQROutlier(Pileup_Depth)[1],
          Pileup_Depth_Outlier_B = getIQROutlier(Pileup_Depth,mult=3)[1],
          Pileup_Depth_Outlier_C = getIQROutlier(Pileup_Depth,mult=6)[1])

isolate_pileup_depth_df <- cfsan_metrics_df %>% 
  select(Species,Analysis_ID,Isolate_ID,Pileup_Depth) %>%
  left_join(pileup_stat_df) %>%
  group_by(Species) %>%
  reframe(Analysis_ID = Analysis_ID,
          Isolate_ID = Isolate_ID,
          Variable = "Pileup_Depth",
          Value = Pileup_Depth,
          Category = ifelse(Pileup_Depth <= Pileup_Depth_Outlier_C,"F3",
                            ifelse(Pileup_Depth <= Pileup_Depth_Outlier_B,"F2",
                                   ifelse(Pileup_Depth <= Pileup_Depth_Outlier_A,"F1",
                                          ifelse(Pileup_Depth >= Pileup_Depth_Q75,"A",
                                                 ifelse(Pileup_Depth >= Pileup_Depth_Q50,"B",
                                                        ifelse(Pileup_Depth >= Pileup_Depth_Q25,"C","D"))))))) %>%
  mutate(Category = factor(Category,levels = c("A","B","C","D","F1",'F2','F3'))) %>%
  group_by(Species) %>%
  mutate(ZScore = modified_z_score_test(Value)) %>%
  ungroup()

cfsan_missing_stat_df <- cfsan_metrics_df %>%
  group_by(Species) %>%
  reframe(CFSAN_Missing_Q25 = quantile(CFSAN_Missing,0.25),
          CFSAN_Missing_Q50 = quantile(CFSAN_Missing,0.50),
          CFSAN_Missing_Q75 = quantile(CFSAN_Missing,0.75),
          CFSAN_Missing_Outlier_A = getIQROutlier(CFSAN_Missing)[2],
          CFSAN_Missing_Outlier_B = getIQROutlier(CFSAN_Missing,mult=3)[2],
          CFSAN_Missing_Outlier_C = getIQROutlier(CFSAN_Missing,mult=6)[2])

isolate_cfsan_missing_df <- cfsan_metrics_df %>% 
  select(Species,Analysis_ID,Isolate_ID,CFSAN_Missing) %>%
  left_join(cfsan_missing_stat_df) %>%
  group_by(Species) %>%
  reframe(Analysis_ID = Analysis_ID,
          Isolate_ID = Isolate_ID,
          Variable = "CFSAN_Missing",
          Value = CFSAN_Missing,
          Category = ifelse(CFSAN_Missing >= CFSAN_Missing_Outlier_C,"F3",
                            ifelse(CFSAN_Missing >= CFSAN_Missing_Outlier_B,"F2",
                                   ifelse(CFSAN_Missing >= CFSAN_Missing_Outlier_A,"F1",
                                          ifelse(CFSAN_Missing <= CFSAN_Missing_Q25,"A",
                                                 ifelse(CFSAN_Missing <= CFSAN_Missing_Q50,"B",
                                                        ifelse(CFSAN_Missing <= CFSAN_Missing_Q75,"C","D"))))))) %>%
  mutate(Category = factor(Category,levels = c("A","B","C","D","F1",'F2','F3'))) %>%
  group_by(Species) %>%
  mutate(ZScore = modified_z_score_test(Value)) %>%
  ungroup()

# Reference screening data

# csp2_ref_screening_df <- lapply(list.files(paste0(data_dir,'Reference_Screening'), full.names = TRUE), function(query) {
#   read_tsv(query,show_col_types = FALSE) %>%
#     mutate(Reference_ID = str_replace(tools::file_path_sans_ext(basename(query)),"_Reference_Screening",""))}) %>%
#   bind_rows()
# 
# csp2_ref_screening_df %>% 
#   left_join(id_df %>% select(Reference_ID = "Isolate_ID",Species,Analysis_ID)) %>%
#   mutate(Comparison = apply(.[, c("Reference_ID", "Query_ID")], 1, function(x) paste(sort(x), collapse = ";"))) %>%
#   select(Species,Analysis_ID,Reference_ID,Query_ID,Comparison,Screen_Category,CSP2_Screen_SNPs,Query_Percent_Aligned,Reference_Percent_Aligned,Query_Contigs,Query_Bases,Reference_Contigs,Reference_Bases,Raw_SNPs,Purged_Length,Purged_Identity,Purged_LengthIdentity,Purged_Invalid,Purged_Indel,Purged_Duplicate,Purged_Het,Purged_Density,Filtered_Query_Edge,Filtered_Ref_Edge,Filtered_Both_Edge,Kmer_Similarity,Shared_Kmers,Query_Unique_Kmers,Reference_Unique_Kmers,MUMmer_gSNPs,MUMmer_gIndels)  %>%

csp2_ref_screening_df <- read_tsv(paste0(data_dir,'CSP2_Reference_Screening.tsv'))
fail_contig_screening_df <- csp2_ref_screening_df %>% filter(Query_ID %in% fail_contig_isolates)

csp2_fail_screening_df <- csp2_ref_screening_df %>% filter(Screen_Category != "Pass")
csp2_fail_species_count <- isolate_df %>% filter(Isolate_ID %in% csp2_fail_screening_df$Query_ID) %>% distinct() %>% count(Species)

# Reference coverage
refcov_stats_df <- csp2_ref_screening_df %>% 
  filter(!Query_ID %in% fail_contig_isolates) %>%
  group_by(Species) %>%
  reframe(Reference_Percent_Aligned_Q25 = quantile(Reference_Percent_Aligned,0.25),
          Reference_Percent_Aligned_Q50 = quantile(Reference_Percent_Aligned,0.50),
          Reference_Percent_Aligned_Q75 = quantile(Reference_Percent_Aligned,0.75),
          Reference_Percent_Aligned_Outlier_A = getIQROutlier(Reference_Percent_Aligned)[1],
          Reference_Percent_Aligned_Outlier_B = getIQROutlier(Reference_Percent_Aligned,mult=3)[1],
          Reference_Percent_Aligned_Outlier_C = getIQROutlier(Reference_Percent_Aligned,mult=6)[1])

comparison_refcov_df <- csp2_ref_screening_df %>%
  left_join(refcov_stats_df) %>%
  rowwise() %>%
  reframe(Species = Species,
          Analysis_ID = Analysis_ID,
          Reference_ID = Reference_ID,
          Seq_1 = Query_ID,
          Seq_2 = Query_ID,
          Variable = "Reference_Coverage",
          Value = Reference_Percent_Aligned,
          Category = ifelse(Reference_Percent_Aligned <= Reference_Percent_Aligned_Outlier_C,"F3",
                            ifelse(Reference_Percent_Aligned <= Reference_Percent_Aligned_Outlier_B,"F2",
                                   ifelse(Reference_Percent_Aligned <= Reference_Percent_Aligned_Outlier_A,"F1",
                                          ifelse(Reference_Percent_Aligned >= Reference_Percent_Aligned_Q75,"A",
                                                 ifelse(Reference_Percent_Aligned >= Reference_Percent_Aligned_Q50,"B",
                                                        ifelse(Reference_Percent_Aligned >= Reference_Percent_Aligned_Q25,"C","D"))))))) %>%
  mutate(Category = factor(Category,levels = c("A","B","C","D","F1",'F2','F3'))) %>%
  group_by(Species) %>%
  mutate(ZScore = modified_z_score_test(Value)) %>%
  ungroup()

# Query coverage
querycov_stats_df <- csp2_ref_screening_df %>%
  filter(!Query_ID %in% fail_contig_isolates) %>%
  group_by(Species) %>%
  reframe(Query_Percent_Aligned_Q25 = quantile(Query_Percent_Aligned,0.25),
          Query_Percent_Aligned_Q50 = quantile(Query_Percent_Aligned,0.50),
          Query_Percent_Aligned_Q75 = quantile(Query_Percent_Aligned,0.75),
          Query_Percent_Aligned_Outlier_A = getIQROutlier(Query_Percent_Aligned)[1],
          Query_Percent_Aligned_Outlier_B = getIQROutlier(Query_Percent_Aligned,mult=3)[1],
          Query_Percent_Aligned_Outlier_C = getIQROutlier(Query_Percent_Aligned,mult=6)[1])

comparison_querycov_df <- csp2_ref_screening_df %>%
  left_join(querycov_stats_df) %>%
  rowwise() %>%
  reframe(Species = Species,
          Analysis_ID = Analysis_ID,
          Reference_ID = Reference_ID,
          Seq_1 = Query_ID,
          Seq_2 = Query_ID,
          Variable = "Query_Coverage",
          Value = Query_Percent_Aligned,
          Category = ifelse(Query_Percent_Aligned <= Query_Percent_Aligned_Outlier_C,"F3",
                            ifelse(Query_Percent_Aligned <= Query_Percent_Aligned_Outlier_B,"F2",
                                   ifelse(Query_Percent_Aligned <= Query_Percent_Aligned_Outlier_A,"F1",
                                          ifelse(Query_Percent_Aligned >= Query_Percent_Aligned_Q75,"A",
                                                 ifelse(Query_Percent_Aligned >= Query_Percent_Aligned_Q50,"B",
                                                        ifelse(Query_Percent_Aligned >= Query_Percent_Aligned_Q25,"C","D"))))))) %>%
  mutate(Category = factor(Category,levels = c("A","B","C","D","F1",'F2','F3'))) %>%
  group_by(Species) %>%
  mutate(ZScore = modified_z_score_test(Value)) %>%
  ungroup()

# Kmer similarity
kmer_similarity_stats_df <- csp2_ref_screening_df %>%
  filter(!Query_ID %in% fail_contig_isolates) %>%
  group_by(Species) %>%
  reframe(Kmer_Similarity_Q25 = quantile(Kmer_Similarity,0.25),
          Kmer_Similarity_Q50 = quantile(Kmer_Similarity,0.50),
          Kmer_Similarity_Q75 = quantile(Kmer_Similarity,0.75),
          Kmer_Similarity_Outlier_A = getIQROutlier(Kmer_Similarity)[1],
          Kmer_Similarity_Outlier_B = getIQROutlier(Kmer_Similarity,mult=3)[1],
          Kmer_Similarity_Outlier_C = getIQROutlier(Kmer_Similarity,mult=6)[1])


comparison_kmer_similarity_df <- csp2_ref_screening_df %>%
  left_join(kmer_similarity_stats_df) %>%
  rowwise() %>%
  reframe(Species = Species,
          Analysis_ID = Analysis_ID,
          Reference_ID = Reference_ID,
          Seq_1 = Query_ID,
          Seq_2 = Query_ID,
          Variable = "Kmer_Similarity",
          Value = Query_Percent_Aligned,
          Category = ifelse(Kmer_Similarity <= Kmer_Similarity_Outlier_C,"F3",
                            ifelse(Kmer_Similarity <= Kmer_Similarity_Outlier_B,"F2",
                                   ifelse(Kmer_Similarity <= Kmer_Similarity_Outlier_A,"F1",
                                          ifelse(Kmer_Similarity >= Kmer_Similarity_Q75,"A",
                                                 ifelse(Kmer_Similarity >= Kmer_Similarity_Q50,"B",
                                                        ifelse(Kmer_Similarity >= Kmer_Similarity_Q25,"C","D"))))))) %>%
  mutate(Category = factor(Category,levels = c("A","B","C","D","F1",'F2','F3'))) %>%
  group_by(Species) %>%
  mutate(ZScore = modified_z_score_test(Value)) %>%
  ungroup()

alignment_quality_df <- bind_rows(comparison_refcov_df,comparison_querycov_df,comparison_kmer_similarity_df) %>%
  select(-ZScore) %>%
  group_by(Species,Analysis_ID,Reference_ID) %>%
  pivot_wider(names_from = Variable,values_from = c("Value","Category")) %>%
  select(Species,Analysis_ID,Reference_ID,Query_ID = Seq_1,
         RefCov = Value_Reference_Coverage,
         QueryCov = Value_Query_Coverage,
         Kmer = Value_Kmer_Similarity,
         Ref_Cat = Category_Reference_Coverage,
         Query_Cat = Category_Query_Coverage,
         Kmer_Cat = Category_Kmer_Similarity)

reference_query_count_df <- csp2_ref_screening_df %>%
  group_by(Species,Analysis_ID,Reference_ID) %>%
  summarize(Query_Count = n())

# Read in CSP2 SNP data

# csp2_raw_pairwise_df <- lapply(list.files(paste0(data_dir,'CSP2_Raw'), full.names = TRUE), function(query) {
#   read_tsv(query,show_col_types = FALSE) %>%
#     mutate(Reference_ID = str_replace(tools::file_path_sans_ext(basename(query)),"_snp_distance_pairwise",""))}) %>%
#   bind_rows() %>%
#   mutate(Comparison = apply(.[, c("Query_1", "Query_2")], 1, function(x) paste(sort(x), collapse = ";"))) %>%
#   rename(CSP2_SNP_Raw = 'SNP_Distance',CSP2_Cocalled_Raw = 'SNPs_Cocalled') %>%
#   select(-Query_1,-Query_2)
# 
# csp2_preserved_pairwise_df <- lapply(list.files(paste0(data_dir,'CSP2_Preserved'), full.names = TRUE), function(query) {
#   read_tsv(query,show_col_types = FALSE) %>%
#     mutate(Reference_ID = str_replace(tools::file_path_sans_ext(basename(query)),"_snp_distance_pairwise_preserved",""))}) %>%
#   bind_rows() %>%
#   mutate(Comparison = apply(.[, c("Query_1", "Query_2")], 1, function(x) paste(sort(x), collapse = ";"))) %>%
#   rename(CSP2_SNP_Preserved = 'SNP_Distance',CSP2_Cocalled_Preserved = 'SNPs_Cocalled') %>%
#   select(-Query_1,-Query_2)
# 
# csp2_pairwise_df <- left_join(csp2_raw_pairwise_df,
#                               csp2_preserved_pairwise_df) %>%
#   left_join(id_df %>% select(Reference_ID = "Isolate_ID",Species,Analysis_ID)) %>%
#   select(Species,Analysis_ID,Reference_ID,Comparison,CSP2_SNP_Raw,CSP2_SNP_Preserved,CSP2_Cocalled_Raw,CSP2_Cocalled_Preserved)
# 
csp2_pairwise_df <- read_tsv(paste0(data_dir,'CSP2_Distances.tsv'))

# csp2_query_coverage_df <- lapply(list.files(paste0(data_dir,'Query_Coverage'), full.names = TRUE), function(query) {
#   read_tsv(query,show_col_types = FALSE) %>%
#     mutate(Reference_ID = str_replace(tools::file_path_sans_ext(basename(query)),"_Query_Coverage",""))}) %>%
#   bind_rows()
# 
# csp2_query_coverage_df %>%
#   left_join(id_df %>% select(Reference_ID = "Isolate_ID",Species,Analysis_ID)) %>%
#   select(Species,Analysis_ID,Reference_ID,Query_ID,SNP,Ref_Base,Percent_Missing,Purged,Uncovered,Rescued_SNP,Purged_Ref_Edge,Purged_Length,Purged_Identity,Purged_Invalid,Purged_Indel,Purged_Heterozygous,Purged_Density) %>%

csp2_query_coverage_df <- read_tsv(paste0(data_dir,'CSP2_Query_Coverage.tsv'))

# csp2_locus_category_df <- lapply(list.files(paste0(data_dir,'Locus_Categories'), full.names = TRUE), function(query) {
#   read_tsv(query,show_col_types = FALSE) %>%
#     mutate(Reference_ID = str_replace(tools::file_path_sans_ext(basename(query)),"_Locus_Categories",""))}) %>%
#   bind_rows()
# 
# csp2_locus_category_df %>%
#   left_join(id_df %>% select(Reference_ID = "Isolate_ID",Species,Analysis_ID)) %>%
#   select(Species,Analysis_ID,Reference_ID,Ref_Loc,SNP_Count,Ref_Base_Count,Uncovered_Count,Purged_Count,Missing_Ratio) %>%

csp2_locus_category_df <- read_tsv(paste0(data_dir,'CSP2_Locus_Categories.tsv'))

# Read in NCBI/CFSAN data

# cfsan_raw_pairwise_df <- lapply(list.files(paste0(data_dir,'CFSAN_Raw'), full.names = TRUE), function(query) {
#   read_tsv(query,show_col_types = FALSE) %>%
#     mutate(Analysis_ID = str_replace(tools::file_path_sans_ext(basename(query)),"_snp_distance_pairwise",""))}) %>%
#   bind_rows() %>%
#   filter(Seq1 != Seq2) %>%
#   rename(CFSAN_Raw = "Distance") %>%
#   mutate(Comparison = apply(.[, c("Seq1", "Seq2")], 1, function(x) paste(sort(x), collapse = ";"))) %>%
#   select(-Seq1,-Seq2)
# 
# cfsan_preserved_pairwise_df <- lapply(list.files(paste0(data_dir,'CFSAN_Preserved'), full.names = TRUE), function(query) {
#   read_tsv(query,show_col_types = FALSE) %>%
#     mutate(Analysis_ID = str_replace(tools::file_path_sans_ext(basename(query)),"_snp_distance_pairwise_preserved",""))}) %>%
#   bind_rows() %>%
#   filter(Seq1 != Seq2) %>%
#   rename(CFSAN_Preserved = "Distance") %>%
#   mutate(Comparison = apply(.[, c("Seq1", "Seq2")], 1, function(x) paste(sort(x), collapse = ";"))) %>%
#   select(-Seq1,-Seq2)
# 
# cfsan_raw_pairwise_df %>% group_by(Comparison) %>% summarize(Count = n(),Distance = n_distinct(CFSAN_Raw)) %>% filter(Distance > 1)
# cfsan_preserved_pairwise_df %>% group_by(Comparison) %>% summarize(Count = n(),Distance = n_distinct(CFSAN_Preserved)) %>% filter(Distance > 1)
# 
# cfsan_distance_df <- left_join(cfsan_raw_pairwise_df %>%
#                                  distinct(),
#                                cfsan_preserved_pairwise_df,multiple = "first") %>%
#   mutate(Species = ifelse(str_detect(Analysis_ID,"Ec_"),"Ecoli",
#                           ifelse(str_detect(Analysis_ID,"Lm_"),"Listeria",
#                                  ifelse(str_detect(Analysis_ID,"Se_"),"Salmonella",NA)))) %>%
#   select(Species,Analysis_ID,Comparison,CFSAN_Raw,CFSAN_Preserved)
# 

cfsan_distance_df <- read_tsv(paste0(data_dir,'CFSAN_Distances.tsv'))

ecoli_clusters <- isolate_df %>% filter(Species == "Ecoli") %>% pull(Analysis_ID) %>% str_replace("Ec_","") %>% unique()
listeria_clusters <- isolate_df %>% filter(Species == "Listeria") %>% pull(Analysis_ID) %>% str_replace("Lm_","") %>% unique()
salmonella_clusters <- isolate_df %>% filter(Species == "Salmonella") %>% pull(Analysis_ID) %>% str_replace("Se_","") %>% unique()

# ecoli_distances <- arrow::read_parquet(paste0(data_dir,'NCBI/Ecoli_SNP_Distances.parquet')) %>%
#   filter(PDS_acc %in% ecoli_clusters) %>%
#   filter(SRR_1 %in% all_isolates & SRR_2 %in% all_isolates) %>%
#   mutate(Comparison = apply(.[, c("SRR_1", "SRR_2")], 1, function(x) paste(sort(x), collapse = ";"))) %>%
#   select(Comparison,NCBI_Distance = delta_positions_unambiguous)
# 
# listeria_distances <- arrow::read_parquet(paste0(data_dir,'NCBI/Listeria_SNP_Distances.parquet')) %>%
#   filter(PDS_acc %in% listeria_clusters) %>%
#   filter(SRR_1 %in% all_isolates & SRR_2 %in% all_isolates) %>%
#   mutate(Comparison = apply(.[, c("SRR_1", "SRR_2")], 1, function(x) paste(sort(x), collapse = ";"))) %>%
#   select(Comparison,NCBI_Distance = delta_positions_unambiguous)
#   
# salmonella_distances <- arrow::read_parquet(paste0(data_dir,'NCBI/Salmonella_SNP_Distances.parquet')) %>%
#   filter(PDS_acc %in% salmonella_clusters) %>%
#   filter(SRR_1 %in% all_isolates & SRR_2 %in% all_isolates) %>%
#   mutate(Comparison = apply(.[, c("SRR_1", "SRR_2")], 1, function(x) paste(sort(x), collapse = ";"))) %>%
#   select(Comparison,NCBI_Distance = delta_positions_unambiguous)
#   
# ncbi_distances <- ecoli_distances %>%
#   bind_rows(listeria_distances) %>%
#   bind_rows(salmonella_distances)

csp2_distance_df <- csp2_pairwise_df %>%
  select(Species,Analysis_ID,Reference_ID,Comparison,CSP2 = CSP2_SNP_Preserved)

cfsan_ncbi_distance_df <- read_tsv(paste0(data_dir,'Single_Reference_Distances.tsv')) %>%
  separate(Comparison,c("Seq_1","Seq_2"),";",remove=FALSE) %>%
  select(Species,Analysis_ID,Comparison,Seq_1,Seq_2,CFSAN = CFSAN_Preserved,NCBI) %>%
  filter(Comparison %in% csp2_distance_df$Comparison)

full_distance_df <- csp2_distance_df %>%
  left_join(cfsan_ncbi_distance_df,by = c("Species","Analysis_ID","Comparison")) %>%
  mutate(Reference_Type = ifelse(Reference_ID %in% cfsan_refs,"CFSAN","CSP2"))

contig_outlier_distance_df <- full_distance_df %>%
  filter(Seq_1 %in% fail_contig_isolates | Seq_2 %in% fail_contig_isolates)

contig_inlier_distance_df <- full_distance_df %>%
  filter(!Seq_1 %in% fail_contig_isolates & !Seq_2 %in% fail_contig_isolates)

quality_distance_df <- full_distance_df %>%
  left_join(alignment_quality_df %>% select(Reference_ID,Seq_1 = Query_ID,Seq_1_Ref_Cat = Ref_Cat)) %>%
  left_join(alignment_quality_df %>% select(Reference_ID,Seq_2 = Query_ID,Seq_2_Ref_Cat = Ref_Cat)) %>%
  left_join(alignment_quality_df %>% select(Reference_ID,Seq_1 = Query_ID,Seq_1_Query_Cat = Query_Cat)) %>%
  left_join(alignment_quality_df %>% select(Reference_ID,Seq_2 = Query_ID,Seq_2_Query_Cat = Query_Cat)) %>%
  left_join(alignment_quality_df %>% select(Reference_ID,Seq_1 = Query_ID,Seq_1_Kmer_Cat = Kmer_Cat)) %>%
  left_join(alignment_quality_df %>% select(Reference_ID,Seq_2 = Query_ID,Seq_2_Kmer_Cat = Kmer_Cat)) %>%
  mutate(Seq_1_Ref_Cat = as.character(Seq_1_Ref_Cat),
         Seq_2_Ref_Cat = as.character(Seq_2_Ref_Cat),
         Seq_1_Query_Cat = as.character(Seq_1_Query_Cat),
         Seq_2_Query_Cat = as.character(Seq_2_Query_Cat),
         Seq_1_Kmer_Cat = as.character(Seq_1_Kmer_Cat),
         Seq_2_Kmer_Cat = as.character(Seq_2_Kmer_Cat),
         Seq_1_Ref_Cat = replace_na(Seq_1_Ref_Cat,"Reference"),
         Seq_2_Ref_Cat = replace_na(Seq_2_Ref_Cat,"Reference"),
         Seq_1_Query_Cat = replace_na(Seq_1_Query_Cat,"Reference"),
         Seq_2_Query_Cat = replace_na(Seq_2_Query_Cat,"Reference"),
         Seq_1_Kmer_Cat = replace_na(Seq_1_Kmer_Cat,"Reference"),
         Seq_2_Kmer_Cat = replace_na(Seq_2_Kmer_Cat,"Reference"))

coverage_distance_df <- quality_distance_df %>%
  mutate(Contig = ifelse(!Seq_1 %in% fail_contig_isolates & !Seq_2 %in% fail_contig_isolates,"A",
                         ifelse((Seq_1 %in% fail_contig_isolates | Seq_2 %in% fail_contig_isolates) & !(Seq_1 %in% fail_contig_isolates & Seq_2 %in% fail_contig_isolates),"B","C")),
         Reference = ifelse(!Seq_1_Ref_Cat %in% c("F1","F2","F3") & !Seq_2_Ref_Cat %in% c("F1","F2","F3"),"A",
                                ifelse((Seq_1_Ref_Cat %in% c("F1","F2","F3")| Seq_2_Ref_Cat %in% c("F1","F2","F3")) & !(Seq_1_Ref_Cat %in% c("F1","F2","F3") & Seq_1_Ref_Cat %in% c("F1","F2","F3")),"B","C")),
         Query = ifelse(!Seq_1_Query_Cat %in% c("F1","F2","F3") & !Seq_2_Query_Cat %in% c("F1","F2","F3"),"A",
                            ifelse((Seq_1_Query_Cat %in% c("F1","F2","F3")| Seq_2_Query_Cat %in% c("F1","F2","F3")) & !(Seq_1_Query_Cat %in% c("F1","F2","F3") & Seq_1_Query_Cat %in% c("F1","F2","F3")),"B","C")),
         Kmer = ifelse(!Seq_1_Kmer_Cat %in% c("F1","F2","F3") & !Seq_2_Kmer_Cat %in% c("F1","F2","F3"),"A",
                            ifelse((Seq_1_Kmer_Cat %in% c("F1","F2","F3")| Seq_2_Kmer_Cat %in% c("F1","F2","F3")) & !(Seq_1_Kmer_Cat %in% c("F1","F2","F3") & Seq_1_Kmer_Cat %in% c("F1","F2","F3")),"B","C"))) %>%
    select(Species,Analysis_ID,Reference_ID,Reference_Type,Comparison,Seq_1,Seq_2,CFSAN,NCBI,CSP2,Contig,Reference,Query,Kmer)

ncbi_cfsan_coverage_distance_df <- coverage_distance_df %>%
  filter(!Seq_1 %in% fail_contig_isolates & !Seq_2 %in% fail_contig_isolates) %>%
  select(Species,Analysis_ID,Comparison,Seq_1,Seq_2,CFSAN,NCBI) %>%
  distinct()

ncbi_cfsan_species_percent_df <- ncbi_cfsan_coverage_distance_df %>%
  mutate(CFSAN_NCBI = CFSAN - NCBI) %>%
  group_by(Species) %>%
  summarize(Percent_Within_One = 100*(sum(abs(CFSAN_NCBI) <= 1)/n()),
            Percent_Within_Three = 100*(sum(abs(CFSAN_NCBI) <= 3)/n()),
            Percent_Above_Three = 100*(sum(CFSAN_NCBI >= 4)/n()),
            Percent_Below_Three = 100*sum(CFSAN_NCBI <= -4)/n())

ncbi_cfsan_species_corr_df <- ncbi_cfsan_coverage_distance_df %>%
  group_by(Species) %>%
  summarize(
    Correlation = cor(CFSAN, NCBI),
    P_value = cor.test(CFSAN, NCBI)$p.value,
    CI_low = cor.test(CFSAN, NCBI)$conf.int[1],
    CI_high = cor.test(CFSAN, NCBI)$conf.int[2])

csp2_ncbi_species_all_percent_df <- coverage_distance_df %>%
  filter(!Seq_1 %in% fail_contig_isolates & !Seq_2 %in% fail_contig_isolates) %>%
  select(Species,CSP2,NCBI,Reference_Type) %>%
  mutate(CSP2_NCBI = CSP2 - NCBI) %>%
  group_by(Species) %>%
  group_by(Species,Reference_Type) %>%
  summarize(Percent_Within_One = 100*(sum(abs(CSP2_NCBI) <= 1)/n()),
            Percent_Within_Three = 100*(sum(abs(CSP2_NCBI) <= 3)/n()),
            Percent_Above_Three = 100*(sum(CSP2_NCBI >= 4)/n()),
            Percent_Below_Three = 100*sum(CSP2_NCBI <= -4)/n()) %>%
  mutate(Comparison = "CSP2_NCBI")

csp2_ncbi_species_all_corr_df <- coverage_distance_df %>%
  filter(!Seq_1 %in% fail_contig_isolates & !Seq_2 %in% fail_contig_isolates) %>%
  group_by(Species) %>%
  group_by(Species,Reference_Type) %>%
  summarize(
    Correlation = cor(CSP2, NCBI),
    P_value = cor.test(CSP2, NCBI)$p.value,
    CI_low = cor.test(CSP2, NCBI)$conf.int[1],
    CI_high = cor.test(CSP2, NCBI)$conf.int[2])

csp2_cfsan_species_all_percent_df <- coverage_distance_df %>%
  filter(!Seq_1 %in% fail_contig_isolates & !Seq_2 %in% fail_contig_isolates) %>%
  select(Species,CSP2,CFSAN,Reference_Type,Reference) %>%
  mutate(CSP2_CFSAN = CSP2 - CFSAN) %>%
  group_by(Species) %>%
  group_by(Species,Reference_Type) %>%
  summarize(Percent_Within_One = 100*(sum(abs(CSP2_CFSAN) <= 1)/n()),
            Percent_Within_Three = 100*(sum(abs(CSP2_CFSAN) <= 3)/n()),
            Percent_Above_Three = 100*(sum(CSP2_CFSAN >= 4)/n()),
            Percent_Below_Three = 100*sum(CSP2_CFSAN <= -4)/n()) %>%
  mutate(Comparison = "CSP2_CFSAN")

csp2_cfsan_species_all_corr_df <- coverage_distance_df %>%
  filter(!Seq_1 %in% fail_contig_isolates & !Seq_2 %in% fail_contig_isolates) %>%
  group_by(Species) %>%
  group_by(Species,Reference_Type) %>%
  summarize(
    Correlation = cor(CSP2, CFSAN),
    P_value = cor.test(CSP2, CFSAN)$p.value,
    CI_low = cor.test(CSP2, CFSAN)$conf.int[1],
    CI_high = cor.test(CSP2, CFSAN)$conf.int[2])

# Analysis/Reference Level
csp2_breakdown_distance_df <- coverage_distance_df %>%
  filter(!Seq_1 %in% fail_contig_isolates & !Seq_2 %in% fail_contig_isolates) %>%
  select(Species,Analysis_ID,Reference_ID,Reference_Type,Comparison,Seq_1,Seq_2,CSP2,CFSAN,NCBI,Reference,Query) %>%
  distinct() %>%
  rowwise() %>%
  mutate(Comparison_Type = paste0(Reference,Query)) %>%
  ungroup() %>%
  mutate(Distance_Type = ifelse(abs(CSP2 - CFSAN) >= 4 & abs(CSP2 - NCBI) >= 4  & abs(NCBI - CFSAN) >= 4,"Outside_Three",
                                ifelse(abs(CSP2 - CFSAN) <= 3 & abs(CSP2 - NCBI) <= 3 & abs(NCBI - CFSAN) <= 3,"Inside_Three",
                                       ifelse(
                                         (abs(CSP2 - CFSAN) <= 3 & abs(CSP2 - NCBI) <= 3) | 
                                         (abs(NCBI - CFSAN) <= 3 & abs(CSP2 - NCBI) <= 3) | 
                                         (abs(CSP2 - CFSAN) <= 3 & abs(CFSAN - NCBI) <= 3),"Daisy_Chain","Other"))))

# Stats
ncbi_cfsan_coverage_distance_df %>%
  mutate(Difference = CFSAN - NCBI) %>%
  group_by(Species) %>%
  summarize(Median = median(Difference,na.rm = TRUE),
            Mean = mean(Difference,na.rm = TRUE),
            SD = sd(Difference,na.rm = TRUE),
            N = n()) %>%
  mutate(SE = SD/sqrt(N),
         CI = map2_chr(Mean,SE,calculate_ci))

coverage_distance_df %>%
  filter(!Seq_1 %in% fail_contig_isolates & !Seq_2 %in% fail_contig_isolates) %>%
  mutate(Difference = CSP2 - CFSAN) %>%
  group_by(Species,Reference_Type) %>%
  summarize(Median = median(Difference,na.rm = TRUE),
            Mean = mean(Difference,na.rm = TRUE),
            SD = sd(Difference,na.rm = TRUE),
            N = n()) %>%
  mutate(SE = SD/sqrt(N),
         CI = map2_chr(Mean,SE,calculate_ci))

coverage_distance_df %>%
  filter(!Seq_1 %in% fail_contig_isolates & !Seq_2 %in% fail_contig_isolates) %>%
  mutate(Difference = CSP2 - NCBI) %>%
  group_by(Species,Reference_Type) %>%
  summarize(Median = median(Difference,na.rm = TRUE),
            Mean = mean(Difference,na.rm = TRUE),
            SD = sd(Difference,na.rm = TRUE),
            N = n()) %>%
  mutate(SE = SD/sqrt(N),
         CI = map2_chr(Mean,SE,calculate_ci))

##### O4: Manuscript Plots ########

fig2_df <- final_sim_isolate_snp_df %>%
  mutate(Dataset = factor(Dataset,levels = c("20X","40X","60X","80X",'100X'))) %>%
  mutate(Distance = Positive + Included_False_Positive) %>%
  mutate(Total_False_Negative = -1*(False_Negative + Missing + Edge_Filtered_Negative)) %>%
  select(Isolate_ID,Assembler,Dataset,Total_False_Negative,Included_False_Positive) %>%
  pivot_longer(cols = c(Total_False_Negative,Included_False_Positive),names_to = "SNP_Category",values_to = "Count") %>%
  mutate(SNP_Category = factor(SNP_Category,levels = c("Total_False_Negative","Included_False_Positive")))


fig2_rough <- ggplot(fig2_df,aes(x=Count,y=Dataset,fill=Assembler,color=Assembler)) +
  geom_density_ridges(alpha=0.65,scale=1)+
  theme_minimal() +
  theme(strip.text = element_text(size=16),
        axis.text.x = element_text(size = 16,face="bold"),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16,face="bold",margin = margin(r = 5)),
        panel.border = element_rect(color = "darkgray", fill = NA)) +
  facet_wrap(~SNP_Category,scales="free_x") +
  scale_color_manual(values = c("SKESA" = "blue","SPAdes" = "orange3")) +
  scale_fill_manual(values = c("SKESA" = "blue","SPAdes" = "orange3")) +
  scale_x_continuous(breaks = seq(min(fig2_df$Count), max(fig2_df$Count), by = 5))

fig2_rough

fig3_df <- coverage_distance_df %>%
  filter(!Seq_1 %in% fail_contig_isolates & !Seq_2 %in% fail_contig_isolates) %>%
  select(Species,Analysis_ID,Reference_ID,Reference_Type,Comparison,Seq_1,Seq_2,CSP2,CFSAN,NCBI,Reference,Query) %>%
  distinct()

cfsan_ncbi_slope_df <- fig3_df %>%
  select(Species,Analysis_ID,Comparison,CFSAN,NCBI) %>%
  distinct() %>%
  group_by(Analysis_ID) %>%
  summarize(Slope = lm(CFSAN ~ NCBI)$coefficients[2],
            Intercept = lm(CFSAN ~ NCBI)$coefficients[1],
            Correlation = cor(NCBI,CFSAN),
            Percent_Within_Three = sum(abs(NCBI - CFSAN) <= 3)/n()) %>%
  mutate(Dataset = "CFSAN_NCBI",Reference_Type = "NCBI")

csp2_cfsan_slope_df <- fig3_df %>%
  select(Species,Analysis_ID,Reference_Type,Comparison,CFSAN,CSP2) %>%
  distinct() %>%
  group_by(Analysis_ID,Reference_Type) %>%
  summarize(Slope = lm(CSP2 ~ CFSAN)$coefficients[2],
            Intercept = lm(CSP2 ~ CFSAN)$coefficients[1],
            Correlation = cor(CFSAN,CSP2),
            Percent_Within_Three = sum(abs(CSP2 - CFSAN) <= 3)/n()) %>%
  mutate(Dataset = "CSP2_CFSAN")

csp2_ncbi_slope_df <- fig3_df %>%
  select(Species,Analysis_ID,Reference_Type,Comparison,CSP2,NCBI) %>%
  distinct() %>%
  group_by(Analysis_ID,Reference_Type) %>%
  summarize(Slope = lm(CSP2 ~ NCBI)$coefficients[2],
            Intercept = lm(CSP2 ~ NCBI)$coefficients[1],
            Correlation = cor(NCBI,CSP2),
            Percent_Within_Three = sum(abs(CSP2 - NCBI) <= 3)/n()) %>%
  mutate(Dataset = "CSP2_NCBI")

cfsan_ncbi_diff_df <- fig3_df %>%
  select(Species,Analysis_ID,Comparison,CFSAN,NCBI) %>%
  distinct() %>%
  mutate(Difference = CFSAN - NCBI,
         Dataset = "CFSAN_NCBI",
         Reference_Type = "NCBI")

csp2_cfsan_diff_df <- fig3_df %>%
  select(Species,Analysis_ID,Reference_Type,Comparison,CFSAN,CSP2) %>%
  distinct() %>%
  mutate(Difference = CSP2 - CFSAN,
         Dataset = "CSP2_CFSAN")

csp2_ncbi_diff_df <- fig3_df %>%
  select(Species,Analysis_ID,Reference_Type,Comparison,CSP2,NCBI) %>%
  distinct() %>%
  mutate(Difference = CSP2 - NCBI,
         Dataset = "CSP2_NCBI")

fig3_df <- bind_rows(cfsan_ncbi_slope_df %>% select(names(csp2_cfsan_slope_df)),
                     csp2_cfsan_slope_df,
                     csp2_ncbi_slope_df) %>%
  left_join(isolate_df %>% select(Species,Analysis_ID),by="Analysis_ID",relationship = "many-to-many") %>%
  mutate(Dataset = as.character(Dataset)) %>%
  select(-Intercept) %>%
  pivot_longer(cols = c("Slope","Correlation","Percent_Within_Three"),names_to = "Measurement",values_to = "Value") %>%
  mutate(Measurement = ifelse(Measurement == "Percent_Within_Three","Percent Within 3 SNPs",Measurement)) %>%
  mutate(Species = ifelse(Species == "Ecoli","E. coli",Species)) %>%
  distinct()

fig3_box_df <- fig3_df %>% filter(Measurement %in% c("Correlation","Slope","Percent Within 3 SNPs")) %>%
  mutate(Measurement = factor(Measurement,levels = c("Percent Within 3 SNPs","Slope","Correlation")),
         Reference_Type <- factor(Reference_Type, levels = c("CSP2", "CFSAN", "NCBI")),
         Dataset <- factor(Dataset, levels = c("CSP2_CFSAN", "CFSAN_NCBI", "CSP2_NCBI")))

ecoli_plot <- ggplot(fig3_box_df %>% filter(Species == "E. coli") ,aes(x=Dataset,y=Value,fill=Reference_Type)) +
  geom_boxplot(color="black",outlier.alpha = 1,outlier.size = 3,outlier.color="black",outlier.shape = 21) +
  facet_grid(Species ~ Measurement) +
  scale_x_discrete(limits = c("CSP2_CFSAN", "CFSAN_NCBI", "CSP2_NCBI"),labels = c("CSP2 vs. CSP1", "CSP1 vs. NCBI", "CSP2 vs. NCBI")) +
  theme_minimal() +
  theme(legend.position = 'none',
        strip.text = element_text(size = 14,face="bold.italic"),
        axis.text.x = element_text(size = 18,face="bold"),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16,face="bold",margin = margin(r = 5)),
        panel.border =  element_rect(color = "black", fill = NA)) +
  scale_color_manual(values = c("CSP2" = "blue","CFSAN" = "orange3","NCBI"="gray45")) +
  scale_fill_manual(values = c("CSP2" = "blue","CFSAN" = "orange3","NCBI"="gray45")) +
  geom_hline(yintercept = 1,linetype = "dashed")

sal_plot <- ggplot(fig3_box_df %>% filter(Species == "Salmonella") ,aes(x=Dataset,y=Value,fill=Reference_Type)) +
  geom_boxplot(color="black",outlier.alpha = 1,outlier.size = 3,outlier.color="black",outlier.shape = 21) +
  facet_grid(Species ~ Measurement) +
  scale_x_discrete(limits = c("CSP2_CFSAN", "CFSAN_NCBI", "CSP2_NCBI"),labels = c("CSP2 vs. CSP1", "CSP1 vs. NCBI", "CSP2 vs. NCBI")) +
  theme_minimal() +
  theme(legend.position = 'none',
        strip.text = element_text(size = 14,face="bold.italic"),
        axis.text.x = element_text(size = 18,face="bold"),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16,face="bold",margin = margin(r = 5)),
        panel.border =  element_rect(color = "black", fill = NA)) +
  scale_color_manual(values = c("CSP2" = "blue","CFSAN" = "orange3","NCBI"="gray45")) +
  scale_fill_manual(values = c("CSP2" = "blue","CFSAN" = "orange3","NCBI"="gray45")) +
  geom_hline(yintercept = 1,linetype = "dashed")

list_plot <- ggplot(fig3_box_df %>% filter(Species == "Listeria") ,aes(x=Dataset,y=Value,fill=Reference_Type)) +
  geom_boxplot(color="black",outlier.alpha = 1,outlier.size = 3,outlier.color="black",outlier.shape = 21) +
  facet_grid(Species ~ Measurement) +
  scale_x_discrete(limits = c("CSP2_CFSAN", "CFSAN_NCBI", "CSP2_NCBI"),labels = c("CSP2 vs. CSP1", "CSP1 vs. NCBI", "CSP2 vs. NCBI")) +
  theme_minimal() +
  theme(legend.position = 'none',
        strip.text = element_text(size = 14,face="bold.italic"),
        axis.text.x = element_text(size = 18,face="bold"),
        axis.text.y = element_text(size = 18),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16,face="bold",margin = margin(r = 5)),
        panel.border =  element_rect(color = "black", fill = NA)) +
  scale_color_manual(values = c("CSP2" = "blue","CFSAN" = "orange3","NCBI"="gray45")) +
  scale_fill_manual(values = c("CSP2" = "blue","CFSAN" = "orange3","NCBI"="gray45")) +
  geom_hline(yintercept = 1,linetype = "dashed")

fig_s1_contig_non20 <- ggplot(sim_isolate_data %>% filter(Isolate_Type != "Reference" & Dataset != "20X"), aes(x = Dataset, color=Assembler,fill = Assembler,y=Contig_Count)) +
  geom_violin(position = position_dodge(width = 1)) +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("SKESA" = "blue","SPAdes" = "orange3")) +
  scale_fill_manual(values = c("SKESA" = "blue","SPAdes" = "orange3")) +
    theme(strip.text = element_text(size=16),
          axis.text.x = element_text(size = 16,face="bold"),
          axis.text.y = element_text(size = 14),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 16,face="bold",margin = margin(r = 5)),
          panel.border = element_rect(color = "darkgray", fill = NA)) +
  ylab("Contig Count")

fig_s1_contig_20 <- ggplot(sim_isolate_data %>% filter(Isolate_Type != "Reference" & Dataset == "20X"), aes(x = Dataset, color=Assembler,fill = Assembler,y=Contig_Count)) +
  geom_violin(position = position_dodge(width = 1)) +
  theme_minimal() +
  scale_color_manual(values = c("SKESA" = "blue","SPAdes" = "orange3")) +
  scale_fill_manual(values = c("SKESA" = "blue","SPAdes" = "orange3")) +
  theme(strip.text = element_text(size=16),
        axis.text.x = element_text(size = 16,face="bold"),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        panel.border = element_rect(color = "darkgray", fill = NA))
  
fig_s1_assembly <- ggplot(sim_isolate_data %>% filter(Isolate_Type != "Reference"), aes(x = Dataset, color=Assembler,fill = Assembler,y=Assembly_Bases)) +
  geom_violin(position = position_dodge(width = 1)) +
  theme_minimal() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("SKESA" = "blue","SPAdes" = "orange3")) +
  scale_fill_manual(values = c("SKESA" = "blue","SPAdes" = "orange3")) +
  theme(strip.text = element_text(size=16),
        axis.text.x = element_text(size = 16,face="bold"),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16,face="bold",margin = margin(r = 5)),
        panel.border = element_rect(color = "darkgray", fill = NA)) +
  ylab("Assembly Bases")

fig_s1_plot <- cowplot::plot_grid(cowplot::plot_grid(fig_s1_contig_non20,fig_s1_contig_20,rel_widths = c(2,1)),fig_s1_assembly,nrow=2)

fig_s2_ref <- sim_ref_coverage_df %>% 
  ggplot(aes(x=Dataset,y=Reference_Percent_Aligned,color=Assembler,fill=Assembler)) +
  geom_violin(alpha = 0.3,position = position_dodge(width = 1)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1,dodge.width = 1)) + 
  theme_minimal() +
  theme(legend.position = "none") +
  scale_color_manual(values = c("SKESA" = "blue","SPAdes" = "orange3")) +
  scale_fill_manual(values = c("SKESA" = "blue","SPAdes" = "orange3")) +
  theme(strip.text = element_text(size=16),
        axis.text.x = element_text(size = 16,face="bold"),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16,face="bold",margin = margin(r = 5)),
        panel.border = element_rect(color = "darkgray", fill = NA)) +
  ylab("Percent Reference Aligned")

sim_query_coverage_df <- bind_rows(skesa_screen_screening_df,spades_screen_screening_df) %>%
  select(Isolate_ID,Dataset,Assembler,Query_Percent_Aligned)

fig_s2_query <- sim_query_coverage_df %>% 
  ggplot(aes(x=Dataset,y=Query_Percent_Aligned,color=Assembler,fill=Assembler)) +
  geom_violin(alpha = 0.3,position = position_dodge(width = 1)) +
  geom_jitter(position = position_jitterdodge(jitter.width = 0.1,dodge.width = 1)) + 
  theme_minimal() +
  scale_color_manual(values = c("SKESA" = "blue","SPAdes" = "orange3")) +
  scale_fill_manual(values = c("SKESA" = "blue","SPAdes" = "orange3")) +
  theme(strip.text = element_text(size=16),
        axis.text.x = element_text(size = 16,face="bold"),
        axis.text.y = element_text(size = 14),
        axis.title.x = element_blank(),
        axis.title.y = element_text(size = 16,face="bold",margin = margin(r = 5)),
        panel.border = element_rect(color = "darkgray", fill = NA)) +
  ylab("Percent Query Aligned")

fig_s2_plot <- cowplot::plot_grid(fig_s2_ref,fig_s2_query,rel_widths = c(1,1.1))

fig_s3 <- ggplot(isolate_df %>% filter(!Isolate_ID %in% fail_contig_isolates),aes(x=Analysis_ID,y=Contig_Count)) +
    geom_boxplot() +
    theme_minimal() +
    theme_minimal() + 
    theme(strip.text = element_text(size=16),
          axis.text.x = element_text(size = 12,face="bold",angle=90),
          axis.text.y = element_text(size = 14),
          axis.title.x = element_blank(),
          axis.title.y = element_text(size = 16,face="bold",margin = margin(r = 5)),
          panel.border = element_rect(color = "darkgray", fill = NA)) +
    facet_wrap(~Species,ncol=1,scales="free") +
    ylab("Contig Count")

fail_qc_isolates <- csp2_ref_screening_df %>% filter(Screen_Category != "Pass") %>% pull(Query_ID) %>% unique()
fail_isolates <- c(fail_contig_isolates,fail_qc_isolates)

fig_s4_df <- comparison_querycov_df %>%
  bind_rows(comparison_refcov_df) %>%
  filter(!Seq_1 %in% fail_isolates & !Seq_2 %in% fail_isolates) %>%
  mutate(Reference_Type = ifelse(Reference_ID %in% csp2_refs,"CSP2","CSP1"))

fig_s4 <- cowplot::plot_grid(ggplot(fig_s4_df %>% filter(Variable == "Reference_Coverage"),aes(x=Analysis_ID,y=Value,color=Reference_Type,fill=Reference_Type)) +
  geom_boxplot() +
  theme_minimal() +
  facet_wrap(~Species,scale="free_x",nrow=1) +
  theme(axis.text.x = element_text(angle=90)) +
  scale_color_manual(values = c("CSP2" = "blue","CSP1" = "orange3")) +
  scale_fill_manual(values = c("CSP2" = "blue","CSP1" = "orange3")),
  ggplot(fig_s4_df %>% filter(Variable == "Query_Coverage"),aes(x=Analysis_ID,y=Value,color=Reference_Type,fill=Reference_Type)) +
    geom_boxplot() +
    theme_minimal() +
    facet_wrap(~Species,scale="free_x",nrow=1) +
    theme(axis.text.x = element_text(angle=90)) +
    scale_color_manual(values = c("CSP2" = "blue","CSP1" = "orange3")) +
    scale_fill_manual(values = c("CSP2" = "blue","CSP1" = "orange3")),nrow=2)
  
