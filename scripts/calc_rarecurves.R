library(dplyr)
library(readr)
library(purrr)
library(tidyr)
library(tibble)
library(ggplot2)
library(vegan)

# read in metadata which maps metagenome sample accession to sample type
# metadata <- read_tsv("inputs/metadata-paired-mgx-mtx.tsv")
metadata <- read_tsv(snakemake@input[['metadata']])


# create a vector of file paths
#files <- list.files(path = "outputs/sourmash_sketch_downsample_filtered_csv", pattern = ".csv$",
#                    full.names = TRUE, recursive = TRUE)
files <- snakemake@input[['csvs']]

# create a dataframe that maps metagenome sample accession to file path
files_df <- data.frame(path = files) %>%
  mutate(mgx_run_accession = gsub("_k[235]1_scaled100k.csv", "", basename(path)))

# use snakemake wild cards to subset sample type 
# this could be run as a for loop instead, but this way snakemake will parallelize the rarefaction curves which are slow to build
# loop over every sample type and create per-sample rarefaction curves and mean slope to estimate convergence of rarefaction curve.
sample_type_i <- snakemake@wildcards[['sample_type']]

#sample_type_i = "activated_sludge"
# filter metadata to metagenomes of a specific sample type
metadata_i <- metadata %>%
  filter(sample_type == sample_type_i)

# filter to files paths for metagenomes of specific sample type
files_df_i <- files_df %>%
  filter(mgx_run_accession %in% metadata_i$mgx_run_accession)

# read in sketches and convert to a wide data frame
sigs_i <- files_df_i$path %>%
  map_dfr(read_csv, col_types = "cdc") %>%
  pivot_wider(id_cols = name, values_from = "abund", names_from = "hash") %>%
  column_to_rownames("name") %>%
  replace(is.na(.), 0)

raremax_i <- min(rowSums(sigs_i))

# return as a tidy data frame with per-sample information to calculate mean slope of each line
rarecurve_i <- rarecurve(sigs_i, step =1, sample = raremax_i, col = "grey", cex = 0.4, tidy = T,
                         xlab = "Number of k-mers sampled", ylab = "Number of distinct k-mers observed")

# rename columns to something that makes sense
colnames(rarecurve_i) <- c("mgx_run_accession", "num_kmers_sampled", "num_kmers_observed")

# plot the results and save the plot as a pdf
pdf(snakemake@output[['pdf']])
#pdf(paste0("rarecurve_", sample_type_i, ".pdf"))
ggplot(rarecurve_i, aes(x = num_kmers_sampled, y = num_kmers_observed, color = mgx_run_accession)) +
  geom_point() +
  theme_minimal()
dev.off()

# output the tsv file recording the per-sample 
write_tsv(rarecurve_i, snakemake@output[['tsv']])

# for each sample, calculate the slope of the rarefaction curve
slope_df <- data.frame()
for(mgx_run_accession_i in files_df_i$mgx_run_accession){
  rarecurve_mgx_run_accession_i <- rarecurve_i %>% 
    filter(mgx_run_accession == mgx_run_accession_i)
  run_i <- diff(rarecurve_mgx_run_accession_i$num_kmers_sampled)
  rise_i <- diff(rarecurve_mgx_run_accession_i$num_kmers_observed)
  slope_i <- rise_i/run_i
  mean_slope_i <- mean(slope_i)
  min_nonzero_slope_i <- min(slope_i[slope_i != 0]) # filter out zero values and then calculate the minimum slope observed
  slope_df_i <- data.frame(mgx_run_accession = mgx_run_accession_i,
                           mean_slope = mean_slope_i,
                           min_nonzero_slope = min_nonzero_slope_i)
  slope_df <- bind_rows(slope_df, slope_df_i)
}


write_tsv(slope_df, snakemake@output[['slopes']])
