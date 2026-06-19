############################################
## CLEAN PEPTIDE PIPELINE (CAN vs NONCAN)
############################################

library(tidyverse)

ROOT <- "C:/Project_R4"

NONCAN_FILE <- file.path(ROOT, "output_non_canonical/combined_peptide.tsv")
CAN_FILE    <- file.path(ROOT, "output_data_4/combined_peptide.tsv")

OUT_FIG <- file.path(ROOT, "results/figures")
OUT_TAB <- file.path(ROOT, "results/tables")

dir.create(OUT_FIG, recursive = TRUE, showWarnings = FALSE)
dir.create(OUT_TAB, recursive = TRUE, showWarnings = FALSE)
noncan <- read_tsv(NONCAN_FILE, show_col_types = FALSE)
canon  <- read_tsv(CAN_FILE, show_col_types = FALSE)

noncan$type <- "noncanonical"
canon$type  <- "canonical"

df <- bind_rows(noncan, canon)

colnames(df) <- trimws(colnames(df))
df <- df %>%
  mutate(
    length = nchar(as.character(Sequence)),
    prob = as.numeric(Probability)
  )
p1 <- ggplot(df, aes(x = length, fill = type)) +
  geom_histogram(bins = 40, alpha = 0.5, position = "identity") +
  theme_minimal() +
  labs(title = "Peptide length distribution")

ggsave(file.path(OUT_FIG, "01_length.png"), p1, width = 7, height = 5)
p2 <- ggplot(df, aes(x = prob, fill = type)) +
  geom_histogram(bins = 40, alpha = 0.5, position = "identity") +
  theme_minimal() +
  labs(title = "Identification probability distribution")

ggsave(file.path(OUT_FIG, "02_probability.png"), p2, width = 7, height = 5)
p3 <- ggplot(df, aes(x = `Charge States`, fill = type)) +
  geom_histogram(position = "dodge", bins = 20) +
  theme_minimal() +
  labs(title = "Charge state distribution")

ggsave(file.path(OUT_FIG, "03_charge.png"), p3, width = 7, height = 5)
p4 <- ggplot(df, aes(x = log10(`GRLCL_1D_ETHCD Spectral Count` + 1), fill = type)) +
  geom_histogram(bins = 40, alpha = 0.5) +
  theme_minimal() +
  labs(title = "Spectral evidence (log scale)")

ggsave(file.path(OUT_FIG, "04_spectral.png"), p4, width = 7, height = 5)
protein_usage <- df %>%
  group_by(Protein, type) %>%
  summarise(n = n(), .groups = "drop")

p5 <- ggplot(protein_usage, aes(x = n, fill = type)) +
  geom_histogram(bins = 40, alpha = 0.5) +
  theme_minimal() +
  labs(title = "Protein-level peptide reuse")

ggsave(file.path(OUT_FIG, "05_protein_reuse.png"), p5, width = 7, height = 5)
summary_df <- df %>%
  group_by(type) %>%
  summarise(mean_prob = mean(prob, na.rm = TRUE))

p6 <- ggplot(summary_df, aes(x = type, y = mean_prob, fill = type)) +
  geom_col() +
  theme_minimal() +
  labs(title = "Mean identification confidence")

ggsave(file.path(OUT_FIG, "06_summary.png"), p6, width = 6, height = 5)
print("DONE: all figures saved")
list.files(OUT_FIG)



