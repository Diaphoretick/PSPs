
R version 4.6.0 (2026-04-24 ucrt) -- "Because it was There"
Copyright (C) 2026 The R Foundation for Statistical Computing
Platform: x86_64-w64-mingw32/x64

R is free software and comes with ABSOLUTELY NO WARRANTY.
You are welcome to redistribute it under certain conditions.
Type 'license()' or 'licence()' for distribution details.

R is a collaborative project with many contributors.
Type 'contributors()' for more information and
'citation()' on how to cite R or R packages in publications.

Type 'demo()' for some demos, 'help()' for on-line help, or
'help.start()' for an HTML browser interface to help.
Type 'q()' to quit R.

[Workspace loaded from ~/.RData]

> load("~/.RData")
> ############################################
> ## COMPARATIVE PROTEOMICS ANALYSIS
  > ############################################
> 
  > library(tidyverse)
── Attaching core tidyverse packages ──────────────────────────── tidyverse 2.0.0 ──
✔ dplyr     1.2.1     ✔ readr     2.2.0
✔ forcats   1.0.1     ✔ stringr   1.6.0
✔ ggplot2   4.0.3     ✔ tibble    3.3.1
✔ lubridate 1.9.5     ✔ tidyr     1.3.2
✔ purrr     1.2.2     
── Conflicts ────────────────────────────────────────────── tidyverse_conflicts() ──
✖ dplyr::filter() masks stats::filter()
✖ dplyr::lag()    masks stats::lag()
ℹ Use the conflicted package to force all conflicts to become errors
> 
  > ############################################
> ## INPUT FILES
  > ############################################
> 
  > NONCANON_FILE <- "C:/Project_R4/output_non_canonical/combined_peptide.tsv"
> CANON_FILE    <- "C:/Project_R4/output_data_4/combined_peptide.tsv"
> 
  > OUTDIR <- "C:/Project_R4/comparative_analysis"
> 
  > dir.create(OUTDIR,
               +            showWarnings = FALSE,
               +            recursive = TRUE)
> 
  > ############################################
> ## LOAD DATA
  > ############################################
> 
  > noncanon <- read_tsv(
    +     NONCANON_FILE,
    +     show_col_types = FALSE
    + )
New names:                                                                          
  • `` -> `...15`
> 
  > canon <- read_tsv(
    +     CANON_FILE,
    +     show_col_types = FALSE
    + )
New names:                                                                          
  • `` -> `...15`
> 
  > noncanon$type <- "Non-canonical"
> canon$type <- "Canonical"
> 
  > df <- bind_rows(noncanon, canon)
> 
  > colnames(df) <- tolower(colnames(df))
> colnames(df) <- gsub(" ", "_", colnames(df))
> 
  > ############################################
> ## FEATURES
  > ############################################
> 
  > df <- df %>%
  +     mutate(
    +         length = nchar(sequence),
    +         probability = as.numeric(probability)
    +     )
> 
  > ############################################
> ## FIGURE 1
  > ## LENGTH DISTRIBUTION
  > ############################################
> 
  > p1 <- ggplot(
    +     df,
    +     aes(length, fill = type)
    + ) +
  +     geom_density(alpha = 0.4) +
  +     theme_minimal() +
  +     labs(
    +         title = "Peptide Length Distribution"
    +     )
> 
  > ggsave(
    +     file.path(OUTDIR,"F1_length_distribution.png"),
    +     p1,
    +     width = 7,
    +     height = 5
    + )
> 
  > ############################################
> ## FIGURE 2
  > ## IDENTIFICATION PROBABILITY
  > ############################################
> 
  > p2 <- ggplot(
    +     df,
    +     aes(probability,
              +         color = type)
    + ) +
  +     geom_density() +
  +     theme_minimal() +
  +     labs(
    +         title = "Identification Probability Distribution"
    +     )
> 
  > ggsave(
    +     file.path(OUTDIR,"F2_probability_distribution.png"),
    +     p2,
    +     width = 7,
    +     height = 5
    + )
> 
  > ############################################
> ## FIGURE 3
  > ## CHARGE STATES
  > ############################################
> 
  > p3 <- ggplot(
    +     df,
    +     aes(
      +         factor(charge_states),
      +         fill = type
      +     )
    + ) +
  +     geom_bar(position="dodge") +
  +     theme_minimal() +
  +     labs(
    +         title = "Charge State Distribution",
    +         x = "Charge State"
    +     )
> 
  > ggsave(
    +     file.path(OUTDIR,"F3_charge_states.png"),
    +     p3,
    +     width = 7,
    +     height = 5
    + )
> 
  > ############################################
> ## FIGURE 4
  > ## PROTEIN REPRESENTATION
  > ############################################
> 
  > protein_rep <- df %>%
  +     group_by(type, protein) %>%
  +     summarise(
    +         peptide_count = n(),
    +         .groups = "drop"
    +     )
> 
  > p4 <- ggplot(
    +     protein_rep,
    +     aes(
      +         peptide_count,
      +         fill = type
      +     )
    + ) +
  +     geom_histogram(
    +         bins = 40,
    +         alpha = 0.6,
    +         position = "identity"
    +     ) +
  +     theme_minimal() +
  +     labs(
    +         title = "Protein Representation"
    +     )
> 
  > ggsave(
    +     file.path(OUTDIR,"F4_protein_representation.png"),
    +     p4,
    +     width = 7,
    +     height = 5
    + )
> 
  > ############################################
> ## FIGURE 5
  > ## RANK-ABUNDANCE CURVE
  > ############################################
> 
  > rank_df <- protein_rep %>%
  +     arrange(type, desc(peptide_count)) %>%
  +     group_by(type) %>%
  +     mutate(rank = row_number())
> 
  > p5 <- ggplot(
    +     rank_df,
    +     aes(rank,
              +         peptide_count,
              +         color = type)
    + ) +
  +     geom_line() +
  +     scale_y_log10() +
  +     theme_minimal() +
  +     labs(
    +         title = "Protein Rank-Abundance Curve"
    +     )
> 
  > ggsave(
    +     file.path(OUTDIR,"F5_rank_abundance.png"),
    +     p5,
    +     width = 7,
    +     height = 5
    + )
> 
  > ############################################
> ## FIGURE 6
  > ## SEQUENCE COMPLEXITY
  > ############################################
> 
  > complexity_fun <- function(x){
    +     
      +     length(unique(
        +         unlist(
          +             strsplit(x,"")
          +         )
        +     ))
    + }
> 
  > df$complexity <- sapply(
    +     df$sequence,
    +     complexity_fun
    + )
> 
  > p6 <- ggplot(
    +     df,
    +     aes(complexity,
              +         fill = type)
    + ) +
  +     geom_density(alpha = 0.4) +
  +     theme_minimal() +
  +     labs(
    +         title = "Sequence Complexity"
    +     )
> 
  > ggsave(
    +     file.path(OUTDIR,"F6_complexity.png"),
    +     p6,
    +     width = 7,
    +     height = 5
    + )
> 
  > ############################################
> ## FIGURE 7
  > ## AMINO ACID COMPOSITION
  > ############################################
> 
  > AA <- c(
    +     "A","C","D","E","F",
    +     "G","H","I","K","L",
    +     "M","N","P","Q","R",
    +     "S","T","V","W","Y"
    + )
> 
  > aa_summary <- map_dfr(
    +     AA,
    +     function(aa){
      +         
        +         df %>%
        +             mutate(
          +                 freq =
            +                     stringr::str_count(
              +                         sequence,
              +                         aa
              +                     ) / length
          +             ) %>%
        +             group_by(type) %>%
        +             summarise(
          +                 mean_freq = mean(freq),
          +                 .groups="drop"
          +             ) %>%
        +             mutate(amino_acid = aa)
      +         
        +     }
    + )
> 
  > p7 <- ggplot(
    +     aa_summary,
    +     aes(
      +         amino_acid,
      +         mean_freq,
      +         fill = type
      +     )
    + ) +
  +     geom_col(position="dodge") +
  +     theme_minimal() +
  +     labs(
    +         title = "Amino Acid Composition"
    +     )
> 
  > ggsave(
    +     file.path(OUTDIR,"F7_amino_acid_composition.png"),
    +     p7,
    +     width = 9,
    +     height = 5
    + )
> 
  > ############################################
> ## FIGURE 8
  > ## PROTEIN DIVERSITY
  > ############################################
> 
  > diversity_df <- protein_rep %>%
  +     group_by(type) %>%
  +     summarise(
    +         proteins = n_distinct(protein),
    +         mean_peptides = mean(peptide_count)
    +     )
> 
  > write_csv(
    +     diversity_df,
    +     file.path(
      +         OUTDIR,
      +         "protein_diversity_summary.csv"
      +     )
    + )
>                                                                                 
  > ############################################
> ## FIGURE 9
  > ## TOP PROTEINS
  > ############################################
> 
  > top_proteins <- protein_rep %>%
  +     group_by(type) %>%
  +     slice_max(
    +         peptide_count,
    +         n = 15
    +     )
> 
  > p9 <- ggplot(
    +     top_proteins,
    +     aes(
      +         reorder(protein,
                        +                 peptide_count),
      +         peptide_count,
      +         fill = type
      +     )
    + ) +
  +     geom_col() +
  +     coord_flip() +
  +     theme_minimal() +
  +     labs(
    +         title = "Top Protein Sources"
    +     )
> 
  > ggsave(
    +     file.path(OUTDIR,"F9_top_proteins.png"),
    +     p9,
    +     width = 10,
    +     height = 6
    + )
> 
  > ############################################
> ## FIGURE 10
  > ## CUMULATIVE COVERAGE
  > ############################################
> 
  > coverage <- protein_rep %>%
  +     arrange(type,
                +             desc(peptide_count)) %>%
  +     group_by(type) %>%
  +     mutate(
    +         cumulative =
      +             cumsum(peptide_count)/
      +             sum(peptide_count),
    +         rank = row_number()
    +     )
> 
  > p10 <- ggplot(
    +     coverage,
    +     aes(
      +         rank,
      +         cumulative,
      +         color = type
      +     )
    + ) +
  +     geom_line(linewidth=1) +
  +     theme_minimal() +
  +     labs(
    +         title = "Cumulative Protein Contribution"
    +     )
> 
  > ggsave(
    +     file.path(OUTDIR,"F10_cumulative_coverage.png"),
    +     p10,
    +     width = 7,
    +     height = 5
    + )
> 
  > ############################################
> ## SUMMARY TABLE
  > ############################################
> 
  > summary_table <- df %>%
  +     group_by(type) %>%
  +     summarise(
    +         peptides = n(),
    +         proteins = n_distinct(protein),
    +         mean_length = mean(length),
    +         median_length = median(length),
    +         mean_probability = mean(probability,
                                      +                                 na.rm = TRUE)
    +     )
> 
  > write_csv(
    +     summary_table,
    +     file.path(
      +         OUTDIR,
      +         "summary_statistics.csv"
      +     )
    + )
>                                                                                 
  > cat("Analysis complete\n")
Analysis complete
> ############################################
> ## 0. PACKAGES (AUTO INSTALL SAFE)
  > ############################################
> 
  > pkgs_cran <- c(
    +     "tidyverse",
    +     "data.table",
    +     "stringr",
    +     "ggplot2",
    +     "scales"
    + )
> 
  > pkgs_bioc <- c(
    +     "limma",
    +     "BiocManager"
    + )
> 
  > install_if_missing <- function(pkgs){
    +     for(p in pkgs){
      +         if(!requireNamespace(p, quietly = TRUE)){
        +             install.packages(p, dependencies = TRUE)
        +         }
      +         library(p, character.only = TRUE)
      +     }
    + }
> 
  > install_if_missing(pkgs_cran)
data.table 1.18.4 using 2 threads (see ?getDTthreads).  Latest news: r-datatable.com

Attaching package: ‘data.table’

The following objects are masked from ‘package:lubridate’:
  
  hour, isoweek, isoyear, mday, minute, month, quarter, second, wday,
week, yday, year

The following objects are masked from ‘package:dplyr’:
  
  between, first, last

The following object is masked from ‘package:purrr’:
  
  transpose

The following object is masked from ‘package:base’:
  
  %notin%
  
  
  Attaching package: ‘scales’

The following object is masked from ‘package:purrr’:
  
  discard

The following object is masked from ‘package:readr’:
  
  col_factor
> 
  > if(!requireNamespace("BiocManager", quietly = TRUE)){
    +     install.packages("BiocManager")
    + }
> 
  > if(!requireNamespace("limma", quietly = TRUE)){
    +     BiocManager::install("limma", ask = FALSE, update = FALSE)
    + }
> 
  > library(limma)
> 
  > ############################################
> ## 1. FILE PATHS
  > ############################################
> 
  > NONCANON_FILE <- "C:/Project_R4/output_non_canonical/combined_peptide.tsv"
> CANON_FILE    <- "C:/Project_R4/output_data_4/combined_peptide.tsv"
> 
  > OUTDIR <- "C:/Project_R4/defense_figures"
> dir.create(OUTDIR, showWarnings = FALSE, recursive = TRUE)
> 
  > ############################################
> ## 2. LOAD DATA
  > ############################################
> 
  > noncanon <- read_tsv(NONCANON_FILE, show_col_types = FALSE)
New names:                                                                          
  • `` -> `...15`
> canon    <- read_tsv(CANON_FILE, show_col_types = FALSE)
New names:                                                                          
  • `` -> `...15`
> 
  > noncanon$type <- "Non-canonical"
> canon$type    <- "Canonical"
> 
  > df <- bind_rows(noncanon, canon)
> 
  > colnames(df) <- tolower(colnames(df))
> colnames(df) <- gsub(" ", "_", colnames(df))
> 
  > ############################################
> ## 3. FEATURE ENGINEERING
  > ############################################
> 
  > df <- df %>%
  +     mutate(
    +         sequence = as.character(sequence),
    +         length = nchar(sequence),
    +         probability = as.numeric(probability),
    +         charge_states = as.factor(charge_states)
    +     )
> 
  > ############################################
> ## 4. FILTER (LIGHT QC ONLY)
  > ############################################
> 
  > df <- df %>%
  +     filter(length >= 6, length <= 40)
> 
  > ############################################
> ## 5. FIGURE 1 — LENGTH DISTRIBUTION
  > ############################################
> 
  > p1 <- ggplot(df, aes(length, fill = type)) +
  +     geom_density(alpha = 0.4) +
  +     theme_minimal() +
  +     labs(title = "Peptide Length Distribution",
             +          x = "Length", y = "Density")
> 
  > ggsave(file.path(OUTDIR,"F1_length.png"), p1, width = 7, height = 5)
> 
  > ############################################
> ## 6. FIGURE 2 — IDENTIFICATION CONFIDENCE
  > ############################################
> 
  > p2 <- ggplot(df, aes(probability, fill = type)) +
  +     geom_density(alpha = 0.4) +
  +     theme_minimal() +
  +     labs(title = "Identification Confidence (Probability)")
> 
  > ggsave(file.path(OUTDIR,"F2_probability.png"), p2, width = 7, height = 5)
> 
  > ############################################
> ## 7. FIGURE 3 — CHARGE STATES
  > ############################################
> 
  > p3 <- ggplot(df, aes(charge_states, fill = type)) +
  +     geom_bar(position = "dodge") +
  +     theme_minimal() +
  +     labs(title = "Charge State Distribution")
> 
  > ggsave(file.path(OUTDIR,"F3_charge.png"), p3, width = 7, height = 5)
> 
  > ############################################
> ## 8. FIGURE 4 — PROTEIN COVERAGE
  > ############################################
> 
  > protein_map <- df %>%
  +     group_by(type, protein) %>%
  +     summarise(n_pep = n(), .groups = "drop")
> 
  > p4 <- ggplot(protein_map, aes(n_pep, fill = type)) +
  +     geom_histogram(bins = 40, alpha = 0.6) +
  +     theme_minimal() +
  +     labs(title = "Protein Representation (Peptides per Protein)")
> 
  > ggsave(file.path(OUTDIR,"F4_protein_coverage.png"), p4, width = 7, height = 5)
> 
  > ############################################
> ## 9. FIGURE 5 — RANK ABUNDANCE (KEY SLIDE)
  > ############################################
> 
  > rank_df <- protein_map %>%
  +     group_by(type) %>%
  +     arrange(desc(n_pep)) %>%
  +     mutate(rank = row_number())
> 
  > p5 <- ggplot(rank_df, aes(rank, n_pep, color = type)) +
  +     geom_line(linewidth = 1) +
  +     scale_y_log10() +
  +     theme_minimal() +
  +     labs(title = "Protein Rank-Abundance (Core vs Rare Proteins)")
> 
  > ggsave(file.path(OUTDIR,"F5_rank_abundance.png"), p5, width = 7, height = 5)
> 
  > ############################################
> ## 10. FIGURE 6 — AMINO ACID SIGNATURE
  > ############################################
> 
  > aa <- c("A","C","D","E","F","G","H","I","K","L",
            +         "M","N","P","Q","R","S","T","V","W","Y")
> 
  > aa_df <- map_dfr(aa, function(a){
    +     df %>%
      +         mutate(freq = str_count(sequence, a) / length) %>%
      +         group_by(type) %>%
      +         summarise(mean_freq = mean(freq, na.rm = TRUE), .groups = "drop") %>%
      +         mutate(aa = a)
    + })
> 
  > p6 <- ggplot(aa_df, aes(aa, mean_freq, fill = type)) +
  +     geom_col(position = "dodge") +
  +     theme_minimal() +
  +     labs(title = "Amino Acid Composition Bias")
> 
  > ggsave(file.path(OUTDIR,"F6_aa_composition.png"), p6, width = 9, height = 5)
> 
  > ############################################
> ## 11. FIGURE 7 — SEQUENCE COMPLEXITY
  > ############################################
> 
  > df$complexity <- sapply(df$sequence, function(x){
    +     length(unique(strsplit(x,"")[[1]]))
    + })
> 
  > p7 <- ggplot(df, aes(complexity, fill = type)) +
  +     geom_density(alpha = 0.4) +
  +     theme_minimal() +
  +     labs(title = "Sequence Complexity Distribution")
> 
  > ggsave(file.path(OUTDIR,"F7_complexity.png"), p7, width = 7, height = 5)
> 
  > ############################################
> ## 12. FIGURE 8 — TOP PROTEINS (BIOLOGICAL DRIVERS)
  > ############################################
> 
  > top_prot <- protein_map %>%
  +     group_by(type) %>%
  +     slice_max(n_pep, n = 12)
> 
  > p8 <- ggplot(top_prot,
                 +              aes(reorder(protein, n_pep), n_pep, fill = type)) +
  +     geom_col() +
  +     coord_flip() +
  +     theme_minimal() +
  +     labs(title = "Top Protein Sources")
> 
  > ggsave(file.path(OUTDIR,"F8_top_proteins.png"), p8, width = 10, height = 6)
> 
  > ############################################
> ## 13. FIGURE 9 — CUMULATIVE CONTRIBUTION
  > ############################################
> 
  > cum <- rank_df %>%
  +     group_by(type) %>%
  +     mutate(cum = cumsum(n_pep)/sum(n_pep),
               +            rank = row_number())
> 
  > p9 <- ggplot(cum, aes(rank, cum, color = type)) +
  +     geom_line(linewidth = 1) +
  +     theme_minimal() +
  +     labs(title = "Cumulative Protein Contribution")
> 
  > ggsave(file.path(OUTDIR,"F9_cumulative.png"), p9, width = 7, height = 5)
> 
  > ############################################
> ## 14. SUMMARY TABLE
  > ############################################
> 
  > summary <- df %>%
  +     group_by(type) %>%
  +     summarise(
    +         peptides = n(),
    +         proteins = n_distinct(protein),
    +         mean_length = mean(length),
    +         mean_prob = mean(probability, na.rm = TRUE)
    +     )
> 
  > write_csv(summary, file.path(OUTDIR,"summary.csv"))
>                                                                                 
  > cat("DONE — DEFENSE FIGURES GENERATED\n")
DONE — DEFENSE FIGURES GENERATED
> pkgs_extra <- c("FactoMineR", "factoextra", "umap", "reshape2")
> 
  > install_if_missing(pkgs_extra)
WARNING: Rtools is required to build R packages but is not currently installed. Please download and install the appropriate version of Rtools before proceeding:
  
  https://cran.rstudio.com/bin/windows/Rtools/
  Installing package into ‘C:/Users/ZHUK/AppData/Local/R/win-library/4.6’
(as ‘lib’ is unspecified)
also installing the dependencies ‘ucminf’, ‘colorspace’, ‘fracdiff’, ‘lmtest’, ‘timeDate’, ‘urca’, ‘zoo’, ‘RcppArmadillo’, ‘ordinal’, ‘cowplot’, ‘Deriv’, ‘forecast’, ‘microbenchmark’, ‘rbibutils’, ‘shape’, ‘pan’, ‘jomo’, ‘doBy’, ‘SparseM’, ‘MatrixModels’, ‘Rdpack’, ‘minqa’, ‘nloptr’, ‘reformulas’, ‘RcppEigen’, ‘lazyeval’, ‘later’, ‘glmnet’, ‘mitml’, ‘iterators’, ‘commonmark’, ‘httpuv’, ‘sourcetools’, ‘xtable’, ‘miniUI’, ‘shinyjs’, ‘carData’, ‘abind’, ‘Formula’, ‘pbkrtest’, ‘quantreg’, ‘lme4’, ‘crosstalk’, ‘promises’, ‘estimability’, ‘mvtnorm’, ‘numDeriv’, ‘mice’, ‘doParallel’, ‘foreach’, ‘shiny’, ‘FactoInvestigate’, ‘colourpicker’, ‘shinydashboard’, ‘shinyjqui’, ‘litedown’, ‘car’, ‘DT’, ‘ellipse’, ‘emmeans’, ‘flashClust’, ‘leaps’, ‘multcompView’, ‘scatterplot3d’, ‘ggrepel’, ‘irlba’, ‘missMDA’, ‘Factoshiny’, ‘markdown’
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/ucminf_1.2.3.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/colorspace_2.1-2.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/fracdiff_1.5-4.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/lmtest_0.9-40.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/timeDate_4052.112.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/urca_1.3-4.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/zoo_1.8-15.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/RcppArmadillo_15.2.7-1.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/ordinal_2025.12-29.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/cowplot_1.2.0.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/Deriv_4.2.0.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/forecast_9.0.2.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/microbenchmark_1.5.0.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/rbibutils_2.4.1.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/shape_1.4.6.1.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/pan_1.9.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/jomo_2.7-6.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/doBy_4.7.1.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/SparseM_1.84-2.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/MatrixModels_0.5-4.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/Rdpack_2.6.6.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/minqa_1.2.8.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/nloptr_2.2.1.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/reformulas_0.4.4.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/RcppEigen_0.3.4.0.2.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/lazyeval_0.2.3.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/later_1.4.8.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/glmnet_5.0.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/mitml_0.4-5.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/iterators_1.0.14.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/commonmark_2.0.0.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/httpuv_1.6.17.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/sourcetools_0.1.7-2.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/xtable_1.8-8.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/miniUI_0.1.2.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/shinyjs_2.1.1.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/carData_3.0-6.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/abind_1.4-8.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/Formula_1.2-5.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/pbkrtest_0.5.5.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/quantreg_6.1.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/lme4_2.0-1.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/crosstalk_1.2.2.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/promises_1.5.0.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/estimability_1.5.1.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/mvtnorm_1.4-1.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/numDeriv_2016.8-1.1.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/mice_3.19.0.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/doParallel_1.0.17.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/foreach_1.5.2.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/shiny_1.13.0.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/FactoInvestigate_1.9.1.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/colourpicker_1.3.0.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/shinydashboard_0.7.3.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/shinyjqui_0.4.1.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/litedown_0.9.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/car_3.1-5.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/DT_0.34.0.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/ellipse_0.5.0.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/emmeans_2.0.3.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/flashClust_1.1-4.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/leaps_3.2.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/multcompView_0.1-11.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/scatterplot3d_0.3-45.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/ggrepel_0.9.8.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/irlba_2.3.7.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/missMDA_1.21.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/Factoshiny_2.7.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/markdown_2.0.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/FactoMineR_2.15.zip'
package ‘ucminf’ successfully unpacked and MD5 sums checked
package ‘colorspace’ successfully unpacked and MD5 sums checked
package ‘fracdiff’ successfully unpacked and MD5 sums checked
package ‘lmtest’ successfully unpacked and MD5 sums checked
package ‘timeDate’ successfully unpacked and MD5 sums checked
package ‘urca’ successfully unpacked and MD5 sums checked
package ‘zoo’ successfully unpacked and MD5 sums checked
package ‘RcppArmadillo’ successfully unpacked and MD5 sums checked
package ‘ordinal’ successfully unpacked and MD5 sums checked
package ‘cowplot’ successfully unpacked and MD5 sums checked
package ‘Deriv’ successfully unpacked and MD5 sums checked
package ‘forecast’ successfully unpacked and MD5 sums checked
package ‘microbenchmark’ successfully unpacked and MD5 sums checked
package ‘rbibutils’ successfully unpacked and MD5 sums checked
package ‘shape’ successfully unpacked and MD5 sums checked
package ‘pan’ successfully unpacked and MD5 sums checked
package ‘jomo’ successfully unpacked and MD5 sums checked
package ‘doBy’ successfully unpacked and MD5 sums checked
package ‘SparseM’ successfully unpacked and MD5 sums checked
package ‘MatrixModels’ successfully unpacked and MD5 sums checked
package ‘Rdpack’ successfully unpacked and MD5 sums checked
package ‘minqa’ successfully unpacked and MD5 sums checked
package ‘nloptr’ successfully unpacked and MD5 sums checked
package ‘reformulas’ successfully unpacked and MD5 sums checked
package ‘RcppEigen’ successfully unpacked and MD5 sums checked
package ‘lazyeval’ successfully unpacked and MD5 sums checked
package ‘later’ successfully unpacked and MD5 sums checked
package ‘glmnet’ successfully unpacked and MD5 sums checked
package ‘mitml’ successfully unpacked and MD5 sums checked
package ‘iterators’ successfully unpacked and MD5 sums checked
package ‘commonmark’ successfully unpacked and MD5 sums checked
package ‘httpuv’ successfully unpacked and MD5 sums checked
package ‘sourcetools’ successfully unpacked and MD5 sums checked
package ‘xtable’ successfully unpacked and MD5 sums checked
package ‘miniUI’ successfully unpacked and MD5 sums checked
package ‘shinyjs’ successfully unpacked and MD5 sums checked
package ‘carData’ successfully unpacked and MD5 sums checked
package ‘abind’ successfully unpacked and MD5 sums checked
package ‘Formula’ successfully unpacked and MD5 sums checked
package ‘pbkrtest’ successfully unpacked and MD5 sums checked
package ‘quantreg’ successfully unpacked and MD5 sums checked
package ‘lme4’ successfully unpacked and MD5 sums checked
package ‘crosstalk’ successfully unpacked and MD5 sums checked
package ‘promises’ successfully unpacked and MD5 sums checked
package ‘estimability’ successfully unpacked and MD5 sums checked
package ‘mvtnorm’ successfully unpacked and MD5 sums checked
package ‘numDeriv’ successfully unpacked and MD5 sums checked
package ‘mice’ successfully unpacked and MD5 sums checked
package ‘doParallel’ successfully unpacked and MD5 sums checked
package ‘foreach’ successfully unpacked and MD5 sums checked
package ‘shiny’ successfully unpacked and MD5 sums checked
package ‘FactoInvestigate’ successfully unpacked and MD5 sums checked
package ‘colourpicker’ successfully unpacked and MD5 sums checked
package ‘shinydashboard’ successfully unpacked and MD5 sums checked
package ‘shinyjqui’ successfully unpacked and MD5 sums checked
package ‘litedown’ successfully unpacked and MD5 sums checked
package ‘car’ successfully unpacked and MD5 sums checked
package ‘DT’ successfully unpacked and MD5 sums checked
package ‘ellipse’ successfully unpacked and MD5 sums checked
package ‘emmeans’ successfully unpacked and MD5 sums checked
package ‘flashClust’ successfully unpacked and MD5 sums checked
package ‘leaps’ successfully unpacked and MD5 sums checked
package ‘multcompView’ successfully unpacked and MD5 sums checked
package ‘scatterplot3d’ successfully unpacked and MD5 sums checked
package ‘ggrepel’ successfully unpacked and MD5 sums checked
package ‘irlba’ successfully unpacked and MD5 sums checked
package ‘missMDA’ successfully unpacked and MD5 sums checked
package ‘Factoshiny’ successfully unpacked and MD5 sums checked
package ‘markdown’ successfully unpacked and MD5 sums checked
package ‘FactoMineR’ successfully unpacked and MD5 sums checked

The downloaded binary packages are in
C:\Users\ZHUK\AppData\Local\Temp\Rtmp4EPb2o\downloaded_packages
WARNING: Rtools is required to build R packages but is not currently installed. Please download and install the appropriate version of Rtools before proceeding:
  
  https://cran.rstudio.com/bin/windows/Rtools/
  Installing package into ‘C:/Users/ZHUK/AppData/Local/R/win-library/4.6’
(as ‘lib’ is unspecified)
also installing the dependencies ‘corrplot’, ‘viridis’, ‘ggsci’, ‘ggsignif’, ‘gridExtra’, ‘polynom’, ‘rstatix’, ‘pixmap’, ‘sp’, ‘dendextend’, ‘ggpubr’, ‘ade4’, ‘ca’, ‘mclust’
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/corrplot_0.95.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/viridis_0.6.5.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/ggsci_5.0.0.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/ggsignif_0.6.4.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/gridExtra_2.3.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/polynom_1.4-1.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/rstatix_0.7.3.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/pixmap_0.4-14.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/sp_2.2-1.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/dendextend_1.19.1.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/ggpubr_0.6.3.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/ade4_1.7-24.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/ca_0.71.1.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/mclust_6.1.2.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/factoextra_2.0.0.zip'
package ‘corrplot’ successfully unpacked and MD5 sums checked
package ‘viridis’ successfully unpacked and MD5 sums checked
package ‘ggsci’ successfully unpacked and MD5 sums checked
package ‘ggsignif’ successfully unpacked and MD5 sums checked
package ‘gridExtra’ successfully unpacked and MD5 sums checked
package ‘polynom’ successfully unpacked and MD5 sums checked
package ‘rstatix’ successfully unpacked and MD5 sums checked
package ‘pixmap’ successfully unpacked and MD5 sums checked
package ‘sp’ successfully unpacked and MD5 sums checked
package ‘dendextend’ successfully unpacked and MD5 sums checked
package ‘ggpubr’ successfully unpacked and MD5 sums checked
package ‘ade4’ successfully unpacked and MD5 sums checked
package ‘ca’ successfully unpacked and MD5 sums checked
package ‘mclust’ successfully unpacked and MD5 sums checked
package ‘factoextra’ successfully unpacked and MD5 sums checked

The downloaded binary packages are in
C:\Users\ZHUK\AppData\Local\Temp\Rtmp4EPb2o\downloaded_packages
Welcome to factoextra!
  Want to learn more? See two factoextra-related books at https://www.datanovia.com/en/product/practical-guide-to-principal-component-methods-in-r/
  WARNING: Rtools is required to build R packages but is not currently installed. Please download and install the appropriate version of Rtools before proceeding:
  
  https://cran.rstudio.com/bin/windows/Rtools/
  Installing package into ‘C:/Users/ZHUK/AppData/Local/R/win-library/4.6’
(as ‘lib’ is unspecified)
also installing the dependencies ‘RcppTOML’, ‘here’, ‘png’, ‘reticulate’, ‘RSpectra’
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/RcppTOML_0.2.3.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/here_1.0.2.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/png_0.1-9.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/reticulate_1.46.0.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/RSpectra_0.16-2.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/umap_0.2.10.0.zip'
package ‘RcppTOML’ successfully unpacked and MD5 sums checked
package ‘here’ successfully unpacked and MD5 sums checked
package ‘png’ successfully unpacked and MD5 sums checked
package ‘reticulate’ successfully unpacked and MD5 sums checked
package ‘RSpectra’ successfully unpacked and MD5 sums checked
package ‘umap’ successfully unpacked and MD5 sums checked

The downloaded binary packages are in
C:\Users\ZHUK\AppData\Local\Temp\Rtmp4EPb2o\downloaded_packages
WARNING: Rtools is required to build R packages but is not currently installed. Please download and install the appropriate version of Rtools before proceeding:
  
  https://cran.rstudio.com/bin/windows/Rtools/
  Installing package into ‘C:/Users/ZHUK/AppData/Local/R/win-library/4.6’
(as ‘lib’ is unspecified)
also installing the dependencies ‘rex’, ‘plyr’, ‘covr’
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/rex_1.2.2.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/plyr_1.8.9.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/covr_3.6.5.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/reshape2_1.4.5.zip'
package ‘rex’ successfully unpacked and MD5 sums checked
package ‘plyr’ successfully unpacked and MD5 sums checked
package ‘covr’ successfully unpacked and MD5 sums checked
package ‘reshape2’ successfully unpacked and MD5 sums checked

The downloaded binary packages are in
C:\Users\ZHUK\AppData\Local\Temp\Rtmp4EPb2o\downloaded_packages

Attaching package: ‘reshape2’

The following objects are masked from ‘package:data.table’:
  
  dcast, melt

The following object is masked from ‘package:tidyr’:
  
  smiths
> library(FactoMineR)
> library(factoextra)
> library(umap)
> library(reshape2)
> df$charge_num <- as.numeric(as.character(df$charge_states))
> 
  > feature_df <- df %>%
  +     select(type, length, probability, charge_num) %>%
  +     na.omit()
> 
  > feature_scaled <- scale(feature_df[, -1])
> pca <- prcomp(feature_scaled)
> 
  > pca_df <- data.frame(
    +     PC1 = pca$x[,1],
    +     PC2 = pca$x[,2],
    +     type = feature_df$type
    + )
> 
  > p_pca <- ggplot(pca_df, aes(PC1, PC2, color = type)) +
  +     geom_point(alpha = 0.4) +
  +     theme_minimal() +
  +     labs(title = "PCA: Global Separation of Peptide Space")
> 
  > ggsave(file.path(OUTDIR,"PCA_peptides.png"), p_pca, width = 7, height = 5)
> set.seed(123)
> 
  > um <- umap(feature_scaled)
Warning message:
  failed creating initial embedding; using random embedding instead 
> 
  > um_df <- data.frame(
    +     UMAP1 = um$layout[,1],
    +     UMAP2 = um$layout[,2],
    +     type = feature_df$type
    + )
> 
  > p_umap <- ggplot(um_df, aes(UMAP1, UMAP2, color = type)) +
  +     geom_point(alpha = 0.5) +
  +     theme_minimal() +
  +     labs(title = "UMAP: Non-linear Structure of Peptide Space")
> 
  > ggsave(file.path(OUTDIR,"UMAP_peptides.png"), p_umap, width = 7, height = 5)
> feat_summary <- feature_df %>%
  +     group_by(type) %>%
  +     summarise(
    +         length = mean(length),
    +         probability = mean(probability),
    +         charge = mean(charge_num)
    +     )
> 
  > heat <- melt(as.matrix(feat_summary[,-1]))
> heat$type <- rep(feat_summary$type, each = 3)
> 
  > p_heat <- ggplot(heat, aes(Var2, type, fill = value)) +
  +     geom_tile() +
  +     theme_minimal() +
  +     labs(title = "Feature Space Shift (Canonical vs Non-canonical)")
> 
  > ggsave(file.path(OUTDIR,"feature_heatmap.png"), p_heat, width = 6, height = 4)
> protein_signal <- df %>%
  +     group_by(type, protein) %>%
  +     summarise(
    +         n_peptides = n(),
    +         mean_prob = mean(probability, na.rm = TRUE),
    +         signal_score = n_peptides * mean_prob,
    +         .groups = "drop"
    +     )
> 
  > top_signal <- protein_signal %>%
  +     group_by(type) %>%
  +     slice_max(signal_score, n = 15)
> 
  > p_signal <- ggplot(top_signal,
                       +                    aes(reorder(protein, signal_score),
                                                +                        signal_score,
                                                +                        fill = type)) +
  +     geom_col() +
  +     coord_flip() +
  +     theme_minimal() +
  +     labs(title = "Protein Signal Strength Index")
> 
  > ggsave(file.path(OUTDIR,"protein_signal_strength.png"), p_signal, width = 10, height = 6)
> protein_signal <- df %>%
  +     group_by(type, protein) %>%
  +     summarise(
    +         n_peptides = n(),
    +         mean_prob = mean(probability, na.rm = TRUE),
    +         signal_score = n_peptides * mean_prob,
    +         .groups = "drop"
    +     )
> 
  > top_signal <- protein_signal %>%
  +     group_by(type) %>%
  +     slice_max(signal_score, n = 15)
> 
  > p_signal <- ggplot(top_signal,
                       +                    aes(reorder(protein, signal_score),
                                                +                        signal_score,
                                                +                        fill = type)) +
  +     geom_col() +
  +     coord_flip() +
  +     theme_minimal() +
  +     labs(title = "Protein Signal Strength Index")
> 
  > ggsave(file.path(OUTDIR,"protein_signal_strength.png"), p_signal, width = 10, height = 6)
> df$bias_score <- (df$length / max(df$length)) *
  +     (df$probability)
> 
  > p_bias <- ggplot(df, aes(bias_score, fill = type)) +
  +     geom_density(alpha = 0.4) +
  +     theme_minimal() +
  +     labs(title = "Composite Bias Score (Length × Confidence)")
> 
  > ggsave(file.path(OUTDIR,"bias_score.png"), p_bias, width = 7, height = 5)
> pkgs_net <- c("igraph", "ggraph", "tidygraph")
> 
  > install_if_missing(pkgs_net)

Attaching package: ‘igraph’

The following objects are masked from ‘package:lubridate’:
  
  %--%, union

The following objects are masked from ‘package:dplyr’:
  
  as_data_frame, groups, union

The following objects are masked from ‘package:purrr’:
  
  compose, simplify

The following object is masked from ‘package:tidyr’:
  
  crossing

The following object is masked from ‘package:tibble’:
  
  as_data_frame

The following objects are masked from ‘package:stats’:
  
  decompose, spectrum

The following object is masked from ‘package:base’:
  
  union

WARNING: Rtools is required to build R packages but is not currently installed. Please download and install the appropriate version of Rtools before proceeding:
  
  https://cran.rstudio.com/bin/windows/Rtools/
  Installing package into ‘C:/Users/ZHUK/AppData/Local/R/win-library/4.6’
(as ‘lib’ is unspecified)
also installing the dependencies ‘proxy’, ‘coda’, ‘lpSolve’, ‘e1071’, ‘wk’, ‘geometries’, ‘tweenr’, ‘polyclip’, ‘statnet.common’, ‘gclus’, ‘qap’, ‘registry’, ‘TSP’, ‘transformr’, ‘classInt’, ‘s2’, ‘units’, ‘lwgeom’, ‘sfheaders’, ‘ggforce’, ‘tidygraph’, ‘graphlayouts’, ‘network’, ‘seriation’, ‘deldir’, ‘gganimate’, ‘sf’, ‘sfnetworks’
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/proxy_0.4-29.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/coda_0.19-4.1.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/lpSolve_5.6.23.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/e1071_1.7-17.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/wk_0.9.5.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/geometries_0.2.5.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/tweenr_2.0.3.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/polyclip_1.10-7.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/statnet.common_4.13.0.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/gclus_1.3.3.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/qap_0.1-2.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/registry_0.5-1.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/TSP_1.2.7.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/transformr_0.1.5.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/classInt_0.4-11.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/s2_1.1.11.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/units_1.0-1.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/lwgeom_0.2-16.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/sfheaders_0.4.5.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/ggforce_0.5.0.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/tidygraph_1.3.1.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/graphlayouts_1.2.3.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/network_1.20.0.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/seriation_1.5.8.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/deldir_2.0-4.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/gganimate_1.0.11.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/sf_1.1-1.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/sfnetworks_0.6.6.zip'
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/ggraph_2.2.2.zip'
package ‘proxy’ successfully unpacked and MD5 sums checked
package ‘coda’ successfully unpacked and MD5 sums checked
package ‘lpSolve’ successfully unpacked and MD5 sums checked
package ‘e1071’ successfully unpacked and MD5 sums checked
package ‘wk’ successfully unpacked and MD5 sums checked
package ‘geometries’ successfully unpacked and MD5 sums checked
package ‘tweenr’ successfully unpacked and MD5 sums checked
package ‘polyclip’ successfully unpacked and MD5 sums checked
package ‘statnet.common’ successfully unpacked and MD5 sums checked
package ‘gclus’ successfully unpacked and MD5 sums checked
package ‘qap’ successfully unpacked and MD5 sums checked
package ‘registry’ successfully unpacked and MD5 sums checked
package ‘TSP’ successfully unpacked and MD5 sums checked
package ‘transformr’ successfully unpacked and MD5 sums checked
package ‘classInt’ successfully unpacked and MD5 sums checked
package ‘s2’ successfully unpacked and MD5 sums checked
package ‘units’ successfully unpacked and MD5 sums checked
package ‘lwgeom’ successfully unpacked and MD5 sums checked
package ‘sfheaders’ successfully unpacked and MD5 sums checked
package ‘ggforce’ successfully unpacked and MD5 sums checked
package ‘tidygraph’ successfully unpacked and MD5 sums checked
package ‘graphlayouts’ successfully unpacked and MD5 sums checked
package ‘network’ successfully unpacked and MD5 sums checked
package ‘seriation’ successfully unpacked and MD5 sums checked
package ‘deldir’ successfully unpacked and MD5 sums checked
package ‘gganimate’ successfully unpacked and MD5 sums checked
package ‘sf’ successfully unpacked and MD5 sums checked
package ‘sfnetworks’ successfully unpacked and MD5 sums checked
package ‘ggraph’ successfully unpacked and MD5 sums checked

The downloaded binary packages are in
C:\Users\ZHUK\AppData\Local\Temp\Rtmp4EPb2o\downloaded_packages

Attaching package: ‘tidygraph’

The following object is masked from ‘package:igraph’:
  
  groups

The following object is masked from ‘package:stats’:
  
  filter
> library(igraph)
> library(ggraph)
> library(tidygraph)
> network_df <- df %>%
  +     select(protein, sequence, type, probability) %>%
  +     filter(!is.na(protein), !is.na(sequence))
> 
  > edges <- network_df %>%
  +     distinct(protein, sequence, type)
> 
  > g <- graph_from_data_frame(edges, directed = FALSE)
> V(g)$type_node <- ifelse(V(g)$name %in% edges$protein, "Protein", "Peptide")
> 
  > V(g)$class <- ifelse(V(g)$name %in% edges$sequence,
                         +                      edges$type[match(V(g)$name, edges$sequence)],
                         +                      "Protein")
> 
  > p_net <- ggraph(g, layout = "fr") +
  +     geom_edge_link(alpha = 0.2) +
  +     geom_node_point(aes(color = class, shape = type_node), size = 2) +
  +     theme_void() +
  +     labs(title = "Protein–Peptide Interaction Network")
> 
  > ggsave(file.path(OUTDIR,"network_protein_peptide.png"),
           +        p_net, width = 9, height = 7)
> pkgs_cluster <- c("cluster")
> 
  > install_if_missing(pkgs_cluster)
> library(cluster)
> cluster_df <- df %>%
  +     select(length, probability, charge_num, type) %>%
  +     na.omit()
> 
  > X <- scale(cluster_df[,1:3])
> set.seed(42)
> km <- kmeans(X, centers = 2, nstart = 25)
> 
  > cluster_df$cluster <- as.factor(km$cluster)
> sil <- silhouette(km$cluster, dist(X))
> 
  > sil_df <- data.frame(
    +     silhouette = sil[,3],
    +     type = cluster_df$type
    + )
> 
  > mean_sil <- mean(sil[,3])
> cat("Mean silhouette score:", mean_sil, "\n")
Mean silhouette score: 0.4111814 
> p_sil <- ggplot(sil_df, aes(silhouette, fill = type)) +
  +     geom_density(alpha = 0.4) +
  +     theme_minimal() +
  +     labs(title = paste0("Cluster Separation Quality (Silhouette = ",
                            +                         round(mean_sil, 3), ")"))
> 
  > ggsave(file.path(OUTDIR,"silhouette_score.png"),
           +        p_sil, width = 7, height = 5)
> list.files("C:/Project_R4/output_non_canonical/results", full.names = TRUE)
[1] "C:/Project_R4/output_non_canonical/results/merged_peptides.tsv"
[2] "C:/Project_R4/output_non_canonical/results/summary_table.tsv"  
> list.files("C:/Project_R4/output_data_4/results", full.names = TRUE)
character(0)
> file.exists("C:/Project_R4/output_data_4/combined_peptide.tsv")
[1] TRUE
> canonical <- read_tsv(
  +     "C:/Project_R4/output_data_4/combined_peptide.tsv",
  +     show_col_types = FALSE
  + )
New names:                                                                          
  • `` -> `...15`
> 
  > dim(canonical)
[1] 8871   15
> dir.create("C:/Project_R4/output_data_4/results", recursive = TRUE)
> 
  > write_tsv(canonical,
              +           "C:/Project_R4/output_data_4/results/merged_peptides.tsv")

> NONCANON_FILE <- "C:/Project_R4/output_non_canonical/combined_peptide.tsv"
> CANON_FILE    <- "C:/Project_R4/output_data_4/combined_peptide.tsv"
> 
  > ROOT <- "C:/Project_R4/output_comparison"
> dir.create(ROOT, recursive = TRUE, showWarnings = FALSE)
> 
  > OUTDIR <- file.path(ROOT, "results")
> FIGDIR <- file.path(ROOT, "figures")
> 
  > dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
> dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)
> pkgs <- c(
  +     "tidyverse", "readr", "stringr",
  +     "ggplot2", "patchwork",
  +     "FactoMineR", "factoextra",
  +     "umap"
  + )
> 
  > install_if_missing <- function(pkgs){
    +     for (p in pkgs){
      +         if (!requireNamespace(p, quietly = TRUE)) {
        +             install.packages(p)
        +         }
      +         library(p, character.only = TRUE)
      +     }
    + }
> 
  > install_if_missing(pkgs)
WARNING: Rtools is required to build R packages but is not currently installed. Please download and install the appropriate version of Rtools before proceeding:
  
  https://cran.rstudio.com/bin/windows/Rtools/
  Installing package into ‘C:/Users/ZHUK/AppData/Local/R/win-library/4.6’
(as ‘lib’ is unspecified)
trying URL 'https://cran.rstudio.com/bin/windows/contrib/4.6/patchwork_1.3.2.zip'
Content type 'application/zip' length 3352556 bytes (3.2 MB)
downloaded 3.2 MB

package ‘patchwork’ successfully unpacked and MD5 sums checked

The downloaded binary packages are in
C:\Users\ZHUK\AppData\Local\Temp\Rtmp4EPb2o\downloaded_packages
> load_peptides <- function(path, type_label){
  +     
    +     df <- read_tsv(path, show_col_types = FALSE)
    +     
      +     colnames(df) <- tolower(colnames(df))
      +     
        +     df$type <- type_label
        +     
          +     df$length <- nchar(df$sequence)
          +     df$probability <- as.numeric(df$probability)
          +     
            +     df$charge <- suppressWarnings(as.numeric(df$`charge states`))
            +     
              +     df <- df %>%
                +         mutate(
                  +             is_long = length > 30,
                  +             is_high_conf = probability > quantile(probability, 0.25, na.rm = TRUE)
                  +         )
              +     
                +     return(df)
              + }
> load_peptides <- function(path, type_label){
  +     
    +     df <- read_tsv(path, show_col_types = FALSE)
    +     
      +     colnames(df) <- tolower(colnames(df))
      +     
        +     df$type <- type_label
        +     
          +     df$length <- nchar(df$sequence)
          +     df$probability <- as.numeric(df$probability)
          +     
            +     df$charge <- suppressWarnings(as.numeric(df$`charge states`))
            +     
              +     df <- df %>%
                +         mutate(
                  +             is_long = length > 30,
                  +             is_high_conf = probability > quantile(probability, 0.25, na.rm = TRUE)
                  +         )
              +     
                +     return(df)
              + }
> noncanon <- load_peptides(NONCANON_FILE, "non_canonical")
New names:                                                                          
  • `` -> `...15`
> canonical <- load_peptides(CANON_FILE, "canonical")
New names:                                                                          
  • `` -> `...15`
> 
  > df <- bind_rows(noncanon, canonical)
> qc_summary <- df %>%
  +     group_by(type) %>%
  +     summarise(
    +         n = n(),
    +         mean_length = mean(length),
    +         mean_prob = mean(probability),
    +         mean_charge = mean(charge, na.rm = TRUE)
    +     )
> 
  > write_tsv(qc_summary, file.path(OUTDIR, "qc_summary.tsv"))

> p_length <- ggplot(df, aes(length, fill = type)) +
  +     geom_density(alpha = 0.4) +
  +     theme_minimal() +
  +     labs(title = "Peptide Length Distribution (Canonical vs Non-canonical)")
> 
  > ggsave(file.path(FIGDIR, "length_dist.png"), p_length)
Saving 7 x 7 in image
> feat <- df %>%
  +     select(length, probability, charge) %>%
  +     mutate(across(everything(), as.numeric)) %>%
  +     na.omit()
> 
  > pca <- prcomp(scale(feat))
> 
  > pca_df <- data.frame(
    +     PC1 = pca$x[,1],
    +     PC2 = pca$x[,2],
    +     type = df$type[as.numeric(rownames(feat))]
    + )
> 
  > p_pca <- ggplot(pca_df, aes(PC1, PC2, color = type)) +
  +     geom_point(alpha = 0.5) +
  +     theme_minimal() +
  +     labs(title = "PCA: Comparative Peptide Space")
> 
  > ggsave(file.path(FIGDIR, "pca.png"), p_pca)
Saving 7 x 7 in image
> set.seed(1)
> um <- umap(scale(feat))
Warning message:
  failed creating initial embedding; using random embedding instead 
> 
  > um_df <- data.frame(
    +     UMAP1 = um$layout[,1],
    +     UMAP2 = um$layout[,2],
    +     type = df$type[as.numeric(rownames(feat))]
    + )
> 
  > p_umap <- ggplot(um_df, aes(UMAP1, UMAP2, color = type)) +
  +     geom_point(alpha = 0.5) +
  +     theme_minimal() +
  +     labs(title = "UMAP: Structural Separation")
> 
  > ggsave(file.path(FIGDIR, "umap.png"), p_umap)
Saving 7 x 7 in image
> protein_signal <- df %>%
  +     group_by(type, protein) %>%
  +     summarise(
    +         n_peptides = n(),
    +         signal = n_peptides * mean(probability, na.rm = TRUE),
    +         .groups = "drop"
    +     )
> 
  > top <- protein_signal %>%
  +     group_by(type) %>%
  +     slice_max(signal, n = 10)
> 
  > p_signal <- ggplot(top,
                       +                    aes(reorder(protein, signal),
                                                +                        signal, fill = type)) +
  +     geom_col() +
  +     coord_flip() +
  +     theme_minimal() +
  +     labs(title = "Top Protein Signal Contribution")
> 
  > ggsave(file.path(FIGDIR, "protein_signal.png"), p_signal)
Saving 7 x 7 in image
> df$bias <- (df$length / max(df$length)) * df$probability
> 
  > p_bias <- ggplot(df, aes(bias, fill = type)) +
  +     geom_density(alpha = 0.4) +
  +     theme_minimal() +
  +     labs(title = "Composite Structural Bias Score")
> 
  > ggsave(file.path(FIGDIR, "bias.png"), p_bias)
Saving 7 x 7 in image
> ks_length <- wilcox.test(length ~ type, data = df)
> ks_prob <- wilcox.test(probability ~ type, data = df)
> 
  > stats <- data.frame(
    +     test = c("length", "probability"),
    +     p_value = c(ks_length$p.value, ks_prob$p.value)
    + )
> 
  > write_tsv(stats, file.path(OUTDIR, "stats_tests.tsv"))

> NONCANON_FILE <- "C:/Project_R4/output_non_canonical/combined_peptide.tsv"
> CANON_FILE    <- "C:/Project_R4/output_data_4/combined_peptide.tsv"
> 
  > ROOT <- "C:/Project_R4/output_comparison_v2"
> OUTDIR <- file.path(ROOT, "results")
> FIGDIR <- file.path(ROOT, "figures")
> 
  > dir.create(OUTDIR, recursive = TRUE, showWarnings = FALSE)
> dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)
> pkgs <- c(
  +     "tidyverse", "readr", "stringr",
  +     "igraph", "ggraph",
  +     "ggplot2", "patchwork"
  + )
> 
  > for (p in pkgs) {
    +     if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
    +     library(p, character.only = TRUE)
    + }
> load_df <- function(path, label){
  +     df <- read_tsv(path, show_col_types = FALSE)
  +     colnames(df) <- tolower(colnames(df))
  +     
    +     df$type <- label
    +     df$length <- nchar(df$sequence)
    +     df$probability <- as.numeric(df$probability)
    +     
      +     return(df)
    + }
> 
  > noncanon <- load_df(NONCANON_FILE, "non_canonical")
New names:                                                                          
  • `` -> `...15`
> canonical <- load_df(CANON_FILE, "canonical")
New names:                                                                          
  • `` -> `...15`
> 
  > df <- bind_rows(noncanon, canonical)
> edges <- df %>%
  +     group_by(type, protein) %>%
  +     summarise(weight = n(), .groups = "drop")
> 
  > g <- graph_from_data_frame(edges)
> 
  > set.seed(1)
> p_net <- ggraph(g, layout = "fr") +
  +     geom_edge_link(aes(width = weight), alpha = 0.4) +
  +     geom_node_point(size = 3) +
  +     theme_void() +
  +     ggtitle("Protein–Peptide Interaction Network")
> 
  > ggsave(file.path(FIGDIR, "network.png"), p_net, width = 8, height = 6)
Warning message:
  The `trans` argument of `continuous_scale()` is deprecated as of ggplot2 3.5.0.
ℹ Please use the `transform` argument instead.
This warning is displayed once per session.
Call lifecycle::last_lifecycle_warnings() to see where this warning was generated. 

> shannon_entropy <- function(seq){
  +     aa <- strsplit(seq, "")[[1]]
  +     p <- table(aa) / length(aa)
  +     -sum(p * log2(p))
  + }
> 
  > df$entropy <- sapply(df$sequence, shannon_entropy)
> 
  > p_entropy <- ggplot(df, aes(entropy, fill = type)) +
  +     geom_density(alpha = 0.4) +
  +     theme_minimal() +
  +     ggtitle("Sequence Complexity (Shannon Entropy)")
> 
  > ggsave(file.path(FIGDIR, "entropy.png"), p_entropy)
Saving 7 x 7 in image
> p_phase <- ggplot(df, aes(length, probability, color = type)) +
  +     geom_point(alpha = 0.5) +
  +     theme_minimal() +
  +     ggtitle("Charge/Length–Probability Phase Space")
> 
  > ggsave(file.path(FIGDIR, "phase_space.png"), p_phase)
Saving 7 x 7 in image
> p_phase <- ggplot(df, aes(length, probability, color = type)) +
  +     geom_point(alpha = 0.5) +
  +     theme_minimal() +
  +     ggtitle("Charge/Length–Probability Phase Space")
> 
  > ggsave(file.path(FIGDIR, "phase_space.png"), p_phase)
Saving 7 x 7 in image
> aa_df <- df %>%
  +     mutate(seq = strsplit(sequence, "")) %>%
  +     unnest(seq) %>%
  +     group_by(type, seq) %>%
  +     summarise(freq = n(), .groups = "drop")
> 
  > top_aa <- aa_df %>%
  +     group_by(type) %>%
  +     slice_max(freq, n = 10)
> 
  > p_aa <- ggplot(top_aa, aes(reorder(seq, freq), freq, fill = type)) +
  +     geom_col() +
  +     coord_flip() +
  +     theme_minimal() +
  +     ggtitle("Amino Acid Usage Bias")
> 
  > ggsave(file.path(FIGDIR, "aa_bias.png"), p_aa)
Saving 7 x 7 in image
> hub <- df %>%
  +     group_by(protein, type) %>%
  +     summarise(n = n(), .groups = "drop")
> 
  > p_hub <- ggplot(hub, aes(n, fill = type)) +
  +     geom_histogram(bins = 30, alpha = 0.5, position = "identity") +
  +     theme_minimal() +
  +     ggtitle("Protein Hub Distribution")
> 
  > ggsave(file.path(FIGDIR, "protein_hubs.png"), p_hub)
Saving 7 x 7 in image
> library(patchwork)
> 
  > master <- (p_entropy | p_phase) /
  +     (p_aa | p_hub) /
  +     (p_net) +
  +     plot_annotation(
    +         title = "Multi-layer Structural Comparison of Peptide Populations"
    +     )
> 
  > ggsave(file.path(FIGDIR, "MASTER_NEXT_LEVEL.png"),
           +        master, width = 12, height = 10)
> ############################################################
> ## 0. SETUP
  > ############################################################
> 
  > NONCAN_FILE <- "C:/Project_R4/output_non_canonical/combined_peptide.tsv"
> CAN_FILE    <- "C:/Project_R4/output_data_4/combined_peptide.tsv"
> 
  > OUTDIR <- "C:/Project_R4/biological_validation_results"
> FIGDIR <- file.path(OUTDIR, "figures")
> TABDIR <- file.path(OUTDIR, "tables")
> 
  > dir.create(FIGDIR, recursive = TRUE, showWarnings = FALSE)
> dir.create(TABDIR, recursive = TRUE, showWarnings = FALSE)
> 
  > ############################################################
> ## 1. PACKAGES
  > ############################################################
> 
  > pkgs <- c("tidyverse", "readr", "stringr")
> 
  > for (p in pkgs) {
    +     if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
    +     library(p, character.only = TRUE)
    + }
> 
  > ############################################################
> ## 2. LOAD DATA
  > ############################################################
> 
  > load_data <- function(path, label){
    +     
      +     df <- read_tsv(path, show_col_types = FALSE)
      +     colnames(df) <- tolower(colnames(df))
      +     
        +     df$type <- label
        +     
          +     # standardize expected columns safely
          +     df$sequence <- as.character(df$sequence)
          +     df$length <- nchar(df$sequence)
          +     
            +     # optional columns
            +     if ("probability" %in% colnames(df)) {
              +         df$probability <- as.numeric(df$probability)
              +     } else df$probability <- NA
              +     
                  +     if ("charge states" %in% colnames(df)) {
                    +         df$charge <- df[["charge states"]]
                    +     } else df$charge <- NA
                    +     
                        +     if ("protein" %in% colnames(df)) {
                          +         df$protein <- as.character(df$protein)
                          +     } else df$protein <- NA
                          +     
                              +     if ("grlcl_1d_ethcd spectral count" %in% colnames(df)) {
                                +         df$spectral <- df[["grlcl_1d_ethcd spectral count"]]
                                +     } else df$spectral <- NA
                                +     
                                    +     return(df)
                                + }
> 
  > can <- load_data(CAN_FILE, "canonical")
New names:                                                                          
  • `` -> `...15`
> noncan <- load_data(NONCAN_FILE, "non_canonical")
New names:                                                                          
  • `` -> `...15`
> 
  > df <- bind_rows(can, noncan)
> 
  > ############################################################
> ## 3. CLEANING
  > ############################################################
> 
  > df <- df %>%
  +     filter(!is.na(sequence)) %>%
  +     mutate(sequence = str_replace_all(sequence, "[^A-Z]", ""))
> 
  > ############################################################
> ## ========================================================
> ## BIOLOGICAL HYPOTHESIS 1:
  > ## detection bias (probability)
  > ## ========================================================
> ############################################################
> 
  > p1 <- ggplot(df, aes(probability, fill = type)) +
  +     geom_density(alpha = 0.4) +
  +     theme_minimal() +
  +     ggtitle("H1: Identification Confidence Distribution")
> 
  > ggsave(file.path(FIGDIR, "H1_probability.png"), p1)
Saving 7 x 7 in image
> 
  > ############################################################
> ## HYPOTHESIS 2:
  > ## peptide length constraint (proteasome/MHC signature)
  > ############################################################
> 
  > p2 <- ggplot(df, aes(type, length, fill = type)) +
  +     geom_boxplot(alpha = 0.6) +
  +     theme_minimal() +
  +     ggtitle("H2: Peptide Length Constraints")
> 
  > ggsave(file.path(FIGDIR, "H2_length.png"), p2)
Saving 7 x 7 in image
> 
  > ############################################################
> ## HYPOTHESIS 3:
  > ## protein origin bias (hub structure)
  > ############################################################
> 
  > hub <- df %>%
  +     group_by(type, protein) %>%
  +     summarise(n = n(), .groups = "drop")
> 
  > p3 <- ggplot(hub, aes(n, fill = type)) +
  +     geom_histogram(bins = 40, alpha = 0.5, position = "identity") +
  +     theme_minimal() +
  +     ggtitle("H3: Protein Origin (Hub Bias)")
> 
  > ggsave(file.path(FIGDIR, "H3_hubs.png"), p3)
Saving 7 x 7 in image
> 
  > ############################################################
> ## HYPOTHESIS 4:
  > ## amino acid composition bias
  > ############################################################
> 
  > aa <- df %>%
  +     mutate(aa = strsplit(sequence, "")) %>%
  +     unnest(aa) %>%
  +     group_by(type, aa) %>%
  +     summarise(n = n(), .groups = "drop")
> 
  > p4 <- ggplot(aa, aes(reorder(aa, n), n, fill = type)) +
  +     geom_col() +
  +     coord_flip() +
  +     theme_minimal() +
  +     ggtitle("H4: Amino Acid Composition Bias")
> 
  > ggsave(file.path(FIGDIR, "H4_amino_acids.png"), p4)
Saving 7 x 7 in image
> 
  > ############################################################
> ## HYPOTHESIS 5:
  > ## charge state distribution
  > ############################################################
> 
  > if (all(!is.na(df$charge))) {
    +     
      +     p5 <- ggplot(df, aes(charge, fill = type)) +
        +         geom_density(alpha = 0.4) +
        +         theme_minimal() +
        +         ggtitle("H5: Charge State Distribution")
      +     
        +     ggsave(file.path(FIGDIR, "H5_charge.png"), p5)
      +     
        + } else {
          +     message("Charge column missing or empty — skipping H5 plot")
          + }
Saving 7 x 7 in image
> 
  > ############################################################
> ## SUMMARY TABLE
  > ############################################################
> 
  > summary_table <- df %>%
  +     group_by(type) %>%
  +     summarise(
    +         n_peptides = n(),
    +         mean_length = mean(length, na.rm = TRUE),
    +         mean_prob = mean(probability, na.rm = TRUE),
    +         median_length = median(length, na.rm = TRUE)
    +     )
> 
  > write_tsv(summary_table, file.path(TABDIR, "summary_table.tsv"))
>                                                                                 
  > ############################################################
> ## FINAL MESSAGE
  > ############################################################
> 
  > cat("DONE: Biological validation pipeline completed\n")
DONE: Biological validation pipeline completed
> cat("Figures saved in:", FIGDIR, "\n")
Figures saved in: C:/Project_R4/biological_validation_results/figures 
> cat("Tables saved in:", TABDIR, "\n")
Tables saved in: C:/Project_R4/biological_validation_results/tables