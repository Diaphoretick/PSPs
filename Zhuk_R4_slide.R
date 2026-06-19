############################################################
# MS/MS REANALYSIS PIPELINE (GitHub CLEAN VERSION)
############################################################

library(data.table)
library(ggplot2)
library(Biostrings)

base_path <- "C:/Project_R4"

out_dir <- file.path(base_path, "out")
fig_dir <- file.path(base_path, "figs")

dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)
dir.create(fig_dir, showWarnings=FALSE, recursive=TRUE)

############################################################
# LOAD DATA
############################################################

df <- fread(file.path(base_path, "out/result.tsv"))

############################################################
# QC
############################################################

df[, is_decoy := grepl("^XXX_", Protein)]

df[, pep := Peptide]
df[, pep := gsub("^[A-Z]\\.", "", pep)]
df[, pep := gsub("\\.[A-Z]$", "", pep)]
df[, pep := gsub("\\+[^A-Z]*", "", pep)]

df[, len := nchar(pep)]

contam <- c("KRT","KERATIN","ALB","TRYPSIN")

df[, contaminant := Reduce(`|`, lapply(contam, function(x){
  grepl(x, Protein, ignore.case=TRUE)
}))]

############################################################
# SUMMARY METRICS
############################################################

summary_dt <- df[, .(
  total_psms = .N,
  decoy_fraction = mean(is_decoy),
  contaminant_fraction = mean(contaminant),
  hla_fraction = mean(len >= 8 & len <= 11)
)]

fwrite(summary_dt, file.path(out_dir, "summary.tsv"))

############################################################
# FIGURE 1 (MAIN SLIDE FIGURE)
############################################################

png(file.path(fig_dir, "Figure1.png"), 1400, 500)
par(mfrow=c(1,3))

hist(df$len, main="Peptide length", col="grey")
hist(df$MSGFScore, main="MSGFScore", col="steelblue")
barplot(table(df$is_decoy), main="Decoy ratio", col="tomato")

dev.off()

############################################################
# GENCODE (FAST VERSION)
############################################################

gencode <- readAAStringSet(
  "https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_44/gencode.v44.pc_translations.fa.gz"
)

proteome <- as.character(gencode)

df[, in_gencode := pep %in% proteome]

fwrite(df[, .(pep, Protein, in_gencode)],
       file.path(out_dir, "gencode_matches.tsv"))

############################################################
# DONE
############################################################

cat("PIPELINE COMPLETE\n")