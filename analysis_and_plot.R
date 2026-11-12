# Dog neuronal somatic-mutation analysis
# =======================================
#
# This script reproduces the statistical analyses and figures associated with
# the manuscript. It is intentionally kept as a single, sequential analysis so
# that the order used for the archived results remains explicit.
#
# Expected repository layout:
#
#   dog_project/
#   |-- code/
#   |   |-- dog_mutation_project.R
#   |   `-- dog_mutation_helpers.R
#   |-- tables/
#   `-- figures/
#
# Human comparison data and the aging-project reference data are expected in
# the relative locations used below. Run the complete analysis with:
#
#   Rscript code/dog_mutation_project.R
#
# An alternative project root can be supplied when the data are elsewhere:
#
#   Rscript code/dog_mutation_project.R --project-root=/path/to/dog_project
#
# The script overwrites manuscript tables, RDS objects, and PDF figures with
# the same names as the archived analysis. Run it in a clean output directory
# when preservation of existing results is required.

# Configuration and dependencies ---------------------------------------------
command_args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", command_args, value = TRUE)
script_dir <- if (length(file_arg) == 1L) {
    dirname(normalizePath(sub("^--file=", "", file_arg), mustWork = TRUE))
} else {
    normalizePath(getwd(), mustWork = TRUE)
}

analysis_args <- commandArgs(trailingOnly = TRUE)
root_arg <- grep("^--project-root=", analysis_args, value = TRUE)
default_root <- if (basename(script_dir) == "code") dirname(script_dir) else script_dir
project_root <- if (length(root_arg) == 1L) sub("^--project-root=", "", root_arg) else default_root
project_root <- normalizePath(project_root, mustWork = TRUE)

if (!dir.exists(file.path(project_root, "tables"))) {
    stop("Project root does not contain a 'tables' directory: ", project_root)
}
setwd(project_root)

source(file.path(script_dir, "dog_mutation_helpers.R"))

required_packages <- c(
    "DESeq2", "epitools", "lmerTest", "MASS", "MutationalPatterns",
    "NMF", "pheatmap", "RColorBrewer"
)
missing_packages <- required_packages[
    !vapply(required_packages, requireNamespace, logical(1), quietly = TRUE)
]
if (length(missing_packages) > 0L) {
    stop(
        "Install the following R packages before running this analysis: ",
        paste(missing_packages, collapse = ", ")
    )
}

suppressPackageStartupMessages({
    library(epitools)
    library(lmerTest)
    library(MASS)
    library(pheatmap)
    library(RColorBrewer)
})

output_directories <- c(
    "figures", "figures/CtoT_in_CpG", "figures/cosmic_signatures",
    "figures/enrichment_analysis", "figures/enrichment_analysis/with_crypt",
    "figures/mutation_profile_human_vs_dog",
    "tables/For_manuscript"
)
invisible(lapply(output_directories, dir.create, recursive = TRUE, showWarnings = FALSE))
set.seed(1)

# Input data and analysis constants -------------------------------------------
# Mutation burdens are divided by 5.833 Gb, the length of the human autosomes,
# to report comparable per-Gb values. Dog burdens were normalized upstream.
# Dog autosomes span 4.408 Gb. Sex-specific analyses are intentionally avoided
# because sex annotations in the preliminary metadata were unreliable; the
# final metadata table is used here.
len_dog_autosome <- 4.408
len_human_autosome <- 5.833  # Autosomes only; unplaced scaffolds are excluded.
samples <- read.table("tables/annotation/metadata.txt", header = TRUE, row.names = 1, comment.char = "",
    sep = "\t", check.names = FALSE)
samples_meta <- read.csv("tables/annotation/scdnaseq.metadata_final.txt", header = TRUE, row.names = 1,
    colClasses = "character", comment.char = "", sep = "\t", check.names = FALSE)
samples_meta$age <- as.numeric(samples_meta$age)
snv_burden <- read.csv("tables/dnaseq/scan2/master.burden_snv.csv", header = T, row.names = 1)
indel_burden <- read.csv("tables/dnaseq/scan2/master.burden_indel.csv", header = T, row.names = 1)
snv.sig <- read.table("tables/dnaseq/scan2/master.snv_profiles", header = T, row.names = 1, check.names = F)
indel.sig <- read.table("tables/dnaseq/scan2/master.indel_profiles", header = T, row.names = 1, check.names = F)
# needs re-do with new data
nmf_res_snv <- readRDS("tables/dnaseq/scan2/signatures.snv.rds")
nmf_res_indel <- readRDS("tables/dnaseq/scan2/signatures.indel.rds")
nmf_res_snv_with_crypt <- readRDS("tables/dnaseq/scan2/signatures.snv_with_crypt.rds")
nmf_res_indel_with_crypt <- readRDS("tables/dnaseq/scan2/signatures.indel_with_crypt.rds")

cs_species <- c(dog = "#e41a1c", human = "#1a81e4")
cs_breed <- c(Belgian_Malinois = "#e41a1c", Boxer = "#377eb8", Golden_Retriever = "#4daf4a", Hound_Cross = "#984ea3",
    Maltese_Cross = "#ff7f00", Mixed_Breed = cs_gg[4], Norwich_Terrier = "#a65628", Pitbull_Cross = "#f781bf",
    Pug_Cross = "#999999", Shih_Tzu = "#8dd3c7", Labrador_Retriever = "#a6cee3", German_Shepard_x_catxho = "#ffd92f",
    Dachshund = "#66c2a5")
cs_sex <- c(Female = cs_gg[1], Male = cs_gg[2], unknown = "grey")

dog_aging.exp <- read.table("tables/rnaseq/dog_aging.gene.tpm", header = T, row.names = 1)[, c(4, 5,
    6, 7, 8, 9, 10, 11, 1, 2, 3)]
dog_aging.count <- read.table("tables/rnaseq/dog_aging.gene.count", header = T, row.names = 1)[, c(4,
    5, 6, 7, 8, 9, 10, 11, 1, 2, 3)]
dog_aging.meta <- read.table("tables/rnaseq/dog_aging.metadata.txt", header = F, row.names = 2)[-1, ]
dog_aging.young <- row.names(dog_aging.meta)[dog_aging.meta[, 3] < 10]
dog_aging.old <- row.names(dog_aging.meta)[dog_aging.meta[, 3] >= 10]
gtex.exp <- read.table("tables/rnaseq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct",
    header = T, row.names = 1, sep = "\t", check.names = F)

tri.context.dog <- read.table("tables/annotation/canFam3.96context.txt", header = F, row.names = 1)
tri.context.human <- read.table("tables/annotation/hs37d5.96context.txt", header = F, row.names = 1)
tri.factors <- read.table("tables/annotation/trinucleotide.factors.txt", header = F, row.names = 1)  # divide dog mutation signature by this factor to correct to human signatures
# correct for trinucleotide context for dog to human
snv.sig <- snv.sig/tri.factors[, 1]

# first remove those suspicious neurons with more than 1000 estimated sSNVs
cells_used <- row.names(snv_burden)[snv_burden$trimmed.genome.burden < 1000]
# also remove the three dogs that have overburden signature A2 mutations (C>T)
# Brad said they are from the same batch, and these overburden might be temperature related
# we have re-did the sequencing for the three dogs now
cells_used <- setdiff(cells_used, c("2206A3231207_1.fastq.gz", "2206A4231207_1.fastq.gz", "2206A5231207_1.fastq.gz",
    "2206A6231207_1.fastq.gz", "2205C3231207_1.fastq.gz", "2205C4231207_1.fastq.gz", "2305D1231207_1.fastq.gz",
    "2305D2231207_1.fastq.gz", "2305D3231207_1.fastq.gz", "2305D4231207_1.fastq.gz"))
# several neurons are problematic
# 2305B4251119 has a very low sSNV burden due to super low coverage
cells_used <- setdiff(cells_used, "2305B4251119")

cells_filtered <- setdiff(row.names(snv_burden), cells_used)
cells_dog <- cells_used
cells_dog_highCov <- cells_dog[as.numeric(samples_meta[cells_dog, ]$coverage) >= 20]
# a big matrix for snv and indel burden and metainfo
out_mat <- cbind(samples_meta[cells_used, ], snv_burden[cells_used, "genome.burden"], indel_burden[cells_used,
    "genome.burden"])
colnames(out_mat)[c(ncol(out_mat) - 1, ncol(out_mat))] <- c("sSNV_burden", "sINDEL_burden")

cs_sig <- c("#56B3D9", "#0A0B09", "#C13D35", "#C6C3C0", "#A4C56E", "#E1C3BC")
col_sig <- c()
for (i in 1:6) {
    col_sig <- c(col_sig, rep(cs_sig[i], 16))
}
col_sig_indel <- c(rep(rgb(244, 192, 122, maxColorValue = 255), 6), rep(rgb(238, 133, 51, maxColorValue = 255),
    6), rep(rgb(182, 217, 145, maxColorValue = 255), 6), rep(rgb(87, 157, 64, maxColorValue = 255), 6),
    rep(rgb(243, 203, 183, maxColorValue = 255), 6), rep(rgb(237, 142, 112, maxColorValue = 255), 6),
    rep(rgb(221, 83, 64, maxColorValue = 255), 6), rep(rgb(172, 45, 38, maxColorValue = 255), 6), rep(rgb(211,
        223, 239, maxColorValue = 255), 6), rep(rgb(156, 193, 222, maxColorValue = 255), 6), rep(rgb(95,
        150, 197, maxColorValue = 255), 6), rep(rgb(46, 98, 165, maxColorValue = 255), 6), rep(rgb(225,
        223, 239, maxColorValue = 255), 1), rep(rgb(181, 184, 213, maxColorValue = 255), 2), rep(rgb(134,
        132, 187, maxColorValue = 255), 3), rep(rgb(91, 67, 148, maxColorValue = 255), 5))
col_sbs <- c("#66c2a5", "#fc8d62", "#80b1d3")
names(col_sbs) <- c("signature A1", "signature A2", "signature A3")
col_id <- c("#8da0cb", "#e78ac3", "#b3de69")
names(col_id) <- c("signature ID-A", "signature ID-B", "signature ID-C")

# out_mat_all$species='dog';out_mat_all[cells_hs,'species']='human';out_mat_all$species=as.factor(out_mat_all$species)
# write.table(out_mat_all,'tables/For_manuscript/summary.dog+human.txt',sep='\\t',quote=F)
out_mat_all <- read.table("tables/For_manuscript/summary.dog+human.txt", header = T, row.names = 1, sep = "\t",
    check.names = F)
out_mat_all[, "species"] <- as.factor(out_mat_all[, "species"])
# aging DEGs from bulk tissue for human and dog
deseq_res <- read.table("tables/rnaseq/dog_aging.deseq.txt", header = T, row.names = 1)
deseq_downgenes <- row.names(deseq_res)[which(deseq_res$log2FoldChange < (-0.585) & deseq_res$padj <
    0.05)]
deseq_upgenes <- row.names(deseq_res)[which(deseq_res$log2FoldChange > (0.585) & deseq_res$padj < 0.05)]
deseq_downgenes2 <- row.names(deseq_res)[which(deseq_res$log2FoldChange < 0 & deseq_res$padj < 0.05)]
deseq_upgenes2 <- row.names(deseq_res)[which(deseq_res$log2FoldChange > 0 & deseq_res$padj < 0.05)]
deseq_gtex <- read.table("tables/rnaseq/GTEx.deseq.txt", header = T, row.names = 1)
deseq_gtex_downgenes <- row.names(deseq_gtex)[which(deseq_gtex$log2FoldChange < (-0.585) & deseq_gtex$padj <
    0.05)]
deseq_gtex_upgenes <- row.names(deseq_gtex)[which(deseq_gtex$log2FoldChange > (0.585) & deseq_gtex$padj <
    0.05)]
# aging DEGs from snRNAseq in our aging paper in excitatory neurons
degs <- read.table("../aging_project/tables/rnaseq/pfc.clean.degs_filterInfants.separtedByCellTypes.txt",
    header = T, row.names = 1)
ct <- "ext"
degs.up.human <- degs[which(degs[, 6] == "elderly" & apply(degs[, 3:4], 1, max) > 0.25 & degs[, 2] >
    0.5 & degs[, 8] == ct), 7]
degs.down.human <- degs[which(degs[, 6] == "elderly" & apply(degs[, 3:4], 1, max) > 0.25 & degs[, 2] <
    (-0.5) & degs[, 8] == ct), 7]
# aging DEGs from snRNAseq from the NSR paper in ext neurons for dog
degs.up.dog <- unique(as.vector(read.table("tables/rnaseq/nwaf388_NSR_dogAging_snRNA/ext_up.txt", header = F,
    row.names = NULL)[, 1]))
degs.down.dog <- unique(as.vector(read.table("tables/rnaseq/nwaf388_NSR_dogAging_snRNA/ext_down.txt",
    header = F, row.names = NULL)[, 1]))

# gene length
gene.len.dog <- read.table("tables/annotation/canFam3.mRNA.genelen", header = T, row.names = 1)
gene.len.human <- read.table("../aging_project/tables/annotation/hs37.gene.len", header = T, row.names = 1)
gene.len.human <- gene.len.human[gene.len.human$gene_type == "protein_coding", ]

# read crypt data
snv.sig.crypt <- read.table("tables/dnaseq/cagan/master.snv_profiles", header = T, row.names = 1)
indel.sig.crypt <- read.table("tables/dnaseq/cagan/master.indel_profiles", header = T, row.names = 1)
meta_crypt <- read.table("tables/annotation/cagan_et_al.metadata.txt", header = T, row.names = 1, sep = "\t")
meta_crypt <- meta_crypt[meta_crypt$Species %in% c("Dog", "Human"), ]
cells_crypt_dog <- row.names(meta_crypt)[meta_crypt$Species == "Dog"]
cells_crypt_dog <- cells_crypt_dog[order(meta_crypt[cells_crypt_dog, "Age..years."])]
cells_crypt_human <- row.names(meta_crypt)[meta_crypt$Species == "Human"]
cells_crypt_human <- cells_crypt_human[order(meta_crypt[cells_crypt_human, "Age..years."])]
cells_crypt <- c(cells_crypt_dog, cells_crypt_human)

# read human mutation data from Jenn's project----
meta_pd <- read.csv("../tables/annotation/scdnaseq.metadata_final.txt", header = T, row.names = 1, check.names = F,
    sep = "\t", comment.char = "")
meta_pd$platform <- "unknown"
meta_pd[c("3142_230403A08", "14279_230403C01", "14279_230403C04"), "platform"] <- "NovaseqX_10B"
meta_pd[c("3142_230403_A2", "14279_230403_A3", "14279_230403_D3"), "platform"] <- "Novaseq_25B"
snv_pd <- read.table("../tables/dnaseq/scan2_rerun/master.burden_snv.csv", header = T, row.names = 1,
    sep = ",")
indel_pd <- read.table("../tables/dnaseq/scan2_rerun/master.burden_indel.csv", header = T, row.names = 1,
    sep = ",")
cells <- row.names(snv_pd)[meta_pd[row.names(snv_pd), "method"] == "FACS" & meta_pd[row.names(snv_pd),
    "cellNumber"] == 1]
cells_Jenn <- row.names(snv_pd)[meta_pd[row.names(snv_pd), "method"] == "FACS" & meta_pd[row.names(snv_pd),
    "cellNumber"] == 1 & meta_pd[row.names(snv_pd), "source"] %in% c("Jenn")]
cells_HMS <- row.names(snv_pd)[meta_pd[row.names(snv_pd), "method"] == "FACS" & meta_pd[row.names(snv_pd),
    "cellNumber"] == 1 & meta_pd[row.names(snv_pd), "source"] %in% c("HMS")]
cells_Jenn2 <- row.names(snv_pd)[meta_pd[row.names(snv_pd), "method"] == "FACS" & meta_pd[row.names(snv_pd),
    "cellNumber"] == 1 & meta_pd[row.names(snv_pd), "source"] %in% c("Jenn", "HMS")]
# now also include HMS normal samples
cells_Jenn <- cells_Jenn2[which(meta_pd[cells_Jenn2, ]$disease != "AD")]
cells_Jenn <- setdiff(cells_Jenn, c("14279_230403_A3", "14279_230403_D3", "3142_230403_A2"))
# remove those HMS neurons that have more than 2000 indels been called
cells_Jenn2 <- cells_Jenn2[which(indel_pd[cells_Jenn2, "genome.burden"] <= 2000)]
cells_Jenn2 <- setdiff(cells_Jenn2, c("14279_230403_A3", "14279_230403_D3", "3142_230403_A2"))
cells_normal <- cells_Jenn[meta_pd[cells_Jenn, "disease"] == "normal"]
cells_pd <- cells_Jenn[meta_pd[cells_Jenn, "disease"] == "PD"]
cells_ad <- row.names(snv_pd)[meta_pd[row.names(snv_pd), "method"] == "FACS" & meta_pd[row.names(snv_pd),
    "cellNumber"] == 1 & meta_pd[row.names(snv_pd), "source"] %in% c("Jenn", "HMS") & meta_pd[row.names(snv_pd),
    "disease"] == "AD"]
cells_ad <- cells_ad[meta_pd[cells_ad, ]$`allele_dropout_rate(%)` < 20]

snv.sig_pd <- read.table("../tables/dnaseq/scan2_rerun/master.snv_profiles", header = T, row.names = 1,
    check.names = F)
indel.sig_pd <- read.table("../tables/dnaseq/scan2_rerun/master.indel_profiles", header = T, row.names = 1,
    check.names = F)
#nmf_res_snv_hs <- readRDS('../tables/dnaseq/scan2_rerun/signatures.snv.rds')
#nmf_res_indel_hs <- readRDS('../tables/dnaseq/scan2_rerun/signatures.indel.rds')
mut_mat_indel_hs <- indel.sig_pd[, cells_Jenn] + 0.01
mut_mat_snv_hs <- snv.sig_pd[, cells_Jenn] + 0.01

samples_hs <- unique(meta_pd[, c("donor", "age", "sex", "disease")])
row.names(samples_hs) <- samples_hs[, 1]
cs_disease <- c(normal = "#2B3990", PD = "#F15A29", AD = cs_gg[1])

cells_hs_young <- cells_normal[meta_pd[cells_normal, "age"] <= 20]
#cells_normal=cells_normal[! meta_pd[cells_normal,'donor'] %in% c('5451','4976')]
cells_hs <- cells_normal

out_mat_hs <- cbind(meta_pd[cells_hs, ], snv_pd[cells_hs, "genome.burden"], indel_pd[cells_hs, "genome.burden"])
colnames(out_mat_hs)[19:20] <- c("sSNV_burden", "sINDEL_burden")
write.table(out_mat_hs, "tables/For_manuscript/summary_human.txt", quote = F, sep = "\t")


# make big metadata matrix; combine human and dog; after doing signature analysis----
nmf_res <- readRDS("tables/dnaseq/scan2/signatures.snv.rds")
t <- apply(nmf_res$contribution * colSums(nmf_res$signatures), 2, function(x) {
    return(x/sum(x))
})
n_sigs <- t(t)
n_sigs[cells_used, ] <- n_sigs[cells_used, ] * out_mat[cells_used, ]$sSNV_burden
n_sigs[cells_normal, ] <- n_sigs[cells_normal, ] * out_mat_hs[cells_normal, ]$sSNV_burden

tm1 <- out_mat[, c(1, 6:ncol(out_mat))]
tm2 <- cbind(out_mat_hs[, c(1, 6)], rep("NA", nrow(out_mat_hs)), out_mat_hs[, c(13, 14, 7, 8, 9, 10,
    11, 12, 19, 20)])
colnames(tm2)[3] <- "breed"
tmp_out <- rbind(tm1, tm2)
tmp_out[, "signature A1"] <- n_sigs[row.names(tmp_out), "signature A1"]
tmp_out[, "signature A2"] <- n_sigs[row.names(tmp_out), "signature A2"]

nmf_res <- readRDS("tables/dnaseq/scan2/signatures.indel.rds")
t <- apply(nmf_res$contribution * colSums(nmf_res$signatures), 2, function(x) {
    return(x/sum(x))
})
n_sigs <- t(t)
n_sigs[cells_used, ] <- n_sigs[cells_used, ] * out_mat[cells_used, ]$sINDEL_burden
n_sigs[cells_normal, ] <- n_sigs[cells_normal, ] * out_mat_hs[cells_normal, ]$sINDEL_burden
tmp_out[, "signature ID-A"] <- n_sigs[row.names(tmp_out), "signature ID-A"]
tmp_out[, "signature ID-B"] <- n_sigs[row.names(tmp_out), "signature ID-B"]
tmp_out[, "species"] <- "dog"
tmp_out[cells_hs, "species"] <- "human"
write.table(tmp_out, "tables/For_manuscript/summary.dog+human.txt", sep = "\t", quote = F)

# data quality----
# overall quality
dog.sc.qua <- read.table("tables/annotation/scdnaseq.metadata_final.txt", header = T, row.names = 1,
    sep = "\t", check.names = F)
dog.bk.qua <- read.table("tables/annotation/bulkdnaseq.metadata2.txt", header = T, row.names = 1, sep = "\t",
    check.names = F)
tcs <- rep(cs_gg[2], dim(dog.sc.qua)[1])
names(tcs) <- row.names(dog.sc.qua)
tcs[cells_filtered] <- "grey"

pdf("figures/data_quality.pdf", width = 12, height = 7, useDingbats = F)
par(mar = c(2, 1, 3, 1), tcl = 0.3, bty = "n", xpd = T, mfrow = c(1, 7))
plot.new()
barplot(dog.sc.qua[, "coverage"], horiz = T, xlim = c(0, 80), main = "coverage (x)", names = row.names(dog.sc.qua),
    las = 1, space = 0, border = F, col = tcs)
barplot(dog.sc.qua[, "mapping_rate(%)"], horiz = T, xlim = c(0, 100), main = "mapping rate (%)", las = 1,
    space = 0, border = F, col = tcs)
barplot(dog.sc.qua[, "duplication_rate(%)"], horiz = T, xlim = c(0, 100), main = "PCR duplication\nrate (%)",
    las = 1, space = 0, border = F, col = tcs)
barplot(dog.sc.qua[, "allele_dropout_rate(%)"], horiz = T, xlim = c(0, 70), main = "allele dropout\nrate (%)",
    las = 1, space = 0, border = F, col = tcs)
barplot(dog.sc.qua[, "locus_dropout_rate(%)"], horiz = T, xlim = c(0, 20), main = "locus dropout\nrate (%)",
    las = 1, space = 0, border = F, col = tcs)
barplot(dog.sc.qua[, "mapd"], horiz = T, xlim = c(0, 2), main = "MAPD", las = 1, space = 0, border = F,
    col = tcs)
dev.off()

# specific to mapd
fun_mapd <- function(cn, mapd, m) {
    cn <- cn[which(cn[, 1] %in% paste("chr", c(1:38, "X", "Y"), sep = "")), ]
    mapd <- mapd[which(mapd[, 1] %in% paste("chr", c(1:38, "X"), sep = "")), ]
    l <- rep(0, 40)
    names(l) <- paste("chr", c(1:38, "X", "Y"), sep = "")
    p <- 0
    for (i in paste("chr", c(1:38, "X", "Y"), sep = "")) {
        l[i] <- p + length(which(cn[, 1] == i))/100
        p <- l[i]
    }
    tmn <- tapply(cn[, 3], (seq_along(cn[, 3]) - 1)%/%100, mean)
    tsd <- tapply(cn[, 3], (seq_along(cn[, 3]) - 1)%/%100, sd)
    plot(tmn, ylim = c(-4, 5), pch = 20, yaxt = "n", xaxt = "n", xlab = "", ylab = "copy number ratio",
        main = m)
    arrows(1:length(tmn), tmn - tsd, 1:length(tmn), tmn + tsd, length = 0)
    axis(2, c(-4, 0, 1, 2, 5), label = c("", 0, 1, 2, ""))
    abline(h = c(0, 1, 2), v = l[1:40], lty = 2, col = "grey")
    axis(1, (l - c(0, l[1:39]))/2 + c(0, l[1:39]), label = names(l), lwd = 0, las = 3, cex.axis = 1)
    text(1, 4.5, label = paste("median MAPD=", round(median(mapd[, 3]), 2)), pos = 4)
    text(1, 4, label = paste("mean MAPD=", round(mean(mapd[, 3]), 2)), pos = 4)
}
pdf("figures/MAPD.pdf", width = 10, height = 4, useDingbats = F)
par(mar = c(6, 4, 3, 1), tcl = 0.3, cex = 5/6, bty = "n")
for (f in list.files("tables/dnaseq/mapd/", pattern = ".mapd")) {
    m <- strsplit(f, ".mapd")[[1]][1]
    cn <- read.table(paste("tables/dnaseq/MAPD/", m, ".cn", sep = ""), header = F, row.names = NULL)
    mapd <- read.table(paste("tables/dnaseq/MAPD/", f, sep = ""), header = F, row.names = NULL)
    fun_mapd(cn, mapd, m)
}
dev.off()

# mutation burden; also add signature burden (after signature analysis)----
# function
# for now, I build linear mixed effect models, and reported regression line and mutation accumulation rates per year
# I also added covariants (coverage, MAPD, and sex) into the model to test robustness, reported p-values
fun_mutburden <- function(out_mat_all, cells_dog, cells_human, colname, xlm, ylab) {
    ylm <- max(out_mat_all[c(cells_dog, cells_human)[which(out_mat_all[c(cells_dog, cells_human), "age"] <
        xlm)], colname]/len_human_autosome)
    ylm <- ceiling(ylm/(10^floor(log10(ylm)))) * (10^floor(log10(ylm)))
    plot(out_mat_all[cells_human, "age"], out_mat_all[cells_human, colname]/len_human_autosome, pch = 21,
        col = "white", bg = cs_species["human"], lwd = 0.5, cex = 1.5, ylab = ylab, ylim = c(0, ylm),
        xlim = c(0, xlm), xlab = "age (yrs)")
    points(out_mat_all[cells_dog, "age"], out_mat_all[cells_dog, colname]/len_human_autosome, bg = cs_species["dog"],
        col = "white", pch = 21, lwd = 0.5, cex = 1.5)
    df <- out_mat_all[c(cells_dog, cells_human), c("age", colname, "species", "donor", "sex", "coverage",
        "mapd")]
    colnames(df)[2] <- "burden"
    df[, 2] <- df[, 2]/len_human_autosome
    model_dog <- lmer(burden ~ age + sex + coverage + mapd + (1 | donor), data = df[cells_dog, ])
    model_human <- lmer(burden ~ age + sex + coverage + mapd + (1 | donor), data = df[cells_human, ])
    model_mixed_coef <- summary(lmer(burden ~ age * species + sex + coverage + mapd + (1 | donor), data = df))$coefficients
    model_mixed_varcor <- as.data.frame(summary(lmer(burden ~ age * species + sex + coverage + mapd +
        (1 | donor), data = df))$varcor)
    model_dog_nocov <- lmer(burden ~ age + (1 | donor), data = df[cells_dog, ])
    model_human_nocov <- lmer(burden ~ age + (1 | donor), data = df[cells_human, ])
    icc <- round(model_mixed_varcor[1, 4]/(model_mixed_varcor[1, 4] + model_mixed_varcor[2, 4]), 2)
    fun_addRegressionLine_lmer(model_human_nocov, "age", 0, max(out_mat_all[cells_human, "age"]), species = "human",
        se = T, lwd = 2, col = cs_species["human"])
    fun_addRegressionLine_lmer(model_dog_nocov, "age", 0, max(out_mat_all[cells_dog, "age"]), species = "dog",
        se = T, lwd = 2, col = cs_species["dog"])
    text(xlm, ylm * 0.1, pos = 2, col = cs_species["human"], cex = 0.7, label = paste(round(summary(model_human_nocov)$coefficients[2,
        1], 1), " mutation/yr/Gb; R^2=", round(cor(df[cells_human, "burden"], df[cells_human, "age"])^2,
        2), "; p=", format(summary(model_human)$coefficients[2, 5], scientific = T, digit = 2), sep = ""))
    text(xlm, ylm * 0.15, pos = 2, col = cs_species["dog"], cex = 0.7, label = paste(round(summary(model_dog_nocov)$coefficients[2,
        1], 1), " mutation/yr/Gb; R^2=", round(cor(df[cells_dog, "burden"], df[cells_dog, "age"])^2,
        2), "; p=", format(summary(model_dog)$coefficients[2, 5], scientific = T, digit = 2), sep = ""))
    text(xlm, ylm * 0.05, pos = 2, col = "black", cex = 0.7, paste("slope diff p=", format(model_mixed_coef["age:specieshuman",
        5], scientific = T, digit = 2), "; intercept diff p=", format(model_mixed_coef["specieshuman",
        5], scientific = T, digit = 2), sep = ""))
    text(xlm, ylm * 0.2, pos = 2, col = "black", cex = 0.7, label = paste("ICC(variance from donor)=",
        icc, sep = ""))
    legend("topleft", pch = c(19, 19), legend = c("dog", "human"), bty = "n", pt.cex = 1.2, col = cs_species)
}
# define young as \342\211\24420 yrs
cells_hs_young <- cells_human[meta_pd[cells_human, "age"] <= 20]
# plot
pdf("figures/mutation_burden.pdf", width = 4.5, height = 4.5, useDingbats = F)
par(mar = c(4, 4, 1, 1), tcl = 0.3, bty = "n")
fun_mutburden(out_mat_all, cells_dog, cells_hs, "sSNV_burden", 110, "sSNV/Gb")
fun_mutburden(out_mat_all, cells_dog, cells_hs, "sSNV_burden", 20, "sSNV/Gb")
fun_mutburden(out_mat_all, cells_dog, cells_hs_young, "sSNV_burden", 20, "sSNV/Gb")
fun_mutburden(out_mat_all, cells_dog, cells_hs, "sINDEL_burden", 110, "sindel/Gb")
fun_mutburden(out_mat_all, cells_dog, cells_hs, "sINDEL_burden", 20, "sindel/Gb")
fun_mutburden(out_mat_all, cells_dog, cells_hs_young, "sINDEL_burden", 20, "sindel/Gb")
dev.off()

pdf("figures/mutation_burden_signatures.pdf", width = 4.5, height = 4.5, useDingbats = F)
par(mar = c(4, 4, 1, 1), tcl = 0.3, bty = "n")
for (s in c("signature A1", "signature A2", "signature ID-A", "signature ID-B")) {
    fun_mutburden(out_mat_all, cells_dog, cells_hs, s, 110, paste(s, "/Gb", sep = ""))
    fun_mutburden(out_mat_all, cells_dog, cells_hs, s, 20, paste(s, "/Gb", sep = ""))
    fun_mutburden(out_mat_all, cells_dog, cells_hs_young, s, 20, paste(s, "/Gb", sep = ""))
}
dev.off()

# mutation burden (only for cells with coverage >= 20); also add signature burden (after signature analysis)----
# function
# for now, I build linear mixed effect models, and reported regression line and mutation accumulation rates per year
# I also added covariants (coverage, MAPD, and sex) into the model to test robustness, reported p-values
cells_human <- cells_hs
fun_mutburden <- function(out_mat_all, cells_dog, cells_human, colname, xlm, ylab) {
    ylm <- max(out_mat_all[c(cells_dog, cells_human)[which(out_mat_all[c(cells_dog, cells_human), "age"] <
        xlm)], colname]/len_human_autosome)
    ylm <- ceiling(ylm/(10^floor(log10(ylm)))) * (10^floor(log10(ylm)))
    plot(out_mat_all[cells_human, "age"], out_mat_all[cells_human, colname]/len_human_autosome, pch = 21,
        col = "white", bg = cs_species["human"], lwd = 0.5, cex = 1.5, ylab = ylab, ylim = c(0, ylm),
        xlim = c(0, xlm), xlab = "age (yrs)")
    points(out_mat_all[cells_dog, "age"], out_mat_all[cells_dog, colname]/len_human_autosome, bg = cs_species["dog"],
        col = "white", pch = 21, lwd = 0.5, cex = 1.5)
    df <- out_mat_all[c(cells_dog, cells_human), c("age", colname, "species", "donor", "sex", "coverage",
        "mapd")]
    colnames(df)[2] <- "burden"
    df[, 2] <- df[, 2]/len_human_autosome
    model_dog <- lmer(burden ~ age + sex + coverage + mapd + (1 | donor), data = df[cells_dog, ])
    model_human <- lmer(burden ~ age + sex + coverage + mapd + (1 | donor), data = df[cells_human, ])
    model_mixed_coef <- summary(lmer(burden ~ age * species + sex + coverage + mapd + (1 | donor), data = df))$coefficients
    model_mixed_varcor <- as.data.frame(summary(lmer(burden ~ age * species + sex + coverage + mapd +
        (1 | donor), data = df))$varcor)
    model_dog_nocov <- lmer(burden ~ age + (1 | donor), data = df[cells_dog, ])
    model_human_nocov <- lmer(burden ~ age + (1 | donor), data = df[cells_human, ])
    icc <- round(model_mixed_varcor[1, 4]/(model_mixed_varcor[1, 4] + model_mixed_varcor[2, 4]), 2)
    fun_addRegressionLine_lmer(model_human_nocov, "age", 0, max(out_mat_all[cells_human, "age"]), species = "human",
        se = T, lwd = 2, col = cs_species["human"])
    fun_addRegressionLine_lmer(model_dog_nocov, "age", 0, max(out_mat_all[cells_dog, "age"]), species = "dog",
        se = T, lwd = 2, col = cs_species["dog"])
    text(xlm, ylm * 0.1, pos = 2, col = cs_species["human"], cex = 0.7, label = paste(round(summary(model_human_nocov)$coefficients[2,
        1], 1), " mutation/yr/Gb; R^2=", round(cor(df[cells_human, "burden"], df[cells_human, "age"])^2,
        2), "; p=", format(summary(model_human)$coefficients[2, 5], scientific = T, digit = 2), sep = ""))
    text(xlm, ylm * 0.15, pos = 2, col = cs_species["dog"], cex = 0.7, label = paste(round(summary(model_dog_nocov)$coefficients[2,
        1], 1), " mutation/yr/Gb; R^2=", round(cor(df[cells_dog, "burden"], df[cells_dog, "age"])^2,
        2), "; p=", format(summary(model_dog)$coefficients[2, 5], scientific = T, digit = 2), sep = ""))
    text(xlm, ylm * 0.05, pos = 2, col = "black", cex = 0.7, paste("slope diff p=", format(model_mixed_coef["age:specieshuman",
        5], scientific = T, digit = 2), "; intercept diff p=", format(model_mixed_coef["specieshuman",
        5], scientific = T, digit = 2), sep = ""))
    text(xlm, ylm * 0.2, pos = 2, col = "black", cex = 0.7, label = paste("ICC(variance from donor)=",
        icc, sep = ""))
    legend("topleft", pch = c(19, 19), legend = c("dog", "human"), bty = "n", pt.cex = 1.2, col = cs_species)
}
# define young as \342\211\24420 yrs
cells_hs_young <- cells_human[meta_pd[cells_human, "age"] <= 20]
cells_hs_highcov <- cells_human[out_mat_all[cells_human, "coverage"] >= 20]
cells_hs_young_highcov <- cells_hs_young[out_mat_all[cells_hs_young, "coverage"] >= 20]
# plot
pdf("figures/highCov_mutation_burden.pdf", width = 4.5, height = 4.5, useDingbats = F)
par(mar = c(4, 4, 1, 1), tcl = 0.3, bty = "n")
fun_mutburden(out_mat_all, cells_dog_highCov, cells_hs_highcov, "sSNV_burden", 110, "sSNV/Gb")
fun_mutburden(out_mat_all, cells_dog_highCov, cells_hs_highcov, "sSNV_burden", 20, "sSNV/Gb")
fun_mutburden(out_mat_all, cells_dog_highCov, cells_hs_young_highcov, "sSNV_burden", 20, "sSNV/Gb")
fun_mutburden(out_mat_all, cells_dog_highCov, cells_hs_highcov, "sINDEL_burden", 110, "sindel/Gb")
fun_mutburden(out_mat_all, cells_dog_highCov, cells_hs_highcov, "sINDEL_burden", 20, "sindel/Gb")
fun_mutburden(out_mat_all, cells_dog_highCov, cells_hs_young_highcov, "sINDEL_burden", 20, "sindel/Gb")
dev.off()

pdf("figures/highCov_mutation_burden_signatures.pdf", width = 4.5, height = 4.5, useDingbats = F)
par(mar = c(4, 4, 1, 1), tcl = 0.3, bty = "n")
for (s in c("signature A1", "signature A2", "signature ID-A", "signature ID-B")) {
    fun_mutburden(out_mat_all, cells_dog_highCov, cells_hs_highcov, s, 110, paste(s, "/Gb", sep = ""))
    fun_mutburden(out_mat_all, cells_dog_highCov, cells_hs_highcov, s, 20, paste(s, "/Gb", sep = ""))
    fun_mutburden(out_mat_all, cells_dog_highCov, cells_hs_young_highcov, s, 20, paste(s, "/Gb", sep = ""))
}
dev.off()

# Mutation burden for the six basic SNV classes -------------------------------
fun_mutburden <- function(out_mat_all, cells_dog, cells_human, colname, xlm, ylab) {
    ylm <- max(out_mat_all[c(cells_dog, cells_human)[which(out_mat_all[c(cells_dog, cells_human), "age"] <
        xlm)], colname]/len_human_autosome)
    ylm <- ceiling(ylm/(10^floor(log10(ylm)))) * (10^floor(log10(ylm)))
    plot(out_mat_all[cells_human, "age"], out_mat_all[cells_human, colname]/len_human_autosome, pch = 21,
        col = "white", bg = cs_species["human"], lwd = 0.5, cex = 1.5, ylab = ylab, ylim = c(0, ylm),
        xlim = c(0, xlm), xlab = "age (yrs)")
    points(out_mat_all[cells_dog, "age"], out_mat_all[cells_dog, colname]/len_human_autosome, bg = cs_species["dog"],
        col = "white", pch = 21, lwd = 0.5, cex = 1.5)
    df <- out_mat_all[c(cells_dog, cells_human), c("age", colname, "species", "donor", "sex", "coverage",
        "mapd")]
    colnames(df)[2] <- "burden"
    df[, 2] <- df[, 2]/len_human_autosome
    model_dog <- lmer(burden ~ age + sex + coverage + mapd + (1 | donor), data = df[cells_dog, ])
    model_human <- lmer(burden ~ age + sex + coverage + mapd + (1 | donor), data = df[cells_human, ])
    model_mixed_coef <- summary(lmer(burden ~ age * species + sex + coverage + mapd + (1 | donor), data = df))$coefficients
    model_mixed_varcor <- as.data.frame(summary(lmer(burden ~ age * species + sex + coverage + mapd +
        (1 | donor), data = df))$varcor)
    model_dog_nocov <- lmer(burden ~ age + (1 | donor), data = df[cells_dog, ])
    model_human_nocov <- lmer(burden ~ age + (1 | donor), data = df[cells_human, ])
    icc <- round(model_mixed_varcor[1, 4]/(model_mixed_varcor[1, 4] + model_mixed_varcor[2, 4]), 2)
    fun_addRegressionLine_lmer(model_human_nocov, "age", 0, max(out_mat_all[cells_human, "age"]), species = "human",
        se = T, lwd = 2, col = cs_species["human"])
    fun_addRegressionLine_lmer(model_dog_nocov, "age", 0, max(out_mat_all[cells_dog, "age"]), species = "dog",
        se = T, lwd = 2, col = cs_species["dog"])
    text(xlm, ylm * 0.1, pos = 2, col = cs_species["human"], cex = 0.7, label = paste(round(summary(model_human_nocov)$coefficients[2,
        1], 1), " mutation/yr/Gb; R^2=", round(cor(df[cells_human, "burden"], df[cells_human, "age"])^2,
        2), "; p=", format(summary(model_human)$coefficients[2, 5], scientific = T, digit = 2), sep = ""))
    text(xlm, ylm * 0.15, pos = 2, col = cs_species["dog"], cex = 0.7, label = paste(round(summary(model_dog_nocov)$coefficients[2,
        1], 1), " mutation/yr/Gb; R^2=", round(cor(df[cells_dog, "burden"], df[cells_dog, "age"])^2,
        2), "; p=", format(summary(model_dog)$coefficients[2, 5], scientific = T, digit = 2), sep = ""))
    text(xlm, ylm * 0.05, pos = 2, col = "black", cex = 0.7, paste("slope diff p=", format(model_mixed_coef["age:specieshuman",
        5], scientific = T, digit = 2), "; intercept diff p=", format(model_mixed_coef["specieshuman",
        5], scientific = T, digit = 2), sep = ""))
    text(xlm, ylm * 0.2, pos = 2, col = "black", cex = 0.7, label = paste("ICC(variance from donor)=",
        icc, sep = ""))
    legend("topleft", pch = c(19, 19), legend = c("dog", "human"), bty = "n", pt.cex = 1.2, col = cs_species)
}
# define young as \342\211\24420 yrs
cells_hs_young <- cells_human[meta_pd[cells_human, "age"] <= 20]
# get 6 basic types of SNV number
tdf <- cbind(snv.sig, snv.sig_pd)
tps <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
for (t in tps) {
    rns <- row.names(tdf)[grep(t, row.names(tdf))]
    cls <- row.names(out_mat_all)
    out_mat_all[, t] <- apply(tdf[rns, cls], 2, sum)/apply(tdf[, cls], 2, sum) * out_mat_all[cls, ]$sSNV_burden
}
pdf("figures/mutation_burden_6types_SNVs.pdf", width = 4.5, height = 4.5, useDingbats = F)
for (t in tps) {
    par(mar = c(4, 4, 1, 1), tcl = 0.3, bty = "n")
    fun_mutburden(out_mat_all, cells_dog, cells_hs, t, 110, paste(t, "/Gb", sep = ""))
    fun_mutburden(out_mat_all, cells_dog, cells_hs_young, t, 20, paste(t, "/Gb", sep = ""))
}
dev.off()

# mutation burden with crypt----
fun_mutburden <- function(mat, cd, ch, colname, xlm, ylab) {
    ylm <- max(mat[c(cd, ch)[which(mat[c(cd, ch), "age"] < xlm)], colname], na.rm = T)
    ylm <- ceiling(ylm/(10^floor(log10(ylm)))) * (10^floor(log10(ylm)))
    plot(mat[ch, "age"], mat[ch, colname], pch = 21, col = "white", bg = cs_species["human"], lwd = 0.5,
        cex = 1.5, ylab = ylab, ylim = c(0, ylm), xlim = c(0, xlm), xlab = "age (yrs)")
    points(mat[cd, "age"], mat[cd, colname], bg = cs_species["dog"], col = "white", pch = 21, lwd = 0.5,
        cex = 1.5)
    df <- mat[c(cd, ch), c("age", colname, "species", "donor")]
    colnames(df)[2] <- "burden"
    df[, 2] <- df[, 2] + 0.1
    model_dog <- lmer(burden ~ age + (1 | donor), data = df[cd, ])
    model_human <- lmer(burden ~ age + (1 | donor), data = df[ch, ])
    model_mixed_coef <- summary(lmer(burden ~ age * species + (1 | donor), data = df))$coefficients
    model_mixed_varcor <- as.data.frame(summary(lmer(burden ~ age * species + (1 | donor), data = df))$varcor)
    icc <- round(model_mixed_varcor[1, 4]/(model_mixed_varcor[1, 4] + model_mixed_varcor[2, 4]), 2)
    fun_addRegressionLine_lmer(model_human, "age", 0, max(mat[ch, "age"]), species = "human", se = T,
        lwd = 2, col = cs_species["human"])
    fun_addRegressionLine_lmer(model_dog, "age", 0, max(mat[cd, "age"]), species = "dog", se = T, lwd = 2,
        col = cs_species["dog"])
    text(xlm, ylm * 0.1, pos = 2, col = cs_species["human"], cex = 0.7, label = paste(round(summary(model_human)$coefficients[2,
        1], 1), " mutation/yr/Gb; R^2=", round(cor(df[ch, "burden"], df[ch, "age"])^2, 2), "; p=", format(summary(model_human)$coefficients[2,
        5], scientific = T, digit = 2), sep = ""))
    text(xlm, ylm * 0.15, pos = 2, col = cs_species["dog"], cex = 0.7, label = paste(round(summary(model_dog)$coefficients[2,
        1], 1), " mutation/yr/Gb; R^2=", round(cor(df[cd, "burden"], df[cd, "age"])^2, 2), "; p=", format(summary(model_dog)$coefficients[2,
        5], scientific = T, digit = 2), sep = ""))
    text(xlm, ylm * 0.05, pos = 2, col = "black", cex = 0.7, paste("slope diff p=", format(model_mixed_coef["age:specieshuman",
        5], scientific = T, digit = 2), "; intercept diff p=", format(model_mixed_coef["specieshuman",
        5], scientific = T, digit = 2), sep = ""))
    text(xlm, ylm * 0.2, pos = 2, col = "black", cex = 0.7, label = paste("ICC(variance from donor)=",
        icc, sep = ""))
    legend("topleft", pch = c(19, 19), legend = c("dog", "human"), bty = "n", pt.cex = 1.2, col = cs_species)
}

meta_crypt[, "SNV_per_Gb"] <- meta_crypt[, "Substitution.burden"]/meta_crypt[, "Total.genome.size..bp."]/2 *
    1e+09
meta_crypt[, "indel_per_Gb"] <- meta_crypt[, "Indel.burden"]/meta_crypt[, "Total.genome.size..bp."]/2 *
    1e+09
meta_crypt[, "age"] <- meta_crypt[, "Age..years."]
meta_crypt[, "species"] <- meta_crypt[, "Species"]
meta_crypt[, "donor"] <- meta_crypt[, "Individual"]

tmp_out <- out_mat_all[, c("donor", "age", "species", "sSNV_burden", "sINDEL_burden")]
tmp_out[, 4:5] <- tmp_out[, 4:5]/len_human_autosome
colnames(tmp_out)[4:5] <- c("SNV_per_Gb", "indel_per_Gb")
tmp_out <- rbind(tmp_out, meta_crypt[, c("donor", "age", "species", "SNV_per_Gb", "indel_per_Gb")])

nmf_res <- readRDS("tables/dnaseq/scan2/signatures.snv_with_crypt.rds")
t <- apply(nmf_res$contribution * colSums(nmf_res$signatures), 2, function(x) {
    return(x/sum(x))
})
n_sigs <- t(t)
t_cells <- c(cells_dog, cells_hs, cells_crypt)
n_sigs[t_cells, ] <- n_sigs[t_cells, ] * tmp_out[t_cells, ]$SNV_per_Gb
tmp_out <- cbind(tmp_out, n_sigs[row.names(tmp_out), ])

nmf_res <- readRDS("tables/dnaseq/scan2/signatures.indel_with_crypt.rds")
t <- apply(nmf_res$contribution * colSums(nmf_res$signatures), 2, function(x) {
    return(x/sum(x))
})
n_sigs <- t(t)
n_sigs[t_cells, ] <- n_sigs[t_cells, ] * tmp_out[t_cells, ]$indel_per_Gb
tmp_out <- cbind(tmp_out, n_sigs[row.names(tmp_out), ])
tmp_out[tmp_out[, "species"] == "Human", "species"] <- "human"
tmp_out[tmp_out[, "species"] == "Dog", "species"] <- "dog"
tmp_out[, "celltype"] <- "neuron"
tmp_out[cells_crypt, "celltype"] <- "crypt"
write.table(tmp_out, "tables/For_manuscript/meta_combined_with_crypt.txt", sep = "\t", quote = F)

pdf("figures/mutation_burden.with_crypt.pdf", width = 4.5, height = 4.5, useDingbats = F)
par(mar = c(4, 4, 1, 1), tcl = 0.3, bty = "n")
fun_mutburden(tmp_out, cells_crypt_dog, cells_crypt_human, "SNV_per_Gb", 110, "sSNV/Gb; crypt")
fun_mutburden(tmp_out, cells_crypt_dog, cells_crypt_human, "indel_per_Gb", 110, "sindel/Gb; crypt")
for (i in colnames(tmp_out)[6:11]) {
    fun_mutburden(tmp_out, cells_crypt_dog, cells_crypt_human, i, 110, paste(i, "/Gb; crypt", sep = ""))
}
for (i in colnames(tmp_out)[6:11]) {
    fun_mutburden(tmp_out, cells_dog, cells_hs, i, 110, paste(i, "/Gb; neuron", sep = ""))
}
dev.off()

# rate of mutation in CpG----
cToTs <- row.names(snv.sig)[33:48]

pdf("figures/CtoT_in_CpG/dog.pdf", width = 3.5, height = 2.2 * 3, useDingbats = F)
par(mfrow = c(3, 1), mar = c(5, 5, 1, 1), tcl = 0.3, bty = "n", las = 3)
tcs <- rep("#D97A74", 16)
tcs[c(3, 7, 11, 15)] <- "#8F2B25"
x <- apply(snv.sig[cToTs, cells_dog], 1, sum)
barplot(x/sum(x) * 100, col = tcs, space = 0, border = "white", ylab = "C>T mutations\nin tri-nucleotide context (%)",
    xaxt = "n", ylim = c(0, 15))
legend("topright", pch = 15, pt.cex = 2, bty = "n", col = tcs[3], legend = paste(round(sum((x/sum(x) *
    100)[c(3, 7, 11, 15)]), 1), "% in CpG", sep = ""))
y <- tri.context.dog[cToTs, 1]
barplot(y/sum(y) * 100, col = tcs, space = 0, border = "white", ylab = "genomic tri-nucleotide context (%)",
    xaxt = "n", ylim = c(0, 15))
legend("topright", pch = 15, pt.cex = 2, bty = "n", col = tcs[3], legend = paste(round(sum((y/sum(y) *
    100)[c(3, 7, 11, 15)]), 1), "% CpG", sep = ""))
barplot(x/sum(x)/(y/sum(y)), col = tcs, space = 0, border = "white", ylim = c(0, 8), ylab = "enrichment of C>T mutations\nin tri-nucleotide context",
    names = cToTs)
legend("topright", pch = 15, pt.cex = 2, bty = "n", col = tcs[3], legend = paste(round(mean((x/sum(x)/(y/sum(y)))[c(3,
    7, 11, 15)]), 1), " average fold in CpG", sep = ""))
abline(h = 1, lty = 2, col = "grey")
dev.off()

pdf("figures/CtoT_in_CpG/human.pdf", width = 3.5, height = 2.2 * 3, useDingbats = F)
par(mfrow = c(3, 1), mar = c(5, 5, 1, 1), tcl = 0.3, bty = "n", las = 3)
tcs <- rep("#D97A74", 16)
tcs[c(3, 7, 11, 15)] <- "#8F2B25"
x <- apply(snv.sig_pd[cToTs, cells_hs], 1, sum)
barplot(x/sum(x) * 100, col = tcs, space = 0, border = "white", ylab = "C>T mutations\nin tri-nucleotide context (%)",
    xaxt = "n", ylim = c(0, 15))
legend("topright", pch = 15, pt.cex = 2, bty = "n", col = tcs[3], legend = paste(round(sum((x/sum(x) *
    100)[c(3, 7, 11, 15)]), 1), "% in CpG", sep = ""))
y <- tri.context.human[cToTs, 1]
barplot(y/sum(y) * 100, col = tcs, space = 0, border = "white", ylab = "genomic tri-nucleotide context (%)",
    xaxt = "n", ylim = c(0, 15))
legend("topright", pch = 15, pt.cex = 2, bty = "n", col = tcs[3], legend = paste(round(sum((y/sum(y) *
    100)[c(3, 7, 11, 15)]), 1), "% CpG", sep = ""))
barplot(x/sum(x)/(y/sum(y)), col = tcs, space = 0, border = "white", ylim = c(0, 8), ylab = "enrichment of C>T mutations\nin tri-nucleotide context",
    names = cToTs)
legend("topright", pch = 15, pt.cex = 2, bty = "n", col = tcs[3], legend = paste(round(mean((x/sum(x)/(y/sum(y)))[c(3,
    7, 11, 15)]), 1), " average fold in CpG", sep = ""))
abline(h = 1, lty = 2, col = "grey")
dev.off()

# exposure of signatures in neurons (barplots)----
pdf("figures/signature_exposure_barplot.pdf", width = 8, height = 4, useDingbats = F)
par(mar = c(8, 4, 3, 1), tcl = 0.3, bty = "n")
tmp_cells <- cells_dog[order(out_mat_all[cells_dog, "age"])]
tdf <- as.matrix(out_mat_all[tmp_cells, c("signature A1", "signature A2")]/len_human_autosome)
barplot(t(tdf), col = c("#66c2a5", "#fc8d62"), las = 3, main = "sSNV signature exposure for dog neurons",
    ylab = "mutation/Gb", names = paste(tmp_cells, "; ", out_mat_all[tmp_cells, "age"], "yrs", sep = ""),
    cex.names = 0.6)
legend("topleft", pch = 15, col = c("#66c2a5", "#fc8d62"), legend = c("A1", "A2"), pt.cex = 2, bty = "n")
tdf <- as.matrix(out_mat_all[tmp_cells, c("signature ID-A", "signature ID-B")]/len_human_autosome)
barplot(t(tdf), col = c("#8da0cb", "#e78ac3"), las = 3, main = "sSNV signature exposure for dog neurons",
    ylab = "mutation/Gb", names = paste(tmp_cells, "; ", out_mat_all[tmp_cells, "age"], "yrs", sep = ""),
    cex.names = 0.6)
legend("topleft", pch = 15, col = c("#8da0cb", "#e78ac3"), legend = c("ID-A", "ID-B"), pt.cex = 2, bty = "n")

tmp_cells <- cells_human[order(out_mat_all[cells_human, "age"])]
tdf <- as.matrix(out_mat_all[tmp_cells, c("signature A1", "signature A2")]/len_human_autosome)
barplot(t(tdf), col = c("#66c2a5", "#fc8d62"), las = 3, main = "sSNV signature exposure for human neurons",
    ylab = "mutation/Gb", names = paste(tmp_cells, "; ", out_mat_all[tmp_cells, "age"], "yrs", sep = ""),
    cex.names = 0.2)
legend("topleft", pch = 15, col = c("#66c2a5", "#fc8d62"), legend = c("A1", "A2"), pt.cex = 2, bty = "n")
tdf <- as.matrix(out_mat_all[tmp_cells, c("signature ID-A", "signature ID-B")]/len_human_autosome)
barplot(t(tdf), col = c("#8da0cb", "#e78ac3"), las = 3, main = "sSNV signature exposure for human neurons",
    ylab = "mutation/Gb", names = paste(tmp_cells, "; ", out_mat_all[tmp_cells, "age"], "yrs", sep = ""),
    cex.names = 0.2)
legend("topleft", pch = 15, col = c("#8da0cb", "#e78ac3"), legend = c("ID-A", "ID-B"), pt.cex = 2, bty = "n")
dev.off()

# signature analysis using COSMIC----
require(NMF)
require(MutationalPatterns)

pdf("figures/cosmic_signatures/cosine_correlation.pdf", width = 8, height = 3, useDingbats = F)
# SNV
mut_mat <- cbind(snv.sig[, cells_used], snv.sig_pd[, cells_normal])
mut_mat <- apply(mut_mat, 2, function(x) {
    return(x/sum(x) * 100 + 0.001)
})
cosmic_snv <- get_known_signatures("snv")
cos_sim_snv <- cos_sim_matrix(nmf_res_snv$signatures, cosmic_snv)
pheatmap(cos_sim_snv, labels_col = sub("^SBS", "", colnames(cos_sim_snv)), breaks = c(0:100)/100)
# get top 3 similar signatures to A1 and A2 respectively
active_cosmic_snv <- unique(c(colnames(cos_sim_snv)[order(cos_sim_snv[1, ], decreasing = T)[1:3]], colnames(cos_sim_snv)[order(cos_sim_snv[2,
    ], decreasing = T)[1:3]]))
refit_snv <- fit_to_signatures(mut_mat, cosmic_snv[, active_cosmic_snv])
t <- apply(refit_snv$contribution, 2, function(x) {
    return(x/sum(x))
})
n_sigs_snv <- t(t)
n_sigs_snv <- n_sigs_snv * out_mat_all[row.names(n_sigs_snv), ]$sSNV_burden/len_human_autosome
# indel
mut_mat <- cbind(indel.sig[, cells_used], indel.sig_pd[, cells_normal])
mut_mat <- apply(mut_mat, 2, function(x) {
    return(x/sum(x) * 100 + 0.001)
})
mut_mat[is.na(mut_mat)] <- 0
cosmic_indel <- read.table("../tables/annotation/COSMIC_v3.4_ID_GRCh37.txt", header = T, row.names = 1,
    check.names = T)
cosmic_indel <- as.matrix(cosmic_indel)
cos_sim_indel <- cos_sim_matrix(nmf_res_indel$signatures, cosmic_indel)
pheatmap(cos_sim_indel, labels_col = sub("^SBS", "", colnames(cos_sim_indel)), breaks = c(0:100)/100)
dev.off()
# get top 3 similar signatures to A1 and A2 respectively
active_cosmic_indel <- unique(c(colnames(cos_sim_indel)[order(cos_sim_indel[1, ], decreasing = T)[1:3]],
    colnames(cos_sim_indel)[order(cos_sim_indel[2, ], decreasing = T)[1:3]]))
refit_indel <- fit_to_signatures(mut_mat, cosmic_indel[, active_cosmic_indel])
t <- apply(refit_indel$contribution, 2, function(x) {
    return(x/sum(x))
})
n_sigs_indel <- t(t)
n_sigs_indel <- n_sigs_indel * out_mat_all[row.names(n_sigs_indel), ]$sINDEL_burden/len_human_autosome

# plot
pdf("figures/cosmic_signatures/signature_exposure_barplot.pdf", width = 8, height = 4, useDingbats = F)
par(mar = c(8, 4, 3, 1), tcl = 0.3, bty = "n")
tmp_cells <- cells_dog[order(out_mat_all[cells_dog, "age"])]
barplot(t(n_sigs_snv[tmp_cells, order(apply(n_sigs_snv, 2, sum), decreasing = T)]), col = brewer.pal(length(active_cosmic_snv),
    "Set3"), ylab = "mutations/Gb", main = "sSNV in dog neurons", names = paste(tmp_cells, "; ", out_mat_all[tmp_cells,
    "age"], "yrs", sep = ""), cex.names = 0.6, las = 3)
legend("topleft", legend = active_cosmic_snv[order(apply(n_sigs_snv, 2, sum), decreasing = T)], col = brewer.pal(length(active_cosmic_snv),
    "Set3"), pch = 15, pt.cex = 2, bty = "n")
barplot(t(n_sigs_indel[tmp_cells, order(apply(n_sigs_indel, 2, sum), decreasing = T)]), col = brewer.pal(length(active_cosmic_indel),
    "Set3"), ylab = "mutations/Gb", main = "sindels in dog neurons", names = paste(tmp_cells, "; ", out_mat_all[tmp_cells,
    "age"], "yrs", sep = ""), cex.names = 0.6, las = 3)
legend("topleft", legend = active_cosmic_indel[order(apply(n_sigs_indel, 2, sum), decreasing = T)], col = brewer.pal(length(active_cosmic_indel),
    "Set3"), pch = 15, pt.cex = 2, bty = "n")
tmp_cells <- cells_human[order(out_mat_all[cells_human, "age"])]
barplot(t(n_sigs_snv[tmp_cells, order(apply(n_sigs_snv, 2, sum), decreasing = T)]), col = brewer.pal(length(active_cosmic_snv),
    "Set3"), ylab = "mutations/Gb", main = "sSNV in human neurons", names = paste(tmp_cells, "; ", out_mat_all[tmp_cells,
    "age"], "yrs", sep = ""), cex.names = 0.2, las = 3)
legend("topleft", legend = active_cosmic_snv[order(apply(n_sigs_snv, 2, sum), decreasing = T)], col = brewer.pal(length(active_cosmic_snv),
    "Set3"), pch = 15, pt.cex = 2, bty = "n")
barplot(t(n_sigs_indel[tmp_cells, order(apply(n_sigs_indel, 2, sum), decreasing = T)]), col = brewer.pal(length(active_cosmic_indel),
    "Set3"), ylab = "mutations/Gb", main = "sindels in human neurons", names = paste(tmp_cells, "; ",
    out_mat_all[tmp_cells, "age"], "yrs", sep = ""), cex.names = 0.2, las = 3)
legend("topleft", legend = active_cosmic_indel[order(apply(n_sigs_indel, 2, sum), decreasing = T)], col = brewer.pal(length(active_cosmic_indel),
    "Set3"), pch = 15, pt.cex = 2, bty = "n")
dev.off()

# Refit all samples to the complete COSMIC signature collections.
cosmic_snv <- get_known_signatures("snv")
refit_snv <- fit_to_signatures(mut_mat, cosmic_snv)

t <- apply(refit_snv$contribution * colSums(cosmic_snv), 2, function(x) {
    return(x/sum(x))
})
n_sigs <- t(t)
n_sigs[cells_used, ] <- n_sigs[cells_used, ] * out_mat[cells_used, ]$sSNV_burden

cosmic_indel <- read.table("tables/annotation/COSMIC_v3.4_ID_GRCh37.txt", header = T, row.names = 1,
    check.names = T)
cosmic_indel <- as.matrix(cosmic_indel)
refit_indel <- fit_to_signatures(mut_mat_indel_hs, cosmic_indel)
contri <- refit_indel$contribution * colSums(cosmic_indel)
tn <- cells_Jenn[order(meta_pd[cells_Jenn, "age"])]
tn <- tn[order(meta_pd[tn, "disease"])]
pdf("figures/COSMIC_contribution_INDEL.pdf", width = 10, height = 5, useDingbats = F)
par(mar = c(4, 4, 1, 1), tcl = 0.3, bty = "n", las = 3, cex.axis = 0.3)
barplot(contri[, tn], col = cs25[1:18], border = F, space = 0, names = meta_pd[tn, "age"], ylab = "n. of detected signatures")
n_normal <- sum(meta_pd[cells_Jenn, "disease"] == "normal")
n_pd <- length(cells_Jenn)[1] - n_normal
abline(v = n_normal)
text(c(n_normal/2, n_normal + n_pd/2), rep(max(contri), 2) * 0.9, label = c("normal", "PD"), cex = 1)
legend("topleft", legend = colnames(cosmic_indel), col = cs25[1:18], pch = 15, ncol = 2)
dev.off()
pdf("figures/COSMIC_exposure_INDEL.pdf", width = 5, height = 4, useDingbats = F)
par(mar = c(4, 4, 3, 1), tcl = 0.3, bty = "n", las = 3, cex.axis = 0.7)
barplot(apply(contri[, cells_Jenn[meta_pd[cells_Jenn, "disease"] == "normal"]], 1, sum)/sum(contri[,
    cells_Jenn[meta_pd[cells_Jenn, "disease"] == "normal"]]), col = cs25[1:18], ylab = "exposure", main = "normal neurons",
    ylim = c(0, 0.4))
abline(h = 0.05)
barplot(apply(contri[, cells_Jenn[meta_pd[cells_Jenn, "disease"] == "PD"]], 1, sum)/sum(contri[, cells_Jenn[meta_pd[cells_Jenn,
    "disease"] == "PD"]]), col = cs25[1:18], ylab = "exposure", main = "PD neurons", ylim = c(0, 0.4))
abline(h = 0.05)
p <- c()
for (i in row.names(contri)) {
    p <- c(p, t.test(contri[i, cells_pd], contri[i, cells_normal])$p.value)
}
p <- p.adjust(p, method = "fdr")
dev.off()

# Power analysis of dog-human mutation-accumulation differences ---------------
set.seed(1)

# Split once
df_h <- subset(df, species == "human")
df_d <- subset(df, species == "dog")

# Fit your current LM to get baseline coefficients and residual SD for reference
m0 <- lm(burden ~ age * species, data = df)
b_hat <- coef(m0)  # '(Intercept)', 'age', 'specieshuman', 'age:specieshuman'
sigma_hat <- sigma(m0)  # overall residual SD (not used below if you fix human y)

# Helper: one simulation + test
sim_once_lm_fix_human <- function(df_h, df_d_sample, delta, sigma_dog) {
    # Keep human y fixed; only simulate dog y
    # Build model matrix for predictions
    Xh <- model.matrix(~age * species, data = df_h)
    Xd <- model.matrix(~age * species, data = df_d_sample)

    # Set the *true* interaction to desired delta; keep other coefs as estimated
    s_h <- b_hat["age"] + b_hat["age:specieshuman"]  # current human slope (keep fixed)
    # Target: human - dog = delta  => dog slope = s_h - delta
    b_true <- b_hat
    b_true["age"] <- s_h - delta  # set dog slope
    b_true["age:specieshuman"] <- s_h - b_true["age"]  # keep human slope = s_h  (this equals delta)


    # Means
    mu_h <- as.numeric(Xh %*% b_true)
    mu_d <- as.numeric(Xd %*% b_true)

    # Generate dog outcomes only; keep human as observed
    y_h <- df_h$burden
    y_d <- mu_d + rnorm(nrow(df_d_sample), 0, sigma_dog)

    dat <- rbind(transform(df_h, burden_sim = y_h), transform(df_d_sample, burden_sim = y_d))

    fit <- lm(burden_sim ~ age * species, data = dat)
    p <- summary(fit)$coef["age:specieshuman", "Pr(>|t|)"]
    as.numeric(!is.na(p) && p < 0.05)
}

#plot(dat$age,dat$burden_sim,col=cs_species[dat$species],pch=19)

# Power curve vs. dog donors (bootstrap donors; keeps your dog age distribution)
power_vs_dog_donors_lm <- function(n_dog_vec, delta, sigma_dog, nsim = 1000) {
    # collapse dogs to donor means if you want donors; otherwise use cells
    # Here: donors (recommended). If you prefer cells, use the raw df_d instead.
    dog_donor <- aggregate(burden ~ donor + species + age, data = df_d, mean)

    sapply(n_dog_vec, function(nD) {
        mean(replicate(nsim, {
            idx <- sample(seq_len(nrow(dog_donor)), size = nD, replace = TRUE)
            df_d_sample <- dog_donor[idx, , drop = FALSE]
            sim_once_lm_fix_human(df_h = aggregate(burden ~ donor + species + age, data = df_h, mean),
                df_d_sample = df_d_sample, delta = delta, sigma_dog = sigma_dog)
        }))
    })
}

# Example usage:
pdf("figures/power_analysis.pdf", width = 5, height = 5, useDingbats = F)
dogs <- seq(5, 80, by = 5)  # number of dog donors to simulate
delta <- (-2)  # target slope difference (human - dog) in mut/yr
sigma_dog <- 17.34  # dog residual SD at donor level (tune to your data)
pow <- power_vs_dog_donors_lm(dogs, delta, sigma_dog, nsim = 100)
plot(dogs, pow, type = "b", xlab = "Dog donors (humans fixed, donor means)", ylab = sprintf("Power to detect \316\224 = %.2f (lm)",
    delta), ylim = c(0, 1), xlim = c(0, 80))
abline(h = 0.8, lty = 2)
abline(v = 34, lty = 2)

delta <- (-1)  # target slope difference (human - dog) in mut/yr
sigma_dog <- 17.34  # dog residual SD at donor level (tune to your data)
pow <- power_vs_dog_donors_lm(dogs, delta, sigma_dog, nsim = 100)
plot(dogs, pow, type = "b", xlab = "Dog donors (humans fixed, donor means)", ylab = sprintf("Power to detect \316\224 = %.2f (lm)",
    delta), ylim = c(0, 1), xlim = c(0, 80))
abline(h = 0.8, lty = 2)
abline(v = 34, lty = 2)
dev.off()

# trinucleotide context for dog and human----
tcs <- c()
for (i in 1:6) {
    tcs <- c(tcs, rep(cs_gg[i], 16))
}
pdf("figures/trinucleotide.background.pdf", width = 12, height = 7.5, useDingbats = F)
par(mfrow = c(2, 1), mar = c(6, 4, 3, 1), tcl = 0.3, bty = "n", cex = 5/6)
barplot(tri.context.dog[, 1], col = tcs, border = "white", space = 0, names = row.names(snv.sig), las = 3,
    main = "dog (canFam3)", ylab = "n. of genomic background")
abline(v = c(1:5) * 16)
text(c(1:6) * 16 - 8, max(tri.context.dog[, 1]) * 19/20, label = c("C>A", "C>G", "C>T", "T>A", "T>C",
    "T>G"))
barplot(tri.context.human[, 1], col = tcs, border = "white", space = 0, names = row.names(snv.sig), las = 3,
    main = "human (hs37d5)", ylab = "n. of genomic background")
abline(v = c(1:5) * 16)
text(c(1:6) * 16 - 8, max(tri.context.human[, 1]) * 19/20, label = c("C>A", "C>G", "C>T", "T>A", "T>C",
    "T>G"))
dev.off()

# make normalization factors to correct dog trinucleotide context to the ones of human
tri.factors <- tri.context.dog[, 1]/sum(tri.context.dog[, 1])/(tri.context.human[, 1]/sum(tri.context.human[,
    1]))
names(tri.factors) <- row.names(tri.context.dog)
write.table(tri.factors, "tables/annotation/trinucleotide.factors.txt", quote = F, sep = "\t", col.names = F)

# overall mutation profiles----
# dog SNV----
rns <- cells_used[order(out_mat[cells_used, "age"])]
smps <- samples_meta[rns, "donor"]
tcs <- col_sig
pdf("figures/snv.profiles.samples.pdf", width = 12, height = 7.5/2, useDingbats = F)
par(tcl = 0.3, bty = "n", cex = 5/6)
for (cn in rns) {
    barplot(snv.sig[, cn], col = tcs, border = "white", space = 0, names = row.names(snv.sig), las = 3,
        main = cn, ylab = "N of identified sSNV")
    abline(v = c(1:5) * 16)
    text(c(1:6) * 16 - 8, max(snv.sig[, cn]) * 19/20, label = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"))
}
dev.off()
pdf("figures/snv.profiles.donors.pdf", width = 12, height = 7.5/2, useDingbats = F)
par(tcl = 0.3, bty = "n", cex = 5/6)
barplot(apply(snv.sig[, rns], 1, sum), col = tcs, border = "white", space = 0, ylim = c(0, 140), names = row.names(snv.sig),
    las = 3, main = "all donors", ylab = "N of identified sSNV")
abline(v = c(1:5) * 16)
text(c(1:6) * 16 - 8, max(apply(snv.sig[, rns], 1, sum)) * 19/20, label = c("C>A", "C>G", "C>T", "T>A",
    "T>C", "T>G"))
for (s in unique(smps)) {
    cns <- which(smps == s)
    if (length(cns) > 1) {
        x <- apply(snv.sig[, rns[cns]], 1, sum)
    } else {
        x <- snv.sig[, rns[cns]]
    }
    barplot(x, col = tcs, border = "white", space = 0, names = row.names(snv.sig), las = 3, main = paste(s,
        "; ", samples[s, 1], "; ", samples[s, 4], "yrs", sep = ""), ylab = "N of identified sSNV")
    abline(v = c(1:5) * 16)
    text(c(1:6) * 16 - 8, max(x) * 19/20, label = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"))
}
dev.off()

# dog indel----
tcs <- c(rep(rgb(244, 192, 122, maxColorValue = 255), 6), rep(rgb(238, 133, 51, maxColorValue = 255),
    6), rep(rgb(182, 217, 145, maxColorValue = 255), 6), rep(rgb(87, 157, 64, maxColorValue = 255), 6),
    rep(rgb(243, 203, 183, maxColorValue = 255), 6), rep(rgb(237, 142, 112, maxColorValue = 255), 6),
    rep(rgb(221, 83, 64, maxColorValue = 255), 6), rep(rgb(172, 45, 38, maxColorValue = 255), 6), rep(rgb(211,
        223, 239, maxColorValue = 255), 6), rep(rgb(156, 193, 222, maxColorValue = 255), 6), rep(rgb(95,
        150, 197, maxColorValue = 255), 6), rep(rgb(46, 98, 165, maxColorValue = 255), 6), rep(rgb(225,
        223, 239, maxColorValue = 255), 1), rep(rgb(181, 184, 213, maxColorValue = 255), 2), rep(rgb(134,
        132, 187, maxColorValue = 255), 3), rep(rgb(91, 67, 148, maxColorValue = 255), 5))
ns <- c(rep(c(1, 2, 3, 4, 5, "6+"), 2), rep(c(0, 1, 2, 3, 4, "5+"), 2), rep(c(1, 2, 3, 4, 5, "6+"), 4),
    rep(c(0, 1, 2, 3, 4, "5+"), 4), c(1, 1, 2, 1, 2, 3, 1, 2, 3, 4, "5+"))
pdf("figures/indel.profiles.samples.pdf", width = 12, height = 7.5/2, useDingbats = F)
par(tcl = 0.3, bty = "n", cex = 5/6)
for (cn in rns) {
    barplot(indel.sig[, cn], col = tcs, border = "white", space = 0, las = 1, names = ns, las = 3, main = cn,
        ylab = "N of identified sindel")
    axis(1, c(6.5, 18.5, 36.5, 60.5, 77.5), label = c("homopolymer\nlength", "homopolymer\nlength", "number of\nrepeat units",
        "number of\nrepeat units", "microhomology\nlength"), lwd = 0, padj = 2, cex.axis = 0.5)
    axis(3, c(6.5, 18.5, 36.5, 60.5, 77.5), label = c("1bp deletion", "1bp insertion", ">1 bp deletion at repeats",
        ">1 bp insertion at repeats", "microhomology"), lwd = 0, padj = 2, cex.axis = 0.5)
    text(c(3.5, 9.5, 15.5, 21.5, 27.5, 33.5, 39.5, 45.5, 51.5, 57.5, 63.5, 69.5, 73, 74.5, 77, 81), max(indel.sig[,
        cn]) * 0.95, label = c(rep(c("C", "T"), 2), rep(c(2, 3, 4, "5+"), 2), 2, 3, 4, "5+"), cex = 0.5,
        col = unique(tcs))
}
dev.off()

pdf("figures/indel.profiles.donors.pdf", width = 12, height = 7.5/2, useDingbats = F)
par(tcl = 0.3, bty = "n", cex = 5/6)
barplot(apply(indel.sig[, rns], 1, sum), col = tcs, border = "white", space = 0, las = 1, names = ns,
    las = 3, main = "all donors", ylab = "N of identified sindel")
axis(1, c(6.5, 18.5, 36.5, 60.5, 77.5), label = c("homopolymer\nlength", "homopolymer\nlength", "number of\nrepeat units",
    "number of\nrepeat units", "microhomology\nlength"), lwd = 0, padj = 2, cex.axis = 0.5)
axis(3, c(6.5, 18.5, 36.5, 60.5, 77.5), label = c("1bp deletion", "1bp insertion", ">1 bp deletion at repeats",
    ">1 bp insertion at repeats", "microhomology"), lwd = 0, padj = 2, cex.axis = 0.5)
text(c(3.5, 9.5, 15.5, 21.5, 27.5, 33.5, 39.5, 45.5, 51.5, 57.5, 63.5, 69.5, 73, 74.5, 77, 81), max(apply(indel.sig[,
    rns], 1, sum)) * 0.95, label = c(rep(c("C", "T"), 2), rep(c(2, 3, 4, "5+"), 2), 2, 3, 4, "5+"), cex = 0.5,
    col = unique(tcs))
for (s in unique(smps)) {
    cns <- which(smps == s)
    if (length(cns) > 1) {
        x <- apply(indel.sig[, rns[cns]], 1, sum)
    } else {
        x <- indel.sig[, rns[cns]]
    }
    barplot(x, col = tcs, border = "white", space = 0, las = 1, names = ns, las = 3, main = paste(s,
        "; ", samples[s, 1], "; ", samples[s, 4], "yrs", sep = ""), ylab = "N of identified sindel")
    axis(1, c(6.5, 18.5, 36.5, 60.5, 77.5), label = c("homopolymer\nlength", "homopolymer\nlength", "number of\nrepeat units",
        "number of\nrepeat units", "microhomology\nlength"), lwd = 0, padj = 2, cex.axis = 0.5)
    axis(3, c(6.5, 18.5, 36.5, 60.5, 77.5), label = c("1bp deletion", "1bp insertion", ">1 bp deletion at repeats",
        ">1 bp insertion at repeats", "microhomology"), lwd = 0, padj = 2, cex.axis = 0.5)
    text(c(3.5, 9.5, 15.5, 21.5, 27.5, 33.5, 39.5, 45.5, 51.5, 57.5, 63.5, 69.5, 73, 74.5, 77, 81), max(x) *
        0.95, label = c(rep(c("C", "T"), 2), rep(c(2, 3, 4, "5+"), 2), 2, 3, 4, "5+"), cex = 0.5, col = unique(tcs))
}
dev.off()

# human SNV----
rns <- cells_normal[order(meta_pd[cells_normal, "age"])]
smps <- meta_pd[rns, "donor"]
tcs <- col_sig
pdf("figures/snv.profiles.samples.human.pdf", width = 12, height = 7.5/2, useDingbats = F)
par(tcl = 0.3, bty = "n", cex = 5/6)
for (cn in rns) {
    barplot(snv.sig_pd[, cn], col = tcs, border = "white", space = 0, names = row.names(snv.sig), las = 3,
        main = cn, ylab = "N of identified sSNV")
    abline(v = c(1:5) * 16)
    text(c(1:6) * 16 - 8, max(snv.sig_pd[, cn]) * 19/20, label = c("C>A", "C>G", "C>T", "T>A", "T>C",
        "T>G"))
}
dev.off()
pdf("figures/snv.profiles.donors.human.pdf", width = 12, height = 7.5/2, useDingbats = F)
par(tcl = 0.3, bty = "n", cex = 5/6)
barplot(apply(snv.sig_pd[, rns], 1, sum), col = tcs, border = "white", space = 0, names = row.names(snv.sig_pd),
    las = 3, main = "all donors", ylab = "N of identified sSNV")
abline(v = c(1:5) * 16)
text(c(1:6) * 16 - 8, max(apply(snv.sig_pd[, rns], 1, sum)) * 19/20, label = c("C>A", "C>G", "C>T", "T>A",
    "T>C", "T>G"))
for (s in unique(smps)) {
    cns <- which(smps == s)
    if (length(cns) > 1) {
        x <- apply(snv.sig_pd[, rns[cns]], 1, sum)
    } else {
        x <- snv.sig_pd[, rns[cns]]
    }
    barplot(x, col = tcs, border = "white", space = 0, names = row.names(snv.sig_pd), las = 3, main = paste(s,
        "; ", meta_pd[rns[cns[1]], "age"], "yrs", sep = ""), ylab = "N of identified sSNV")
    abline(v = c(1:5) * 16)
    text(c(1:6) * 16 - 8, max(x) * 19/20, label = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"))
}
dev.off()

# human indel----
tcs <- c(rep(rgb(244, 192, 122, maxColorValue = 255), 6), rep(rgb(238, 133, 51, maxColorValue = 255),
    6), rep(rgb(182, 217, 145, maxColorValue = 255), 6), rep(rgb(87, 157, 64, maxColorValue = 255), 6),
    rep(rgb(243, 203, 183, maxColorValue = 255), 6), rep(rgb(237, 142, 112, maxColorValue = 255), 6),
    rep(rgb(221, 83, 64, maxColorValue = 255), 6), rep(rgb(172, 45, 38, maxColorValue = 255), 6), rep(rgb(211,
        223, 239, maxColorValue = 255), 6), rep(rgb(156, 193, 222, maxColorValue = 255), 6), rep(rgb(95,
        150, 197, maxColorValue = 255), 6), rep(rgb(46, 98, 165, maxColorValue = 255), 6), rep(rgb(225,
        223, 239, maxColorValue = 255), 1), rep(rgb(181, 184, 213, maxColorValue = 255), 2), rep(rgb(134,
        132, 187, maxColorValue = 255), 3), rep(rgb(91, 67, 148, maxColorValue = 255), 5))
ns <- c(rep(c(1, 2, 3, 4, 5, "6+"), 2), rep(c(0, 1, 2, 3, 4, "5+"), 2), rep(c(1, 2, 3, 4, 5, "6+"), 4),
    rep(c(0, 1, 2, 3, 4, "5+"), 4), c(1, 1, 2, 1, 2, 3, 1, 2, 3, 4, "5+"))
pdf("figures/indel.profiles.samples.human.pdf", width = 12, height = 7.5/2, useDingbats = F)
par(tcl = 0.3, bty = "n", cex = 5/6)
for (cn in rns) {
    barplot(indel.sig_pd[, cn], col = tcs, border = "white", space = 0, las = 1, names = ns, las = 3,
        main = cn, ylab = "N of identified sindel")
    axis(1, c(6.5, 18.5, 36.5, 60.5, 77.5), label = c("homopolymer\nlength", "homopolymer\nlength", "number of\nrepeat units",
        "number of\nrepeat units", "microhomology\nlength"), lwd = 0, padj = 2, cex.axis = 0.5)
    axis(3, c(6.5, 18.5, 36.5, 60.5, 77.5), label = c("1bp deletion", "1bp insertion", ">1 bp deletion at repeats",
        ">1 bp insertion at repeats", "microhomology"), lwd = 0, padj = 2, cex.axis = 0.5)
    text(c(3.5, 9.5, 15.5, 21.5, 27.5, 33.5, 39.5, 45.5, 51.5, 57.5, 63.5, 69.5, 73, 74.5, 77, 81), max(indel.sig_pd[,
        cn]) * 0.95, label = c(rep(c("C", "T"), 2), rep(c(2, 3, 4, "5+"), 2), 2, 3, 4, "5+"), cex = 0.5,
        col = unique(tcs))
}
dev.off()

pdf("figures/indel.profiles.donors.human.pdf", width = 12, height = 7.5/2, useDingbats = F)
par(tcl = 0.3, bty = "n", cex = 5/6)
barplot(apply(indel.sig_pd[, rns], 1, sum), col = tcs, border = "white", space = 0, las = 1, names = ns,
    las = 3, main = "all donors", ylab = "N of identified sindel")
axis(1, c(6.5, 18.5, 36.5, 60.5, 77.5), label = c("homopolymer\nlength", "homopolymer\nlength", "number of\nrepeat units",
    "number of\nrepeat units", "microhomology\nlength"), lwd = 0, padj = 2, cex.axis = 0.5)
axis(3, c(6.5, 18.5, 36.5, 60.5, 77.5), label = c("1bp deletion", "1bp insertion", ">1 bp deletion at repeats",
    ">1 bp insertion at repeats", "microhomology"), lwd = 0, padj = 2, cex.axis = 0.5)
text(c(3.5, 9.5, 15.5, 21.5, 27.5, 33.5, 39.5, 45.5, 51.5, 57.5, 63.5, 69.5, 73, 74.5, 77, 81), max(apply(indel.sig_pd[,
    rns], 1, sum)) * 0.95, label = c(rep(c("C", "T"), 2), rep(c(2, 3, 4, "5+"), 2), 2, 3, 4, "5+"), cex = 0.5,
    col = unique(tcs))
for (s in unique(smps)) {
    cns <- which(smps == s)
    if (length(cns) > 1) {
        x <- apply(indel.sig_pd[, rns[cns]], 1, sum)
    } else {
        x <- indel.sig_pd[, rns[cns]]
    }
    barplot(x, col = tcs, border = "white", space = 0, las = 1, names = ns, las = 3, main = paste(s,
        "; ", meta_pd[rns[cns[1]], "age"], "yrs", sep = ""), ylab = "N of identified sindel")
    axis(1, c(6.5, 18.5, 36.5, 60.5, 77.5), label = c("homopolymer\nlength", "homopolymer\nlength", "number of\nrepeat units",
        "number of\nrepeat units", "microhomology\nlength"), lwd = 0, padj = 2, cex.axis = 0.5)
    axis(3, c(6.5, 18.5, 36.5, 60.5, 77.5), label = c("1bp deletion", "1bp insertion", ">1 bp deletion at repeats",
        ">1 bp insertion at repeats", "microhomology"), lwd = 0, padj = 2, cex.axis = 0.5)
    text(c(3.5, 9.5, 15.5, 21.5, 27.5, 33.5, 39.5, 45.5, 51.5, 57.5, 63.5, 69.5, 73, 74.5, 77, 81), max(x) *
        0.95, label = c(rep(c("C", "T"), 2), rep(c(2, 3, 4, "5+"), 2), 2, 3, 4, "5+"), cex = 0.5, col = unique(tcs))
}
dev.off()

# cosine similarity between human and dog SNV and indel----
library(MutationalPatterns)
cos_mat <- cos_sim_matrix(snv.sig[, cells_dog], snv.sig_pd[, cells_hs])
row.names(cos_mat) <- paste(row.names(cos_mat), "; ", out_mat_all[row.names(cos_mat), "age"], "yrs",
    sep = "")
colnames(cos_mat) <- paste(colnames(cos_mat), "; ", out_mat_all[colnames(cos_mat), "age"], "yrs", sep = "")
o1 <- order(out_mat_all[cells_dog, "age"])
o2 <- order(out_mat_all[cells_hs, "age"])
pheatmap(cos_mat[o1, o2], cluster_rows = F, cluster_cols = F, cellwidth = 8, cellheight = 8, fontsize = 8,
    legend = T, filename = "figures/mutation_profile_human_vs_dog/snv.cell.cos_sim.pdf")

pdf("figures/mutation_profile_human_vs_dog/overall.snv.profile.pdf", width = 12, height = 7.5 * 3/2,
    useDingbats = F)
p1 <- apply(snv.sig[, cells_dog], 1, sum)
p2 <- apply(snv.sig_pd[, cells_hs[out_mat_all[cells_hs, "age"] < 20]], 1, sum)
p3 <- apply(snv.sig_pd[, cells_hs], 1, sum)
par(mfrow = c(3, 1), tcl = 0.3, bty = "n", cex = 5/6)
fun_snv_profile(p1, main = "dog sSNVs", ylab = "%")
fun_snv_profile(p2, main = paste("human (<=20yrs) sSNVs; cosine similarity = ", round(cos_sim(p1, p2),
    2), sep = ""), ylab = "%")
fun_snv_profile(p3, main = paste("human sSNVs; cosine similarity = ", round(cos_sim(p1, p3), 2), sep = ""),
    ylab = "%")
dev.off()

cos_mat <- cos_sim_matrix(indel.sig[, cells_dog], indel.sig_pd[, cells_hs])
row.names(cos_mat) <- paste(row.names(cos_mat), "; ", out_mat_all[row.names(cos_mat), "age"], "yrs",
    sep = "")
colnames(cos_mat) <- paste(colnames(cos_mat), "; ", out_mat_all[colnames(cos_mat), "age"], "yrs", sep = "")
o1 <- order(out_mat_all[cells_dog, "age"])
o2 <- order(out_mat_all[cells_hs, "age"])
pheatmap(cos_mat[o1, o2], cluster_rows = F, cluster_cols = F, cellwidth = 8, cellheight = 8, fontsize = 8,
    legend = T, filename = "figures/mutation_profile_human_vs_dog/indel.cell.cos_sim.pdf")

pdf("figures/mutation_profile_human_vs_dog/overall.indel.profile.pdf", width = 12, height = 7.5 * 3/2,
    useDingbats = F)
p1 <- apply(indel.sig[, cells_dog], 1, sum)
p2 <- apply(indel.sig_pd[, cells_hs[out_mat_all[cells_hs, "age"] < 20]], 1, sum)
p3 <- apply(indel.sig_pd[, cells_hs], 1, sum)
par(mfrow = c(3, 1), tcl = 0.3, bty = "n", cex = 5/6)
fun_indel_profile(p1, main = "dog sindels", ylab = "n. sindel")
fun_indel_profile(p2, main = paste("human (<20yrs) sindels; cosine similarity = ", round(cos_sim(p1,
    p2), 2), sep = ""), ylab = "n. sindel")
fun_indel_profile(p3, main = paste("human sindels; cosine similarity = ", round(cos_sim(p1, p3), 2),
    sep = ""), ylab = "n. sindel")
dev.off()

p3 <- apply(indel.sig_pd[, setdiff(cells_hs, cells_hs_young)], 1, sum)
cos_sim(p1, p3)

# overall mutation profile for species/celltypes----
fun_plot_snv_profile <- function(sig.vector, ylab = "mutation", main = "", percent = T) {
    if (percent) {
        sig.vector <- sig.vector/sum(sig.vector) * 100
    }
    barplot(sig.vector, col = col_sig, border = "white", space = 0, names = names(sig.vector), las = 3,
        main = main, ylab = ylab)
    abline(v = c(1:5) * 16)
    text(c(1:6) * 16 - 8, max(sig.vector) * 19/20, label = c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G"))
}
fun_plot_indel_profile <- function(sig.vector, ylab = "mutation", main = "", percent = T) {
    if (percent) {
        sig.vector <- sig.vector/sum(sig.vector) * 100
    }
    tcs <- c(rep(rgb(244, 192, 122, maxColorValue = 255), 6), rep(rgb(238, 133, 51, maxColorValue = 255),
        6), rep(rgb(182, 217, 145, maxColorValue = 255), 6), rep(rgb(87, 157, 64, maxColorValue = 255),
        6), rep(rgb(243, 203, 183, maxColorValue = 255), 6), rep(rgb(237, 142, 112, maxColorValue = 255),
        6), rep(rgb(221, 83, 64, maxColorValue = 255), 6), rep(rgb(172, 45, 38, maxColorValue = 255),
        6), rep(rgb(211, 223, 239, maxColorValue = 255), 6), rep(rgb(156, 193, 222, maxColorValue = 255),
        6), rep(rgb(95, 150, 197, maxColorValue = 255), 6), rep(rgb(46, 98, 165, maxColorValue = 255),
        6), rep(rgb(225, 223, 239, maxColorValue = 255), 1), rep(rgb(181, 184, 213, maxColorValue = 255),
        2), rep(rgb(134, 132, 187, maxColorValue = 255), 3), rep(rgb(91, 67, 148, maxColorValue = 255),
        5))
    ns <- c(rep(c(1, 2, 3, 4, 5, "6+"), 2), rep(c(0, 1, 2, 3, 4, "5+"), 2), rep(c(1, 2, 3, 4, 5, "6+"),
        4), rep(c(0, 1, 2, 3, 4, "5+"), 4), c(1, 1, 2, 1, 2, 3, 1, 2, 3, 4, "5+"))
    barplot(sig.vector, col = tcs, border = "white", space = 0, las = 1, names = ns, las = 3, main = main,
        ylab = ylab)
    axis(1, c(6.5, 18.5, 36.5, 60.5, 77.5), label = c("homopolymer\nlength", "homopolymer\nlength", "number of\nrepeat units",
        "number of\nrepeat units", "microhomology\nlength"), lwd = 0, padj = 2, cex.axis = 0.5)
    axis(3, c(6.5, 18.5, 36.5, 60.5, 77.5), label = c("1bp deletion", "1bp insertion", ">1 bp deletion at repeats",
        ">1 bp insertion at repeats", "microhomology"), lwd = 0, padj = 2, cex.axis = 0.5)
    text(c(3.5, 9.5, 15.5, 21.5, 27.5, 33.5, 39.5, 45.5, 51.5, 57.5, 63.5, 69.5, 73, 74.5, 77, 81), max(sig.vector) *
        0.95, label = c(rep(c("C", "T"), 2), rep(c(2, 3, 4, "5+"), 2), 2, 3, 4, "5+"), cex = 0.5, col = unique(tcs))
}

pdf("figures/overall_mutation_spectrum.pdf", width = 8, height = 3, useDingbats = F)
par(mar = c(6, 4, 3, 1), tcl = 0.3, cex = 5/6)
fun_plot_snv_profile(apply(snv.sig[, cells_dog], 1, sum), main = "dog neurons", ylab = "%")
fun_plot_snv_profile(apply(snv.sig_pd[, cells_hs], 1, sum), main = "human neurons", ylab = "%")
fun_plot_snv_profile(apply(snv.sig_pd[, cells_hs_young], 1, sum), main = "human neurons (\342\211\24420yrs)",
    ylab = "%")
fun_plot_snv_profile(apply(snv.sig_pd[, setdiff(cells_hs, cells_hs_young)], 1, sum), main = "human neurons (>20yrs)",
    ylab = "%")
fun_plot_snv_profile(apply(snv.sig.crypt[, cells_crypt_dog], 1, sum), main = "dog crypt", ylab = "%")
fun_plot_snv_profile(apply(snv.sig.crypt[, cells_crypt_human], 1, sum), main = "human crypt", ylab = "%")
fun_plot_indel_profile(apply(indel.sig[, cells_dog], 1, sum), main = "dog neurons", ylab = "%")
fun_plot_indel_profile(apply(indel.sig_pd[, cells_hs], 1, sum), main = "human neurons", ylab = "%")
fun_plot_indel_profile(apply(indel.sig_pd[, cells_hs_young], 1, sum), main = "human neurons (\342\211\24420yrs)",
    ylab = "%")
fun_plot_indel_profile(apply(indel.sig_pd[, setdiff(cells_hs, cells_hs_young)], 1, sum), main = "human neurons (>20yrs)",
    ylab = "%")
fun_plot_indel_profile(apply(indel.sig.crypt[, cells_crypt_dog], 1, sum), main = "dog crypt", ylab = "%")
fun_plot_indel_profile(apply(indel.sig.crypt[, cells_crypt_human], 1, sum), main = "human crypt", ylab = "%")
dev.off()

# signature factorization for snv----
library(NMF)
library(MutationalPatterns)
mut_mat <- cbind(snv.sig[, cells_used], snv.sig_pd[, cells_normal])
mut_mat <- apply(mut_mat, 2, function(x) {
    return(x/sum(x) * 100 + 0.001)
})
# mut_mat=cbind(snv.sig[,cells_used],snv.sig_pd[,cells_normal])+0.01
estimate <- nmf(mut_mat, rank = 1:8, method = "brunet", nrun = 20, seed = 123456, .opt = "v-p")
pdf("figures/signature_NMF_rankNumber_SNV.pdf", width = 7.5, height = 7.5, useDingbats = F)
plot(estimate)
dev.off()
nmf_res <- extract_signatures(mut_mat, rank = 2, nrun = 10)
colnames(nmf_res$signatures) <- c("signature A2", "signature A1")
row.names(nmf_res$contribution) <- c("signature A2", "signature A1")
saveRDS(nmf_res, "tables/dnaseq/scan2/signatures.snv.rds")
#nmf_res_snv=nmf_res
nmf_res <- nmf_res_snv

pdf("figures/signature_NMF_SNV.pdf", width = 12, height = 7.5, useDingbats = F)
par(mfrow = c(ncol(nmf_res$signatures), 1), tcl = 0.3, bty = "n", cex = 5/6)
tcs <- col_sig
for (i in colnames(nmf_res$signatures)) {
    barplot(nmf_res_snv$signatures[, i]/sum(nmf_res_snv$signatures[, i]) * 100, col = tcs, border = "white",
        space = 0, names = row.names(snv.sig), las = 3, main = i, ylab = "%")
    abline(v = c(1:5) * 16)
    text(c(1:6) * 16 - 8, max(nmf_res_snv$signatures[, i]) * 19/20, label = c("C>A", "C>G", "C>T", "T>A",
        "T>C", "T>G"))
}
par(mar = c(10, 4, 1, 1), mfrow = c(1, 1), tcl = 0.3, bty = "n", cex = 5/6)
o1 <- cells_used[order(samples_meta[cells_used, "age"])]
o2 <- cells_normal[order(meta_pd[cells_normal, "age"])]
t <- apply(nmf_res$contribution[, c(o1, o2)] * apply(nmf_res$signatures, 2, sum), 2, sum)
barplot(nmf_res$contribution[, c(o1, o2)] * apply(nmf_res$signatures, 2, sum)/rbind(t, t), col = cs_gg[1:2],
    ylab = "proportion of identified sSNVs", las = 3, names = c(paste(samples_meta[o1, "donor"], samples_meta[o1,
        c("age")], sep = "; "), paste(meta_pd[o2, "donor"], meta_pd[o2, c("age")], sep = "; ")), border = cs_gg[1:2],
    space = c(rep(0, length(o1)), 2, rep(0, length(o2) - 1)))
legend("topleft", pch = 15, col = cs_gg[1:2], legend = colnames(nmf_res_snv$signatures), bty = "n")
dev.off()

# signature factorization for indel----
require(NMF)
require(MutationalPatterns)
mut_mat <- cbind(indel.sig[, cells_used], indel.sig_pd[, cells_normal])
mut_mat <- apply(mut_mat, 2, function(x) {
    return(x/sum(x) * 100 + 0.001)
})
mut_mat[is.na(mut_mat)] <- 0.001
estimate <- nmf(mut_mat, rank = 1:8, method = "brunet", nrun = 20, seed = 123456, .opt = "v-p")
pdf("figures/signature_NMF_rankNumber_INDEL.pdf", width = 7.5, height = 7.5, useDingbats = F)
plot(estimate)
dev.off()

nmf_res <- extract_signatures(mut_mat, rank = 2, nrun = 10)
colnames(nmf_res$signatures) <- c("signature ID-A", "signature ID-B")
row.names(nmf_res$contribution) <- c("signature ID-A", "signature ID-B")
saveRDS(nmf_res, "tables/dnaseq/scan2/signatures.indel.rds")
nmf_res <- nmf_res_indel
pdf("figures/signature_NMF_INDEL.pdf", width = 12, height = 7.5, useDingbats = F)
par(mfrow = c(ncol(nmf_res_indel$signatures), 1), tcl = 0.3, bty = "n", cex = 5/6)
tcs <- c(rep(rgb(244, 192, 122, maxColorValue = 255), 6), rep(rgb(238, 133, 51, maxColorValue = 255),
    6), rep(rgb(182, 217, 145, maxColorValue = 255), 6), rep(rgb(87, 157, 64, maxColorValue = 255), 6),
    rep(rgb(243, 203, 183, maxColorValue = 255), 6), rep(rgb(237, 142, 112, maxColorValue = 255), 6),
    rep(rgb(221, 83, 64, maxColorValue = 255), 6), rep(rgb(172, 45, 38, maxColorValue = 255), 6), rep(rgb(211,
        223, 239, maxColorValue = 255), 6), rep(rgb(156, 193, 222, maxColorValue = 255), 6), rep(rgb(95,
        150, 197, maxColorValue = 255), 6), rep(rgb(46, 98, 165, maxColorValue = 255), 6), rep(rgb(225,
        223, 239, maxColorValue = 255), 1), rep(rgb(181, 184, 213, maxColorValue = 255), 2), rep(rgb(134,
        132, 187, maxColorValue = 255), 3), rep(rgb(91, 67, 148, maxColorValue = 255), 5))
ns <- c(rep(c(1, 2, 3, 4, 5, "6+"), 2), rep(c(0, 1, 2, 3, 4, "5+"), 2), rep(c(1, 2, 3, 4, 5, "6+"), 4),
    rep(c(0, 1, 2, 3, 4, "5+"), 4), c(1, 1, 2, 1, 2, 3, 1, 2, 3, 4, "5+"))
for (i in 1:ncol(nmf_res_indel$signatures)) {
    barplot(nmf_res_indel$signatures[, i]/sum(nmf_res_indel$signatures[, i]) * 100, col = tcs, border = "white",
        space = 0, las = 1, names = ns, las = 3, main = colnames(nmf_res_indel$signatures)[i], ylab = "%")
    axis(1, c(6.5, 18.5, 36.5, 60.5, 77.5), label = c("homopolymer\nlength", "homopolymer\nlength", "number of\nrepeat units",
        "number of\nrepeat units", "microhomology\nlength"), lwd = 0, padj = 2, cex.axis = 0.5)
    axis(3, c(6.5, 18.5, 36.5, 60.5, 77.5), label = c("1bp deletion", "1bp insertion", ">1 bp deletion at repeats",
        ">1 bp insertion at repeats", "microhomology"), lwd = 0, padj = 2, cex.axis = 0.5)
    text(c(3.5, 9.5, 15.5, 21.5, 27.5, 33.5, 39.5, 45.5, 51.5, 57.5, 63.5, 69.5, 73, 74.5, 77, 81), max(nmf_res_indel$signatures[,
        i]) * 0.95, label = c(rep(c("C", "T"), 2), rep(c(2, 3, 4, "5+"), 2), 2, 3, 4, "5+"), cex = 0.5,
        col = unique(tcs))
}
par(mar = c(10, 4, 1, 1), mfrow = c(1, 1), tcl = 0.3, bty = "n", cex = 5/6)
o1 <- cells_used[order(samples_meta[cells_used, "age"])]
o2 <- cells_normal[order(meta_pd[cells_normal, "age"])]
t <- apply(nmf_res_indel$contribution[, c(o1, o2)] * apply(nmf_res_indel$signatures, 2, sum), 2, sum)
barplot(nmf_res_indel$contribution[, c(o1, o2)] * apply(nmf_res_indel$signatures, 2, sum)/rbind(t, t),
    col = cs_gg[1:2], ylab = "proportion of identified sSNVs", las = 3, names = c(paste(samples_meta[o1,
        "donor"], samples_meta[o1, c("age")], sep = "; "), paste(meta_pd[o2, "donor"], meta_pd[o2, c("age")],
        sep = "; ")), border = cs_gg[1:2], space = c(rep(0, length(o1)), 2, rep(0, length(o2) - 1)))
dev.off()

# signature factorization for snv and indel after combining crypt and neuron data----
# do SNV first----
# downsample number of neurons to 59 for human, as we don't have a lot of crypt data
tmp_cells_human <- sample(cells_normal, 58)
mut_mat <- cbind(snv.sig[, cells_used], snv.sig_pd[, tmp_cells_human], snv.sig.crypt)
mut_mat <- apply(mut_mat, 2, function(x) {
    return(x/sum(x) * 100 + 0.001)
})
mut_mat[is.na(mut_mat)] <- 0.001
mut_mat_all <- cbind(snv.sig[, cells_used], snv.sig_pd[, cells_human], snv.sig.crypt)
estimate <- nmf(mut_mat, rank = 1:8, method = "brunet", nrun = 20, seed = 123456, .opt = "v-p")
pdf("figures/signature_NMF_rankNumber_SNV_with_crypt.pdf", width = 7.5, height = 7.5, useDingbats = F)
plot(estimate)
dev.off()

nmf_res_snv_with_crypt <- extract_signatures(mut_mat, rank = 3, nrun = 10)
colnames(nmf_res_snv_with_crypt$signatures) <- c("signature A2", "signature A3", "signature A1")
row.names(nmf_res_snv_with_crypt$contribution) <- c("signature A2", "signature A3", "signature A1")
nmf_res_snv_with_crypt$contribution <- fit_to_signatures(mut_mat_all, nmf_res_snv_with_crypt$signatures)$contribution
saveRDS(nmf_res_snv_with_crypt, "tables/dnaseq/scan2/signatures.snv_with_crypt.rds")

pdf("figures/signature_NMF_SNV_with_scrypt.pdf", width = 12, height = 7.5 * 3/2, useDingbats = F)
par(mfrow = c(ncol(nmf_res_snv_with_crypt$signatures), 1), tcl = 0.3, bty = "n", cex = 5/6)
tcs <- col_sig
for (i in colnames(nmf_res_snv_with_crypt$signatures)) {
    barplot(nmf_res_snv_with_crypt$signatures[, i]/sum(nmf_res_snv_with_crypt$signatures[, i]) * 100,
        col = tcs, border = "white", space = 0, names = row.names(snv.sig), las = 3, main = i, ylab = "%")
    abline(v = c(1:5) * 16)
    text(c(1:6) * 16 - 8, max(nmf_res_snv_with_crypt$signatures[, i]) * 19/20, label = c("C>A", "C>G",
        "C>T", "T>A", "T>C", "T>G"))
}
par(mar = c(10, 4, 1, 1), mfrow = c(2, 1), tcl = 0.3, bty = "n", cex = 5/6)
o1 <- cells_used[order(samples_meta[cells_used, "age"])]
o2 <- cells_normal[order(meta_pd[cells_normal, "age"])]
tmp_cells <- c(o1, o2, cells_crypt_dog, cells_crypt_human)
t <- nmf_res_snv_with_crypt$contribution[, tmp_cells] * apply(nmf_res_snv_with_crypt$signatures, 2, sum)
sig_expo <- t(t)/apply(t, 2, sum) * c(out_mat_all[c(o1, o2), "sSNV_burden"]/len_human_autosome, meta_crypt[cells_crypt,
    "Substitution.burden"]/meta_crypt[cells_crypt, "Total.genome.size..bp."]/2 * 1e+09)
sig_expo_snv <- sig_expo
barplot(t(sig_expo[, names(col_sbs)]), col = col_sbs, ylab = "signature exposure / Gb", las = 3, names = c(paste(samples_meta[o1,
    "donor"], samples_meta[o1, c("age")], sep = "; "), paste(meta_pd[o2, "donor"], meta_pd[o2, c("age")],
    sep = "; "), paste(meta_crypt[cells_crypt, "Individual"], round(meta_crypt[cells_crypt, c("Age..years.")],
    1), sep = "; ")), border = col_sbs, space = c(rep(0, length(o1)), 2, rep(0, length(o2) - 1), 2, rep(0,
    length(cells_crypt_dog) - 1), 2, rep(0, length(cells_crypt_human) - 1)), cex.names = 0.3)
legend("topleft", pch = 15, col = col_sbs, legend = names(col_sbs), bty = "n")
dev.off()

# then do indel----
# downsample number of neurons to 59 for human, as we don't have a lot of crypt data
tmp_cells_human <- sample(cells_normal, 58)
mut_mat <- cbind(indel.sig[, cells_used], indel.sig_pd[, tmp_cells_human], indel.sig.crypt)
mut_mat <- apply(mut_mat, 2, function(x) {
    return(x/sum(x) * 100 + 0.001)
})
mut_mat[is.na(mut_mat)] <- 0.001
mut_mat_all <- cbind(indel.sig[, cells_used], indel.sig_pd[, cells_human], indel.sig.crypt)
estimate <- nmf(mut_mat, rank = 1:8, method = "brunet", nrun = 20, seed = 123456, .opt = "v-p")
pdf("figures/signature_NMF_rankNumber_indel_with_crypt.pdf", width = 7.5, height = 7.5, useDingbats = F)
plot(estimate)
dev.off()

nmf_res_indel_with_crypt <- extract_signatures(mut_mat, rank = 3, nrun = 10)
colnames(nmf_res_indel_with_crypt$signatures) <- c("signature ID-C", "signature ID-A", "signature ID-B")
row.names(nmf_res_indel_with_crypt$contribution) <- c("signature ID-C", "signature ID-A", "signature ID-B")
nmf_res_indel_with_crypt$contribution <- fit_to_signatures(mut_mat_all, nmf_res_indel_with_crypt$signatures)$contribution
saveRDS(nmf_res_indel_with_crypt, "tables/dnaseq/scan2/signatures.indel_with_crypt.rds")

pdf("figures/signature_NMF_indel_with_scrypt.pdf", width = 12, height = 7.5 * 3/2, useDingbats = F)
par(mfrow = c(ncol(nmf_res_indel_with_crypt$signatures), 1), tcl = 0.3, bty = "n", cex = 5/6)
tcs <- c(rep(rgb(244, 192, 122, maxColorValue = 255), 6), rep(rgb(238, 133, 51, maxColorValue = 255),
    6), rep(rgb(182, 217, 145, maxColorValue = 255), 6), rep(rgb(87, 157, 64, maxColorValue = 255), 6),
    rep(rgb(243, 203, 183, maxColorValue = 255), 6), rep(rgb(237, 142, 112, maxColorValue = 255), 6),
    rep(rgb(221, 83, 64, maxColorValue = 255), 6), rep(rgb(172, 45, 38, maxColorValue = 255), 6), rep(rgb(211,
        223, 239, maxColorValue = 255), 6), rep(rgb(156, 193, 222, maxColorValue = 255), 6), rep(rgb(95,
        150, 197, maxColorValue = 255), 6), rep(rgb(46, 98, 165, maxColorValue = 255), 6), rep(rgb(225,
        223, 239, maxColorValue = 255), 1), rep(rgb(181, 184, 213, maxColorValue = 255), 2), rep(rgb(134,
        132, 187, maxColorValue = 255), 3), rep(rgb(91, 67, 148, maxColorValue = 255), 5))
ns <- c(rep(c(1, 2, 3, 4, 5, "6+"), 2), rep(c(0, 1, 2, 3, 4, "5+"), 2), rep(c(1, 2, 3, 4, 5, "6+"), 4),
    rep(c(0, 1, 2, 3, 4, "5+"), 4), c(1, 1, 2, 1, 2, 3, 1, 2, 3, 4, "5+"))
for (i in 1:ncol(nmf_res_indel_with_crypt$signatures)) {
    barplot(nmf_res_indel_with_crypt$signatures[, i]/sum(nmf_res_indel_with_crypt$signatures[, i]) *
        100, col = tcs, border = "white", space = 0, las = 1, names = ns, las = 3, main = colnames(nmf_res_indel_with_crypt$signatures)[i],
        ylab = "%")
    axis(1, c(6.5, 18.5, 36.5, 60.5, 77.5), label = c("homopolymer\nlength", "homopolymer\nlength", "number of\nrepeat units",
        "number of\nrepeat units", "microhomology\nlength"), lwd = 0, padj = 2, cex.axis = 0.5)
    axis(3, c(6.5, 18.5, 36.5, 60.5, 77.5), label = c("1bp deletion", "1bp insertion", ">1 bp deletion at repeats",
        ">1 bp insertion at repeats", "microhomology"), lwd = 0, padj = 2, cex.axis = 0.5)
    text(c(3.5, 9.5, 15.5, 21.5, 27.5, 33.5, 39.5, 45.5, 51.5, 57.5, 63.5, 69.5, 73, 74.5, 77, 81), max(nmf_res_indel_with_crypt$signatures[,
        i]) * 0.95, label = c(rep(c("C", "T"), 2), rep(c(2, 3, 4, "5+"), 2), 2, 3, 4, "5+"), cex = 0.5,
        col = unique(tcs))
}
par(mar = c(10, 4, 1, 1), mfrow = c(2, 1), tcl = 0.3, bty = "n", cex = 5/6)
o1 <- cells_used[order(samples_meta[cells_used, "age"])]
o2 <- cells_normal[order(meta_pd[cells_normal, "age"])]
tmp_cells <- c(o1, o2, cells_crypt_dog, cells_crypt_human)
t <- nmf_res_indel_with_crypt$contribution[, tmp_cells] * apply(nmf_res_indel_with_crypt$signatures,
    2, sum)
sig_expo <- t(t)/apply(t, 2, sum) * c(out_mat_all[c(o1, o2), "sINDEL_burden"]/len_human_autosome, meta_crypt[cells_crypt,
    "Indel.burden"]/meta_crypt[cells_crypt, "Total.genome.size..bp."]/2 * 1e+09)
barplot(t(sig_expo[, names(col_id)]), col = col_id, ylab = "signature exposure / Gb", las = 3, names = c(paste(samples_meta[o1,
    "donor"], samples_meta[o1, c("age")], sep = "; "), paste(meta_pd[o2, "donor"], meta_pd[o2, c("age")],
    sep = "; "), paste(meta_crypt[cells_crypt, "Individual"], round(meta_crypt[cells_crypt, c("Age..years.")],
    1), sep = "; ")), border = col_id, space = c(rep(0, length(o1)), 2, rep(0, length(o2) - 1), 2, rep(0,
    length(cells_crypt_dog) - 1), 2, rep(0, length(cells_crypt_human) - 1)), cex.names = 0.3)
legend("topleft", pch = 15, col = col_id, legend = names(col_id), bty = "n")
dev.off()

write.table(cbind(sig_expo_snv, sig_expo), "tables/For_manuscript/signature_exposure_perGb_with_crypt.txt",
    sep = "\t", quote = F)

# test signature robustness----
library(NMF)
library(MutationalPatterns)
# SNV signatures
mut_mat <- cbind(snv.sig[, cells_used], snv.sig_pd[, cells_normal])
mut_mat <- apply(mut_mat, 2, function(x) {
    return(x/sum(x) * 100 + 0.001)
})
# mut_mat=cbind(snv.sig[,cells_used],snv.sig_pd[,cells_normal])+0.01
res <- list()
for (i in 1:10) {
    res[[i]] <- extract_signatures(mut_mat, rank = 2, nrun = 10, seed = i)$signatures
}
sig_snv_robust <- matrix(0, 10, 2)
colnames(sig_snv_robust) <- colnames(nmf_res_snv$signatures)
for (i in 1:10) {
    t <- cos_sim_matrix(nmf_res_snv$signatures, res[[i]])
    sig_snv_robust[i, ] <- c(max(t[1, ]), max(t[2, ]))
}
# indel signatures
mut_mat <- cbind(indel.sig[, cells_used], indel.sig_pd[, cells_normal])
mut_mat <- apply(mut_mat, 2, function(x) {
    return(x/sum(x) * 100 + 0.001)
})
mut_mat[is.na(mut_mat)] <- 0.001
# mut_mat=cbind(snv.sig[,cells_used],snv.sig_pd[,cells_normal])+0.01
res <- list()
for (i in 1:10) {
    res[[i]] <- extract_signatures(mut_mat, rank = 2, nrun = 10, seed = i)$signatures
}
sig_indel_robust <- matrix(0, 10, 2)
colnames(sig_indel_robust) <- colnames(nmf_res_indel$signatures)
for (i in 1:10) {
    t <- cos_sim_matrix(nmf_res_indel$signatures, res[[i]])
    sig_indel_robust[i, ] <- c(max(t[1, ]), max(t[2, ]))
}
# plot
tdf <- cbind(sig_snv_robust[, 2:1], sig_indel_robust)
tl <- list()
for (c in colnames(tdf)) {
    tl[[c]] <- tdf[, c]
}
pdf("figures/signature_robustness.pdf", width = 3.5, height = 3, useDingbats = F)
par(mar = c(2, 4, 3, 1), tcl = 0.3, bty = "n")
fun_boxplot(tl, ylab = "cosine similarity", main = "signature robustness", ylim = c(0.7, 1), points = T,
    cex = 1.2, col = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3"))
axis(1, 1:4, label = c("A1", "A2", "ID-A", "ID-B"), lwd = 0)
dev.off()

# test dog-human normalization robustness in signature analysis----
# directly get signature from dog data doesn't performed well: cos_sim = 0.89 and 0.85 for A1 and A1
# I think it is because of the sample size.
# now try using down-sample method
library(NMF)
library(MutationalPatterns)
# SNV signatures
mut_mat <- cbind(snv.sig[, cells_used], snv.sig_pd[, cells_normal])
mut_mat <- apply(mut_mat, 2, function(x) {
    return(x/sum(x) * 100 + 0.001)
})
res <- list()
for (i in 1:10) {
    resampled_cells <- c(cells_dog, sample(cells_normal, length(cells_dog)))
    res[[i]] <- extract_signatures(mut_mat[, resampled_cells], rank = 2, nrun = 10)$signatures
}
sig_snv_downsample <- matrix(0, 10, 2)
colnames(sig_snv_downsample) <- colnames(nmf_res_snv$signatures)
for (i in 1:10) {
    t <- cos_sim_matrix(nmf_res_snv$signatures, res[[i]])
    sig_snv_downsample[i, ] <- c(max(t[1, ]), max(t[2, ]))
}
# indel signatures
mut_mat <- cbind(indel.sig[, cells_used], indel.sig_pd[, cells_normal])
mut_mat <- apply(mut_mat, 2, function(x) {
    return(x/sum(x) * 100 + 0.001)
})
mut_mat[is.na(mut_mat)] <- 0.001
res <- list()
for (i in 1:10) {
    resampled_cells <- c(cells_dog, sample(cells_normal, length(cells_dog)))
    res[[i]] <- extract_signatures(mut_mat[, resampled_cells], rank = 2, nrun = 10)$signatures
}
sig_indel_downsample <- matrix(0, 10, 2)
colnames(sig_indel_downsample) <- colnames(nmf_res_indel$signatures)
for (i in 1:10) {
    t <- cos_sim_matrix(nmf_res_indel$signatures, res[[i]])
    sig_indel_downsample[i, ] <- c(max(t[1, ]), max(t[2, ]))
}
# plot
tdf <- cbind(sig_snv_downsample[, 2:1], sig_indel_downsample)
tl <- list()
for (c in colnames(tdf)) {
    tl[[c]] <- tdf[, c]
}
pdf("figures/signature_downsample_robustness.pdf", width = 3.5, height = 3, useDingbats = F)
par(mar = c(2, 4, 3, 1), tcl = 0.3, bty = "n")
fun_boxplot(tl, ylab = "cosine similarity", main = "signature robustness after downsampling", ylim = c(0.7,
    1), points = T, cex = 1.2, col = c("#66c2a5", "#fc8d62", "#8da0cb", "#e78ac3"))
axis(1, 1:4, label = c("A1", "A2", "ID-A", "ID-B"), lwd = 0)
dev.off()

# snv and indel burden for signatures; after factorization----
pdf("figures/mutation_burden_signatures.pdf", width = 6, height = 4.5, useDingbats = F)
par(mar = c(4, 4, 1, 1), tcl = 0.3, bty = "n")
for (s in c(colnames(nmf_res_snv$signatures), colnames(nmf_res_indel$signatures))) {
    r <- round(cor(out_mat_all[cells_normal, "age"], out_mat_all[cells_normal, s]), 2)
    plot(out_mat_all[cells_normal, "age"], out_mat_all[cells_normal, s], bg = cs_species["human"], pch = 21,
        col = "white", lwd = 0.3, ylab = "sSNVs per neuron", main = s, xlab = "age (yrs)", cex = 1.5,
        xlim = c(0, 110))
    legend("topleft", pch = 19, legend = names(cs_species), bty = "n", col = cs_species)
    lreg <- lm(out_mat_all[cells_normal, s] ~ out_mat_all[cells_normal, "age"])
    abline(lreg, col = cs_species["human"])
    text(110, 100, pos = 2, col = cs_species["human"], label = paste(round(lreg$coefficients[2], 1),
        " sSNVs / year; r = ", round(r, 2), sep = ""))
    r <- round(cor(out_mat_all[cells_used, "age"], out_mat_all[cells_used, s]), 2)
    points(out_mat_all[cells_used, "age"], out_mat_all[cells_used, s], bg = cs_species["dog"], pch = 21,
        col = "white", lwd = 0.3, cex = 1.5)
    lreg <- lm(out_mat_all[cells_used, s] ~ out_mat_all[cells_used, "age"])
    abline(lreg, col = cs_species["dog"])
    text(110, 200, pos = 2, col = cs_species["dog"], label = paste(round(lreg$coefficients[2], 1), " sSNVs / year; r = ",
        round(r, 2), sep = ""))
    # zoom-in version
    r <- round(cor(out_mat_all[cells_hs_young, "age"], out_mat_all[cells_hs_young, s]), 2)
    plot(out_mat_all[cells_hs_young, "age"], out_mat_all[cells_hs_young, s], bg = cs_species["human"],
        pch = 21, col = "white", lwd = 0.3, ylab = "sSNVs per neuron", main = s, xlab = "age (yrs)",
        cex = 1.5, xlim = c(0, 30))
    legend("topleft", pch = 19, legend = names(cs_species), bty = "n", col = cs_species)
    lreg <- lm(out_mat_all[cells_hs_young, s] ~ out_mat_all[cells_hs_young, "age"])
    abline(lreg, col = cs_species["human"])
    text(30, 100, pos = 2, col = cs_species["human"], label = paste(round(lreg$coefficients[2], 1), " sSNVs / year; r = ",
        round(r, 2), sep = ""))
    r <- round(cor(out_mat_all[cells_used, "age"], out_mat_all[cells_used, s]), 2)
    points(out_mat_all[cells_used, "age"], out_mat_all[cells_used, s], bg = cs_species["dog"], pch = 21,
        col = "white", lwd = 0.3, cex = 1.5)
    lreg <- lm(out_mat_all[cells_used, s] ~ out_mat_all[cells_used, "age"])
    abline(lreg, col = cs_species["dog"])
    text(30, 200, pos = 2, col = cs_species["dog"], label = paste(round(lreg$coefficients[2], 1), " sSNVs / year; r = ",
        round(r, 2), sep = ""))
}
dev.off()

# genes differentially expressed in dog during aging----
library(DESeq2)
library(pheatmap)
# deseq_res=fun_deseq(dog_aging.count,which(colnames(dog_aging.count)%in%dog_aging.old),
#                     which(colnames(dog_aging.count)%in%dog_aging.young))
# write.table(deseq_res,'tables/rnaseq/dog_aging.deseq.txt',sep='\\t',quote=F)
deseq_res <- read.table("tables/rnaseq/dog_aging.deseq.txt", row.names = 1, header = T)
deseq_downgenes <- row.names(deseq_res)[which(deseq_res$log2FoldChange < (-0.585) & deseq_res$padj <
    0.05)]
deseq_upgenes <- row.names(deseq_res)[which(deseq_res$log2FoldChange > (0.585) & deseq_res$padj < 0.05)]
deseq_downgenes2 <- row.names(deseq_res)[which(deseq_res$log2FoldChange < 0 & deseq_res$padj < 0.05)]
deseq_upgenes2 <- row.names(deseq_res)[which(deseq_res$log2FoldChange > 0 & deseq_res$padj < 0.05)]

pheatmap(dog_aging.exp[deseq_downgenes2, ], scale = "row", cluster_rows = F, cluster_cols = F, filename = "figures/gene_expression_downregulated.pdf",
    cellwidth = 16, cellheight = 1, fontsize_row = 1, fontsize_col = 16)
pheatmap(dog_aging.exp[deseq_upgenes2, ], scale = "row", cluster_rows = F, cluster_cols = F, filename = "figures/gene_expression_upregulated.pdf",
    cellwidth = 16, cellheight = 1, fontsize_row = 1, fontsize_col = 16)
pheatmap(log10(dog_aging.exp[deseq_downgenes2, ] + 1), cluster_rows = T, cluster_cols = F, filename = "figures/gene_expression_downregulated.unnormalized.pdf",
    cellwidth = 16, cellheight = 1, fontsize_row = 1, fontsize_col = 16)
pheatmap(log10(dog_aging.exp[deseq_upgenes2, ] + 1), cluster_rows = T, cluster_cols = F, filename = "figures/gene_expression_upregulated.unnormalized.pdf",
    cellwidth = 16, cellheight = 1, fontsize_row = 1, fontsize_col = 16)

# genes differentially expressed in human during aging; analyzed from GTEx data----
gtex.count <- read.table("tables/rnaseq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_reads.Frontal_Cortext.gct",
    header = T, row.names = 1, check.names = F)
gtex.meta <- read.table("tables/rnaseq/metadata3.Frontal_Cortex.txt", header = F, row.names = 1, sep = "\t")
# which(gtex.meta[colnames(gtex.count),3]<=40)
# deseq_gtex=fun_deseq(gtex.count,which(gtex.meta[colnames(gtex.count),3]>=60),
#           which(gtex.meta[colnames(gtex.count),3]<=40))
# write.table(deseq_gtex,'tables/rnaseq/GTEx.deseq.txt',sep='\\t',quote=F)
deseq_gtex <- read.table("tables/rnaseq/GTEx.deseq.txt", header = T, row.names = 1)
deseq_gtex_downgenes <- row.names(deseq_gtex)[which(deseq_gtex$log2FoldChange < (-0.585) & deseq_gtex$padj <
    0.05)]
deseq_gtex_upgenes <- row.names(deseq_gtex)[which(deseq_gtex$log2FoldChange > (0.585) & deseq_gtex$padj <
    0.05)]

gtex.tpm <- read.table("tables/rnaseq/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_tpm.Frontal_Cortext.gct",
    header = T, row.names = 1, check.names = F)
gtex.tpm <- gtex.tpm[, order(gtex.meta[colnames(gtex.tpm), 3])]
pheatmap(gtex.tpm[deseq_gtex_downgenes, ], scale = "row", cluster_rows = F, cluster_cols = F, filename = "figures/gene_expression_GTEx_downregulated.pdf",
    cellwidth = 2, cellheight = 1, fontsize_row = 1, fontsize_col = 2, labels_col = gtex.meta[colnames(gtex.tpm),
        3])
pheatmap(gtex.tpm[deseq_gtex_upgenes, ], scale = "row", cluster_rows = F, cluster_cols = F, filename = "figures/gene_expression_GTEx_upregulated.pdf",
    cellwidth = 2, cellheight = 1, fontsize_row = 1, fontsize_col = 2, labels_col = gtex.meta[colnames(gtex.tpm),
        3])


# gene length and aging genes----
pdf("figures/gene_length_aging.pdf", width = 5, height = 5, useDingbats = F)
par(bty = "n", tcl = 0.3)
plot(density(log10(gene.len.dog[rowMeans(dog_aging.exp[row.names(gene.len.dog), ]) > 1, 1]), na.rm = T),
    xaxt = "n", xlim = c(1, 7), main = "dog bulk data", xlab = "", ylim = c(0, 0.8), col = "grey", lwd = 2)
lines(density(log10(gene.len.dog[deseq_upgenes2, 1]), na.rm = T), col = cs_gg[1], lwd = 2)
lines(density(log10(gene.len.dog[deseq_downgenes2, 1]), na.rm = T), col = cs_gg[2], lwd = 2)
legend("topleft", lty = 1, col = c("grey", cs_gg[1:2]), lwd = 3, legend = c("expressed genes", "up during aging",
    "down during aging"), bty = "n")
add_axis(1)
plot(density(log10(gene.len.dog[rowMeans(dog_aging.exp[row.names(gene.len.dog), ]) > 1, 1]), na.rm = T),
    xaxt = "n", xlim = c(1, 7), main = "dog snRNA (excitatory neuron) data", xlab = "", ylim = c(0, 1.2),
    col = "grey", lwd = 2)
lines(density(log10(gene.len.dog[degs.up.dog, 1]), na.rm = T), col = cs_gg[1], lwd = 2)
lines(density(log10(gene.len.dog[degs.down.dog, 1]), na.rm = T), col = cs_gg[2], lwd = 2)
legend("topleft", lty = 1, col = c("grey", cs_gg[1:2]), lwd = 3, legend = c("expressed genes", "up during aging",
    "down during aging"), bty = "n")
add_axis(1)
plot(density(log10(gene.len.human[rowMeans(gtex.tpm[row.names(gene.len.human), ]) > 1, 1]), na.rm = T),
    xaxt = "n", xlim = c(1, 7), main = "human bulk (GTEx) data", xlab = "", ylim = c(0, 0.8), col = "grey",
    lwd = 2)
lines(density(log10(gene.len.dog[deseq_gtex_upgenes, 1]), na.rm = T), col = cs_gg[1], lwd = 2)
lines(density(log10(gene.len.dog[deseq_gtex_downgenes, 1]), na.rm = T), col = cs_gg[2], lwd = 2)
legend("topleft", lty = 1, col = c("grey", cs_gg[1:2]), lwd = 3, legend = c("expressed genes", "up during aging",
    "down during aging"), bty = "n")
add_axis(1)
plot(density(log10(gene.len.human[rowMeans(gtex.tpm[row.names(gene.len.human), ]) > 1, 1]), na.rm = T),
    xaxt = "n", xlim = c(1, 7), main = "human snRNA (excitatory neuron) data", xlab = "", ylim = c(0,
        1), col = "grey", lwd = 2)
lines(density(log10(gene.len.dog[degs.up.human, 1]), na.rm = T), col = cs_gg[1], lwd = 2)
lines(density(log10(gene.len.dog[degs.down.human, 1]), na.rm = T), col = cs_gg[2], lwd = 2)
legend("topleft", lty = 1, col = c("grey", cs_gg[1:2]), lwd = 3, legend = c("expressed genes", "up during aging",
    "down during aging"), bty = "n")
add_axis(1)
dev.off()

boxplot(log10(rowMeans(dog_aging.exp[deseq_upgenes2, ]) + 0.1), log10(rowMeans(dog_aging.exp[deseq_downgenes2,
    ]) + 0.1), log10(rowMeans(gtex.tpm[deseq_gtex_upgenes, ]) + 0.1), log10(rowMeans(gtex.tpm[deseq_gtex_downgenes,
    ]) + 0.1), log10(rowMeans(gtex.tpm[degs.up.human, ]) + 0.1), log10(rowMeans(gtex.tpm[degs.up.human,
    ]) + 0.1), staplewex = 0, outline = F, lty = 1)

# enrichment analysis----
# read data
require(MutationalPatterns)
enrich.mut <- read.table("tables/enrichment/merged.sMutation.all.overlap.reducedCol.tab", header = T,
    row.names = 1)
meta.mut <- read.table("tables/enrichment/mutation.metadata.tab", header = F, row.names = 1)
enrich.ctrl <- read.table("tables/enrichment/control.sMutation.all.overlap.reducedCol.tab", header = T,
    row.names = 1)
meta.ctrl <- read.table("tables/enrichment/control.metadata.tab", header = F, row.names = 1)

t <- read.table("tables/annotation/canFam3.merged.elements.bed", header = F, row.names = NULL)
inter_n <- t[t[, 5] == "intergenic", 4]
gene_n <- t[t[, 5] == "mRNA.gene", 4]

dog.exp.extN <- read.table("tables/rnaseq/dog.snRNA.ext_neuron.exp.txt", header = T, row.names = 1)
dog.exp.extN <- apply(dog.exp.extN, 1, mean)

tcs <- c(rgb(41, 52, 115, maxColorValue = 255), rgb(58, 106, 161, maxColorValue = 255), rgb(84, 158,
    200, maxColorValue = 255), rgb(164, 205, 227, maxColorValue = 255), rgb(199, 57, 54, maxColorValue = 255))
# read expression data the dog aging paper
dog_aging.exp.avg <- apply(dog_aging.exp[gene_n, ], 1, mean)
dog_aging.exp.avg.young <- apply(dog_aging.exp[gene_n, 1:5], 1, mean)

dog.exp.extN <- dog.exp.extN[intersect(gene_n, names(dog.exp.extN))]
#list_tpm=chunk(names(dog_aging.exp.avg[order(dog_aging.exp.avg)]),5)
list_tpm <- chunk(names(dog.exp.extN[order(dog.exp.extN)]), 5)
list_len <- chunk(row.names(gene.len.dog)[order(gene.len.dog[, 1])], 5)

pdf("figures/boxplot_gene_exp_aging_genes_dog_bulk.pdf", width = 3, height = 3, useDingbats = F)
par(mar = c(4, 4, 1, 1), tcl = 0.3, bty = "n")
boxplot(log10(dog_aging.exp.avg[deseq_downgenes2] + 1), log10(dog_aging.exp.avg[deseq_upgenes2] + 1),
    log10(dog_aging.exp.avg + 1), pch = 20, outline = F, staplewex = 0, lty = 1, ylim = c(0, 4), yaxt = "n",
    ylab = "TPM in PFC", xaxt = "n")
add_axis(2)
axis(1, 1:3, label = c("Aging\nDown", "Aging\nUp", "All\nGenes"), lwd = 0)
dev.off()

boxplot(log10(dog_aging.exp.avg[degs.down.dog] + 1), log10(dog_aging.exp.avg[degs.up.dog] + 1), log10(dog_aging.exp.avg +
    1))

boxplot(log10(gtex.exp[deseq_gtex_downgenes, "Brain - Frontal Cortex (BA9)"] + 1), log10(gtex.exp[deseq_gtex_upgenes,
    "Brain - Frontal Cortex (BA9)"] + 1), log10(gtex.exp[, "Brain - Frontal Cortex (BA9)"] + 1))

boxplot(log10(gtex.exp[degs.down.human, "Brain - Frontal Cortex (BA9)"] + 1), log10(gtex.exp[degs.up.human,
    "Brain - Frontal Cortex (BA9)"] + 1), log10(gtex.exp[, "Brain - Frontal Cortex (BA9)"] + 1))


#list_tpm=chunk(names(dog_aging.exp.avg[order(dog_aging.exp.avg)]),5)
tl <- list()
for (i in 1:5) {
    tl[[i]] <- log10(dog_aging.exp.avg[list_tpm[[i]]] + 1)
}
tl_len <- list()
for (i in 1:5) {
    tl_len[[i]] <- log10(gene.len.dog[list_tpm[[i]], 1] + 1)
}
tl_gc <- list()
for (i in 1:5) {
    tl_gc[[i]] <- gene.len.dog[list_tpm[[i]], "GC"]
}
# tl_tmp=tl
# tl_tmp[['down_in_aging']]=log10(dog_aging.exp.avg[deseq_downgenes2]+1)
# tl_tmp[['up_in_aging']]=log10(dog_aging.exp.avg[deseq_upgenes2]+1)
pdf("figures/gene_sets_definition.pdf", width = 6, height = 5, useDingbats = F)
par(mfrow = c(1, 3), mar = c(10, 4, 3, 1), tcl = 0.3, bty = "n", cex = 5/6)
boxplot(tl, col = c(tcs), ylab = "log10(TPM+1)", xlab = "", pch = 20, lty = 1, staplewex = 0, outline = F,
    xaxt = "n")
axis(1, 3, "gene sets by expression\nlow >>> high", lwd = 0)
axis(1, c(7:8), c("down in aging", "up in aging"), las = 3, lwd = 0)
boxplot(tl_len, col = c(tcs), ylab = "gene length", xlab = "", pch = 20, lty = 1, staplewex = 0, outline = F,
    xaxt = "n", yaxt = "n")
add_axis(2)
axis(1, 3, "gene sets by expression\nlow >>> high", lwd = 0)
axis(1, c(7:8), c("down in aging", "up in aging"), las = 3, lwd = 0)
boxplot(tl_gc, col = c(tcs), ylab = "GC content", xlab = "", pch = 20, lty = 1, staplewex = 0, outline = F,
    xaxt = "n")
axis(1, 3, "gene sets by expression\nlow >>> high", lwd = 0)
axis(1, c(7:8), c("down in aging", "up in aging"), las = 3, lwd = 0)
dev.off()

list_elements <- list(mRNA.gene.1st_quantile = "mRNA.gene", mRNA.gene.2nd_quantile = "mRNA.gene", mRNA.gene.3rd_quantile = "mRNA.gene",
    mRNA.gene.4th_quantile = "mRNA.gene", mRNA.gene.5th_quantile = "mRNA.gene", mRNA.gene = "mRNA.gene",
    intergenic = "intergenic", aging_down.gene = "mRNA.gene", aging_up.gene = "mRNA.gene", aging_down.snRNA = "mRNA.gene",
    aging_up.snRNA = "mRNA.gene", mRNA.gene.len_1 = "mRNA.gene", mRNA.gene.len_2 = "mRNA.gene", mRNA.gene.len_3 = "mRNA.gene",
    mRNA.gene.len_4 = "mRNA.gene", mRNA.gene.len_5 = "mRNA.gene")
list_genes <- list(mRNA.gene.1st_quantile = list_tpm[[1]], mRNA.gene.2nd_quantile = list_tpm[[2]], mRNA.gene.3rd_quantile = list_tpm[[3]],
    mRNA.gene.4th_quantile = list_tpm[[4]], mRNA.gene.5th_quantile = list_tpm[[5]], mRNA.gene = gene_n,
    intergenic = inter_n, aging_down.gene = deseq_downgenes2, aging_up.gene = deseq_upgenes2, aging_down.snRNA = degs.down.dog,
    aging_up.snRNA = degs.up.dog, mRNA.gene.len_1 = list_len[[1]], mRNA.gene.len_2 = list_len[[2]], mRNA.gene.len_3 = list_len[[3]],
    mRNA.gene.len_4 = list_len[[4]], mRNA.gene.len_5 = list_len[[5]])

# function----
sd <- 1
fun_enrich <- function(enrich_mut, enrich_ctrl, list_elements, list_genes, column_for_elements, column_for_genes) {
    ce <- column_for_elements
    cg <- column_for_genes
    vector_out <- rep(NA, length(list_elements))
    names(vector_out) <- names(list_elements)
    vector_out_upper <- vector_out
    vector_out_lower <- vector_out
    for (n in names(vector_out)) {
        tn_m <- which(enrich_mut[, ce] %in% list_elements[[n]] & enrich_mut[, cg] %in% list_genes[[n]])
        tn_c <- which(enrich_ctrl[, ce] %in% list_elements[[n]] & enrich_ctrl[, cg] %in% list_genes[[n]])
        vector_out[n] <- (length(tn_m) + sd)/((length(tn_c))/dim(enrich_ctrl)[1] * dim(enrich_mut)[1] +
            sd)
        vector_out_lower[n] <- (pois.exact(length(tn_m))[1, 4] + sd)/((length(tn_c))/dim(enrich_ctrl)[1] *
            dim(enrich_mut)[1] + sd)
        vector_out_upper[n] <- (pois.exact(length(tn_m))[1, 5] + sd)/((length(tn_c))/dim(enrich_ctrl)[1] *
            dim(enrich_mut)[1] + sd)
    }
    return(rbind(vector_out, vector_out_lower, vector_out_upper))
}
fun_enrich_pv_num <- function(enrich_mut, enrich_ctrl, list_elements, list_genes, column_for_elements,
    column_for_genes) {
    ce <- column_for_elements
    cg <- column_for_genes
    vector_out1 <- rep(NA, length(list_elements))
    names(vector_out1) <- names(list_elements)
    vector_out2 <- rep(NA, length(list_elements))
    names(vector_out2) <- names(list_elements)
    vector_out3 <- rep(NA, length(list_elements))
    names(vector_out3) <- names(list_elements)
    for (n in names(vector_out1)) {
        tn_m <- which(enrich_mut[, ce] %in% list_elements[[n]] & enrich_mut[, cg] %in% list_genes[[n]])
        tn_c <- which(enrich_ctrl[, ce] %in% list_elements[[n]] & enrich_ctrl[, cg] %in% list_genes[[n]])
        n1 <- length(tn_m)
        n2 <- (length(tn_c))/dim(enrich_ctrl)[1] * dim(enrich_mut)[1]
        vector_out1[n] <- if (n2 > 0)
            poisson.test(n1, T = 1, r = n2)$p.value else NA_real_
        vector_out2[n] <- n1
        vector_out3[n] <- n2
    }
    return(rbind(vector_out1, vector_out2, vector_out3))
}
# when decomposite mutations into signatures, also correct for trinucleotide composition
fun_mut_profile <- function(enrich_mut, enrich_ctrl, meta_mut, list_elements, list_genes, column_for_elements,
    column_for_genes, nmf.sig) {
    ce <- column_for_elements
    cg <- column_for_genes
    mat_out <- matrix(NA, nrow(nmf.sig), length(list_elements))
    colnames(mat_out) <- names(list_elements)
    row.names(mat_out) <- row.names(nmf.sig)
    for (n in colnames(mat_out)) {
        tn_m <- which(enrich_mut[, ce] %in% list_elements[[n]] & enrich_mut[, cg] %in% list_genes[[n]])
        mut_profile <- data.frame(context1 = rep(0, dim(nmf.sig)[1]), context_all = rep(0, dim(nmf.sig)[1]))
        row.names(mut_profile) <- row.names(nmf.sig)
        tmp1 <- table(meta_mut[tn_m, 3])
        tmp2 <- table(meta_mut[, 3])
        mut_profile[intersect(row.names(mut_profile), names(tmp1)), 1] <- tmp1[intersect(row.names(mut_profile),
            names(tmp1))]
        mut_profile[intersect(row.names(mut_profile), names(tmp2)), 2] <- tmp2[intersect(row.names(mut_profile),
            names(tmp2))]
        if (dim(mat_out)[1] == 96) {
            mat_out[, n] <- mut_profile[, 1]/tri.factors[, 1]
        } else {
            mat_out[, n] <- mut_profile[, 1]
        }
    }
    return(mat_out)
}
fun_enrich_strand <- function(enrich_mut, list_elements, list_genes, column_for_elements, column_for_genes,
    column_for_strand) {
    ce <- column_for_elements
    cg <- column_for_genes
    cstr <- column_for_strand
    vector_out <- rep(NA, length(list_elements))
    names(vector_out) <- names(list_elements)
    vector_out2 <- rep(NA, length(list_elements))
    names(vector_out2) <- names(list_elements)
    for (n in names(vector_out)) {
        tn_w <- which(enrich_mut[, ce] %in% list_elements[[n]] & enrich_mut[, cg] %in% list_genes[[n]] &
            enrich_mut[, cstr] == "+")
        tn_c <- which(enrich_mut[, ce] %in% list_elements[[n]] & enrich_mut[, cg] %in% list_genes[[n]] &
            enrich_mut[, cstr] == "-")
        nw <- length(tn_w)
        nc <- length(tn_c)
        vector_out[n] <- (nw + 1)/(nw + nc + 2)
        vector_out2[n] <- if (nw + nc > 0)
            binom.test(nw, nw + nc, p = 0.5)$p.value else NA_real_
    }
    return(cbind(vector_out, vector_out2))
}
fun_num_sig <- function(enrich_mut, enrich_ctrl, meta_mut, list_elements, list_genes, column_for_elements,
    column_for_genes, nmf.sig) {
    ce <- column_for_elements
    cg <- column_for_genes
    mat_out <- matrix(NA, ncol(nmf.sig) * 2, length(list_elements))
    colnames(mat_out) <- names(list_elements)
    row.names(mat_out) <- c(colnames(nmf.sig), unlist(lapply(colnames(nmf.sig), paste, "expected")))
    for (n in colnames(mat_out)) {
        tn_m <- which(enrich_mut[, ce] %in% list_elements[[n]] & enrich_mut[, cg] %in% list_genes[[n]])
        tn_c <- which(enrich_ctrl[, ce] %in% list_elements[[n]] & enrich_ctrl[, cg] %in% list_genes[[n]])
        mut_profile <- data.frame(context1 = rep(0, dim(nmf.sig)[1]), context_all = rep(0, dim(nmf.sig)[1]))
        row.names(mut_profile) <- row.names(nmf.sig)
        tmp1 <- table(meta_mut[tn_m, 3])
        tmp2 <- table(meta_mut[, 3])
        mut_profile[intersect(row.names(mut_profile), names(tmp1)), 1] <- tmp1[intersect(row.names(mut_profile),
            names(tmp1))]
        mut_profile[intersect(row.names(mut_profile), names(tmp2)), 2] <- tmp2[intersect(row.names(mut_profile),
            names(tmp2))]
        if (length(mut_profile) == 96) {
            mut_profile/tri.factors[, 1]
        }
        tmp_res <- fit_to_signatures(mut_profile, as.matrix(nmf.sig))
        tmp <- tmp_res$contribution * colSums(nmf.sig)
        row.names(tmp) <- colnames(nmf.sig)
        for (s in row.names(tmp)) {
            mat_out[s, n] <- tmp[s, 1]
            mat_out[paste(s, "expected"), n] <- (length(tn_c) + sd)/dim(enrich_ctrl)[1] * tmp[s, 2]
        }
    }
    return(mat_out)
}
fun_enrich_sig <- function(enrich_mut, enrich_ctrl, meta_mut, list_elements, list_genes, column_for_elements,
    column_for_genes, nmf.sig, correct = T) {
    ce <- column_for_elements
    cg <- column_for_genes
    mat_out <- matrix(NA, ncol(nmf.sig) * 3, length(list_elements))
    colnames(mat_out) <- names(list_elements)
    row.names(mat_out) <- c(colnames(nmf.sig), paste(colnames(nmf.sig), "lower", sep = "_"), paste(colnames(nmf.sig),
        "upper", sep = "_"))
    for (n in colnames(mat_out)) {
        tn_m <- which(enrich_mut[, ce] %in% list_elements[[n]] & enrich_mut[, cg] %in% list_genes[[n]])
        tn_c <- which(enrich_ctrl[, ce] %in% list_elements[[n]] & enrich_ctrl[, cg] %in% list_genes[[n]])
        mut_profile <- data.frame(context1 = rep(0, dim(nmf.sig)[1]), context_all = rep(0, dim(nmf.sig)[1]))
        row.names(mut_profile) <- row.names(nmf.sig)
        tmp1 <- table(meta_mut[tn_m, 3])
        tmp2 <- table(meta_mut[, 3])
        mut_profile[intersect(row.names(mut_profile), names(tmp1)), 1] <- tmp1[intersect(row.names(mut_profile),
            names(tmp1))]
        mut_profile[intersect(row.names(mut_profile), names(tmp2)), 2] <- tmp2[intersect(row.names(mut_profile),
            names(tmp2))]
        if (length(mut_profile) == 96 & correct) {
            mut_profile/tri.factors[, 1]
        }
        tmp_res <- fit_to_signatures(mut_profile, as.matrix(nmf.sig))
        tmp <- tmp_res$contribution * colSums(nmf.sig)
        row.names(tmp) <- colnames(nmf.sig)
        for (s in row.names(tmp)) {
            mat_out[s, n] <- (tmp[s, 1] + sd)/(length(tn_c) + sd) * dim(enrich_ctrl)[1]/tmp[s, 2]
            mat_out[paste(s, "lower", sep = "_"), n] <- (pois.exact(tmp[s, 1])[1, 4] + sd)/(length(tn_c) +
                sd) * dim(enrich_ctrl)[1]/tmp[s, 2]
            mat_out[paste(s, "upper", sep = "_"), n] <- (pois.exact(tmp[s, 1])[1, 5] + sd)/(length(tn_c) +
                sd) * dim(enrich_ctrl)[1]/tmp[s, 2]
        }
    }
    return(mat_out)
}
fun_enrich_sig_pv <- function(enrich_mut, enrich_ctrl, meta_mut, list_elements, list_genes, column_for_elements,
    column_for_genes, nmf.sig, correct = T) {
    ce <- column_for_elements
    cg <- column_for_genes
    mat_out <- matrix(NA, ncol(nmf.sig), length(list_elements))
    colnames(mat_out) <- names(list_elements)
    row.names(mat_out) <- colnames(nmf.sig)
    for (n in colnames(mat_out)) {
        tn_m <- which(enrich_mut[, ce] %in% list_elements[[n]] & enrich_mut[, cg] %in% list_genes[[n]])
        tn_c <- which(enrich_ctrl[, ce] %in% list_elements[[n]] & enrich_ctrl[, cg] %in% list_genes[[n]])
        mut_profile <- data.frame(context1 = rep(0, dim(nmf.sig)[1]), context_all = rep(0, dim(nmf.sig)[1]))
        row.names(mut_profile) <- row.names(nmf.sig)
        tmp1 <- table(meta_mut[tn_m, 3])
        tmp2 <- table(meta_mut[, 3])
        mut_profile[intersect(row.names(mut_profile), names(tmp1)), 1] <- tmp1[intersect(row.names(mut_profile),
            names(tmp1))]
        mut_profile[intersect(row.names(mut_profile), names(tmp2)), 2] <- tmp2[intersect(row.names(mut_profile),
            names(tmp2))]
        if (length(mut_profile) == 96 & correct) {
            mut_profile/tri.factors[, 1]
        }
        tmp_res <- fit_to_signatures(mut_profile, as.matrix(nmf.sig))
        tmp <- tmp_res$contribution * colSums(nmf.sig)
        row.names(tmp) <- colnames(nmf.sig)
        for (s in row.names(tmp)) {
            n1 <- tmp[s, 1]
            n2 <- (length(tn_c))/dim(enrich_ctrl)[1] * tmp[s, 2]
            mat_out[s, n] <- if (n2 > 0)
                poisson.test(n1, T = 1, r = n2)$p.value else NA_real_
        }
    }
    return(mat_out)
}
plot_enrichment_with_ci <- function(mat_enrich, mat_enrich_pv = NULL, main = "sSNV", ylim = c(0, 2),
    tcs = NULL, ylab = "Observed / Expected", star_y = NULL, cex.names = 1, cex.star = 0.6, las = 3) {

    ord <- c("mRNA.gene", "intergenic", "mRNA.gene.1st_quantile", "mRNA.gene.2nd_quantile", "mRNA.gene.3rd_quantile",
        "mRNA.gene.4th_quantile", "mRNA.gene.5th_quantile", "aging_down.gene", "aging_up.gene", "aging_down.snRNA",
        "aging_up.snRNA", "mRNA.gene.len_1", "mRNA.gene.len_2", "mRNA.gene.len_3", "mRNA.gene.len_4",
        "mRNA.gene.len_5")

    ord <- ord[ord %in% colnames(mat_enrich)]

    vals <- as.numeric(mat_enrich[1, ord])
    lower <- as.numeric(mat_enrich[2, ord])
    upper <- as.numeric(mat_enrich[3, ord])

    labs <- c(mRNA.gene = "genic", intergenic = "intergenic", mRNA.gene.1st_quantile = "1st quintile",
        mRNA.gene.2nd_quantile = "2nd quintile", mRNA.gene.3rd_quantile = "3rd quintile", mRNA.gene.4th_quantile = "4th quintile",
        mRNA.gene.5th_quantile = "5th quintile", aging_down.gene = "down; bulk", aging_up.gene = "up; bulk",
        aging_down.snRNA = "down; snRNA", aging_up.snRNA = "up; snRNA", mRNA.gene.len_1 = "length 1st quintile",
        mRNA.gene.len_2 = "length 2nd quintile", mRNA.gene.len_3 = "length 3rd quintile", mRNA.gene.len_4 = "length 4th quintile",
        mRNA.gene.len_5 = "length 5th quintile")[ord]

    if (is.null(tcs)) {
        tcs <- colorRampPalette(c("#313695", "#74add1", "#d73027"))(5)
    }

    cols <- c(mRNA.gene = "black", intergenic = "grey", mRNA.gene.1st_quantile = tcs[1], mRNA.gene.2nd_quantile = tcs[2],
        mRNA.gene.3rd_quantile = tcs[3], mRNA.gene.4th_quantile = tcs[4], mRNA.gene.5th_quantile = tcs[5],
        aging_down.gene = "black", aging_up.gene = "black", aging_down.snRNA = "black", aging_up.snRNA = "black",
        mRNA.gene.len_1 = tcs[1], mRNA.gene.len_2 = tcs[2], mRNA.gene.len_3 = tcs[3], mRNA.gene.len_4 = tcs[4],
        mRNA.gene.len_5 = tcs[5])[ord]

    space <- c(mRNA.gene = 0.5, intergenic = 0.5, mRNA.gene.1st_quantile = 2.5, mRNA.gene.2nd_quantile = 0.5,
        mRNA.gene.3rd_quantile = 0.5, mRNA.gene.4th_quantile = 0.5, mRNA.gene.5th_quantile = 0.5, aging_down.gene = 2.5,
        aging_up.gene = 0.5, aging_down.snRNA = 0.5, aging_up.snRNA = 0.5, mRNA.gene.len_1 = 2.5, mRNA.gene.len_2 = 0.5,
        mRNA.gene.len_3 = 0.5, mRNA.gene.len_4 = 0.5, mRNA.gene.len_5 = 0.5)[ord]

    bp <- barplot(vals, names.arg = labs, space = space, border = FALSE, col = cols, las = las, main = main,
        ylim = ylim, ylab = ylab, cex.names = cex.names)

    abline(h = 1, lty = 2)

    arrows(x0 = bp, y0 = lower, x1 = bp, y1 = upper, angle = 90, code = 3, length = 0.05, lwd = 1.2,
        col = "darkgrey")

    if (!is.null(mat_enrich_pv)) {
        pvals <- as.numeric(mat_enrich_pv[ord])

        # stars <- rep('', length(pvals))
        # stars[pvals < 0.05] <- '*'
        # stars[pvals < 0.01] <- '**'
        # stars[pvals < 0.001] <- '***'
        stars <- format(pvals, scientific = T, digit = 3)

        if (is.null(star_y)) {
            star_y <- ylim[2] * 0.85
        }

        text(x = bp, y = rep(star_y, length(bp)), labels = stars, cex = cex.star, srt = 90)
    }

    invisible(bp)
}

# on genomic elements; non-seperated by signatures----
mat_enrich <- fun_enrich(enrich.mut[meta.mut[, 2] == "SNV", ], enrich.ctrl, list_elements, list_genes,
    1, 7)
mat_enrich_pv <- fun_enrich_pv_num(enrich.mut[meta.mut[, 2] == "SNV", ], enrich.ctrl, list_elements,
    list_genes, 1, 7)
for (i in 1) {
    mat_enrich_pv[i, c(1:7, 10, 11)] <- p.adjust(mat_enrich_pv[i, c(1:7, 10, 11)], method = "fdr")
}
mat_enrich2 <- fun_enrich(enrich.mut[meta.mut[, 2] != "SNV", ], enrich.ctrl, list_elements, list_genes,
    1, 7)
mat_enrich_pv2 <- fun_enrich_pv_num(enrich.mut[meta.mut[, 2] != "SNV", ], enrich.ctrl, list_elements,
    list_genes, 1, 7)
for (i in 1) {
    mat_enrich_pv2[i, c(1:7, 10, 11)] <- p.adjust(mat_enrich_pv2[i, c(1:7, 10, 11)], method = "fdr")
}
pdf("figures/enrichment_analysis/enrichment_genes.dog.pdf", width = 5, height = 4, useDingbats = F)
par(mar = c(10, 4, 3, 1), tcl = 0.3, bty = "n", cex = 5/6)
plot_enrichment_with_ci(mat_enrich, mat_enrich_pv[1, ], main = "sSNV", ylim = c(0, 2))
plot_enrichment_with_ci(mat_enrich2, mat_enrich_pv2[1, ], main = "sindel", ylim = c(0, 3))
dev.off()

# on genomic elements----
mat_enrich_sig <- fun_enrich_sig(enrich.mut[meta.mut[, 2] == "SNV", ], enrich.ctrl, meta.mut[meta.mut[,
    2] == "SNV", ], list_elements, list_genes, 1, 7, nmf_res_snv$signatures)
mat_enrich_sig_pv <- fun_enrich_sig_pv(enrich.mut[meta.mut[, 2] == "SNV", ], enrich.ctrl, meta.mut[meta.mut[,
    2] == "SNV", ], list_elements, list_genes, 1, 7, nmf_res_snv$signatures)
mat_enrich_sig_num <- fun_num_sig(enrich.mut[meta.mut[, 2] == "SNV", ], enrich.ctrl, meta.mut[meta.mut[,
    2] == "SNV", ], list_elements, list_genes, 1, 7, nmf_res_snv$signatures)
mat_profile <- fun_mut_profile(enrich.mut[meta.mut[, 2] == "SNV", ], enrich.ctrl, meta.mut[meta.mut[,
    2] == "SNV", ], list_elements, list_genes, 1, 7, nmf_res_snv$signatures)
pdf("figures/enrichment_analysis/profiles.snv_dog.pdf", width = 12, height = 7.5/2, useDingbats = F)
par(tcl = 0.3, bty = "n", cex = 5/6)
for (n in colnames(mat_profile)) {
    barplot(mat_profile[, n], col = col_sig, border = "white", space = 0, names = row.names(mat_profile),
        las = 3, main = n, ylab = "N of identified sSNV")
    abline(v = c(1:5) * 16)
    text(c(1:6) * 16 - 8, max(mat_profile[, n]) * 19/20, label = c("C>A", "C>G", "C>T", "T>A", "T>C",
        "T>G"))
}
dev.off()
#for(i in 1:2){mat_enrich_sig_pv[i,]=p.adjust(mat_enrich_sig_pv[i,],method='fdr')}
mat_enrich_sig2 <- fun_enrich_sig(enrich.mut[meta.mut[, 2] != "SNV", ], enrich.ctrl, meta.mut[meta.mut[,
    2] != "SNV", ], list_elements, list_genes, 1, 7, nmf_res_indel$signatures)
mat_enrich_sig_pv2 <- fun_enrich_sig_pv(enrich.mut[meta.mut[, 2] != "SNV", ], enrich.ctrl, meta.mut[meta.mut[,
    2] != "SNV", ], list_elements, list_genes, 1, 7, nmf_res_indel$signatures)
mat_enrich_sig_num2 <- fun_num_sig(enrich.mut[meta.mut[, 2] != "SNV", ], enrich.ctrl, meta.mut[meta.mut[,
    2] != "SNV", ], list_elements, list_genes, 1, 7, nmf_res_indel$signatures)
mat_profile2 <- fun_mut_profile(enrich.mut[meta.mut[, 2] != "SNV", ], enrich.ctrl, meta.mut[meta.mut[,
    2] != "SNV", ], list_elements, list_genes, 1, 7, nmf_res_indel$signatures)
pdf("figures/enrichment_analysis/profiles.indel_dog.pdf", width = 12, height = 7.5/2, useDingbats = F)
par(tcl = 0.3, bty = "n", cex = 5/6)
for (n in colnames(mat_profile2)) {
    ns <- c(rep(c(1, 2, 3, 4, 5, "6+"), 2), rep(c(0, 1, 2, 3, 4, "5+"), 2), rep(c(1, 2, 3, 4, 5, "6+"),
        4), rep(c(0, 1, 2, 3, 4, "5+"), 4), c(1, 1, 2, 1, 2, 3, 1, 2, 3, 4, "5+"))
    barplot(mat_profile2[, n], col = col_sig_indel, border = "white", space = 0, names = ns, las = 3,
        main = n, ylab = "N of identified sSNV")
    axis(1, c(6.5, 18.5, 36.5, 60.5, 77.5), label = c("homopolymer\nlength", "homopolymer\nlength", "number of\nrepeat units",
        "number of\nrepeat units", "microhomology\nlength"), lwd = 0, padj = 2, cex.axis = 0.5)
    axis(3, c(6.5, 18.5, 36.5, 60.5, 77.5), label = c("1bp deletion", "1bp insertion", ">1 bp deletion at repeats",
        ">1 bp insertion at repeats", "microhomology"), lwd = 0, padj = 2, cex.axis = 0.5)
    text(c(3.5, 9.5, 15.5, 21.5, 27.5, 33.5, 39.5, 45.5, 51.5, 57.5, 63.5, 69.5, 73, 74.5, 77, 81), max(mat_profile2[,
        n]) * 0.95, label = c(rep(c("C", "T"), 2), rep(c(2, 3, 4, "5+"), 2), 2, 3, 4, "5+"), cex = 0.5,
        col = unique(tcs))
}
dev.off()
#for(i in 1:2){mat_enrich_sig_pv2[i,]=p.adjust(mat_enrich_sig_pv2[i,],method='fdr')}
pdf("figures/enrichment_analysis/enrichment_genes_signature.pdf", width = 4, height = 4, useDingbats = F)
par(mar = c(10, 4, 3, 1), tcl = 0.3, bty = "n", cex = 5/6)
plot_enrichment_with_ci(mat_enrich_sig[c(2, 4, 6), ], mat_enrich_sig_pv[2, ], main = "signature A1",
    ylim = c(0, 3))
plot_enrichment_with_ci(mat_enrich_sig[c(1, 3, 5), ], mat_enrich_sig_pv[1, ], main = "signature A2",
    ylim = c(0, 8))
plot_enrichment_with_ci(mat_enrich_sig2[c(1, 3, 5), ], mat_enrich_sig_pv2[1, ], main = "signature ID-A",
    ylim = c(0, 20))
plot_enrichment_with_ci(mat_enrich_sig2[c(2, 4, 6), ], mat_enrich_sig_pv2[2, ], main = "signature ID-B",
    ylim = c(0, 20))
dev.off()

# enrichment of strand bias----
list_elements <- list(mRNA.gene.1st_quantile = "mRNA.gene", mRNA.gene.2nd_quantile = "mRNA.gene", mRNA.gene.3rd_quantile = "mRNA.gene",
    mRNA.gene.4th_quantile = "mRNA.gene", mRNA.gene.5th_quantile = "mRNA.gene")
list_genes <- list(mRNA.gene.1st_quantile = list_tpm[[1]], mRNA.gene.2nd_quantile = list_tpm[[2]], mRNA.gene.3rd_quantile = list_tpm[[3]],
    mRNA.gene.4th_quantile = list_tpm[[4]], mRNA.gene.5th_quantile = list_tpm[[5]])
tmat <- matrix(NA, 6, 5)
row.names(tmat) <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
colnames(tmat) <- names(list_elements)
tmat2 <- matrix(NA, 6, 5)
row.names(tmat2) <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
colnames(tmat2) <- names(list_elements)
for (i in row.names(tmat)) {
    tn <- which(substring(meta.mut[, 3], 3, 5) == i)
    tmp <- fun_enrich_strand(enrich.mut[tn, ], list_elements, list_genes, 1, 7, 4)
    tmat[i, ] <- tmp[, 1]
    tmat2[i, ] <- tmp[, 2]
}

tv1 <- c()
for (i in 1:5) {
    tv1 <- c(tv1, tmat[, i])
}
tv2 <- c()
for (i in 1:5) {
    tv2 <- c(tv2, tmat2[, i])
}
tv3 <- rep("", 30)
tv3[p.adjust(tv2, method = "fdr") < 0.05] <- "*"
at <- c(1:6, 8:13, 15:20, 22:27, 29:34)
pdf("figures/enrichment_analysis/strand_analysis.dog.pdf", width = 3.5, height = 6, useDingbats = F)
par(mar = c(4, 1, 1, 1), tcl = 0.3, bty = "n", cex = 5/6)
plot(100 - tv1 * 100, at, pch = 19, yaxt = "n", xlab = "%sSNV in transcribed strand", xlim = c(20, 80),
    col = c("#80b1d3", "black", "#e41a1c", "#d9d9d9", "#b3de69", "#fccde5"))
abline(h = c(7, 14, 21, 28))
text(20, at, label = tv3, cex = 3)
text(20, 1:6, label = row.names(tmat), pos = 4)
dev.off()

# regression model to see if correlation between T>C (A1 mutation), transcription, gene length, and GC----
tc_mut <- row.names(meta.mut)[substring(meta.mut[, 3], 3, 5) == "T>C"]
t <- table(enrich.mut[tc_mut[enrich.mut[tc_mut, 1] == "mRNA.gene"], ]$name1)
gene.len.dog[, "ct_mut"] <- 0
tn <- intersect(names(t), row.names(gene.len.dog))
gene.len.dog[tn, "ct_mut"] <- t[tn]
gene.len.dog[, "expression"] <- NA
# tn=intersect(row.names(gene.len.dog),names(dog_aging.exp.avg))
# gene.len.dog[tn,'expression']=dog_aging.exp.avg[tn]
# now use exicitatory neuron expression instead of PFC bulk
tn <- intersect(row.names(gene.len.dog), names(dog.exp.extN))
gene.len.dog[tn, "expression"] <- dog.exp.extN[tn]
fit1 <- glm.nb(ct_mut ~ expression + GC + offset(log(gene_length)), data = gene.len.dog)
fit2 <- glm.nb(ct_mut ~ expression + GC + log10(gene_length) + offset(log(gene_length)), data = gene.len.dog)
summary(fit2)
tdf <- summary(fit2)$coefficients
row.names(tdf) <- c("Intercept", "ln(Transcription)", "GC", "log10GeneLength")
write.table(tdf, "tables/For_manuscript/model_negativeBinomial_TC_mutation_in_genes_dog.txt", sep = "\t",
    quote = F)

tc_mut <- row.names(meta.mut)[meta.mut[, 2] == "SNV"]  # also use all sSNVs
t <- table(enrich.mut[tc_mut[enrich.mut[tc_mut, 1] == "mRNA.gene"], ]$name1)
gene.len.dog[, "ct_mut"] <- 0
tn <- intersect(names(t), row.names(gene.len.dog))
gene.len.dog[tn, "ct_mut"] <- t[tn]
gene.len.dog[, "expression"] <- NA
# tn=intersect(row.names(gene.len.dog),names(dog_aging.exp.avg))
# gene.len.dog[tn,'expression']=dog_aging.exp.avg[tn]
# now use exicitatory neuron expression instead of PFC bulk
tn <- intersect(row.names(gene.len.dog), names(dog.exp.extN))
gene.len.dog[tn, "expression"] <- dog.exp.extN[tn]
fit1 <- glm.nb(ct_mut ~ expression + GC + offset(log(gene_length)), data = gene.len.dog)
fit2 <- glm.nb(ct_mut ~ expression + GC + log10(gene_length) + offset(log(gene_length)), data = gene.len.dog)
summary(fit2)
tdf <- summary(fit2)$coefficients
row.names(tdf) <- c("Intercept", "ln(Transcription)", "GC", "log10GeneLength")
write.table(tdf, "tables/For_manuscript/model_negativeBinomial_all_SNV_in_genes_dog.txt", sep = "\t",
    quote = F)

tc_mut <- row.names(meta.mut)[meta.mut[, 2] != "SNV"]  # also use all sSNVs
t <- table(enrich.mut[tc_mut[enrich.mut[tc_mut, 1] == "mRNA.gene"], ]$name1)
gene.len.dog[, "ct_mut"] <- 0
tn <- intersect(names(t), row.names(gene.len.dog))
gene.len.dog[tn, "ct_mut"] <- t[tn]
gene.len.dog[, "expression"] <- NA
# tn=intersect(row.names(gene.len.dog),names(dog_aging.exp.avg))
# gene.len.dog[tn,'expression']=dog_aging.exp.avg[tn]
# now use exicitatory neuron expression instead of PFC bulk
tn <- intersect(row.names(gene.len.dog), names(dog.exp.extN))
gene.len.dog[tn, "expression"] <- dog.exp.extN[tn]
fit1 <- glm.nb(ct_mut ~ expression + GC + offset(log(gene_length)), data = gene.len.dog)
fit2 <- glm.nb(ct_mut ~ expression + GC + log10(gene_length) + offset(log(gene_length)), data = gene.len.dog)
summary(fit2)
tdf <- summary(fit2)$coefficients
row.names(tdf) <- c("Intercept", "ln(Transcription)", "GC", "log10GeneLength")
write.table(tdf, "tables/For_manuscript/model_negativeBinomial_all_sindel_in_genes_dog.txt", sep = "\t",
    quote = F)

# enrichment analysis for human----
# read data
require(MutationalPatterns)
enrich.mut <- read.table("../tables/enrichment/merged.sMutation.all.overlap.reducedCol.tab", header = T,
    row.names = 1)
meta.mut <- read.table("../tables/enrichment/mutation.metadata.tab", header = F, row.names = 1)
enrich.ctrl <- read.table("../tables/enrichment/control.sMutation.all.overlap.reducedCol.tab", header = T,
    row.names = 1)
meta.ctrl <- read.table("../tables/enrichment/control.metadata.tab", header = F, row.names = 1)

neuron.cpm <- read.table("../aging_project/tables/rnaseq/neurons.CPM.mRNA.txt", header = F, row.names = 1,
    check.names = F)
colnames(neuron.cpm) <- c("all_neuron", "ext", "inb")
neuron.cpm <- neuron.cpm[order(neuron.cpm[, 2]), ]
inter_n <- as.vector(read.table("../aging_project/tables/annotation/hs37.intergenic.bed", header = F,
    row.names = NULL)[, 4])
list_tpm <- chunk(row.names(neuron.cpm), 5)
gene_n <- row.names(neuron.cpm)
tn <- intersect(row.names(gtex.exp), gene_n)
list_tpm_crypt <- chunk(tn[order(gtex.exp[tn, "Colon - Sigmoid"])], 5)
list_tpm_bulk <- chunk(tn[order(gtex.exp[tn, "Brain - Frontal Cortex (BA9)"])], 5)

tcs <- c(rgb(41, 52, 115, maxColorValue = 255), rgb(58, 106, 161, maxColorValue = 255), rgb(84, 158,
    200, maxColorValue = 255), rgb(164, 205, 227, maxColorValue = 255), rgb(199, 57, 54, maxColorValue = 255))

colnames(neuron.cpm) <- c("all_neuron", "ext", "inb")
gene.len.human <- gene.len.human[order(gene.len.human[, 1]), ]
list_len <- chunk(row.names(gene.len.human), 5)
boxplot(neuron.cpm[list_len[[1]], 2], neuron.cpm[list_len[[2]], 2], neuron.cpm[list_len[[3]], 2], neuron.cpm[list_len[[4]],
    2], neuron.cpm[list_len[[5]], 2])

tl <- list()
for (i in 1:5) {
    tl[[i]] <- neuron.cpm[list_tpm[[i]], 2]
}
tl_len <- list()
for (i in 1:5) {
    tl_len[[i]] <- log10(gene.len.human[list_tpm[[i]], 1] + 1)
}
tl_gc <- list()
for (i in 1:5) {
    tl_gc[[i]] <- gene.len.human[list_tpm[[i]], "GC"]
}
pdf("figures/gene_sets_definition_human.pdf", width = 6, height = 5, useDingbats = F)
par(mfrow = c(1, 3), mar = c(10, 4, 3, 1), tcl = 0.3, bty = "n", cex = 5/6)
boxplot(tl, col = c(tcs), ylab = "log(CPM+1)", xlab = "", pch = 20, lty = 1, staplewex = 0, outline = F,
    xaxt = "n")
axis(1, 3, "gene sets by expression\nlow >>> high", lwd = 0)
axis(1, c(7:8), c("down in aging", "up in aging"), las = 3, lwd = 0)
boxplot(tl_len, col = c(tcs), ylab = "gene length", xlab = "", pch = 20, lty = 1, staplewex = 0, outline = F,
    xaxt = "n", yaxt = "n")
add_axis(2)
axis(1, 3, "gene sets by expression\nlow >>> high", lwd = 0)
axis(1, c(7:8), c("down in aging", "up in aging"), las = 3, lwd = 0)
boxplot(tl_gc, col = c(tcs), ylab = "GC content", xlab = "", pch = 20, lty = 1, staplewex = 0, outline = F,
    xaxt = "n")
axis(1, 3, "gene sets by expression\nlow >>> high", lwd = 0)
axis(1, c(7:8), c("down in aging", "up in aging"), las = 3, lwd = 0)
dev.off()
list_elements <- list(mRNA.gene.1st_quantile = "mRNA.gene", mRNA.gene.2nd_quantile = "mRNA.gene", mRNA.gene.3rd_quantile = "mRNA.gene",
    mRNA.gene.4th_quantile = "mRNA.gene", mRNA.gene.5th_quantile = "mRNA.gene", mRNA.gene = "mRNA.gene",
    intergenic = "intergenic", aging_down.gene = "mRNA.gene", aging_up.gene = "mRNA.gene", aging_down.snRNA = "mRNA.gene",
    aging_up.snRNA = "mRNA.gene", mRNA.gene.len_1 = "mRNA.gene", mRNA.gene.len_2 = "mRNA.gene", mRNA.gene.len_3 = "mRNA.gene",
    mRNA.gene.len_4 = "mRNA.gene", mRNA.gene.len_5 = "mRNA.gene")
list_genes <- list(mRNA.gene.1st_quantile = list_tpm[[1]], mRNA.gene.2nd_quantile = list_tpm[[2]], mRNA.gene.3rd_quantile = list_tpm[[3]],
    mRNA.gene.4th_quantile = list_tpm[[4]], mRNA.gene.5th_quantile = list_tpm[[5]], mRNA.gene = gene_n,
    intergenic = inter_n, aging_down.gene = deseq_gtex_downgenes, aging_up.gene = deseq_gtex_upgenes,
    aging_down.snRNA = degs.down.human, aging_up.snRNA = degs.up.human, mRNA.gene.len_1 = list_len[[1]],
    mRNA.gene.len_2 = list_len[[2]], mRNA.gene.len_3 = list_len[[3]], mRNA.gene.len_4 = list_len[[4]],
    mRNA.gene.len_5 = list_len[[5]])

# function----
sd <- 1
fun_enrich <- function(enrich_mut, enrich_ctrl, list_elements, list_genes, column_for_elements, column_for_genes) {
    ce <- column_for_elements
    cg <- column_for_genes
    vector_out <- rep(NA, length(list_elements))
    names(vector_out) <- names(list_elements)
    vector_out_upper <- vector_out
    vector_out_lower <- vector_out
    for (n in names(vector_out)) {
        tn_m <- which(enrich_mut[, ce] %in% list_elements[[n]] & enrich_mut[, cg] %in% list_genes[[n]])
        tn_c <- which(enrich_ctrl[, ce] %in% list_elements[[n]] & enrich_ctrl[, cg] %in% list_genes[[n]])
        vector_out[n] <- (length(tn_m) + sd)/((length(tn_c))/dim(enrich_ctrl)[1] * dim(enrich_mut)[1] +
            sd)
        vector_out_lower[n] <- (pois.exact(length(tn_m))[1, 4] + sd)/((length(tn_c))/dim(enrich_ctrl)[1] *
            dim(enrich_mut)[1] + sd)
        vector_out_upper[n] <- (pois.exact(length(tn_m))[1, 5] + sd)/((length(tn_c))/dim(enrich_ctrl)[1] *
            dim(enrich_mut)[1] + sd)
    }
    return(rbind(vector_out, vector_out_lower, vector_out_upper))
}
fun_enrich_pv_num <- function(enrich_mut, enrich_ctrl, list_elements, list_genes, column_for_elements,
    column_for_genes) {
    ce <- column_for_elements
    cg <- column_for_genes
    vector_out1 <- rep(NA, length(list_elements))
    names(vector_out1) <- names(list_elements)
    vector_out2 <- rep(NA, length(list_elements))
    names(vector_out2) <- names(list_elements)
    vector_out3 <- rep(NA, length(list_elements))
    names(vector_out3) <- names(list_elements)
    for (n in names(vector_out1)) {
        tn_m <- which(enrich_mut[, ce] %in% list_elements[[n]] & enrich_mut[, cg] %in% list_genes[[n]])
        tn_c <- which(enrich_ctrl[, ce] %in% list_elements[[n]] & enrich_ctrl[, cg] %in% list_genes[[n]])
        n1 <- length(tn_m)
        n2 <- (length(tn_c))/dim(enrich_ctrl)[1] * dim(enrich_mut)[1]
        vector_out1[n] <- if (n2 > 0)
            poisson.test(n1, T = 1, r = n2)$p.value else NA_real_
        vector_out2[n] <- n1
        vector_out3[n] <- n2
    }
    return(rbind(vector_out1, vector_out2, vector_out3))
}
# for human, doesn't need to correct for trinucleotide composition
fun_mut_profile <- function(enrich_mut, enrich_ctrl, meta_mut, list_elements, list_genes, column_for_elements,
    column_for_genes, nmf.sig) {
    ce <- column_for_elements
    cg <- column_for_genes
    mat_out <- matrix(NA, nrow(nmf.sig), length(list_elements))
    colnames(mat_out) <- names(list_elements)
    row.names(mat_out) <- row.names(nmf.sig)
    for (n in colnames(mat_out)) {
        tn_m <- which(enrich_mut[, ce] %in% list_elements[[n]] & enrich_mut[, cg] %in% list_genes[[n]])
        mut_profile <- data.frame(context1 = rep(0, dim(nmf.sig)[1]), context_all = rep(0, dim(nmf.sig)[1]))
        row.names(mut_profile) <- row.names(nmf.sig)
        tmp1 <- table(meta_mut[tn_m, 3])
        tmp2 <- table(meta_mut[, 3])
        mut_profile[intersect(row.names(mut_profile), names(tmp1)), 1] <- tmp1[intersect(row.names(mut_profile),
            names(tmp1))]
        mut_profile[intersect(row.names(mut_profile), names(tmp2)), 2] <- tmp2[intersect(row.names(mut_profile),
            names(tmp2))]
        mat_out[, n] <- mut_profile[, 1]
    }
    return(mat_out)
}
fun_enrich_strand <- function(enrich_mut, list_elements, list_genes, column_for_elements, column_for_genes,
    column_for_strand) {
    ce <- column_for_elements
    cg <- column_for_genes
    cstr <- column_for_strand
    vector_out <- rep(NA, length(list_elements))
    names(vector_out) <- names(list_elements)
    vector_out2 <- rep(NA, length(list_elements))
    names(vector_out2) <- names(list_elements)
    for (n in names(vector_out)) {
        tn_w <- which(enrich_mut[, ce] %in% list_elements[[n]] & enrich_mut[, cg] %in% list_genes[[n]] &
            enrich_mut[, cstr] == "+")
        tn_c <- which(enrich_mut[, ce] %in% list_elements[[n]] & enrich_mut[, cg] %in% list_genes[[n]] &
            enrich_mut[, cstr] == "-")
        nw <- length(tn_w)
        nc <- length(tn_c)
        vector_out[n] <- (nw + 1)/(nw + nc + 2)
        vector_out2[n] <- if (nw + nc > 0)
            binom.test(nw, nw + nc, p = 0.5)$p.value else NA_real_
    }
    return(cbind(vector_out, vector_out2))
}
fun_num_sig <- function(enrich_mut, enrich_ctrl, meta_mut, list_elements, list_genes, column_for_elements,
    column_for_genes, nmf.sig) {
    ce <- column_for_elements
    cg <- column_for_genes
    mat_out <- matrix(NA, ncol(nmf.sig) * 2, length(list_elements))
    colnames(mat_out) <- names(list_elements)
    row.names(mat_out) <- c(colnames(nmf.sig), unlist(lapply(colnames(nmf.sig), paste, "expected")))
    for (n in colnames(mat_out)) {
        tn_m <- which(enrich_mut[, ce] %in% list_elements[[n]] & enrich_mut[, cg] %in% list_genes[[n]])
        tn_c <- which(enrich_ctrl[, ce] %in% list_elements[[n]] & enrich_ctrl[, cg] %in% list_genes[[n]])
        mut_profile <- data.frame(context1 = rep(0, dim(nmf.sig)[1]), context_all = rep(0, dim(nmf.sig)[1]))
        row.names(mut_profile) <- row.names(nmf.sig)
        tmp1 <- table(meta_mut[tn_m, 3])
        tmp2 <- table(meta_mut[, 3])
        mut_profile[intersect(row.names(mut_profile), names(tmp1)), 1] <- tmp1[intersect(row.names(mut_profile),
            names(tmp1))]
        mut_profile[intersect(row.names(mut_profile), names(tmp2)), 2] <- tmp2[intersect(row.names(mut_profile),
            names(tmp2))]
        tmp_res <- fit_to_signatures(mut_profile, as.matrix(nmf.sig))
        tmp <- tmp_res$contribution * colSums(nmf.sig)
        row.names(tmp) <- colnames(nmf.sig)
        for (s in row.names(tmp)) {
            mat_out[s, n] <- tmp[s, 1]
            mat_out[paste(s, "expected"), n] <- (length(tn_c) + sd)/dim(enrich_ctrl)[1] * tmp[s, 2]
        }
    }
    return(mat_out)
}

# on genomic elements; non-seperated by signatures----
mat_enrich <- fun_enrich(enrich.mut[meta.mut[, 2] == "SNV" & meta.mut[, 1] %in% cells_normal, ], enrich.ctrl,
    list_elements, list_genes, 1, 7)
mat_enrich_pv <- fun_enrich_pv_num(enrich.mut[meta.mut[, 2] == "SNV" & meta.mut[, 1] %in% cells_normal,
    ], enrich.ctrl, list_elements, list_genes, 1, 7)
for (i in 1) {
    mat_enrich_pv[i, c(1:7, 10, 11)] <- p.adjust(mat_enrich_pv[i, c(1:7, 10, 11)], method = "fdr")
}
mat_enrich2 <- fun_enrich(enrich.mut[meta.mut[, 2] != "SNV" & meta.mut[, 1] %in% cells_normal, ], enrich.ctrl,
    list_elements, list_genes, 1, 7)
mat_enrich_pv2 <- fun_enrich_pv_num(enrich.mut[meta.mut[, 2] != "SNV" & meta.mut[, 1] %in% cells_normal,
    ], enrich.ctrl, list_elements, list_genes, 1, 7)
for (i in 1) {
    mat_enrich_pv2[i, c(1:7, 10, 11)] <- p.adjust(mat_enrich_pv2[i, c(1:7, 10, 11)], method = "fdr")
}
pdf("figures/enrichment_analysis/enrichment_genes_human.pdf", width = 5, height = 4, useDingbats = F)
par(mar = c(10, 4, 3, 1), tcl = 0.3, bty = "n", cex = 5/6)
plot_enrichment_with_ci(mat_enrich, mat_enrich_pv[1, ], main = "sSNV", ylim = c(0, 2))
plot_enrichment_with_ci(mat_enrich2, mat_enrich_pv2[1, ], main = "sindel", ylim = c(0, 3))
dev.off()

list_genes <- list(mRNA.gene.1st_quantile = list_tpm_bulk[[1]], mRNA.gene.2nd_quantile = list_tpm_bulk[[2]],
    mRNA.gene.3rd_quantile = list_tpm_bulk[[3]], mRNA.gene.4th_quantile = list_tpm_bulk[[4]], mRNA.gene.5th_quantile = list_tpm_bulk[[5]],
    mRNA.gene = gene_n, intergenic = inter_n, aging_down.gene = deseq_gtex_downgenes, aging_up.gene = deseq_gtex_upgenes,
    aging_down.snRNA = degs.down.human, aging_up.snRNA = degs.up.human, mRNA.gene.len_1 = list_len[[1]],
    mRNA.gene.len_2 = list_len[[2]], mRNA.gene.len_3 = list_len[[3]], mRNA.gene.len_4 = list_len[[4]],
    mRNA.gene.len_5 = list_len[[5]])
mat_enrich <- fun_enrich(enrich.mut[meta.mut[, 2] == "SNV" & meta.mut[, 1] %in% cells_normal, ], enrich.ctrl,
    list_elements, list_genes, 1, 7)
mat_enrich_pv <- fun_enrich_pv_num(enrich.mut[meta.mut[, 2] == "SNV" & meta.mut[, 1] %in% cells_normal,
    ], enrich.ctrl, list_elements, list_genes, 1, 7)
#for(i in 1){mat_enrich_pv[i,]=p.adjust(mat_enrich_pv[i,],method='fdr')}
mat_enrich2 <- fun_enrich(enrich.mut[meta.mut[, 2] != "SNV" & meta.mut[, 1] %in% cells_normal, ], enrich.ctrl,
    list_elements, list_genes, 1, 7)
mat_enrich_pv2 <- fun_enrich_pv_num(enrich.mut[meta.mut[, 2] != "SNV" & meta.mut[, 1] %in% cells_normal,
    ], enrich.ctrl, list_elements, list_genes, 1, 7)
#for(i in 1:2){mat_enrich_sig_pv2[i,]=p.adjust(mat_enrich_sig_pv2[i,],method='fdr')}
pdf("figures/enrichment_analysis/enrichment_genes_human_expQuintile_by_bulk.pdf", width = 5, height = 4,
    useDingbats = F)
par(mar = c(10, 4, 3, 1), tcl = 0.3, bty = "n", cex = 5/6)
plot_enrichment_with_ci(mat_enrich, mat_enrich_pv, main = "sSNV", ylim = c(0, 2))
plot_enrichment_with_ci(mat_enrich2, mat_enrich_pv2, main = "sindel", ylim = c(0, 3))
dev.off()

# on genomic elements----
list_genes <- list(mRNA.gene.1st_quantile = list_tpm[[1]], mRNA.gene.2nd_quantile = list_tpm[[2]], mRNA.gene.3rd_quantile = list_tpm[[3]],
    mRNA.gene.4th_quantile = list_tpm[[4]], mRNA.gene.5th_quantile = list_tpm[[5]], mRNA.gene = gene_n,
    intergenic = inter_n, aging_down.gene = deseq_gtex_downgenes, aging_up.gene = deseq_gtex_upgenes,
    aging_down.snRNA = degs.down.human, aging_up.snRNA = degs.up.human, mRNA.gene.len_1 = list_len[[1]],
    mRNA.gene.len_2 = list_len[[2]], mRNA.gene.len_3 = list_len[[3]], mRNA.gene.len_4 = list_len[[4]],
    mRNA.gene.len_5 = list_len[[5]])

mat_enrich_sig <- fun_enrich_sig(enrich.mut[meta.mut[, 2] == "SNV" & meta.mut[, 1] %in% cells_normal,
    ], enrich.ctrl, meta.mut[meta.mut[, 2] == "SNV" & meta.mut[, 1] %in% cells_normal, ], list_elements,
    list_genes, 1, 7, nmf_res_snv$signatures, correct = F)
mat_enrich_sig_pv <- fun_enrich_sig_pv(enrich.mut[meta.mut[, 2] == "SNV" & meta.mut[, 1] %in% cells_normal,
    ], enrich.ctrl, meta.mut[meta.mut[, 2] == "SNV" & meta.mut[, 1] %in% cells_normal, ], list_elements,
    list_genes, 1, 7, nmf_res_snv$signatures, correct = F)
mat_enrich_sig_num <- fun_num_sig(enrich.mut[meta.mut[, 2] == "SNV" & meta.mut[, 1] %in% cells_normal,
    ], enrich.ctrl, meta.mut[meta.mut[, 2] == "SNV" & meta.mut[, 1] %in% cells_normal, ], list_elements,
    list_genes, 1, 7, nmf_res_snv$signatures)
mat_profile <- fun_mut_profile(enrich.mut[meta.mut[, 2] == "SNV" & meta.mut[, 1] %in% cells_normal, ],
    enrich.ctrl, meta.mut[meta.mut[, 2] == "SNV" & meta.mut[, 1] %in% cells_normal, ], list_elements,
    list_genes, 1, 7, nmf_res_snv$signatures)
pdf("figures/enrichment_analysis/profiles.snv_human.pdf", width = 12, height = 7.5/2, useDingbats = F)
par(tcl = 0.3, bty = "n", cex = 5/6)
for (n in colnames(mat_profile)) {
    barplot(mat_profile[, n], col = col_sig, border = "white", space = 0, names = row.names(mat_profile),
        las = 3, main = n, ylab = "N of identified sSNV")
    abline(v = c(1:5) * 16)
    text(c(1:6) * 16 - 8, max(mat_profile[, n]) * 19/20, label = c("C>A", "C>G", "C>T", "T>A", "T>C",
        "T>G"))
}
dev.off()
#for(i in 1:2){mat_enrich_sig_pv[i,]=p.adjust(mat_enrich_sig_pv[i,],method='fdr')}
mat_enrich_sig2 <- fun_enrich_sig(enrich.mut[meta.mut[, 2] != "SNV" & meta.mut[, 1] %in% cells_normal,
    ], enrich.ctrl, meta.mut[meta.mut[, 2] != "SNV" & meta.mut[, 1] %in% cells_normal, ], list_elements,
    list_genes, 1, 7, nmf_res_indel$signatures, correct = F)
mat_enrich_sig_pv2 <- fun_enrich_sig_pv(enrich.mut[meta.mut[, 2] != "SNV" & meta.mut[, 1] %in% cells_normal,
    ], enrich.ctrl, meta.mut[meta.mut[, 2] != "SNV" & meta.mut[, 1] %in% cells_normal, ], list_elements,
    list_genes, 1, 7, nmf_res_indel$signatures, correct = F)
mat_enrich_sig_num2 <- fun_num_sig(enrich.mut[meta.mut[, 2] != "SNV" & meta.mut[, 1] %in% cells_normal,
    ], enrich.ctrl, meta.mut[meta.mut[, 2] != "SNV" & meta.mut[, 1] %in% cells_normal, ], list_elements,
    list_genes, 1, 7, nmf_res_indel$signatures)
mat_profile2 <- fun_mut_profile(enrich.mut[meta.mut[, 2] != "SNV" & meta.mut[, 1] %in% cells_normal,
    ], enrich.ctrl, meta.mut[meta.mut[, 2] != "SNV" & meta.mut[, 1] %in% cells_normal, ], list_elements,
    list_genes, 1, 7, nmf_res_indel$signatures)
pdf("figures/enrichment_analysis/profiles.indel_human.pdf", width = 12, height = 7.5/2, useDingbats = F)
par(tcl = 0.3, bty = "n", cex = 5/6)
for (n in colnames(mat_profile2)) {
    ns <- c(rep(c(1, 2, 3, 4, 5, "6+"), 2), rep(c(0, 1, 2, 3, 4, "5+"), 2), rep(c(1, 2, 3, 4, 5, "6+"),
        4), rep(c(0, 1, 2, 3, 4, "5+"), 4), c(1, 1, 2, 1, 2, 3, 1, 2, 3, 4, "5+"))
    barplot(mat_profile2[, n], col = col_sig_indel, border = "white", space = 0, names = ns, las = 3,
        main = n, ylab = "N of identified sSNV")
    axis(1, c(6.5, 18.5, 36.5, 60.5, 77.5), label = c("homopolymer\nlength", "homopolymer\nlength", "number of\nrepeat units",
        "number of\nrepeat units", "microhomology\nlength"), lwd = 0, padj = 2, cex.axis = 0.5)
    axis(3, c(6.5, 18.5, 36.5, 60.5, 77.5), label = c("1bp deletion", "1bp insertion", ">1 bp deletion at repeats",
        ">1 bp insertion at repeats", "microhomology"), lwd = 0, padj = 2, cex.axis = 0.5)
    text(c(3.5, 9.5, 15.5, 21.5, 27.5, 33.5, 39.5, 45.5, 51.5, 57.5, 63.5, 69.5, 73, 74.5, 77, 81), max(mat_profile2[,
        n]) * 0.95, label = c(rep(c("C", "T"), 2), rep(c(2, 3, 4, "5+"), 2), 2, 3, 4, "5+"), cex = 0.5,
        col = unique(tcs))
}
dev.off()
#for(i in 1:2){mat_enrich_sig_pv2[i,]=p.adjust(mat_enrich_sig_pv2[i,],method='fdr')}
pdf("figures/enrichment_analysis/enrichment_genes_signature_human.pdf", width = 4, height = 4, useDingbats = F)
par(mar = c(10, 4, 3, 1), tcl = 0.3, bty = "n", cex = 5/6)
plot_enrichment_with_ci(mat_enrich_sig[c(2, 4, 6), ], mat_enrich_sig_pv[2, ], main = "singature A1",
    ylim = c(0, 3))
plot_enrichment_with_ci(mat_enrich_sig[c(1, 3, 5), ], mat_enrich_sig_pv[1, ], main = "singature A2",
    ylim = c(0, 3))
plot_enrichment_with_ci(mat_enrich_sig2[c(1, 3, 5), ], mat_enrich_sig_pv2[1, ], main = "singature ID-A",
    ylim = c(0, 3))
plot_enrichment_with_ci(mat_enrich_sig2[c(2, 4, 6), ], mat_enrich_sig_pv2[2, ], main = "singature ID-B",
    ylim = c(0, 10))
dev.off()

list_genes <- list(mRNA.gene.1st_quantile = list_tpm_bulk[[1]], mRNA.gene.2nd_quantile = list_tpm_bulk[[2]],
    mRNA.gene.3rd_quantile = list_tpm_bulk[[3]], mRNA.gene.4th_quantile = list_tpm_bulk[[4]], mRNA.gene.5th_quantile = list_tpm_bulk[[5]],
    mRNA.gene = gene_n, intergenic = inter_n, aging_down.gene = deseq_gtex_downgenes, aging_up.gene = deseq_gtex_upgenes,
    aging_down.snRNA = degs.down.human, aging_up.snRNA = degs.up.human, mRNA.gene.len_1 = list_len[[1]],
    mRNA.gene.len_2 = list_len[[2]], mRNA.gene.len_3 = list_len[[3]], mRNA.gene.len_4 = list_len[[4]],
    mRNA.gene.len_5 = list_len[[5]])
mat_enrich_sig <- fun_enrich_sig(enrich.mut[meta.mut[, 2] == "SNV" & meta.mut[, 1] %in% cells_normal,
    ], enrich.ctrl, meta.mut[meta.mut[, 2] == "SNV" & meta.mut[, 1] %in% cells_normal, ], list_elements,
    list_genes, 1, 7, nmf_res_snv$signatures, correct = F)
mat_enrich_sig_pv <- fun_enrich_sig_pv(enrich.mut[meta.mut[, 2] == "SNV" & meta.mut[, 1] %in% cells_normal,
    ], enrich.ctrl, meta.mut[meta.mut[, 2] == "SNV" & meta.mut[, 1] %in% cells_normal, ], list_elements,
    list_genes, 1, 7, nmf_res_snv$signatures, correct = F)
mat_enrich_sig2 <- fun_enrich_sig(enrich.mut[meta.mut[, 2] != "SNV" & meta.mut[, 1] %in% cells_normal,
    ], enrich.ctrl, meta.mut[meta.mut[, 2] != "SNV" & meta.mut[, 1] %in% cells_normal, ], list_elements,
    list_genes, 1, 7, nmf_res_indel$signatures, correct = F)
mat_enrich_sig_pv2 <- fun_enrich_sig_pv(enrich.mut[meta.mut[, 2] != "SNV" & meta.mut[, 1] %in% cells_normal,
    ], enrich.ctrl, meta.mut[meta.mut[, 2] != "SNV" & meta.mut[, 1] %in% cells_normal, ], list_elements,
    list_genes, 1, 7, nmf_res_indel$signatures, correct = F)
pdf("figures/enrichment_analysis/enrichment_genes_signature_human_expQuintile_by_bulk.pdf", width = 4,
    height = 4, useDingbats = F)
par(mar = c(10, 4, 3, 1), tcl = 0.3, bty = "n", cex = 5/6)
plot_enrichment_with_ci(mat_enrich_sig[c(2, 4, 6), ], mat_enrich_sig_pv[2, ], main = "singature A1",
    ylim = c(0, 3))
plot_enrichment_with_ci(mat_enrich_sig[c(1, 3, 5), ], mat_enrich_sig_pv[1, ], main = "singature A2",
    ylim = c(0, 3))
plot_enrichment_with_ci(mat_enrich_sig2[c(1, 3, 5), ], mat_enrich_sig_pv2[1, ], main = "singature ID-A",
    ylim = c(0, 3))
plot_enrichment_with_ci(mat_enrich_sig2[c(2, 4, 6), ], mat_enrich_sig_pv2[2, ], main = "singature ID-B",
    ylim = c(0, 10))
dev.off()

# enrichment of strand bias----
list_elements <- list(mRNA.gene.1st_quantile = "mRNA.gene", mRNA.gene.2nd_quantile = "mRNA.gene", mRNA.gene.3rd_quantile = "mRNA.gene",
    mRNA.gene.4th_quantile = "mRNA.gene", mRNA.gene.5th_quantile = "mRNA.gene")
list_genes <- list(mRNA.gene.1st_quantile = list_tpm[[1]], mRNA.gene.2nd_quantile = list_tpm[[2]], mRNA.gene.3rd_quantile = list_tpm[[3]],
    mRNA.gene.4th_quantile = list_tpm[[4]], mRNA.gene.5th_quantile = list_tpm[[5]])
tmat <- matrix(NA, 6, 5)
row.names(tmat) <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
colnames(tmat) <- names(list_elements)
tmat2 <- matrix(NA, 6, 5)
row.names(tmat2) <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
colnames(tmat2) <- names(list_elements)
for (i in row.names(tmat)) {
    tn <- which(substring(meta.mut[, 3], 3, 5) == i)
    tmp <- fun_enrich_strand(enrich.mut[tn, ], list_elements, list_genes, 1, 7, 4)
    tmat[i, ] <- tmp[, 1]
    tmat2[i, ] <- tmp[, 2]
}
tv1 <- c()
for (i in 1:5) {
    tv1 <- c(tv1, tmat[, i])
}
tv2 <- c()
for (i in 1:5) {
    tv2 <- c(tv2, tmat2[, i])
}
tv3 <- rep("", 30)
tv3[p.adjust(tv2, method = "fdr") < 0.05] <- "*"
at <- c(1:6, 8:13, 15:20, 22:27, 29:34)
pdf("figures/enrichment_analysis/strand_analysis_human.pdf", width = 3.5, height = 6, useDingbats = F)
par(mar = c(4, 1, 1, 1), tcl = 0.3, bty = "n", cex = 5/6)
plot(100 - tv1 * 100, at, pch = 19, yaxt = "n", xlab = "%sSNV in transcribed strand", xlim = c(30, 70),
    col = c("#80b1d3", "black", "#e41a1c", "#d9d9d9", "#b3de69", "#fccde5"))
abline(h = c(7, 14, 21, 28))
text(30, at, label = tv3, cex = 3)
text(30, 1:6, label = row.names(tmat), pos = 4)
dev.off()

# regression model to see if correlation between T>C (A1 mutation), transcription, gene length, and GC----
tc_mut <- row.names(meta.mut)[substring(meta.mut[, 3], 3, 5) == "T>C"]
t <- table(enrich.mut[tc_mut[enrich.mut[tc_mut, 1] == "mRNA.gene"], ]$name1)
gene.len.human[, "ct_mut"] <- 0
tn <- intersect(names(t), row.names(gene.len.human))
gene.len.human[tn, "ct_mut"] <- t[tn]
gene.len.human[, "expression"] <- 0
tn <- intersect(row.names(gene.len.human), row.names(neuron.cpm))
gene.len.human[tn, "expression"] <- neuron.cpm[tn, 2]
fit1 <- glm.nb(ct_mut ~ expression + GC + offset(log(gene_len)), data = gene.len.human)
fit2 <- glm.nb(ct_mut ~ expression + GC + log10(gene_len) + offset(log(gene_len)), data = gene.len.human)
summary(fit2)
tdf <- summary(fit2)$coefficients
row.names(tdf) <- c("Intercept", "ln(Transcription)", "GC", "log10GeneLength")
write.table(tdf, "tables/For_manuscript/model_negativeBinomial_TC_mutation_in_genes_human.txt", sep = "\t",
    quote = F)

tc_mut <- row.names(meta.mut)[meta.mut[, 2] == "SNV"]
t <- table(enrich.mut[tc_mut[enrich.mut[tc_mut, 1] == "mRNA.gene"], ]$name1)
gene.len.human[, "ct_mut"] <- 0
tn <- intersect(names(t), row.names(gene.len.human))
gene.len.human[tn, "ct_mut"] <- t[tn]
gene.len.human[, "expression"] <- 0
tn <- intersect(row.names(gene.len.human), row.names(neuron.cpm))
gene.len.human[tn, "expression"] <- neuron.cpm[tn, 2]
fit1 <- glm.nb(ct_mut ~ expression + GC + offset(log(gene_len)), data = gene.len.human)
fit2 <- glm.nb(ct_mut ~ expression + GC + log10(gene_len) + offset(log(gene_len)), data = gene.len.human)
summary(fit2)
tdf <- summary(fit2)$coefficients
row.names(tdf) <- c("Intercept", "ln(Transcription)", "GC", "log10GeneLength")
write.table(tdf, "tables/For_manuscript/model_negativeBinomial_all_sSNV_in_genes_human.txt", sep = "\t",
    quote = F)

tc_mut <- row.names(meta.mut)[meta.mut[, 2] != "SNV"]
t <- table(enrich.mut[tc_mut[enrich.mut[tc_mut, 1] == "mRNA.gene"], ]$name1)
gene.len.human[, "ct_mut"] <- 0
tn <- intersect(names(t), row.names(gene.len.human))
gene.len.human[tn, "ct_mut"] <- t[tn]
gene.len.human[, "expression"] <- 0
tn <- intersect(row.names(gene.len.human), row.names(neuron.cpm))
gene.len.human[tn, "expression"] <- neuron.cpm[tn, 2]
fit1 <- glm.nb(ct_mut ~ expression + GC + offset(log(gene_len)), data = gene.len.human)
fit2 <- glm.nb(ct_mut ~ expression + GC + log10(gene_len) + offset(log(gene_len)), data = gene.len.human)
summary(fit2)
tdf <- summary(fit2)$coefficients
row.names(tdf) <- c("Intercept", "ln(Transcription)", "GC", "log10GeneLength")
write.table(tdf, "tables/For_manuscript/model_negativeBinomial_all_sindel_in_genes_human.txt", sep = "\t",
    quote = F)

# enrichment analysis with crypt----
# read data
require(MutationalPatterns)
enrich.mut <- read.table("tables/enrichment/merged.sMutation.all.overlap.reducedCol.tab", header = T,
    row.names = 1)
meta.mut <- read.table("tables/enrichment/mutation.metadata.tab", header = F, row.names = 1)
enrich.ctrl <- read.table("tables/enrichment/control.sMutation.all.overlap.reducedCol.tab", header = T,
    row.names = 1)
meta.ctrl <- read.table("tables/enrichment/control.metadata.tab", header = F, row.names = 1)
enrich.mut.crypt <- read.table("tables/enrichment/cagan_dog.sMutation.all.overlap.reducedCol.tab", header = T,
    row.names = 1)
meta.mut.crypt <- read.table("tables/enrichment/cagan_dog.metadata.tab", header = F, row.names = 1)
t <- read.table("tables/annotation/canFam3.merged.elements.bed", header = F, row.names = NULL)
inter_n <- t[t[, 5] == "intergenic", 4]
gene_n <- t[t[, 5] == "mRNA.gene", 4]

dog.exp.colon <- read.table("tables/rnaseq/dog.salmon.tximport.tpm.txt", header = T, row.names = 1)
dog.exp.colon <- dog.exp.colon[, grep("Colon", colnames(dog.exp.colon))]
dog.id_name.map <- read.table("tables/annotation/dog.id_name.map", header = F, row.names = 1)
tn <- intersect(row.names(dog.id_name.map), row.names(dog.exp.colon))
tn <- tn[!duplicated(dog.id_name.map[tn, 1])]
dog.exp.colon <- dog.exp.colon[tn, ]
row.names(dog.exp.colon) <- dog.id_name.map[tn, 1]

trn <- intersect(row.names(dog_aging.exp), row.names(dog.exp.colon))
tdf <- cbind(apply(dog_aging.exp[trn, ], 1, mean), apply(dog.exp.colon[trn, ], 1, mean))
colnames(tdf) <- c("PFC", "Colon")
write.table(tdf, "tables/For_manuscript/dog.gene_exp.pfc_colon.tpm.txt", quote = F, sep = "\t")

tcs <- c(rgb(41, 52, 115, maxColorValue = 255), rgb(58, 106, 161, maxColorValue = 255), rgb(84, 158,
    200, maxColorValue = 255), rgb(164, 205, 227, maxColorValue = 255), rgb(199, 57, 54, maxColorValue = 255))
# read expression data the dog aging paper
dog_aging.exp.avg <- apply(dog_aging.exp[gene_n, ], 1, mean)
dog_aging.exp.avg.young <- apply(dog_aging.exp[gene_n, 1:5], 1, mean)
dog_aging.exp.avg.colon <- apply(dog.exp.colon[intersect(gene_n, row.names(dog.exp.colon)), ], 1, mean)
list_tpm <- chunk(names(dog_aging.exp.avg[order(dog_aging.exp.avg)]), 5)
list_len <- chunk(row.names(gene.len.dog)[order(gene.len.dog[, 1])], 5)
list_tpm_crypt <- chunk(names(dog_aging.exp.avg.colon[order(dog_aging.exp.avg.colon)]), 5)

# function----
sd <- 1
# when decomposite mutations into signatures, also correct for trinucleotide composition
fun_enrich <- function(enrich_mut, enrich_ctrl, list_elements, list_genes, column_for_elements, column_for_genes) {
    ce <- column_for_elements
    cg <- column_for_genes
    vector_out <- rep(NA, length(list_elements))
    names(vector_out) <- names(list_elements)
    vector_out_upper <- vector_out
    vector_out_lower <- vector_out
    for (n in names(vector_out)) {
        tn_m <- which(enrich_mut[, ce] %in% list_elements[[n]] & enrich_mut[, cg] %in% list_genes[[n]])
        tn_c <- which(enrich_ctrl[, ce] %in% list_elements[[n]] & enrich_ctrl[, cg] %in% list_genes[[n]])
        vector_out[n] <- (length(tn_m) + sd)/((length(tn_c))/dim(enrich_ctrl)[1] * dim(enrich_mut)[1] +
            sd)
        vector_out_lower[n] <- (pois.exact(length(tn_m))[1, 4] + sd)/((length(tn_c))/dim(enrich_ctrl)[1] *
            dim(enrich_mut)[1] + sd)
        vector_out_upper[n] <- (pois.exact(length(tn_m))[1, 5] + sd)/((length(tn_c))/dim(enrich_ctrl)[1] *
            dim(enrich_mut)[1] + sd)
    }
    return(rbind(vector_out, vector_out_lower, vector_out_upper))
}
fun_enrich_pv_num <- function(enrich_mut, enrich_ctrl, list_elements, list_genes, column_for_elements,
    column_for_genes) {
    ce <- column_for_elements
    cg <- column_for_genes
    vector_out1 <- rep(NA, length(list_elements))
    names(vector_out1) <- names(list_elements)
    vector_out2 <- rep(NA, length(list_elements))
    names(vector_out2) <- names(list_elements)
    vector_out3 <- rep(NA, length(list_elements))
    names(vector_out3) <- names(list_elements)
    for (n in names(vector_out1)) {
        tn_m <- which(enrich_mut[, ce] %in% list_elements[[n]] & enrich_mut[, cg] %in% list_genes[[n]])
        tn_c <- which(enrich_ctrl[, ce] %in% list_elements[[n]] & enrich_ctrl[, cg] %in% list_genes[[n]])
        n1 <- length(tn_m) + sd
        n2 <- (length(tn_c))/dim(enrich_ctrl)[1] * dim(enrich_mut)[1] + sd
        vector_out1[n] <- chisq.test(matrix(c(n1, n2, mean(n1, n2), mean(n1, n2))))$p.value
        vector_out2[n] <- n1
        vector_out3[n] <- n2
    }
    return(rbind(vector_out1, vector_out2, vector_out3))
}
fun_enrich_sig <- function(enrich_mut, enrich_ctrl, meta_mut, list_elements, list_genes, column_for_elements,
    column_for_genes, nmf.sig, correct = T) {
    ce <- column_for_elements
    cg <- column_for_genes
    mat_out <- matrix(NA, ncol(nmf.sig) * 3, length(list_elements))
    colnames(mat_out) <- names(list_elements)
    row.names(mat_out) <- c(colnames(nmf.sig), paste(colnames(nmf.sig), "lower", sep = "_"), paste(colnames(nmf.sig),
        "upper", sep = "_"))
    for (n in colnames(mat_out)) {
        tn_m <- which(enrich_mut[, ce] %in% list_elements[[n]] & enrich_mut[, cg] %in% list_genes[[n]])
        tn_c <- which(enrich_ctrl[, ce] %in% list_elements[[n]] & enrich_ctrl[, cg] %in% list_genes[[n]])
        mut_profile <- data.frame(context1 = rep(0, dim(nmf.sig)[1]), context_all = rep(0, dim(nmf.sig)[1]))
        row.names(mut_profile) <- row.names(nmf.sig)
        tmp1 <- table(meta_mut[tn_m, 3])
        tmp2 <- table(meta_mut[, 3])
        mut_profile[intersect(row.names(mut_profile), names(tmp1)), 1] <- tmp1[intersect(row.names(mut_profile),
            names(tmp1))]
        mut_profile[intersect(row.names(mut_profile), names(tmp2)), 2] <- tmp2[intersect(row.names(mut_profile),
            names(tmp2))]
        if (length(mut_profile) == 96 & correct) {
            mut_profile/tri.factors[, 1]
        }
        tmp_res <- fit_to_signatures(mut_profile, as.matrix(nmf.sig))
        tmp <- tmp_res$contribution * colSums(nmf.sig)
        row.names(tmp) <- colnames(nmf.sig)
        for (s in row.names(tmp)) {
            mat_out[s, n] <- (tmp[s, 1] + sd)/(length(tn_c) + sd) * dim(enrich_ctrl)[1]/tmp[s, 2]
            mat_out[paste(s, "lower", sep = "_"), n] <- (pois.exact(tmp[s, 1])[1, 4] + sd)/(length(tn_c) +
                sd) * dim(enrich_ctrl)[1]/tmp[s, 2]
            mat_out[paste(s, "upper", sep = "_"), n] <- (pois.exact(tmp[s, 1])[1, 5] + sd)/(length(tn_c) +
                sd) * dim(enrich_ctrl)[1]/tmp[s, 2]
        }
    }
    return(mat_out)
}
fun_enrich_sig_pv <- function(enrich_mut, enrich_ctrl, meta_mut, list_elements, list_genes, column_for_elements,
    column_for_genes, nmf.sig, correct = T) {
    ce <- column_for_elements
    cg <- column_for_genes
    mat_out <- matrix(NA, ncol(nmf.sig), length(list_elements))
    colnames(mat_out) <- names(list_elements)
    row.names(mat_out) <- colnames(nmf.sig)
    for (n in colnames(mat_out)) {
        tn_m <- which(enrich_mut[, ce] %in% list_elements[[n]] & enrich_mut[, cg] %in% list_genes[[n]])
        tn_c <- which(enrich_ctrl[, ce] %in% list_elements[[n]] & enrich_ctrl[, cg] %in% list_genes[[n]])
        mut_profile <- data.frame(context1 = rep(0, dim(nmf.sig)[1]), context_all = rep(0, dim(nmf.sig)[1]))
        row.names(mut_profile) <- row.names(nmf.sig)
        tmp1 <- table(meta_mut[tn_m, 3])
        tmp2 <- table(meta_mut[, 3])
        mut_profile[intersect(row.names(mut_profile), names(tmp1)), 1] <- tmp1[intersect(row.names(mut_profile),
            names(tmp1))]
        mut_profile[intersect(row.names(mut_profile), names(tmp2)), 2] <- tmp2[intersect(row.names(mut_profile),
            names(tmp2))]
        if (length(mut_profile) == 96 & correct) {
            mut_profile/tri.factors[, 1]
        }
        tmp_res <- fit_to_signatures(mut_profile, as.matrix(nmf.sig))
        tmp <- tmp_res$contribution * colSums(nmf.sig)
        row.names(tmp) <- colnames(nmf.sig)
        for (s in row.names(tmp)) {
            n1 <- tmp[s, 1] + sd
            n2 <- (length(tn_c) + sd)/dim(enrich_ctrl)[1] * tmp[s, 2]
            mat_out[s, n] <- chisq.test(matrix(c(n1, n2, mean(n1, n2), mean(n1, n2))))$p.value
        }
    }
    return(mat_out)
}
plot_enrichment_with_ci <- function(mat_enrich, mat_enrich_pv = NULL, main = "sSNV", ylim = c(0, 2),
    tcs = NULL, ylab = "Observed / Expected", star_y = NULL, cex.names = 1, cex.star = 0.6, las = 3) {

    ord <- c("mRNA.gene", "intergenic", "mRNA.gene.1st_quantile", "mRNA.gene.2nd_quantile", "mRNA.gene.3rd_quantile",
        "mRNA.gene.4th_quantile", "mRNA.gene.5th_quantile", "aging_down.gene", "aging_up.gene", "aging_down.snRNA",
        "aging_up.snRNA", "mRNA.gene.len_1", "mRNA.gene.len_2", "mRNA.gene.len_3", "mRNA.gene.len_4",
        "mRNA.gene.len_5")

    ord <- ord[ord %in% colnames(mat_enrich)]

    vals <- as.numeric(mat_enrich[1, ord])
    lower <- as.numeric(mat_enrich[2, ord])
    upper <- as.numeric(mat_enrich[3, ord])

    labs <- c(mRNA.gene = "genic", intergenic = "intergenic", mRNA.gene.1st_quantile = "1st quintile",
        mRNA.gene.2nd_quantile = "2nd quintile", mRNA.gene.3rd_quantile = "3rd quintile", mRNA.gene.4th_quantile = "4th quintile",
        mRNA.gene.5th_quantile = "5th quintile", aging_down.gene = "down; bulk", aging_up.gene = "up; bulk",
        aging_down.snRNA = "down; snRNA", aging_up.snRNA = "up; snRNA", mRNA.gene.len_1 = "length 1st quintile",
        mRNA.gene.len_2 = "length 2nd quintile", mRNA.gene.len_3 = "length 3rd quintile", mRNA.gene.len_4 = "length 4th quintile",
        mRNA.gene.len_5 = "length 5th quintile")[ord]

    if (is.null(tcs)) {
        tcs <- colorRampPalette(c("#313695", "#74add1", "#d73027"))(5)
    }

    cols <- c(mRNA.gene = "black", intergenic = "grey", mRNA.gene.1st_quantile = tcs[1], mRNA.gene.2nd_quantile = tcs[2],
        mRNA.gene.3rd_quantile = tcs[3], mRNA.gene.4th_quantile = tcs[4], mRNA.gene.5th_quantile = tcs[5],
        aging_down.gene = "black", aging_up.gene = "black", aging_down.snRNA = "black", aging_up.snRNA = "black",
        mRNA.gene.len_1 = tcs[1], mRNA.gene.len_2 = tcs[2], mRNA.gene.len_3 = tcs[3], mRNA.gene.len_4 = tcs[4],
        mRNA.gene.len_5 = tcs[5])[ord]

    space <- c(mRNA.gene = 0.5, intergenic = 0.5, mRNA.gene.1st_quantile = 2.5, mRNA.gene.2nd_quantile = 0.5,
        mRNA.gene.3rd_quantile = 0.5, mRNA.gene.4th_quantile = 0.5, mRNA.gene.5th_quantile = 0.5, aging_down.gene = 2.5,
        aging_up.gene = 0.5, aging_down.snRNA = 0.5, aging_up.snRNA = 0.5, mRNA.gene.len_1 = 2.5, mRNA.gene.len_2 = 0.5,
        mRNA.gene.len_3 = 0.5, mRNA.gene.len_4 = 0.5, mRNA.gene.len_5 = 0.5)[ord]

    bp <- barplot(vals, names.arg = labs, space = space, border = FALSE, col = cols, las = las, main = main,
        ylim = ylim, ylab = ylab, cex.names = cex.names)

    abline(h = 1, lty = 2)

    arrows(x0 = bp, y0 = lower, x1 = bp, y1 = upper, angle = 90, code = 3, length = 0.05, lwd = 1.2,
        col = "darkgrey")

    if (!is.null(mat_enrich_pv)) {
        pvals <- as.numeric(mat_enrich_pv[ord])

        # stars <- rep('', length(pvals))
        # stars[pvals < 0.05] <- '*'
        # stars[pvals < 0.01] <- '**'
        # stars[pvals < 0.001] <- '***'
        stars <- format(pvals, scientific = T, digit = 3)

        if (is.null(star_y)) {
            star_y <- ylim[2] * 0.85
        }

        text(x = bp, y = rep(star_y, length(bp)), labels = stars, cex = cex.star, srt = 90)
    }

    invisible(bp)
}

# on genomic elements----
list_elements <- list(mRNA.gene.1st_quantile = "mRNA.gene", mRNA.gene.2nd_quantile = "mRNA.gene", mRNA.gene.3rd_quantile = "mRNA.gene",
    mRNA.gene.4th_quantile = "mRNA.gene", mRNA.gene.5th_quantile = "mRNA.gene", mRNA.gene = "mRNA.gene",
    intergenic = "intergenic", aging_down.gene = "mRNA.gene", aging_up.gene = "mRNA.gene", aging_down.snRNA = "mRNA.gene",
    aging_up.snRNA = "mRNA.gene", mRNA.gene.len_1 = "mRNA.gene", mRNA.gene.len_2 = "mRNA.gene", mRNA.gene.len_3 = "mRNA.gene",
    mRNA.gene.len_4 = "mRNA.gene", mRNA.gene.len_5 = "mRNA.gene")
list_genes <- list(mRNA.gene.1st_quantile = list_tpm[[1]], mRNA.gene.2nd_quantile = list_tpm[[2]], mRNA.gene.3rd_quantile = list_tpm[[3]],
    mRNA.gene.4th_quantile = list_tpm[[4]], mRNA.gene.5th_quantile = list_tpm[[5]], mRNA.gene = gene_n,
    intergenic = inter_n, aging_down.gene = deseq_downgenes2, aging_up.gene = deseq_upgenes2, aging_down.snRNA = degs.down.dog,
    aging_up.snRNA = degs.up.dog, mRNA.gene.len_1 = list_len[[1]], mRNA.gene.len_2 = list_len[[2]], mRNA.gene.len_3 = list_len[[3]],
    mRNA.gene.len_4 = list_len[[4]], mRNA.gene.len_5 = list_len[[5]])

# neuron
mat_enrich_sig <- fun_enrich_sig(enrich.mut[meta.mut[, 2] == "SNV", ], enrich.ctrl, meta.mut[meta.mut[,
    2] == "SNV", ], list_elements, list_genes, 1, 7, nmf_res_snv_with_crypt$signatures)
mat_enrich_sig_pv <- fun_enrich_sig_pv(enrich.mut[meta.mut[, 2] == "SNV", ], enrich.ctrl, meta.mut[meta.mut[,
    2] == "SNV", ], list_elements, list_genes, 1, 7, nmf_res_snv_with_crypt$signatures)
mat_enrich_sig2 <- fun_enrich_sig(enrich.mut[meta.mut[, 2] != "SNV", ], enrich.ctrl, meta.mut[meta.mut[,
    2] != "SNV", ], list_elements, list_genes, 1, 7, nmf_res_indel_with_crypt$signatures)
mat_enrich_sig_pv2 <- fun_enrich_sig_pv(enrich.mut[meta.mut[, 2] != "SNV", ], enrich.ctrl, meta.mut[meta.mut[,
    2] != "SNV", ], list_elements, list_genes, 1, 7, nmf_res_indel_with_crypt$signatures)
pdf("figures/enrichment_analysis/with_crypt/dog.neuron.enrichment_genes_signature.pdf", width = 4, height = 4,
    useDingbats = F)
par(mar = c(10, 4, 3, 1), tcl = 0.3, bty = "n", cex = 5/6)
n_sigs <- nrow(mat_enrich_sig)/3
for (i in 1:n_sigs) {
    plot_enrichment_with_ci(mat_enrich_sig[c(i, i + n_sigs, i + n_sigs * 2), ], mat_enrich_sig_pv[i,
        ], main = row.names(mat_enrich_sig)[i], ylim = c(0, 3))
}
n_sigs <- nrow(mat_enrich_sig2)/3
for (i in 1:n_sigs) {
    plot_enrichment_with_ci(mat_enrich_sig2[c(i, i + n_sigs, i + n_sigs * 2), ], mat_enrich_sig_pv2[i,
        ], main = row.names(mat_enrich_sig2)[i], ylim = c(0, 10))
}
dev.off()

pdf("figures/enrichment_analysis/with_crypt/dog.neuron.enrichment_genes_signature_autoScaled.pdf", width = 4,
    height = 4, useDingbats = F)
par(mar = c(10, 4, 3, 1), tcl = 0.3, bty = "n", cex = 5/6)
n_sigs <- nrow(mat_enrich_sig)/3
for (i in 1:n_sigs) {
    plot_enrichment_with_ci(mat_enrich_sig[c(i, i + n_sigs, i + n_sigs * 2), ], mat_enrich_sig_pv[i,
        ], main = row.names(mat_enrich_sig)[i], ylim = c(0, max(mat_enrich_sig[c(i, i + n_sigs, i + n_sigs *
        2), ])))
}
n_sigs <- nrow(mat_enrich_sig2)/3
for (i in 1:n_sigs) {
    plot_enrichment_with_ci(mat_enrich_sig2[c(i, i + n_sigs, i + n_sigs * 2), ], mat_enrich_sig_pv2[i,
        ], main = row.names(mat_enrich_sig2)[i], ylim = c(0, max(mat_enrich_sig2[c(i, i + n_sigs, i +
        n_sigs * 2), ])))
}
dev.off()

# crypt
mat_enrich_sig <- fun_enrich_sig(enrich.mut.crypt[meta.mut.crypt[, 2] == "SNV", ], enrich.ctrl, meta.mut.crypt[meta.mut.crypt[,
    2] == "SNV", ], list_elements, list_genes, 1, 7, nmf_res_snv_with_crypt$signatures)
mat_enrich_sig_pv <- fun_enrich_sig_pv(enrich.mut.crypt[meta.mut.crypt[, 2] == "SNV", ], enrich.ctrl,
    meta.mut.crypt[meta.mut.crypt[, 2] == "SNV", ], list_elements, list_genes, 1, 7, nmf_res_snv_with_crypt$signatures)
mat_enrich_sig2 <- fun_enrich_sig(enrich.mut.crypt[meta.mut.crypt[, 2] != "SNV", ], enrich.ctrl, meta.mut.crypt[meta.mut.crypt[,
    2] != "SNV", ], list_elements, list_genes, 1, 7, nmf_res_indel_with_crypt$signatures)
mat_enrich_sig_pv2 <- fun_enrich_sig_pv(enrich.mut.crypt[meta.mut.crypt[, 2] != "SNV", ], enrich.ctrl,
    meta.mut.crypt[meta.mut.crypt[, 2] != "SNV", ], list_elements, list_genes, 1, 7, nmf_res_indel_with_crypt$signatures)
pdf("figures/enrichment_analysis/with_crypt/dog.crypt.enrichment_genes_signature.pdf", width = 4, height = 4,
    useDingbats = F)
par(mar = c(10, 4, 3, 1), tcl = 0.3, bty = "n", cex = 5/6)
n_sigs <- nrow(mat_enrich_sig)/3
for (i in 1:n_sigs) {
    plot_enrichment_with_ci(mat_enrich_sig[c(i, i + n_sigs, i + n_sigs * 2), ], mat_enrich_sig_pv[i,
        ], main = row.names(mat_enrich_sig)[i], ylim = c(0, 3))
}
n_sigs <- nrow(mat_enrich_sig2)/3
mat_enrich_sig_pv2[2, ] <- 1
for (i in 1:n_sigs) {
    plot_enrichment_with_ci(mat_enrich_sig2[c(i, i + n_sigs, i + n_sigs * 2), ], mat_enrich_sig_pv2[i,
        ], main = row.names(mat_enrich_sig2)[i], ylim = c(0, 10))
}
dev.off()

pdf("figures/enrichment_analysis/with_crypt/dog.crypt.enrichment_genes_signature_autoScaled.pdf", width = 4,
    height = 4, useDingbats = F)
par(mar = c(10, 4, 3, 1), tcl = 0.3, bty = "n", cex = 5/6)
n_sigs <- nrow(mat_enrich_sig)/3
for (i in 1:n_sigs) {
    plot_enrichment_with_ci(mat_enrich_sig[c(i, i + n_sigs, i + n_sigs * 2), ], mat_enrich_sig_pv[i,
        ], main = row.names(mat_enrich_sig)[i], ylim = c(0, max(mat_enrich_sig[c(i, i + n_sigs, i + n_sigs *
        2), ])))
}
n_sigs <- nrow(mat_enrich_sig2)/3
mat_enrich_sig_pv2[2, ] <- 1
for (i in 1:n_sigs) {
    plot_enrichment_with_ci(mat_enrich_sig2[c(i, i + n_sigs, i + n_sigs * 2), ], mat_enrich_sig_pv2[i,
        ], main = row.names(mat_enrich_sig2)[i], ylim = c(0, min(max(mat_enrich_sig2[c(i, i + n_sigs,
        i + n_sigs * 2), ]), 100)))
}
dev.off()

# on genomic elements; non-seperated by signatures----
mat_enrich <- fun_enrich(enrich.mut.crypt[meta.mut.crypt[, 2] == "SNV", ], enrich.ctrl, list_elements,
    list_genes, 1, 7)
mat_enrich_pv <- fun_enrich_pv_num(enrich.mut.crypt[meta.mut.crypt[, 2] == "SNV", ], enrich.ctrl, list_elements,
    list_genes, 1, 7)
#for(i in 1){mat_enrich_pv[i,]=p.adjust(mat_enrich_pv[i,],method='fdr')}
mat_enrich2 <- fun_enrich(enrich.mut.crypt[enrich.mut.crypt[, 2] != "SNV", ], enrich.ctrl, list_elements,
    list_genes, 1, 7)
mat_enrich_pv2 <- fun_enrich_pv_num(enrich.mut.crypt[meta.mut.crypt[, 2] != "SNV", ], enrich.ctrl, list_elements,
    list_genes, 1, 7)
#for(i in 1:2){mat_enrich_sig_pv2[i,]=p.adjust(mat_enrich_sig_pv2[i,],method='fdr')}
pdf("figures/enrichment_analysis/with_crypt/dog.crypt.enrichment_genes.pdf", width = 5, height = 4, useDingbats = F)
par(mar = c(10, 4, 3, 1), tcl = 0.3, bty = "n", cex = 5/6)
plot_enrichment_with_ci(mat_enrich, mat_enrich_pv[1, ], main = "sSNV", ylim = c(0, 2))
plot_enrichment_with_ci(mat_enrich2, mat_enrich_pv2[1, ], main = "sindel", ylim = c(0, 3))
dev.off()

# enrichment analysis for human----
# read data
require(MutationalPatterns)
enrich.mut <- read.table("../tables/enrichment/merged.sMutation.all.overlap.reducedCol.tab", header = T,
    row.names = 1)
meta.mut <- read.table("../tables/enrichment/mutation.metadata.tab", header = F, row.names = 1)
enrich.ctrl <- read.table("../tables/enrichment/control.sMutation.all.overlap.reducedCol.tab", header = T,
    row.names = 1)
meta.ctrl <- read.table("../tables/enrichment/control.metadata.tab", header = F, row.names = 1)
enrich.mut.crypt <- read.table("tables/enrichment/cagan_human.sMutation.all.overlap.reducedCol.tab",
    header = T, row.names = 1)
meta.mut.crypt <- read.table("tables/enrichment/cagan_human.metadata.tab", header = F, row.names = 1)

neuron.cpm <- read.table("../aging_project/tables/rnaseq/neurons.CPM.mRNA.txt", header = F, row.names = 1,
    check.names = F)
colnames(neuron.cpm) <- c("all_neuron", "ext", "inb")
neuron.cpm <- neuron.cpm[order(neuron.cpm[, 2]), ]
inter_n <- as.vector(read.table("../aging_project/tables/annotation/hs37.intergenic.bed", header = F,
    row.names = NULL)[, 4])
list_tpm <- chunk(row.names(neuron.cpm), 5)
gene_n <- row.names(neuron.cpm)
tn <- intersect(row.names(gtex.exp), gene_n)
list_tpm_crypt <- chunk(tn[order(gtex.exp[tn, "Colon - Sigmoid"])], 5)
list_tpm_bulk <- chunk(tn[order(gtex.exp[tn, "Brain - Frontal Cortex (BA9)"])], 5)

# on genomic elements----
list_elements <- list(mRNA.gene.1st_quantile = "mRNA.gene", mRNA.gene.2nd_quantile = "mRNA.gene", mRNA.gene.3rd_quantile = "mRNA.gene",
    mRNA.gene.4th_quantile = "mRNA.gene", mRNA.gene.5th_quantile = "mRNA.gene", mRNA.gene = "mRNA.gene",
    intergenic = "intergenic", aging_down.gene = "mRNA.gene", aging_up.gene = "mRNA.gene", aging_down.snRNA = "mRNA.gene",
    aging_up.snRNA = "mRNA.gene", mRNA.gene.len_1 = "mRNA.gene", mRNA.gene.len_2 = "mRNA.gene", mRNA.gene.len_3 = "mRNA.gene",
    mRNA.gene.len_4 = "mRNA.gene", mRNA.gene.len_5 = "mRNA.gene")
list_genes <- list(mRNA.gene.1st_quantile = list_tpm[[1]], mRNA.gene.2nd_quantile = list_tpm[[2]], mRNA.gene.3rd_quantile = list_tpm[[3]],
    mRNA.gene.4th_quantile = list_tpm[[4]], mRNA.gene.5th_quantile = list_tpm[[5]], mRNA.gene = gene_n,
    intergenic = inter_n, aging_down.gene = deseq_gtex_downgenes, aging_up.gene = deseq_gtex_upgenes,
    aging_down.snRNA = degs.down.human, aging_up.snRNA = degs.up.human, mRNA.gene.len_1 = list_len[[1]],
    mRNA.gene.len_2 = list_len[[2]], mRNA.gene.len_3 = list_len[[3]], mRNA.gene.len_4 = list_len[[4]],
    mRNA.gene.len_5 = list_len[[5]])

# neuron
mat_enrich_sig <- fun_enrich_sig(enrich.mut[meta.mut[, 2] == "SNV" & meta.mut[, 1] %in% cells_normal,
    ], enrich.ctrl, meta.mut[meta.mut[, 2] == "SNV" & meta.mut[, 1] %in% cells_normal, ], list_elements,
    list_genes, 1, 7, nmf_res_snv_with_crypt$signatures, correct = F)
mat_enrich_sig_pv <- fun_enrich_sig_pv(enrich.mut[meta.mut[, 2] == "SNV" & meta.mut[, 1] %in% cells_normal,
    ], enrich.ctrl, meta.mut[meta.mut[, 2] == "SNV" & meta.mut[, 1] %in% cells_normal, ], list_elements,
    list_genes, 1, 7, nmf_res_snv_with_crypt$signatures, correct = F)
mat_enrich_sig2 <- fun_enrich_sig(enrich.mut[meta.mut[, 2] != "SNV" & meta.mut[, 1] %in% cells_normal,
    ], enrich.ctrl, meta.mut[meta.mut[, 2] != "SNV" & meta.mut[, 1] %in% cells_normal, ], list_elements,
    list_genes, 1, 7, nmf_res_indel_with_crypt$signatures, correct = F)
mat_enrich_sig_pv2 <- fun_enrich_sig_pv(enrich.mut[meta.mut[, 2] != "SNV" & meta.mut[, 1] %in% cells_normal,
    ], enrich.ctrl, meta.mut[meta.mut[, 2] != "SNV" & meta.mut[, 1] %in% cells_normal, ], list_elements,
    list_genes, 1, 7, nmf_res_indel_with_crypt$signatures, correct = F)

pdf("figures/enrichment_analysis/with_crypt/human.neuron.enrichment_genes_signature.pdf", width = 4,
    height = 4, useDingbats = F)
par(mar = c(10, 4, 3, 1), tcl = 0.3, bty = "n", cex = 5/6)
n_sigs <- nrow(mat_enrich_sig)/3
for (i in 1:n_sigs) {
    plot_enrichment_with_ci(mat_enrich_sig[c(i, i + n_sigs, i + n_sigs * 2), ], mat_enrich_sig_pv[i,
        ], main = row.names(mat_enrich_sig)[i], ylim = c(0, 3))
}
n_sigs <- nrow(mat_enrich_sig2)/3
for (i in 1:n_sigs) {
    plot_enrichment_with_ci(mat_enrich_sig2[c(i, i + n_sigs, i + n_sigs * 2), ], mat_enrich_sig_pv2[i,
        ], main = row.names(mat_enrich_sig2)[i], ylim = c(0, 10))
}
dev.off()

# neuron; expression quintile by bulk data
list_genes <- list(mRNA.gene.1st_quantile = list_tpm_bulk[[1]], mRNA.gene.2nd_quantile = list_tpm_bulk[[2]],
    mRNA.gene.3rd_quantile = list_tpm_bulk[[3]], mRNA.gene.4th_quantile = list_tpm_bulk[[4]], mRNA.gene.5th_quantile = list_tpm_bulk[[5]],
    mRNA.gene = gene_n, intergenic = inter_n, aging_down.gene = deseq_gtex_downgenes, aging_up.gene = deseq_gtex_upgenes,
    aging_down.snRNA = degs.down.human, aging_up.snRNA = degs.up.human, mRNA.gene.len_1 = list_len[[1]],
    mRNA.gene.len_2 = list_len[[2]], mRNA.gene.len_3 = list_len[[3]], mRNA.gene.len_4 = list_len[[4]],
    mRNA.gene.len_5 = list_len[[5]])

mat_enrich_sig <- fun_enrich_sig(enrich.mut[meta.mut[, 2] == "SNV" & meta.mut[, 1] %in% cells_normal,
    ], enrich.ctrl, meta.mut[meta.mut[, 2] == "SNV" & meta.mut[, 1] %in% cells_normal, ], list_elements,
    list_genes, 1, 7, nmf_res_snv_with_crypt$signatures, correct = F)
mat_enrich_sig_pv <- fun_enrich_sig_pv(enrich.mut[meta.mut[, 2] == "SNV" & meta.mut[, 1] %in% cells_normal,
    ], enrich.ctrl, meta.mut[meta.mut[, 2] == "SNV" & meta.mut[, 1] %in% cells_normal, ], list_elements,
    list_genes, 1, 7, nmf_res_snv_with_crypt$signatures, correct = F)
mat_enrich_sig2 <- fun_enrich_sig(enrich.mut[meta.mut[, 2] != "SNV" & meta.mut[, 1] %in% cells_normal,
    ], enrich.ctrl, meta.mut[meta.mut[, 2] != "SNV" & meta.mut[, 1] %in% cells_normal, ], list_elements,
    list_genes, 1, 7, nmf_res_indel_with_crypt$signatures, correct = F)
mat_enrich_sig_pv2 <- fun_enrich_sig_pv(enrich.mut[meta.mut[, 2] != "SNV" & meta.mut[, 1] %in% cells_normal,
    ], enrich.ctrl, meta.mut[meta.mut[, 2] != "SNV" & meta.mut[, 1] %in% cells_normal, ], list_elements,
    list_genes, 1, 7, nmf_res_indel_with_crypt$signatures, correct = F)

pdf("figures/enrichment_analysis/with_crypt/human.neuron.enrichment_genes_signature_expQuintile_by_bulk.pdf",
    width = 4, height = 4, useDingbats = F)
par(mar = c(10, 4, 3, 1), tcl = 0.3, bty = "n", cex = 5/6)
n_sigs <- nrow(mat_enrich_sig)/3
for (i in 1:n_sigs) {
    plot_enrichment_with_ci(mat_enrich_sig[c(i, i + n_sigs, i + n_sigs * 2), ], mat_enrich_sig_pv[i,
        ], main = row.names(mat_enrich_sig)[i], ylim = c(0, 3))
}
n_sigs <- nrow(mat_enrich_sig2)/3
for (i in 1:n_sigs) {
    plot_enrichment_with_ci(mat_enrich_sig2[c(i, i + n_sigs, i + n_sigs * 2), ], mat_enrich_sig_pv2[i,
        ], main = row.names(mat_enrich_sig2)[i], ylim = c(0, 10))
}
dev.off()

# crypt
mat_enrich_sig <- fun_enrich_sig(enrich.mut.crypt[meta.mut.crypt[, 2] == "SNV", ], enrich.ctrl, meta.mut.crypt[meta.mut.crypt[,
    2] == "SNV", ], list_elements, list_genes, 1, 7, nmf_res_snv_with_crypt$signatures, correct = F)
mat_enrich_sig_pv <- fun_enrich_sig_pv(enrich.mut.crypt[meta.mut.crypt[, 2] == "SNV", ], enrich.ctrl,
    meta.mut.crypt[meta.mut.crypt[, 2] == "SNV", ], list_elements, list_genes, 1, 7, nmf_res_snv_with_crypt$signatures,
    correct = F)
mat_enrich_sig2 <- fun_enrich_sig(enrich.mut.crypt[meta.mut.crypt[, 2] != "SNV", ], enrich.ctrl, meta.mut.crypt[meta.mut.crypt[,
    2] != "SNV", ], list_elements, list_genes, 1, 7, nmf_res_indel_with_crypt$signatures, correct = F)
mat_enrich_sig_pv2 <- fun_enrich_sig_pv(enrich.mut.crypt[meta.mut.crypt[, 2] != "SNV", ], enrich.ctrl,
    meta.mut.crypt[meta.mut.crypt[, 2] != "SNV", ], list_elements, list_genes, 1, 7, nmf_res_indel_with_crypt$signatures,
    correct = F)
pdf("figures/enrichment_analysis/with_crypt/human.crypt.enrichment_genes_signature.pdf", width = 4, height = 4,
    useDingbats = F)
par(mar = c(10, 4, 3, 1), tcl = 0.3, bty = "n", cex = 5/6)
n_sigs <- nrow(mat_enrich_sig)/3
for (i in 1:n_sigs) {
    plot_enrichment_with_ci(mat_enrich_sig[c(i, i + n_sigs, i + n_sigs * 2), ], mat_enrich_sig_pv[i,
        ], main = row.names(mat_enrich_sig)[i], ylim = c(0, 3))
}
n_sigs <- nrow(mat_enrich_sig2)/3
for (i in 1:n_sigs) {
    plot_enrichment_with_ci(mat_enrich_sig2[c(i, i + n_sigs, i + n_sigs * 2), ], mat_enrich_sig_pv2[i,
        ], main = row.names(mat_enrich_sig2)[i], ylim = c(0, 10))
}
dev.off()

# on genomic elements; non-seperated by signatures----
mat_enrich <- fun_enrich(enrich.mut.crypt[meta.mut.crypt[, 2] == "SNV", ], enrich.ctrl, list_elements,
    list_genes, 1, 7)
mat_enrich_pv <- fun_enrich_pv_num(enrich.mut.crypt[meta.mut.crypt[, 2] == "SNV", ], enrich.ctrl, list_elements,
    list_genes, 1, 7)
#for(i in 1){mat_enrich_pv[i,]=p.adjust(mat_enrich_pv[i,],method='fdr')}
mat_enrich2 <- fun_enrich(enrich.mut.crypt[enrich.mut.crypt[, 2] != "SNV", ], enrich.ctrl, list_elements,
    list_genes, 1, 7)
mat_enrich_pv2 <- fun_enrich_pv_num(enrich.mut.crypt[meta.mut.crypt[, 2] != "SNV", ], enrich.ctrl, list_elements,
    list_genes, 1, 7)
#for(i in 1:2){mat_enrich_sig_pv2[i,]=p.adjust(mat_enrich_sig_pv2[i,],method='fdr')}
pdf("figures/enrichment_analysis/with_crypt/human.crypt.enrichment_genes.pdf", width = 5, height = 4,
    useDingbats = F)
par(mar = c(10, 4, 3, 1), tcl = 0.3, bty = "n", cex = 5/6)
plot_enrichment_with_ci(mat_enrich, mat_enrich_pv[1, ], main = "sSNV", ylim = c(0, 2))
plot_enrichment_with_ci(mat_enrich2, mat_enrich_pv2[1, ], main = "sindel", ylim = c(0, 3))
dev.off()

# enrichment of strand bias----
list_elements <- list(mRNA.gene.1st_quantile = "mRNA.gene", mRNA.gene.2nd_quantile = "mRNA.gene", mRNA.gene.3rd_quantile = "mRNA.gene",
    mRNA.gene.4th_quantile = "mRNA.gene", mRNA.gene.5th_quantile = "mRNA.gene")
list_genes <- list(mRNA.gene.1st_quantile = list_tpm[[1]], mRNA.gene.2nd_quantile = list_tpm[[2]], mRNA.gene.3rd_quantile = list_tpm[[3]],
    mRNA.gene.4th_quantile = list_tpm[[4]], mRNA.gene.5th_quantile = list_tpm[[5]])
tmat <- matrix(NA, 6, 5)
row.names(tmat) <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
colnames(tmat) <- names(list_elements)
tmat2 <- matrix(NA, 6, 5)
row.names(tmat2) <- c("C>A", "C>G", "C>T", "T>A", "T>C", "T>G")
colnames(tmat2) <- names(list_elements)
for (i in row.names(tmat)) {
    tn <- which(substring(meta.mut[, 3], 3, 5) == i)
    tmp <- fun_enrich_strand(enrich.mut[tn, ], list_elements, list_genes, 1, 7, 4)
    tmat[i, ] <- tmp[, 1]
    tmat2[i, ] <- tmp[, 2]
}
tv1 <- c()
for (i in 1:5) {
    tv1 <- c(tv1, tmat[, i])
}
tv2 <- c()
for (i in 1:5) {
    tv2 <- c(tv2, tmat2[, i])
}
tv3 <- rep("", 30)
tv3[p.adjust(tv2, method = "fdr") < 0.05] <- "*"
at <- c(1:6, 8:13, 15:20, 22:27, 29:34)
pdf("figures/enrichment_analysis/strand_analysis_human.pdf", width = 3.5, height = 6, useDingbats = F)
par(mar = c(4, 1, 1, 1), tcl = 0.3, bty = "n", cex = 5/6)
plot(100 - tv1 * 100, at, pch = 19, yaxt = "n", xlab = "%sSNV in transcribed strand", xlim = c(30, 70),
    col = c("#80b1d3", "black", "#e41a1c", "#d9d9d9", "#b3de69", "#fccde5"))
abline(h = c(7, 14, 21, 28))
text(30, at, label = tv3, cex = 3)
text(30, 1:6, label = row.names(tmat), pos = 4)
dev.off()

# power analysis for enrichment statistics----
min_O_for_sig <- function(OE, alpha = 0.05) {
    O <- 1
    while (TRUE) {
        E <- O/OE
        tab <- matrix(c(O, E, (O + E)/2, (O + E)/2), nrow = 2)
        p <- suppressWarnings(chisq.test(tab)$p.value)
        if (p < alpha)
            return(O)
        O <- O + 1
    }
}
tv <- c()
for (oe in 1 + (2:20)/10) {
    tv <- c(tv, min_O_for_sig(oe))
}
pdf("figures/power_analysis_enrichment.pdf", width = 4, height = 4, useDingbats = F)
par(mar = c(4, 4, 3, 1), tcl = 0.3, bty = "n")
plot(tv, ylim = c(0, 600), ylab = "min mutations needed to reach p < 0.05", main = "power analysis for enrichment",
    xaxt = "n", xlab = "observed / expected")
axis(1, 1:length(tv), label = 1 + (2:20)/10, lwd = 0, las = 3)
dev.off()

# Reproducibility record -------------------------------------------------------
writeLines(
    capture.output(sessionInfo()),
    con = "tables/For_manuscript/sessionInfo.txt"
)
message("Analysis complete. Outputs were written under: ", project_root)
