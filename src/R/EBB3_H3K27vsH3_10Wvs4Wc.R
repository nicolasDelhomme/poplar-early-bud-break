#' ---
#' title: "Abundance of H3K27me3 in EBB3 gene"
#' author: "Katja Stojkovič"
#' date: "`r Sys.Date()`"
#' output:
#'  html_document:
#'    toc: true
#'    number_sections: true
#' ---


#' # Setup  
#' Load the libraries
suppressPackageStartupMessages(library(GenomicAlignments))
suppressPackageStartupMessages(library(RColorBrewer))
suppressPackageStartupMessages(library(genomeIntervals))
suppressPackageStartupMessages(library(stringr))

#' Read in the file names...
bam_files <- list.files("/mnt/picea/home/miskolczi/histone-methylation_d/star/subsamples/",
                        pattern = ".out.bam$",
                        full.names = TRUE)

#' ...and name them
names(bam_files) <- sub("_sortmerna_trimmomatic_10_millionAligned.out.bam","",
                        sapply(strsplit(bam_files, "/"), .subset, 10))

#' Select the samples
bam_files_H3K27 <- bam_files[grepl("K", names(bam_files))]
bam_files_H3 <- bam_files[grepl("H", names(bam_files))]

#' Read in sample information
samples_ChIP_H3K27 <- read.csv("~/Git/UPSCb/projects/histone-methylation/src/R/figures/sample-information_ChIP-Seq.csv", row.names = 1)
samples_ChIP_H3K27 <- droplevels(samples_ChIP_H3K27[grepl("10", samples_ChIP_H3K27$samples), ])
samples_ChIP_H3 <- samples_ChIP_H3K27
samples_ChIP_H3$samples <- names(bam_files_H3)[grep("10", names(bam_files_H3))]


#' Read in gene information file
Potrx01 <- read.delim("/mnt/picea/home/miskolczi/projects/aspseq/rbhalerao/histone-methylation/Reference/Potrx01-gene_16-11-18.gff3", header = FALSE)

Potrx01 <- Potrx01[-1, -c(2:3,6,8)]
rownames(Potrx01) <- substr(Potrx01$V9, start = 4, stop = 20)
Potrx01 <- Potrx01[ , 1:4]
colnames(Potrx01) <- c("chr", "start", "end", "strand")

#' Numbers for correcting coverage for abundance differences defined earlier in the analysis
mat.basepair_10Million <- c(13342925, 10920785, 13145706, 19498055 ,19897640, 20539094, 22537301, 21745033, 23896687, 18119787, 17863719, 15323753) 
names(mat.basepair_10Million) <- paste0(rep(c("0","6","10","104"), each = 3), rep(c("A","B","C"), 4)) 
abund.allsamp.offset<- 1/(mat.basepair_10Million/mean(mat.basepair_10Million))
# choose only those important for the two time points
abund.allsamp.offset <- abund.allsamp.offset[grepl("10", names(abund.allsamp.offset))]
# Correction will be calculated this way: cov_gA_corr <- mapply("/", cov_gA, abund.allsamp.offset)

#' Choose a colour palette
pal <- brewer.pal(12, "Paired")

#'---------------------------------------------------------------------  
#'
#' # EBB3

#' Define gene and region of interest  
# write name of the gene
gene_name <- "EBB3"
# write gene accession number
gene <- "Potrx003580g02767"
# define length of upstream and downstream regions
upstream <- 2000
downstream <- 2000

# get chromosome, start and end of the gene
chr <- as.character(Potrx01[gene, "chr"])
TSS <- Potrx01[gene, "start"]
TTS <- Potrx01[gene, "end"]

# get the region of interest
region <- paste0(chr,":", TSS-upstream, "-", TTS+downstream)

#' ## H3K27
#' Read in region of interest from .bam files  
gA_H3K27 <- lapply(as.list(bam_files_H3K27), function(bam) {
  readGAlignmentPairs(file = bam, 
                      index = bam,
                      param = ScanBamParam(flag = scanBamFlag(isUnmappedQuery = FALSE, isPaired = TRUE, hasUnmappedMate = FALSE),
                                           which= GRanges(seqnames = region))
  )
})

#' Compute coverage  
cov_gA_H3K27 <- lapply(gA_H3K27, FUN = coverage)
# extract covarage in region of interest
cov_gA_H3K27 <- lapply(cov_gA_H3K27, function(x) {
  x[GRanges(region)]
})

#' Correct coverage for abundance differences defined earlier in the analysis
cov_gA_corr_H3K27 <- mapply("/", cov_gA_H3K27, abund.allsamp.offset)

df_cov_H3K27 <- as.data.frame(lapply(cov_gA_corr_H3K27, function(x) as.numeric(unlist(x, use.names = FALSE))))


#' ## H3  
#' Read in region of interest from .bam files  
gA_H3 <- lapply(as.list(bam_files_H3), function(bam) {
  readGAlignmentPairs(file = bam, 
                      index = bam,
                      param = ScanBamParam(flag = scanBamFlag(isUnmappedQuery = FALSE, isPaired = TRUE, hasUnmappedMate = FALSE),
                                           which= GRanges(seqnames = region))
  )
})

#' Compute coverage  
cov_gA_H3 <- lapply(gA_H3, FUN = coverage)
# extract covarage in region of interest
cov_gA_H3 <- lapply(cov_gA_H3, function(x) {
  x[GRanges(region)]
})

#' Correct coverage for abundance differences defined earlier in the analysis
cov_gA_corr_H3 <- mapply("/", cov_gA_H3, abund.allsamp.offset)

# make a data frame
df_cov_H3 <- as.data.frame(lapply(cov_gA_corr_H3, function(x) as.numeric(unlist(x, use.names = FALSE))))

#' ## H3K27/H3  
#' Divide abundance of H3K27 by H3  
#' Calculate log2 (and add 1), so there is no problems with 0 
df_log_H3K27 <- log2(df_cov_H3K27 + 1)
df_log_H3 <- log2(df_cov_H3 + 1)

#' subtract H3 from H3K27 (subtract, not divide, because thevalues are logarithmic)
df_cov_ratio <- df_log_H3K27 - df_log_H3

#' Calculate mean of three replicates in each time point  
df_cov_mean_ratio <- do.call(
  cbind,
  lapply(split.data.frame(t(df_cov_ratio),
                          samples_ChIP_H3K27$time_point),
         colMeans))
df_cov_mean_ratio <- as.data.frame(df_cov_mean_ratio)


#' ## correct H3K27/H3 by the offset H3-mean H3  
#' Calculate mean of H3 replicates in each time point  
df_log_mean_H3 <- do.call(
  cbind,
  lapply(split.data.frame(t(df_log_H3),
                          samples_ChIP_H3$time_point),
         colMeans))
df_log_mean_H3 <- as.data.frame(df_log_mean_H3)

#' Calculate mean of H3 of all time points
log_mean_H3_timep <- apply(df_log_mean_H3, 1, mean)
#' Calculate offset from the mean
df_log_mean_H3_offset <- log_mean_H3_timep - df_log_mean_H3

#' Subtract offset from H3K27 values
df_log_mean_ratio_offset <- df_cov_mean_ratio - df_log_mean_H3_offset

#' # Bar plot  

# split region in bins, show mean +/- SD  
length_region <- nrow(df_log_mean_ratio_offset)
# Length of the region is 8565 nt. For the purpose of plotting, remove 1 nt, because the length cannot be divided in bins.
df_cut <- df_log_mean_ratio_offset[1:(nrow(df_log_mean_ratio_offset)-15), ]
length_region <- nrow(df_cut)
# bin_size <- nrow(df_cut)/23 #size of the bin = length/nr of bins
bin_size <- 50
split_region <- split(df_cut, ceiling((1:length_region)/bin_size))

# mean
mean_bin <- sapply(split_region, function(bin) {
  apply(bin, 2, mean)
})
# sd
sd_bin <- sapply(split_region, function(bin) {
  apply(bin, 2, sd)
})

stat_10w <- t(rbind(mean = data.frame(mean_bin)["SD-10W", ], sd = data.frame(sd_bin)["SD-10W", ]))
stat_4wc <- t(rbind(mean = data.frame(mean_bin)["SD-10Wc4W", ], sd = data.frame(sd_bin)["SD-10Wc4W", ]))

plot_mean <- barplot(stat_10w[ , "mean"], 
                     col = rgb(177, 89, 40, maxColorValue = 255, alpha = 70), 
                     ylim = c(-1, 4.0),
                     xaxt = "n",
                     ylab = "normalised H3K27me3 abundance")
arrows(x0=plot_mean,
       y0=stat_10w[ , "mean"]+stat_10w[ , "sd"],
       y1=stat_10w[ , "mean"]-stat_10w[ , "sd"], 
       angle=90,
       code=3,
       length=0.03,
       col = rgb(177, 89, 40, maxColorValue = 255, alpha = 255))

plot_mean2 <- barplot(stat_4wc[ , "mean"], 
                      col = rgb(31, 120, 180, maxColorValue = 255, alpha = 50), 
                      xaxt = "n",
                      add = TRUE)
arrows(x0=plot_mean2,
       y0=stat_4wc[ , "mean"]+stat_4wc[ , "sd"],
       y1=stat_4wc[ , "mean"]-stat_4wc[ , "sd"], 
       angle=90,
       code=3,
       length=0.03,
       col = rgb(31, 120, 180, maxColorValue = 255, alpha = 200))

abline(v = c((upstream/bin_size)+1, ((upstream+(TTS-TSS))/bin_size)+1), col = "dark grey", lty = 2, lwd = 1)
text((upstream/bin_size)+1, -0.5, "TSS", col = "dark grey", adj = c(0, -.1))
text(((upstream+(TTS-TSS))/bin_size)+1, -0.5, "TTS", col = "grey", adj = c(0, -.1))

legend("topleft",
       legend = c("10WSD", "4WC"),
       fill = c(rgb(177, 89, 40, maxColorValue = 255, alpha = 70), rgb(31, 120, 180, maxColorValue = 255, alpha = 50)),
       border = "black",
       bty = "n")

# Window size = 100

#' # Bar plot  

# split region in bins, show mean +/- SD  
length_region <- nrow(df_log_mean_ratio_offset)
# Length of the region is 8565 nt. For the purpose of plotting, remove 1 nt, because the length cannot be divided in bins.
df_cut <- df_log_mean_ratio_offset[1:(nrow(df_log_mean_ratio_offset)-65), ]
length_region <- nrow(df_cut)
# bin_size <- nrow(df_cut)/23 #size of the bin = length/nr of bins
bin_size <- 100
split_region <- split(df_cut, ceiling((1:length_region)/bin_size))

# mean
mean_bin <- sapply(split_region, function(bin) {
  apply(bin, 2, mean)
})
# sd
sd_bin <- sapply(split_region, function(bin) {
  apply(bin, 2, sd)
})

stat_10w <- t(rbind(mean = data.frame(mean_bin)["SD-10W", ], sd = data.frame(sd_bin)["SD-10W", ]))
stat_4wc <- t(rbind(mean = data.frame(mean_bin)["SD-10Wc4W", ], sd = data.frame(sd_bin)["SD-10Wc4W", ]))

plot_mean <- barplot(stat_10w[ , "mean"], 
                     col = rgb(177, 89, 40, maxColorValue = 255, alpha = 70), 
                     ylim = c(-1, 4.0),
                     xaxt = "n",
                     ylab = "normalised H3K27me3 abundance")
arrows(x0=plot_mean,
       y0=stat_10w[ , "mean"]+stat_10w[ , "sd"],
       y1=stat_10w[ , "mean"]-stat_10w[ , "sd"], 
       angle=90,
       code=3,
       length=0.03,
       col = rgb(177, 89, 40, maxColorValue = 255, alpha = 255))

plot_mean2 <- barplot(stat_4wc[ , "mean"], 
                      col = rgb(31, 120, 180, maxColorValue = 255, alpha = 50), 
                      xaxt = "n",
                      add = TRUE)
arrows(x0=plot_mean2,
       y0=stat_4wc[ , "mean"]+stat_4wc[ , "sd"],
       y1=stat_4wc[ , "mean"]-stat_4wc[ , "sd"], 
       angle=90,
       code=3,
       length=0.03,
       col = rgb(31, 120, 180, maxColorValue = 255, alpha = 200))

abline(v = c((upstream/bin_size)+1, ((upstream+(TTS-TSS))/bin_size)+1), col = "dark grey", lty = 2, lwd = 1)
text((upstream/bin_size)+1, -0.5, "TSS", col = "dark grey", adj = c(0, -.1))
text(((upstream+(TTS-TSS))/bin_size)+1, -0.5, "TTS", col = "grey", adj = c(0, -.1))

legend("topleft",
       legend = c("10WSD", "4WC"),
       fill = c(rgb(177, 89, 40, maxColorValue = 255, alpha = 70), rgb(31, 120, 180, maxColorValue = 255, alpha = 50)),
       border = "black",
       bty = "n")
# For the whole region:
boxplot(df_log_mean_ratio_offset,
        names = c("10WSD", "4WC"),
        ylim = c(-1, 4),
        ylab = "normalised H3K27me3 abundance",
        xlab = "time point")