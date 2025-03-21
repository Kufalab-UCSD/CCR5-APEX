if (!exists("l_baselineCorrect")) { l_baselineCorrect <- TRUE }
if (!exists("maxnPerPlex")) { maxnPerPlex <- 15 }
if (!exists("maxnPerCond")) { maxnPerCond <- 3 }
if (!exists("r_density.adjust")) { r_density.adjust <- 1. }

##################### TRANFSORM and FILTER THE DATA ##########################
densityPlots <- c("none", "gui", "file")[3] ### a flag indicating where and whether intermediate density plots should be built

for (colnm in semantics$sample) { data[, colnm] <- as.numeric(data[,colnm]) ; data[(data[, colnm] == 0 & !is.na(data[, colnm])), colnm] <- NA } ### make into numerical form, replace 0's by NAs. Exclude NAs to prevent error.
data[, semantics$sample] <- log(data[, semantics$sample]) ### log-transform

### delete the so-called "blanks" if there are any
colnms <- colnames(data)[!(colnames(data) %in% semantics$sample)] ; semantics <- semantics[!grepl("blank", semantics$cell_line), ] ; data <- data[, c(colnms, semantics$sample)]

### if needed, generate density plots for individual cell lines and plexes ######
colnms <- c("plex", "cell_line", "replicate") ; for (colnm in colnms) { if (colnm %in% colnames(semantics)) { break } } ### find a column from a predefined dictionary that is in semantics
if (densityPlots!="none") { oufile <- "" ; if (densityPlots=="file") { oufile=paste0(dirOut, "/images/density_raw.png") } ; plot_density(data, semantics, ymax=0.9, adjust=r_density.adjust, semacol=colnm, oufile=oufile) ; }

if (TRUE) { ### identify samples with abnormal distributions - these may be compromised
  semantics$mode <- 0 ; semantics$mode_density <- 0 ; semantics$log_mode <- 0 ; semantics$log_mode_density <- 0 ; semantics$median <- 0 ; semantics$sigma <- 0
  for (i in 1:nrow(semantics)) {
    d <- density(data[, semantics$sample[i]], adjust=r_density.adjust, na.rm=TRUE) ; d <- data.frame(x=d$x, y=d$y) ; semantics$log_mode[i] <- d$x[order(d$y, decreasing=TRUE)] [1] ; semantics$log_mode_density[i] <- max(d$y)
    d <- density(exp(data[, semantics$sample[i]]), adjust=r_density.adjust, na.rm=TRUE) ; d <- data.frame(x=d$x, y=d$y) ;  semantics$mode[i] <- d$x[order(d$y, decreasing=TRUE)] [1] ; semantics$mode_density[i] <- max(d$y)
    semantics$median[i] <- median(data[, semantics$sample[i]], na.rm=TRUE) ; semantics$sigma[i] <- sd(data[, semantics$sample[i]], na.rm=TRUE)
  }

  ### 2023-06-01 re-implemented the chunk below
  colnms <- c("log_mode", "median", "sigma") ; if ("OD" %in% colnames(semantics)) colnms <- c(colnms, "OD", "concentration")
  data.stats <- as.data.frame(t(semantics[, colnms]), stringsAsFactors = FALSE)
  colnames(data.stats) <- semantics$sample ; data.stats$AC <- rownames(data.stats) ; rownames(data.stats) <- 1:nrow(data.stats)
  data.stats[, semantics$sample] <- as.data.frame(lapply(data.stats[, semantics$sample], function(x) as.numeric(x))) ### t() coerces numbers to strings, this is to fix it
  data.stats.plot <- plotPoint(data.stats, semantics, ncol=1) ; print(data.stats.plot) ### this print does not work by some reason, but data.stats.plot is created alright

  ### older version of the same
  # modes <- as.data.frame(t(semantics[, 'log_mode', drop=FALSE]), stringsAsFactors = FALSE) ; colnames(modes) <- semantics$sample ; rownames(modes) <- 1:1
  # modes <- as.data.frame(lapply(modes, function(x) as.numeric(x))) ### t() coerces numbers to strings, this is to fix it
  # medians <- as.data.frame(t(semantics[, 'median', drop=FALSE]), stringsAsFactors = FALSE) ; colnames(medians) <- semantics$sample ; rownames(medians) <- 1:1
  # medians <- as.data.frame(lapply(medians, function(x) as.numeric(x))) ### t() coerces numbers to strings, this is to fix it
  # sds <- as.data.frame(t(semantics[, 'sigma', drop=FALSE]), stringsAsFactors = FALSE) ; colnames(sds) <- semantics$sample ; rownames(sds) <- 1:1
  # sds <- as.data.frame(lapply(sds, function(x) as.numeric(x))) ### t() coerces numbers to strings, this is to fix it
  # modes$AC <- "modes" ; medians$AC <- "medians" ; sds$AC <- "SDs"
  # plotPoint(modes, semantics, 1) ; plotPoint(medians, semantics, 1) ; plotPoint(sds, semantics, 1)
  # if ("OD" %in% colnames(semantics)) {
  #   ODs <- as.data.frame(t(semantics[, 'OD', drop=FALSE]), stringsAsFactors = FALSE) ; colnames(ODs) <- semantics$sample ; rownames(ODs) <- 1:1
  #   ODs <- as.data.frame(lapply(ODs, function(x) as.numeric(x))) ### t() coerces numbers to strings, this is to fix it
  #   concentrations <- as.data.frame(t(semantics[, 'concentration', drop=FALSE]), stringsAsFactors = FALSE) ; colnames(concentrations) <- semantics$sample ; rownames(concentrations) <- 1:1
  #   concentrations <- as.data.frame(lapply(concentrations, function(x) as.numeric(x))) ### t() coerces numbers to strings, this is to fix it
  #   plotPoint(ODs, semantics, 1) ; plotPoint(concentrations, semantics, 1) ; rm(concentrations, ODs)
  # }
  # #colnms <- colnames(data)[!(colnames(data) %in% semantics$sample)] ; semantics <- semantics[semantics$mode > 80., ] ; data <- data[, c(colnms, semantics$sample)] ### 80 seemed like a good cutoff for protein-level quants on a linear scale; must be adjusted for other versions of the data
  # rm(medians, modes, sds)

  semantics <- within(semantics, rm(plex_channel, mode, mode_density, log_mode, log_mode_density, median, sigma))
  if ("OD" %in% colnames(semantics)) semantics <- within(semantics, rm(OD, concentration))
}

#### evaluate and plot peptide/protein coverage by plexes and conditions #####
##### then remove peptide/condition combos with insufficient replicates, #####
########################### and do it again ##################################

known_cell_lines <- c("CCR5", "CCR2wG", "CCR2noG",
                      "PM", "Endo", "Cyto", "WT",
                      "dbArr", "normalG", "deltaG", "noAPEX",
                      "BioID2-Gai", "BioID2-Gai-QL", "BioID2-Caax",
                      "PP", "AA", "Vec",
                      "Gai.119.FA", "Gai.92.FA", "NT.FA.Gai",
                      "B2AR", "DOR")

if (FALSE) { # filter peptide/condition combinations with insufficient replicates
  # plot peptide coverage by plexes and conditions pre-filtering
  print(paste0(sum(rowSums(!is.na(data[, semantics$sample]))), " non-NA values pre filtering insufficient replicates"))
  semanticsX <- semantics[semantics$cell_line %in% known_cell_lines, ] ; semanticsX <- semanticsX[order(semanticsX$sample), ] ; dataX <- data[, semanticsX$sample]
  cvg_plex <- plot_coverage(dataX, semanticsX, semacol='plex',      maxn=maxnPerPlex, paste(dirOut, "images", paste0(kwd, '_coverage_plex_preFilt.png'), sep='/'))
  cvg_cond <- plot_coverage(dataX, semanticsX, semacol='condition', maxn=maxnPerCond, paste(dirOut, "images", paste0(kwd, '_coverage_condition_preFilt.png'), sep='/'))
  # filter
  for (cond in unique(semantics$condition)) {
    colnms <- semantics$sample[semantics$condition == cond]
    data[rowSums(!is.na(data[,colnms, drop=FALSE])) < 2, colnms] <- NA
  }
  rm(cond, colnms)
}

print(paste0(sum(rowSums(!is.na(data[, semantics$sample]))), " non-NA values in the dataset"))
semanticsX <- semantics[semantics$cell_line %in% known_cell_lines, ] ; semanticsX <- semanticsX[order(semanticsX$sample), ] ; dataX <- data[, semanticsX$sample]
if (length(unique(semanticsX$plex))>1) { cvg_plex <- plot_coverage(dataX, semanticsX, semacol='plex', maxn=maxnPerPlex, paste0(dirOut, "/images/coverage_plex.png")) }
cvg_cond <- plot_coverage(dataX, semanticsX, semacol='condition', maxn=maxnPerCond, paste0(dirOut, "/images/coverage_condition.png"))
rm(semanticsX, dataX, cvg_plex, cvg_cond, known_cell_lines)

####### mean-center columns within plexes (H2O2 and no-H2O2 separately) #######
centerMethod <- c("mean", "median", "mode")[2]
colnms <- semantics$sample ### only the data columns
means <- colMeans(data[, colnms], na.rm=TRUE) ### by default, use arithmetic means
if (centerMethod == "median") { for (i in 1:nrow(semantics)) { colnm <- semantics$sample[i] ; means[i] <- median(data[, colnm], na.rm=TRUE) }}
if (centerMethod == "mode") { ### use modes deduced from the kernel density function
  for (i in 1:nrow(semantics)) {
    d <- density(data[, semantics$sample[i]], na.rm=TRUE) ; d <- data.frame(x=d$x, y=d$y)
    semantics$mode[i] <- d$x[order(d$y, decreasing=TRUE)] [1] ; semantics$mode_density[i] <- max(d$y)
  }
  means <- semantics$mode ### modes of the data columns
}

plexes <- unique(semantics$plex) ; plexes <- plexes[order(plexes)]
for (p in plexes) {
  ixs <- ((semantics$plex==p) & (!grepl("noH2O2|noAPEX", semantics$cell_line))) ### only H2O2-treated data columns from plex p
  if (sum(ixs) > 0) { means[ixs] <- means[ixs] - mean(means[ixs], na.rm=TRUE) }
  ixs <- ((semantics$plex==p) & (grepl("noH2O2|noAPEX", semantics$cell_line))) ### only non-H2O2-treated data columns from plex p
  if (sum(ixs) > 0) { means[ixs] <- means[ixs] - mean(means[ixs], na.rm=TRUE) }
}
means <- t(as.data.frame(means))[rep(1, nrow(data)),]
means[,colSums(is.na(means)) == nrow(means)] <- NA
data[, semantics$sample] <- data[, semantics$sample]-means #; plot_all(data, semantics)
rm(means, p)
print(paste0(sum(rowSums(!is.na(data[, semantics$sample]))), " non-NA values post plex-specific mean-centering "))
if (densityPlots!="none") { oufile <- "" ; if (densityPlots=="file") { oufile=paste0(dirOut, "/images/density_post_", centerMethod, "Centering.png") } ; plot_density(data, semantics, ymax=0.9, adjust=r_density.adjust, semacol='plex', oufile=oufile) ; }

############# subtract no-H2O2 quants, remove the no-H2O2 columns #############
### DISADVANTAGE of no-H2O2 subtraction: it clears off all peptides for which noH2O2|noAPEX quant in the same plex is NA
### 2022-05-09: this can be avoided by keeping those peptides whose noH2O2|noAPEX quant in the same plex is NA *unchanged*
### ADVANTAGE of no-H2O2 subtraction: we can distinguish specific preys from non-specific ones by monitoring how positive the difference with no-H2O2 is
### BUT: we never do it so maybe there is no point in no-H2O2 subtraction. Or otherwise need to implement the check.
if (l_baselineCorrect) {
  for (p in plexes) {
    colnms <- semantics$sample[(semantics$plex==p) & (grepl("noH2O2|noAPEX", semantics$cell_line))] ### only non-H2O2-treated data columns from plex p
    rmeans <- rowMeans(data[, colnms, drop=FALSE], na.rm=TRUE)
    colnms <- semantics$sample[semantics$plex==p] ### all data columns from plex p
    data[is.na(rmeans), colnms] <- NA ; ### clear off those peptides whose noH2O2|noAPEX quant in the same plex is NA. When commented out, these peptides remain unchanged
    rmeans[is.na(rmeans)] <- 0.
    for (colnm in colnms) { data[, colnm] <- data[, colnm] - rmeans }
  }
  rm(rmeans, p)
  print(paste0(sum(rowSums(!is.na(data[, semantics$sample]))), " non-NA values post no-H2O2 quant subtraction"))
  colnms <- semantics$sample[!grepl("noH2O2|noAPEX", semantics$cell_line)] ; tmp <- data[, colnms] ; tmp[is.na(tmp)] <- 999.;
  data[rje::rowMins(tmp)<0, semantics$sample] <- NA ; rm(tmp) ### clear off everything that has at least one negative quant after noH2O2|noAPEX subtraction
  print(paste0(sum(rowSums(!is.na(data[, semantics$sample]))), " non-NA values post removal of rows with negative quants"))
  if (densityPlots!="none") { oufile <- "" ; if (densityPlots=="file") { oufile=paste(dirOut, "images/density_post_noH2O2_subtraction.png", sep='/') } ; plot_density(data, semantics, ymax=2.5, adjust=r_density.adjust, semacol='plex', oufile=oufile) ; }
}
colnms <- colnames(data)[!(colnames(data) %in% semantics$sample)] ; semantics <- semantics[!grepl("noH2O2|noAPEX", semantics$cell_line), ] ; data <- data[, c(colnms, semantics$sample)]
if (densityPlots!="none") { oufile <- "" ; if (densityPlots=="file") { oufile=paste(dirOut, "images/density_post_noH2O2_removal.png", sep='/') } ; plot_density(data, semantics, ymax=1.25, adjust=r_density.adjust, semacol='plex', oufile=oufile) ; }

#################### mean-center columns across all plexes ####################

print(paste0(sum(rowSums(!is.na(data[, semantics$sample]))), " non-NA values post no-H2O2 column removal"))
colnms <- semantics$sample ### only the data columns
means <- colMeans(data[, colnms], na.rm=TRUE) ### by default, use arithmetic means
if (centerMethod == "median") { for (i in 1:nrow(semantics)) { colnm <- semantics$sample[i] ; means[i] <- median(data[, colnm], na.rm=TRUE) }}
if (centerMethod == "mode") { ### use modes deduced from the kernel density function
  for (i in 1:nrow(semantics)) {
    d <- density(data[, semantics$sample[i]], na.rm=TRUE) ; d <- data.frame(x=d$x, y=d$y)
    semantics$mode[i] <- d$x[order(d$y, decreasing=TRUE)] [1] ; semantics$mode_density[i] <- max(d$y)
  }
  means <- semantics$mode ### modes of the data columns
}
means <- means-mean(means, na.rm=TRUE)
means <- t(as.data.frame(means))[rep(1, nrow(data)),]
means[,colSums(is.na(means)) == nrow(means)] <- NA
data[, semantics$sample] <- data[, semantics$sample]-means
rm(means)

print(paste0(sum(rowSums(!is.na(data[, semantics$sample]))), " non-NA values post across-plexes mean-centering "))

#colnms <- c("plex", "cell_line", "replicate") ; for (colnm in colnms) { if (colnm %in% colnames(semantics)) { break } } ### find a column from a predefined dictionary that is in semantics
#if (densityPlots!="none") { oufile <- "" ; if (densityPlots=="file") { oufile=paste0(dirOut, "/images/density_post_", centerMethod, "Centering.png") } ; plot_density(data, semantics, ymax=0.9, adjust=r_density.adjust, semacol=colnm, oufile=oufile) ; }

save(list=c("data", "semantics", "l_baselineCorrect", "data.stats.plot", "data.stats"), file=paste(dirOut, '03_transformFilterBaselineCorrect.Rdata', sep="/"))
