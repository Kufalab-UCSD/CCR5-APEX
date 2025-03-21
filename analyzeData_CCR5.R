#library(stats)
#library(reshape2)
#library(Rcpp) ### for melt?
library(ggplot2)
library(ProteomicsToolkit)

ggplot2::theme_update(text=ggplot2::element_text(size=16)) ### this is a ggplot command; use larger font for all plots

#################################################################################
################## set up directories and file names ############################
#################################################################################

## Kufalab data is stored in a shared Google Drive.
if (dir.exists('G:/Shared drives')) { Gdrive <- 'G:/Shared drives/' } else { Gdrive <- '/Volumes/GoogleDrive/Shared drives/' }
if (dir.exists(paste0("C:/Users/", Sys.info()[["user"]], "/Documents/GitHub/"))) { github <- paste0("C:/Users/", Sys.info()[["user"]], "/Documents/GitHub/")} else { github <- '~/GitHub/'}

### code directory
dirCode <- paste0(github, "/CCR5_APEX/") ##New home for analysis scripts

kwd <-  c("peptides", "proteins", "PSMs")[3] #Default to PSMs

### output directory and uniprot directory need to be local for speed and to prevent unending crashes of GDrive for Windows
dirOut <- paste0(gsub("\\", "/", Sys.getenv("HOME"), fixed=TRUE), "/CCR5_APEX_LOCAL/")
dirUniprot <- paste0(gsub("\\", "/", Sys.getenv("HOME"), fixed=TRUE), "/uniprot/")

### if the output dir, the image dir in it, or the uniprot dir do not exist, create them
if (!dir.exists(file.path(dirOut))) { dir.create(file.path(dirOut), recursive=TRUE) }
if (!dir.exists(file.path(paste0(dirOut, "images")))) { dir.create(file.path(paste0(dirOut, "images")), recursive=TRUE) }
if (!dir.exists(file.path(dirUniprot))) { dir.create(file.path(dirUniprot), recursive=TRUE) }

#################################################################################
########################## read and rename files ################################
#################################################################################
# source calls are commented out to avoid accidentlly overwritting cached pipeline steps
# Uncomment lines to run.

#source(paste0(dirCode, "01_readAndRename_2021_CCR5_APEX.R"))
load(file=paste(dirOut, '01_readAndRename.Rdata', sep="/"))

#################################################################################
########################### re-map UniProt names ################################
#################################################################################

#source(paste0(dirCode, "02_mapUniprot.R"))
load(file=paste0(dirOut, '02_mapUniprot.Rdata'))

#################################################################################
############################ remove contaminants ################################
#################################################################################

if (FALSE) { ### remove high-average-residual (calculated on a full set below) peptides
  ###############################################################################
  load(file=paste(dirOut, 'resi.Rdata', sep="/")) ### this reads the resi data frame
  resi$grp <- paste0(resi$AC, " ", resi$sequence) ; data$grp <- paste0(data$AC, " ", data$sequence)
  resi <- resi[, c("grp", "maxAbsResi")] ; data <- plyr::join(data, resi, by='grp', type="left")

  if (FALSE) { ### sanity check visualization
    tmp <- data
    for (colnm in semantics$sample) { tmp[, colnm] <- as.numeric(tmp[,colnm]) } ; tmp[tmp == 0] <- NA ### make into numerical form, replace 0's by NAs
    tmp[, semantics$sample] <- log(tmp[, semantics$sample]) ### log-transform
    tmp$avgQuant <- rowMeans(tmp[, semantics$sample], na.rm=TRUE)
    print(ggplot2::ggplot(data=tmp, ggplot2::aes(x=avgQuant, y=maxAbsResi)) + ggplot2::geom_point() + ggplot2::labs(x="average log quant", y="max absolute residual"))
    rm(tmp)
  }

  data[(!is.na(data$maxAbsResi)) & (data$maxAbsResi > 1.), semantics$sample] <- NA ### clear off the peptides with high residuals
  data <- within(data, rm(grp, maxAbsResi)) ; rm(resi)
  ###############################################################################
} else { ### remove random peptides in accordance with the exclusion list
  ### 2022-08-12 added the quants in the no-H2O2 condition as a proxy for badness, too

  ### append the log-transformed, normalized, and averaged quants without H2O2
  colnms <- semantics$sample[grepl("noH2O2", semantics$cell_line)] ; tmp <- data.frame(AC=data$AC)
  for (colnm in colnms) {
    tmp[, colnm] <- log(data[, colnm]) ; tmp[!is.na(data[, colnm]) & data[, colnm] == 0, colnm] <- NA
    tmp[, colnm] <- 0.75*(tmp[, colnm]-stats::quantile(tmp[, colnm], probs=0.5, na.rm=TRUE)) / sd(tmp[, colnm], na.rm=TRUE)
  } ### like z-score but using median not mean
  #plot_density(tmp, semantics[semantics$cell_line == "noH2O2", ])
  data <- cbind(data.frame(avgNoH2O2=rowMeans(tmp[, colnms], na.rm=TRUE)), data)
  # ... maybe this should have been at the level of proteins not peptides ... but currently not implemented

  ### load and append the CRAPome scores
  load(file=paste0(dirOut, "/../APEX_exclusionList_LOCAL/proteins_with_perc_excluded.Rdata"))
  tmp <- plyr::join(data, proteins, by='AC', type="left") ; tmp$perc_excluded[is.na(tmp$perc_excluded)] <- 0 ; tmp$badness_score[is.na(tmp$badness_score)] <- 0 ; tmp$status[is.na(tmp$status)] <- "regular"
  data <- cbind(tmp[, c("status", "badness_score", "perc_excluded")], data)

  ### combine CRAPome scores with no-H2O2 quants
  #print(ggplot2::ggplot(data=data, ggplot2::aes(x=avgNoH2O2, y=badness_score)) + ggplot2::geom_point() + ggplot2::labs(x="average log quant without H2O2", y="badness"))
  data$avgNoH2O2[is.na(data$avgNoH2O2)] <- -1001. ; data$badness_score[is.na(data$badness_score)] <- 0.
  data$badness_score <- rje::rowMaxs(data[, c("avgNoH2O2", "badness_score")]) ; data$avgNoH2O2[data$avgNoH2O2 < -1000.] <- NA ;

  ### re-calculate perc_excluded - make it more stringent
  trim_percentile <- 0.01
  badness_max <- stats::quantile(proteins$badness_score, probs=1.-trim_percentile, na.rm=TRUE)
  badness_min <- stats::quantile(proteins$badness_score, probs=0.25, na.rm=TRUE)
  data$perc_excluded <- 100.*data$badness_score/(badness_max-badness_min) ### 0-100
  data$perc_excluded[is.na(data$perc_excluded)] <- 0. ; data$perc_excluded[data$perc_excluded > 100.] <- 100.
  data$perc_excluded <- 0.98*data$perc_excluded ; ixs <- (data$status == "protected") ; data$perc_excluded[ixs] <- data$perc_excluded[ixs]*1.5-55.
  data$perc_excluded[data$perc_excluded < 0.] = 0. ; data$perc_excluded[data$perc_excluded > 100.] = 100.

  badness <- data ; for (colnm in semantics$sample) { ixs <- !is.na(badness[, colnm]) ; badness[ixs, colnm] <- badness$badness_score[ixs] }
  if (FALSE) { oufile=paste(dirOut, "images/badness_preExclusionList.png", sep='/') ; plot_density(badness, semantics, ymax=1.5, adjust=1., semacol='ligand', oufile=oufile) ; rm(oufile) }

  tmp <- data ; tmp$exclude <- 100*pracma::rand(1, nrow(tmp))[1, ] ; tmp$exclude <- tmp$exclude - tmp$perc_excluded ; tmp$exclude <- tmp$exclude < 0.
  if (FALSE) { ### 2022-08-15 contaminant removal moved to the peptide-to-protein projection section
    data <- data[!tmp$exclude, ]
    badness <- data ; for (colnm in semantics$sample) { ixs <- !is.na(badness[, colnm]) ; badness[ixs, colnm] <- badness$badness_score[ixs] }
    oufile=paste(dirOut, "images/badness_postExclusionList.png", sep='/') ; plot_density(badness, semantics, ymax=1.5, adjust=1., semacol='ligand', oufile=oufile)
    rm(oufile)
  }
  rm(tmp, badness, badness_max, badness_min, proteins, trim_percentile, colnm, colnms, ixs)
  #print(ggplot2::ggplot(data=tmp, ggplot2::aes(x=log(badness_score, base=10), y=resi)) + ggplot2::geom_point() + ggplot2::labs(x="badness score", y="max abs residual"))
}

#################################################################################
########## merge (or not) peptides that map to more than one protein ############
################ effectively, create "protein groups" again #####################
#################################################################################

if (FALSE) { ### sanity check on quants-per-peptide, assuming them to be exactly the same when the peptide is the same)
  #dataX <- aggregate(data[, semantics$sample], by=list(seq_mod=data$seq_mod), FUN="sd", na.rm=TRUE)
  dataX <- aggregate(data[, semantics$sample], by=list(sequence=data$sequence), FUN="sd", na.rm=TRUE)
  dataX <- dataX[rowSums(!is.na(dataX[, semantics$sample])) > 0, ]
  for (colnm in semantics$sample) { dataX[is.na(dataX[, colnm]), colnm] <- -999. }
  max(rje::rowMaxs(dataX[, semantics$sample]))
}

l_uniquePeptides <- TRUE ; if (l_uniquePeptides) {
  #source(paste0(dirCode, "025_findProteinGroups.R"))
  load(file=paste0(dirOut, '025_findProteinGroups.Rdata'))
}

#################################################################################
###### log-transform, filter, baseline-correct (or not), and center data ########
#################################################################################

l_baselineCorrect <- FALSE ; source(paste0(dirCode, "03_transformFilterBaselineCorrect.R"))
load(file=paste(dirOut, '03_transformFilterBaselineCorrect.Rdata', sep="/"))

if (TRUE) { ## Clean Up Sigma Plot for Supplemental Figures
  semanticsX <- semantics[semantics$cell_line == "CCR5", ]
  data.stats.list <- lplot(data = data.stats, semanticsX)
  g <- data.stats.list[[3]] ##Double check, ggplot objects really complicated to make non-index method of retrieval. Doable though
  g <- g + ggtitle(ggplot2::element_blank()) + labs(x = "Sample", y = "Log10 Quant Standard Deviation")
  ggsave(plot = g, filename = "~/Desktop/figure_images/231030_supFig2_SD.pdf", width = 9, height = 6)
}

#################################################################################
########################### remove batch effects ################################
#################################################################################

if (FALSE) { ### plot pairwise correlation heatmap, clustering dendrogram, and ARI pre-batch effect removal
  ggplot2::theme_update(text=ggplot2::element_text(size=16))
  semanticsX <- semantics[!(semantics$cell_line %in% c("bridge", "noH2O2")), ] ; semanticsX <- semanticsX[order(semanticsX$sample), ] ; dataX <- data[, semanticsX$sample]
  plotHeatmap(   dataX, semanticsX, oufile=paste0(dirOut, "/images/preBER_corrMat.png"))
  plotDendrogram(dataX, semanticsX, oufile=paste0(dirOut, "/images/preBER_dendrogram.pdf"))
  plotARI(       dataX, semanticsX, oufile=paste0(dirOut, "/images/preBER_ARI.png"))
  rm(dataX, semanticsX)
}

methodBER <- c("EB", "LM", "both")[3]
if (FALSE) {
  #source(paste0(dirCode, "04_removeBatchEffects.R"))
  load(file=paste0(dirOut, '/04_removeBatchEffects_', methodBER, '.Rdata'))
  if (TRUE) {plot_BERcvg(cvg, mode=2, oufile=paste0(dirOut, "/images/CCR5_BERcoverage.pdf"))}#, oufile = paste0(dirOut, "images/CCR5_BERcoverage.pdf"))}
  rm(cvg)
} else {
  #source(paste0(dirCode, "04_removeBatchEffects.R")); rm(cvg)
  load(file=paste0(dirOut, '/04_removeBatchEffects_', methodBER, '.Rdata'))
}


if (FALSE) { ### plot pairwise correlation heatmap, clustering dendrogram, and ARI post-batch effect removal
  ggplot2::theme_update(text=ggplot2::element_text(size=16))
  semanticsX <- semantics[!(semantics$cell_line %in% c("bridge", "noH2O2")), ] ; semanticsX <- semanticsX[order(semanticsX$sample), ] ; dataX <- data[, semanticsX$sample]
  plotHeatmap(   dataX, semanticsX, oufile=paste0(dirOut, "/images/postBER_corrMat_", methodBER, ".png"))
  plotDendrogram(dataX, semanticsX, oufile=paste0(dirOut, "/images/postBER_dendrogram_", methodBER, ".pdf"), semacolnms= c("cell_line","condition","plex"))
  plotARI(       dataX, semanticsX, oufile=paste0(dirOut, "/images/postBER_ARI_", methodBER, ".png"))
  rm(dataX, semanticsX)
}

if (FALSE) { ### plot pairwise correlation heatmap, clustering dendrogram, and ARI post-batch effect removal
  ggplot2::theme_update(text=ggplot2::element_text(size=16))
  semanticsX <- semantics[!(semantics$cell_line %in% c("bridge", "noH2O2")), ] ; semanticsX <- semanticsX[order(semanticsX$sample), ] ; dataX <- data[, semanticsX$sample]

  pdf(file = paste0(dirOut, "/images/postBER_corrMat_", methodBER, ".pdf"))
  plotHeatmap(   dataX, semanticsX, oufile="")
  dev.off()

  pdf(file = paste0(dirOut, "/images/postBER_ARI_", methodBER, ".pdf"))
  plotARI(       dataX, semanticsX, oufile="", semacolnms = c("cell_line", "condition", "plex"))
  dev.off()

  rm(dataX, semanticsX)
}

##############################################################################################
##### calculate high-average-residual peptides, potentially to be used as a filter above #####
################# added 2022-05-06 - trying to identify "outliers" ###########################
############ the idea is to identify peptides whose residuals are too high/too low ###########
##################### (>0.75 or <-0.75 according to my estimate - IK) ########################
########## then  delete them from either respective samples or from the entire set ###########
############### then go back to baseline correct and remove batch effects... #################
##############################################################################################

if (FALSE) { ### currently not used in any form
  resi <- calc_residuals(data, semantics) ; plot_density(resi, semantics, ymax=2.1, adjust=1., semacol='plex')
  tmp <- resi[, semantics$sample] ; tmp[is.na(tmp)] <-  999. ; resi$minResi <- rje::rowMins(tmp) ; resi$minResi[resi$minResi > 998.] <- NA
  resi$maxAbsResi <- -resi$minResi
  tmp <- resi[, c(semantics$sample, "maxAbsResi")] ; tmp[is.na(tmp)] <- -999. ; resi$maxResi <- rje::rowMaxs(tmp) ; resi$maxResi[resi$maxResi < -998.] <- NA
  resi$maxAbsResi <- resi$maxResi
  rm(tmp)
  resi <- resi[order(resi$maxAbsResi, decreasing=TRUE), ]
  #save(list=c("resi"), file=paste(dirOut, 'resi.Rdata', sep="/"))
}

#################################################################################
#################### project (or not) peptides onto proteins ####################
#################################################################################

### revised 2022-08-16
### 05_ now has options to project not only all peptides but also top X ordered by average quants
### by p-val across all conditions in semantics, or a random selection of peptides per protein
### capitalizing on the latter to build a tighter projection and also to calculate more statistically robust p-vals

if (FALSE){
  length(unique(data$AC))
  cat("Number of rows in data BEFORE contaminant filtering:", nrow(data))
  ### remove contaminants
  data_all <- data ;
  tmp <- data_all ; tmp$exclude <- 100*pracma::rand(1, nrow(tmp))[1, ] ; tmp$exclude <- tmp$exclude - tmp$perc_excluded ; tmp$exclude <- tmp$exclude < 0.
  data <- data_all[!tmp$exclude, ]# ; data <- data[(!grepl("^RPL", data$GN)) & (!grepl("^RPS", data$GN)), ] ### don't know what else to do with this ribosome information
  rownames(data) <- 1:nrow(data) ; rm(tmp)
  cat("Number of rows in data AFTER contaminant filtering:", nrow(data))
  length(unique(data$AC))
}


l_projectPeptidesOntoProteins <- TRUE

if (l_projectPeptidesOntoProteins) { ### really project

  ### first perform ntrials trial projections with 3 random peptides per protein in each, calculate the significance and average
  ntrials <- 5 ; prj_stats <- data.frame(AC=unique(data$AC), avg.p.pVal.adj=0, nof.non.NA=0)
  for (tri in 1:ntrials) {
    prjmode <- "3 random" ; source(paste0(dirCode, "05_projectPeptidesToProteins.R"))
    data_var <- data_wide_norm ; semantics_var <- semantics_wide ; source(paste0(dirCode, "06_ANOVA.R")) ;
    tmp <- plyr::join(prj_stats, data_var[, c("AC", "df.residual", "F", "p.pVal", "p.pVal.adj")], by="AC", type="left")
    ixs <- !is.na(tmp$p.pVal.adj) ; prj_stats$avg.p.pVal.adj[ixs] <- prj_stats$avg.p.pVal.adj[ixs] + tmp$p.pVal.adj[ixs] ; prj_stats$nof.non.NA[ixs] <- prj_stats$nof.non.NA[ixs] + 1
  }
  prj_stats$avg.p.pVal.adj <- prj_stats$avg.p.pVal.adj / (prj_stats$nof.non.NA+0.000001) ### average over trials

  ### next remove contaminants and project the maximum of 10 random peptides of what remains
  prjmode <- "top 10 by p.pVal.adj" ; source(paste0(dirCode, "05_projectPeptidesToProteins.R"))

  ### finally append prj_stats from the mini-trials
  tmp <- plyr::join(data_wide_norm, prj_stats, by="AC", type="left")
  data_wide_norm <- cbind(data.frame(ceil.p.pVal.adj=tmp$avg.p.pVal.adj), data_wide_norm)
  rm(tmp, prj_stats)

  save(list=c("data_wide_norm", "semantics_wide"), file=paste(dirOut, '05_projectPeptidesToProteins.Rdata', sep="/"))
  #load(file=paste(dirOut, '05_projectPeptidesToProteins.Rdata', sep="/"))

  rm(ANOVA_result, mo, ntrials, tri)

} else { ### run this if the analysis needs to be done at the level of peptides (not proteins) instead
  data_wide_norm <- data ; semantics_wide <- semantics
  rownames(data_wide_norm) <- 1:nrow(data_wide_norm)
}
data_wide_norm <- data_wide_norm[(!grepl("^RPL", data_wide_norm$GN)) & (!grepl("^RPS", data_wide_norm$GN)), ] ### don't know what else to do with this ribosome information


#################################################################################
########################## run answer questions script ##########################
#################################################################################

