library(pracma) ### for randi
library(limma)
set.seed(1)

if (!exists("methodBER")){methodBER <- "both"}

if (FALSE) {
  semanticsX <- semantics[semantics$cell_line %in% c("Cyto", "PM", "Endo"), ] ; dataX <- data[, semanticsX$sample] ; plotHeatMapHeavy(dataX, semanticsX, "compartment markers")
  semanticsX <- semantics[semantics$ligand=="5P14", ] ; dataX <- data[, semanticsX$sample] ; plotHeatMapHeavy(dataX, semanticsX, "ligand=5P14")
  semanticsX <- semantics[semantics$ligand=="6P4", ] ; dataX <- data[, semanticsX$sample] ; plotHeatMapHeavy(dataX, semanticsX, "ligand=6P4")
  semanticsX <- semantics[semantics$ligand=="RANTES", ] ; dataX <- data[, semanticsX$sample] ; plotHeatMapHeavy(dataX, semanticsX, "ligand=RANTES")
  semanticsX <- semantics[grepl("bridge", semantics$cell_line), ] ; dataX <- data[, semanticsX$sample] ; plotHeatMapHeavy(dataX, semanticsX, "bridges")
  rm(dataX, semanticsX)
}

####### Correct data for batch (plex) effects using experimental bridge ################
data1 <- data #Preservation of input df, used later for method analysis
cat(paste0("\n", sum(rowSums(!is.na(data[, semantics$sample]))), " non-NA values pre experimental-bridge-based BER")) #Number of values
cat(paste0("\nBatch effect removal using experimental bridges..."))

plexes <- unique(semantics$plex) ; plexes <- plexes[order(plexes)] ; bridges <- data.frame(AC=data$AC, has_data=TRUE) ; #initializes a new dataframe with a row for each AC. Also gathers all available plexes
bridges$has_data[rowSums(!is.na(data[, semantics$sample[semantics$cell_line!="bridge"]])) == 0] <- FALSE #If the AC has no data, then it is marked
bridges$EB_resolvable <- TRUE ### experimental-bridge-resolvable, initializes new column
for (p in plexes) {
  colnms <- semantics$sample[semantics$plex==p & semantics$cell_line=="bridge"] #Gathers the names of all bridge samples from the particular plex. Usually just 1 per plex
  bridges <- cbind(bridges, data.frame( rowMeans( data[, colnms, drop=FALSE], na.rm=TRUE)) ) #Find the average bridge value of a protein for the particular plex
}
colnms <- paste0("bridge", plexes) ; colnames(bridges) <- c("AC", "has_data", "EB_resolvable", colnms) #Rename bridge columns
bridge <- rowMeans(bridges[, colnms, drop=FALSE], na.rm=TRUE) ; bridge[is.na(bridge)] <- 0. #Gets the average of all bridge values for a protein, then replaces NAs with 0

## Apply the bridge-based batch effect removal
dataEB <- data
for (p in plexes) {
  bcolnm <- paste0("bridge", p) ; bridges[,bcolnm] <- bridges[,bcolnm] - bridge; bridges[is.na(bridges[,bcolnm]),bcolnm] <- 0 #Mean centers bridge, unless it is na then it remains 0
  nofNonNAwBridge <- sum(rowSums(!is.na(   data[ semantics$sample[semantics$plex==p & semantics$cell_line!="bridge"]]  )) > 0 & rowSums(!is.na(data[semantics$sample[semantics$plex==p & semantics$cell_line=="bridge"]])) > 0) #Must have normal data and bridge data
  nofNonNA <- sum(rowSums(!is.na(data[semantics$sample[semantics$plex==p & semantics$cell_line!="bridge"]])) > 0) #Only non bridge samples
  #cat(paste0("\nplex ", p, ": ", nofNonNAwBridge, " of ", nofNonNA, " peptides have bridge (of ", nrow(data), " total)"))
  colnms <- semantics$sample[semantics$plex==p] ; for (colnm in colnms) { dataEB[,colnm] <- data[,colnm] - bridges[,bcolnm] } #Iterates through each sample from the plex, subtracting the bridge value from the data value
  bridges$EB_resolvable[rowSums(!is.na(dataEB[semantics$sample[semantics$plex==p & semantics$cell_line!="bridge"]])) > 0 & rowSums(!is.na(dataEB[semantics$sample[semantics$plex==p & semantics$cell_line=="bridge"]])) == 0] <- FALSE
}
cat(paste0("\n", sum(rowSums(!is.na(dataEB[, semantics$sample]))), " non-NA values post experimental-bridge-based BER"))
cat(paste0("\n", sum(bridges$has_data &  bridges$EB_resolvable), " rows have data and are experimental-bridge-resolvable, ",
             sum(bridges$has_data & !bridges$EB_resolvable), " have data but are NOT experimental-bridge-resolvable"))
rm(bridge, colnm, bcolnm, colnms, p, plexes, nofNonNA, nofNonNAwBridge) ; data2 <- data
if (methodBER != "LM") {data <- dataEB}; rm(dataEB)

if (FALSE) {
  cat("\n", methodBER, " method(s) used for batch effect removal.")
  semanticsX <- semantics[semantics$cell_line %in% c("Cyto", "PM", "Endo"), ] ; dataX <- data[, semanticsX$sample] ; plotHeatMapHeavy(dataX, semanticsX, "compartment markers")
  semanticsX <- semantics[semantics$ligand=="5P14", ] ; dataX <- data[, semanticsX$sample] ; plotHeatMapHeavy(dataX, semanticsX, "ligand=5P14")
  semanticsX <- semantics[semantics$ligand=="6P4", ] ; dataX <- data[, semanticsX$sample] ; plotHeatMapHeavy(dataX, semanticsX, "ligand=6P4")
  semanticsX <- semantics[semantics$ligand=="RANTES", ] ; dataX <- data[, semanticsX$sample] ; plotHeatMapHeavy(dataX, semanticsX, "ligand=RANTES")
  semanticsX <- semantics[grepl("bridge", semantics$cell_line), ] ; dataX <- data[, semanticsX$sample] ; plotHeatMapHeavy(dataX, semanticsX, "bridges")
  rm(dataX, semanticsX)
}

#################### LINEAR MODELING PART ############################

### delete experimental bridges from data (or not (current as of 221229))
#colnms <- colnames(data)[!(colnames(data) %in% semantics$sample)] ; semantics <- semantics[semantics$cell_line!="bridge", ] ; data <- data[, c(colnms, semantics$sample)]

### figure out the coverage of different peptides by plexes and conditions, plot while we are at it
cvg_plex <- plot_coverage(data = data, semantics = semantics, semacol = 'plex',      oufile = paste(dirOut, "images", paste0(kwd, '_coverage_plex_preBER.png'), sep='/'))
cvg_cond <- plot_coverage(data = data, semantics = semantics, semacol = 'condition', oufile = paste(dirOut, "images", paste0(kwd, '_coverage_condition_preBER.png'), sep='/'))

### build a data frame called cvg for figuring out suitable base plexes and base conditions for every row
cvg <- data.frame(has_data=bridges$has_data, EB_resolvable=bridges$EB_resolvable, LM_resolvable=TRUE, nof_plex=cvg_plex$nof, nof_cond=cvg_cond$nof) #Pulls values calculated earlier, from the "bridges" df
cvg$base_plex <- NA ; cvg$base_cond <- NA ; cvg$plex_max <- 0 ; cvg$cond_max <- 0 #More initialization of columns
cvg <- cbind(cvg, cvg_plex[,colnames(cvg_plex) != "nof"]) ; cvg$plex_max <- apply(cvg, 1, function(x) as.numeric(max(x[as.character(unique(semantics$plex))]))) #Finds the plex that has the most values for a given protein
cvg <- cbind(cvg, cvg_cond[,colnames(cvg_cond) != "nof"]) ; cvg$cond_max <- apply(cvg, 1, function(x) as.numeric(max(x[as.character(unique(semantics$cond))]))) #Finds the condition that has the most values for a given protein
rm(cvg_cond, cvg_plex, bridges) #all of these dataframes have been aggregated in one way or another into "cvg"

if (methodBER != "EB") {cat("\nLinear modeling method starting.")

  ### add and fill in columns for the different condition/plex pairs and whether they can be used as base for LM
  cond_x_plex <- expand.grid(unique(semantics$condition), unique(semantics$plex)) ; colnms <- paste(cond_x_plex$Var1, cond_x_plex$Var2, sep="@") ; rm(cond_x_plex) #Makes vector of strings with all combinations of cond and plex
  for (colnm in colnms) {
    cvg[,colnm] = 0 #init new column for cond/plex combo
    cond_ <- unlist(strsplit(colnm, "@"))[1] ; plex_ <- as.character(unlist(strsplit(colnm, "@"))[2]) #Method for getting current plex and cond that are being dealt with
    cvg[(cvg[, cond_] > 0) & (cvg[, plex_] > 0) & (cvg[, plex_] > cvg[,'plex_max']-3), colnm] <- 1 #If prot has data in plex & cond, and is within 3 difference of plex max, mark with 1
  }

  ### find a close-to-minimal subset of condition@plex pairs that covers all rows
  cat("\n\nSearching for a close-to-minimal subset of condition@plex pairs that covers all rows:")
  cvg_orig <- cvg ; if (exists('subcolnms')) { rm(subcolnms) }
  #The colnms here are still the cond plex combos
  repeat {
    tmp <- cbind(as.data.frame(colnms), colSums(cvg[,colnms])) #Column 1 is are the combos, column 2 are the number of proteins the combos are good for. Rownames are also the combos
    tmp <- tmp[(tmp[,2] == max(tmp[,2])), ] #Finds column that has the most coverage
    colnm <- tmp[1,1] ; if (tmp[1,2] == 0) { break } #Stores current combo name; if the last remaining combo is 0 max, then break the loop
    cat(paste("\n", colnm, tmp[1,2]))
    cvg[cvg[, colnm] == 1, colnms] <- 0 #Back in the original 'cvg' df, look up the proteins that have been marked by this particular combo. Remove this mark. This also in effect removes the combo from the loop
    if (!exists('subcolnms')) { subcolnms <- colnm } else { subcolnms <- c(subcolnms, colnm) } #Adds the combo that has just been scrubbed clean of marks in 'cvg' to a list
  }
  cvg <- cvg_orig; cvg1 <- cvg[, !(colnames(cvg) %in% colnms)]; cvg2 <- cvg[, subcolnms] ; cvg <- cbind(cvg1, cvg2) ; rm(cvg1, cvg2, cvg_orig, tmp) #Preserves only the combo columns that were collected in the previous loop

  ### find a semi-balanced condition@plex partition consisting of the same pairs as the coverage above
  cat("\n"); ndata <- nrow(cvg) ; pb <- txtProgressBar(min=1, max=ndata, style=3)
  for (i in 1:ndata) {
    nofOnes <- sum(cvg[i,subcolnms]) ; if (nofOnes > 0) { #If there is coverage (Maybe can be > 1?)
      ix <- pracma::randi(nofOnes)[1] ### will keep this 1 in the current row, replace the rest with 0's; randomly selects one of the available combos
      colnm <- subcolnms[cvg[i,subcolnms] == 1] [ix] #colnm is the name of the combo column randomly selected
      cvg[i,subcolnms] <- 0 ; cvg[i,colnm] <- 1 #Removes all marks except for one
      cvg[i, 'base_cond'] <- unlist(strsplit(colnm, "@"))[1] ; cvg[i, 'base_plex'] <- as.character(unlist(strsplit(colnm, "@"))[2]) #Labels the base cond/plex from selected combo
    }
    setTxtProgressBar(pb, i)
  }
  close(pb) ; rm(pb, ndata, nofOnes, ix)

  ### check that the partition is a partition (pairwise intersections are empty) and that its total covers all rows
  cat(paste('\nThe identified partition covers', sum(cvg[,subcolnms]), 'rows with the following intersections:'))
  for (colnm in subcolnms) { cat("\n"); cat(colSums(cvg[, subcolnms] * cvg[, colnm, drop=TRUE])) } ; rm(subcolnms, colnms, colnm)

  ### calculate LM "resolvability" using the identified base plexes and base conditions
  cat("\n\nCalculating LM-resolvability of rows:\n")
  cvg$LM_resolvable[cvg$nof_cond < 2] <- FALSE # rows with less than 2 conditions are not resolvable
  # cvg$LM_resolvable[cvg$plex_max == 2] <- FALSE # rows with 2 or less points in the maximally populated plex are not resolvable?
  cvg$LM_resolvable[(rowSums(!is.na(data)) < 4) & (cvg$nof_plex > 1)] <- FALSE # rows with less than 4 non-NA values in more than 1 plex are not resolvable (I think this non-NA call is missing a semantics$sample subsetter?)
  ndata <- nrow(data) ; pb <- txtProgressBar(min=1, max=ndata, style=3)
  for (ix in 1:ndata) {
    if (cvg$LM_resolvable[ix]) { #LM_resolvable is initialized as TRUE for all rows
      res <- testPoint(data, semantics, ix, cvg$base_plex[ix], cvg$base_cond[ix]) ; cvg$LM_resolvable[ix] <- res$resolvable
      if (res$resolvable) { cvg$base_plex[ix] <- res$base_plex ; cvg$base_cond[ix] <- res$base_cond }
    }
    setTxtProgressBar(pb, ix)
  }
  close(pb) ; rm(pb, res, ndata)

  ### Build linear model; correct data for batch effects using the obtained coefficients
  ### https://genomicsclass.github.io/book/pages/expressing_design_formula.html
  ### https://www.rdocumentation.org/packages/limma/versions/3.28.14/topics/06.LinearModels
  ### https://www.rdocumentation.org/packages/limma/versions/3.28.14
  ### http://www.bioconductor.org/packages/devel/bioc/vignettes/limma/inst/doc/usersguide.pdf

  cat("\n\nLM-based batch-effect removal:")
  cvg_cnt <- plyr::count(cvg[cvg$LM_resolvable, c('base_plex','base_cond')]) #Count for each condition and plex combo, how many prot are resolvable
  cvg_cnt <- cvg_cnt[order(cvg_cnt$freq, decreasing=TRUE),] #Move maximum coverage to top
  cvg_cnt <- cvg_cnt[!is.na(cvg_cnt$base_plex) & !is.na(cvg_cnt$base_cond),] #Remove any NAs, nonresolvable?
  dataNormRBE <- data #copy the data before any BER is preformed
  for (i in 1:nrow(cvg_cnt)) {
    plex_ <- relevel(factor(semantics$plex), cvg_cnt$base_plex[i]) #Moves target plex to top of list
    cond_ <- relevel(factor(semantics[,'condition']), cvg_cnt$base_cond[i]) #Moves target condition to top of list
    ixsX  <- cvg$LM_resolvable & (cvg$base_plex == cvg_cnt$base_plex[i]) & (cvg$base_cond == cvg_cnt$base_cond[i]) #Gets ixs of proteins that match the target plex & condition
    cat(c("\n", i, sum(ixsX)))

    ### extract dataX and normalize using removeBatchEffect
    dataX <- data[ixsX, semantics$sample, drop=FALSE]
    design <- stats::model.matrix(~0+cond_)
    dataX <- as.data.frame(removeBatchEffect(dataX, batch=plex_, design=design))
    dataNormRBE[ixsX, semantics$sample] <- dataX
  }
  data <- dataNormRBE
  rm(dataX, dataNormRBE, ixsX, cond_, plex_, design) #; rm(cvg_cnt, cvg)
}

cat(paste0("\n\n", sum(rowSums(!is.na(data[, semantics$sample]))), " non-NA values post linear-modeling-based BER"))

cat(paste0("\n", sum( cvg$has_data &  cvg$EB_resolvable &  cvg$LM_resolvable), " rows have data and are resolvable by both methods"))
cat(paste0("\n", sum( cvg$has_data &  cvg$EB_resolvable & !cvg$LM_resolvable), " rows have data and are resolvable only by experimental bridge"))
cat(paste0("\n", sum( cvg$has_data & !cvg$EB_resolvable &  cvg$LM_resolvable), " rows have data and are resolvable only by linear modeling"))
cat(paste0("\n", sum( cvg$has_data & (cvg$EB_resolvable | cvg$LM_resolvable)), " rows have data and are resolvable somehow"))
cat(paste0("\n", sum( cvg$has_data & !cvg$EB_resolvable & !cvg$LM_resolvable), " rows have data and are not resolvable at all"))
cat(paste0("\n", sum(!cvg$has_data), " rows have no data"))
if (methodBER == "both"){cat("\nboth methods used for batch effect removal.\n")} else {cat("\n", methodBER, "method used for batch effect removal.\n")}

#tmp <- data[cvg$has_data &  cvg$EB_resolvable &  cvg$LM_resolvable, ]

if (FALSE) {
  semanticsX <- semantics[semantics$cell_line %in% c("Cyto", "PM", "Endo"), ] ; dataX <- data[, semanticsX$sample] ; plotHeatMapHeavy(dataX, semanticsX, "compartment markers")
  semanticsX <- semantics[semantics$ligand=="5P14", ] ; dataX <- data[, semanticsX$sample] ; plotHeatMapHeavy(dataX, semanticsX, "ligand=5P14")
  semanticsX <- semantics[semantics$ligand=="6P4", ] ; dataX <- data[, semanticsX$sample] ; plotHeatMapHeavy(dataX, semanticsX, "ligand=6P4")
  semanticsX <- semantics[semantics$ligand=="RANTES", ] ; dataX <- data[, semanticsX$sample] ; plotHeatMapHeavy(dataX, semanticsX, "ligand=RANTES")
  semanticsX <- semantics[grepl("bridge", semantics$cell_line), ] ; dataX <- data[, semanticsX$sample] ; plotHeatMapHeavy(dataX, semanticsX, "bridges")
  rm(dataX, semanticsX)
}

if (FALSE) { ### sanity check for the batch effect removal routines: the within-plexes relative quants should have stayed the same
  plexes <- unique(semantics$plex) ; plexes <- plexes[order(plexes)] ; ndata <- nrow(data)
  pb <- txtProgressBar(min=1, max=ndata, style=3)
  for (i in 1:ndata) {
    for (p in plexes) {
      colnms <- semantics$sample[semantics$plex==p]
      v1 <- data1[i, colnms] ; v <- data[i, colnms] ; d <- v-v1 ; d <-d[!is.na(d)]
      if (sum(is.na(v1)) - sum(is.na(v)) != 0) { cat(paste0("\nplex ", p, ", peptide ", i)) }
      if (length(d) > 0) { if (max(d)-min(d)>0.00001) { cat(paste0("\nplex ", p, ", peptide ", i, ": d = <", paste0(d, collapse=" "), ">")) } }
    }
    setTxtProgressBar(pb, i)
  }
  close(pb) ; rm(pb, plexes, ndata, colnms, v, v1, p, d, i)
}
rm(data1, data2)

### here delete experimental bridges from data for sure
colnms <- colnames(data)[!(colnames(data) %in% semantics$sample)] ; semantics <- semantics[!grepl("bridge", semantics$cell_line), ] ; data <- data[, c(colnms, semantics$sample)]

data[cvg$has_data & !cvg$EB_resolvable & !cvg$LM_resolvable, semantics$sample] <- NA ### clear out the data that is not alignable
rm(cvg_cnt)

cvg_plex <- plot_coverage(data = data, semantics = semantics, semacol = 'plex',      oufile = paste(dirOut, "images", paste0(kwd, '_coverage_plex_postBER.png'), sep='/'))
cvg_cond <- plot_coverage(data = data, semantics = semantics, semacol = 'condition', oufile = paste(dirOut, "images", paste0(kwd, '_coverage_condition_postBER.png'), sep='/'))
rm(cvg_plex, cvg_cond)

if (exists("methodBER")) {
  if (exists("extraExt")) {
    save(list=c("data", "semantics", "cvg"), file=paste0(dirOut, '/04_removeBatchEffects_', methodBER , '_', extraExt,'.Rdata'))
  } else {
    save(list=c("data", "semantics", "cvg"), file=paste0(dirOut, '/04_removeBatchEffects_', methodBER ,'.Rdata'))
  }
}else{
  save(list=c("data", "semantics", "cvg"), file=paste0(dirOut, '/04_removeBatchEffects.Rdata'))
}
