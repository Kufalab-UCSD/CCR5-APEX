#library(limma)

### 2022-08-15: introduced options for peptide-to-protein projection:
# * all
# * top # ordered by average quants
# * top # ordered by p-val across all conditions in semantics
# * top # ordered by any other column ? - not yet implemented
# * random #

#prjmode <- "all"
#prjmode <- "top 3 by avg_quant"
#prjmode <- "top 3 by p.pVal.adj"
#prjmode <- "3 random"

if (!exists("prjmode")) { prjmode <- "all" } ### "all" is the default
if (!exists("verbose")) { verbose <- TRUE }

prjmodes <- c("all", "top [0-9]+ by avg_quant", "top [0-9]+ by p.pVal.adj", "[0-9]+ random", "invalid")
for (mo in prjmodes) { if (grepl(mo, prjmode)) { break } }
if (mo == "invalid") { stop(paste0("Error: invalid mode ", prjmode)) }

maxpepcnt <- as.numeric(stringr::str_extract(prjmode, "[0-9]+"))
if (verbose) cat(paste0("\n  Projection mode: ", prjmode, " (", mo, "); max ", maxpepcnt, " peptides per protein\n"))

########### remove empty rows, sort according to the projection mode (prjmode), calculate stats of peptides per protein #############
data <- data[rowSums(!is.na(data[, semantics$sample])) > 1, ] ### remove entries with less than 2 quants

### may need: avg_quant, p.pVal.adj, completeness
if ((mo == "all") | (mo == "top [0-9]+ by avg_quant")) {
  data$avg_quant <- rowMeans(data[, semantics$sample], na.rm=TRUE)
  median_quant <- stats::quantile(data$avg_quant, probs=0.68, na.rm=TRUE) ; data$avg_quant[data$avg_quant > median_quant] <- median_quant ### trim at 1 sigma
  colnm <- "avg_quant"
}
if (mo == "top [0-9]+ by p.pVal.adj") {
  if (!("p.pVal.adj" %in% colnames(data))) { data_var <- data ; semantics_var <- semantics ; source(paste0(dirCode, "06_ANOVA.R")) ; data$p.pVal.adj <- data_var$p.pVal.adj ; rm(data_var, semantics_var) }
  data$p.pVal.adj[data$p.pVal.adj > 3.] <- 3. ### trim at 3.
  colnm <- "p.pVal.adj"
}
if (mo == "[0-9]+ random") { data$rnd <- base::sample(1:nrow(data)) ; data$rnd <- floor(data$rnd / length(unique(semantics$sample))) ; colnm <- "rnd" }
#used to be: if (mo == "[0-9]+ random") { data$rnd <- base::sample(1:nrow(data)) ; data$rnd <- floor(data$rnd / 10.) ; colnm <- "rnd" }

data$nof_quants <- rowSums(!is.na(data[, semantics$sample]))
data$nof_conds <- 0 ; for (cond in unique(semantics$condition)) { data$nof_conds <- data$nof_conds + as.numeric(rowSums(!is.na(data[ , semantics$sample[semantics$condition == cond], drop=FALSE])) > 0) }

### order data by colnm while prioritizing most complete peptides
### NEED TO ADD: all other conditions equal, prioritize those peptides that have more "friends" (same protein, same conditions)
data <- data[order(data[, colnm], decreasing=TRUE), ]
data <- data[order(data$nof_quants, decreasing=TRUE), ]
data <- data[order(data$nof_conds, decreasing=TRUE), ]
rownames(data) <- 1:nrow(data) ; data[!duplicated(data[, colnm]), colnm] <- max(data[, colnm])+1
data <- data[order(data[, colnm], decreasing=TRUE), ]
data <- data[, !(colnames(data) %in% c(colnm, "avg_quant", "p.pVal.adj", "nof_quants", "nof_conds", "rnd"))]
rownames(data) <- 1:nrow(data)

if (!grepl("all", prjmode)) { ### remove everything except the first maxpepcnt rows for each AC
  toProject <- data[!duplicated(data$AC), ] ; tmp <- data[duplicated(data$AC), ]
  if (maxpepcnt > 1) { for (i in 2:maxpepcnt) { toProject <- rbind(toProject, tmp[!duplicated(tmp$AC), ]) ; tmp <- tmp[duplicated(tmp$AC), ] }}
  rm(tmp)
} else {
  toProject <- data
}

prot_cnt <- plyr::count(toProject[, 'AC']) ; colnames(prot_cnt)[1] <- "AC" ; prot_cnt <- prot_cnt[!is.na(prot_cnt$AC), ]
toProject$pepcnt <- 1
for (i in 1:nrow(prot_cnt)) { if (prot_cnt$freq[i] > 1) { ixs <- toProject$AC == prot_cnt$AC[i] ; toProject$pepcnt[ixs] <- 1:sum(ixs) ; } } ### append peptide counts per protein
rm(prot_cnt)

#################### PROJECT PEPTIDES ONTO PROTEINS #####################

n.char <- ceiling(log10(max(toProject$pepcnt)+1)) ### an accessory variable for forming peptide-dependent column names with "001" or similar

colnms <- colnames(toProject)[!(colnames(toProject) %in% semantics$sample)] ; toProject <- toProject[, c(colnms, semantics$sample)]
semantics_wide <- semantics ; peplex <- paste0('pep', substr(as.character(10^n.char+1), 2, 1+n.char)) ; semantics_wide$plex <- peplex ; semantics_wide$sample <- paste0(peplex, "_", semantics$sample)
data_wide <- toProject[toProject$pepcnt == 1, ] ; colnames(data_wide) <- c(colnms, semantics_wide$sample)

prot_cnt <- plyr::count(toProject[, 'pepcnt']) ; colnames(prot_cnt)[1] <- 'pepcnt'
# build a bar graph here
rm(prot_cnt)

if (maxpepcnt>1) {
  for (i in 2:max(toProject$pepcnt)) {
    if (sum(toProject$pepcnt == i) > 0) {
      tmps <- semantics ; peplex <- paste0('pep', substr(as.character(10^n.char+i), 2, 1+n.char)) ; tmps$plex <- peplex
      tmps$sample <- paste0(peplex, "_", semantics$sample)
      semantics_wide <- rbind(semantics_wide, tmps)
      tmp <- toProject[toProject$pepcnt == i, c("AC", "sequence", semantics$sample)]
      colnames(tmp) <- c('AC', "tmpseq", tmps$sample) ; data_wide <- plyr::join(data_wide, tmp, by='AC', type="left")
      ixs <- !is.na(data_wide$tmpseq) ; data_wide$pepcnt[ixs] <- data_wide$pepcnt[ixs]+1
      data_wide$sequence[ixs] <- paste(data_wide$sequence[ixs], data_wide$tmpseq[ixs], sep=", ") ; data_wide <- within(data_wide, rm(tmpseq))
    }
  }
  rm(tmp, tmps, ixs)
}
rm(peplex)

colnms <- colnames(data_wide)[!(colnames(data_wide) %in% semantics_wide$sample)] ; AC_pepcnt <- data_wide[, colnms] ;

### define a mini-function for getting the completeness and resolvability stats on a dataset
### this is redundant with ProteomicsToolkit::testPoint() but works much faster
resolvable <- function(data, semantics, design.mat=NULL) { ### returns a list with n.non.NA, n.rows.w.data, n.rows.resolvable, and rows.resolvable (the latter a list of logicals)
  if (is.null(design.mat)) { ### build a list of design matrices ~cond_+plex_ with all possible conds as base condition
    plex_ <- factor(semantics$plex) ; design.mat <- list() ; conds <- unique(semantics$condition) ; nconds <- length(conds)
    for (i in 1:nconds) { cond_ <- factor(semantics$condition, levels=conds) ; design.mat[[conds[1]]] <- model.matrix(~cond_+plex_) ; conds <- c(conds[2:nconds], conds[1]) }
  }
  if (is(design.mat, "matrix")) { d.m <- design.mat ; design.mat <- list() ; design.mat[[1]] <- d.m } ### if a single matrix, make it into a 1-element list for the loop below
  n <- nrow(data) ; result <- list(n.non.NA=NA, n.rows.w.data=NA, n.rows.resolvable=NA, rows.resolvable=rep(FALSE, n)) ; names(result$rows.resolvable) <- rownames(data)
  defined <- !is.na(data[, semantics$sample])
  for (d.m in design.mat) {
    for (i in 1:n) {
      if (result$rows.resolvable[i]) next ;
      if (sum(defined[i, ])<1) next ; d <- d.m[unlist(defined[i, ]), ]
      if (sum(colSums(d)>0)<1) next ; d <- d[, colSums(d)>0, drop=FALSE]
      result$rows.resolvable[i] <- (Matrix::rankMatrix(d)==ncol(d))
    }
  }
  result$n.rows.resolvable <- sum(result$rows.resolvable) ; result$n.non.NA <- sum(defined) ; result$n.rows.w.data <- sum(rowSums(defined)>0) ; return(result)
} ### end definition for the resolvable(data, semantics) function

### make sure all data is resolvable
if (maxpepcnt == 1) { data_wide <- data_wide[, semantics_wide$sample] ; data_wide_norm <- data_wide[, semantics_wide$sample] } else {
  reso <- resolvable(data_wide, semantics_wide)
  if (verbose) cat(paste0("\n", reso$n.non.NA, " non-NA value(s), ", reso$n.rows.w.data, " protein(s) with data, and ", reso$n.rows.resolvable, " resolvable protein(s) before the first chew-back cycle"))
  while (TRUE) {
    pepcnt.new <- AC_pepcnt$pepcnt
    for (p in unique(AC_pepcnt$pepcnt)) {
      peplex <- paste0('pep', substr(as.character(10^n.char+p), 2, 1+n.char))
      ixs <- AC_pepcnt$pepcnt==p & !reso$rows.resolvable ### find all unresolvable proteins with exactly p peptides
      data_wide[ixs, semantics_wide$sample[semantics_wide$plex==peplex]] <- NA ### "chew back", i.e. remove the last peptide in these proteins
      pepcnt.new[ixs] <- AC_pepcnt$pepcnt[ixs]-1 ### decrease the peptide count for the chewed-back proteins
    }
    AC_pepcnt$pepcnt <- pepcnt.new
    reso <- resolvable(data_wide, semantics_wide)
    if (verbose) cat(paste0("\n", reso$n.non.NA, " non-NA value(s), ", reso$n.rows.w.data, " protein(s) with data, and ", reso$n.rows.resolvable, " resolvable protein(s) after the current chew-back cycle"))
    if (sum(reso$n.rows.resolvable)==reso$n.rows.w.data) break
  }
  if (verbose) cat("\n\n")

  #################################################################################################
  ####### Correct data for "batch" (peptide) effects using linear modeling (limma) ################
  #################################################################################################

  ### remove the non-data columns (they are saved in AC_pepcnt); order the data columns as in semantics_wide
  data_wide <- data_wide[, semantics_wide$sample]
  ### build a list of design matrices ~cond_+_plex with all possible conds as base condition
  plex_ <- factor(semantics_wide$plex) ; design.mat <- list() ; conds <- unique(semantics_wide$condition) ; nconds <- length(conds)
  for (i in 1:nconds) { cond_ <- factor(semantics_wide$condition, levels=conds) ; design.mat[[conds[1]]] <- model.matrix(~cond_+plex_) ; conds <- c(conds[2:nconds], conds[1]) }

  data_wide_norm <- data_wide ; data_wide_norm[!reso$rows.resolvable, semantics_wide$sample] <- NA
  conds <- unique(semantics_wide$condition) ; nconds <- length(conds)
  for (i in 1:nconds) {
    reso <- resolvable(data_wide, semantics_wide, design.mat=design.mat[[conds[i]]]) ; ixs <- reso$rows.resolvable & AC_pepcnt$pepcnt > 1
    cond_ <- factor(semantics_wide$condition, levels=conds)
    data_wide_norm[ixs, semantics_wide$sample] <- as.data.frame(limma::removeBatchEffect(data_wide_norm[ixs, semantics_wide$sample], batch=plex_, design=model.matrix(~0+cond_)))
    conds <- c(conds[2:nconds], conds[1]) ### circular permutation
  }
  rm(plex_, cond_, cond, conds, i, reso, nconds, design.mat, pepcnt.new, peplex, p, ixs)

  # ### test and remove those quants that could not be aligned with RBE
  # ### 2023-04-08 this is incorrect and is not working, replaced by the above
  # res <- limma::lmFit(data_wide[ , semantics_wide$sample], model.matrix(~cond_+plex_)) ; res <- limma::eBayes(res)
  # coeff.col <- colnames(res$design)[grepl("plex_|Intercept", colnames(res$design))]
  # res <- limma::topTable(res, adjust="BH", coef=coeff.col, number=nrow(data_wide), sort.by="none")
  # l_done <- TRUE ; peplexes <- unique(semantics_wide$plex)
  # for (p in 2:length(peplexes)) {
  #   ixs <- is.na(res[, paste0("plex_", peplexes[p])]) ; colnms <- semantics_wide$sample[semantics_wide$plex==peplexes[p]]
  #   ixs <- ixs & (rowSums(!is.na(data_wide_norm[, colnms])) > 0) ; if (sum(ixs) < 1) next
  #   l_done <- FALSE
  #   if (verbose) cat(paste0("\n  * ", peplexes[p], ": clearing ", sum(!is.na(data_wide_norm[ixs, colnms])), " quant(s) for ", sum(ixs), " protein(s) because batch effects are not estimable for them"))
  #   data_wide_norm[ixs, colnms] <- NA ### added 2023-04-08; by some reason was not here before although this is the actual *clearing* command
  # }
  # if (verbose) cat("\n\n")
  # rm(ixs, coeff.col, res, colnms)
} ### endif (maxpepcnt > 1)

### reattach non-quant columns
data_wide <- cbind(AC_pepcnt, data_wide) ; data_wide_norm <- cbind(AC_pepcnt, data_wide_norm)
rm(AC_pepcnt, toProject, verbose)
rm(prjmode, prjmodes, maxpepcnt, colnm)

### no longer saves, save in the calling script instead
#save(list=c("data_wide", "data_wide_norm", "semantics_wide"), file=paste(dirOut, '05_projectPeptidesToProteins.Rdata', sep="/"))

