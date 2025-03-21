#########################################################################################
####### Calculate statistical distances between proteins in data_XY, ####################
####### assuming variation across conditions in semantics_XY ############################
#########################################################################################

# spline_df = 0
# if (grepl("CCR2_phospho", dirCode)) {
#   step06_filepath <- paste0(strsplit(dirCode, "//CCR2_phospho/")[[1]], "//CCR5_APEX/06_ANOVA_multi.R")
#   source(step06_filepath)
#   if (Q == "Q_R1R2_no720") { spline_df = NULL }
# } else { source(paste0(dirCode, "06_ANOVA_multi.R")) } ### define the function
if (file.exists(paste0(dirCode, "06_ANOVA_multi.R"))){source(paste0(dirCode, "06_ANOVA_multi.R"))
} else if (!exists("step06_ANOVA", envir = .GlobalEnv)){stop("Cannot find step06_ANOVA function.")}

#########################################################################################
### Inputs:
### * data_XY - a df of *proteins* (not protein pairs), semantics_XY
### * model.on, stats.on - how to assess the significance of these proteins
### * splineDF - an initial/default number degrees of freedom for spline modeling of numeric columns:
###   - if 0, use knots at experimentally measured points
###   - if -1, place knots between exp points, except for the boundaries
###   - if NULL, do not use splines at all but will instead treat numerical variables as categorical
### * these will be arguments for the call of step06_ANOVA(data_XY, semantics_XY, model.on=model.on, stats.on=stats.on, verbose=FALSE, splineDF=splineDF)
#########################################################################################
### Result: protein_pairs, a dataframe with the following columns:
### * ACx and ACy for the protein pair in question
### * p.pVal.x and p.pVal.y, the significances of the individual proteins in the pair in accordance with model.on and stats.on
### * x.or.y: max of p.pVal.x and p.pVal.y
### * x.and.y: collated significance, also in accordance with model.on and stats.on
### * x.vs.y.11: x.or.y-x.and.y (script 11 philosophy)
### * x.vs.y.int: significance of the interaction coefficients
### * x.vs.y.F: significance of the Prism-style F-test
### * x.vs.y.AIC: significance of the Prism-style Akaike Information Criterion (AIC)
#########################################################################################
### All significances are on -log10 scale, larger number indicates bigger difference
#########################################################################################

step07_distanceMatrix <- function(data_XY, semantics_XY, model.on="condition", stats.on="condition", splineDF=0, nproc=1) {

  ### accessory function to perform calculations on chunk #chunkix
  ### this is necessary for the parallel implementation
  step07_distanceMatrix_onChunk <- function(chunkix) {
    if (chunkix > nrow(chunks)) return()

    cat(paste0("\n\n****************** Chunk ", chunks$from[chunkix], "-", chunks$to[chunkix], ": calculating distances ******************\n"))

    ### form a name of the output, if the file exists, remove it
    ou07fnm <- paste0(substr(as.character(10^n.char+as.numeric(chunks$from[chunkix])), 2, 1+n.char), "_", substr(as.character(10^n.char+as.numeric(chunks$to[chunkix])), 2, 1+n.char))
    ou07fnm <- paste0(tmpdirnm, "/", ou07fnm, ".Rdata") ; if (file.exists(ou07fnm)) { tmp <- do.call(file.remove, list(ou07fnm)) ; rm(tmp) }

    ### determine the indices of protein pairs in the current chunk
    ixs <- chunks$from[chunkix]:chunks$to[chunkix]

    ### build the current chunk of protein pairs
    data_var <- protein_pairs[ixs, ]
    tmp <- data_XY[, semantics_XY$sample] ; colnames(tmp) <- paste0("x_", colnames(tmp)) ; tmp <- cbind(data.frame(ACx=data_XY$AC), tmp)
    data_var <- plyr::join(data_var, tmp, by='ACx', type="left")
    tmp <- data_XY[, semantics_XY$sample] ; colnames(tmp) <- paste0("y_", colnames(tmp)) ; tmp <- cbind(data.frame(ACy=data_XY$AC), tmp)
    data_var <- plyr::join(data_var, tmp, by='ACy', type="left")

    ### build a "separate" model
    model.on.plex <- paste0("plex*(", model.on, ")") ; stats.on.plex <- paste0("plex_y:", gsub("-", "-plex_y:", stats.on))
    t <- step06_ANOVA(data_var, semantics_var, model.on=model.on.plex, stats.on=stats.on.plex, verbose=TRUE, splineDF=splineDF)
    separate <- t$fit.properties

    ### get the interaction coefficients
    t.intCoeff <- t[[names(t)[grepl("^plex_y:", names(t))][1]]] ### first of the coefficient/pval tables that contains the interaction coeff; hopefully there is only one

    ### build a "collated" model
    t <- step06_ANOVA(data_var, semantics_var, model.on=paste0(model.on, "+plex"), stats.on=stats.on, verbose=FALSE, splineDF=splineDF)
    collated <- t$fit.properties
    nm <- "" ; for (ix in length(names(t)):1) { if (sum(grepl("p.pVal", colnames(t[[ix]]), fixed=TRUE))>0) { nm <- names(t)[ix] } }
    t.and <- t[[nm]] ### first of the coefficient/pval tables, hope this is the right one

    rm(t, data_var)

    save(list=c("collated", "separate", "t.intCoeff", "t.and", "model.on.plex", "stats.on.plex"), file=ou07fnm)
  } ### end function to perform calculations on chunk #chunkix

  ###########################################################
  ### main activities of step07_distanceMatrix start here ###
  ###########################################################

  stats.on <- stats.on[1] ### take the first one if there are many

  if (sum(duplicated(data_XY$AC)) > 0) {
    warning("The AC column of data_XY contains duplicates, make it unique before proceeding")
    return(NULL)
  }

  ### first, calculate sigCol on proteins in accordance with model.on and stats.on
  ### at this step, no interactions or AICs or whatever are considered, even if data_XY allows for it
  ### we are just assessing the individual profiles here before we start comparing profile pairs.
  cat("\nInfo (step07_distanceMatrix): initializing protein_pairs...\n")
  t <- step06_ANOVA(data_XY, semantics_XY, model.on=model.on, stats.on=stats.on, verbose=FALSE, splineDF=splineDF)
  nm <- "" ; cnt <- 0
  for (ix in length(names(t)):1) { if (sum(grepl("p.pVal", colnames(t[[ix]]), fixed=TRUE))>0) { nm <- names(t)[ix] ; cnt<-cnt+1 } }
  if (cnt>1) { ### if stats.on is ambiguous (i.e. generates more than one p-val), give a warning and take the first p-val
    warning(paste0("stats.on=(", paste0(stats.on, collapse=", "), ") generated more than one p-value-containing data frame in the output, using the first one, ", nm))
    }
  t <- t[[nm]] ; sigCol <- "p.pVal.adj" ; data_XY[, sigCol] <- t[, sigCol] ; rm(cnt, nm, t)

  #data_XY <- data_XY[order(data_XY[, sigCol], decreasing=TRUE)[1:100], ]
  #data_XY <- data_XY[!is.na(data_XY[, sigCol]), ]

  ### generate a data frame of pairs using AC as the key
  protein_pairs <- base::expand.grid(data_XY$AC, data_XY$AC) ; colnames(protein_pairs) <- c("AC", "Var2")
  protein_pairs <- plyr::join(protein_pairs, data_XY[, c("AC", sigCol)], by='AC', type="left") ; colnames(protein_pairs) <- c("ACx", "AC", "p.pVal.x")
  protein_pairs <- plyr::join(protein_pairs, data_XY[, c("AC", sigCol)], by='AC', type="left") ; colnames(protein_pairs) <- c("ACx", "ACy", "p.pVal.x", "p.pVal.y")
  protein_pairs$ACx <- as.character(protein_pairs$ACx) ; protein_pairs$ACy <- as.character(protein_pairs$ACy)

  protein_pairs$x.or.y <- rje::rowMaxs(protein_pairs[, c("p.pVal.x", "p.pVal.y")])
  protein_pairs$x.and.y <- NA ### collated significance on stats.on
  protein_pairs$x.vs.y.11  <- NA ### difference x.or.y-x.and.y (will be filled at the bottom)
  protein_pairs$x.vs.y.int <- NA ### interaction coefficients
  protein_pairs$x.vs.y.F   <- NA ### Prism-style F-test
  protein_pairs$x.vs.y.F.stat   <- NA ### Prism-style F-test statistic
  protein_pairs$x.vs.y.AIC <- NA ### Prism-style AIC

  ### generate semantics for protein pairs
  semantics_x <- semantics_XY ; semantics_x$sample <- paste0("x_", semantics_x$sample) ; semantics_x$plex <- "x"
  semantics_y <- semantics_XY ; semantics_y$sample <- paste0("y_", semantics_y$sample) ; semantics_y$plex <- "y"
  semantics_var <- rbind(semantics_x, semantics_y) ; rm(semantics_x, semantics_y)

  ### create a temporary directory
  letters.digits <- c(letters, unlist(stringr::str_split("123456789", "")))
  s_tempDir <- dirname(normalizePath(tempdir()))
  while (TRUE) {
    tmpdirnm <- paste0(c("step07_DM_tmp", letters[unlist(ramify::randi(length(letters), 1)[1, 1])], letters.digits[unlist(ramify::randi(length(letters.digits), 5)[, 1])]), collapse="")
    tmpdirnm <- paste0(s_tempDir, "/", tmpdirnm)
    if (!dir.exists(tmpdirnm)) { dir.create(tmpdirnm, recursive=TRUE) ; break ; }
  }
  n.char <- ceiling(log10(nrow(protein_pairs)+1)) ### an accessory variable for forming chunk-dependent file names with "0001" or similar
  cat("\nInfo (step07_distanceMatrix): temporary directory set to", tmpdirnm, "\n")

  ### will go over protein pairs in chunks, because attempts to do the same in one go fail due to memory limits
  chunk <- 5000 ; nchunks <- nrow(protein_pairs) %/% chunk
  rem <- nrow(protein_pairs) %% chunk ; chunk <- chunk + (rem %/% nchunks) + 1
  nchunks <- ceiling(nrow(protein_pairs) / chunk) ; pb <- txtProgressBar(min=1, max=nchunks, style=3)

  ### build a two-column dataframe chunks(from, to)
  chunks <- data.frame(from=rep(NA, nchunks), to=rep(NA, nchunks)) ; from <- 1 ; i <- 1 ;
  while (from < nrow(protein_pairs)) { chunks$from[i] <- from ; chunks$to[i] <- min(nrow(protein_pairs), from+chunk-1) ; from <- from+chunk ; i <- i+1 ; }
  ti <- Sys.time() ;

  if (!Sys.info()["user"] == "anl70"){nproc <- min(floor(parallel::detectCores()/2), nproc)}

  cat("\nInfo (step07_distanceMatrix): calculating stats on protein_pairs in", nchunks, "chunk(s) of", chunk, "pair(s) each; using", nproc, "thread(s)\n")
  if (nproc < 2) { ### 2023-03-31 the usual boring single-threaded loop
    for (chunkix in 1:nrow(chunks)) step07_distanceMatrix_onChunk(chunkix=chunkix)
  } else { ### try to parallelize chunk processing (maybe this is nonsensical, memory limits are still a thing)
    virtualCluster.oufnm <- paste0(tmpdirnm, "/virtualCluster.ou") ; if (file.exists(virtualCluster.oufnm)) { tmp <- do.call(file.remove, list(virtualCluster.oufnm)) ; rm(tmp) }
    virtualCluster <- parallel::makeCluster(nproc, outfile=virtualCluster.oufnm)
    parallel::clusterExport(cl=virtualCluster,
                            varlist=c("dirOut", "dirCode", "step07_distanceMatrix_onChunk", "step06_ANOVA",
                                      "protein_pairs", "chunks", "n.char", "tmpdirnm",
                                      "semantics_XY", "model.on", "stats.on"),
                            envir=environment())
    tmp <- parallel::clusterApplyLB(cl=virtualCluster, x=1:nrow(chunks), fun='step07_distanceMatrix_onChunk')
    #tmp <- parallel::parLapply(cl=virtualCluster, X=1:nrow(chunks), fun='step07_distanceMatrix_onChunk')
    parallel::stopCluster(virtualCluster)

    if (Sys.info()["sysname"] == "Windows") {
      nrow(installr::get_tasklist()) ### how many processes (all kinds of) currently?
      installr::kill_all_Rscript_s()
      nrow(installr::get_tasklist()) ### kill all Rscript processes and count total again
    }

    #if (file.exists(virtualCluster.oufnm)) { tmp <- do.call(file.remove, list(virtualCluster.oufnm)) ; rm(tmp) } ### remove the  virtual cluster output file
  }
  cat("\nInfo (step07_distanceMatrix): stats calculation completed in", format(difftime(Sys.time(), ti, units="mins"), digits=4), "\n")

  ### now that chunk data is calculated, fill in the columns of protein_pairs
  cat("\nInfo (step07_distanceMatrix): tranferring the calculated stats into protein_pairs in a single thread...\n"); pb <- txtProgressBar(min=1, max=nchunks, style=3)
  for (chunkix in 1:nrow(chunks)) {

    ### determine the indices of protein pairs in the current chunk
    ixs <- chunks$from[chunkix]:chunks$to[chunkix]

    ### clear the environment from the previous iteration data
    nms <- c("collated", "separate", "t.intCoeff", "t.and", "model.on.plex", "stats.on.plex") ; nms <- nms[nms %in% ls()] ; if (length(nms)>0) rm(list=nms) ; rm(nms)

    ### form a name of the output, if the does not exists, skip
    ou07fnm <- paste0(substr(as.character(10^n.char+as.numeric(chunks$from[chunkix])), 2, 1+n.char), "_", substr(as.character(10^n.char+as.numeric(chunks$to[chunkix])), 2, 1+n.char))
    ou07fnm <- paste0(tmpdirnm, "/", ou07fnm, ".Rdata") ; if (!file.exists(ou07fnm)) { cat("\nError (step07_distanceMatrix): file", ou07fnm, "does not exist; skipping\n\n") ; next ; }

    ### load the next chunk, if incomplete, skip
    load(file=ou07fnm) ; if (!all(c("collated", "separate", "t.intCoeff", "t.and") %in% ls())) { cat("\nError (step07_distanceMatrix): file", ou07fnm, "lacks one or more of collated, separate, t.and, and t.intCoeff; skipping\n\n") ; next ; }

    ### fill in collated significance
    protein_pairs$x.and.y[ixs] <- t.and$p.pVal.adj

    ### fill in interaction coefficients
    protein_pairs$x.vs.y.int[ixs] <- t.intCoeff$p.pVal.adj

    ### fill in Prism-style F-test
    F <- ((collated$ss-separate$ss)/(collated$df.residual-separate$df.residual)) / (separate$ss/separate$df.residual)
    dfn <- collated$df.residual-separate$df.residual ; dfd <- separate$df.residual ;
    protein_pairs$x.vs.y.F[ixs] <- -log10(pf(F, dfn, dfd, lower.tail=FALSE))
    protein_pairs$x.vs.y.F.stat[ixs] <- F ## Adding stat that does not account for df
    rm(F, dfn, dfd)

    ### fill in Prism-style AIC
    protein_pairs$x.vs.y.AIC[ixs] <- -(separate$AICc-collated$AICc)/2*log10(exp(1))

    setTxtProgressBar(pb, chunkix)
  } ### end filling in the columns of protein_pairs
  close(pb)

  protein_pairs$x.vs.y.11 <- protein_pairs$x.or.y - protein_pairs$x.and.y

  errcode <- unlink(tmpdirnm, recursive=TRUE) ### this supposed to remove the temporary directory but does not by whatever reason
  if (errcode == 0) { cat("\nInfo (step07_distanceMatrix): temporary directory", tmpdirnm, "deleted\n")
  } else { cat("\nError (step07_distanceMatrix): cannot delete the temporary directory", tmpdirnm, ", error code", errcode, "\n") }
  rm(ixs, semantics_var, tmpdirnm, pb)

  return(protein_pairs)
}
