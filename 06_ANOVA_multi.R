#library(limma) - now specified explicitly with :: at every limma function call
#library(splines) - now specified explicitly with :: at every splines function call
#set.seed(7)

### 2023-01-24 IK thinking aloud ###
### maybe we should revise this and finally make it into a function
### there are two things that ANOVA analysis needs to know:
###  * what to run the analysis on (and it can be an arbitrarily complicated formula).
###  * for what coefficients to return stat significance
###  * Examples:
###    - time (run on)
###      - time=3min (return pval for)
###      - time=3min and time=10min (return pval for)
###      - etc.
###    - cell_line (run on)
###      - cell_line (return pval for, there are no other options here)
###    - time+cell_line (run on)
###      - cell_line (return pval on - this is supposed to prioritize offsets but responses may be asynchronous or synchronous, and HQ or crappy)
###      - time (return pval on - synchronous responses with or without a cell-line-mediated offset)
###      - all coeff (return pval on - not sure what this will prioritize)
###    - time*cell_line (run on; there will be interaction coeff)
###      - time (return pval on - prioritize HQ responses in the 1st condition only)
###      - cell_line (return pval on - this is supposed to prioritize offsets but responses may be asynchronous or synchronous, and HQ or crappy)
###      - time*cell_line (return pval on - prioritizes distinct responses regardless of offset)
###      - etc.
###    - ... and similar for CCR5...
###  * for the coefficients themselves, there is nothing wrong with returning them all.
### And what is the function supposed to return? res? ANOVA_result?
###
### Something along the lines of:
### 06_ANOVA(data_var, semantics_var, model.on="condition", stats.on="condition", verbose=TRUE, splineDF=0)
### where model.on is a single string (ignore anything past the 1st element if it is a list)
### and stats.on can be a list of either semantics column names or design matrix column names,
### or greppable substrings for the latter
###
### R formula documentation: https://www.rdocumentation.org/packages/stats/versions/3.6.2/topics/formula

################################################################################################################
### the function returns a list with the following elements:
### * fit.properties: a df with columns
###   - df.residual from limma::lmFit()
###   - sigma       from limma::lmFit()
###   - ss - sum of squares of the residuals
###   - AIC  - Akaike information criterion
###   - AICc - Akaike information criterion with a correction for low DFs
### * spline.basis: a matrix where
###   - the first column is X
###   - the remaining columns are Y's of the basis functions
###   - NOTE! - only the last basis is saved, so the function is not going
###     to work properly with two or more orthogonal numerical columns
### * one data frame for each  element of stats.on, containing
###   - intercept
###   - coefficients for the relevant design matrix columns
###   - avg.quant
###   - F from limma::topTable()
###   - p.pVal from limma::topTable()
###   - p.pVal.adj from limma::topTable()
################################################################################################################

step06_ANOVA <- function(data_var, semantics_var, model.on="condition", stats.on="condition", verbose=TRUE, splineDF=0, l_sqrtTime=FALSE) {

  if (missing(verbose) || is.null(verbose)) { verbose <- TRUE }

  ################################################################
  ### The model.on variable specifies on which column(s) #########
  ### of the semantics data frame the model will be built. #######
  ### It can have a syntax of a formula. #########################
  ################################################################
  if (missing(model.on) || is.null(model.on) || model.on=="") { model.on <- "condition" } ; model.on <- model.on[1] ### if a list, use only the 1st element
  if (!grepl("^~", model.on, fixed=FALSE)) { model.on <- paste0("~", model.on) }
  semacol <- all.vars(as.formula(model.on))

  if (verbose) {cat(paste0("\n", length(semacol), "-way ANOVA using semantics_var column(s): ", paste0(semacol, collapse=", "), "\n"))}

  ### initial/default number degrees of freedom for spline modeling of numeric columns
  ### if 0, will use experimentally measured values of numeric semantics_var columns as knots.
  ### if -1, will place knots between exp points, except for the boundaries; this option may be removed in the future because it is not so useful
  ### if NULL, will not use splines at all but will instead treat numerical variables as categorical
  if (missing(splineDF)) { splineDF <- 0 }

  ### determine which columns in semacol are numerical
  semacolnum <- (semacol == "") ### initialize with all FALSE's
  for (colix in 1:length(semacol)) { semacolnum[colix] <- (sum(is.na(as.numeric(semantics_var[, semacol[colix]]))) == 0) }

  if (verbose & (sum(semacolnum)==0)) { cat("\nNo numerical columns in model.on, all variables are categorical\n") }
  if (verbose & (sum(semacolnum)>0) &  is.null(splineDF)) { cat("\nNumerical columns will be treated as categorical, no splines will be used\n") }
  if (verbose & (sum(semacolnum)>0) & !is.null(splineDF)) {
    cat(paste0("\nColumn(s) ", paste0(semacol[semacolnum], collapse=", "), " is/are numeric and will be modeled using splines"))
    cat(paste0("\nInitial # of DFs for the splines is set to ", splineDF, "\n"))
  }

  ### for multi-way ANOVA, build a list of factors or spline bases, one for each column in semacol
  ways <- semacol ### initialize with just semacol values
  if (verbose) { cat("\nCalculating factors/spline bases for limma:") }
  for (colix in 1:length(semacol)) {
    colVals <- unique(semantics_var[, semacol[colix]])
    ways[colix] <- paste0(substr(semacol[colix], 1, 4), "_")
    if (is.null(splineDF) | !semacolnum[colix]) { ### either we are not using splines at all, or this is a non-numeric column: use factors
      if (verbose) {cat(paste0("\n * column ", semacol[colix], ", factor levels: ", paste0(colVals, collapse=", ")))}
      assign(ways[colix], factor(semantics_var[, semacol[colix]], levels=colVals))
    } else { ### this is a numeric column, use splines
      if (l_sqrtTime) {colVals <- sqrt(as.numeric(colVals)); colVals[colVals == "-Inf"] <- 0} else {colVals <- as.numeric(colVals)}
      minVal <- min(colVals) ; maxVal <- max(as.numeric(colVals))
      ### set df_colix. df_colix cannot be more than the # of unique values in semantics_var[, semacol[colix]] minus 1.
      df_colix <- splineDF ; if (df_colix > length(colVals)-1) { df_colix <- length(colVals)-1 }
      ### calculate the natural spline basis for the entire range of semacol[colix] values, not just the experimentally measured ones
      ### make sure the experimentally measured points are included
      splineBasisX <- unique(c(seq(minVal, maxVal, length.out=61), colVals)) ; splineBasisX <- splineBasisX[order(splineBasisX)]
      ### build the spline knots
      if (df_colix < 1) { ### knots guided by the data
        splineKnotsX <- colVals[order(colVals)]
        if (df_colix < 0) { splineKnotsX <- (splineKnotsX[2:(length(splineKnotsX))] + splineKnotsX[1:(length(splineKnotsX)-1)])/2. }
      } else { ### equally spaced knots
        splineKnotsX <- (0:df_colix) * (maxVal-minVal) / df_colix
      }
      splineBasis <- splines::ns(splineBasisX, knots=splineKnotsX[2:(length(splineKnotsX)-1)], Boundary.knots=splineKnotsX[c(1, length(splineKnotsX))]) ; df_colix <- length(splineKnotsX)-1
      if (verbose) {cat(paste0("\n * column ", semacol[colix], ", splines with ", df_colix, " degrees of freedom and knots at ", paste0(splineKnotsX, collapse="/")))}
      if (l_sqrtTime){
        tmp <- plyr::join(data.frame(x=sqrt(as.numeric(semantics_var[, semacol[colix]])) ), cbind(data.frame(x=splineBasisX), splineBasis), by="x", type="left")
      } else{
        tmp <- plyr::join(data.frame(x=as.numeric(semantics_var[, semacol[colix]])), cbind(data.frame(x=splineBasisX), splineBasis), by="x", type="left")
      }
      assign(ways[colix], as.matrix(tmp[, colnames(splineBasis)])) ; rm(tmp)
    }
  } ### at the completion of this loop, have "ways". Also can have splineBasisX and splineBasis for *each* numerical column in semacol. But need to save them somehow...
  if (verbose) { cat("\n") } ; if ((sum(semacolnum)>0) & !is.null(splineDF)) { rm(df_colix, maxVal, minVal) } ; rm(colix, colVals, splineDF, semacolnum)

  ### Now use "ways" to create design matrix for limma: this is supposed to be a universal way for all projects and model.on values
  ### Side benefit: the names of coefficients are more descriptive now
  ### To generate the formula for the design matrix, simply replace elements of semacol with elements of ways in model.on
  fml <- model.on ; for (colix in 1:length(semacol)) { fml <- gsub(semacol[colix], ways[colix], fml, fixed=TRUE) } ; rm(colix)
  #Used to be:
  #fml <- paste0("~", paste0(ways, collapse="*")) #This generates a string to be passed to as.formula()
  #fml <- paste0("~", paste0(ways, collapse="+")) #This generates a string to be passed to as.formula() but with a +
  design <- model.matrix(as.formula(fml)) #Creates the design matrix input for limma. Based off "ways"

  ### Build the linear model
  if (verbose) { cat(paste0("\nBuilding a linear model for ", fml, "\n")) }
  res0 <- limma::lmFit(data_var[ , semantics_var$sample], design) ; res <- limma::eBayes(res0)
  ANOVA_result <- data.frame(df.residual=res$df.residual, sigma=res$sigma)

  ### R-squared values
  #sst <- rowSums(data_var[ , semantics_var$sample] ^ 2, na.rm = TRUE)
  #ssr <- sst - res$df.residual*res$sigma^2
  #ANOVA_result$R_squared <- ssr/sst; rm(ssr, sst)

  ################################################################################################
  #### Calculate data for subsequent application of Prism-like approach to curve comparisons #####
  #### For this, ANOVA_result needs the following: ###############################################
  #### * df.residual (which it already has) ######################################################
  #### * ss, only for the given model ############################################################
  #### * AIC and AICc, although these are kind of a cherry on top, they are easy enough ##########
  ####   to calculate from df and ss but still ###################################################
  ################################################################################################

  ### calculate residuals and sums-of-squares; append to ANOVA_result
  resi <- limma::residuals.MArrayLM(res, data_var[ , semantics_var$sample])
  ANOVA_result$ss <- rowSums(resi*resi, na.rm=TRUE) ; ANOVA_result$ss[rowSums(!is.na(resi))<1] <- NA

  ### calculate AIC (Akaike Information Criterion, see Prism manual)
  npoints <- rowSums(!is.na(data_var[ , semantics_var$sample]))
  ANOVA_result$AIC <- npoints*log(ANOVA_result$ss/npoints)+2*(ncol(res$design)+1)

  ### when # of points is small compared to ncol(res$design)+1, AIC correction is needed:
  nparam <- ncol(res$design)+1
  ixs <- npoints > nparam+1 ; correction <- rep(NA, nrow(ANOVA_result))
  correction[ixs] <- 2*nparam*(nparam+1)/(npoints[ixs]-(nparam+1))
  ANOVA_result$AICc <- ANOVA_result$AIC + correction ;

  rm(resi, nparam, npoints, ixs, correction)

  ### in the calling script, use the following to calculate probability that the "separate" model is correct
  ### dAICc <- AICcSeparate-AICcCollated
  ### prob <- exp(-dAICc/2.)/(1.+exp(-dAICc/2))

  ### Initialize result - a list of ANOVA_results style data frames
  result <- list() ; result[['fit.properties']] <- ANOVA_result ;
  if (exists("splineBasisX") && exists("splineBasis")) { result[['spline.basis']] <- cbind(x=splineBasisX, splineBasis) }

  ##########################################################################
  ### Calculate stats for the sets of coefficients described in stats.on ###
  ### stats.on can be a list of: ###########################################
  ### * semantics_var column names or their combinations via ":" ###########
  ### * design matrix column names or their combinations via ":" ###########
  ### * any greppable substrings for design matrix column names ############
  ##########################################################################
  ### Examples of valid elements in stats.on: ##############################
  ### "ligand" #############################################################
  ### "liga_" (from ways) ##################################################
  ### "liga_6P4" (from design matrix) ######################################
  ### "ligand:time" ########################################################
  ### "liga_6P4:time_" #####################################################
  ### "*" (this means all columns) #########################################
  ### Special values for stats.on: #########################################
  ### "_all_" - same as "*", call columns of the design matrix #############
  ### "_none_" - same as NULL, "", or not provided at all, no stats ########
  ### "_terms_" - terms of the model.on formula ############################
  ### whenever there is a "-" in stats.on, trying to treat it as a contrast#
  ### etc. #################################################################
  ### An ANOVA_result table will be built for each stats.on element ########
  ##########################################################################

  if (missing(stats.on) || is.null(stats.on) || paste0(stats.on, collapse="")=="" || length(stats.on)<1 || sum(stats.on=="_none_")>0) { ### stats.on is not available
    if (verbose) { cat("\nstats.on is unavailable, no stats will be calculated\n\n") }
    rm(list=ways) ; rm(verbose, res, design, ANOVA_result, fml, ways, semacol)
    return(result)
  }

  if (sum(stats.on ==  "_all_")>0) { stats.on <- "*" }
  if (sum(stats.on == "_term_")>0) { stats.on <- labels(terms(as.formula(model.on))) }

  ### if stats.on entries have a "-" in them, that indicates that we need to use contrast, w.r.to will specify contrasts with respect to what
  w.r.to   <- unlist(lapply(strsplit(paste0(stats.on, "- "), "-"), "[[", 2)) ; w.r.to <- gsub(" ", "", w.r.to)
  stats.on <- unlist(lapply(strsplit(paste0(stats.on, "- "), "-"), "[[", 1))

  ### build greppable versions of stats.on entries and w.r.to entries
  stats.on.greppable <- stats.on ; w.r.to.greppable <- w.r.to
  for (colix in 1:length(semacol)) {
    stats.on.greppable <- gsub(semacol[colix], ways[colix], stats.on.greppable, fixed=TRUE)
    w.r.to.greppable   <- gsub(semacol[colix], ways[colix], w.r.to.greppable,   fixed=TRUE)
  } ; rm(colix)
  stats.on.greppable <- paste0("^", stats.on.greppable) ; stats.on.greppable <- gsub("^^", "^", stats.on.greppable, fixed=TRUE)
  stats.on.greppable <- paste0(stats.on.greppable, "$") ; stats.on.greppable <- gsub("$$", "$", stats.on.greppable, fixed=TRUE)
  stats.on.greppable <- gsub("__", "_", stats.on.greppable)
  stats.on.greppable <- gsub("_:", "_[^:]{1,}:", stats.on.greppable, fixed=TRUE)
  stats.on.greppable <- gsub("_$", "_[^:]{1,}$", stats.on.greppable, fixed=TRUE)
  w.r.to.greppable <- paste0("^", w.r.to.greppable) ; w.r.to.greppable <- gsub("^^", "^", w.r.to.greppable, fixed=TRUE)
  w.r.to.greppable <- paste0(w.r.to.greppable, "$") ; w.r.to.greppable <- gsub("$$", "$", w.r.to.greppable, fixed=TRUE)
  w.r.to.greppable <- gsub("__", "_", w.r.to.greppable)
  w.r.to.greppable <- gsub("_:", "_[^:]{1,}:", w.r.to.greppable, fixed=TRUE)
  w.r.to.greppable <- gsub("_$", "_[^:]{1,}$", w.r.to.greppable, fixed=TRUE)

  ### The section above used to be:
  ### If we want to  select the coef columns manually, ie only the interaction coef or only change with treatment w/o respect to time
  #if (!exists("coefCols")) {coefCols <- NULL} #Manually select which coefficient to run statistical analysis on
  #if (is.null(coefCols)) {coefCols <- colnames(design)[grep(":", colnames(design), fixed=TRUE)]}
  #if (length(semacol) == 1) {coefCols <- colnames(design)[colnames(design) != "(Intercept)"]} #If only one-way, use all coeficients, no intercept
  #coefCols <- colnames(design)[!grepl("Intercept", colnames(design), fixed=TRUE) & !grepl(":", colnames(design), fixed=TRUE)]

  ### now finally calculate stats, with or without contrasts
  if (verbose) { cat("\nCalculating stats:") } ; t <- data.frame()
  for (colix in 1:length(stats.on)) { ### calculate the significances

    coeff.col <- colnames(design)[grepl(stats.on.greppable[colix], colnames(design)) & colnames(design)!="(Intercept)"]
    if (length(coeff.col) < 1) {
      if (verbose) { cat(paste0("\n * \"", stats.on[colix], "\" (grep \"", stats.on.greppable[colix], "\") does not match any design columns (available: ", paste0(colnames(design), collapse=", "),"), skipping")) }
      next()
    }

    w.r.to.col <- colnames(design)[grepl(w.r.to.greppable[colix], colnames(design)) & colnames(design)!="(Intercept)"]
    if (length(w.r.to.col) > 1) { ### this is an ambiguous situation
      if (verbose) { cat(paste0("\n * contrast \"", stats.on[colix], "-", w.r.to[colix], "\" (grep \"", stats.on.greppable[colix], "\" w.r.t. grep \"", w.r.to.greppable[colix],"\") is ambiguous, skipping")) }
      next()
    }

    if (verbose) {
      if (length(w.r.to.col) < 1) {
        cat(paste0("\n * on \"", stats.on[colix], "\" (grep \"", stats.on.greppable[colix], "\"): design matrix column(s) ", paste0(coeff.col, collapse=", "), ""))
      } else {
        cat(paste0("\n * on \"", stats.on[colix], "-", w.r.to[colix], "\" contrast (grep \"", stats.on.greppable[colix], "\" w.r.t. grep \"", w.r.to.greppable[colix],"\"): design matrix column(s) ", paste0(coeff.col, collapse=", "), " vs ", paste0(w.r.to.col, collapse=", ")))
      }
    }

    if (length(w.r.to.col) == 0) { resC <- res } else { ### go back to the pre-eBayes res0 model and calculate the fit on contrasts, not just on coefficients
      resC <- res0 ### a version with contrasts
      contrasts <- paste0(coeff.col, "-", w.r.to.col)
      design.colnames <- colnames(resC$design)       ; colnames(resC$design)       <- make.names(colnames(resC$design))       ### save the existing design colnames, then make them "syntactically correct"
      coeff.colnames  <- colnames(resC$coefficients) ; colnames(resC$coefficients) <- make.names(colnames(resC$coefficients)) ### save the existing coefficients colnames, then make them "syntactically correct"
      contrast.matrix <- limma::makeContrasts(contrasts=contrasts, levels=resC$design)
      resC <- limma::contrasts.fit(resC, contrast.matrix) ; resC <- limma::eBayes(resC)
      coeff.col <- colnames(resC$coefficients)
    }

    t <- limma::topTable(resC, adjust="BH", coef=coeff.col, number=nrow(data_var), sort.by="none")
    ### other p.adjust.methods: c("holm", "hochberg", "hommel", "bonferroni", "BH", "BY", "fdr", "none")
    colnames(t) <- gsub("AveExpr", "avg.quant", colnames(t))
    if (sum(grepl("logFC", colnames(t))) > 0) { ### this is a two-condition, t-test scenario
      colnames(t) <- gsub("logFC", coeff.col[1], colnames(t))
      colnames(t) <- gsub("^t$", "F", colnames(t)) ; t$F <- t$F*t$F ### F=t^2 by definition
      t <- within(t, rm(B))
    }
    t$p.pVal <- -log(t$P.Value, base=10.) ; t$p.pVal.adj <- -log(t$adj.P.Val, base=10.) ; t <- within(t, rm(P.Value, adj.P.Val))

    ### add intercept or 0 if it is not available
    result[[stats.on[colix]]] <- cbind(data.frame(intercept=0), t)
    if (sum(grepl("(Intercept)", colnames(resC$coefficients)))>0) { result[[stats.on[colix]]]$intercept <- resC$coefficients[, "(Intercept)"] }
  }
  if (verbose) { cat("\n\n") } ; rm(colix, coeff.col)

  ### something Alexis implemented, Alexis please clarify what this is
  ## This was implemented for the original, non function version. It's meant to create the ANOVA_results df. Does not seem necessary anymore with new function architecture
  if (FALSE) {
    if (TRUE) {data_var <- cbind(data.frame(df.residual=ANOVA_result$df.residual, R_squared = ANOVA_result$R_squared, F=ANOVA_result$F, p.pVal=ANOVA_result$p.pVal, p.pVal.adj=ANOVA_result$p.pVal.adj), data_var)}
    if (!("ACx" %in% colnames(data_var))) {ANOVA_result <- cbind(AC = data_var$AC, protein = data_var$protein, ANOVA_result)
    } else {ANOVA_result <- cbind(data.frame(ACx = data_var$ACx, ACy = data_var$ACy), ANOVA_result)}
  }

  rm(list=ways) ; rm(verbose, res, design, ANOVA_result, t, fml, ways, semacol, stats.on.greppable)

  return(result)
}
