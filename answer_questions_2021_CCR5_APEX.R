### the call below defines a function, step06_ANOVA(data_var, semantics_var, model.on="condition", stats.on="condition", verbose=TRUE, splineDF=0)
source(paste0(dirCode, "06_ANOVA_multi.R"))

### and the call below defines a function, step07_distanceMatrix(data_XY, semantics_XY, model.on="condition", stats.on="condition", nproc=1), in its PARALLEL implementation
source(paste0(dirCode, "07_distanceMatrix_parallel.R"))

if (!exists("l_clustering")) {l_clustering <- FALSE} ### just a simple way to avoid any clustering, it takes a while to do
if (!exists("l_plots")){l_plots <- TRUE} ### if we only need statistics, use this to skip plotting
if (!exists("var_from")){var_from_ <- 1}
nplot <- 20; ggplot2::theme_update(text=ggplot2::element_text(size=12))
ncol <- 5
splineDF <- 0 ### will use this value whenever possible when modeling numerical data
l_sqrtTime <- TRUE

##################################################################################################################

### Q_CCL5: How does the proximity environment of CCR5 vary over time when cells are stimulated with RANTES(=CCL5)?
### Q_5P14: How does it vary over time when cells are stimulated with [5P14]CCL5?
### Q_6P4 : How does it vary over time when cells are stimulated with [6P4]CCL5?
### Q3 "In the basal state..." is irrelevant for the CCR5 experiment
### Q_CCL5_vs_5P14: How do temporally-resolved responses to CCL5 differ from those to 5P14?
### Q_CCL5_vs_6P4:  How do temporally-resolved responses to CCL5 differ from those to 6P4?
### Q_5P14_vs_6P4:  How do temporally-resolved responses to 5P14 differ from those to 6P4? ... Will this be a concern ?

if (!exists("Qs")){Qs <- c("Q_CCL5", "Q_5P14", "Q_6P4", "Q_CCL5_vs_5P14", "Q_CCL5_vs_6P4")} #If Qs is assigned externally, it will not be overwritten

for (Q in Qs) {
  cat(paste0("\n\n************************* Question ", Q, " *************************\n"))

  semantics_var <- semantics_wide[semantics_wide$cell_line=="CCR5", ]
  ### triplicate the ligand="none" rows of the semantics_var data frame, to represent t=0 for each separate ligand.
  ### the structure of the data_var data frame is unchanged, these "none" columns/sample names are each just used three times in semantics_var
  ### this is necessary for the multi-way analysis in ANOVA_multi
  tmp <- semantics_var[semantics_var$ligand == "none", ] ; semantics_var <- semantics_var[!(semantics_var$sample %in% tmp$sample), ]
  tmp$ligand <- "RANTES" ; semantics_var <- rbind(tmp, semantics_var)
  tmp$ligand <- "5P14"   ; semantics_var <- rbind(tmp, semantics_var)
  tmp$ligand <- "6P4"    ; semantics_var <- rbind(tmp, semantics_var)
  semantics_var <- semantics_var[order(semantics_var$time), ]
  rm(tmp)

  ### seantics_var can be reordered so that the "time" coef corresponds to specific ligand
  if (Q=="Q_CCL5") semantics_var <- rbind(semantics_var[grepl("RANTES", semantics_var$ligand), ],
                                          semantics_var[grepl("5P14"  , semantics_var$ligand), ],
                                          semantics_var[grepl("6P4"   , semantics_var$ligand), ])
  if (Q=="Q_5P14") semantics_var <- rbind(semantics_var[grepl("5P14"  , semantics_var$ligand), ],
                                          semantics_var[grepl("RANTES", semantics_var$ligand), ],
                                          semantics_var[grepl("6P4"   , semantics_var$ligand), ])
  if (Q=="Q_6P4")  semantics_var <- rbind(semantics_var[grepl("6P4"   , semantics_var$ligand), ],
                                          semantics_var[grepl("RANTES", semantics_var$ligand), ],
                                          semantics_var[grepl("5P14"  , semantics_var$ligand), ])

  if (Q=="Q_CCL5_vs_5P14") semantics_var <- rbind(semantics_var[grepl("RANTES", semantics_var$ligand), ],
                                                  semantics_var[grepl("5P14"  , semantics_var$ligand), ])
  if (Q=="Q_CCL5_vs_6P4")  semantics_var <- rbind(semantics_var[grepl("RANTES", semantics_var$ligand), ],
                                                  semantics_var[grepl("6P4"   , semantics_var$ligand), ])

  ### Filter 10 min samples if analyzing questions with 5P14
  if (TRUE && grepl("5P14", Q)){
    semantics_var <- semantics_var[ semantics_var$time != "10", ]
    cat("\n####################\nREMOVING 10 MIN 5P14\n####################\n")
  }


  ### calculate a "separate" fit and save model coefficients for plotting
  separate <- step06_ANOVA(data_wide_norm, semantics_var, model.on="ligand*time", stats.on="_term_", verbose=TRUE, splineDF=splineDF)
  coeff_separate <- list()
  t_ligand   <- separate[["ligand"]]      ; t_ligand   <- t_ligand  [ , !(colnames(t_ligand)   %in% c("avg.quant", "F", "p.pVal", "p.pVal.adj")), drop=FALSE]
  t_time     <- separate[["time"]]        ; t_time     <- t_time    [ , !(colnames(t_time)     %in% c("intercept", "avg.quant", "F", "p.pVal", "p.pVal.adj")), drop=FALSE]
  t_interact <- separate[["ligand:time"]] ; t_interact <- t_interact[ , !(colnames(t_interact) %in% c("intercept", "avg.quant", "F", "p.pVal", "p.pVal.adj")), drop=FALSE]
  if (Q=="Q_CCL5") coeff_separate[["RANTES"]] <- as.matrix(cbind(data.frame(ligand=t_ligand$intercept), t_time))
  if (Q=="Q_5P14") coeff_separate[["5P14"]]   <- as.matrix(cbind(data.frame(ligand=t_ligand$intercept), t_time))
  if (Q=="Q_6P4")  coeff_separate[["6P4"]]    <- as.matrix(cbind(data.frame(ligand=t_ligand$intercept), t_time))
  if (grepl("_vs_", Q))     coeff_separate[["RANTES"]] <- as.matrix(cbind(data.frame(ligand=t_ligand$intercept), t_time))
  if (grepl("_vs_5P14", Q)) coeff_separate[["5P14"]]   <- as.matrix(cbind(data.frame(ligand=t_ligand$intercept+t_ligand$liga_5P14), t_time+t_interact))
  if (grepl("_vs_6P4" , Q)) coeff_separate[["6P4"]]    <- as.matrix(cbind(data.frame(ligand=t_ligand$intercept+t_ligand$liga_6P4) , t_time+t_interact))
  rm(t_ligand, t_time, t_interact)

  ### for the interaction questions, calculate a "collated" fit and save model coefficients for plotting
  if (grepl("_vs_", Q)) {
    collated <- step06_ANOVA(data_wide_norm, semantics_var, model.on="ligand+time", stats.on="_term_", verbose=TRUE, splineDF=splineDF)
    coeff_collated <- list()
    t_ligand <- collated[["ligand"]] ; t_ligand <- t_ligand[ , !(colnames(t_ligand) %in% c("avg.quant", "F", "p.pVal", "p.pVal.adj")), drop=FALSE]
    t_time   <- collated[["time"]]   ; t_time   <- t_time  [ , !(colnames(t_time)   %in% c("intercept", "avg.quant", "F", "p.pVal", "p.pVal.adj")), drop=FALSE]
    coeff_collated[["RANTES"]]  <- as.matrix(cbind(data.frame(ligand=t_ligand$intercept), t_time))
    if (grepl("_vs_5P14", Q)) coeff_collated[["5P14"]] <- as.matrix(cbind(data.frame(ligand=t_ligand$intercept+t_ligand$liga_5P14), t_time))
    if (grepl("_vs_6P4" , Q)) coeff_collated[["6P4"]]  <- as.matrix(cbind(data.frame(ligand=t_ligand$intercept+t_ligand$liga_6P4) , t_time))
    rm(t_ligand, t_time)
  }

  ### for the interaction questions, compare response profiles using F-test & p-val approach
  if (grepl("_vs_", Q)) {
    separatefp <- separate$fit.properties ; collatedfp <- collated$fit.properties ;
    F <- ((collatedfp$ss-separatefp$ss)/(collatedfp$df.residual-separatefp$df.residual)) / (separatefp$ss/separatefp$df.residual)
    dfn <- collatedfp$df.residual-separatefp$df.residual ; dfd <- separatefp$df.residual ;
    prism <- data.frame(F=F, p.pVal.F=-log10(pf(F, dfn, dfd, lower.tail=FALSE)))
    rm(F, dfn, dfd, collatedfp, separatefp)
  }

  ### for the interaction questions, also compare response profiles using AICc approach:
  if (grepl("_vs_", Q)) {
    separatefp <- separate$fit.properties ; collatedfp <- collated$fit.properties ;
    dAIC <- separatefp$AICc - collatedfp$AICc ### corrected - for the peptide-level scenario; for the projected scenario the correction factor is negligible i.e. AICc still works
    prism <- cbind(prism, data.frame(p.prob.ratio=-dAIC/2*log10(exp(1)), prob.different=exp(-dAIC/2.)/(1.+exp(-dAIC/2)))) ### probability that the separate model is correct
    rm(dAIC, collatedfp, separatefp)
  }

  ### define semanticsX for plotting
  if (Q == "Q_CCL5") semanticsX <- semantics_var[semantics_var$ligand == "RANTES", ]
  if (Q == "Q_5P14") semanticsX <- semantics_var[semantics_var$ligand == "5P14"  , ]
  if (Q == "Q_6P4")  semanticsX <- semantics_var[semantics_var$ligand == "6P4"   , ]
  if (!grepl("_vs_", Q))   semanticsX <- within(semanticsX, rm(plex))
  if ( grepl("_vs_", Q)) { semanticsX <- semantics_var ; semanticsX$plex <- semanticsX$ligand }
  semanticsX$condition <- semanticsX$time

  if (exists("stat.env")){ ## This will assign the "separate" variable with names unique per Q
    if (l_projectPeptidesOntoProteins) {projLevel <- "protein"} else {projLevel <- "peptide"}
    assign(paste(projLevel, Q, "separate", sep = "_"), separate, envir = stat.env)
    if (grepl("vs", Q)) {assign(paste(projLevel, Q, "prism", sep = "_"), prism, envir = stat.env)}
  }

  ############################## --PLOTTING-- #############################

  ### plot highest- and lowest-ranking peptides/proteins for each question
  if (l_plots) {
    for (l_decreasing in c(TRUE, FALSE)) {

      ### options for ranking hits
      if (Q %in% c("Q_CCL5", "Q_5P14", "Q_6P4")) opts <- "r.p"
      if (grepl("_vs_", Q)) { opts <- c("r.int", "r.F", "r.AIC") }

      for (opt in opts) {
        if (opt=="r.p")   { t <- separate[["time"]]        ; ixs <- order(t$p.pVal.adj,       decreasing=l_decreasing)[var_from:(var_from+nplot-1)] }
        if (opt=="r.int") { t <- separate[["ligand:time"]] ; ixs <- order(t$p.pVal.adj,       decreasing=l_decreasing)[var_from:(var_from+nplot-1)] }
        if (opt=="r.F")   {                                  ixs <- order(prism$p.pVal.F,     decreasing=l_decreasing)[var_from:(var_from+nplot-1)] }
        if (opt=="r.AIC") {                                  ixs <- order(prism$p.prob.ratio, decreasing=l_decreasing)[var_from:(var_from+nplot-1)] }

        ### get the spline coefficients, trim them to only indices that are going to be plotted
        coeff_ixs <- coeff_separate ; basis <- separate$spline.basis
        if (!is.null(coeff_ixs)) for (ix in 1:length(coeff_ixs)) coeff_ixs[[ix]] <- coeff_ixs[[ix]][ixs, , drop=FALSE]

        ### plot into a png file with a question- and opt-dependent name
        if (l_decreasing) cat(paste0("\n", Q, ", opt=", opt, ": ", sum(!is.na(t$p.pVal.adj) & t$p.pVal.adj>-log10(0.05)), " significant proteins, plotting\n"))
        g <- ProteomicsToolkit::plotPoint(data_wide_norm[ixs, ], semanticsX, maxPoints=nplot, ncol=ncol, basis=basis, coeff=coeff_ixs, legendPos="bottom")
        oufnm <- paste0(dirOut, "/peptides_", Q, "/", Q, "_bottom.png")
        if (l_projectPeptidesOntoProteins) { oufnm <- gsub("peptides_", "proteins_", oufnm) } ; if (l_decreasing) { oufnm <- gsub("_bottom", "_top", oufnm) }
        if (opt=="r.p")    oufnm <- gsub(".png", "_response.png"      , oufnm, fixed=TRUE)
        if (opt=="r.int" ) oufnm <- gsub(".png", "_interactCoeff.png" , oufnm, fixed=TRUE)
        if (opt=="r.F"   ) oufnm <- gsub(".png", "_Ftest.png"         , oufnm, fixed=TRUE)
        if (opt=="r.AIC" ) oufnm <- gsub(".png", "_AIC.png"           , oufnm, fixed=TRUE)
        if (!dir.exists(dirname(oufnm))) dir.create(dirname(oufnm), recursive=TRUE)

         if (exists("e_plotEnvironment") && is.environment(e_plotEnvironment)) {
          cat("Placing plots into e_plotEnvironment.\n")

          oufnm <- paste0("peptides_", Q, "_bottom.png")
          if (l_projectPeptidesOntoProteins) { oufnm <- gsub("peptides_", "proteins_", oufnm) } ; if (l_decreasing) { oufnm <- gsub("_bottom", "_top", oufnm) }
          oufnm <- paste0(from, "-", (from+nplot-1), "_", gsub(".png", paste0("_", opt), oufnm))
          assign(oufnm, g, envir=e_plotEnvironment)
        }else{
          ggplot2::ggsave(oufnm, g, width=ncol, height=ceiling(nplot/ncol), units="in", scale=3, limitsize=FALSE)
        }
      }
    }
  }

  if (l_clustering) { ### calculate or load the matrices for the different questions, then run clustering and EnAndOR analyses

      cutType <- "h_8" ### or not, need to test

    if (!grepl("_vs_", Q)) { ### Q_CCL5, Q_5P14, or Q_6P4 time course
      oufnm <- paste0(dirOut, "/07_DM_", Q, "_time_time.Rdata") ### this is dependent on model.on & stats.on, will change accordingly after testing them
      data_XY <- data_wide_norm ; semantics_XY <- semantics_var[semantics_var$ligand==gsub("CCL5", "RANTES", gsub("Q_", "", Q)), ]
      model.on <- "time" ; stats.on <- "time_" ### not at all sure about this, need to test
    } else { ### this is a response difference question, Q_CCL5_vs_5P14 or Q_CCL5_vs_6P4. This should NOT need any correction for the basal differences but who knows?
      oufnm <- paste0(dirOut, "/07_DM_", Q, "_condition_condition_blc.Rdata") ### blc stands for baseline-corrected
      data_XY <- data_wide_norm ; semantics_XY <- semantics_var
      ### remove "batch" (receptor) effects in data_XY before calculating the DM
      tmp <- step06_ANOVA(data_wide_norm, semantics_var, model.on="ligand*time", stats.on="_term_", verbose=TRUE, splineDF=NULL) ### run a two-way separate model without splines
      t_ligand <- tmp$ligand ; ixs <- !is.na(t_ligand[, gsub("Q_CCL5_vs_", "liga_", Q)]) ### offsets for 5P14 or 6P4 vs WT CCL5
      for (colix in 1:nrow(semantics_XY)) {
        if (semantics_XY$ligand[colix] != gsub("Q_CCL5_vs_", "", Q)) next ;
        data_XY[ixs, semantics_XY$sample[colix]] <- data_XY[ixs, semantics_XY$sample[colix]] - t_ligand[ixs, gsub("Q_CCL5_vs_", "liga_", Q)]
      }
      model.on <- "condition" ; stats.on <- "condition"
    }

    if (file.exists(oufnm))   { ### retrieve a previously calculated DM, together with data_XY and semantics_XY
      load(file=oufnm); } else {### calculate the DM and save, together with data_XY and semantics_XY
      ### the call below takes quite a while, do not re-run if there is no reason for it
      protein_pairs <- step07_distanceMatrix(data_XY, semantics_XY, model.on=model.on, stats.on=stats.on, splineDF=splineDF, nproc=8)
      save(list=c("protein_pairs", "data_XY", "semantics_XY", "model.on", "stats.on"), file=oufnm)
    }

    ### append the right pval calculated above to data_XY so that 07_hclust call does not do this
    if (grepl("basal", Q)) {
      data_XY$p.pVal.adj <- onewaycat$cond_deltaG_00$p.pVal.adj
    } else if (grepl("normalG|deltaG", Q)) { ### Q1_normalG or Q2_deltaG
      if (grepl("normalG", Q)) data_XY$p.pVal.adj <- separate_no_intercept[["cell_normalG:time"]]$p.pVal.adj
      if (grepl("deltaG", Q))  data_XY$p.pVal.adj <- separate_no_intercept[["cell_deltaG:time"]]$p.pVal.adj
    } else { ### Q4_responseDiff
      opt <- "r.int"
      if (opt=="r.int") data_XY$p.pVal.adj <- separate[["ligand:time"]]$p.pVal.adj
      if (opt=="r.F")   data_XY$p.pVal.adj <- prism$p.pVal.F
      if (opt=="r.AIC") data_XY$p.pVal.adj <- prism$p.prob.ratio
    }

    ### here run hclust and ontology enrichment
    source(paste0(dirCode, "07_hclustTimeCourses.R"))
    source(paste0(dirCode, "08_GSEA_Cluster_Batches.R"))

  } ### end calculate or load the matrices + run clustering and EnAndOR analyses

}

### clear the environment
nms <- c("basis", "coeff_collated", "coeff_separate_no_intercept", "coeff_separate")
nms <- c(nms, "semanticsX", "oufnm", "coeff_ixs", "g", "ix", "opt", "opts", "var_from", "ixs", "l_decreasing", "nplot", "basis", "dAIC", "t_cell_line")
nms <- c(nms, "colix", "tmp", "ncol")

if (!exists("l_exportVars") || !l_exportVars) {nms <- c(nms, "separate", "collated", "separate_no_intercept", "prism", "t", "coeff_separate")}
nms <- nms[nms %in% ls()] ; if (length(nms)>0) rm(list=nms) ; rm(nms)
