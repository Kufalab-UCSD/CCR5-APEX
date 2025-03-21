suppressMessages(require(clusterProfiler))
suppressMessages(require(ReactomePA))
### genome-wide annotations for human, mouse, and rat require the following:
suppressMessages(library(org.Hs.eg.db)) # to install: BiocManager::install("org.Hs.eg.db")
suppressMessages(library(org.Mm.eg.db)) # to install: BiocManager::install("org.Mm.eg.db")
suppressMessages(library(org.Rn.eg.db)) # to install: BiocManager::install("org.Rn.eg.db")
suppressMessages(library(org.Ss.eg.db)) # to install: BiocManager::install("org.Ss.eg.db")
suppressMessages(require(Kmisc)) ### if not installed, run devtools::install_github("kevinushey/Kmisc")
suppressMessages(require(stats))
#set.seed(7)

ID_to_ID <- function(list_ID, fromType='UNIPROT', toType='ENTREZID') {
  result <- vector("list", length(list_ID)) ; result <- setNames(result, list_ID)
  if (!(fromType %in% c('UNIPROT', 'ENSEMBL', 'SYMBOL', 'ENTREZID'))) { cat(paste0("\nError: cannot yet map from type '", fromType, "'\n")) ; return(result) }
  if (!(toType   %in% c('UNIPROT', 'ENSEMBL', 'SYMBOL', 'ENTREZID'))) { cat(paste0("\nError: cannot yet map to type '", toType, "'\n")) ; return(result) }

  ### guess the organism
  for (org in c('HUMAN', 'MOUSE', 'RAT', NA)) {
    if (org=='HUMAN') fromType.allowed.vals <- AnnotationDbi::keys(org.Hs.eg.db, keytype=fromType)
    if (org=='MOUSE') fromType.allowed.vals <- AnnotationDbi::keys(org.Mm.eg.db, keytype=fromType)
    if (org=='RAT')   fromType.allowed.vals <- AnnotationDbi::keys(org.Rn.eg.db, keytype=fromType)
    if (sum(unique(list_ID) %in% fromType.allowed.vals) > 0.75*length(unique(list_ID))) break
  }
  if (is.na(org)) { cat(paste0("\nError: could not map ", fromType, " IDs ", paste0(list_ID[1:min(10, length(list_ID))], collapse=", "), "...\n")) ; return(result) }

  tc <- tryCatch(
    expr = {
     if (org=='HUMAN') ou <- capture.output(newIDs <- AnnotationDbi::mapIds(org.Hs.eg.db, keys=unique(list_ID), column=toType, keytype=fromType), type="message")
     if (org=='MOUSE') ou <- capture.output(newIDs <- AnnotationDbi::mapIds(org.Mm.eg.db, keys=unique(list_ID), column=toType, keytype=fromType), type="message")
     if (org=='RAT')   ou <- capture.output(newIDs <- AnnotationDbi::mapIds(org.Rn.eg.db, keys=unique(list_ID), column=toType, keytype=fromType), type="message")
     if (!exists('newIDs')) { cat(paste0("\nError: could not map ", fromType, " IDs ", paste0(list_ID[1:min(10, length(list_ID))], collapse=", "), "... to ", toType, "\n")) ; return(result) }
     return(plyr::mapvalues(list_ID, from=names(newIDs), to=newIDs, warn_missing=FALSE))
    }, error = function(e)   { cat(paste0("\nError: ", e)) ; return("error") }, warning = function(w) { cat(paste0("\nWarning: ", w)) ; return("warning") }, finally = { }
  )
  if (tc=="error") { cat(paste0("\nError: could not map ", fromType, " IDs ", paste0(list_ID[1:min(10, length(list_ID))], collapse=", "), "...\n")) ; return(result) }
}

########################################################################################
################################ DEFINE A FUNCTION #####################################
####### that calculates Reactome and GO term enrichment and overrepresentation #########
############### using the clusterProfiler and associated libraries #####################
######## https://yulab-smu.top/biomedical-knowledge-mining-book/index.html #############
########################################################################################

step08_ontologyEnAndOR <- function(named.scores, score.cutoff, term_pValCutoff, l_prioritizeEnrichment=FALSE, l_multiACcheck=FALSE, verbose=TRUE, s_ont="Reactome", fromType="UNIPROT") {
  ### named.scores are, well, a named list of "scores"
  ### the names MUST be UniProt AC and they MUST be unique; not sure what to do if names (AC) are not unique
      ### 2024-09-19: well maybe no they don't have to be unique
  ### higher "scores" signify more interesting proteins
  ### scores can be -log10(p-val) with p-val describing variation across an arbitrary set of conditions
  ### scores can also be cluster membership, or differences in variation between conditions
  ### score.cutoff is where to trim for over-representation; is only used if enrichment did not work
  ### score.cutoff scale must be consistent with that of named.scores!
  ### term_pValCutoff <- 0.1   ### term_pValCutoff is active but needs to be provided explicitly as an argument
  ### l_prioritizeEnrichment specifies whether enrichment more valuable than over-representation and should it be tried first. In not, go straight for overrepresentation
  ### l_multiACcheck  <- FALSE is AL's argument controlling whether to split rows with |-separated AC's into multiple rows. clusterProfiler will only take the greatest value if this creates duplicate rows.
  ### l_printPlots    <- FALSE ### l_printPlots is retired from the function, plots are built and returned but not printed
  ### fromType        <- "Uniprot" ### Default to starting with uniprot ACs as IDs. Can use another. Use same arguments as bitr library.
  ### s_ont This variable passes which ontology database to use. ie c("Reactome", "GO.CC", "GO.BP", "GO.MF")

  if (is.data.frame(named.scores)) df08 <- named.scores else df08 <- data.frame(AC=names(named.scores), score=named.scores)

  result <-list(df08=df08, Enrichment=NA, Overrepresentation=NA, Enplot=NA, ORplot=NA)
  # 2023-04-03 used to be: result <-list(df08=df08, GO_enrichment=NA, GO_overrepresentation=NA, Reactome_enrichment=NA, Reactome_overrepresentation=NA, goEnplot=NA, goORplot=NA, ReEnplot=NA, ReORplot=NA)

  if (sum(c("AC", "score") %in% colnames(df08))<2)  { warning(paste0("Error (step08_ontologyEnAndOR): if named.scores is a data.frame, it needs to have columns 'AC' and 'score'; currently at least one is absent; skipping\n")) ; return(result) }
  s_ont.input <- s_ont ; s_ont <- s_ont.input[s_ont.input %in% c("Reactome", "GO.CC", "GO.BP", "GO.MF")]
  if (length(s_ont)< 1) { warning(paste0("Error (step08_ontologyEnAndOR): Unknown ontology", paste0(s_ont.input, collapse="|"), "; skipping")) ; return(result) ; }
  s_ont <- s_ont[1] ### 2023-04-03 for now, will be doing one ontology at a time - IK

  nms <- df08$AC ; df08$AC <- unlist(lapply(stringr::str_split(df08$AC, " "), "[[", 1)) 
  df08 <- cbind(data.frame(name=nms, ENTREZID=ID_to_ID(df08$AC, toType="ENTREZID"), SYMBOL=ID_to_ID(df08$AC, toType="SYMBOL")), df08) ; rm(nms)
  df08$score[is.na(df08$score)] <- min(df08$score[!is.na(df08$score)])-0.0001 ### prior to 2024-09-19, used to be df08 <- df08[!is.na(df08$score), ]
  df08 <- df08[order(df08$score, decreasing=TRUE), ] ; result$df08  <- df08
  tmp <- df08 |> dplyr::select(-name, -SYMBOL, -AC) |> dplyr::filter(!is.na(ENTREZID)) |> dplyr::group_by(ENTREZID) |> dplyr::slice_max(order_by=score, n=1) |> dplyr::distinct()  # 241030 for unique ENTREZID inputs to enrichment/overrepresentation analysis
  scores <- as.numeric(tmp$score) ; names(scores) <- as.character(tmp$ENTREZID) ; scores <- scores[order(scores, decreasing=TRUE)]  # 241030 sort scores
  scoreType <- c("std", "pos")[1]

### DEAD BLOCK: re-written on 2024-09-21. New code needs testing
#  if (FALSE) {  
#    cat("Pretending that the ENSEMBL ids are actually ACs.")
#    toType <- c("ENTREZID", "SYMBOL")[1]
#    tmp <- clusterProfiler::bitr(df08$AC, fromType="ENSEMBL", toType=c("ENTREZID"), OrgDb="org.Hs.eg.db")# ; colnames(tmp) <- c("AC", toType)
#    tmp <- tmp[!duplicated(tmp$ENSEMBL), ]
#    # if (verbose) {
#    #   cat(paste0("\nWarning (step08_ontologyEnAndOR()): ", sum(!(df08$AC %in% tmp$AC)), " AC numbers did not match any ", toType, ", removing")) ; tmp <- tmp[!is.na(tmp[, toType]), ]
#    #   cat(paste0("\nWarning (step08_ontologyEnAndOR()): ", length(unique(tmp$AC[duplicated(tmp$AC)])), " AC numbers matched two or more ", toType, "s"))
#    # }
#    df08$ENSEMBL <- df08$AC
#    df08 <- plyr::join(tmp, df08, by='ENSEMBL', type="left"); rm(tmp)
#    df08 <- df08[order(df08$score, decreasing=TRUE), ] ; df08 <- df08[!duplicated(df08[, "ENTREZID"]), ]
#    
#    if (TRUE) { ### append gene symbols so that we have a way to match AC to them
#      tmp <- clusterProfiler::bitr(df08[, toType], fromType="ENTREZID", toType="SYMBOL", OrgDb="org.Hs.eg.db") ; colnames(tmp) <- c(toType, "SYMBOL")
#      df08 <- plyr::join(df08, tmp, by=toType, type="left"); rm(tmp)
#    }
#  }else{
#    
#    #if (sum(duplicated(df08$AC)) > 0) { warning("Error (step08_ontologyEnAndOR()): the names of named.scores contain duplicates, make them unique before proceeding") ; return(result) ; } ### 2024-09-19 removed this requirement
#    if (l_multiACcheck) { df08 <- tidyr::separate_rows(df08, AC, sep = "\\|") ; if (verbose) cat("\nInfo (step08_ontologyEnAndOR()): rows with \"|\" in AC separated\n\n") }
#  
#    ### Match Uniprot AC to something that can be used by the overrepresentation and enrichment functions.
#    ### For ReactomePA::gsePathway, was only successful with ENTREZID so far
#    ### For clusterProfiler::gseGO, same thing, plus
#    ### Select the ID to map to from c("ACCNUM", "ALIAS", "ENSEMBL", "ENSEMBLPROT", "ENSEMBLTRANS", "ENTREZID", "ENZYME", "EVIDENCE", "EVIDENCEALL", "GENENAME", "GENETYPE", "GO",
#    ###                                "GOALL", "IPI", "MAP", "OMIM", "ONTOLOGY", "ONTOLOGYALL", "PATH", "PFAM", "PMID", "PROSITE", "REFSEQ", "SYMBOL", "UCSCKG", "UNIPROT")
#    ### This is the output of idType("org.Hs.eg.db") by the way
#    toType <- c("ENTREZID", "SYMBOL")[1]
#  
#    ### bitr function from clusterProfiler returns with columns named for the fromType and toType verbatim
#    tmp <- clusterProfiler::bitr(df08$AC, fromType=fromType, toType=toType, OrgDb="org.Hs.eg.db") ; colnames(tmp) <- c("AC", toType)
#    if (verbose) {
#      cat(paste0("\nWarning (step08_ontologyEnAndOR()): ", sum(!(df08$AC %in% tmp$AC)), " AC numbers did not match any ", toType, ", removing")) ; tmp <- tmp[!is.na(tmp[, toType]), ]
#      cat(paste0("\nWarning (step08_ontologyEnAndOR()): ", length(unique(tmp$AC[duplicated(tmp$AC)])), " AC numbers matched two or more ", toType, "s"))
#    }
#    if (sum(df08$AC %in% tmp$AC)<1) { warning("Error (step08_ontologyEnAndOR()): the names of named.scores do not look like UniProt ACs, make sure they are UniProt ACs before proceeding") ; return(result) ; }
#    df08 <- plyr::join(tmp, df08, by='AC', type="left"); rm(tmp)
#
#    ### For toType's that matched multiple ACs, retain the one with the highest score; remove the rest
#    ### 2024-09-19 removed this requirement
#    #if (verbose) cat(paste0("\nWarning (step08_ontologyEnAndOR()): ", length(unique(df08[duplicated(df08[, toType]), toType])), " ", toType, "s matched two or more ACs, the one with the higher score will be retained\n"))
#    #df08 <- df08[order(df08$score, decreasing=TRUE), ] ; df08 <- df08[!duplicated(df08[, toType]), ]
#  
#    if (toType == "ENTREZID") { ### append gene symbols so that we have a way to match AC to them
#      tmp <- clusterProfiler::bitr(df08[, toType], fromType=toType, toType="SYMBOL", OrgDb="org.Hs.eg.db") ; colnames(tmp) <- c(toType, "SYMBOL")
#      df08 <- plyr::join(df08, tmp, by=toType, type="left"); rm(tmp)
#    }
#  }

  ##################################################################################
  ######## Different kinds of overrepresentation and enrichment analyses ###########
  ##################################################################################

  minGSSize <- 5 ; maxGSSize <- 500 #Changed 230106
  ### clear the environment from the possible results of earlier runs
  nms <- c("GO_overrepresentation", "GO_enrichment", "Reactome_overrepresentation", "Reactome_enrichment", "goEnplot", "goORplot", "ReEnplot", "ReORplot") ; nms <- nms[nms %in% ls()] ; if (length(nms)>0) rm(list=nms) ; rm(nms)
  GO_enrichment <- NA ; GO_overrepresentation <- NA ; Reactome_enrichment <- NA ; Reactome_overrepresentation <- NA
  s_error <- ""

  if ("Reactome" %in% s_ont) { ### attempt Reactome enrichment and over-representation analyses
    ### ReactomePA::enrichPathway() is an overrepresentation using a *hypergeometric* model
    # Reactome_enrichment <- ReactomePA::gsePathway(scores, scoreType=scoreType, pvalueCutoff=0.1, by="DOSE", nPerm=3) ### this does not find any enriched terms
    if (l_prioritizeEnrichment) { Reactome_enrichment <- tryCatch(
      ReactomePA::gsePathway(scores, scoreType=scoreType, by="fgsea", pvalueCutoff=term_pValCutoff, minGSSize=minGSSize, maxGSSize=maxGSSize),
      error=function(cond) { s_error <<- paste0(s_error, "\nError (ReactomePA::gsePathway() for Reactome_enrichment): ", gsub("Error: ", "", cond)) ; return(NA) }
    ) } else { Reactome_enrichment <- NA }
    if (sum(grepl("enrichResult|gseaResult", is(Reactome_enrichment)))>0) {
      Reactome_enrichment <- clusterProfiler::setReadable(Reactome_enrichment, OrgDb=org.Hs.eg.db)
      Reactome_overrepresentation <- NA
    } else {
      if (l_prioritizeEnrichment & verbose) cat("\nInfo (step08_ontologyEnAndOR()): Reactome enrichment analysis unsuccessful; trying overrepresentation analysis next\n")
      Reactome_enrichment <- NA ; Reactome_overrepresentation <- NA
      while ((score.cutoff > -log10(0.05)) && (sum(scores>score.cutoff)<2)) { score.cutoff <- score.cutoff-log10(sqrt(5)) }
      if (sum(scores>score.cutoff)>0) { Reactome_overrepresentation <- tryCatch(
        ReactomePA::enrichPathway(names(scores)[scores>score.cutoff], organism="human", universe=names(scores), pvalueCutoff=term_pValCutoff, minGSSize=minGSSize, maxGSSize=maxGSSize),
        error=function(cond) { s_error <<- paste0(s_error, "\nError (ReactomePA::enrichPathway() for Reactome_overrepresentation): ", gsub("Error: ", "", cond)) ; return(NA) }
      ) }
      if (sum(grepl("enrichResult|gseaResult", is(Reactome_overrepresentation)))>0) {
        Reactome_overrepresentation <- clusterProfiler::setReadable(Reactome_overrepresentation, OrgDb=org.Hs.eg.db)
        Reactome_overrepresentation@result <- Reactome_overrepresentation@result[Reactome_overrepresentation@result$p.adjust<term_pValCutoff, ]
        if (nrow(Reactome_overrepresentation@result)<1) Reactome_overrepresentation <- NA
      }
      if (sum(grepl("enrichResult|gseaResult", is(Reactome_overrepresentation)))>0) {
        message(paste0("Info (step08_ontologyEnAndOR()): Reactome overrepresentation analysis successful with ", nrow(Reactome_overrepresentation@result), " terms identified"))
      } else { message("Info (step08_ontologyEnAndOR()): Reactome overrepresentation analysis failed") ; Reactome_overrepresentation <- NA }
    }

    # Other arguments of ReactomePA::gsePathway:
    # * organism="human"
    # * exponent=1 # "weight of each step" - whatever this means
    # * minGSSize=10, maxGSSize=500 # "minimal/maximal size of each geneSet for analyzing"
    # * nPerm # permutation numbers, described in documentation but absent from the actual function. Apparently required when by="DOSE"
    # * eps=1e-10 # the opposite, not mentioned in the documentation but listed in the function
    # * pvalueCutoff=0.05 # for adjusted significance of term enrichment
    # * pAdjustMethod="BH"
    # * verbose=TRUE
    # * seed=FALSE
    # * by="fgsea" # "Users can use the GSEA algorithm implemented in DOSE or fgsea by specifying the parameter by="DOSE" or by="fgsea". By default, the fgsea method will be used since it is much more faster."

    # Reactome_enrichment@result is a data frame of terms. Among other columns, each term has:
    # * leading_edge "...reports:
    #   - Tags to indicate the percentage of genes contributing to the enrichment score,
    #   - List to indicate where in the list the enrichment score is attained and
    #   - Signal for enrichment signal strength.")
    # * core_enrichment ("core enriched genes... that contribute to the enrichment."

    ### Reactome enrichment/overrepresentation structures cannot be simplified despite having a hierarchical structure...
    result$Enrichment <- Reactome_enrichment
    result$Overrepresentation <- Reactome_overrepresentation
  } ### end Reactome enrichment and over-representation analyses

  if (sum(c("GO.CC", "GO.BP", "GO.MF") %in% s_ont)>0) { ### attempt GO enrichment and over-representation analyses
    #GO_enrichment <- clusterProfiler::gseGO(scores, ont="CC"/"BP", OrgDb=org.Hs.eg.db, by="DOSE", nPerm=5, keyType="ENTREZID") ### this does not find any enriched terms either
    s_GO <- unlist(lapply(stringr::str_split(s_ont, "[.]"), "[[", 2)) ; s_GO <- s_GO[1]
    if (l_prioritizeEnrichment) { GO_enrichment <- tryCatch(
      clusterProfiler::gseGO(scores, ont=s_GO, OrgDb=org.Hs.eg.db, by="fgsea", keyType="ENTREZID", pvalueCutoff=term_pValCutoff, minGSSize=minGSSize, maxGSSize=maxGSSize),
      error=function(cond) { s_error <<- paste0(s_error, "\nError (clusterProfiler::gseGO() for GO_enrichment): ", gsub("Error: ", "", cond)) ; return(NA) }
    ) } else { GO_enrichment <- NA }
    if (sum(grepl("enrichResult|gseaResult", is(GO_enrichment)))>0) {
      GO_enrichment <- clusterProfiler::setReadable(GO_enrichment, OrgDb=org.Hs.eg.db)
      GO_overrepresentation <- NA
    } else {
      if (l_prioritizeEnrichment & verbose) cat(paste0("\nInfo (step08_ontologyEnAndOR()): GO enrichment analysis (", s_GO, ") unsuccessful; trying overrepresentation analysis next\n"))
      GO_enrichment <- NA ; GO_overrepresentation <- NA
      while ((score.cutoff > -log10(0.05)) && (sum(scores>score.cutoff)<2)) { score.cutoff <- score.cutoff-log10(sqrt(5)) }
      if (sum(scores>score.cutoff)>0) { GO_overrepresentation <- tryCatch(
        clusterProfiler::enrichGO(names(scores)[scores>score.cutoff], ont=s_GO, OrgDb=org.Hs.eg.db, universe=names(scores), pvalueCutoff=term_pValCutoff, minGSSize=minGSSize, maxGSSize=maxGSSize),
        error=function(cond) { s_error <<- paste0(s_error, "\nError (clusterProfiler::enrichGO() for GO_overrepresentation): ", gsub("Error: ", "", cond)) ; return(NA) }
      ) }
      if (sum(grepl("enrichResult|gseaResult", is(GO_overrepresentation)))>0) {
        GO_overrepresentation <- clusterProfiler::setReadable(GO_overrepresentation, OrgDb=org.Hs.eg.db)
        GO_overrepresentation@result <- GO_overrepresentation@result[GO_overrepresentation@result$p.adjust<term_pValCutoff, ] ; if (nrow(GO_overrepresentation@result)<1) GO_overrepresentation <- NA
      }
      if (sum(grepl("enrichResult|gseaResult", is(GO_overrepresentation)))>0) {
        message(paste0("Info (step08_ontologyEnAndOR()): GO overrepresentation analysis (", s_GO, ") successful with ", nrow(GO_overrepresentation@result), " terms identified"))
      } else { message(paste0("Info (step08_ontologyEnAndOR()): GO overrepresentation analysis (", s_GO, ") failed")) ; GO_overrepresentation <- NA }
    }
    ### simplify GO enrichment/overrepresentation results
    term_sim_method <- c("JC", "Resnik", "Lin", "Rel", "Jiang")[1] ### for everything except JC, GOSemSim::godata() call is needed
    if (sum(grepl("enrichResult|gseaResult", is(GO_enrichment)))>0 && sum(GO_enrichment@result$p.adjust<term_pValCutoff)>0) { ### simplify GO enrichment structures
      term_sim <- enrichplot::pairwise_termsim(GO_enrichment, method=term_sim_method)
      GO_enrichment <- clusterProfiler::simplify(GO_enrichment, cutoff=0.7, by="p.adjust", select_fun=min) ### not sure about the min/max function here
    }
    if (sum(grepl("enrichResult|gseaResult", is(GO_overrepresentation)))>0 && sum(GO_overrepresentation@result$p.adjust<term_pValCutoff)>0) { ### simplify GO overrepresentation structures
      term_sim <- enrichplot::pairwise_termsim(GO_overrepresentation, method=term_sim_method)
      GO_overrepresentation <- clusterProfiler::simplify(GO_overrepresentation, cutoff=0.7, by="p.adjust", select_fun=min) ### not sure about the min/max function here
      GO_overrepresentation@result <- GO_overrepresentation@result[GO_overrepresentation@result$p.adjust<term_pValCutoff, ] ; if (nrow(GO_overrepresentation@result)<1) GO_overrepresentation <- NA
    }
    result$Enrichment <- GO_enrichment
    result$Overrepresentation <- GO_overrepresentation
  } ### done GO enrichment and over-representation analyses

  ### 2023-04-03 used to be (but is no longer necessary):
  ### update the result
  # result$df08 <- df08
  # result$GO_enrichment <- GO_enrichment
  # result$GO_overrepresentation <- GO_overrepresentation
  # result$Reactome_enrichment <- Reactome_enrichment
  # result$Reactome_overrepresentation <- Reactome_overrepresentation

  #return(result) ### debug return

  ##############################################################################################
  ############### Different ways of visualizing the results of the analysis ####################
  ## see https://bioconductor.org/packages/release/bioc/manuals/enrichplot/man/enrichplot.pdf ##
  ############################ can also output the graphs to a file ############################
  ##############################################################################################

  #oufnm <- paste0(dirOut, "/images/Reactome_enrichment.png")
  #cat(paste0("\nSaving plot to ", oufnm)) ; png(filename=oufnm, width=1200, height=1200)
  #
  # print(enrichplot::dotplot(Reactome_enrichment))
  # print(enrichplot::heatplot(Reactome_enrichment))
  # print(enrichplot::cnetplot(Reactome_enrichment, categorySize="pvalue", foldChange=scores, max.overlaps=10000)) ### this does not show edges in the R session but outputs them to afile
  # print(enrichplot::ridgeplot(Reactome_enrichment))
  # print(enrichplot::upsetplot(Reactome_enrichment))
  #
  #dev.off()

  if (FALSE) { ### term  enrichment / running score, several options
    enrichplot::gseaplot (Reactome_enrichment, geneSetID=1, by="runningScore", title=Reactome_enrichment$Description[1])
    enrichplot::gseaplot (Reactome_enrichment, geneSetID=1, by="preranked", title=Reactome_enrichment$Description[1])
    enrichplot::gseaplot (Reactome_enrichment, geneSetID=1, title=Reactome_enrichment$Description[1])
    enrichplot::gseaplot2(Reactome_enrichment, geneSetID=1:5, title="first five terms", pvalue_table=TRUE)
  }

  ### PMC stats of enrichment terms
  # p1 <- enrichplot::pmcplot(Reactome_enrichment$Description[1:5], 2010:2020, proportion=TRUE)
  # p2 <- enrichplot::pmcplot(Reactome_enrichment$Description[1:5], 2010:2020, proportion=FALSE)
  # cowplot::plot_grid(p1, p2, ncol=2)

  ### A graph of term relatedness for things that did not fail to calculate
  ### used to be constructed by enrichMap() but "enrichMap() is deprecated and is no longer part of the DOSE package"
  ### The way to achieve the same now is this:

  ##############################################
  ### generate plots and append terms to df08
  ##############################################

  term_sim_method <- c("JC", "Resnik", "Lin", "Rel", "Jiang")[1] ### for everything except JC, GOSemSim::godata() call is needed

  ### GO_enrichment
  if (sum(grepl("enrichResult|gseaResult", is(GO_enrichment)))>0 && sum(GO_enrichment@result$p.adjust<term_pValCutoff)>0) {
    cat("\nInfo (step08_ontologyEnAndOR()): generating plots for GO_enrichment...\n")
    #goEnplot <- enrichplot::emapplot(enrichplot::pairwise_termsim(GO_enrichment, method=term_sim_method), max.overlaps=10000) ### default max.overlaps=10
    goEnplot <- tryCatch(enrichplot::emapplot(enrichplot::pairwise_termsim(GO_enrichment, method=term_sim_method)), error=function(cond) { s_error <<- paste0(s_error, "\nError (enrichplot::emapplot() for GO_enrichment): ", gsub("Error: ", "", cond)) ; return(NA) } )
    result$Enplot <- goEnplot
    cat("Info (step08_ontologyEnAndOR()): assigning terms to proteins...\n") ; termcolnm <- paste0(s_ont, ".terms") ; df08[, termcolnm] <- ""
    for (i in 1:nrow(GO_enrichment@result)) {
      proteins_in_term <- unlist(stringr::str_split(GO_enrichment@result$core_enrichment[i], "/")) ; ixs <- df08$SYMBOL %in% proteins_in_term
      df08[ixs, termcolnm] <- paste0(df08[ixs, termcolnm], "; ", GO_enrichment@result$Description[i])
    }
    df08[, termcolnm] <- gsub("^;", "", df08[, termcolnm]) ; df08[, termcolnm] <- gsub("^ ", "", df08[, termcolnm])
  }

  ### GO_overrepresentation
  if (sum(grepl("enrichResult|gseaResult", is(GO_overrepresentation)))>0 && sum(GO_overrepresentation@result$p.adjust<term_pValCutoff)>0) {
    cat("\nInfo (step08_ontologyEnAndOR()): generating plots for GO_overrepresentation...\n")
    #goORplot <- enrichplot::emapplot(enrichplot::pairwise_termsim(GO_overrepresentation, method=term_sim_method), max.overlaps=10000) ### default max.overlaps=10
    goORplot <- tryCatch(enrichplot::emapplot(enrichplot::pairwise_termsim(GO_overrepresentation, method=term_sim_method)), error=function(cond) { s_error <<- paste0(s_error, "\nError (enrichplot::emapplot() for GO_overrepresentation): ", gsub("Error: ", "", cond)) ; return(NA) } )
    result$ORplot <- goORplot
    cat("Info (step08_ontologyEnAndOR()): assigning terms to proteins...\n") ; termcolnm <- paste0(s_ont, ".terms") ; df08[, termcolnm] <- ""
    for (i in 1:nrow(GO_overrepresentation@result)) {
      proteins_in_term <- unlist(stringr::str_split(GO_overrepresentation@result$geneID[i], "/")) ; ixs <- df08$SYMBOL %in% proteins_in_term
      df08[ixs, termcolnm] <- paste0(df08[ixs, termcolnm], "; ", GO_overrepresentation@result$Description[i])
    }
    df08[, termcolnm] <- gsub("^;", "", df08[, termcolnm]) ; df08[, termcolnm] <- gsub("^ ", "", df08[, termcolnm])
  }

  ### Reactome_enrichment
  if (sum(grepl("enrichResult|gseaResult", is(Reactome_enrichment)))>0 && sum(Reactome_enrichment@result$p.adjust<term_pValCutoff)>0) {
    cat("\nInfo (step08_ontologyEnAndOR()): generating plots for Reactome_enrichment...\n")
    #ReEnplot <- enrichplot::emapplot(enrichplot::pairwise_termsim(Reactome_enrichment, method=term_sim_method), max.overlaps=10000) ### default max.overlaps=10
    ReEnplot <- tryCatch(enrichplot::emapplot(enrichplot::pairwise_termsim(Reactome_enrichment, method=term_sim_method)), error=function(cond) { s_error <<- paste0(s_error, "\nError (enrichplot::emapplot() for Reactome_enrichment): ", gsub("Error: ", "", cond)) ; return(NA) } )
    result$Enplot <- ReEnplot
    cat("Info (step08_ontologyEnAndOR()): assigning terms to proteins...\n") ; termcolnm <- paste0(s_ont, ".terms") ; df08[, termcolnm] <- ""
    for (i in 1:nrow(Reactome_enrichment@result)) {
      proteins_in_term <- unlist(stringr::str_split(Reactome_enrichment@result$core_enrichment[i], "/")) ; ixs <- df08$SYMBOL %in% proteins_in_term
      df08[ixs, termcolnm] <- paste0(df08[ixs, termcolnm], "; ", Reactome_enrichment@result$Description[i])
    }
    df08[, termcolnm] <- gsub("^;", "", df08[, termcolnm]) ; df08[, termcolnm] <- gsub("^ ", "", df08[, termcolnm])
  }

  ### Reactome_overrepresentation
  if (sum(grepl("enrichResult|gseaResult", is(Reactome_overrepresentation)))>0 && sum(Reactome_overrepresentation@result$p.adjust<term_pValCutoff)>0) {
    cat("\nInfo (step08_ontologyEnAndOR()): generating plots for Reactome_overrepresentation...\n")
    #ReORplot <- enrichplot::emapplot(enrichplot::pairwise_termsim(Reactome_overrepresentation, method=term_sim_method), max.overlaps=10000) ### default max.overlaps=10
    ReORplot <- tryCatch(enrichplot::emapplot(enrichplot::pairwise_termsim(Reactome_overrepresentation, method=term_sim_method)), error=function(cond) { s_error <<- paste0(s_error, "\nError (enrichplot::emapplot() for Reactome_overrepresentation): ", gsub("Error: ", "", cond)) ; return(NA) } )
    result$ORplot <- ReORplot
    cat("Info (step08_ontologyEnAndOR()): assigning terms to proteins...\n") ; termcolnm <- paste0(s_ont, ".terms") ; df08[, termcolnm] <- ""
    for (i in 1:nrow(Reactome_overrepresentation@result)) {
      proteins_in_term <- unlist(stringr::str_split(Reactome_overrepresentation@result$geneID[i], "/")) ; ixs <- df08$SYMBOL %in% proteins_in_term
      df08[ixs, termcolnm] <- paste0(df08[ixs, termcolnm], "; ", Reactome_overrepresentation@result$Description[i])
    }
    df08[, termcolnm] <- gsub("^;", "", df08[, termcolnm]) ; df08[, termcolnm] <- gsub("^ ", "", df08[, termcolnm])
  }

  ### update and return the result
  result$df08 <- df08
  result$error <- s_error
  #if (exists("goEnplot")) result$goEnplot <- goEnplot
  #if (exists("goORplot")) result$goORplot <- goORplot
  #if (exists("ReEnplot")) result$ReEnplot <- ReEnplot
  #if (exists("ReORplot")) result$ReORplot <- ReORplot

  return(result)
}

##################################################################################
##### create a temporary directory ###############################################
##################################################################################

# letters.digits <- c(letters, stringr::str_split("123456789", ""))
# s_tempDir <- dirname(normalizePath(tempdir()))
# while (TRUE) {
#   tmpdirnm <- paste0(c("step08_ontEnAndOR_tmp", letters[unlist(ramify::randi(length(letters), 1)[1, 1])], letters.digits[unlist(ramify::randi(length(letters.digits), 5)[, 1])]), collapse="")
#   tmpdirnm <- paste0(s_tempDir, "/", tmpdirnm)
#   if (!dir.exists(tmpdirnm)) { dir.create(tmpdirnm, recursive=TRUE) ; break ; }
# }
# cat("\nInfo (step08_ontologyEnAndOR): temporary directory set to", tmpdirnm, "\n")

# write("TMP = '<your-desired-tempdir>'", file=file.path(Sys.getenv('R_USER'), '.Renviron'))

##################################################################################
##### print the current system temporary directory ###############################
##################################################################################

# cat("\nInfo (step08_ontologyEnAndOR): Rtmp directory is", normalizePath(tempdir()), "\n")
