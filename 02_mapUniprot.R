if (!exists("uniprot_HUMAN") | !exists("sec_ac") | !exists("iso")) {
  if (!exists("dirUniprot")) {dirUniprot <- "~/Documents/uniprot/"}
  speciesList = c("HUMAN", "MOUSE"); if (grepl("mansoni", tolower(dirOut))){speciesList = c("HUMAN", "SCHMA")}
  .GlobalEnv <- list2env(getUniprot(dirUniprot, querydate = format(Sys.time(), "%d-%b-%Y"), speciesList=speciesList), .GlobalEnv)
  rm(uniprot_MOUSE, speciesList) ;
}
uniprot_SPECIES <- uniprot_HUMAN;
if (grepl("mansoni", tolower(dirOut))){
  uniprot_SPECIES <- uniprot_SCHMA; unreviewed <- FALSE
  iso <- iso[grepl("_SCHMA", iso$Entry.name), ]
}else{
  iso <- iso[grepl("_HUMAN", iso$Entry.name), ]
}

if (!exists("experiment_list")) {experiment_list <- read.csv(paste0(dirCode, "kufalab_experiments.csv")); experiment_list <- experiment_list[experiment_list$experiment_type == "phospho", ]}
if (!exists("save_file")) { save_file <- TRUE }
if (!exists("rm_unreviewed")) { rm_unreviewed <- TRUE}

if (exists("kwd") && grepl("mansoni", kwd)) {
  cat("Not including APEX constructs.")
}else{
  APEXnPTX <- data.frame( AC=      c("PMAPEX",
                                     "EndoAPEX",
                                     "CytoAPEX",
                                     "BORPTX",
                                     "aa.rGai3.92.FlagAPEX2.C351I",
                                     "aa.rGai3.119.FlagAPEX2.C351I",
                                     "aa.NT.FlagAPEX2.FLCC.rGai3.C351I",
                                     "aa.NT.FLCC.rGai3.C351I"
                                     ),
                          Entry.name="", Status="reviewed", Protein.names="", Gene.names="",
                          Sequence=c("MGCIKSKGKDSADSAGSAGSAGVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYKSGLRTRAQASNSAAATMDYKDDDDKGKSYPTVSADYQDAVEKAKKKLRGFIAEKRCAPLMLRLAFHSAGTFDKGTKTGGPFGTIKHPAELAHSANNGLDIAVRLLEPLKAEFPILSYADFYQLAGVVAVEVTGGPKVPFHPGREDKPEPPPEGRLPDPTKGSDHLRDVFGKAMGLTDQDIVALSGGHTIGAAHKERSGFEGPWTSNPLIFDNSYFTELLSGEKEGLLQLPSDKALLSDPVFRPLVDKYAADEDAFFADYAEAHQKLSELGFADALQLPPLERLTLD",
                                     "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYKYSDLELKLRIEFESDAMFAAERAPDWVDAEECHRCRVQFGVVTRKHHCRACGQIFCGKCSSKYSTIPKFGIEKEVRVCEPCYEQLNKKAQGQGSESDAMFAAERAPDWVDAEECHRCRVQFGVVTRKHHCRACGQIFCGKCSSKYSTIPKFGIEKEVRVCEPCYEQLNKKAYVDGTAGPGSTGSAAATMDYKDDDDKGKSYPTVSADYQDAVEKAKKKLRGFIAEKRCAPLMLRLAFHSAGTFDKGTKTGGPFGTIKHPAELAHSANNGLDIAVRLLEPLKAEFPILSYADFYQLAGVVAVEVTGGPKVPFHPGREDKPEPPPEGRLPDPTKGSDHLRDVFGKAMGLTDQDIVALSGGHTIGAAHKERSGFEGPWTSNPLIFDNSYFTELLSGEKEGLLQLPSDKALLSDPVFRPLVDKYAADEDAFFADYAEAHQKLSELGFADALQLPPLERLTLD",
                                     "MVSKGEELFTGVVPILVELDGDVNGHKFSVSGEGEGDATYGKLTLKFICTTGKLPVPWPTLVTTLTYGVQCFSRYPDHMKQHDFFKSAMPEGYVQERTIFFKDDGNYKTRAEVKFEGDTLVNRIELKGIDFKEDGNILGHKLEYNYNSHNVYIMADKQKNGIKVNFKIRHNIEDGSVQLADHYQQNTPIGDGPVLLPDNHYLSTQSALSKDPNEKRDHMVLLEFVTAAGITLGMDELYKSGLRTRAQASNSAAATMDYKDDDDKGKSYPTVSADYQDAVEKAKKKLRGFIAEKRCAPLMLRLAFHSAGTFDKGTKTGGPFGTIKHPAELAHSANNGLDIAVRLLEPLKAEFPILSYADFYQLAGVVAVEVTGGPKVPFHPGREDKPEPPPEGRLPDPTKGSDHLRDVFGKAMGLTDQDIVALSGGHTIGAAHKERSGFEGPWTSNPLIFDNSYFTELLSGEKEGLLQLPSDKALLSDPVFRPLVDKYAADEDAFFADYAEAHQKLSELGFADALQLPPLERLTLD",
                                     "DDPPATVYRYDSRPPEDVFQNGFTAWGNNDNVLDHLTGRSCQVGSSNSAFVSTSSSRRYTEVYLEHRMQEAVEAERAGRGTGHFIGYIYEVRADNNFYGAASSYFEYVDTYGDNAGRILAGALATYQSEYLAHRRIPPENIRRVTRVYHNGITGETTTTEYSNARYVSQQTRANPNPYTSRRSVASIVGTLVRMAPVIGACMARQAESSEAMAAWSERAGEAMVLVYYESIAYSF",
                                     "MGCTLSAGDKAAVERSKMIDRNLREDGEKAAKEVKLLLLGAGESGKSTIVKQMKIIHEDGYSEDECKQYKVVVYSNTIQSIIAIIRAMGRLKSGSGSGDYKDDDDKGSGGKSYPTVSADYQDAVEKAKKKLRGFIAEKRCAPLMLRLAFHSAGTFDKGTKTGGPFGTIKHPAELAHSANNGLDIAVRLLEPLKAEFPILSYADFYQLAGVVAVEVTGGPKVPFHPGREDKPEPPPEGRLPDPTKGSDHLRDVFGKAMGLTDQDIVALSGGHTIGAAHKERSGFEGPWTSNPLIFDNSYFTELLSGEKEGLLQLPSDKALLSDPVFRPLVDKYAADEDAFFADYAEAHQKLSELGFADASGSIDFGEAARADDARQLFVLAGSAEEGVMTSELAGVIKRLWRDGGVQACFSRSREYQLNDSASYYLNDLDRISQTNYIPTQQDVLRTRVKTTGIVETHFTFKELYFKMFDVGGQRSERKKWIHCFEGVTAIIFCVALSDYDLVLAEDEEMNRMHESMKLFDSICNNKWFTDTSIILFLNKKDLFEEKIKRSPLTICYPEYTGSNTYEEAAAYIQCQFEDLNRRKDTKEVYTHFTCATDTKNVQFVFDAVTDVIIKNNLKEIGLY",
                                     "MGCTLSAGDKAAVERSKMIDRNLREDGEKAAKEVKLLLLGAGESGKSTIVKQMKIIHEDGYSEDECKQYKVVVYSNTIQSIIAIIRAMGRLKIDFGEAARADDARQLFVLAGSAEEGVMSGSGSGDYKDDDDKGSGGKSYPTVSADYQDAVEKAKKKLRGFIAEKRCAPLMLRLAFHSAGTFDKGTKTGGPFGTIKHPAELAHSANNGLDIAVRLLEPLKAEFPILSYADFYQLAGVVAVEVTGGPKVPFHPGREDKPEPPPEGRLPDPTKGSDHLRDVFGKAMGLTDQDIVALSGGHTIGAAHKERSGFEGPWTSNPLIFDNSYFTELLSGEKEGLLQLPSDKALLSDPVFRPLVDKYAADEDAFFADYAEAHQKLSELGFADASGSTSELAGVIKRLWRDGGVQACFSRSREYQLNDSASYYLNDLDRISQTNYIPTQQDVLRTRVKTTGIVETHFTFKELYFKMFDVGGQRSERKKWIHCFEGVTAIIFCVALSDYDLVLAEDEEMNRMHESMKLFDSICNNKWFTDTSIILFLNKKDLFEEKIKRSPLTICYPEYTGSNTYEEAAAYIQCQFEDLNRRKDTKEVYTHFTCATDTKNVQFVFDAVTDVIIKNNLKEIGLY",
                                     "MGCTLSAGDKAAVERSKMIDSGSGSGDYKDDDDKGSGGKSYPTVSADYQDAVEKAKKKLRGFIAEKRCAPLMLRLAFHSAGTFDKGTKTGGPFGTIKHPAELAHSANNGLDIAVRLLEPLKAEFPILSYADFYQLAGVVAVEVTGGPKVPFHPGREDKPEPPPEGRLPDPTKGSDHLRDVFGKAMGLTDQDIVALSGGHTIGAAHKERSGFEGPWTSNPLIFDNSYFTELLSGEKEGLLQLPSDKALLSDPVFRPLVDKYAADEDAFFADYAEAHQKLSELGFADASGSGSDYKDDDDKGSSGLCCTLSAGDKAAVERSKMIDRNLREDGEKAAKEVKLLLLGAGESGKSTIVKQMKIIHEDGYSEDECKQYKVVVYSNTIQSIIAIIRAMGRLKIDFGEAARADDARQLFVLAGSAEEGVMTSELAGVIKRLWRDGGVQACFSRSREYQLNDSASYYLNDLDRISQTNYIPTQQDVLRTRVKTTGIVETHFTFKELYFKMFDVGGQRSERKKWIHCFEGVTAIIFCVALSDYDLVLAEDEEMNRMHESMKLFDSICNNKWFTDTSIILFLNKKDLFEEKIKRSPLTICYPEYTGSNTYEEAAAYIQCQFEDLNRRKDTKEVYTHFTCATDTKNVQFVFDAVTDVIIKNNLKEIGLY",
                                     "MGCTLSAGDKAAVERSKMIDSGSGSGSDYKDDDDKGSSGLCCTLSAGDKAAVERSKMIDRNLREDGEKAAKEVKLLLLGAGESGKSTIVKQMKIIHEDGYSEDECKQYKVVVYSNTIQSIIAIIRAMGRLKIDFGEAARADDARQLFVLAGSAEEGVMTSELAGVIKRLWRDGGVQACFSRSREYQLNDSASYYLNDLDRISQTNYIPTQQDVLRTRVKTTGIVETHFTFKELYFKMFDVGGQRSERKKWIHCFEGVTAIIFCVALSDYDLVLAEDEEMNRMHESMKLFDSICNNKWFTDTSIILFLNKKDLFEEKIKRSPLTICYPEYTGSNTYEEAAAYIQCQFEDLNRRKDTKEVYTHFTCATDTKNVQFVFDAVTDVIIKNNLKEIGLY"
                                     )
                        )
  APEXnPTX$Entry.name <- APEXnPTX$AC ; APEXnPTX$Protein.names <- APEXnPTX$AC ; APEXnPTX$Gene.names <- APEXnPTX$AC
  uniprot_SPECIES <- rbind(uniprot_SPECIES, APEXnPTX) ; uniprot_SPECIES <- uniprot_SPECIES[!duplicated(uniprot_SPECIES$AC), ]
  APEXnPTX <- data.frame(AC=APEXnPTX$AC, primaryAC=APEXnPTX$AC) ; sec_ac <- rbind(sec_ac, APEXnPTX) ; rm(APEXnPTX) ; sec_ac <- sec_ac[!duplicated(sec_ac$AC), ]
}

################# split the multi-sequence rows (for phospho project) #################
if (experiment_list[experiment_list$kwd == kwd, "experiment_type"] == "phospho") {
  cat('The "kwd" variable matches a phosphoproteomics project from the experiment_list\nUsing phospho specific algorithm\n')
  cat("\n\nSplitting multi-sequence rows")
  npros <- data.frame(N=unlist(lapply(strsplit(gsub(" ", "", data$sequence), ";"), length))) ### how many sequences in each row?
  if (max(npros) == 1) {
    cat("\nEach row only has one AC. Skipping.")
  }else{
    data_orig <- data
    for (npro in 1:max(npros$N)) {
      cat(paste("\n  processing rows with >=", npro, "sequence values", sep=" "))
      ixs <- npros$N >= npro ; data_tmp <- data_orig[ixs, ]
      data_tmp$sequence <- unlist(lapply(strsplit(data_tmp$sequence  , ";"), "[[", npro))
      if (npro == 1) { data <- data_tmp } else { data <- rbind(data, data_tmp) }
      rm(data_tmp)
    }
    rownames(data) <- 1:nrow(data) ; rm(npros)
  }
}
################# split the multi-AC rows - first attempt ##########################

cat("\n\nSplitting multi-AC rows - first attempt")
if (experiment_list[experiment_list$kwd == kwd, "experiment_type"] != "phospho") {
  ### implemented 2022-12-18: circularly permute AC's in the data$AC column so that reviewed ACs come first - this will speed up the subsequent processing
  ### but this should only be done for non-phospho experiments because for phospho, pPos would need to be mutated similarly
  reviewed_primary_AC <- uniprot_SPECIES$AC[uniprot_SPECIES$Status == "reviewed"]
  npros <- data.frame(N=unlist(lapply(strsplit(gsub(" ", "", data$AC), ";"), length))) ### how many fields in each row?
  rot <- 0 ; ixs <- !is.na(data$AC) ### initialize with all TRUE's
  while (sum(ixs)>0) {
    rot <- rot+1
    AC1 <- unlist(lapply(strsplit(data$AC, ";"), "[[", 1)) ### the first AC in each row
    ixs <- !(AC1 %in% reviewed_primary_AC) ### for these row, the first AC is unreviewed; rotate them
    ixs <- ixs & (npros$N > rot)
    for (i in 1:nrow(data)) { if (ixs[i]) { data$AC[i] <- paste0(gsub(paste0(AC1[i], ";"), "", data$AC[i]), ";", AC1[i]) } }
  }
}

npros <- data.frame(N=unlist(lapply(strsplit(gsub(" ", "", data$AC), ";"), length))) ### how many fields in each row?
data_orig <- data
for (npro in 1:max(npros$N)) {
  cat(paste("\n  processing rows with >=", npro, "AC values", sep=" "))
  ixs <- npros$N >= npro ; data_tmp <- data_orig[ixs, ]
  data_tmp$AC <- unlist(lapply(strsplit(data_tmp$AC, ";"), "[[", npro))

  # for phospho, split multi-pPos rows as well
  if (experiment_list[experiment_list$kwd == kwd, "experiment_type"] == "phospho") data_tmp$pPos <- unlist(lapply(strsplit(as.character(data_tmp$pPos), ";"), "[[", npro))

  # replace secondary accession numbers by primary ones ; try matching to UniProt
  tmp <- plyr::join(data_tmp, sec_ac, by='AC', type="left") ; ixs <- !is.na(tmp$primaryAC) ; data_tmp$AC[ixs] <- tmp$primaryAC[ixs] ; ### substitute the primary accession numbers
  tmp <- plyr::join(data_tmp, uniprot_SPECIES, by='AC', type="left") ; data_tmp$reviewed <- tmp[, "Status"] ; ### try matching to UniProt and get status

  if (npro == 1) { data <- data_tmp } else {
    data_tmp <- data_tmp[!is.na(data_tmp$reviewed), ] ### on repeats, remove unmatched entries (preserve them on the first iteration though)
    data <- rbind(data, data_tmp)
  }
  rm(data_tmp)
}
rownames(data) <- 1:nrow(data) ; rm(tmp, npros) ;
if (exists("semantics")) { colnms <- colnames(data)[!(colnames(data) %in% semantics$sample)] ; data <- data[, c(colnms, semantics$sample)] }
data$reviewed[is.na(data$reviewed)] <- "no match"
cat(paste0("\n", sum(data$reviewed != "reviewed"), " rows remained unreviewed or unmatched"))

######### try to rescue unmatched and unreviewed peptides ###########################
uniprot_grep <- uniprot_SPECIES[uniprot_SPECIES$Status=="reviewed" | unreviewed, ] ; uniprot_grep$Sequence <- paste0("-", uniprot_grep$Sequence, "-")
iso_grep <- iso ; iso_grep$sequence <- paste0("-", iso_grep$sequence, "-")

data$grepSeq <- gsub("<([A-Z,-])>", "\\1", gsub(".", "", gsub("[", "<", gsub("]", ">", gsub("_", "", data$sequence), fixed=TRUE), fixed=TRUE), fixed=TRUE))
ixs <- grep("<[A-Z]{1,}-[A-Z]{1,}>", data$grepSeq) ; strsFrom <- stringr::str_extract(data$grepSeq[ixs], "<[A-Z]{1,}-[A-Z]{1,}>") ; strsTo <- gsub(">", "->", gsub("-", "", strsFrom))
if (length(ixs)>0) { ### this is to replace things like "<Y-R>" with "<YR->"
  for (i in 1:length(ixs)) { data$grepSeq[ixs[i]] <- gsub(strsFrom[i], strsTo[i], data$grepSeq[ixs[i]], fixed=TRUE) }
  rm (ixs, strsFrom, strsTo)
}
data$grepSeq <- gsub(">", "]", gsub("<", "[", data$grepSeq, fixed=TRUE), fixed=TRUE)
data$grepSeq <- gsub("; ", "|", data$grepSeq)
data$grepSeq <- toupper(data$grepSeq)
data$grepFixed <- !grepl("[", data$grepSeq, fixed=TRUE) & !grepl("|", data$grepSeq, fixed=TRUE) ### grepFixed indicates if the search can be performed w/o regex

cat("\n\nsearching for unmatched peptides in primary SwissProt sequences:")
ixs <- (1:nrow(data))[data$reviewed != "reviewed"] ; pb <- txtProgressBar(min=1, max=nrow(data), style=3)
for (i in ixs) {
  setTxtProgressBar(pb, i) ; ACs <- uniprot_grep$AC[grepl(data$grepSeq[i], uniprot_grep$Sequence, fixed=data$grepFixed[i])]
  if (length(ACs) > 0) { data$reviewed[i] <- "reviewed" ; data$AC[i] <- paste(ACs, collapse=";") }
}
close(pb) ; cat(paste0("\n", sum(data$reviewed == "unreviewed"), " peptides remained unreviewed, ", sum(data$reviewed == "no match"), " unmatched"))

cat("\n\nsearching for unmatched peptides in SwissProt isoform sequences:")
ixs <- (1:nrow(data))[data$reviewed != "reviewed"] ; pb <- txtProgressBar(min=1, max=nrow(data), style=3)
for (i in ixs) {
  setTxtProgressBar(pb, i) ; ACs <- unique(unlist(lapply(strsplit(iso_grep$AC[grepl(data$grepSeq[i], iso_grep$sequence, fixed=data$grepFixed[i])], "-"), "[[", 1)))
  if (length(ACs) > 0) { data$reviewed[i] <- "reviewed" ; data$AC[i] <- paste(ACs, collapse=";") }
}
close(pb) ; cat(paste0("\n", sum(data$reviewed == "unreviewed"), " peptides remained unreviewed, ", sum(data$reviewed == "no match"), " unmatched"))

rm(pb) ; rm(uniprot_grep, iso_grep) ; data <- within(data, rm(grepSeq, grepFixed))

############ split the multi-AC rows - second attempt after rematching ###############
######### this is needed because rematching has generated new multi-AC rows ##########
############## this time, also bring additional columns from UniProt #################

cat("\n\nSplitting multi-AC rows - second attempt after rematching")
npros <- data.frame(N=unlist(lapply(strsplit(gsub(" ", "", data$AC), ";"), length))) ### now, how many fields in each row?
data_orig <- data
for (npro in 1:max(npros$N)) {
  cat(paste0("\nprocessing rows with >= ", npro, " AC values "))
  ixs <- npros$N >= npro ; data_tmp <- data_orig[ixs, ]
  data_tmp$AC <- unlist(lapply(strsplit(data_tmp$AC, ";"), "[[", npro))

  # bring the protein names and gene names from UniProt
  tmp <- plyr::join(data_tmp, uniprot_SPECIES, by='AC', type="left") ; data_tmp$reviewed <- tmp[, "Status"] ;
  data_tmp$GN <- tmp[, "Gene.names"] ; data_tmp$protein <- tmp[, "Entry.name"] ; data_tmp$protein <- unlist(lapply(strsplit(data_tmp$protein, "_"), "[[", 1))

  if (npro == 1) { data <- data_tmp } else {
    data_tmp <- data_tmp[!is.na(data_tmp$reviewed), ] ### on repeats, remove unmatched entries (preserve them on the first iteration though)
    data <- rbind(data, data_tmp)
  }
  rm(data_tmp)
}
rownames(data) <- 1:nrow(data) ; rm(tmp, npros) ;
data$prot_seq <- paste0("unknown ", data$sequence) ; ixs <- !is.na(data$protein) ; data$prot_seq[ixs] <- paste0(data$protein[ixs], " ", data$sequence[ixs])
if ("m2z" %in% colnames(data)) { ### this is for Conce's exclusion list samples
  data$prot_seq <- paste0("unknown ", as.character(data$nMC), " ", as.character(data$m2z), " ", as.character(data$mh), " ", as.character(data$mhTheoretical), " ", data$sequence)
  ixs <- !is.na(data$protein) ;
  data$prot_seq[ixs] <- paste0(data$protein[ixs], " ", as.character(data$nMC[ixs]), " ", as.character(data$m2z[ixs]), " ", as.character(data$mh[ixs]), " ", as.character(data$mhTheoretical[ixs]), " ", data$sequence[ixs])
}
if ("multiplicity" %in% colnames(data)) data$prot_seq <- paste0(data$prot_seq, data$multiplicity)

if (exists("semantics")) { colnms <- colnames(data)[!(colnames(data) %in% semantics$sample)] ; data <- data[, c(colnms, semantics$sample)] }
if (rm_unreviewed){
  data$reviewed[is.na(data$reviewed)] <- "no match"
} else {
  data$reviewed[is.na(data$reviewed)] <- "unreviewed"
}
rm(data_orig)

################ match to identify phospho positions if need to #################
if (experiment_list[experiment_list$kwd == kwd, "experiment_type"] == "phospho") {
  cat("\n\nMatching to identify phospho positions\n")
  data$sequence <- stringr::str_to_upper(data$sequence) ; data$sequence <- gsub("[^A-Z]", "", data$sequence)

  data$prot_site <- NA ; ndata <- nrow(data) ; pb <- txtProgressBar(min=1, max=ndata, style=3)
  for (i in 1:ndata) {
    setTxtProgressBar(pb, i)
    s <- data$sequence[i] ; ixs <- grep(data$AC[i], uniprot_SPECIES$AC, fixed=TRUE)
    if (length(ixs) == 0) {
      cat(paste0("\n", data$AC[i], "(", data$protein[i], ") not found in human UniProt"))
    } else {
      S <- uniprot_SPECIES$Sequence[ixs[1]] ; aa <- substr(S, data$pPos[i], data$pPos[i]) ; prot <- strsplit(uniprot_SPECIES$Entry.name[ixs[1]], "_")[[1]][1]
      #if (grepl(s, S, fixed=TRUE) & (aa %in% c("S", "T", "Y"))) { ### strict criterion: ensure that the peptide sequence is identical to the UniProt sequence AND aa at position pPos is STY
      if (aa %in% c("S", "T", "Y")) {                              ### 2024-09-14 relaxed criterion: don't require that the peptide sequence is identical to the UniProt sequence, just ensure that aa at position pPos is STY
        data$prot_site[i] <- paste0(prot, "_", aa, data$pPos[i])
      } else {
        ixs <- grepl(data$AC[i], iso$AC, fixed=TRUE) & (substr(iso$sequence, data$pPos[i], data$pPos[i]) %in% c("S","T","Y")) & grepl(s, iso$sequence, fixed=TRUE)
        if (sum(ixs) == 0) {
          cat(paste0("\n", data$AC[i], "(", data$protein[i], "): sequence is unmatched or site is not STY in Uniprot and isoforms")) }
        else {
          S <- iso$sequence[ixs] [1] ; aa <- substr(S, data$pPos[i], data$pPos[i]) ; prot <- strsplit(iso$Entry.name[ixs] [1], "_")[[1]][1]
          data$prot_site[i] <- paste0(prot, "_", aa, data$pPos[i])
        } ### found in iso
      } ### sequence not matching or site not STY, search in iso
    } ### AC found in human UniProt
  } ### for i in 1:ndata
  close(pb)
  data <- within(data, rm(pPos))
}

############ remove duplicates of prot_seq and save ###############
data <- data[order(data$prot_seq, data$reviewed, data$AC), ] ; rownames(data) <- 1:nrow(data)
data <- data[!duplicated(data$prot_seq), ] ; data <- within(data, rm(prot_seq))

str <- ""
if (experiment_list[experiment_list$kwd == kwd, "experiment_type"] == "phospho") { ixs <- is.na(data$prot_site) } else { ixs <- data$reviewed=="no match" }
if (sum(ixs)<1000) { str <- paste0(": ", paste0(unique(paste0(data$AC[ixs], " (", data$protein[ixs], ")")), collapse=", ")) }
cat(paste0("\ndeleting ", sum(ixs), " unmatched peptide(s) from ", length(unique(data$AC[ixs])), " unique proteins", str))
data <- data[!ixs, ]

if (save_file) {
  if (exists("semantics")) {
    save(list=c("data", "semantics"), file=paste0(dirOut, "/02_mapUniprot.Rdata"))
  } else {
    save(list=c("data"), file=paste0(dirOut, "/02_mapUniprot.Rdata"))
  }
}
rm(save_file)
