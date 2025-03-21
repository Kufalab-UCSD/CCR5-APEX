#library(openxlsx)
library(reshape2) ### for melt
#library(stringr) ### for str_extract and str_to_upper
library(plyr)

dirIn <- paste0(Gdrive, "kufalab_projects/active/CCR5_APEX/210830_TMT_data/") ;

########### read and standardize the semantics dataframe ###########

filenm = paste0(dirIn, "Annotation file.xlsx") ; semantics <- openxlsx::read.xlsx(filenm, sheet=1, colNames = TRUE)
colnames(semantics) <- gsub(".0", "", colnames(semantics), fixed=TRUE) ### 2021-09-02 read.xlsx converts 126 into 126.0, correcting for that
semantics <- reshape2::melt(t(semantics), as.is = TRUE); colnames(semantics) <- c("channel", "plex", "sample"); semantics <- semantics[semantics$channel != "Experiment", ]
semantics$condition <- "" ; semantics$replicate <- as.numeric(gsub("_", "", stringr::str_extract(semantics$sample, "_[0-9]$"))); semantics$replicate[grep("5P14_10", semantics$sample)] <- 1:3
semantics$sample <- gsub("_[0-9]$", "", semantics$sample)

semantics$plex[semantics$plex == 1] <- 0 ; semantics$plex[semantics$plex == 2] <- 1 ; semantics$plex[semantics$plex == 0] <- 2 ### 2021-09-26 tentatively swap plexes

semantics$cell_line <- "CCR5"
semantics$cell_line[semantics$sample == "Cyto_APEX"] <- "Cyto"
semantics$cell_line[semantics$sample == "PM_APEX"] <- "PM"
semantics$cell_line[semantics$sample == "Endo_APEX"] <- "Endo"
semantics$cell_line[semantics$sample == "bridge"] <- "bridge"
semantics$cell_line[semantics$sample == "Blank"] <- "blank"
semantics$cell_line[grep("No_H2O2", semantics$sample)] <- "noH2O2"
semantics$ligand <- "none" ; semantics$time <- "00"

ixs <- grep("_[0-9]{1,2}m$", semantics$sample)
semantics$ligand[ixs] <- unlist(lapply(strsplit(semantics$sample [ixs], "_"), "[[", 1))
semantics$time[ixs] <- unlist(lapply(strsplit(semantics$sample [ixs], "_"), "[[", 2))
semantics$time <- gsub("m", "", semantics$time) ; semantics$time[semantics$time == "3"] <- "03"

semantics$sample <- paste(semantics$cell_line, semantics$ligand, semantics$time, semantics$replicate, sep="_")
semantics$condition <- paste(semantics$cell_line, semantics$ligand, semantics$time, sep="_")
semantics$plex_channel <- paste(semantics$plex, semantics$channel, sep = '_')

################### READ AND PARSE THE DATA ########################

if (kwd == "peptides") { filenm <- paste0(dirIn, "Handel_APEX_TMT16_PeptideGroups.txt") }
if (kwd == "proteins") { filenm <- paste0(dirIn, "Handel_APEX_TMT16_Proteins.txt") }
if (kwd == "PSMs") {  filenm <- paste0(dirIn, "Handel_APEX_TMT16_PSMs_Excel.xlsx") ; df <- openxlsx::read.xlsx(filenm, sheet=1, colNames = TRUE) ;
} else { df <- read.table(filenm, sep="\t", header=TRUE) }

### extract quant columns and 1-2 annotation columns only
if (kwd == "peptides") { colnms <- c("Master.Protein.Accessions", "Annotated.Sequence", colnames(df)[grep("Abundances..Grouped...", colnames(df), fixed=TRUE)]) }
if (kwd == "proteins") { colnms <- c("Accession", colnames(df)[grep("Abundances..Grouped...", colnames(df), fixed=TRUE)]) } ### can also save Description and parse out GN from it
if (kwd == "PSMs")     { colnms <- c("Master.Protein.Accessions", "Annotated.Sequence", "File.ID", colnames(df)[grep("Abundance:.", colnames(df), fixed=TRUE)]) }
data <- df[, colnms]

if (kwd == "PSMs") { colnms <- gsub("..", "_", gsub("Abundance:.", "ch", colnames(data)), fixed = TRUE) ;
} else { colnms <- gsub("..", "_", gsub("Abundances..Grouped...F", "", colnames(data)), fixed = TRUE) }
colnms[1] <- "AC" ; if (kwd != "proteins") { colnms[2] <- "sequence" } ; if (kwd == "PSMs") { colnms[3] <- "plex" ; colnames(data) <- colnms ; }

### fake a PSM table pivot by duplicating columns, then sum up PSMs
if (kwd == "PSMs") {
  ### fake pivot
  data <- data[rowSums(!is.na(data[, grep('ch1', colnames(data))])) > 0, ]
  chnms <- colnames(data)[grep('ch1', colnames(data))] ;
  plexes <- unique(data$plex) ; plexes <- plexes[order(plexes)]
  for (plex in plexes) { tmp <- data[, chnms] ; tmp[data$plex != plex, ] <- NA ; colnames(tmp) <- gsub("ch", paste0(gsub("F", "", plex), "_"), colnames(tmp)) ; data <- cbind(data, tmp) ; rm(tmp) }
  data <- data[, !(colnames(data) %in% chnms)] ; data <- within(data, rm(plex)) ;
  data <- cbind(data.frame(pepID=stringr::str_to_upper(paste0(data$AC, "_", data$sequence))), data) ### assuming that in annotated sequences upper case and lower case letters mean the same
  if (TRUE) { ### add up PSMs to the peptide level
    data <- aggregate(data[, semantics$plex_channel], by=list(pepID=data$pepID), FUN="sum", na.rm=TRUE) #https://stackoverflow.com/questions/1660124/how-to-sum-a-variable-by-group
    data <- cbind(data.frame(AC = unlist(lapply(strsplit(data$pepID, "_"), "[[", 1)), sequence = unlist(lapply(strsplit(data$pepID, "_"), "[[", 2))), data[, semantics$plex_channel])
  } else { data <- within(data, rm(pepID))}
  colnms <- colnames(data) ; data[data == 0] <- NA
}

### rename data columns so that they agree with the semantics dataframe
colnms <- plyr::mapvalues(colnms, from = semantics$plex_channel, to = semantics$sample) ### this won't do anything to PSM colnms
colnames(data) <- colnms
data$AC <- gsub(" ", "", data$AC)

rm(df)

save(list=c("data", "semantics"), file=paste(dirOut, '01_readAndRename.Rdata', sep="/"))
