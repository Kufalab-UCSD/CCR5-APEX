# Pipeline for CCR5-APEX proximity labeling mass spectrometry data analysis



**General Workflow**: Data analysis is performed using R code written into multiple steps, as described below. Each step is run through sourced R scripts, organized into a master R script referred to as the "analyzeData" script, which aims to walk through data preparation from parsing to statistical analysis. The framework  of the "answer_questions" script is to build models and perform statistical analysis for ranking top hits/proteins to answer specific biological questions, then builds clusters for each questions.

# Part 1: Pre-processing
## Step 0: R packages and other prerequisites

A list of necessary libraries is below. 

- R libraries: `igraph`, `tidyr`, `dplyr`, `Rcpp`, `reshape2`, `ggplot2`, `httr`, `installr`, `readr`, `parallel`, `BiocManager`. All installed using `install.packages("<package>"))`.
- `limma` is also required and is installed using `BiocManager::install("<package>")`
- `ProteomicsToolkit` is in-house R library of functions for plotting and statistical analysis. It is required for running the analysis pipeline. The GitHub project needs to be cloned, compiled, and installed. It is built using the `roxygen2` library and on Windows machines an Rtools install is necessary for compiling R libraries. (https://github.com/Kufalab-UCSD/ProteomicsToolkit)
- The Contaminant Repository for Affinity Purification (CRAPome, https://reprint-apms.org/) must be preprocessed using tools in the `ProteomicsToolkit` library. The resulting processed contaminant list for proximity dependent biotinylation is stored in the `dirUniprot` local directory.

**Additional Notes**
- Paths for commonly used filesystem directories are assigned as strings to multiple variables. This includes paths for main system GitHub directory, GitHub repositories, UniProtKB cache, and a directory for writing output files to.
- Each step will write to disk, in `dirOut`, [step].Rdata checkpoint files. 

In case any are missed, please let us know through the repository so they can be included!

## Step 1: Parsing Data `01_readAndRename_2021_CCR5_APEX.R`
Reads the input spreadsheet provided by mass spectrometry (MS) platform and creates two standard data frames to assist in organizing subsequent steps. The `kwd` object defines the which version of output from MS software to use. Peptide spectrum matches (PSMs) are the lowest level outputted form of data, which can then be mapped to peptides and subsequently proteins. 

- `data` is the data frame that contains peptide metadata, samples as individual columns, and MS quants. 
- `semantics` is the data frame that outlines the experiment design. This is assists in easily accessing sample columns from the `data`. 
- Multiple PSMs for a peptide are summed. 


## Step 2: Mapping Peptides to UniProtKB `02_mapUniprot.R'

Peptides are matched to UniProtKB. The `get_uniprot` function from `ProteomicsToolkit` manages which version to use and can cache a species specific UniProt version locally. Peptide sequences are matched to primary accession numbers (ACs) if unique and unambiguous. Otherwise, sequences are mapped to isoform sequences, secondary ACs, or multiple proteins/ACs. 

#### Intermediate Step 1: Label proteins as potential contaminants.
- Requires preprocessed and compiled CRAPome data in the `dirUniprot` local directory. 
- Proteins are not immediately removed and instead are managed using a "badness" score. The number of peptides representing a protein is limited based on this score.

#### Intermediate Step 2: Non-specific peptides are grouped using a `|` symbol. 
- Can be turned on/off by setting `l_uniquePeptides` to true/false.


## Step 3: Transform MS Values `03_transformFilterBaselineCorrect.R`

MS values are transformed to a natural log scale. For an experimental batch, log quants are median centered. Baseline correction can be applied if the `l_baselineCorrect` is set to true. 
- By default, experiment sample distributions before and after re-scaling are plotted and saved to `dirOut`.



## Step 4: Remove Batch Effects `04_removeBatchEffects.R`

In multiplexed MS experiments, sample values correlate strongly within technical runs of the instruments. By default, plots before and after corrections are saved to disk. This script can utilize one or both of the following strategies to correct these batch effects:
- Experiment bridge based batch effect removal. Requires a "bridge" sample be included within each technical batch. 
- Linear model based batch effect removal. This approach utilizes the `limma` R library for removing batch effects. Does not require "bridge" samples.  

**Note:** Not required if there are no technical batches. Analyze distributions and correlations of samples according to how the experiment was performed and instruments were used.


## Step 5: Project Peptides onto Proteins `05_projectPeptidesToProteins.R`
This script will map peptides to proteins. The result will create two new objects:
- Data frames `data_wide` and `data_wide_norm` which contain unique proteins as rows. Individual peptides are assigned as sample columns for a protein. The `data_wide_norm` will normalize MS log quants within a protein row and is the recommended form to use. 
- `semantics_wide` data frame is the projected equivalent of `semantics`. Samples are expanded to include peptides. This form can be used to easily access sample columns from `data_wide` or `data_wide_norm`. 

# Part 2: Statistical Analysis
## Answer Biological Questions `answer_questions_2021_CCR5_APEX.R`
Once data is processed, it can now be analyzed statistically and with clustering and ontologies. Beyond Step 5, the workflow continues in the "answer_questions" script which is organized into specific and separate *biological questions* to guide statistical analysis. 
- First build multiple models of the data depending on the *question*. Time resolved data utilizes a fitted spline model. When applicable, interaction coefficients are used for multi-variable ANOVA in subsequent steps. 
- Run statistical analysis of proteins with the fitted model. `limma` functions are used for its implementation of moderated ANVOA. Additionally, Akaike Information Criterion (AIC) and F-test based metrics can be used as an alternative to multi-variable ANOVA. 
- Plot top and bottom ranked proteins by significance.
- Generate a distance matrix for hierarchical clustering using a *question* based model. 
- Run individual clusters through gene set enrichment analysis (GSEA) and/or overrepresentation analysis (ORA).

## Step 6: ANOVA `06_ANOVA_multi.R`
While this step was previously a script that could be sourced like earlier stages, the ANOVA script was later written as a function. When sourced, it will load the `step06_ANOVA` function into the R global environment. This function is used in test code after step 5 in the "analyzeData" script and in the "answer_questions" script. 

## Step 7: Protein Clustering
A series of scripts to accomplish the following:
- Load into the R global environment custom functions for generating a distance matrix for hierarchical clustering. This process is time and resource intensive relative to the previous steps and can utilize multiple system CPUs through the `parallel` library. 
- Generate clusters with a set cut point using the distance matrix. The cut point is set to the `cutType` object as a string to either a height *h* or number of clusters *k*.

## Step 8: GSEA & ORA for Clusters
This script will plot and output to `dirOut` a series of proteins and ontology analysis results. When run within the "answer_questions" script, the results are organized into specific *questions*. 
- By default, GSEA is prioritized.




