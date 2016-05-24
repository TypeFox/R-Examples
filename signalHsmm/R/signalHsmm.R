#' signalHsmm - prediction of signal peptides
#'
#' @details
#' Secretory signal peptides are short (20-30 residues) N-terminal amino acid sequences 
#' tagging among others tag among others hormons, immune system proteins, structural 
#' proteins, and metabolic enzymes. They direct a protein to the endomembrane system and 
#' next to the extracellular localization. All signal peptides possess three distinct 
#' domains with variable length 
#' and characteristic amino acid composition. Despite their variability, signal peptides 
#' are universal enough to direct properly proteins in different secretory systems. 
#' For example, artifically introduced bacterial signal peptides can guide proteins in 
#' mammals and plants. 
#' 
#' The development of signalHsmm was funded by National Science Center 
#' (2015/17/N/NZ2/01845).
#' 
#' @description Using hidden semi-Markov models as a probabilistic framework, signalHsmm 
#' is new, highly accurate signal peptide predictor for eukaryotic proteins. 
#' @importFrom seqinr read.fasta a
#' @importFrom shiny runApp
#' @importFrom utils tail
#' @importFrom stats sd median na.omit
#' @importFrom graphics plot axis text lines legend 
#' @docType package
#' @name signalHsmm
#' @examples
#' few_predictions <- run_signalHsmm(benchmark_dat[1:3])
#' #see all predictions
#' pred2df(few_predictions)
#' #summary one prediction
#' summary(few_predictions[[1]])
#' #plot one prediction
#' plot(few_predictions[[1]])
#' 
#' #have fun with GUI
#' \dontrun{
#' gui_signalHsmm()
#' }
NULL
