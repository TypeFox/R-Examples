#' Read in an AGS or FGS test file
#' 
#' Imports a TAB-delimited AGS or FGS file.
#' Checks the format of the given FGS file. Reads in unique gene/protein IDs from the given file \code{tbl} so that the created object can be directly submitted to \code{\link{nea.render}} and \code{\link{gsea.render}}. The function \code{import.gs} will be automatically called from \code{\link{nea.render}} if parameters \code{AGS} and/or \code{FGS} of the latter correspond to text files rather than lists.
#' 
#' @param tbl Input TAB-delimited FGS file.
#' @param Lowercase render gene/protein IDs lower-case. 
#' @param col.gene Number of the column with gene/protein IDs in the input file.
#' @param col.set Number of the column with set (sample, pathway etc.) IDs in the input file. 
#' Special options: 
#' \code{col.set = 0}: all gene/protein IDs in the file become members of a single FGS; 
#' \code{col.set < 0} each single gene/protein IDs becomes a separate FGS (ID of which equals the gene/protein ID).  
#' @param gs.type GS type, 'a' produces AGS, 'f' produces FGS. 
#' @return A list with entry names that correspond AGS/FGS IDs in the input file and gene/protein IDs as elements of each entry.
#' @seealso \code{\link{mutations2ags}}, \link{samples2ags}
#' @examples
#' data(can.sig.go)
#' fpath <- can.sig.go
#' fgs.list <- import.gs(fpath)
#' summary(fgs.list)
#'
#' @importFrom utils read.table
#' @export



import.gs <- function(tbl, Lowercase = 1, col.gene = 2, col.set = 3, gs.type = '') {
#if col.set = 0, then a single FGS list is created
#if col.set < 0, then each single gene becomes an FGS
if (is.null(tbl)) {stop("No FGS file name given...");}
if (is.data.frame(tbl)){ 
  f1 <- tbl
} else {
f1 <- read.table(tbl, row.names = NULL, header = FALSE, sep = "\t", quote = "", dec = ".", na.strings = "", skip = 0, colClasses="character");
}
if (is.null(f1)) {stop(paste("Not a proper geneset file: ", tbl, "...", sep=" "));}
if (nrow(f1) < 1) {stop(paste("Not a proper geneset file: ", tbl, "...", sep=" "));}
if (length(unique(f1[,col.gene])) < 2) {stop(paste("Multiple gene IDs are not found: check parameter 'col.gene' ", tbl, "...", sep=" "));}

for (i in 1:ncol(f1)) {
if (Lowercase > 0 ) { 
f1[,i] <- tolower(f1[,i]);
}}
gs.list <- NULL
if (col.set > 0) {
for (f2  in unique(f1[,col.set])) {
gs.list[[f2]] <- as.vector(unique(f1[which(f1[,col.set] == f2),col.gene]));
}} else {
if (col.set < 0) {
for (f2  in unique(f1[,col.gene])) {
gs.list[[f2]] <- as.vector(c(f2));
}} else {
gs.list[[paste('users_single_', gs.type, 'gs',sep="")]] <- as.vector(unique(f1[,col.gene]));
}}
gs.list = as.list(gs.list)
 return(gs.list);
}
