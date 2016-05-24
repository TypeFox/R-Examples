#' Create AGS from a raw data matrix.
#'
#' This function creates a number of new AGSs from a given dataset, such as gene copy number, gene/protein expression, gene methylation etc. Such matrices M typically have size Ngenes x Nsamples, so that the current function returns a list of \code{length=ncol(M)}. The AGSs for each of the Nsamples are created with one of the five available methods (see parameter \code{method}).

#' @param method Method to select sample-specific genes. One of
#' \itemize{
#' \item "significant" : The default. Using one-sided z-test, selects genes values of which in the given sample \code{i} deviate from the mean over all the samples, requiring q-value (Benjamini-Hochberg FDR) be below  \code{cutoff.q}.
#' \item "topnorm"     : similar to "significant", i.e. calculates \code{(x[i] - mean(x)) / sd(x)} but does not evaluate significance. Instead, top N ranked genes \code{Ntop} are taken into AGS.
#' \item "top"         : similar to "topnorm", but \code{x[i] - mean(x)} is not divided with \code{sd(x)}. This might help to prioritize genes with higher \code{mean(x)} and ignore ones  with low signal. Consider also that AGSs from "top" overlap much more with each other than those from "topnorm", i.e. would be less sample-specific.
#' \item "toppos"     : similar to "top", but retrives only genes with positive values of \code{x[i] - mean(x)}. This might be useful when the gene expression values are small counts (such as in sincle-cell RNA sequencing), so that considering the left part of the distribution would not bring high-quality AGS.
#' \item "toprandom"     : generates lists of Ntop random genes for each AGS.
#' }
#' @param m0 input matrix.
#' @param col.mask To include only columns with IDs that contain the specified mask. Follows the regular expression synthax.
#' @param namesFromColumn Number of the column (if any) that contains the gene/protein names. Note that it is only necessary of the latter  are NOT the unique rownames of the matrix. This could be sometimes needed to be able to process redundant expression etc. profiles.
#' @param Lowercase Render gene/protein IDs lower-case.
#' @param cutoff.q cutoff value. Default 0.05. Mutually exclusive with "Ntop".
#' @param Ntop Number of top ranking genes to include into each sample-specific AGS. Mutually exclusive with "cutoff.q". A practically recommended value of Ntop could be in the range 30...300. Ntop>1000 might decrease the analysis specificity.
#'
#' @examples
#' data("fantom5.43samples", package="NEArender")
#' ags.list <- samples2ags(fantom5.43samples, cutoff.q = 0.01, method="significant")
#' @export



samples2ags <- function(m0, Ntop=NA, col.mask=NA, namesFromColumn=NA, method=c("significant", "top", "toppos", "topnorm", "toprandom"), Lowercase = 1, cutoff.q = 0.05)
{
if (!method %in% c("significant", "top", "toppos", "topnorm", "toprandom")) {
print(paste("Parameter 'method' should be one of: ", paste(c("significant", "top", "topnorm", "toprandom"), collapse=", "), ". The default is 'method = significant'...", sep=""));
}
if (length(method) < 1 ) {
print("Parameter 'method' is undefined. The default 'method = significant' will be used.");
method="significant";
}
if (method =="significant" & !is.na(Ntop)) {stop("Parameter 'Ntop' is irrelevant when 'method = significant'. Terminated...");}
if (is.null(method)) {stop("Parameter 'method' is missing...");}
if (grepl("top", method, ignore.case = T) & is.na(Ntop)) {stop("Parameter 'Ntop' is missing...");}

ags.list <- NULL;
if (is.na(namesFromColumn)) {m1 <- m0;} else {m1 <- m0[,(namesFromColumn+1):ncol(m0)];}
if (!is.na(col.mask)) {m1 <- m1[,colnames(m1)[grep(col.mask,colnames(m1))]];}
uc <- sweep(m1, 1, rowMeans(m1), FUN="-");
if (method=="significant" | method=="topnorm") {
SD <- apply(m1, 1, sd, na.rm=T);
uc <- sweep(uc, 1, SD, FUN="/");
}
if (method=="top" | method=="toppos" | method=="topnorm") {
for (label in colnames(uc)) {
if (method=="toppos") {
x = uc[,label];
} else {
x = abs(uc[,label]);
}
ags.list[[label]] <- tolower(names(x))[order(x, decreasing=T)][1:Ntop];
}
}

if (method=="significant") {
p1 <- 2*pnorm(abs(uc), lower.tail = FALSE);
q1 <- apply(p1, 2, function (x) p.adjust(x, method="BH"));
ags.list <- apply(q1, 2, function (x) tolower(names(x))[which(x < cutoff.q)]);
}
if (method=="toprandom") {
for (label in colnames(uc)) {
x = uc[,label];
ags.list[[label]] <- sample(tolower(names(x)), Ntop);
}
}
return(ags.list);
}
