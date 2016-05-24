#' Plot a heatmap of NEA/GSEA output
#'
#' Plots a heatmap where width and height of each element reflect respective sizes of input FGS and AGS.

#' @param List1 AGS or FGS object that lists members of each individual AGS/FGS.
#' @param List2 FGS or AGS object that lists members of each individual FGS/AGS.
#' @param Z matrix (output of \code{\link{nea.render}} with Z-scores that will define coloring of the heatmap.
#' @param Log If \code{TRUE}, then the Z values will be log-transformed (the default).

#' @seealso \code{\link{nea.render}}

#' @examples
#' ags.list <- samples2ags(fantom5.43samples, Ntop=20, method="topnorm")
#' data(can.sig.go)
#' fpath <- can.sig.go
#' fgs.list <- import.gs(fpath)
#' data(net.kegg)
#' netpath <- net.kegg
#' net <- import.net(netpath)
#' n1 <- nea.render(AGS=ags.list, FGS=fgs.list, NET=net, graph=TRUE)
#' set.heat(ags.list, fgs.list, n1$z, Log=FALSE)
#' @keywords internal
#' @importFrom grDevices topo.colors
#' @importFrom graphics hist image
#' @export


set.heat <- function (
List1, List2, # <- gene set lists created e.g. with import.gs()
Z, # <-  similarity/dissimilarity matrix of gene lists
Log = TRUE # <- the heatmap will be on the log scale
) {
Borders <-  NULL;
t1 <- names(table(c(length(List1), length(List2)) %in% dim(Z)))
if (length(t1) == 1 & t1 == "TRUE") {
for (s in c("x", "y")) {
S <- c(0);
if (s == "x") {
if (length(List1) == dim(Z)[2]) {
GS <- List1;
Xlab = "List1";
} else {
GS <- List2;
Xlab = "List2";
}}
if (s == "y") {
if (length(List2) == dim(Z)[1]) {
GS <- List2;
Ylab = "List2";
} else {
GS <- List1;
Ylab = "List1";
}}
for (i in names(GS)) {
S <- c(S, (S[length(S)] + length(GS[[i]])));
}
Borders[[s]]  <- S;
}
if (Log) {
Scores <- log(Z + abs(min(Z, na.rm=T)) + 0.1);
} else {
Scores <- Z;
}
Breaks=hist(sort(Scores,decreasing=T), breaks=100, plot=FALSE)$breaks;
image(Borders$x, Borders$y, t(Scores), breaks=Breaks, col = topo.colors(length(Breaks)-1), xaxs = "i", yaxs = "i", xlab=Xlab, ylab=Ylab)
return(NULL);
} else {
stop("Dimensions of the input elements do not match each other...");
}
}
