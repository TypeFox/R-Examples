#' ROC for NEA benchmarks
#' 
#' Plot ROC curve(s) for benchmarked network(s)
#' @param tpvsfp a list contaning data for the ROC curve(s) from \code{\link{benchmark}}: \code{cutoffs}, \code{fp}, \code{tp}, and the point \code{cross.z} that matches \code{coff.fdr}.
#' @param coff.z since this is meaningless to consider statistically insignificant results as "successful", the ROC curve should rather terminate at a certain minimally justifiable level of confidence. This is set by the Z-score cutoff which follows normal distribution parameters, e.g. z=1.96 ~ p-value=0.05 or z=2.57 ~ p-value=0.01 in one-sided tests.
#' @param coff.fdr to make significance levels comparable between different curves, the point where FDR=\code{coff.fdr} will be labeled with a circle (think of TP/FP ratio at this level).
#' @param main title name for the plot
#' @param cex.leg font size in plot

#' @seealso \link{benchmark}
#' @references \url{http://www.biomedcentral.com/1471-2105/15/308}
#'
#' @keywords benchmark

#' @examples
#' # Benchmark and plot one networks on the whole set of test GSs, using no mask:
#' data(can.sig.go);
#' fpath <- can.sig.go
#' gs.list <- import.gs(fpath, Lowercase = 1, col.gene = 2, col.set = 3);
#' data(net.kegg)
#' netpath <- net.kegg
#' net <- import.net(netpath)
#' \donttest{
#' b0 <- benchmark (NET = net, 
#'  GS = gs.list[c("kegg_04270_vascular_smooth_muscle_contraction")],
#'  echo=1, graph=TRUE, na.replace = 0, mask = ".", minN = 0, 
#'  coff.z = 1.965, coff.fdr = 0.1, Parallelize=2);
#'  roc(b0, coff.z = 1.64, coff.fdr = 0.05);
#' }
#' \dontrun{
#' ## Benchmark and plot a number of networks on GO terms and KEGG pathways separately, using masks
#' b1 <- NULL;
#' for (mask in c("kegg_", "go_")) {
#' b1[[mask]] <- NULL;
#' for (file.net in c("netpath")) { 
#' # a series of networks can be put here: c("netpath1", "netpath2", "netpath3")
#' net <- import.net(netpath, col.1 = 1, col.2 = 2, Lowercase = 1, echo = 1)
#' b1[[mask]][[file.net]] <- benchmark (NET = net, GS = gs.list, 
#' gs.gene.col = 2, gs.group.col = 3, net.gene1.col = 1, net.gene2.col = 2, 
#' echo=1, graph=FALSE, na.replace = 0, mask = mask, minN = 0, Parallelize=2);
#' }}
#' par(mfrow=c(2,1));
#' roc(b1[["kegg_"]], coff.z = 2.57, coff.fdr = 0.01);
#' roc(b1[["go_"]], coff.z = 2.57, coff.fdr = 0.01);
#' }
#'
#' @importFrom grDevices rainbow
#' @importFrom graphics abline legend points
#' @export


roc <- function(tpvsfp, coff.z = 1.965, coff.fdr = 0.1, cex.leg = 0.75, main=NA) {
if (mode(tpvsfp) != "list") {
print('Submitted first agrument is wrong. Should be either a single list with entries c("cross.z", "cutoffs", "fp", "tp", "ne", "nv") or a list of such lists...');
return();
}

if (suppressWarnings(all(sort(names(tpvsfp)) == c("cross.z", "cutoffs", "fp", "ne", "nv", "tp")))) {
objs <- NULL; objs[["1"]]=tpvsfp;
} else {
objs=tpvsfp;
} 
Nets <- sort(names(objs));

if (length(Nets) > 1) {
Col <- rainbow(length(Nets)); names(Col) <- Nets;
} else {
Col <- NULL; Col["1"] <- c("black");
}
Max.tp <- 0;
for (p1 in Nets) {
pred <- objs[[p1]];
Max.signif = max(pred$tp[[1]][which(pred$cutoffs[[1]] > coff.z)], na.rm=T);
if (Max.tp < Max.signif ) {
Max.tp <- Max.signif;
}
}

Xlim = c(0,  Max.tp/1); Ylim = c(0, Max.tp);
plot(1, 1, xlim=Xlim, ylim=Ylim,  
cex.main = 1.0, type="n", xlab="False predictions", ylab="True predictions", main=main);
if (length(Nets) > 1) {
Le = paste(substr(Nets, 1, 15), "; n(V)=", sapply(objs, function (x) x$nv), "; n(E)=", sapply(objs, function (x) x$ne), sep="");
legend(x="bottomright", legend=Le, title=paste( "", sep=""), col=Col, bty="n", pch=15, cex=cex.leg);
}

for (p1 in Nets) {
pred <- objs[[p1]];
X = pred$fp[[1]][which(pred$cutoffs[[1]] > coff.z)];
Y = pred$tp[[1]][which(pred$cutoffs[[1]] > coff.z)];
points(X,Y, pch=".", type="l", col=Col[p1]);
abline(0,1,col = "gray60",lwd=0.5, lty=2);
Pos <- pred$cutoffs[[1]][which(pred$cutoffs[[1]] > pred$cross.z)];
points(
pred$fp[[1]][which((Pos - pred$cross.z) == min(Pos - pred$cross.z, na.rm=T))],
pred$tp[[1]][which((Pos - pred$cross.z) == min(Pos - pred$cross.z, na.rm=T))], 
pch="O", type="p", col=Col[p1]);
}
}
