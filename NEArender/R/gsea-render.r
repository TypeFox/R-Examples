#' Gene Set Enrichment Analysis (GSEA)
#' 
#' A binomial version of GSEA, unified as much as possible with \code{\link{nea.render}}. Given the altered gene sets (AGS) and functional gene sets (FGS), calculates no. of members (genes/protein IDs) shared by each AGS-FGS pair as well as respective enrichment statistics. Returns matrices of size \code{length(FGS)} x \code{length(AGS)} (see "Value"). Each of these two parameters can be submitted as either a text file or as an R list which have been preloaded with \code{\link{import.gs}}.

#' @param AGS Either a text file or a list of members of each AGS (see Details). Group IDs should be found in \code{ags.group.col} and gene IDs would be found in \code{ags.gene.col}. Identical to AGS needed for in \code{\link{nea.render}} - see also details there. 
#' @param FGS Either a text file or a list of members of each FGS (see Details). Group IDs should be found in \code{fgs.group.col} and gene IDs would be found in \code{fgs.gene.col}. Alsmost identical to FGS needed for in \code{\link{nea.render}}. 
#' @param Lowercase render node and group IDs lower-case (Default:1, i.e. 'yes').
#' @param ags.gene.col number of the column containing AGS genes (only needed if AGS is submitted as a text file).
#' @param ags.group.col number of the column containing group IDs (only needed if AGS is submitted as a text file).
#' @param fgs.gene.col number of the column containing FGS genes (only needed if FGS is submitted as a text file).
#' @param fgs.group.col number of the column containing group IDs (only needed if FGS is submitted as a text file).
#' @param Ntotal The important parameter for precise calculation of the Fisher's statistics: how big is the whole gene universe? Defaults to 20000 but should be changed depending on the hypothesis and genome/proteome size. 
#' @param echo if messages about execution progress should appear.
#' @param Parallelize The number of CPU cores to be used for calculating the gene set overlap. The other steps are sufficiently fast. The option is not supported in Windows. 
#' @return A list of entries \code{estimate}, \code{p}, \code{q}, and \code{n}, each of which is a matrix of size \code{length(FGS)} x \code{length(AGS)}. The two former ones contain respective output of \code{\link{fisher.test}}: \code{p.value} and \code{estimate}, whereas \code{q} is produced by \code{p.adjust(p.value, method="BH")} and \code{n} is the no. of shared  members. Input to \code{\link{fisher.test}} is \code{matrix(c(<no. of shared  members>, <no. of solely FGS  members>, <no. of solely AGS  members>, <no. of non-members>), nrow=2)}.
#' @seealso \code{\link{nea.render}}, \code{\link{import.gs}}
#' @keywords AGS
#' @keywords FGS
#' @examples
#' ags.list <- samples2ags(fantom5.43samples, Ntop=20, method="topnorm")
#' data(can.sig.go)
#' fpath <- can.sig.go
#' fgs.list <- import.gs(fpath)
#' g1 <- gsea.render(AGS=ags.list, FGS=fgs.list, Lowercase = 1, 
#' ags.gene.col = 2, ags.group.col = 3, fgs.gene.col = 2, fgs.group.col = 3, 
#' echo=1, Ntotal = 20000, Parallelize=1)
#' hist(g1$estimate, breaks=100)
#' hist(g1$n, breaks=100)
#' hist(g1$p, breaks=100)
#' hist(g1$q, breaks=100)
#'@export

gsea.render <- function (AGS, FGS , Lowercase = 1, ags.gene.col = 2, ags.group.col = 3, fgs.gene.col = 2, fgs.group.col = 3, echo=1, Ntotal = 20000, Parallelize=1) {
if (echo>0) {print("Preparing input datasets:");}

if (is.list(FGS)) {
fgs.list <- FGS;
} else {
fgs.list <- import.gs(FGS, Lowercase=Lowercase, fgs.gene.col, fgs.group.col, gs.type = 'f') 
if (echo>0) {print(paste("FGS: ", length(unique(unlist(fgs.list))), " genes in ", length(fgs.list), " groups.", sep = ""));} 
}

if (is.list(AGS)) {ags.list <- AGS
} else {
ags.list <- import.gs(AGS, Lowercase=Lowercase, ags.gene.col, ags.group.col, gs.type = 'a')
}

#we make here an important assumption: the full set size of the cluster members ('the gene universe') is the number of distinct members of the  union {AGS, FGS}. Otherwise, one can submit a specified value of Ntotal.
if (is.na(Ntotal) | is.null(Ntotal)) {
Ntotal = length(unique(c(unlist(AGS), unlist(FGS))));
} 
if (echo>0) {
print(paste("AGS: ",  length(unique(unlist(ags.list))), " genes in ", length(ags.list), " groups.", sep = ""));		
}
if (echo>0) {print("Calculating overlap statistics...");}
Nov.par <- function(i, ags) {
return(
lapply(ags, 
function (x) {
y = fgs.list[[i]];
Nol <- length(which(unlist(x) %in% unlist(y)));
f1 <- fisher.test(matrix(c(Nol, length(x) - Nol, length(y) - Nol, Ntotal - length(x) - length(y) + Nol), nrow = 2));
return(c(f1$p, Nol, f1$estimate));
}
));
}
if (Parallelize > 1) {
l1 <- mclapply(names(fgs.list), Nov.par, ags.list, mc.cores = Parallelize);
} else {
l1 <- lapply(names(fgs.list), Nov.par, ags.list);
}
GSEA.ags_fgs <- array(unlist(l1), dim=c(3, length(ags.list), length(fgs.list)), dimnames=list(Stats=c("P", "N", "OR"), ags.names=names(ags.list), fgs.names=names(fgs.list)));
stats <- NULL;
stats$n <-  t(GSEA.ags_fgs["N",,]);
stats$estimate <-  t(GSEA.ags_fgs["OR",,]);
va = 'estimate'; 
# Min = min(stats[[va]][which(!is.infinite(stats[[va]]))], na.rm=T) - 0;
Max = max(stats[[va]][which(!is.infinite(stats[[va]]))], na.rm=T) + 1;
# stats[[va]][which(is.infinite(stats[[va]]) & (stats[[va]] < 0))] = Min;
stats[[va]][which(is.infinite(stats[[va]]) & (stats[[va]] > 0))] = Max;

stats$p <- t(GSEA.ags_fgs["P",,]);
stats$q <- matrix(
p.adjust(stats$p, method="BH"), 
nrow = nrow(stats$p), ncol = ncol(stats$p), 
dimnames = list(rownames(stats$p), colnames(stats$p)), 
byrow=FALSE);
stats$q[which(stats$estimate < 1)] = 1
# stats$q <- t(p.adjust(stats$p, method="BH"));
return(stats);
}
