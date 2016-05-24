metaMA<-function(data, varname, moderated = c("limma", "SMVar", "t")[1], 
    BHth = 0.05, which = c("pval", "ES")[1]){
if (!(which %in%  c("pval", "ES"))) stop("which must be one of: pval, ES")
esets <- GEDM(data)
classes.f <- selectClass(data, varname, "factor")
if (which == "pval") res<-pvalcombination(esets, classes.f, moderated = moderated, BHth = BHth) 
if (which == "ES") res<-EScombination(esets, classes.f, moderated = moderated, BHth = BHth) 
if (!all(sapply(1:(length(GEDM(data))-1), function(x) all(rownames(GEDM(data)[[x]])==rownames(GEDM(data)[[x+1]]))))) stop("The gene expression data matrices have not equal rownames")
res[["gene.names"]]<-rownames(GEDM(data)[[1]])
class(res)<-"metaMA.res"
return(res)
}

ES.GeneMeta<-function(data, varname, useREM = TRUE, CombineExp = 1:length(esets), nperm = 1000)  
{
esets <- GEDM(data)
classes <- selectClass(data, varname, "binary")
theScores <- zScores(esets, classes, useREM = useREM, CombineExp = CombineExp)
ScoresFDR <- zScoreFDR(esets, classes, useREM = useREM, nperm = nperm, CombineExp = CombineExp) 
res<-list(theScores = theScores, ScoresFDR = ScoresFDR)
class(res)<- "ES.GeneMeta.res"
return(res)
}

RankProduct<-function(data, varname,  num.perm = 100, logged = TRUE, na.rm = FALSE, 
    gene.names = NULL, plot = FALSE, rand = NULL, cutoff= 0.05)
{
if (!require(RankProd)) 
        stop("library RankProd is missing")
if (is.null(gene.names)) gene.names=rownames(GEDM(data)[[1]])
 rankdata <- mergedata(data, varname)
 RP.out <- RPadvance(rankdata$dat, rankdata$cl-1, rankdata$origin, 
  num.perm = num.perm, logged = logged, na.rm = na.rm, 
    gene.names = gene.names, plot = plot, rand = rand)
  RankRes <- topGene(RP.out, cutoff = cutoff)  
  class(RankRes)<-"RankProduct.res"
 return(RankRes)
}


VennMapper<-function(data, varname, cutoff)
{
fc <- fold.change(data, varname)
list <- gene.select.FC(fc, cutoff)
ct<-conting.tab(list)
z<-Z(list, length(fc))
gl<-gene.list(list)
res<- list( conting.tab = ct, z.score = z, genes = gl)
class(res)<-"VennMapper.res"
return(res)
}

METRADISC<-function(data, varname, nperm = 1000)
{
metra <- meta.test(data, varname)
RANK <- rank.genes.adv(metra)
RQ <- compute.RQ(RANK)
MC <- MCtest(RANK, RQ, nper = nperm)
res<-list(ranks = RANK, RQ = RQ, MCtest = MC)
class(res)<-"METRADISC.res"
return(res)
}

MAP.Matches<-function(data, varname, t.cutoff = "98.00%", multiple = TRUE, 
	perm=c("both", "columns", "labels")[1], nperm = 1000, test = c("t", "t.equalvar")[1], sig.col, sig.cutoff = 0.05)
{
cat("Examinig the data...\n")
stat.real <- meta.test(data, varname, stat = test)$test
stat <- c(stat.real)
quan <- quantile(abs(stat), seq(0.00, 1.00, 0.0001))
T.default <- quan[t.cutoff]
value.dis <- apply(stat.real, MARGIN = c(1, 2),  function(x) ifelse(abs(x) > T.default, 1, 0))
results <- ratio(value.dis)

if (multiple) {
MAPmat <- MAPmatrix(value.dis)
MAPmat <- MAPmat[MAPmat$n.sig > 1, ]
} else MAPmat <- MAPmatrix(value.dis)

unique.pat <- as.character(MAPmat[, 1])

cat("Statistical analysis...\n")
if (perm == "columns") {
  p1 <- MAPsig1(unique.pat, value.dis, iter = nperm)
  resx <- cbind(MAPmat, p1)
 colnames(resx) <- c(colnames(MAPmat), "p.col.strong", "p.col.weak")
}

if (perm == "labels") {
 out<-test.group.shuffle(data, varname, B = nperm)
 p2 <- MAPsig2(out, value.dis, unique.pat, B = nperm)
 resx <- cbind(MAPmat, p2)
 colnames(resx) <- c(colnames(MAPmat), "p.lab.strong", "p.lab.weak")
}

if (perm == "both") {
p1 <- MAPsig1(unique.pat, value.dis, iter = nperm)
out <- test.group.shuffle(data, varname, B= nperm)
p2 <- MAPsig2(out, value.dis, unique.pat, B = nperm)
resx <- cbind(MAPmat, p1, p2)
colnames(resx) <- c(colnames(MAPmat), "p.col.strong", "p.col.weak", "p.lab.strong", "p.lab.weak")
}

intx <- as.data.frame(t(resx[which(resx[, sig.col] <= sig.cutoff),  ]))
probs <- MAP.genes(resx, value.dis, files = FALSE)
names(probs) <- rownames(resx)

all.genes<-rownames(stat.real)
out<-list(tests = stat.real, bin.matrix = value.dis, sumarization = results, MAP = MAPmat, stat.analysis = intx, genes = probs, all.genes=all.genes)
class(out) <- "MAP.Matches.res"
return(out)
}
