# code by Stephen Turner (with slight modification "p" -> "Fdr" in ylab)
# http://StephenTurner.us/
# http://GettingGeneticsDone.blogspot.com/
# See license at http://gettinggeneticsdone.blogspot.com/p/copyright.html

manhattan <- function(dataframe, colors=c("gray10", "gray50"), ymax="max", limitchromosomes=1:23, suggestiveline=-log10(1e-5), genomewideline=-log10(5e-8), annotate=NULL, ...) {
  
  CHR <- P <- NULL
  d=dataframe
  if (!("CHR" %in% names(d) & "BP" %in% names(d) & "P" %in% names(d))) stop("Make sure your data frame contains columns CHR, BP, and P")
  
  if (any(limitchromosomes)) d=d[d$CHR %in% limitchromosomes, ]
  d=subset(na.omit(d[order(d$CHR, d$BP), ]), (P>0 & P<=1)) # remove na's, sort, and keep only 0<P<=1
  d$logp = -log10(d$P)
  d$pos=NA
  ticks=NULL
  lastbase=0
  colors <- rep(colors,max(d$CHR))[1:max(d$CHR)]
  if (ymax=="max") ymax<-ceiling(max(d$logp))
  if (ymax<8) ymax<-8
  
  numchroms=length(unique(d$CHR))
  if (numchroms==1) {
    d$pos=d$BP
    ticks=floor(length(d$pos))/2+1
  } else {
    for (i in unique(d$CHR)) {
      if (i==1) {
        d[d$CHR==i, ]$pos=d[d$CHR==i, ]$BP
      } else {
        lastbase=lastbase+tail(subset(d,CHR==i-1)$BP, 1)
        d[d$CHR==i, ]$pos=d[d$CHR==i, ]$BP+lastbase
      }
      ticks=c(ticks, d[d$CHR==i, ]$pos[floor(length(d[d$CHR==i, ]$pos)/2)+1])
    }
  }
  
  if (numchroms==1) {
    with(d, plot(pos, logp, ylim=c(0,ymax), ylab=expression(-log[10](italic(Fdr))), xlab=paste("Chromosome",unique(d$CHR),"position"), ...))
  }	else {
    with(d, plot(pos, logp, ylim=c(0,ymax), ylab=expression(-log[10](italic(Fdr))), xlab="Chromosome", xaxt="n", type="n", ...))
    axis(1, at=ticks, labels=unique(d$CHR), ...)
    icol=1
    for (i in unique(d$CHR)) {
      with(d[d$CHR==i, ],points(pos, logp, col=colors[icol], ...))
      icol=icol+1
    }
  }
  
  if (!is.null(annotate)) {
    d.annotate=d[which(d$SNP %in% annotate), ]
    with(d.annotate, points(pos, logp, col="green3", ...)) 
  }
  
  if (suggestiveline) abline(h=suggestiveline, col="blue")
  if (genomewideline) abline(h=genomewideline, col="red")
}
