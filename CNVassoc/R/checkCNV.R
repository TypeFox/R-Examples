checkCNV <-
function(obj, alpha=0.05) {
  out <- TRUE
  meanRatio <- attr(obj,"meanRatio")
  if (is.null(meanRatio))
    stop("CNV must contain intensity data")
  else{
    mins <- by(meanRatio,obj,min)
    maxs <- by(meanRatio,obj,max)
    out <- cor(mins,maxs,method="spearman")==1
  }
  pvalue <- attr(obj,"mixture")$P
  if (!is.null(pvalue))
    out <- !is.nan(pvalue) && pvalue >= alpha && out
  return(out)
}