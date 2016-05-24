selfcontained.test<-function(pvalue,weight=NA,p_permu=NA) 
{ 
  handle<-function(x)
  {
    x[x>=0.999999]<-0.999999
    x[x< 0.000001]<-0.000001
    return( x )
  }
  pvalue=as.numeric(pvalue)
  pvalue<-handle(pvalue)  
  if (sum(is.na(pvalue))>0) {stop ("Error: pvalue is NA.")}
  if (min(pvalue)<0) {stop ("Error: pvalue is negative.")}
  if (max(pvalue)>1) {stop ("Error: pvalue is >1.")}
  if (is.na(weight)[1]) {weight=rep(2, length(pvalue))}
  if (length(pvalue)!=length(weight)){stop ("Error: length(pvalue)!=length(weight).")}
  if (min(weight)<0) {stop ("Error: weight is negative.")}
  if (!is.na(p_permu)[1]) {if (nrow(p_permu)!=length(pvalue)) {stop ("Error: nrow(p_permu)!=length(pvalue)")}}
  if (!is.na(p_permu)[1]) {p_permu<-handle(p_permu)}
  
  et<-sum(weight)
  fun2<-function(x){qchisq(1-x,df=weight)}
  if (is.na(p_permu)[1]) {vart<-2*sum(weight)} else {vart<-sum(cov(t(apply(p_permu,2,fun2))))}
  v<-2*et^2/vart
  c<-vart/2/et
  T<-sum(qchisq(1-pvalue,weight)) 
  if (is.na(p_permu)[1]) {cor_p<-diag(length(pvalue))} else {cor_p<-cor(t(p_permu))}
  return(list("significance level for combining pvalues"=pchisq(T/c, v, lower.tail=FALSE), "correlation among pvalues"=cor_p))
}

competitive.test<-function(Pvalue,Weight,p_random=NA)
{
  if (is.na(p_random)[1]) {p_random<-matrix(runif(length(Pvalue)*100000),nrow=length(Pvalue))}
  if (nrow(p_random)!=length(Pvalue)) {stop ("Error: nrow(p_random)!=length(Pvalue)")}
  
  random<-NULL
  for (i in 1:ncol(p_random))
  {
    random[i]<-selfcontained.test(pvalue=p_random[,i],weight=Weight,p_permu=NA)[[1]]
  }
  return(list("significance level for combining pvalues"=mean(selfcontained.test(Pvalue,Weight,p_permu=NA)[[1]]>random))) 
}