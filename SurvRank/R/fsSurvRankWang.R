#'Random Cox ranking function
#'
#' This function ranks the input features according to their mean significance in univariate cox models for a randomly selected subset of the data
#' @param x x matrix or data.frame
#' @param y response y as a survival object, generated with \code{Surv()}
#' @param bb repeatements of the bootstrap
#' @param ns sample size per bootstrap
#' @param ... other arguments, not used now
#' @keywords SurvRank
#' @export
#'
fsSurvRankWang = function(x,y,bb=400,ns=80,...){
  res=list()
  scorei = matrix(NA,ncol=bb,nrow=ncol(x));rownames(scorei) = colnames(x)
  for (j in 1:bb){
    subi = sample(x=1:nrow(y),size=ns,replace=T)
    scorei[,j] = unlist(lapply(1:ncol(x),function(q) 1/summary(survival::coxph(y[subi,]~.,data=data.frame(x[subi,q])))$coefficients[,"Pr(>|z|)"])) #univariate cox models, variables with smaller pvalue get higher scorei
  }
  scorei = t(apply(scorei,1,function(x)sort(x)))
  scorei = scorei[,(ceiling(bb*0.05):ceiling(bb*0.95))]
  m_scor = apply(scorei,1,mean)
  res$rank = names(sort(m_scor,decreasing=T))
  res$nrcoef = length(res$rank)
  return(res)
}
