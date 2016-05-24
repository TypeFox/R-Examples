#'Random Cox ranking function
#'
#' This function ranks the input features according to their median significance in cox models with randomly selected nmax variables
#' @param x x matrix or data.frame
#' @param y response y as a survival object, generated with \code{Surv()}
#' @param bb set to 500. Number of random subsamples
#' @param ... other arguments, not used now
#' @keywords SurvRank
#' @export
fsSurvRankRandCox = function(x,y,bb=500,...){
  res=list()
  ranki= fsSurvRankGlmnet(x=x,y=y)
  r100 = ranki$rank[1:ranki$nrcoef]
  nmax = ceiling(nrow(y)/3)
  scorei = matrix(NA,ncol=bb,nrow=ranki$nrcoef);rownames(scorei) = r100
  #scorei = matrix(NA,ncol=bb,nrow=ncol(x));rownames(scorei) = colnames(x)
  for(m in 1:bb){
    subi = sample(x = r100,size = nmax,replace = F)
    #subi = sample(x = colnames(x),size = nmax,replace = F)
    coxi = tryCatch(summary(survival::coxph(y~.,data=data.frame(x[,subi])))$coefficients[,"z"], error = function(e) NA)
    scorei[match(names(coxi),rownames(scorei)),m] = abs(coxi)
  }
  m_scor = apply(scorei,1,median,na.rm=T)
  res$rank = names(sort(m_scor,decreasing=T))
  res$nrcoef = length(res$rank)
  return(res)
}
