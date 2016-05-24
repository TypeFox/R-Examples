#' Weigthing function for selected features
#'
#' Weigthing function for selected features, according to the performance in the outer loop
#' @param out.mat matrix of dimension t.times*folds times features according to training set selections. 0/1 matrix: 1 (0) if the feature was (not) selected in the run
#' @param t.times number of times the cross-validation should be repeated
#' @param cv.out number of folds in outer cross validation loop (for estimation of the predictive accuracy)
#' @param out.sc matrix of t.times (of repeated cross validation) times folds (number of outer CV folds) of prediction accuracy evaluations
#' @param n1 parameter for functions of two different weighting functions. Deviation from average performance leading to double weight of selected coefficients.
#' @keywords SurvRank
#' @export

weighting_fct = function(out.mat,t.times,cv.out,out.sc,n1){
  res=list()
  nv=1/(t.times*cv.out)
  auc_vec = c(unlist(lapply(out.sc,max)))         #maximum auc for each round
  av.rel = (auc_vec-mean(auc_vec,na.rm=T))/mean(auc_vec,na.rm=T)
  wv = nv*exp(log(2)/n1*av.rel)      #weighting factor
  wv[is.na(wv)]=0
  wv[which(auc_vec<0.5)]=0       #set the weighting factor for those auc smaller than 0.5 to
  wv = wv/sum(wv)
  tl1 = data.frame(rel.freq.w=out.mat%*%wv,rel.freq=apply(out.mat,1,sum)/ncol(out.mat))
  res$tl1 = tl1[order(tl1[,1],decreasing=T),]
  res$wv = wv
  return(res)
}
