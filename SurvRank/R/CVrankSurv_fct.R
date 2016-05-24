#' Main function of SurvRank.
#'
#' Main input function for SurvRank.
#' @param data input of data as a list in the format: list.name$x data.frame of covariates. list.name$y response as a survival object, derived from \code{\link{Surv}}.
#' @param t.times number of times the cross-validation should be repeated
#' @param cv.out number of folds in outer cross validation loop (for estimation of the predictive accuracy)
#' @param cv.in number of folds in inner cross validation loop (for model selection on the training set)
#' @param fs.method Defaults to "lasso.rank". Ranking method to be applied. One of c("lasso.rank","conc.rank","rf.rank","boost.rank","cox.rank","rpart.rank","randcox.rank","wang.rank")
#' @param nr.var Number of variables up to which stepwise selection should be carried out. Has to be smaller than n number of observations.
#' @param sd1 factor to which sparser solutions should be chosen. Not maximum Survival AUC in inner loop is used in stepwise selection, instead \code{max(survAUC)*sd1} leading to sparser solutions
#' @param ncl Defaults to 1. Number of clusters for parallel execution.
#' @param weig.t Defaults to TRUE. Should a weighting of features be performed.
#' @param n1 used in weighting function if weig.t=T. Find details in \code{\link{weighting_fct}}
#' @param c.time as defined in package \code{survAUC} time; a positive number restricting the upper limit of the time range under consideration.
#' @param ... arguments that can be passed to underlying functions, not used now
#' @keywords SurvRank
#' @export
#' @details details to follow
#' @return Output of the \code{CVrankSurv_fct}, basically a list containing the following elements
#' \item{\code{method}}{ranking method}
#' \item{\code{accuracy$ranking}}{full ranking of all model estimations}
#' \item{\code{accuracy$pred.in}}{averaged inner AUCs of stepwise selection}
#' \item{\code{accuracy$pred.out}}{predictions of testset}
#' \item{\code{accuracy$used.rank}}{only used features according to stepwise selection}
#' \item{\code{accuracy$used.rank1se}}{only used features according to stepwise selection with factor \code{sd1}}
#' \item{\code{accuracy$auc.out}}{matrix of dimension \code{cv.out} times \code{t.times} of survival AUCs.}
#' \item{\code{accuracy$auc.out1se}}{matrix of dimension \code{cv.out} times \code{t.times} of survival AUCs with factor \code{sd1}.}
#' \item{\code{rank.mat}}{matrix of ranks per feature. If not selected, it is set to number of features.}
#' \item{\code{out.mat}}{0/1 matrix for features selected}
#' \item{\code{out.mat1se}}{0/1 matrix for features selected with factor \code{sd1} application}
#' \item{\code{top1se}}{unweighted toplist with factor \code{sd1}}
#' \item{\code{toplist}}{unweighted toplist}
#' \item{\code{weighted}}{weighted toplist with applied weighting function}
#' \item{\code{rank}}{toplist of ranked features according to ranks}
#' @examples
#' ## Simulating a survival data set
#' N=100; p=10; n=4
#' x=data.frame(matrix(rnorm(N*p),nrow=N,p))
#' beta=rnorm(n)
#' mx=matrix(rnorm(N*n),N,n)
#' fx=mx[,seq(n)]%*%beta/3
#' hx=exp(fx)
#' ty=rexp(N,hx)
#' tcens=1-rbinom(n=N,prob=.3,size=1)
#' y=Surv(ty,tcens)
#' data=list()
#' data$x<-x; data$y<-y
#' ## Ranking the features according to their significance in the univariate cox models
#' out.cox<-CVrankSurv_fct(data,2,3,3,fs.method="cox.rank")
#' ## Ranking the features according to the LASSO algorithm
#' \dontrun{
#' out.lasso<-CVrankSurv_fct(data,2,5,5,fs.method="lasso.rank")}

CVrankSurv_fct = function(data, t.times, cv.out, cv.in, fs.method="lasso.rank",
                            nr.var=10, sd1=0.95,ncl=1,weig.t = T,n1=0.1,c.time=10,...){
  pb = utils::txtProgressBar(style=3,min=1,max=t.times*cv.out)
  res = ranking = pred.in = pred.out = used.rank = used.rank1se = list()
  auc.out = auc.out1 = NULL
  stopifnot(nrow(data$x)==nrow(data$y),
            is.list(data),
            length(data)==2,
            survival::is.Surv(data$y))
  if(is.null(colnames(data$x))){
    colnames(data$x) = paste("X",1:ncol(data$x),sep="")
    warning("No column names in data$x, set to X1-Xncol")
  }
  # prediction accuracy
  for (t in 1:t.times){
    if(cv.out!=length(data$y[,2])){             # create stratified cross validation folds
      fold.outer = crossvalFolds(data$y[,2],k=cv.out)  # using the function crossvalFolds
    }
    utils::setTxtProgressBar(pb,t*cv.out-1)
    output = msSurv_fct(data=data,fold=fold.outer,fs.method=fs.method,
                         ncl=ncl,cv.out=cv.out,cv.in=cv.in,nr.var=nr.var,t=t,sd1=sd1,c.time=c.time,...)  #go to the ranking and inner_AUC functions, creates the output
    for (i in 1:cv.out){
      l=(t-1)*cv.out+i
      ranking[[l]] = output[[i]][[1]];pred.in[[l]] = output[[i]][[2]];pred.out[[l]]=output[[i]][[3]];used.rank[[l]]=output[[i]][[4]];auc.out[l]=output[[i]][[5]];used.rank1se[[l]]=output[[i]][[6]];auc.out1[l]=output[[i]][[7]]
    }
  }       #creating an output list
  auc.out[which(auc.out==0)]=NA    #those models that have an AUC of 0, set it to NA
  auc.o = matrix(auc.out,nrow=cv.out,ncol=t.times)  #create an cv.out*t.times AUC matrix
  res$method = c(fs.method)        #the ranking method used
  res$accuracy = list(ranking = ranking, pred.in = pred.in, pred.out=pred.out,used.rank=used.rank, auc.out = auc.o) #create list with the outcomes for each of the rounds
  out.mat=matrix(0,nrow=ncol(data$x),ncol=t.times*cv.out,dimnames = list(colnames(data$x))) #matrix of dimensions number of variables times number of rounds
  rank.mat=matrix(ncol(data$x),nrow=ncol(data$x),ncol=t.times*cv.out,dimnames = list(colnames(data$x)))
  for (i in 1:length(res$accuracy$used.rank)){
    out.mat[match(res$accuracy$used.rank[[i]],rownames(out.mat)),i]=1    #for each of the variables, set 0 or 1 on whether they were used in the model for the considered round
    rank.mat[match(res$accuracy$used.rank[[i]],rownames(out.mat)),i] = 1:length(res$accuracy$used.rank[[i]])   #gives the order in which the variables were used
  }
  res$rank.mat=rank.mat

 res$out.mat=out.mat
  auc.o1 = matrix(auc.out1,nrow=cv.out,ncol=t.times)   #create an cv.out*t.times AUC matrix for sparsed solutions
  res$accuracy = list(ranking = ranking, pred.in = pred.in, pred.out=pred.out,used.rank=used.rank,used.rank1se=used.rank1se, auc.out = auc.o,auc.out1se = auc.o1) #update the list with the outcomes
 out.mat1=matrix(0,nrow=ncol(data$x),ncol=t.times*cv.out,dimnames = list(colnames(data$x)))
  for (i in 1:length(res$accuracy$used.rank1se)){                           #create a 0/1 matrix for the variables used in the models
    out.mat1[match(res$accuracy$used.rank1se[[i]],rownames(out.mat1)),i]=1
  }
  res$out.mat1se = out.mat1
  tl1 = data.frame(table(unlist(res$accuracy$used.rank1se)))
  tl1 = data.frame(name=tl1[,1],abs.freq=tl1[,2],rel.freq=tl1[,2]/(cv.out*t.times))
  tl1 = tl1[order(tl1[,2],decreasing=T),]                    #a list with the absolute frequency for each of the variables in the sparser models (in how many models it was included)
  res$top1se = tl1
  tl = data.frame(table(unlist(res$accuracy$used.rank)))
  tl = data.frame(name=tl[,1],abs.freq=tl[,2],rel.freq=tl[,2]/(cv.out*t.times))   #the same for the non sparse solutions
  tl = tl[order(tl[,2],decreasing=T),]
  res$toplist = tl
  if(weig.t==T){          #if we use a weighting function
    ww_f = weighting_fct(out.mat=out.mat,t.times=t.times,cv.out=cv.out,out.sc=res$accuracy$auc.out,n1=n1)
    res$weighted = ww_f$tl1
    rr = rank.mat%*%ww_f$wv
    rr1 = rr[order(rr),]
    rr2 = sum(apply(rank.mat,2,function(x) mean(unique(x)))*ww_f$wv)
    rr3=rep(0,length(rr1))
    rr3[which(rr1<rr2)]=1
    res$rank = data.frame(rank=rr1,selected=rr3)
  }
  return(res)
}
