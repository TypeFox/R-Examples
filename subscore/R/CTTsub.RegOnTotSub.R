#' This function Computes subscores based on both the observed subscore and the observed total in classical test theory (CTT)
#' score
#' @description This function computes CTT subscores based on both the observed subscore and the observed total
#' score, using the method introduced in Haberman (2008), Puhan, Sinharay, Haberman, and Larkin(2008), and Sinharay (2010), which return:\cr
#' 	(1) Original observed subscore; \cr
#' 	(2) The true subscore is estimated based on 
#' 	    both the observed subscore and the observed  total score.\cr
#' @param test.data A list that contains datasets of all subtests and the whole test,
#'                  which can be obtained using function 'data.prep'.
#' @return A list of objects that include both test information and subscores. \cr
#'         (1) "subscore.information" - It contains test information of both subtests
#'                                      and the total test, such as mean of subscores and total score, reliability, and PRMSE.\cr
#'         (2) "subscore.original" - It contains original subscores and total score. \cr
#'         (3) "subscore.RegOnTotSub" - It contains subscores that are estimated based on both the observed subscore and the
#'                                      observed total score. 
#' @import CTT
#' @import stats
#' @examples 
#'         CTTsub.RegOnSub(test.data)
#' @export

CTTsub.RegOnTotSub<-function (test.data) {
  
  n.tests<-length(test.data)
  n.subtests<-n.tests-1
  n.items<-rep(NA,n.tests)
  n.cases<-rep(NA,n.tests)

  for (t in 1:n.tests) {
    n.items[t]<-dim(test.data[[t]])[2] 
    n.cases[t]<-dim(test.data[[t]])[1]  
    } 
  
  n.items.total<-n.items[n.tests]
  reliability.alpha<-rep(NA, (n.tests))  
  mylist.names <- c(paste ('Subscore.',rep(1:n.subtests),sep=''),'Score.Total')
  subscore.list <- as.list(rep(NA, length(mylist.names)))
  names(subscore.list) <- mylist.names
  
  for (t in 1 : (n.tests))  {
    subscore.list[[t]]<- rowSums(test.data[[t]],na.rm = TRUE)
  } 
  subscore.original.matrix<-do.call(cbind, subscore.list) 
  
  for (r in 1:(n.tests)) {
    reliability.alpha[r]<-reliability(test.data[[r]],itemal=TRUE,NA.Delete=T)[[3]]
  } 
 
  sigma.obs<-rep(NA,n.tests)
  
  for (t in 1:n.tests) {
    sigma.obs[t]<-sd(subscore.list[[t]],na.rm = TRUE)
  }
  var.obs<-sigma.obs^2
  corr<-cor(subscore.original.matrix)
  CovMat.Obs<-cov(subscore.original.matrix)
  var.true<-var.obs*reliability.alpha
  sigma.true<-sqrt(var.true)
  CovMat.true<-CovMat.Obs
  for (t in 1:n.tests) {
    CovMat.true[t,t]<-var.true[t]
  }
  
  mean<-rep(NA,n.tests)
  for (t in 1:n.tests) {
    mean[t]<-mean(subscore.list[[t]],na.rm = TRUE)
  }
  subscore.dataframe<-as.data.frame(subscore.original.matrix)
  
  cov.rowsum<-rowSums(CovMat.true[,1:n.subtests],na.rm = TRUE)
  
  r.StXt<-rep(NA,n.tests)
  for (t in 1:(n.tests-1)) {     
    r.StXt[t]<-cov.rowsum[t]^2/(var.true[t]*var.true[n.tests])
  } 
  
  tao<-rep(NA,n.tests)
  beta<-rep(NA,n.tests)
  gamma<-rep(NA,n.tests)

  for (t in 1:n.tests) {
    tao[t]<-(sqrt(reliability.alpha[n.tests])*sqrt(r.StXt[t])-corr[t,n.tests]*sqrt(reliability.alpha[t]))/(1-corr[t,n.tests]^2)
    beta[t]<- sqrt(reliability.alpha[t])*(sqrt(reliability.alpha[t])-corr[t,n.tests]*tao[t])
    gamma[t]<-sqrt(reliability.alpha[t])*tao[t]*(sigma.obs[t]/sigma.obs[n.tests])
  } 
  
  PRMSE.RegOnTotSub<-rep(NA, n.tests)
  
  for (t in 1:n.subtests) { 
    PRMSE.RegOnTotSub[t]<-reliability.alpha[t]+tao[t]^2*(1-corr[t,n.tests]^2)
  } 

    mylist.names <- c(paste ('RegOnTotSub.Score.',rep(1:n.subtests),sep=''))
  subscore.list.RegOnTotSub <- as.list(rep(NA, length(mylist.names)))
  names(subscore.list.RegOnTotSub) <- mylist.names
  
  for (t in 1: n.subtests) { 
    subscore.list.RegOnTotSub[[t]]<-mean[t]+beta[t]*(subscore.dataframe[,t]-mean[t])+gamma[t]*(subscore.dataframe[,n.tests]-mean[n.tests])
  } 
  
  subscore.information.list<-list(mean=mean, reliability=reliability.alpha, 
                                  PRMSE.RegOnTotSub=PRMSE.RegOnTotSub)
  subscore.information<-do.call(cbind,subscore.information.list)
  rownames.list<-c(paste('Subtest.',rep(1:n.subtests),sep=''),'Total.test')
  rownames(subscore.information)<-rownames.list
  subscore.original<-do.call(cbind,subscore.list)
  subscore.RegOnTotSub<-do.call(cbind,subscore.list.RegOnTotSub)
  subscores<-list(subscore.original=subscore.original, subscore.RegOnTotSub=subscore.RegOnTotSub)
  
  return (list(subscore.information=subscore.information, 
               subscore.original=subscore.original,
               subscore.RegOnTotSub=subscore.RegOnTotSub))
} 