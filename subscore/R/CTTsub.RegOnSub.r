#' This function computes subscores based on the observed subscore in classical test theory (CTT).
#' @description This function computes CTT subscores based on the observed subscore, using the methods introduced in 
#' Haberman (2008), Puhan, Sinharay, Haberman, and Larkin(2008), and Sinharay (2010), which returns:\cr
#' 	(1) Original observed subscore; \cr
#' 	(2) The true subscore is estimated based on the observed subscore;\cr
#' @param test.data A list that contains datasets of all subtests and the whole test,
#'                  which can be obtained using function 'data.prep'.
#' @return A list of objects that include both test information and subscores: \cr
#'         (1) "subscore.information" - It contains test information of both subtests
#'                                      and the total test, such as mean of subscores and total score, reliability, and PRMSE.\cr
#'         (2) "subscore.original" - It contains original subscores and total score.\cr
#'         (3) "subscore.RegOnSub" - It contains subscores that are estimated based on the observed subscore. 
#' @import CTT
#' @import stats
#' @examples 
#'         CTTsub.RegOnSub(test.data)
#' @export


CTTsub.RegOnSub<-function (test.data) {
  
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
  
  mylist.names <- c(paste ('RegOnSub.Score.',rep(1:n.subtests),sep=''))
  subscore.list.RegOnSub <- as.list(rep(NA, length(mylist.names)))
  names(subscore.list.RegOnSub) <- mylist.names
  
  subscore.dataframe<-as.data.frame(subscore.original.matrix)
  
  for (t in 1: n.subtests) {
    subscore.list.RegOnSub[[t]]<-mean[t]+reliability.alpha[t]*(subscore.dataframe[,t]-mean[t])
  } 
 
  PRMSE.RegOnSub<-rep(NA,n.tests)
  PRMSE.RegOnSub[1:n.subtests]<-reliability.alpha[1:n.subtests]
  
  subscore.information.list<-list(mean=mean, reliability=reliability.alpha, 
                                  PRMSE.RegOnSub=PRMSE.RegOnSub)
  subscore.information<-do.call(cbind,subscore.information.list)
  rownames.list<-c(paste('Subtest.',rep(1:n.subtests),sep=''),'Total.test')
  rownames(subscore.information)<-rownames.list
  
  subscore.original<-do.call(cbind,subscore.list)
  subscore.RegOnSub<-do.call(cbind,subscore.list.RegOnSub)
  
  subscores<-list(subscore.original=subscore.original, subscore.RegOnSub=subscore.RegOnSub)
  
  return (list(subscore.information=subscore.information, 
               subscore.original=subscore.original,
               subscore.RegOnSub=subscore.RegOnSub))
} 