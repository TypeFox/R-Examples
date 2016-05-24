##' For one constituent study and one simulation of its outcome, test if p<type-1 error rate (not for end user).
##' @param r_events_control Number of events in untreated
##' @param n_sample_size_control Sample size in untreated group
##' @param n_sample_size_treated Sample size in treated group
##' @param alpha Type-1 error rate.
##' @param OR_hat Summary odds ratio from meta-analysis
##' @return Number, 1 if positive, 0 if negative.
test_one <- function(r_events_control, n_sample_size_control, n_sample_size_treated, OR_hat, alpha) {
  RR_hat<-OR_hat/(1-((sum(r_events_control))/(sum(n_sample_size_control)))+OR_hat*((sum(r_events_control))/(sum(n_sample_size_control))))
  pi_control <- r_events_control/n_sample_size_control
  pi_exposed <-1/ (1+((n_sample_size_control-r_events_control)/(RR_hat*r_events_control)))
  sim_control_event <- rbinom(1,n_sample_size_control,pi_control)
  sim_exposed_event <-rbinom(1, n_sample_size_treated, pi_exposed)
  sim_control_no_events <-n_sample_size_control-sim_control_event
  sim_exposed_no_events <- n_sample_size_treated-sim_exposed_event
  matrix1 <-matrix(c(sim_exposed_event,sim_exposed_no_events,sim_control_event, sim_control_no_events), nrow=2,ncol=2) 
  test<-fisher.test(matrix1)
  positive<-ifelse(test$p.value<alpha & test$estimate<1,1,0)
  return(as.numeric(positive))
}
NULL

##' This runs test_one function for many replicates, usually 10,000, for example, to generate the type-2 error rate for one constituent study (not for end user).
##' @param r_events_control Number of events in untreated
##' @param n_sample_size_control Sample size in untreated group
##' @param n_sample_size_treated Sample size in treated group
##' @param alpha Type-1 error rate.
##' @param n Number of iterations used to generate constituent study power; suggest use 10,000.
##' @param OR_hat Summary odds ratio from meta-analysis
##' @return Type-2 error for one constituent study
test.n <- function(r_events_control, n_sample_size_control, n_sample_size_treated, OR_hat, n, alpha) {
  x <- as.vector(1:n)
  for (i in 1:n) x[i] <- test_one(r_events_control, n_sample_size_control, n_sample_size_treated, OR_hat,alpha)
  y<- sum(x)/length(x)
  return(y)
}

NULL
##' calculate the number of expected events from a study included in a meta-analysis assuming summary effect estimate (OR_hat) is true (not for end user).
##'  @param vec_r_events_control an ordered vector of number of events in the untreated group of constituent studies from a meta-analysis.
##'  @param vec_n_sample_size_control an ordered vector of the number of participants in the untreated group.
##'  @param vec_n_sample_size_treated an ordered vector of the number of participants in the treated group.
##'  @param alpha Type-1 error rate.
##'  @param n Number of iterations used to generate constituent study power; suggest use 10,000.
##'  @param OR_hat Summary odds ratio from meta-analysis
##'  @return Expected number of events at given alpha level.
expected_events <- function(vec_r_events_control, vec_n_sample_size_control, vec_n_sample_size_treated, OR_hat, n, alpha){
  power_ind<-as.vector(1:length(vec_r_events_control))
  for (i in 1:length(power_ind)) {
    power_ind[i] <- test.n(vec_r_events_control[i], vec_n_sample_size_control[i], vec_n_sample_size_treated[i], OR_hat, n, alpha) 
  }
  expected<-sum(power_ind)
  return(expected)
}
NULL

##' Chi-square test to test for significant difference between observed and expected number of positive studies (not for end user).
##'  @param vec_r_events_control an ordered vector of number of events in the untreated group of constituent studies from a meta-analysis.
##'  @param vec_n_sample_size_control an ordered vector of the number of participants in the untreated group.
##'  @param vec_n_sample_size_treated an ordered vector of the number of participants in the treated group.
##'  @param alpha Type-1 error rate.
##'  @param OR_hat Summary odds ratio from meta-analysis
##'  @param n Number of iterations used to generate constituent study power; suggest use 10,000.
##'  @param vec_pos Vector of positive results from constituent studies, returned by test.n.treated funtion.
##'  @return Vector of p-values for difference between observed and expected number of positive studies from meta-analysis, along with vector of expected values.
ChisqTest_expect <- function(vec_r_events_control, vec_n_sample_size_control, vec_n_sample_size_treated, OR_hat, n, alpha, vec_pos)  {
  obs<-sum(as.vector(vec_pos))
  expect <- expected_events(vec_r_events_control, vec_n_sample_size_control, vec_n_sample_size_treated, OR_hat, n, alpha)
  teststat<- ((obs-expect)^2)/expect + (obs-expect)^2/(length(vec_pos)-expect)
  a<-1-pchisq(teststat, df=1)
  d<- data.frame(cbind(a,expect))
  return(d)
}
NULL

##' Is one constituent study observed significant, in favour of treatment, at a given alpha level? (Not intended for end user)
##' @param r_events_control Number of events in untreated
##' @param r_events_treated Number of events in treated group
##' @param n_sample_size_control Sample size in untreated group
##' @param n_sample_size_treated Sample size in treated group
##' @param alpha Type-1 error rate.
##' @return Number, 1 if positive, 0 if negative.
test.one.treated <- function(r_events_control,r_events_treated, n_sample_size_control,n_sample_size_treated, alpha=0.05){
  treated_nonevents<-n_sample_size_treated-r_events_treated    
  control_nonevents<-n_sample_size_control-r_events_control   
  cd<- matrix(c(r_events_treated,treated_nonevents,
                r_events_control, control_nonevents), nrow=c(2),ncol=c(2))
  test1<-fisher.test(cd)
  positive1<-ifelse(test1$p.value<alpha & test1$estimate<1,1,0)  
  return(positive1)
}
NULL

##' Observed number of positive by testing for significance from observed findings at given alpha level (not intended for end user).
##'  @param vec_r_events_control an ordered vector of number of events in the untreated group of constituent studies from a meta-analysis.
##'  @param vec_r_events_treated an ordered vector of number of events in the treated group of constituent studies from a meta-analysis.
##'  @param vec_n_sample_size_control an ordered vector of the number of participants in the untreated group.
##'  @param vec_n_sample_size_treated an ordered vector of the number of participants in the treated group.
##'  @param alpha Type-1 error rate.
##'  @return Vector of results, 1 if positive, 0 if negative.
test.n.treated <- function(vec_r_events_control, vec_r_events_treated, 
                           vec_n_sample_size_control, vec_n_sample_size_treated, alpha) {
  xx <- as.vector(1:length(vec_r_events_control))
  for (i in 1:length(xx)) xx[i] <- test.one.treated(vec_r_events_control[i], vec_r_events_treated[i],
                                                    vec_n_sample_size_control[i], vec_n_sample_size_treated[i], alpha)
  return(xx)
}
NULL
##' @title BMort
##' @name BMort
##' @docType data
##' @format Brugts meta-analysis study of effect of statins on overall mortality (BMJ 2009).
##' @source Brugts JJ, Yetgin T, Hoeks SE, Gotto AM, Shepherd J, Westendorp RGJ, et al. The benefits of statins in people without established cardiovascular disease but with cardiovascular risk factors: meta-analysis of randomised controlled trials. BMJ 2009;338
NULL



##' From a meta-analysis, analyse for publication bias.
##'  Calculates observed and expected number of positive studies and P for difference.
##'  @param vec_r_events_control an ordered vector of number of events in the untreated group of constituent studies from a meta-analysis.
##'  @param vec_r_events_treated an ordered vector of number of events in the treated group of constituent studies from a meta-analysis.
##'  @param vec_n_sample_size_control an ordered vector of the number of participants in the untreated group.
##'  @param vec_n_sample_size_treated an ordered vector of the number of participants in the treated group.
##'  @param n Number of iterations used to generate constituent study power; suggest use 10,000.
##'  @param low.alpha Lower limit of type-1 error rate used to judge whether constituent studies are positive; suggest 0.001.
##'  @param high.alpha Upper limit of type-1 error rate used to judge whether constituent studies are positive; suggest 0.3.
##'  @param by.alpha Interval of type-2 error rate at which observed and expected values and P for difference evaluated. 
##'  @return a dataframe with columns which include alpha level, observed number of positive studies, expected number, and P for difference, OR_hat (summary measure of effect for meta-analysis) with varying levels of significance for constituent studies.
##'  @export
##'  @examples
##'   data("BMort") ## Meta-analysis of statin use (Brugts 2009, BMJ)
##'  Btmort<-with(BMort, plot_chase_observed_expected(r_events_control, 
##'    r_events_treated, n_sample_size_control, n_sample_size_treated, n=10, 
##'    low.alpha=.001, high.alpha=0.3, by.alpha=0.01)) 
##' plot(Btmort$alpha, Btmort$observed,  type="l", las=1, lwd=2, xlim=c(.0001,0.3),
##'    xlab=c("Significance level"),  #### Brugts study mortality outcome; n set low for speed.
##'    ylab=c(""), main=c("(a) Brugts; all-cause mortality."))
##' lines(Btmort$alpha,Btmort$observed)
##' lines(Btmort$alpha,Btmort$expected, lty=3)
##' abline(v=0.05, lty=2)
##' par(new=TRUE)
##' plot(Btmort$alpha, Btmort$p.value, type="l", xlab="",lty=4,lwd=2, 
##' col="grey", axes=FALSE, ylab="")
##' abline(h=0.1, lty=2)
##' axis(4,las=1)
##' mtext(side=4,line=2.5,"P for difference")
plot_chase_observed_expected <-function(vec_r_events_control, vec_r_events_treated, 
                                          vec_n_sample_size_control, vec_n_sample_size_treated, 
                                          n, 
                                          low.alpha, high.alpha, by.alpha)  {
              metaOR<-NULL
              OR_hat<-NULL
              metaOR<-summary(meta.MH(vec_n_sample_size_treated,vec_n_sample_size_control,vec_r_events_treated,vec_r_events_control))
               OR_hat<-metaOR$MHci[[2]]
  alpha_list <- seq(low.alpha,high.alpha, by=by.alpha)
  b<- as.list(1:length(alpha_list))
  pb <- txtProgressBar(min = 0, max = length(b), style=3)
  for (i in 1:length(alpha_list)){
   
    b[[i]]<-test.n.treated(vec_r_events_control, vec_r_events_treated, 
                           vec_n_sample_size_control, vec_n_sample_size_treated, alpha=alpha_list[i])
  }
  
  a<-as.list(1:length(alpha_list))
  e<-as.vector(1:length(alpha_list))
    
  for (i in 1:length(a)) {
    a[[i]]<-ChisqTest_expect(vec_r_events_control, vec_n_sample_size_control,
                    vec_n_sample_size_treated, OR_hat, n, alpha=alpha_list[i], vec_pos=as.vector(b[[i]]))
    
  setTxtProgressBar(pb, i)
  }
  close(pb)
  #browser()
  a<-do.call(rbind,a)
  b<-do.call(rbind,b)
  f<-as.vector(1:nrow(a))
  for (i in 1:nrow(a)) {
    f[i]<-sum(b[i,])
  } 
  o<-rep(OR_hat,length(a))            
  d<-data.frame(cbind(alpha_list,a,f,o))
  names(d)<-c("alpha","p.value","expected","observed","OR_hat")
  return(d)
}
