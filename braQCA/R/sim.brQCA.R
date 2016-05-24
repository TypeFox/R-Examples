#' Simulation Application
#' 
#' Internal function to calculate the Bootstrapped Recommendation.
#' @import QCAGUI bootstrap
#' @importFrom graphics hist
#' @importFrom stats glm plogis predict quantile
#' @importFrom utils flush.console
#' @param qca.data the data frame of causal conditions.
#' @param outcome the outcome variable (object name) in the QCA data frame of causal conditions; \code{"OUT"} is the outcome variable for an application of QCA. Default set to \code{outcome="OUT"}.
#' @param conditions a set of causal conditions. Default set to \code{conditions=c("")}
#' @param sim number of simulations to run. Default set to \code{sim=10}.
#' @param ncut configurational n levels for inclusion. Default set to \code{ncut=2}.
#' @param type type of QCA application, \code{"crisp"} or \code{"fuzzy"} sets. Default set to \code{type = "crisp"}.
#' @param inclcut minimum sufficiency score for inclusion. Default set to \code{inclcut=""}.
#' @param neg.out [from QCAGUI package] ``Logical, use negation of outcome (ignored if data is a truth table object).'' Default set to \code{neg.out=F}.
#' @param verbose prints the system time used to run the simulation and the percent complete. Default set to \code{verbose=T}.
#' @return Simulation information later passed on to conf.table.
#' @export
sim.brQCA<-function(qca.data, outcome="OUT", conditions=c(""), sim=10, ncut=2, type="crisp", inclcut = "", neg.out=F, verbose=T){
ptm <- proc.time()

if (all(conditions == c(""))) {
  conditions <- names(qca.data[,!(names(qca.data) %in% outcome)]) #use all coniditions that are not the outcome, if no conditions are specified
}

if (inclcut == ""){
  inclcut<-seq(from=.5, to = 1, by=.01) #use precise consistency thresholds
}

pop<-dim(qca.data)[1]  #sample size
rows<-sim*length(ncut)*length(inclcut)*length(pop)*2 #total number of rows
out<-qca.data[,outcome] #outcome vector
qca.data<-qca.data[,!(names(qca.data) %in% outcome)] #matrix of causal conditions
data<-data.frame(CTH=0,CNTH=0,CPI=0,NTH=0,OUT=rep(NA,rows)) #empty data set to simulate into

kk<-0 #set counter to 0

for (j in 1:sim) {
  
  for (k in 1:length(inclcut)){
    
    for (n in ncut){
      
      
      
      kk<-kk+1
      
      data[kk,1]<-inclcut[k]
      data[kk,2]<-n
      data[kk,4]<-pop
      
      
      s.qca.data<-qca.data
      
      if (type=="crisp"){
        
        for (i in 1:length(qca.data)){ #simulate random causal conditions
          prob<-c(sum(qca.data[,i]==0)/(dim(qca.data)[1]),sum(qca.data[,i]==1)/dim(qca.data)[1]) #match distributions of data set
          s.qca.data[,i]<-sample(c(0,1),pop,prob=prob,replace=T)} 
        
        #simulate random outcome variable
        prob<-c(sum(out==0)/(length(out)),sum(out==1)/length(out)) #match distributions of data set
        s.qca.data$OUT<-sample(c(0,1),pop,prob=prob,replace=T)
      }
      
      if (type == "fuzzy"){
        
        for (i in 1:length(qca.data)){ #simulate random causal conditions
          ranges<-seq(from=0.1, to=1, by=.1) #better way to do this? could do a for loop
          prob<-hist(qca.data[,i])[[2]]/dim(qca.data)[1]
          s.qca.data[,i]<-sample(ranges,pop,prob=prob,replace=T)
        } 
        
        #simulate random outcome variable
        prob<-hist(out)[[2]]/length(out)
        s.qca.data$OUT<-sample(ranges,pop,prob=prob,replace=T)
      }
      
      ##########parsimonious
      
      parsimonious <- tryCatch( #trap error
        eqmcc(s.qca.data,  outcome=c("OUT"),  n.cut=n, incl.cut1=inclcut[k], include = "?", neg.out=neg.out,
              conditions= c(names(s.qca.data[,!(names(s.qca.data) %in% 'OUT')])),details = TRUE, show.cases = TRUE),
        error=function(e) e
      )
      
      if(!inherits(parsimonious, "error")){
        
        data[kk,5]<-1 #1 = it returned a configuration with random data!
        data[kk,3]<-0 # parsimonious solution
      }
      
      #if(grepl("None of the values",parsimonious)[1] | grepl("All combinations have been included into analysis", parsimonious)[1]){
      if(inherits(parsimonious, "error")){
        #REAL WORK
        data[kk,5]<-0 #0 = it can't find the pattern that isn't there!
        data[kk,3]<-0 # parsimonious solution
      }
      
      
      ########complex
      
      kk<-kk+1 # increment row
      
      data[kk,1]<-inclcut[k]
      data[kk,2]<-n
      data[kk,4]<-pop
      
      complex <- tryCatch( #trap error
        eqmcc(s.qca.data,  outcome=c("OUT"),  n.cut=n, incl.cut1=inclcut[k], neg.out=neg.out,
              conditions = c(names(s.qca.data[,!(names(s.qca.data) %in% 'OUT')])), details = TRUE, show.cases = TRUE),
        error=function(e) e
      )
      
      if(!inherits(complex, "error")){  
        data[kk,5]<-1 #1 = it returned a configuration with random data!
        data[kk,3]<-1 # complex solution
      }
      
      # if(grepl("None of the values",complex)[1] | grepl("All combinations have been included into analysis", complex)[1]){
      if(inherits(complex, "error")){  
        data[kk,5]<-0 #0 = it can't find the pattern that isn't there!
        data[kk,3]<-1 # complex solution
      }
      
      
      captureError<-tryCatch(truthTable(s.qca.data,  outcome=c("OUT"),  n.cut=10, incl.cut1=inclcut[k], include = "?", neg.out=neg.out,
                                        conditions= c(names(s.qca.data[,!(names(s.qca.data) %in% 'OUT')])),details = TRUE, show.cases = TRUE)[[1]][,1], error=function(e) e)
      
      # if (grepl("replacement has 0 items",captureError)){data$OUT[kk]<-NA}
      
      #if (length(captureError)<=2){data$OUT[kk]<-NA}
      
      
      
      
      if (verbose == T){
        print(paste(round(100*kk/rows, digits=2),"% done", sep=""))
        print(proc.time()-ptm)
        flush.console()}
    }}}
return(data)
}

