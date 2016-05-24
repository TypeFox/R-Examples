TPidm <-
function(data, s, t="last", CI=FALSE, B=199, level=0.95, ci.transformation="linear", method="NM"){
  ## conditions:
  if(missing(data)) stop("Argument 'data' is missing with no default")
  if(!is.data.frame(data)) stop("Argument 'data' must be a data.frame")
  data<-data[,c("time1","event1","Stime","event")]
  if(length(names(data))!=4) stop("'data' must have 4 variables")
  if(sum(c("time1","event1","Stime","event")==names(data))!=4) stop("'data' must contain the rigth variables")
  if(!is.numeric(data$time1)) {stop ("'time1' must be a numeric vector, which describes the time to firt event")}
  if(!is.numeric(data$Stime)) {stop ("'Stime' must be a numeric vector, which describes the survival or censoring time")}
  n=length(data[,1])
  for (i in 1:n){
    if(data$event1[i]==0 & data$event[i]==1) {stop("Argument 'event' must be 0 or argument 'event1' must be 1")}
  }
  
  mint<-min(data$time1[data$event1==1])
  
  if(s<0){stop("Argument 's' must be 0")} 
  if(s>0 & s<=mint){stop("Argument 's' must be 0 or greater than min(data$time1)=",mint)}
  if(t<=s){stop("Argument 't' must be greater than 's'")}
  
  # choose method
  ci.transf.type<-c("linear","log","log-log","cloglog")
  q<-charmatch(ci.transformation,ci.transf.type,nomatch=0)
  # charmatch function returns:
  # value 0 (if ci.transformation != c("linear","log","log-log","cloglog") 
  if (q==0)
    stop("ci.transformation should be 'linear', 'log', 'log-log' or 'cloglog'")
  
  
  # states:
  states<-c("1","2","3")
  
  # ns== number of states
  ns<-length(states)
  
  # non-absorbing states
  tr.states <- states[!states=="3"]
  
  
  # choose method
  method.type<-c("NM","AJ")
  m<-charmatch(method,method.type,nomatch=0)
  # charmatch function returns:
  # value 0 (if method != c(NM,AJ)) 
  # value 1 (if method = NM) or value 2 (if method = AJ)
  if (m==0)
    stop("method should be 'NM' or 'AJ'")
  if (m==1){
    # the method is 'NM'
    # prepare data set and NM method:
    NM.est<- fun.NM(data, states, s, t, tr.states)
    
    if(CI==TRUE){
      # var using B Monte Carlo replicates:
      variances<- var.NM(data, ns, states, tr.states, s, t,  B,  NM.est$e.times)   
      
      # ci using B Monte Carlo replicates:
      ci<- ci.NM(s, t, NM.est$all.probs, variances$var.b, level, ci.transformation, NM.est$e.times)
      
      # results:
      res <- list(
        # states information:
        method=method,s=s,t=NM.est$t,states=states, ns=ns, tr.states=tr.states, 
        ci.transformation=ci.transformation,
        # event times:
        times=NM.est$e.times,
        # occupation or transition probabilities:
        #probs=NM.est$probs, #all.est=new.est$all.probs,
        # confidence intervals:
        probs=ci$CI, all.probs=ci$all.CI,
        # posible transitions:
        p.trans=NM.est$p.trans,CI=CI)
    }else{
      # results:
      res <- list(
        # states information:
        method=method,s=s,t=NM.est$t,states=states, ns=ns, tr.states=tr.states,
        ci.transformation=ci.transformation,
        # event times:
        times=NM.est$e.times,
        # occupation or transition probabilities:
        probs=round(NM.est$probs,7), all.probs=round(NM.est$all.est,7),
        # posible transitions:
        p.trans=NM.est$p.trans,CI=CI)
    }
  }
  if (m==2){
    # the method is 'AJ'
    # start time ==0
    data$start.time<-0
    
    # initial probabilities for each initial state
    i.state<-integer(length(data[,1]))
    for(i in 1:length(data[,1])){
      if(data$start.time[i]==0) i.state[i]=1
      if(data$start.time[i]==data$time1[i]) i.state[i]=2
      if(data$start.time[i]==data$Stime[i]) i.state[i]=3
    }
    i.state <- factor(i.state,levels=states,labels=states)
    
    initial.probs <- prop.table(table(i.state))
    
    
    # prepare data set to compute AJ method:
    ds.prep.AJ <- prep.data.AJ(data, states, tr.states)
    
    # reduces to event times:
    ds.event.AJ <- prep.data.event.AJ(ds.prep.AJ$dNs, ds.prep.AJ$Ys, ds.prep.AJ$sum_dNs, states, tr.states)
    event.times <- as.numeric(as.character(rownames(ds.event.AJ$dNs)))
    
    # AJ estimator:
    AJ.est <- fun.AJ(ns,states, ds.event.AJ$dNs, ds.event.AJ$Ys, ds.event.AJ$sum_dNs,
                     s, t,  event.times, initial.probs)
    
    
    if(CI==TRUE){
      # variances of AJ estimator:
      variances <- var.AJ(ns, states, AJ.est$dNs.id_tr, AJ.est$Ys.id_tr, AJ.est$sum_dNs.id_tr,
                          AJ.est$TP.AJs, AJ.est$all.I.dA, tr.states)
      
      # CI of AJ estimator:
      ci <- ci.AJ(s,t, level, ci.transformation, AJ.est$dNs.id_tr, AJ.est$TP.AJs,
                  variances$cov.AJs, AJ.est$e.times.id_tr)
      
      
      # results:
      res <- list(
        # states information:
        method=method,s=s,t=AJ.est$t,states=states, ns=ns, tr.states=tr.states,
        ci.transformation=ci.transformation,
        # event times:
        times=AJ.est$e.times.id_tr,
        # occupation or transition probabilities:
        #est=AJ.est$probs, #all.est=AJ.est$TP.AJs,
        # confidence intervals:
        probs=ci$CI, all.probs=ci$all.CI,
        # posible transitions:
        p.trans=AJ.est$p.trans,CI=CI)
    }else{      
      # results:
      res <- list(
        # states information:
        method=method,s=s,t=AJ.est$t,states=states, ns=ns, tr.states=tr.states, 
        ci.transformation=ci.transformation,
        # event times:
        times=AJ.est$e.times.id_tr,
        # occupation or transition probabilities:
        probs=round(AJ.est$probs,7),all.probs=round(AJ.est$all.est,7),
        # posible transitions:
        p.trans=AJ.est$p.trans,CI=CI)
    }
  }
  
  res$call<-match.call()
  class(res) = "TPidm"
  res
  
}
