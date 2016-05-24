prep.data.AJ <-
function(data, states, tr.states){
  ttime<-c(data$time1,data$Stime)
  times<-sort(unique(ttime)) # T_ik* indicates the instants of the k transitions
  
  
  # matrix with posible transitions:
  mat_w<-matrix(FALSE,nrow=3,ncol=3)
  mat_w[1,2:3]<-TRUE
  mat_w[2,3]<-TRUE
  
  # matrix with all transitions (including censoring transitions):
  mat_c<-matrix(TRUE,nrow=3,ncol=3)
  mat_c[2,1]<-FALSE
  mat_c[3,1:3]<-FALSE
  
  colnames(mat_w)<-rownames(mat_w)<-states
  colnames(mat_c)<-rownames(mat_c)<-states
  
  # output states whit contempling censoring (tr 11, tr 22)
  output.states.c<-lapply(1:dim(mat_c)[2],function(i){
    rownames(mat_c)[mat_c[i,]==TRUE]})
  
  # into states whitout censoring (tr 11, tr 22)
  into.states<-lapply(1:dim(mat_w)[2],function(i){
    colnames(mat_w)[mat_w[,i]==TRUE]
  })
  
  
  # possible transitions (including censoring)
  
  to<-c(output.states.c[[1]],output.states.c[[2]])
  from<-c(rep(tr.states[[1]],length(output.states.c[[1]])),rep(tr.states[[2]],length(output.states.c[[2]])))
  transitions<-paste("tr",from,to)
  
  
  # risk sets for each non-absorbing state (names)
  ys <- paste("Y", tr.states) 
  
  n=length(data$time1)
  m=length(times)
  
  # number of patientes with time1 = time k-transition (if event1=0 & event=0 --> 1 cens)
  g11<-integer(m)
  for(i in 1:m){
    ss11<-subset(data,data$event1==0 & data$event==0)
    g11[i]=sum(ss11$time1==times[i])
  }
  
  # number of patients with time1<time k-transition (for transition 1->3)
  g13<-integer(m)
  for (i in 1:m){
    ss13<-subset(data,data$event1==1 & data$event==1 & data$time1==data$Stime)
    g13[i]=sum(ss13$time1==times[i])
  }
  
  # number of patients with time1<time k-transition (for transition 1->2)
  g12<-integer(m)
  for (i in 1:m){
    ss12<-subset(data,data$time1<data$Stime)
    g12[i]=sum(ss12$time1==times[i])
  }
  
  # number of patients with time1<time k-transition (for transition 2->2)
  g22<-integer(m)
  for (i in 1:m){
    ss22<-subset(data,data$event1==1 & data$event==0 & data$time1<data$Stime)
    g22[i]=sum(ss22$Stime==times[i])
  }
  
  # number of patients with time1<time k-transition (for transition 2->3)
  g23<-integer(m)
  for (i in 1:m){
    ss23<-subset(data,data$event1==1 & data$event==1 & data$time1<data$Stime)
    g23[i]=sum(ss23$Stime==times[i])
  }
  
  ##  matrix of number of transitions:
  dNs<-cbind(g11,g12,g13,g22,g23)
  colnames(dNs)<-transitions
  rownames(dNs)<-times
  
  ## matrix of total number of transitions from each non-absorbing state:
  sum_n1<-rowSums(dNs[,2:3])
  sum_n2<-dNs[,5]
  sum_dNs<-cbind(sum_n1,sum_n2)
  
  
  
  data$start.time<-0
  
  initial.state<-integer(length(data[,1]))
  for(i in 1:length(data[,1])){
    if(data$start.time[i]==0) initial.state[i]=1
    else initial.state[i]=2
  }
  initial.state<-factor(initial.state,levels=tr.states,labels=tr.states)
  
  initial_risk_set<-table(initial.state)
  
  
  
  
  
  i=1
  n<-initial_risk_set[i]
  name<-paste("y",i)
  
  if(length(into.states[[i]])>0)
    in.st<-paste("tr",into.states[[i]],i) else in.st<-NULL
  
  if(length(output.states.c[[i]])>0)
    out.st<-paste("tr",i,output.states.c[[i]]) else out.st<-NULL
  
  Y1<-c(n,n+cumsum(rowSums(dNs[,in.st,drop=FALSE]))-cumsum(rowSums(dNs[,out.st,drop=FALSE])))
  
  i=2
  n<-initial_risk_set[i]
  name<-paste("y",i)
  
  if(length(into.states[[i]])>0)
    in.st<-paste("tr",into.states[[i]],i) else in.st<-NULL
  
  if(length(output.states.c[[i]])>0)
    out.st<-paste("tr",i,output.states.c[[i]]) else out.st<-NULL
  
  Y2<-c(n,n+cumsum(rowSums(dNs[,in.st,drop=FALSE]))-cumsum(rowSums(dNs[,out.st,drop=FALSE])))
  
  Ys<-cbind(Y1[1:m],Y2[1:m])
  
  
  rownames(dNs) <- rownames(sum_dNs) <- rownames(Ys) <- times
  colnames(dNs) <- transitions
  colnames(Ys) <- ys
  colnames(sum_dNs) <- paste("from", tr.states)
  
  split_dNs <- strsplit(colnames(dNs), "  ") ## string splits names
  split_Ys <- strsplit(colnames(Ys), " ")  ## string split names of Ys
  split_sum_dNs <- strsplit(colnames(sum_dNs), " ") ## string splits names of dNs
  
  ## censored columns in dNs, Ys and sum_dNs counting matrix, needed for D-S est
  cens.dNs.id <- which(sapply(split_dNs, function(x) x[1]=="tr 1 1"|x[1]=="tr 2 2"))
  cens.Ys.id <- which(sapply(split_Ys, function(x) x[2]%in%tr.states))
  cens.sum_dNs.id <- which(sapply(split_sum_dNs, function(x) x[2]%in%tr.states))
  
  
  K <- vector(length=nrow(dNs))
  dN.cens <- rowSums(dNs[, cens.dNs.id, drop=FALSE])
  Y.cens <- rowSums(Ys[, cens.Ys.id, drop=FALSE]) ## those at risk of being censored
  N.Y.cens <- ifelse(dN.cens/Y.cens=="NaN", 0, dN.cens/Y.cens)
  colnames(N.Y.cens) <- NULL
  H.t <- cumsum(N.Y.cens) ## calculating the hazard
  k <- exp(-H.t)
  K <- c(1, k[-length(k)])
  
  # weighted for right censoring
  dNs.w <- dNs/K  ## D-S dNs
  Ys.w <- Ys/K  ## D-S Ys
  sum_dNs.w <- sum_dNs/K
  
  res <- list(dNs=dNs.w, Ys=Ys.w, sum_dNs=sum_dNs.w, transitions=transitions)
  return(res)
}
