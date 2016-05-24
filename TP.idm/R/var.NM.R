var.NM <-
function(data, ns, states, tr.states, s, t, B, e.times){  
  
  data_or<-data[order(data$time1),]
  
  n=nrow(data)
  s1<-c(rep(0,nrow(data)))
  for (i in 1:n){
    if(data_or$time1[i]>s) s1[i]=1
    else s1[i]=0
  }
  data_or$s1<-s1
  
  s2<-c(rep(0,nrow(data)))
  for (i in 1:n){
    if(data_or$time1[i]<=s & s<data_or$Stime[i]) s2[i]=1
    else s2[i]=0
  }
  data_or$s2<-s2
  
  data<-data_or
  
  split_ev <- strsplit(colnames(data), "  ") ## string splits names    
  s1.col <- which(sapply(split_ev, function(x) x[1]=="s1"))
  s2.col <- which(sapply(split_ev, function(x) x[1]=="s2"))
  ev1.col <- which(sapply(split_ev, function(x) x[1]=="event1"))
  ev.col <- which(sapply(split_ev, function(x) x[1]=="event"))
  t1.col <- which(sapply(split_ev, function(x) x[1]=="time1"))
  St.col <- which(sapply(split_ev, function(x) x[1]=="Stime"))
  
  ##  looks at noncensored individuals
  event.id <- which(sapply(split_ev, function(x) x[1]=="s1" | x[1]=="s2"))
  event.rows<- which(apply(data[, event.id, drop=FALSE], 1, function(x) any(x>0)))
  data.event <- data[event.rows, , drop=FALSE] 
  
  event.id.s <- which(sapply(split_ev, function(x) x[1]=="event1"))
  event.rows.s<- which(apply(data.event[, event.id.s, drop=FALSE], 1, function(x) any(x>0)))
  data.event.s <- data.event[event.rows.s, , drop=FALSE] 
  
  time1.s<-data.event.s[,t1.col][data.event.s[,ev1.col]==1]
  Stime.s<-data.event.s[,St.col][data.event.s[,ev.col]==1]
  
  ttime<-c(time1.s,Stime.s)
  times<- sort(unique(ttime))
  
  if (t=="last") t <- times[length(times)]
  
  
  # MC replicates:
  if (s==0){    
    # subset 1
    set_s1<-subset(data,data$s1==1)
    
    n_s1<-nrow(set_s1)
    
    ## subset s1 from initial state
    status1<-integer(n_s1)
    for(i in 1:n_s1){
      if(set_s1[i,ev1.col]==0 & set_s1[i,ev.col]==0) status1[i]=0
      if(set_s1[i,t1.col]<set_s1[i,St.col]) status1[i]=1
      if(set_s1[i,ev1.col]==1 & set_s1[i,ev.col]==1 & set_s1[i,t1.col]==set_s1[i,St.col]) status1[i]=1  
    }
    
    ## subset s1 from state 2:
    status<-integer(n_s1)
    for(i in 1:n_s1){
      if(set_s1[i,t1.col]<set_s1[i,St.col] & set_s1[i,ev1.col]==1 & set_s1[i,ev.col]==0) status[i]=0
      if(set_s1[i,t1.col]<set_s1[i,St.col] & set_s1[i,ev1.col]==1 & set_s1[i,ev.col]==1) status[i]=1
      if(set_s1[i,t1.col]==set_s1[i,St.col] & set_s1[i,ev1.col]==1 & set_s1[i,ev.col]==1) status[i]=1
    }
    
    res.ci.B.p11 <- res.ci.B.p12 <- res.ci.B.p13 <- matrix(0, nrow=length(e.times), ncol=B)
    
    for (j in 1:B){
      n <- dim(set_s1)[1]
      xx <- sample.int(n, size=n, replace=TRUE)
      
      ndata <- set_s1[xx,]
      
      
      ##########   S.Z (survival function for subset s1 in time1)  
      
      M1.s1.B <- cbind(ndata[,t1.col], status1[xx])
      colnames(M1.s1.B) <- c("time1", "status1")
      
      pes.s1.B <- KMW(M1.s1.B[,1], M1.s1.B[,2])
      
      ##########   S.T1 (survival function for subset s1 in Stime)  
      
      M1.s1.2.B <- cbind(ndata[,St.col], status[xx])
      colnames(M1.s1.2.B) <- c("Stime","status")
      
      pes.s1.2.B <- KMW(M1.s1.2.B[,1], M1.s1.2.B[,2])
      
      
      for(k in 1:length(e.times)){
        ti<-e.times[k]
        tts1<-which(M1.s1.B[,1]<=ti)
        
        tts1.2<-which(M1.s1.2.B[,1]<=ti)        
        
        #S1#
        res.ci.B.p11[k,j]<- 1 - sum(pes.s1.B[tts1])
        res.ci.B.p13[k,j]<- sum(pes.s1.2.B[tts1.2])
        res.ci.B.p12[k,j]<- 1- res.ci.B.p11[k,j] - res.ci.B.p13[k,j]         
      }      
    }
        
    var.b.p11<-apply(res.ci.B.p11,1,var)
    var.b.p12<-apply(res.ci.B.p12,1,var)
    var.b.p13<-apply(res.ci.B.p13,1,var)
    
    p.transitions<-c("1 1", "1 2", "1 3")
    var.b<-matrix(NA,nrow=length(e.times),ncol=length(p.transitions))   
    var.b<-cbind(var.b.p11, var.b.p12, var.b.p13)
    rownames(var.b)<-e.times
    colnames(var.b)<-paste(p.transitions)           
    
  } else {
    # s>0:    
        
    # subset 1
    set_s1<-subset(data,data$s1==1)
    # status 1:
    n_s1<-nrow(set_s1)
    status1<-integer(n_s1)
    for(i in 1:n_s1){
      if(set_s1[i,ev1.col]==0 & set_s1[i,ev.col]==0) status1[i]=0
      if(set_s1[i,t1.col]<set_s1[i,St.col]) status1[i]=1
      if(set_s1[i,ev1.col]==1 & set_s1[i,ev.col]==1 & set_s1[i,t1.col]==set_s1[i,St.col]) status1[i]=1  
    }
    
    ## subset s1 from state 2:
    status<-integer(n_s1)
    for(i in 1:n_s1){
      if(set_s1[i,t1.col]<set_s1[i,St.col] & set_s1[i,ev1.col]==1 & set_s1[i,ev.col]==0) status[i]=0
      if(set_s1[i,t1.col]<set_s1[i,St.col] & set_s1[i,ev1.col]==1 & set_s1[i,ev.col]==1) status[i]=1
      if(set_s1[i,t1.col]==set_s1[i,St.col] & set_s1[i,ev1.col]==1 & set_s1[i,ev.col]==1) status[i]=1
    }
    
    # subset s2:
    set_s2<-subset(data,data$s2==1)
    # status:
    n_s2<-nrow(set_s2)
    status2<-integer(n_s2)
    for(i in 1:n_s2){
      if(set_s2[i,t1.col]<set_s2[i,St.col] & set_s2[i,ev1.col]==1 & set_s2[i,ev.col]==0) status2[i]=0
      if(set_s2[i,t1.col]<set_s2[i,St.col] & set_s2[i,ev1.col]==1 & set_s2[i,ev.col]==1) status2[i]=1
    }
    
    res.ci.B.p11 <- res.ci.B.p12 <- res.ci.B.p13 <- matrix(0, nrow=length(e.times), ncol=B)
    res.ci.B.p22 <- res.ci.B.p23 <- matrix(0, nrow=length(e.times), ncol=B)
    
    
    for (j in 1:B){
      n <- dim(set_s1)[1]
      xx <- sample.int(n, size=n, replace=TRUE)
      
      ndata <- set_s1[xx,]
      
      ##########   S.Z (survival function for subset s1 in time1)  
      
      M1.s1.B <- cbind(ndata[,t1.col], status1[xx])
      colnames(M1.s1.B) <- c("time1", "status1")
      
      pes.s1.B <- KMW(M1.s1.B[,1], M1.s1.B[,2])      
      
      ##########   S.T1 (survival function for subset s1 in Stime)    
      
      M1.s1.2.B <- cbind(ndata[,St.col], status[xx])
      colnames(M1.s1.2.B) <- c("Stime","status")    
      
      pes.s1.2.B <- KMW(M1.s1.2.B[,1], M1.s1.2.B[,2])
      
      
      ##########   S.T2 (survival function for subset s2 in Stime)    
      
      n2 <- dim(set_s2)[1]
      xx2 <- sample.int(n2, size=n2, replace=TRUE)
      
      ndata2 <- set_s2[xx2,]
      
      M1.s2.B <- cbind(ndata2[,St.col], status2[xx2])
      colnames(M1.s2.B) <- c("Stime","status")    
      
      pes.s2.B <- KMW(M1.s2.B[,1], M1.s2.B[,2])

      
      #############   Transiton probabilities:
      
      
      for(k in 1:length(e.times)){
        
        ti<-e.times[k]
        tts1<-which(M1.s1.B[,1]<=ti)
        
        tts1.2<-which(M1.s1.2.B[,1]<=ti)
        
        tts2<-which(M1.s2.B[,1]<=ti)
        
        #S1#        
        res.ci.B.p11[k,j]<- 1 - sum(pes.s1.B[tts1])
        res.ci.B.p13[k,j]<- sum(pes.s1.2.B[tts1.2])
        res.ci.B.p12[k,j]<- 1 - res.ci.B.p11[k,j]- res.ci.B.p13[k,j]
        
        #S2#        
        #p22 y p23#
        res.ci.B.p23[k,j]<- sum(pes.s2.B[tts2])
        res.ci.B.p22[k,j]<- 1 - sum(pes.s2.B[tts2])
      }      
    }
    
    var.b.p11<-apply(res.ci.B.p11,1,var)
    var.b.p12<-apply(res.ci.B.p12,1,var)
    var.b.p13<-apply(res.ci.B.p13,1,var)
    var.b.p22<-apply(res.ci.B.p22,1,var)
    var.b.p23<-apply(res.ci.B.p23,1,var)
    
    p.transitions<-c("1 1", "1 2", "1 3", "2 2", "2 3")
    var.b<-matrix(NA,nrow=length(e.times),ncol=length(p.transitions))   
    var.b<-cbind(var.b.p11, var.b.p12, var.b.p13, var.b.p22, var.b.p23)
    rownames(var.b)<-e.times
    colnames(var.b)<-paste(p.transitions)   
    
  }
  
  res <- list(var.b=var.b)
  return(res)
  
}
