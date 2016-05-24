fun.NM <-
function(data, states, s, t, tr.states){
  
  
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
  
  id_tr.s <- which(s<=times & times<=t) ## location of those [s, t]
  l.id_tr.s <- length(id_tr.s)
  
  e.times<-times[id_tr.s]
  
  
  if (s==0){
    p.transitions<-c("1 1", "1 2", "1 3")
    
    ##########   S.Z (survival function for subset s1 in time1)
    
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
    
    M1.s1<-cbind(set_s1[,t1.col], status1)
    colnames(M1.s1)<-c("time1","status1")
    
    
    pes.s1=KMW(M1.s1[,1],M1.s1[,2])
    
    
    ##########   S.T1 (survival function for subset s1 in Stime)
    
    ## subset s1 from state 2:
    status<-integer(n_s1)
    for(i in 1:n_s1){
      if(set_s1[i,t1.col]<set_s1[i,St.col] & set_s1[i,ev1.col]==1 & set_s1[i,ev.col]==0) status[i]=0
      if(set_s1[i,t1.col]<set_s1[i,St.col] & set_s1[i,ev1.col]==1 & set_s1[i,ev.col]==1) status[i]=1
      if(set_s1[i,t1.col]==set_s1[i,St.col] & set_s1[i,ev1.col]==1 & set_s1[i,ev.col]==1) status[i]=1
    }
    
    M1.s1.2<-cbind(set_s1[,St.col], status)
    colnames(M1.s1.2)<-c("Stime","status")
    
    pes.s1.2=KMW(M1.s1.2[,1],M1.s1.2[,2])
    
    
    
    #############   Occupation Probabilities:
    
    all.OP.new<-array(dim=c(1, length(states), length(id_tr.s)),
                      dimnames=list(rows="1",cols=states, time=e.times))
    
    OP.new<-c()
    
    for(i in 1:l.id_tr.s){
      ti<-e.times[i]
      tts1<-which(M1.s1[,1]<=ti)
      
      tts1.2<-which(M1.s1.2[,1]<=ti)
      
      
      #S1#
 
      all.OP.new[1,1,i]<- 1 - sum(pes.s1[tts1])
      all.OP.new[1,3,i]<- sum(pes.s1.2[tts1.2])
      all.OP.new[1,2,i]<- 1- all.OP.new[1,1,i] - all.OP.new[1,3,i]
      
    }
    
    OP.new<-all.OP.new[,,i]
    
    
  } else {
    
    p.transitions<-c("1 1", "1 2", "1 3", "2 2", "2 3")
    
    ##########   S.Z (survival function for subset s1 in time1)
    
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
    
    M1.s1<-cbind(set_s1[,t1.col], status1)
    colnames(M1.s1)<-c("time1","status1")    
    
    pes.s1=KMW(M1.s1[,1],M1.s1[,2])
    
    
    
    ##########   S.T1 (survival function for subset s1 in Stime)  
    
    ## subset s1 from state 2:
    status<-integer(n_s1)
    for(i in 1:n_s1){
      if(set_s1[i,t1.col]<set_s1[i,St.col] & set_s1[i,ev1.col]==1 & set_s1[i,ev.col]==0) status[i]=0
      if(set_s1[i,t1.col]<set_s1[i,St.col] & set_s1[i,ev1.col]==1 & set_s1[i,ev.col]==1) status[i]=1
      if(set_s1[i,t1.col]==set_s1[i,St.col] & set_s1[i,ev1.col]==1 & set_s1[i,ev.col]==1) status[i]=1
    }
    
    M1.s1.2<-cbind(set_s1[,St.col], status)
    colnames(M1.s1.2)<-c("Stime","status")    
    
    pes.s1.2=KMW(M1.s1.2[,1],M1.s1.2[,2])
        
    
    ##########   S.T2 (survival function for subset s2 in Stime)  
    
    # subset s2:
    set_s2<-subset(data,data$s2==1)
    # status:
    n_s2<-nrow(set_s2)
    status<-integer(n_s2)
    for(i in 1:n_s2){
      if(set_s2[i,t1.col]<set_s2[i,St.col] & set_s2[i,ev1.col]==1 & set_s2[i,ev.col]==0) status[i]=0
      if(set_s2[i,t1.col]<set_s2[i,St.col] & set_s2[i,ev1.col]==1 & set_s2[i,ev.col]==1) status[i]=1
    }
    
    M1.s2<-cbind(set_s2[,St.col], status)
    colnames(M1.s2)<-c("Stime","status")    
    
    pes.s2=KMW(M1.s2[,1],M1.s2[,2])
          
    
    #############   Transiton probabilities:
    
    all.TP.new<-array(dim=c(length(states), length(states), length(id_tr.s)),
                      dimnames=list(rows=states,cols=states, time=e.times))
    
    ns<-length(states)
    TP.new<-diag(ns)   
    
    for(i in 1:l.id_tr.s){
      ti<-e.times[i]
      tts1<-which(M1.s1[,1]<=ti)
      
      tts1.2<-which(M1.s1.2[,1]<=ti)
      
      tts2<-which(M1.s2[,1]<=ti)
      
      #S1#
      
      all.TP.new[1,1,i]<- 1 - sum(pes.s1[tts1])
      all.TP.new[1,3,i]<- sum(pes.s1.2[tts1.2])
      all.TP.new[1,2,i]<- 1- all.TP.new[1,1,i] - all.TP.new[1,3,i]
      
      #S2#
      all.TP.new[2,1,i]<-0
      
      #p22 y p23#
      all.TP.new[2,3,i]<- sum(pes.s2[tts2])
      all.TP.new[2,2,i]<- 1 - sum(pes.s2[tts2])
      
      #p31, p32#
      all.TP.new[3,1,i]<-all.TP.new[3,2,i]<-0
      #p33#
      all.TP.new[3,3,i]<-1
    }
    
    TP.new<-all.TP.new[,,i]
    
  }
  
  if (s==0){
    all.est<-array(0,dim=c(l.id_tr.s,1,length(p.transitions)),
                   dimnames=list(rows=e.times,cols="probs",trans=p.transitions))
    
    for(j in 1:length(p.transitions)){
      idx<-unlist(strsplit(p.transitions[j]," "))
      all.est[,1,j]<-all.OP.new[idx[1],idx[2],]
    }
    
    
    res <- list(probs=OP.new, all.est=all.est, all.probs=all.OP.new, e.times=e.times, p.trans=p.transitions,t=t)
    return(res)
  } else{
    all.est<-array(0,dim=c(l.id_tr.s,1,length(p.transitions)),
                   dimnames=list(rows=e.times,cols="probs",trans=p.transitions))
    
    for(j in 1:length(p.transitions)){
      idx<-unlist(strsplit(p.transitions[j]," "))
      all.est[,1,j]<-all.TP.new[idx[1],idx[2],]
    }
    
    
    res <- list(probs=TP.new, all.est=all.est, all.probs=all.TP.new, e.times=e.times, p.trans=p.transitions,t=t)
    return(res)
  }
  
}
