fun.AJ <-
function(ns,states,dNs,Ys,sum_dNs,s,t, event.times,initial.probs){
  
  if (t=="last") t <- event.times[length(event.times)]
  
  id_tr <- which(s<event.times& event.times<=t) ## location of those (s, t]
  
  l.id_tr <- length(id_tr)
  
  # to use in variance function:
  e.times.id_tr<-event.times[id_tr]
  
  ## Occupation Probability matrix
  OP <- matrix(NA, nrow=l.id_tr, ncol=length(states))
  rownames(OP) <- rownames(dNs)[id_tr]; colnames(OP) <- states
  all.dA <- all.I.dA <- array(dim=c(length(states), length(states), l.id_tr),
                              dimnames=list(rows=states,cols=states, time=rownames(dNs)[id_tr]))
  
  ## Transition Probability matrix
  cum.prod <- diag(ns)
  rownames(cum.prod) <- states
  
  TP.AJs <- array(dim=c(ns, ns, l.id_tr),dimnames=list(rows=states, cols=states,
                                                       dim=rownames(dNs)[id_tr]))
  
  dNs.id_tr<-dNs[id_tr,]
  Ys.id_tr<-Ys[id_tr,]
  sum_dNs.id_tr<-sum_dNs[id_tr,]
  
  
  for(i in 1:l.id_tr){
    
    I.dA <- diag(ns) ## creates trans matrix for current time
    
    dA <- matrix(0, nrow=ns, ncol=ns)
    
    colnames(I.dA) <- rownames(I.dA) <- colnames(dA) <- rownames(dA) <- states
    
    i_tr <- which(dNs.id_tr[i, , drop=FALSE]>0)  ## indicator for kind of transition at time i
    dNs.event.names <- colnames(dNs.id_tr)[i_tr] ## gets names of transitions (ie:  dN##)
    split_dNs.event <- strsplit(dNs.event.names, " ")    ## splits title of dN##
    st.start <- sapply(split_dNs.event, function(x) x[2])
    st.end <- sapply(split_dNs.event, function(x) x[3])  ## start & stop states as character strings
    i_tr.s <- matrix(as.character(c(st.start, st.end)), ncol=2)
    i_tr.s2 <- matrix(as.character(c(st.start, st.start)), ncol=2)
    
    dA[i_tr.s] <- dNs.id_tr[i, i_tr]/Ys.id_tr[i, paste("Y", st.start)]
    if (length(i_tr)==1) {
      dA[st.start, st.start] <- -dNs.id_tr[i, i_tr]/Ys.id_tr[i, paste("Y", st.start)]
    } else {
      dA[i_tr.s2] <- -rowSums(dA[st.start, ])
    }
    
    I.dA <- I.dA + dA ## I+dA (transition) matrix
    
    all.dA[, , i] <- dA     ## stores all dA matrices
    all.I.dA[, , i] <- I.dA ## array for storing all tran matrices
    
    cum.prod <- cum.prod %*% I.dA
    TP.AJs[,,i] <- cum.prod
  } ## end loop
  
  if (s==0){    
    OP[i, ] <- initial.probs%*%TP.AJs[, , i] ## state occupation probabilities 
    op<-OP[i,]   
    
    p.transitions<-c("1 1", "1 2", "1 3")
    
    all.est<-array(0,dim=c(l.id_tr,1,length(p.transitions)),
                   dimnames=list(rows=rownames(dNs)[id_tr],cols="probs",trans=p.transitions))
    
    for(j in 1:length(p.transitions)){
      idx<-unlist(strsplit(p.transitions[j]," "))
      all.est[,1,j]<-TP.AJs[idx[1],idx[2],]
    }
    
    
    res <- list(probs=op,all.est=all.est,TP.AJs=TP.AJs,all.I.dA=all.I.dA, dNs.id_tr=dNs.id_tr, 
                Ys.id_tr=Ys.id_tr, sum_dNs.id_tr=sum_dNs.id_tr, e.times.id_tr=e.times.id_tr,p.trans=p.transitions,t=t)
    return(res)
  } else{
    p.transitions<-c("1 1", "1 2", "1 3", "2 2", "2 3") 
    
    all.est<-array(0,dim=c(l.id_tr,1,length(p.transitions)),
                   dimnames=list(rows=rownames(dNs)[id_tr],cols="probs",trans=p.transitions))
    
    for(j in 1:length(p.transitions)){
      idx<-unlist(strsplit(p.transitions[j]," "))
      all.est[,1,j]<-TP.AJs[idx[1],idx[2],]
    }
    
    
    res <- list(probs=cum.prod,all.est=all.est,TP.AJs=TP.AJs,all.I.dA=all.I.dA, dNs.id_tr=dNs.id_tr, 
                Ys.id_tr=Ys.id_tr, sum_dNs.id_tr=sum_dNs.id_tr, e.times.id_tr=e.times.id_tr,p.trans=p.transitions,t=t)
    return(res)
  }
  
  
}
