ctc <-
function(pmat,parvec){  
          #check treatment cardinality condition; parvec contains the cardinality restrictions of resp. a1 and a2
          #a1 cardinality T=1 in a node ;   #a2 cardinality T=2 in a node
          cond<-sapply(1:nrow(pmat),function(kk,pmat,parvec){ifelse(sum(pmat[kk,1:2]>=parvec)==2,1,0)},pmat=pmat, parvec=parvec)  
          condvec<-ifelse(sum(cond)==nrow(pmat),1,0)
          #if condvec==1 then the cardinality conditions for each node after the split are met 
         return(condvec)}
