dsrwon<-function(x,y) {
  # x1 et y1 doivent être des matrices pour assurer le type numérique dans inters()
  x1<-as.matrix(x$DempsterRule) 
  y1<-as.matrix(y$DempsterRule)
  nc<-ncol(x1)-1
  N12<-inters(x1,y1)
  V12<-combmasses(x1,y1)
  ## transformer tableau : MxN lignes by K hypothesis
  N12<-aperm(N12,c(2,1,3))
  N12<-array(c(N12),c(dim(N12)[1],prod(dim(N12)[-1]))) ##tab k x (MxN)
  N12<-aperm(N12,c(2,1)) 
  ## remove duplicates from the table
  W1<- doubles(N12)
  ## idendify  contributions to each subset
  I12<-dotprod(W1,aperm(N12,c(2,1)),g="&",f="==")
  ## calcul de la masse totale de chaque ss-ensemble 
  MAC<-apply(I12*t(array(V12,dim(t(I12)))),1,sum)
  ## order the subsets to check if and where the empty subset is there
  TRI<-order(apply(W1,1,sum))
  z<- sum(W1[TRI[1],])
  ## identify the indice of the empty set 
  if (z==0) IVIDE<-TRI[1] else IVIDE<-0 
 ## if (z==0) W2<-rbind(rep(0,times=nc),matrix(W1[TRI[-1],],ncol=nc))  else W2<-matrix(W1,ncol=nc)
  W2<-matrix(W1,ncol=nc)
  ## calcul de MVIDE
  if (z==0) MVIDE<-MAC[IVIDE] else MVIDE<-0
  con12<-1-(1-x$con)*(1-y$con)
  con<-1-(1-con12)*(1-MVIDE)
   ## result
  W2<-cbind(matrix(MAC,ncol=1),W2)
  colnames(W2)<-colnames(x1)
 # result: list m12 of four elements:
 # 1: un-normalized Dempster's rule of combination
 # 2: table of intersections
 # 3: indices for the sort
 # 4: measure of conflict between beliefs
  list(combination=W2,I12=I12,TRI=TRI,con=con)
  }
 