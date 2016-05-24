demosimhybrid <-
function(x,M=matrix(1,ncol=6, nrow=6),F=c(1,1,1,1,1,1))
  {
  options(digits=4)
    if (sum(x)!=1)
    stop ("The initial proportions must sum 1")
    ono<-function(x){
  L<- x%*%t(x)
  Crosses<- L*M
  Z<-  c((Crosses[1,1]+Crosses[1,4]+Crosses[4,1]),(Crosses[2,2]+Crosses[2,5]+Crosses[5,2]),(Crosses[1,2]+Crosses[2,1]), (Crosses[1,3]+Crosses[3,1]),(Crosses[2,3]+Crosses[3,2]), (Crosses[1,5]+Crosses[1,6]+Crosses[2,4]+Crosses[2,6]+Crosses[3,3]+Crosses[3,4]+Crosses[3,5]+Crosses[3,6]+Crosses[4,2]+Crosses[4,3]+Crosses[4,4]+Crosses[4,5]+Crosses[4,6]+ Crosses[5,1]+Crosses[5,3]+Crosses[5,4]+Crosses[5,5]+Crosses[5,6]+Crosses[6,1]+Crosses[6,2]+Crosses[6,3]+Crosses[6,4]+Crosses[6,5]+Crosses[6,6]))
  x<- ((Z * F)/sum(Z*F))
  print(x)
  if (x[1]> 0.999 | x[2] >0.999 |x[3] >0.999 |x[4] >0.999 |x[5] >0.999 |x[6] >0.999) return(x)
  else Recall(x)
  } 
  # Save output into a matrix 
  out<-capture.output(ono(x))
  matrixout<-t(sapply(strsplit(gsub(" +"," ",out)," "),function(l)as.numeric(l[2:7])))
  #Add initial generation
  matrixout2<-rbind(x,matrixout)
  l<-dim(matrixout2)[1]
  g<-rep("G",l)
  num<-(0:(l-1))
  gen<-paste(g,num, sep="")
  colnames(matrixout2)<-c("PA","PB","F1","BPA","BPB","Fx")
  rownames(matrixout2)<-gen
  matrixout3<-matrixout2[-l,]
  class( matrixout3) <- "demosim.hybrid"
  matrixout3
  }
