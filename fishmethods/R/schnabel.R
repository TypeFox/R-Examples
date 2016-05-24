schnabel<-function(catch=NULL, recaps=NULL, newmarks=NULL,alpha=0.05){
   marks<-as.data.frame(cbind(catch,recaps,newmarks))
   if(ncol(marks)<3) stop("Only 2 columns of data.")
   if(any(length(marks[,1])!=c(length(marks[,2]),length(marks[,3]))))
         stop("Length of columns unequal")
 
  #Schnabel Method Using Reciprocal Confidence Interval
  S<-length(marks[,1])
  marks[1,4]<-0
  marks[2:S,4]<-cumsum(marks[1:as.numeric(S-1),3])
  invvar<-sum(marks[,2])/(sum(marks[,1]*marks[,4])^2)
  Ns<-sum(marks[,4]*marks[,1])/sum(marks[,2])
  if(sum(marks[,2])<50){
    UCI<-sum(marks[,4]*marks[,1])/qpois(alpha/2,sum(marks[,2]))
    LCI<-sum(marks[,4]*marks[,1])/qpois(1-(alpha/2),sum(marks[,2]))
    label<-"Poisson"
   }
  if(sum(marks[,2])>=50){
    UCI<-1/((1/Ns)+qt(alpha/2,df=S-1)*sqrt(invvar))
    LCI<-1/((1/Ns)+qt(1-alpha/2,df=S-1)*sqrt(invvar))
    label<-"t"
    }
  N1<-data.frame(N=Ns,invSE=sqrt(invvar),LCI=LCI,UCI=UCI,CI_Distribution=label)
  rownames(N1)<-c("Schnabel")

#Schumacher and Eschmeyer regression method

  Nse<-sum(marks[,1]*marks[,4]^2)/sum(marks[,2]*marks[,4])
  invvarse<-(sum(marks[,2]^2/marks[,1])-((sum(marks[,2]*marks[,4]))^2)/sum(marks[,1]*marks[,4]^2))/(S-2)
  invSEse <-sqrt(invvarse/sum(marks[,1]*marks[,4]^2))
  LCIse<-1/((1/Nse)+qt(1-alpha/2,df=S-2)*invSEse)
  UCIse<-1/((1/Nse)-qt(1-alpha/2,df=S-2)*invSEse)
  N2<-data.frame(N=Nse,invSE=invSEse,LCI=LCIse,UCI=UCIse,CI_Distribution="t")
  rownames(N2)<-c("Schumacher-Eschmeyer") 
  return(rbind(N1,N2))
}


