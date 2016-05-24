`orderChromosome` <-function(x)
 {
   temp<-rep(NA,13)
   for (i in 10:22)
     temp[i-9]<-grep(i,x)
   temp2<-rep(NA,9)
   for (i in 1:9){
      aux<-grep(i,x)
      temp2[i]<-aux[!aux%in%temp]
   }
   temp3<-grep("X",x)
   res<-c(temp2,temp,temp3)
   res
  }

