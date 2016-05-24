#data<-read.csv("P:/Rwork/book/jensen2.csv")     

agesurvcl<-function(age=NULL,group=NULL,full=NULL,last=NULL){
 if(is.null(age)) 
         stop ("age vector does not exist")
 if(is.null(group)) 
         stop ("group vector does not exist")
 if(is.null(full)) 
         stop ("fully-recruited age not specified") 
 if(length(age)!=length(group)) 
         stop ("length of age and group vectors are different")

   if(is.null(last)) last<-max(age) else last<-last 
   d<-as.data.frame(cbind(age,group))
   d<-subset(d,d$age>=full &d$age<=last)
   sage<-as.data.frame(table(d$age))
   if(length(sage[,1])<=2){
        print(paste("warning: only", length(sage[,1]),"ages!!!"))
   }

   names(sage)<-c("age","freq")
   a<-sage$freq[1]/sum(sage$freq)				
   S<-1-a;Z<--1*log(S)						  

   get1<-as.data.frame(table(d$group,d$age))
   names(get1)<-c("group","age","number")
   get2<-as.data.frame(table(d$group))
   names(get2)<-c("group","total")

  sub<-merge(get1,get2,by.x="group",by.y="group")
  sub<-sub[sub$age==full,]

  s2<-sum((sub$number-a*sub$total)^2)/(nlevels(get1$group)-1)
  SEa<-sqrt(s2/(nlevels(get1$group)*mean(sub$total)^2))
  SEZ<-sqrt(SEa^2/(1-a)^2)
         ans<-NULL  
         ans<-matrix(NA,1L,7L)
         ans<-rbind(cbind(round(a,2),round(SEa,3)),
                            cbind(round(S,2),round(SEa,3)),
                            cbind(round(Z,2),round(SEZ,3)))   
         dimnames(ans)<-list(cbind("a","S","Z"),c("Estimate","SE")) 
 
  return(ans)
}    

#dd<-agesurvcl(age=data$age,group=data$group,full=0)


