###############################################################################
#
#                  Survival Rate Estimators
#        Catch Curve, Chapman-Robson, and Heinecke Methods
#        No plus groups; see Chapman and Robson (1961) for plus-groups
#                    Data on rock bass (C&R 1961)
###############################################################################

agesurv<-function(age=NULL,full=NULL,last=NULL,estimate=c("s","z"),method=c("cc","he","cr","crcb","mil")){
     if(is.null(age)) 
         stop ("vector does not exist")
      if(!is.numeric(age)) 
         stop ("vector is not numeric")
      if(is.null(full)) 
         stop ("fully-recruited age not specified") 
      age<-age[!is.na(age)]
      st_obj<-as.data.frame(table(age))
      st_obj[,1]<-as.numeric(levels(st_obj[,1])[as.integer(st_obj[,1])])
      if(is.null(last)) last<-max(age) else last<-last 
      d<<-subset(st_obj,st_obj[,1]>=full & st_obj[,1]<=last)
        names(d)<-c("age","number")
      if(d[1,1]!=full) 
         stop ("Age specified as fully-recruited does not exist.")
      cnt<-1
     if(length(d[,1])<=2){
        print(paste("warning: only", length(d[,1]),"ages!!!"))
        rown<-length(method)*length(estimate)-1
      }
     if(length(d[,1])>2) 
        rown<-length(method)*length(estimate)  
     results<-data.frame(Method=rep("NA",rown),
                    Parameter=rep("NA",rown),
                    Estimate=rep(as.numeric(NA),rown),
                    SE=rep(as.numeric(NA),rown),stringsAsFactors=F)


    if(any(method=="cc")){ 
      if(length(d[,1])>2){   
	   cc<-lm(log(d[,2])~d[,1])
         Zcc<-coef(summary(cc))[2,1]*-1
         SEZcc<-round(coef(summary(cc))[2,2],3)
         Scc<-exp(coef(summary(cc))[2,1])
         SEScc<-Scc*SEZcc

      if(any(estimate=="s")){
         results[cnt,1]<-"Catch Curve"
         results[cnt,2]<-"S"
         results[cnt,3]<-round(Scc,2)
         results[cnt,4]<-round(SEScc,3)
         cnt<-cnt+1
       }
      if(any(estimate=="z")){
         results[cnt,1]<-"Catch Curve"
         results[cnt,2]<-"Z"
         results[cnt,3]<-round(Zcc,2)
         results[cnt,4]<-round(SEZcc,3)
         cnt<-cnt+1
       }
     }
   }
   if(any(method=="he")){
     	Sh<-1-d[,2][d[,1]==full]/sum(d[,2])
	SESh<-sqrt(Sh*(1-Sh)/sum(d[,2]))
	Zh<--log(1-d[,2][d[,1]==full]/sum(d[,2]))
	SEZh<-sqrt(((1-Sh)^2)/(sum(d[,2])*Sh))

      if(any(estimate=="s")){
         results[cnt,1]<-"Heinecke"
         results[cnt,2]<-"S"
         results[cnt,3]<-round(Sh,2)
         results[cnt,4]<-round(SESh,3)
         cnt<-cnt+1
       }
      if(any(estimate=="z")){
         results[cnt,1]<-"Heinecke"
         results[cnt,2]<-"Z"
         results[cnt,3]<-round(Zh,2)
         results[cnt,4]<-round(SEZh,3)
         cnt<-cnt+1
       }
    }
   if(any(method=="cr")){
       j<-seq(0,length(d[,1])-1,1)
       Scr<-sum(j*d[,2])/(sum(d[,2])+sum(j*d[,2])-1)
       SEScr<-sqrt(Scr*(Scr-(sum(j*d[,2])-1)/(sum(d[,2])+
                 sum(j*d[,2])-2)))
 	      Zcr<--log(Scr)
	      SEZcr<-sqrt(((1-Scr)^2)/(sum(d[,2])*Scr))          
       if(any(estimate=="s")){
          results[cnt,1]<-"Chapman-Robson"
          results[cnt,2]<-"S"
          results[cnt,3]<-round(Scr,2)
          results[cnt,4]<-round(SEScr,3)
          cnt<-cnt+1
        }
       if(any(estimate=="z")){
          results[cnt,1]<-"Chapman-Robson"
          results[cnt,2]<-"Z"
          results[cnt,3]<-round(Zcr,2)
          results[cnt,4]<-round(SEZcr,3)
          cnt<-cnt+1
        }
    }
     if(any(method=="crcb")){
       j<-seq(0,length(d[,1])-1,1)
       Scr<-sum(j*d[,2])/(sum(d[,2])+sum(j*d[,2])-1)
       SEScr<-sqrt(Scr*(Scr-(sum(j*d[,2])-1)/(sum(d[,2])+
            sum(j*d[,2])-2)))
       X<-sum(j*d[,2])
       n<-sum(d[,2])
       bc<-((n-1)*(n-2))/(n*(n+X-1)*(X+1))
       Zcr<--log(Scr)
       pN<-d[1,2]*Scr^j
       chi<-sum((d[,2]-pN)^2/pN)
       chat<-sqrt(chi/(length(d[,2])-1))
       SEZcr<-sqrt(((1-Scr)^2)/(sum(d[,2])*Scr))          
       if(any(estimate=="s")){
         results[cnt,1]<-"Chapman-Robson CB"
         results[cnt,2]<-"S"
         results[cnt,3]<-round(Scr,2)
         results[cnt,4]<-round(SEScr,3)
         cnt<-cnt+1
       }
       if(any(estimate=="z")){
         results[cnt,1]<-"Chapman-Robson CB"
         results[cnt,2]<-"Z"
         results[cnt,3]<-round(-log(Scr)-bc,2)
         results[cnt,4]<-round(SEZcr*chat,3)
         cnt<-cnt+1
       }
     }
  if(any(method=="mil")){
      max.age<-max(d[,1])
    extd<-rbind(d,cbind(age=(max.age+1):(2*max.age),number=rep(0,max.age)))
    wer<-glmer(number~age+(1|age),family=poisson,data=extd)
    if(any(estimate=="s")){
      results[cnt,1]<-"Millar glmer"
      results[cnt,2]<-"S"
      tf<-bt.log(summary(wer)$coefficients[2,1],summary(wer)$coefficients[2,2],sum(d[,2]))
      results[cnt,3]<-round(as.numeric(tf[1]),2)
      results[cnt,4]<-round(as.numeric(tf[4]),3)
      cnt<-cnt+1
    }
    if(any(estimate=="z")){
     results[cnt,1]<-"Millar glmer"
     results[cnt,2]<-"Z"
     results[cnt,3]<-round(abs(summary(wer)$coefficients[2,1]),2)
     results[cnt,4]<-round(summary(wer)$coefficients[2,2],3) 
    }
  
  }
   out<-list(results,d);names(out)<-c("results","data")
   return(out)
}

