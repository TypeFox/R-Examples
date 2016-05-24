out.stats<-function(stat,model.type){
  dec=3

  if (model.type=='grri'){
    R1.m1=stat[[1]][1,]
    R2.m1=stat[[1]][2,]
    Im.m1=stat[[1]][3,]
    R1.m2=stat[[2]][1,]
    R2.m2=stat[[2]][2,]
    Im.m2=stat[[2]][3,]
    DV.m1=stat[[3]]
    DV.m2=stat[[4]]

    sum1=matrix(c(
      exp(R1.m1[1]),exp(R1.m1[1]-1.96*R1.m1[2]),exp(R1.m1[1]+1.96*R1.m1[2]),pchisq(R1.m1[3]^2,1,lower.tail=FALSE),
      exp(R2.m1[1]),exp(R2.m1[1]-1.96*R2.m1[2]),exp(R2.m1[1]+1.96*R2.m1[2]),pchisq(R2.m1[3]^2,1,lower.tail=FALSE),
      exp(Im.m1[1]),exp(Im.m1[1]-1.96*Im.m1[2]),exp(Im.m1[1]+1.96*Im.m1[2]),pchisq(Im.m1[3]^2,1,lower.tail=FALSE),
      pchisq(DV.m1,3,lower.tail=FALSE),
      exp(R1.m2[1]),exp(R1.m2[1]-1.96*R1.m2[2]),exp(R1.m2[1]+1.96*R1.m2[2]),pchisq(R1.m2[3]^2,1,lower.tail=FALSE),
      exp(R2.m2[1]),exp(R2.m2[1]-1.96*R2.m2[2]),exp(R2.m2[1]+1.96*R2.m2[2]),pchisq(R2.m2[3]^2,1,lower.tail=FALSE),
      exp(Im.m2[1]),exp(Im.m2[1]-1.96*Im.m2[2]),exp(Im.m2[1]+1.96*Im.m2[2]),pchisq(Im.m2[3]^2,1,lower.tail=FALSE),
      pchisq(DV.m2,3,lower.tail=FALSE)),
      ncol=13, byrow=TRUE)

    sum1=cbind(c(1,2),data.frame(sum1))
    sum1=data.frame(sum1)
    names(sum1)=c("Model",
                  "R1","LCI","UCI","p-value",
                  "R2","LCI","UCI","p-value",
                  "T","LCI","UCI","p-value",
                  "LRT p-value")

    sum2=matrix(c(
      sum1[1,1],
      paste(signif(sum1[1,2],dec),"(",signif(sum1[1,3],dec),",",signif(sum1[1,4],dec),")",sep=""),
      signif(sum1[1,5],dec),
      paste(signif(sum1[1,6],dec),"(",signif(sum1[1,7],dec),",",signif(sum1[1,8],dec),")",sep=""),
      signif(sum1[1,9],dec),
      paste(signif(sum1[1,10],dec),"(",signif(sum1[1,11],dec),",",signif(sum1[1,12],dec),")",sep=""),
      signif(sum1[1,13],dec),
      signif(sum1[1,14],dec),
      sum1[2,1],
      paste(signif(sum1[2,2],dec),"(",signif(sum1[2,3],dec),",",signif(sum1[2,4],dec),")",sep=""),
      signif(sum1[2,5],dec),
      paste(signif(sum1[2,6],dec),"(",signif(sum1[2,7],dec),",",signif(sum1[2,8],dec),")",sep=""),
      signif(sum1[2,9],dec),
      paste(signif(sum1[2,10],dec),"(",signif(sum1[2,11],dec),",",signif(sum1[2,12],dec),")",sep=""),
      signif(sum1[2,13],dec),
      signif(sum1[2,14],dec)),
      byrow=T,ncol=8)

    sum2=data.frame(sum2)
    names(sum2)=c("Model",
                  "R1(95%CI)","p-value",
                  "R2(95%CI)","p-value",
                  "T(95%CI)","p-value",
                  "LRT p-value")

  }else if (model.type=='gdi'|model.type=='domi'){
    R1.m1=stat[[1]][1,]
    Im.m1=stat[[1]][2,]
    R1.m2=stat[[2]][1,]
    Im.m2=stat[[2]][2,]
    DV.m1=stat[[3]]
    DV.m2=stat[[4]]

    sum1=matrix(c(
      exp(R1.m1[1]),exp(R1.m1[1]-1.96*R1.m1[2]),exp(R1.m1[1]+1.96*R1.m1[2]),pchisq(R1.m1[3]^2,1,lower.tail=FALSE),
      exp(Im.m1[1]),exp(Im.m1[1]-1.96*Im.m1[2]),exp(Im.m1[1]+1.96*Im.m1[2]),pchisq(Im.m1[3]^2,1,lower.tail=FALSE),
      pchisq(DV.m1,2,lower.tail=FALSE),
      exp(R1.m2[1]),exp(R1.m2[1]-1.96*R1.m2[2]),exp(R1.m2[1]+1.96*R1.m2[2]),pchisq(R1.m2[3]^2,1,lower.tail=FALSE),
      exp(Im.m2[1]),exp(Im.m2[1]-1.96*Im.m2[2]),exp(Im.m2[1]+1.96*Im.m2[2]),pchisq(Im.m2[3]^2,1,lower.tail=FALSE),
      pchisq(DV.m2,2,lower.tail=FALSE)),
      ncol=9, byrow=TRUE)

    sum1=cbind(c(1,2),data.frame(sum1))
    sum1=data.frame(sum1)
    names(sum1)=c("Model",
                  "R","LCI","UCI","p-value",
                  "T","LCI","UCI","p-value",
                  "LRT p-value")

    sum2=matrix(c(
      sum1[1,1],
      paste(signif(sum1[1,2],dec),"(",signif(sum1[1,3],dec),",",signif(sum1[1,4],dec),")",sep=""),
      signif(sum1[1,5],dec),
      paste(signif(sum1[1,6],dec),"(",signif(sum1[1,7],dec),",",signif(sum1[1,8],dec),")",sep=""),
      signif(sum1[1,9],dec),
      signif(sum1[1,10],dec),
      sum1[2,1],
      paste(signif(sum1[2,2],dec),"(",signif(sum1[2,3],dec),",",signif(sum1[2,4],dec),")",sep=""),
      signif(sum1[2,5],dec),
      paste(signif(sum1[2,6],dec),"(",signif(sum1[2,7],dec),",",signif(sum1[2,8],dec),")",sep=""),
      signif(sum1[2,9],dec),
      signif(sum1[2,10],dec)),
      byrow=T,ncol=6)

    sum2=data.frame(sum2)
    names(sum2)=c("Model",
                  "R(95%CI)","p-value",
                  "T(95%CI)","p-value",
                  "LRT p-value")

    }else if (model.type=='grr'){
      R1.m1=stat[[1]][1,]
      R2.m1=stat[[1]][2,]
      R1.m2=stat[[2]][1,]
      R2.m2=stat[[2]][2,]
      DV.m1=stat[[3]]
      DV.m2=stat[[4]]

      sum1=matrix(c(
        exp(R1.m1[1]),exp(R1.m1[1]-1.96*R1.m1[2]),exp(R1.m1[1]+1.96*R1.m1[2]),pchisq(R1.m1[3]^2,1,lower.tail=FALSE),
        exp(R2.m1[1]),exp(R2.m1[1]-1.96*R2.m1[2]),exp(R2.m1[1]+1.96*R2.m1[2]),pchisq(R2.m1[3]^2,1,lower.tail=FALSE),
        pchisq(DV.m1,2,lower.tail=FALSE),
        exp(R1.m2[1]),exp(R1.m2[1]-1.96*R1.m2[2]),exp(R1.m2[1]+1.96*R1.m2[2]),pchisq(R1.m2[3]^2,1,lower.tail=FALSE),
        exp(R2.m2[1]),exp(R2.m2[1]-1.96*R2.m2[2]),exp(R2.m2[1]+1.96*R2.m2[2]),pchisq(R2.m2[3]^2,1,lower.tail=FALSE),
        pchisq(DV.m2,2,lower.tail=FALSE)),
        ncol=9, byrow=TRUE)

      sum1=cbind(c(1,2),data.frame(sum1))
      sum1=data.frame(sum1)
      names(sum1)=c("Model",
                    "R1","LCI","UCI","p-value",
                    "R2","LCI","UCI","p-value",
                    "LRT p-value")

      sum2=matrix(c(
        sum1[1,1],
        paste(signif(sum1[1,2],dec),"(",signif(sum1[1,3],dec),",",signif(sum1[1,4],dec),")",sep=""),
        signif(sum1[1,5],dec),
        paste(signif(sum1[1,6],dec),"(",signif(sum1[1,7],dec),",",signif(sum1[1,8],dec),")",sep=""),
        signif(sum1[1,9],dec),
        signif(sum1[1,10],dec),
        sum1[2,1],
        paste(signif(sum1[2,2],dec),"(",signif(sum1[2,3],dec),",",signif(sum1[2,4],dec),")",sep=""),
        signif(sum1[2,5],dec),
        paste(signif(sum1[2,6],dec),"(",signif(sum1[2,7],dec),",",signif(sum1[2,8],dec),")",sep=""),
        signif(sum1[2,9],dec),
        signif(sum1[2,10],dec)),
        byrow=T,ncol=6)

      sum2=data.frame(sum2)
      names(sum2)=c("Model",
                    "R1(95%CI)","p-value",
                    "R2(95%CI)","p-value",
                    "LRT p-value")

    }else if (model.type=='gd'|model.type=='dom'){
      R1.m1=stat[[1]]
      R1.m2=stat[[2]]
      DV.m1=stat[[3]]
      DV.m2=stat[[4]]

      sum1=matrix(c(
        exp(R1.m1[1]),exp(R1.m1[1]-1.96*R1.m1[2]),exp(R1.m1[1]+1.96*R1.m1[2]),pchisq(R1.m1[3]^2,1,lower.tail=FALSE),
        pchisq(DV.m1,1,lower.tail=FALSE),
        exp(R1.m2[1]),exp(R1.m2[1]-1.96*R1.m2[2]),exp(R1.m2[1]+1.96*R1.m2[2]),pchisq(R1.m2[3]^2,1,lower.tail=FALSE),
        pchisq(DV.m2,1,lower.tail=FALSE)),
        ncol=5, byrow=TRUE)

      sum1=cbind(c(1,2),data.frame(sum1))
      sum1=data.frame(sum1)
      names(sum1)=c("Model",
                    "R","LCI","UCI","p-value",
                    "LRT p-value")

      sum2=matrix(c(
        sum1[1,1],
        paste(signif(sum1[1,2],dec),"(",signif(sum1[1,3],dec),",",signif(sum1[1,4],dec),")",sep=""),
        signif(sum1[1,5],dec),
        signif(sum1[1,6],dec),
        sum1[2,1],
        paste(signif(sum1[2,2],dec),"(",signif(sum1[2,3],dec),",",signif(sum1[2,4],dec),")",sep=""),
        signif(sum1[2,5],dec),
        signif(sum1[2,6],dec)),
        byrow=T,ncol=4)

      sum2=data.frame(sum2)
      names(sum2)=c("Model",
                    "R(95%CI)","p-value",
                    "LRT p-value")
    }


  sum2

}
