SLD<-function(fac, lev){
 if(fac >12) 
  stop("This function only handles 12 or less factors.","\n")
 if(lev >5)
  stop("This function only works for 5 levels or less")
             
cnames<-paste("x",1:fac, sep="")
   while (fac==2) {
SL<-expand.grid(seq(0,1,by=(1/lev)), seq(0,1,by=(1/lev)))
names(SL)<-cnames
SL<-SL[SL$x1+SL$x2==1, ]
dimS<-dim(SL)
rows<-dimS[1]
rnames<-paste(1:rows)
rownames(SL)<-rnames
return(SL)
                 }
  while (fac==3)   {
SL<-expand.grid(seq(0,1,by=(1/lev)), seq(0,1,by=(1/lev)),
      seq(0,1,by=(1/lev)))
names(SL)<-cnames
SL<-SL[SL$x1+SL$x2+SL$x3==1, ]
dimS<-dim(SL)
rows<-dimS[1]
rnames<-paste(1:rows)
rownames(SL)<-rnames
return(SL)
                  }
  while (fac==4)   {
SL<-expand.grid(seq(0,1,by=(1/lev)), seq(0,1,by=(1/lev)),
      seq(0,1,by=(1/lev)), seq(0,1,by=(1/lev)))
names(SL)<-cnames
SL<-SL[SL$x1+SL$x2+SL$x3+SL$x4==1, ]
dimS<-dim(SL)
rows<-dimS[1]
rnames<-paste(1:rows)
rownames(SL)<-rnames
return(SL)
                  }
  while (fac==5)   {
SL<-expand.grid(seq(0,1,by=(1/lev)), seq(0,1,by=(1/lev)),
      seq(0,1,by=(1/lev)), seq(0,1,by=(1/lev)),
      seq(0,1,by=(1/lev)))
names(SL)<-cnames
SL<-SL[SL$x1+SL$x2+SL$x3+SL$x4+SL$x5==1, ]
dimS<-dim(SL)
rows<-dimS[1]
rnames<-paste(1:rows)
rownames(SL)<-rnames
return(SL)
                  }
  while (fac==6)   {
SL<-expand.grid(seq(0,1,by=(1/lev)), seq(0,1,by=(1/lev)),
      seq(0,1,by=(1/lev)), seq(0,1,by=(1/lev)),
      seq(0,1,by=(1/lev)), seq(0,1,by=(1/lev)))
names(SL)<-cnames
SL<-SL[SL$x1+SL$x2+SL$x3+SL$x4+SL$x5+SL$x6==1, ]
dimS<-dim(SL)
rows<-dimS[1]
rnames<-paste(1:rows)
rownames(SL)<-rnames
return(SL)
                  }
  while (fac==7)   {
SL<-expand.grid(seq(0,1,by=(1/lev)), seq(0,1,by=(1/lev)),
      seq(0,1,by=(1/lev)), seq(0,1,by=(1/lev)),
      seq(0,1,by=(1/lev)), seq(0,1,by=(1/lev)),
      seq(0,1,by=(1/lev)))

names(SL)<-cnames
SL<-SL[SL$x1+SL$x2+SL$x3+SL$x4+SL$x5+SL$x6+SL$x7==1, ]
dimS<-dim(SL)
rows<-dimS[1]
rnames<-paste(1:rows)
rownames(SL)<-rnames
SL<-data.frame(SL)
return(SL)
                  }
  while (fac==8)   {
SL<-expand.grid(seq(0,1,by=(1/lev)), seq(0,1,by=(1/lev)),
      seq(0,1,by=(1/lev)), seq(0,1,by=(1/lev)),
      seq(0,1,by=(1/lev)), seq(0,1,by=(1/lev)),
      seq(0,1,by=(1/lev)),seq(0,1,by=(1/lev)))

names(SL)<-cnames
SL<-SL[SL$x1+SL$x2+SL$x3+SL$x4+SL$x5+SL$x6+SL$x7+SL$x8==1, ]
dimS<-dim(SL)
rows<-dimS[1]
rnames<-paste(1:rows)
rownames(SL)<-rnames
SL<-data.frame(SL)
return(SL)
                  }
  while (fac==9)   {
SL<-expand.grid(seq(0,1,by=(1/lev)), seq(0,1,by=(1/lev)),
      seq(0,1,by=(1/lev)), seq(0,1,by=(1/lev)),
      seq(0,1,by=(1/lev)), seq(0,1,by=(1/lev)),
      seq(0,1,by=(1/lev)), seq(0,1,by=(1/lev)),
      seq(0,1,by=(1/lev)))

names(SL)<-cnames
SL<-SL[SL$x1+SL$x2+SL$x3+SL$x4+SL$x5+SL$x6+SL$x7+SL$x8+SL$x9==1, ]
dimS<-dim(SL)
rows<-dimS[1]
rnames<-paste(1:rows)
rownames(SL)<-rnames
SL<-data.frame(SL)
return(SL)
                  }
  while (fac==10)   {
SL<-expand.grid(seq(0,1,by=(1/lev)), seq(0,1,by=(1/lev)),
      seq(0,1,by=(1/lev)), seq(0,1,by=(1/lev)),
      seq(0,1,by=(1/lev)), seq(0,1,by=(1/lev)),
      seq(0,1,by=(1/lev)), seq(0,1,by=(1/lev)),
      seq(0,1,by=(1/lev)), seq(0,1,by=(1/lev)))

names(SL)<-cnames
SL<-SL[SL$x1+SL$x2+SL$x3+SL$x4+SL$x5+SL$x6+SL$x7+SL$x8+SL$x9+SL$x10==1, ]
dimS<-dim(SL)
rows<-dimS[1]
rnames<-paste(1:rows)
rownames(SL)<-rnames
SL<-data.frame(SL)
return(SL)
                  }

  while (fac==11)   {
SL<-expand.grid(seq(0,1,by=(1/lev)), seq(0,1,by=(1/lev)),
      seq(0,1,by=(1/lev)), seq(0,1,by=(1/lev)),
      seq(0,1,by=(1/lev)), seq(0,1,by=(1/lev)),
      seq(0,1,by=(1/lev)), seq(0,1,by=(1/lev)),
      seq(0,1,by=(1/lev)), seq(0,1,by=(1/lev)), 
      seq(0,1,by=(1/lev)))

names(SL)<-cnames
SL<-SL[SL$x1+SL$x2+SL$x3+SL$x4+SL$x5+SL$x6+SL$x7+SL$x8+SL$x9+SL$x10+SL$x11==1, ]
dimS<-dim(SL)
rows<-dimS[1]
rnames<-paste(1:rows)
rownames(SL)<-rnames
SL<-data.frame(SL)
return(SL)
                  }
  while (fac==12)   {
SL<-expand.grid(seq(0,1,by=(1/lev)), seq(0,1,by=(1/lev)),
      seq(0,1,by=(1/lev)), seq(0,1,by=(1/lev)),
      seq(0,1,by=(1/lev)), seq(0,1,by=(1/lev)),
      seq(0,1,by=(1/lev)), seq(0,1,by=(1/lev)),
      seq(0,1,by=(1/lev)), seq(0,1,by=(1/lev)), 
      seq(0,1,by=(1/lev)), seq(0,1,by=(1/lev)))

names(SL)<-cnames
SL<-SL[SL$x1+SL$x2+SL$x3+SL$x4+SL$x5+SL$x6+SL$x7+SL$x8+SL$x9+SL$x10+SL$x11+SL$x12==1, ]
dimS<-dim(SL)
rows<-dimS[1]
rnames<-paste(1:rows)
rownames(SL)<-rnames
SL<-data.frame(SL)
return(SL)
                  }
                        }

