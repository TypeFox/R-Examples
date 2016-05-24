## ----B,echo=FALSE--------------------------------------------------------
B=10000

## ----loadpkg-------------------------------------------------------------
library(survSNP)

## ----ex1-----------------------------------------------------------------
res1<-sim.snp.expsurv.power(1.25, n=500, raf=0.1, erate=0.75, 
                            pilm=0.5, lm=1, B=0,
                            model="additive",test="additive",alpha=0.05)


## ----printex1run---------------------------------------------------------
res1[,c("n","GHR","erate","raf","B","alpha","pow0","pow","powB")]

## ----ex1b----------------------------------------------------------------
set.seed(123)
res1b<-sim.snp.expsurv.power(1.25, n=500, raf=0.1, erate=0.75, 
                             pilm=0.5, lm=1, 
                             exactvar=TRUE,B=B,
                             model="additive",test="additive",alpha=0.05)


## ----printex1brun--------------------------------------------------------
res1b[,c("n","GHR","erate","raf","B","alpha","pow0","pow","powB")]

## ----ex1setup------------------------------------------------------------
GHRs<-seq(1.05,1.5,by=0.05)
ns<-c(100,500,700)
rafs<-c(0.1,0.3,0.5)
erates<-c(0.5,0.7,0.9)

## ----ex2run--------------------------------------------------------------
res2<-survSNP.power.table(GHRs,ns,rafs,erates,
                         pilm=0.5,lm=1,
                         model="additive",test="additive",
                         alpha=0.05)

## ----printex2------------------------------------------------------------
res2[1:3,c("n","GHR","erate","raf","pow0","pow","powB")]

## ----res3----------------------------------------------------------------
cols<-c("n","GHR","erate","raf","pow0")
res3<-subset(res2,GHR==1.25&raf==0.3&n==500,select=cols)
res3

## ----tab,results='asis'--------------------------------------------------
print(xtable(res3,digits=c(0,0,1,1,1,3)),
      include.rownames=FALSE,floating=FALSE)

## ----ex1plot-------------------------------------------------------------
KEY=paste("q=",levels(factor(res2$raf)),sep="")
KEY<-list(lines=list(col=1:length(KEY),lty=1:length(KEY)),
          text=list(labels=paste("q=",levels(factor(res2$raf)),sep="")),
          column=3)
print(xyplot(pow0~GHR|factor(erate)*factor(n),group=factor(raf),
             data=res2,type="l",lty=KEY$lines$lty,col=KEY$lines$col,
             key=KEY,
             xlab="Genotype Hazard Ratio",ylab="Power"))

## ----ex2plot-------------------------------------------------------------
print(xyplot(pow0~GHR|factor(erate),group=factor(raf),
             data=subset(res2,n==ns[1]),
             type="l",lty=KEY$lines$lty,col=KEY$lines$col,
             key=KEY,
             xlab="Genotype Hazard Ratio",ylab="Power",
             sub=paste("n=",ns[1],", alpha=",round(unique(res2$alpha),2))))

## ----sessioninfo,results='asis',echo=FALSE-------------------------------
print(toLatex(sessionInfo()))

