"phxvpc.basicplot"<-function(
  model.name=NULL,
  xlab=NULL,
  ylab=NULL,
  xlab.cex=1.5,
  ylab.cex=1.5,
  x.cex=1.5,
  y.cex=1.5,
  main.title=NULL,
  main.cex=1.7,
  xlim=NULL,
  ylim=NULL,
  obs.pt=FALSE,
  obs.pch=16,
  logY=FALSE,
  Q.obs.line=TRUE,
  Q.pred.line=FALSE,

  CI.Q.pred="area",
  CI.Q.pred.area1="pink",
  CI.Q.pred.area2="grey",
  data.obs=data.obs,
  data.Q.obs=data.Q.obs,
  data.Q.pred=data.Q.pred,
  data.Q.CI.pred=data.Q.CI.pred,
  showkey=TRUE,
  key.position=c(0.5,-0.35),
  main.just=c(0.5,2)
  ){
  
  if(is.null(xlim)) xlim=c(NA,NA)                  
  if(is.null(ylim)) ylim=c(NA,NA)
  if(is.null(xlab)) {xlab=list(label="IVAR", cex=xlab.cex)} else {xlab=list(label=xlab, cex=xlab.cex)}
  if(is.null(ylab)) {ylab=list(label="DV", cex=ylab.cex, just=c(0.5,1) )} else {ylab=list(label=ylab, cex=ylab.cex, just=c(0.5,1))}
  if(logY)
    {  DV.all<-c(data.obs$DV, data.Q.obs$DV0, data.Q.pred$DV, data.Q.CI.pred$DV)
       DV.all<-DV.all[!is.na(DV.all)]
      if(any(DV.all<0)) {cat ("Error: negative DV exist. Data can not be log transformed...")
                         setwd("..")
                         return(NULL)
      }else{yscale=list(log = 10, cex=y.cex)}       
    }else{yscale=list(cex=y.cex)}
  
  xscale=list(cex=x.cex)
  if(!is.null(main.title)) {main=list(label=main.title, cex=main.cex, just=main.just)} else{
   main=list(label=paste(model.name, "Visual Predictive Check", sep=" "), cex=main.cex, just=main.just)}
  
  
  listk=list()
  if (obs.pt==T) listk=c(listk, list(text=list("obs"), points=list(pch=obs.pch,col="blue")))
  if(Q.obs.line==T)  { lty.Q.obs<-vector("numeric", length(unique(data.Q.obs$QI)))
                       for(i in 1:length(unique(data.Q.obs$QI))){
                         lty.Q.obs[i]<-ifelse((i==1|i==length(unique(data.Q.obs$QI))),5,1)
                       }
                       listk=c(listk, list( text=list(paste(as.character(unique(data.Q.obs$QI)), "obs", sep=" ")), 
                                            lines=list( lty=lty.Q.obs, col=c("red" ) )))
                       
  }
  if(Q.pred.line==T){
    lty.Q.pred<-vector("numeric", length(unique(data.Q.pred$Q)))
    for(i in 1:length(unique(data.Q.pred$Q))){
       lty.Q.pred<-ifelse((i==1|i==length(unique(data.Q.pred$Q))),5,1)   
    }
    listk=c(listk, list( text=list(paste(as.character(unique(data.Q.pred$Q)), "pred", sep="")), 
                         lines=list( lty=lty.Q.obs, col=c("grey10" ) )))
    
  }
  if(CI.Q.pred=="lines")  listk=c(listk, list( text=list(paste(paste(as.character(unique(data.Q.CI.pred$QE)[-1]), collapse=","), "CI.pred", sep=" ")), 
                                              lines=list( lty=3, col=c("black") )))
 
  if(CI.Q.pred=="area") { n=length(unique(data.Q.CI.pred$QE))
                          listk=c(listk, list( text=list(paste(as.character(unique(data.Q.CI.pred$QE)[2]), 
                                                         "-",as.character(unique(data.Q.CI.pred$QE)[n]),"CI.pred", sep=" ")), 
                                               rectangles = list(height=c(0.8, 0.8), col=c(CI.Q.pred.area1, CI.Q.pred.area2 ) ) ))
                                  
  }
  
  listk=c(listk, list(corner=key.position, rep=F ) )#space="bottom" )
  
  if(showkey==F) listk=list(text=list(c(" ")))

  pl<-xyplot(DV~IVAR,data.obs, xlim=xlim, ylim=ylim, xlab=xlab, ylab=ylab, scales=list(x=xscale, y=yscale), main=main,
       panel=function(x,y){
         if (obs.pt==T) panel.xyplot(x,y,pch=obs.pch, col="blue")
                         
         
         #show lines of observed quantiles
         if(Q.obs.line==T){
         lty.Q.obs<-vector("numeric", length(unique(data.Q.obs$QI)))
         col.Q.obs<-rep("red", length(unique(data.Q.obs$QI)))
         for(i in 1:length(unique(data.Q.obs$QI))){
           dat<-data.Q.obs[which(data.Q.obs$QI==unique(data.Q.obs$QI)[i]),]              
           panel.xyplot(dat$IVAR, dat$DV0, 
                        type="l", lty=ifelse((i==1|i==length(unique(data.Q.obs$QI))),5,1),lwd=2,col=col.Q.obs[i])             
           lty.Q.obs[i]<-ifelse((i==1|i==length(unique(data.Q.obs$QI))),5,1)
         }
         }
         
         #show lines of predicted quantiles
         if(Q.pred.line==T){
           lty.Q.pred<-vector("numeric", length(unique(data.Q.pred$Q)))
           col.Q.pred<-rep("grey10", length(unique(data.Q.pred$Q)))
           for(i in 1:length(unique(data.Q.pred$Q))){
             dat<-data.Q.pred[which(data.Q.pred$Q==unique(data.Q.pred$Q)[i]),]              
             panel.xyplot(dat$IVAR, dat$DV, 
                          type="l", lty=ifelse((i==1|i==length(unique(data.Q.pred$Q))),5,1),lwd=2,col=col.Q.pred[i])
             lty.Q.pred<-ifelse((i==1|i==length(unique(data.Q.pred$Q))),5,1)   
           }
         }
         
         #use lines to show confidence intervials of predicted quantiles 
         if(!is.null(CI.Q.pred))
         {
           if(CI.Q.pred=="lines")
           {
            for(i in 1:length(unique(data.Q.CI.pred$QI)))
            {dat0<-data.Q.CI.pred[which(data.Q.CI.pred$QI==unique(data.Q.CI.pred$QI)[i]),]
             for (j in 1:length(unique(dat0$QE)))
             {
               dat<-dat0[which(dat0$QE==unique(dat0$QE)[j]),]
               panel.xyplot(dat$IVAR, dat$DV,
                            type="l", lty=3, lwd=1, col="black")
             }            
            }
           }else        #use shaded areas to show confidence intervials of predicted quantiles  
           {if(CI.Q.pred=="area")
           {for(i in 1:length(unique(data.Q.CI.pred$QI)))
           {color=ifelse((i==1|i==length(unique(data.Q.CI.pred$QI))),CI.Q.pred.area2,CI.Q.pred.area1)
            col.fill=adjustcolor(color, alpha.f=0.3)
            dat<-data.Q.CI.pred[which(data.Q.CI.pred$QI==unique(data.Q.CI.pred$QI)[i]),]
            n=length(unique(dat$QE))
            
            CI.lower<-unique(data.Q.CI.pred$QE)[2]
            CI.upper<-unique(data.Q.CI.pred$QE)[n]            
            CI.lower.DV<-dat[which(dat$QE==CI.lower),]$DV
            CI.lower.IVAR<-dat[which(dat$QE==CI.lower),]$IVAR
            CI.upper.DV<-dat[which(dat$QE==CI.upper),]$DV
            CI.upper.IVAR<-dat[which(dat$QE==CI.upper),]$IVAR
            if(logY){
              panel.polygon(c(CI.lower.IVAR, rev(CI.upper.IVAR)), c(log(CI.lower.DV,10), rev(log(CI.upper.DV,10) )), col = col.fill, border = FALSE)
            }else{panel.polygon(c(CI.lower.IVAR, rev(CI.upper.IVAR)), c(CI.lower.DV, rev(CI.upper.DV) ), col = col.fill, border = FALSE)
            }
            }                      
           }         
           }
         }
         
         
       }
       ,key=listk
       #       ,key=list(text=list(c("DV1", "DV2")),lines=list(pch=c(19,0), col=c("blue", "red")), columns=2)
)
return(pl) 
  
}
