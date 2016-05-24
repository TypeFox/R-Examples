#' @export
"phxvpc.plot"<-function(
  vpcpath="",
  xlab=NULL,
  ylab=NULL,
  xlab.cex=1.3,
  ylab.cex=1.3,
  x.cex=1.3,
  y.cex=1.3,
  main.title=NULL,
  main.cex=1.3,
  xlim=NULL,
  ylim=NULL,
  obs.pt=FALSE,
  obs.pch=16,
  logY=FALSE,
  Q.obs.line=TRUE,
  Q.pred.line=TRUE,
  CI.Q.pred="area",
  CI.Q.pred.area1="pink",
  CI.Q.pred.area2="grey",
  ppp=4,
  legend=T,
  result.path=NULL,
  pred.corr=FALSE,
  data.obs=NULL,
  data.Q.obs=NULL,
  data.Q.pred=NULL,
  data.Q.CI.pred=NULL
){
  setwd(vpcpath)
  setwd('..')
  rootwd=getwd()
  if(is.null(result.path)) result.path=paste(rootwd,"Results",sep="/")
  if(!file.exists(result.path)) dir.create(result.path)
  
  setwd(vpcpath)
    
  if(is.null(data.obs)) data.obs<-read.csv("predcheck0.csv")
  if(is.null(data.Q.obs)) data.Q.obs<-read.csv("predcheck2.csv")  
  if(is.null(data.Q.pred)) data.Q.pred<-read.csv("predcheck1.csv")
  if(is.null(data.Q.CI.pred)) data.Q.CI.pred<-read.csv("predcheck2.csv")
  
  if(is.null(xlab)) {xlab=readLines("ivar.txt")
   if(xlab=="t") xlab="Time" 
  }
  model.name=basename(getwd())
  if(is.null(main.title)) main.title=paste(model.name, "Visual Predictive Check", sep=" ")
  
  #get stratification info
  if(file.exists("strata.txt")) 
  {strata.names=read.table("strata.txt")
   strata.lable=paste(strata.names[,1], collapse=",")
   
   strata<-function(data){
     if(length(strata.names[,1])==1) data$strata=data$Strat1
     if(length(strata.names[,1])==2) data$strata=paste(data$Strat1, data$Strat2, sep=",")
     if(length(strata.names[,1])==3) data$strata=paste(data$Strat1, data$Strat2, data$Strat3, sep=",")
     return(data)
   }
   
   data.obs<-strata(data.obs)
   data.Q.obs<-strata(data.Q.obs) 
   data.Q.pred<-strata(data.Q.pred)
   data.Q.CI.pred<-strata(data.Q.CI.pred)
   
   strata.level<-unique(c(data.obs$strata,data.Q.obs$strata,data.Q.pred$strata, data.Q.CI.pred$strata))
  }else{
   strata.level=0
  }
  
  pdfname=paste(result.path, paste(model.name,"pdf",sep="."), sep="/")
  showkey=TRUE
  if(length(strata.level)<=1) 
    { if(legend==F) showkey=F
      pdf(file=pdfname)
      VPC<-phxvpc.basicplot(   model.name=model.name,
      xlab=xlab,
      ylab=ylab,
      xlab.cex=xlab.cex,
      ylab.cex=ylab.cex,
      x.cex=x.cex,
      y.cex=y.cex,
      main.title=main.title,
      main.cex=main.cex,
      xlim=xlim,
      ylim=ylim,
      obs.pt=obs.pt,
      obs.pch=obs.pch,
      logY=logY,
      Q.obs.line=Q.obs.line,
      Q.pred.line=Q.pred.line,                                          
      CI.Q.pred=CI.Q.pred,
      CI.Q.pred.area1=CI.Q.pred.area1,
      CI.Q.pred.area2=CI.Q.pred.area2,
      data.obs=data.obs,
      data.Q.obs=data.Q.obs,
      data.Q.pred=data.Q.pred,
      data.Q.CI.pred=data.Q.CI.pred,
      showkey=showkey
       )
     print(VPC, position=c(0.05,0.1,0.9,0.95))
     dev.off()
    }
  
  if(length(strata.level)>1)
  { VPC<-list()
    pdf(file=pdfname)
    par(mar=c(0, 0, 0, 0))
    par(xaxs='i', yaxs='i' )
    plot.new()
    text(0.5, 0.95, main.title, offset=0, cex=main.cex)
    
    if(ppp!=1&ppp!=4) 
    {ppp=1
     write("Warning: ppp is an integer number either 1 or 4. ppp set to 1... ")
    } 
    
    if(ppp==1) {
      if(legend==F) showkey=F
      for(i in 1:length(strata.level)){
      main.title.i=paste(strata.lable, strata.level[i], sep=":")
      data.obs.i<-data.obs[which(data.obs$strata==strata.level[i]),]
      data.Q.obs.i<-data.Q.obs[which(data.Q.obs$strata==strata.level[i]),]
      data.Q.pred.i<-data.Q.pred[which(data.Q.pred$strata==strata.level[i]),]
      data.Q.CI.pred.i<-data.Q.CI.pred[which(data.Q.CI.pred$strata==strata.level[i]),]
      VPC[[i]]<-phxvpc.basicplot(model.name=model.name,
                                 xlab=xlab,
                                 ylab=ylab,
                                 xlab.cex=xlab.cex,
                                 ylab.cex=ylab.cex,
                                 x.cex=x.cex,
                                 y.cex=y.cex,
                                 main.title=main.title.i,
                                 main.cex=main.cex,
                                 xlim=xlim,
                                 ylim=ylim,
                                 obs.pt=obs.pt,
                                 obs.pch=obs.pch,
                                 logY=logY,
                                 Q.obs.line=Q.obs.line,
                                 Q.pred.line=Q.pred.line,                                          
                                 CI.Q.pred=CI.Q.pred,
                                 CI.Q.pred.area1=CI.Q.pred.area1,
                                 CI.Q.pred.area2=CI.Q.pred.area2,
                                 data.obs=data.obs.i,
                                 data.Q.obs=data.Q.obs.i,
                                 data.Q.pred=data.Q.pred.i,
                                 data.Q.CI.pred=data.Q.CI.pred.i,
                                 showkey=showkey
      )
      print(VPC[[i]], position=c(0.05,0.1,0.9,0.95), newpage=ifelse(i==1,F,T))
      }
      dev.off()
    }
    
      
    ### adjust cex size for 4 panels per page
    if(ppp==4){
      xlab.cex=xlab.cex*0.67
      ylab.cex=ylab.cex*0.67
      x.cex=x.cex*0.67
      y.cex=y.cex*0.67
      main.cex=main.cex*0.67     
    

    for(i in 1:length(strata.level))
    { if(i%%4==1) {newpage=ifelse(i==1,F,T) 
                   if((i==length(strata.level)|i==length(strata.level)-1)&(legend==T)) {
                     showkey=T
                     position=c(0.05,0.55,0.5,0.925)
                     key.position=c(0.1,-1.3)
                     main.just=c(0.5,1.5)}else{
                     position=c(0.05,0.55,0.5,0.95)
                     key.position=c(0,-0.35)
                     main.just=c(0.5,3)
                     showkey=F   
                   }
                    }
      if(i%%4==2) {position=c(0.5,0.55,0.95,0.95)
                   key.position=c(0,-0.35)
                   main.just=c(0.5,3)
                   newpage=F
                   showkey=F}
      if(i%%4==3) {if(legend==T) {
                    position=c(0.05,0.15,0.5,0.525)
                    newpage=F
                    showkey=T
                    key.position=c(0.1,-1.3)
                    main.just=c(0.5,1.5)
                    }else{
                      position=c(0.05,0.15,0.5,0.55)
                      newpage=F
                      showkey=F
                      key.position=c(0,-0.35)
                      main.just=c(0.5,3)                      
                    }
                  }
      if(i%%4==0) {position=c(0.5,0.15,0.95,0.55)
                   newpage=F
                   showkey=F
                   key.position=c(0,-0.35)
                   main.just=c(0.5,3)}
      
      
      main.title.i=paste(strata.lable, strata.level[i], sep=":")
      data.obs.i<-data.obs[which(data.obs$strata==strata.level[i]),]
      data.Q.obs.i<-data.Q.obs[which(data.Q.obs$strata==strata.level[i]),]
      data.Q.pred.i<-data.Q.pred[which(data.Q.pred$strata==strata.level[i]),]
      data.Q.CI.pred.i<-data.Q.CI.pred[which(data.Q.CI.pred$strata==strata.level[i]),]
      VPC[[i]]<-phxvpc.basicplot(model.name=model.name,
                               xlab=xlab,
                               ylab=ylab,
                               xlab.cex=xlab.cex,
                               ylab.cex=ylab.cex,
                               x.cex=x.cex,
                               y.cex=y.cex,
                               main.title=main.title.i,
                               main.cex=main.cex,
                               xlim=xlim,
                               ylim=ylim,
                               obs.pt=obs.pt,
                               obs.pch=obs.pch,
                               logY=logY,
                               Q.obs.line=Q.obs.line,
                               Q.pred.line=Q.pred.line,                                          
                               CI.Q.pred=CI.Q.pred,
                               CI.Q.pred.area1=CI.Q.pred.area1,
                               CI.Q.pred.area2=CI.Q.pred.area2,
                               data.obs=data.obs.i,
                               data.Q.obs=data.Q.obs.i,
                               data.Q.pred=data.Q.pred.i,
                               data.Q.CI.pred=data.Q.CI.pred.i,
                               showkey=showkey,
                               key.position=key.position,
                               main.just=main.just
        )
    
    
    print(VPC[[i]], position=position, newpage=newpage)      
    }
    dev.off()
    
    }
      
      
       
   
  }
  
  # open created file
  if (file.exists(pdfname)) {
    if (Sys.info()['sysname'] == 'Windows') { shell.exec(pdfname) }  # windows
    else if (Sys.info()['sysname'] == 'Darwin') { system(paste ("open ",pdfname, sep="")) } # mac
    else { system(paste ("xdg-open ",pdfname, sep="")) } # linux
  }
  
  setwd(rootwd)
}
  
  
 