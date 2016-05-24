# several auxiliary plotting functions:
library(lattice)
library(plotrix)

plotH<-function(X,Y,INT,DISP=0.3,horiz=FALSE,...)
{
 if(horiz)
 { segments( (X-INT),Y,(X+INT),Y,...)
   segments( (X-INT),Y-DISP,(X-INT),Y+DISP,...)
   segments( (X+INT),Y-DISP,(X+INT),Y+DISP,...)
 }
 else { segments( X,Y-INT,X,Y+INT,...)
        segments( X-DISP,Y-INT,X+DISP,Y-INT,...)
        segments( X-DISP,Y+INT,X+DISP,Y+INT,...)
      }
}
#----------------------------------------
mhistogram=function(x,Levels=NULL,PDF="",main="",xlab="Values",ylab="Frequency",xaxt="n",yaxt="n",xaxis=FALSE,size=NULL,density=FALSE,...)
{
 if(!is.null(size)){ width=size[1];height=size[2];} 
 else {width=5;height=5;}
 if(PDF!="") { file=paste(PDF,".pdf",sep="");pdf(file,width=width,height=height); par(mar=c(3.0,4.0,0.2,0.2)) }
 if(is.factor(x)) 
 { L=length(levels(x))
   if(is.null(Levels)) Levels=levels(x)
   #print(Levels)
   hist(as.numeric(x),breaks=seq(0.5,(L+0.5),by=1),xlim=c(0.5,L+0.5),xlab=xlab,ylab=ylab,main=main,xaxt="n",yaxt=yaxt,...) 
   #hist(as.numeric(x),breaks=seq(0.5,(L+0.5),by=1),xlim=c(1,L),xlab=xlab,ylab=ylab,main=main,yaxt=yaxt,...) 
   #if(density) lines(density(as.numeric(x)))
   if(xaxt!="n") axis(1,1:L,Levels)
 }
 else { hist(x,xlab=xlab,ylab="Frequency",main=main,xaxt="n",...)
        #if(density) lines(density(as.numeric(x)))
      }
 if(PDF!="") dev.off()
}


# todo: add xlab and ylab options?
# mgraph
# graph options:
# "REC"
# "ROC"
# "LIFT"
# "IMP"
# "VEC"
# "RSC" - regression scatter
# "REP" - regression error plot
# "REG" - regression plot
# "DLC" - distance line comparation, for comparing errors
# y - list of minings or just 1 mining or =Y (test set)
# y can be just one vector
# x - predicions if y is not mining(s)
# intbar - if confidence interval bars should be plotted
# TC - target class or if "IMP" and TC=1 no above axis is used
# leg - NULL - not used, -1 for ROC and LIFT class level name is used, else can be list(pos=c(x,y),leg=text) or simple legend with text
# xval - domain of x in "REG", if xval!=-1
# Grid - if > 1, then a grey grid is plotted
mgraph=function(y,x=NULL,graph,leg=NULL,xval=-1,PDF="",PTS=-1,size=c(5,5),sort=TRUE,ranges=NULL,data=NULL,digits=NULL,TC=-1,
                intbar=TRUE,lty=1,col="black",main="",metric="MAE",baseline=FALSE,Grid=0,axis=NULL,cex=1.0)
{
#y=L;graph="REC";Grid=10;leg=c("mlpe","mr");main="REC curve";xval=0.1
#x=NULL;PDF="";PTS=-1;size=c(5,5);sort=TRUE;ranges=NULL;data=NULL;digits=NULL;TC=-1;intbar=TRUE;lty=1;col="black";metric="MAE";baseline=FALSE;axis=NULL;cex=1.0
#x=NULL;xval=-1;PDF="";PTS=-1;size=c(5,5);sort=TRUE;ranges=NULL;data=NULL;digits=NULL;TC=-1;intbar=TRUE;lty=1;col="black";main="";metric="MAE";baseline=FALSE;axis=NULL;cex=1.0
 if(!is.null(x)) { N=1;runs=1; 
                   TT=vector("list",1); PP=TT;
                   TT[[1]]=y;PP[[1]]=x; 
                   y=vector("list",1);
                   y[[1]]=list(test=TT,pred=PP,runs=runs)
                   MC=vector("list",1);
                 }
 else if(graph=="REC" || graph=="ROC" || graph=="LIFT" || graph=="RSC" || graph=="REP" || graph=="DLC" || graph=="REG") 
                                  { 
                                    if(is.null(y$runs)) # list with more than 1 mining
                                    {
                                     N=length(y)
                                     runs=y[[1]]$runs
                                     MC=vector("list",N)
                                    }
                                    else{ N=1; runs=y$runs; # transform 1 mining into list of mining 
                                          y2=vector("list",1);y2[[1]]=y;y=y2
                                          MC=vector("list",1);
                                        }
                                  }
 else if(graph=="IMP") { N=1
                         if(is.list(y)) {runs=y$runs}
                         else {sen=matrix(nrow=1,ncol=length(y));sen[1,]=y;y=list(sen=sen);runs=1}
                       }
 else if(graph=="VEC") { runs=y$runs
                         if(xval==-1) { N=length(y$sresponses)
                                        MC=vector("list",N)
                                      }
                         else {N=1;MC=vector("list",1)}
                         xlab=""; ylab="";
                       }
 #------------------------------------------------------------
 if(PTS==-1 && (graph=="ROC"||graph=="REC"||graph=="LIFT")) PTS=11
 else if(PTS==-1 && graph=="VEC") PTS=6

 if(graph=="ROC") {xlab="FPR";ylab="TPR";
                   if(xval==-1) xval=1;
                   ymax=1;DISP=(0.1/mean(PTS)) 
                   NX=NCOL(y[[1]]$pred[[1]])
                   if(NX>2 && TC==-1) TC=NX # last class?
                  }
 if(graph=="LIFT") {xlab="Sample size";ylab="Responses";
                   if(xval==-1) xval=1;
                   ymax=1;DISP=(0.1/mean(PTS)) 
                   NX=NCOL(y[[1]]$pred[[1]])
                   if(NX>2 && TC==-1) TC=NX # last class?
                  }
 else if(graph=="REC") { ymax=1;xlab="Absolute deviation";
                         if(xval==-1) xval= max(y[[1]]$test[[1]]) # print("Error: need to define xval\n"); return (-1);}
                         ylab="Accuracy";DISP=0.1*(xval/mean(PTS))
                       }

 if(graph=="VEC")
 {
   NX=NCOL(y$pred[[1]])
   if(NX>2 && TC==-1) TC=1
   C=length(y$sresponses); 
   if(is.null(ranges))
   {
    ranges=matrix(ncol=2,nrow=C) 
    if(xval==-1) {
                  for(i in 1:C)
                  { 
                   if(!is.null(y$sresponses[[i]]))
                   {
                    if(is.factor(y$sresponses[[i]]$x)) ranges[i,]=c(1,y$sresponses[[i]]$l)
                    else ranges[i,]=c(min(y$sresponses[[i]]$x),max(y$sresponses[[i]]$x))
                   }
                  }
                 }
     else if(is.factor(y$sresponses[[xval]]$x)) ranges[xval,]=c(1,y$sresponses[[xval]]$l)
     else ranges[xval,]=c(min(y$sresponses[[xval]]$x),max(y$sresponses[[xval]]$x))
   }
   if(is.null(leg)) { leg=vector(length=C)
                      if(xval==-1){ for(i in 1:C) if(!is.null(y$sresponses[[i]])) leg[i]=y$sresponses[[i]]$n }
                      else leg[xval]=y$sresponses[[xval]]$n
                    }
 }

 if(graph=="LIFT" || graph=="ROC") 
  { if(TC==-1 & NX<=2) TC=2 
    if(is.factor(y[[1]]$test[[1]])) POSITIVE=levels(y[[1]]$test[[1]])[TC]
  }

 if(graph=="LIFT" || graph=="ROC" || graph=="REC" || graph=="IMP" || graph=="VEC") 
 {
  for(j in 1:N)
  {
   #cat("j:",j,"\n")
   if(graph=="ROC" || graph=="REC" || graph=="LIFT")
   {
    C=vector("list",runs)
    for(i in 1:runs) 
     { 
     #cat("j:",j,"i:",i,"\n")
     if(graph=="ROC") { # warning? 
                       C[[i]]=ROCcurve(y[[j]]$test[[i]],y[[j]]$pred[[i]],TC=TC)$roc
                     }
     else if(graph=="LIFT") C[[i]]=twoclassLift(y[[j]]$test[[i]],y[[j]]$pred[[i]][,TC],Positive=POSITIVE,type=3,STEPS=(PTS-1))
     else if(graph=="REC") C[[i]]=RECcurve(y[[j]]$test[[i]]-y[[j]]$pred[[i]])
     }
   }
   #cat("j:",j,"\n")
   if(graph=="ROC"||graph=="REC"||graph=="LIFT") 
   { 
    if(length(PTS)>1) MC[[j]]=vaveraging(PTS[[j]],C,0,xval)
    else MC[[j]]=vaveraging(PTS,C,0,xval)
   }
   else if(graph=="IMP") MC=meanint(y$sen)
   else if(graph=="VEC") 
   { 
    if(xval!=-1 && !is.null(y$sresponses[[xval]]))
      {
       AUX=resp_to_list(y$sresponses[[xval]],TC)
       MC[[j]]=vaveraging(y$sresponses[[xval]]$l,AUX,ranges[xval,1],ranges[xval,2])
      }
    else if(!is.null(y$sresponses[[j]]))
      { 
       AUX=resp_to_list(y$sresponses[[j]],TC) 
       MC[[j]]=vaveraging(y$sresponses[[j]]$l,AUX,ranges[j,1],ranges[j,2])
      }
   }
  }
 } # if(graph== ...

 if(PDF!="" && graph!="VEC") 
  { file=paste(PDF,".pdf",sep="");pdf(file,width=size[1],height=size[2]); 
    # mar=c(bottom, left, top, right)
    if(graph=="ROC"||graph=="REC"||graph=="LIFT") 
      { if(main!="") par(mar=c(4.0,4.0,1.0,0.7)) 
        else par(mar=c(4.0,4.0,0.3,0.7)) 
      }
    else if(graph=="IMP") 
      { if(main!="") par(mar=c(2.0,1.0,2.0,0.0)) 
        else par(mar=c(2.0,1.0,2.0,0.0))
      }
    else if(graph=="RSC") par(mar=c(2.0,2.0,2.0,2.0))
    else if(graph=="DLC") par(mar=c(2.0,1.0,0.0,1.0))
   }

 if(graph=="ROC"||graph=="REC"||graph=="LIFT") # multi mining graphs
 {
  plot(MC[[1]][,1],MC[[1]][,2],type="n",xlab=xlab,ylab=ylab,ylim=c(0,ymax),lwd=2,main=main,panel.first=grid(Grid,Grid))
  if(length(col)==1) col=rep(col,N)
  if(length(lty)==1) lty=1:N
  for(j in 1:N)
  {
   lines(MC[[j]][,1],MC[[j]][,2],lty=lty[j],col=col[j],lwd=2)
   if(runs>1 && intbar) plotH(MC[[j]][,1],MC[[j]][,2],MC[[j]][,3],DISP=DISP,horiz=FALSE)
  }
  if(baseline && graph=="REC") baseline=FALSE
  if(baseline) { BSL=matrix(ncol=2,nrow=2);BSL[1,]=c(0,0);BSL[2,]=c(1,1)
                 lines(BSL,col="gray40",lwd=2)
               }
  if(is.numeric(leg) && leg[1]==-1 && (graph=="ROC"||graph=="LIFT")) leg=POSITIVE
  if(!is.null(leg)) { if(is.list(leg)) { cleg=leg$pos; leg=leg$leg; } 
                      else{
                           if(graph=="ROC"||graph=="LIFT") cleg=c(0.65,0.45)
                           else if(graph=="REC") cleg=c((0.6*xval),0.3) 
                          }
                      if(baseline) { leg=c(leg,"baseline");col=c(col,"gray40");lty=c(lty,1);}
                      if(length(cleg)==2) legend(cleg[1],cleg[2],leg,lty=lty,col=col,lwd=2) else legend(cleg,leg,lty=lty,col=col,lwd=2)
                    }
 }
 else if(graph=="VEC")
 {
  # b,l,t,r
  for(i in 1:N)
  {
   if(xval==-1) XATT=i
   else XATT=xval
   #cat("XATT:",XATT,"\n")
   if(!is.null(y$sresponses[[XATT]])) 
   {
    if(PDF!="") { file=paste(PDF,"-",XATT,".pdf",sep=""); 
                  pdf(file,width=size[1],height=size[2]); 
                  par(mar=c(2.0,2.0,2.0,2.0)) # b,l,t,r
                }
    if(runs>1 && intbar) intbar=MC[[i]][,3] else intbar=NULL
    resp=vector("list",length=1)
    x=y$sresponses[[XATT]]$x
    if(is.factor(x)) x=levels(x) else if(!is.character(x)) x=MC[[i]][,1]
    resp[[1]]=list(n="",l=nrow(MC[[i]]),x=x,y=MC[[i]][,2])
    IMC=list(imp=c(1),val=c(1),sresponses=resp)
    vecplot(I=IMC,graph="VEC",leg=leg,xval=1,sort=sort,data=data[,c(xval,xval)],digits=digits,TC=TC,intbar=intbar,lty=lty,col=col,datacol="gray90",main=main,Grid=Grid,xlab=xlab,ylab=ylab)
    if(PDF!="") dev.off()
   }
  }
 }
 else if(graph=="IMP") # single mining graphs
 {
  if(is.null(axis)) axis=c(1,3)
  # zero influence is not show in the graph? => yes!
  #MMC<<-MC;MC=MMC;leg=NULL;sort=TRUE;Grid=10;main="";
  if(xval==-1) xval=0.01
  XLAB=xval
  CI=which(MC$mean!=0) # & MC$int!=0) 
  #CI=which(MC$mean!=0 & MC$int!=0) 
  #CI=1:11

  # deal with output label if not included in label: clean this ugly code latter...
  LL=length(leg);LMC=length(MC$mean)

  if(LL>0 && LL<LMC)
   {
    if(!is.null(y$attributes)) yind=y$attributes[[1]][length(y$attributes[[1]])]
    else { # some heuristics that will not work in every situation:
          if(MC$mean[LMC]==0 && MC$int[LMC]==0) yind=LMC
          else if(MC$mean[1]==0 && MC$int[1]==0) yind=1
          else yind=CI[1]
        }
    leg=append(leg,yind,(yind-1))
   }
  # ---
  MM=MC$mean[CI]
  II=MC$int[CI]
  leg=leg[CI]
  IM=which.max(MM); xval=1.07*(max(MM)+II[IM])
  LMM=length(MM)
  if(sort==TRUE)
   {
    #cat("sort:",sort,"\n")
    S=sort(MM,index.return =TRUE) 
    MM=S$x
    II=II[S$ix]
    if(!is.null(leg)) leg=leg[S$ix]
   }
  else
   {
    INVI=LMM:1   
    MM=MM[INVI]
    II=II[INVI]
    if(!is.null(leg)) leg=leg[INVI]
   }

  #cat("MM:",MM,"\n")
  #cat("leg:",leg,"\n")
  if(col=="black") col="white"
  Yx=barplot(MM,names.arg=NULL,axes=FALSE,col=col,xlim=c(0,xval),horiz=TRUE,main=main,panel.first=grid(Grid,Grid)) #,ylab="inputs",xlab="relative importance")
  Yx=Yx[,1]
  #D1=0.7; D2=2.1;
  #Yx=seq(D1,LMM+D2,length=LMM)
  #cat("mean:",MM,"\n")
  #cat("text:",leg,"\n")
  DISP=0.3
  if(runs>1 && intbar) plotH(MM,Yx,II,DISP=DISP,horiz=TRUE)
  Xt=rep(XLAB,LMM)
  if(!is.null(leg)){ YDisp=0.1
                     text(Xt,Yx+YDisp,leg,pos=4,cex=cex)
                   }
  #cat("TC",TC,"\n")
  if(sum(axis==1)==1) axis(1,lty = 1,lwd=1,cex.axis=cex)
  if(sum(axis==3)==1) axis(3,lty = 1,lwd=1,cex.axis=cex)
 } 
 else if(graph=="RSC")
 {
  MIN=Inf;MAX=-Inf
  for(j in 1:N)
   for(i in 1:runs) 
    {
     MIN=min(MIN,y[[j]]$test[[i]]);MAX=max(MAX,y[[j]]$test[[i]])
     MIN=min(MIN,y[[j]]$pred[[i]]);MAX=max(MAX,y[[j]]$pred[[i]])
    }
  xlab="Observed"; ylab="Predicted"
  plot(y$test[[1]],y$pred[[1]],type="n",xlim=c(MIN,MAX),ylim=c(MIN,MAX),main=main,xlab=xlab,ylab=ylab,panel.first=grid(Grid,Grid))
  if(length(col)==1) col=rep(col,N*runs)
  for(j in 1:N) for(i in 1:runs) 
  {
     points(y[[j]]$test[[i]],y[[j]]$pred[[i]],pch=19,cex=0.8,col=col[j+i-1])
  }
  abline(0,1)
  if(!is.null(leg)) { if(is.list(leg)) { cleg=leg$pos; leg=leg$leg;}
                      else{
                           if(length(leg)==1 && leg==-1) leg=c("Predictions")
                           cleg="bottomright"
                          }
                     if(length(cleg)==2) legend(cleg[1],cleg[2],leg,lty=lty,col=col,lwd=2) else legend(cleg,leg,pch=19,lty=0,col=col,lwd=2)
                    }
 }
 else if(graph=="DLC")
 {
  if(runs>1 && intbar) AUX=matrix(ncol=2,nrow=N)
  else AUX=matrix(ncol=1,nrow=N)
  for(j in 1:N)
    {
     mi=meanint(mmetric(y[[j]],metric=metric,val=xval))
     AUX[j,1]=mi$mean; if(runs>1 && intbar) AUX[j,2]=mi$int
    }
  if(worst(metric)==Inf) {llab="better";rlab="worse";} 
  else {llab="worse";rlab="better";} 
  dlplot(AUX,main=main,labels=leg,llab=llab,rlab=rlab,size=size)
 }
 else if(graph=="REG") # to do: sort data?
 {
  MIN=Inf;MAX=-Inf
  if(length(xval)>1) XX=xval 
  else if(length(xval)==1 && xval>0) { LL=length(y[[1]]$test[[1]]); XX=xval:LL; }
  else { LL=length(y[[1]]$test[[1]]); XX=1:LL; }
  for(j in 1:N)
   for(i in 1:runs) 
    {
     MIN=min(MIN,y[[j]]$test[[i]][XX]);MAX=max(MAX,y[[j]]$test[[i]][XX])
     MIN=min(MIN,y[[j]]$pred[[i]][XX]);MAX=max(MAX,y[[j]]$pred[[i]][XX])
    }
  xlab="Examples"; ylab="Values"
  #plot(XX,y$test[[1]][XX],type="n",xlim=c(XX[1],XX[length(XX)]),ylim=c(MIN,MAX),main=main,xlab=xlab,ylab=ylab,panel.first=grid(Grid,Grid))
  plot(XX,y[[1]]$test[[1]][XX],type="n",xlim=c(XX[1],XX[length(XX)]),ylim=c(MIN,MAX),main=main,xlab=xlab,ylab=ylab,panel.first=grid(Grid,Grid))
  if(length(lty)==1) lty=rep(lty,2)
  if(length(col)==1) col=rep(col,runs*N)
  for(j in 1:N) for(i in 1:runs) 
  {
     lines(XX,y[[j]]$test[[i]][XX],type="b",pch=1,cex=0.8,lty=lty[1],lwd=2,col=col[1])
     lines(XX,y[[j]]$pred[[i]][XX],type="b",pch=19,cex=0.6,lty=lty[2],lwd=2,col=col[i+j])
  }
  if(!is.null(leg)) { if(is.list(leg)) { cleg=leg$pos; leg=leg$leg;}
                      else{
                           if(length(leg)==1 && leg==-1) leg=c("Target","Predictions")
                           cleg="bottomright"
                          }
                     if(length(cleg)==2) legend(cleg[1],cleg[2],leg,lty=lty,col=col,lwd=2) else legend(cleg,leg,lty=lty,col=col,lwd=2)
                    }
 }
 else if(graph=="REP") # to do: sort data?
 {
  AUX=NULL
  for(j in 1:N)
  for(i in 1:runs) 
   {
     AUX=c(AUX, (y[[j]]$test[[i]]-y[[j]]$pred[[i]]) )
   }
  if(sort) AUX=sort(AUX) 
  MIN=min(AUX);MAX=max(AUX)
  xlab="Examples"; ylab="Residuals"
  plot(AUX,type="n",ylim=c(MIN,MAX),main=leg,xlab=xlab,ylab=ylab,panel.first=grid(Grid,Grid))

  for(j in 1:N) for(i in 1:runs) points(AUX,pch=19,cex=0.8)
  abline(0,0)
 }
 if(PDF!="" && graph!="VEC") dev.off()
}

# distance line comparation plot (for comparing errors)
dlplot=function(m,main="",yval=0.5,pch=15,cex=2,labels=NULL,llab="worse",rlab="better",size=c(4,2))
{
 DISP=0.01; YDISP=DISP*0.7; 
 NR=nrow(m); NC=ncol(m);
 m2=matrix(ncol=2,nrow=NR)
 XMIN=Inf;XMAX=-Inf
 for(i in 1:NR) { m2[i,1]=m[i,1]; m2[i,2]=yval;
                  if(XMIN>m[i,1]) XMIN=m[i,1]
                  if(XMAX<m[i,1]) XMAX=m[i,1]
                  if(NC>1 && XMIN>(m[i,1]-m[i,2])) XMIN=m[i,1]-m[i,2]
                  if(NC>1 && XMAX<(m[i,1]+m[i,2])) XMAX=m[i,1]+m[i,2]
                }
 XDISP=0.01*max(abs(c(XMIN,XMAX)))
 XLIM=c(XMIN-XDISP,XMAX+XDISP)
 #op=par(din=size)
 plot(m2,type="n",xlim=XLIM,ylim=c((yval-2*DISP),(yval+2*DISP)),xlab="",ylab="",axes=FALSE,frame.plot =FALSE)
 axis(1)
 #axis(2)
 m3=matrix(ncol=2,nrow=2)
 m3[,1]=XLIM; m3[,2]=yval 
 lines(m3)
 points(m2,pch=pch,cex=cex)
 if(NC>1)
 {
  plotH(m2[,1],m2[,2],m[,2],DISP=DISP*0.2,horiz=TRUE)
 }
 if(!is.null(labels))
 {
  text(m2[,1],m2[,2]+YDISP,labels=labels)
 }
 if(llab!="") text(XMIN-XDISP,yval-1.8*DISP,llab,pos=4) 
 if(rlab!="") text(XMAX+XDISP,yval-1.8*DISP,rlab,pos=2) 
 if(main!="") text(mean(c(XMIN,XMAX)),yval+1.5*DISP,cex=1.25,labels=bquote(bold(.(main))))
 #par(op)
}

# plot.rpart uniform=TRUE or FALSE, branch = 0 to 1, compresss=TRUE or FALSE, pretty = 0 (no label abbreviation) or NULL (label abbreviation)
# drawtree extra: nodeinfo=TRUE
modelplot=function(FIT,main="",drawtree=FALSE,FANCY=FALSE,data=NULL,...)
{
 if(FIT@model=="dt")
 {
   #T=try(library(maptree),silent=TRUE)
   #if(class(T)=="try-error" || drawtree==FALSE) 
   # {
      plot(FIT@object,main=main,uniform=TRUE,branch=0,compress=TRUE)
      if(FANCY) text(FIT@object,pretty=0,xpd=TRUE,fancy=TRUE,fwidth=0.2,fheight=0.2)
      else text(FIT@object,pretty=0,xpd=TRUE)
   # }
   #else
   #{
   #  if(main!="") 
   #  {
   #     plot(0,type="n",main=main,frame.plot=FALSE,xaxt="n",yaxt="n",ylab="",xlab="")
   #     draw.tree(FIT@object,new=FALSE,...)
   #  }
   #  else draw.tree(FIT@object,...)
   #}
 }
 else if(FIT@model=="svm" && FIT@task=="prob" && length(FIT@levels)==2 ) # svm only supports binary classification
 {
  # in development...
  svm=FIT@object$svm
  print(svm)
 }
}

# to do: check plot VEC3 and VECC, improve these, VEC - pch in points?
# importanceplot
# xval - attribute, for "VEC", xval can be a vector of attributes to included (multi-vec plot!)
# sort="increasing", TRUE, "decreasing", "increasing2",  "decreasing2", or FALSE 
# TC= target class
# screeb=list(z = 40, x = -60)
# screen=list(z = 20, x = -70, y = 0)
# screen=list(z=0,x=-90,y= 0)  # X, Z
# screen=list(z = 0, x = -90, y = -20)) # from X to Y
# screen=list(z = 0, x = -90, y = -45)) # from X to Y
# screen=list(z = 0, x = -90, y = -20)) # from X to Y
# screen=list(z=0,x=-90,y=-90) # Y, Z
# screen=list(z=-90,x=-90,y=-90)) #  
# screen=list(z=0,x=0,y=0)) #  X vs Y, plane Z

# showlevels = FALSE, TRUE, c(TRUE,FALSE,TRUE), ...
# col = "black"
# col = "grayrange" - for "VEC3", "VECC"
# col = "white"

#mgraph=function(LM,X=NULL,graph,leg=NULL,xval=-1,PDF="",PTS=-1,size=c(5,5),sort=TRUE,ranges=NULL,data=NULL,digits=NULL,TC=-1,
#                intbar=TRUE,lty=1,col="black",main="",metric="MAD",baseline=FALSE,Grid=FALSE)

# TC used also for metric: min, average or max
vecplot=function(I,graph="VEC",leg=NULL,xval=1,sort=FALSE,data=NULL,digits=c(1,1),TC=1,
                 intbar=NULL,lty=1,pch=19,col=NULL,datacol=NULL,main="",main2="",Grid=0,
                 xlab="",ylab="",zlab="",levels=NULL,levels2=NULL,showlevels=FALSE,screen=list(z=40,x=-60),zoom=1,cex=1.0)
{
#if(length(TC)>1) { LTC=length(TC);
#cat("LTC:",LTC,"\n")
#                   par(mar=c(2.0,2.0,0.1,0.1))
#                   par(mfrow=c(1,LTC))
#                   for(t in 1:LTC) vecplot(I=I,graph=graph,leg=leg,xval=xval,sort=sort,data=data,digits=digits,TC=t,
#                                           intbar=intbar,lty=lty,pch=pch,col=col,datacol=datacol,main=main,main2=main2,Grid=Grid,
#                                           xlab=xlab,ylab=ylab,zlab=zlab,levels=levels,levels2=levels2,showlevels=showlevels,screen=screen,zoom=zoom)
#                 }
#else
#{ 

#I=i_1_2;graph="VEC3";xlab=n1[1];ylab=n1[2];TC=TC;zoom=1.1;screen=list(z=0,x=0,y=0)
#leg=NULL;xval=1;sort=FALSE;data=NULL;digits=c(1,1)
#intbar=NULL;lty=1;pch=19;col=NULL;datacol=NULL;main="";main2="";Grid=0
#zlab="";levels=NULL;levels2=NULL;showlevels=FALSE;cex=1.0
 sort2=FALSE;decreasing=FALSE;
 if(sort=="increasing2"){sort=FALSE;sort2=TRUE}
 else if(sort=="decreasing2"){decreasing=TRUE;sort=FALSE;sort2=TRUE}
 else if(sort=="decreasing"){decreasing=TRUE;sort=TRUE} else if(sort=="increasing"){sort=TRUE}
 LS=length(showlevels)
 if(LS==1 && showlevels==FALSE) showlevels=c(FALSE,FALSE,FALSE)
 else if(LS==1 && showlevels==TRUE) showlevels=c(TRUE,TRUE,FALSE)
 else if(LS==2) showlevels=c(showlevels,FALSE)

 if(!is.list(screen)) screen=switch(screen,x=list(z=0,x=-90,y=0),X=list(x=-75),y=list(z=0,x=-90,y=-90),Y=list(z=10,x=-90,y=-90),z=list(z=0,x=0,y=0),xy=list(z=10,x=-90,y=-45))
#cat("--- ord:",sort,"ord2:",sort2,"dec:",decreasing,"showlevels:",showlevels,"col:",col,"\n")

 if(!is.null(I$method))
 { if( I$method=="GSA") { 
                        if(length(I$interactions)==1)
                        {
                         if(length(xval)>1) xval=xval[1]
                        }
                        else
                        {
                         I=aggregate_imp(I,AT=xval,measure=I$measure,Aggregation=I$agg,method=I$method,L=I$Llevels)
                         if(length(xval)>1) xval=xval[1]
                        }
                       }
   else if( (graph!="VECB" && graph!="VEC") && (I$method=="DSA"||I$method=="MSA"))
                       {
                        # xval with 2 attributes:
                        I=aggregate_imp(I,AT=xval,measure=I$measure,Aggregation=I$agg,method=I$method,L=I$Llevels)
                        if(length(xval)>1) xval=xval[1]
                       }
 }

 if(graph=="VEC")
 {
  if(is.null(col)) col="black"
  if(length(xval)==1)
  { x=I$sresponses[[xval]]$x;
    y=I$sresponses[[xval]]$y;
    if(is.factor(y)){ylevels=levels(y);y=as.numeric(y);yaxt="n"}else{ylevels=NULL;yaxt="s";}
    if(class(y)=="matrix"||class(y)=="data.frame"){ if(I$agg==3) y=y[2,] else y=y[,TC]}
    n=I$sresponses[[xval]]$n
    if(!is.null(digits)) M=paste(n,", var: ",round(I$val[xval],digits=digits[1]),sep="") else M=paste(n,", var: ",I$val[xval],sep="")
    if(!is.null(intbar)) YLIM=c(min(y-intbar),max(y+intbar)) else YLIM=range(y)

    if( (class(x)!="character" && !is.factor(x)) && !is.null(data))
    { DIF=diff(x)/2;LX=length(x)
      BREAKS=rep(NA,LX+1)
      BREAKS[1]=x[1]-DIF[1]
      for(i in 1:(LX-1)) BREAKS[i+1]=x[i]+DIF[i]
      BREAKS[LX+1]=x[LX]+DIF[LX-1]
      XLIM=range(BREAKS)
    }
    else if(class(x)!="character" && !is.factor(x)){XLIM=range(x);DIF=(x[2]-x[1])/2;}
    else if(class(x)=="character" || is.factor(x)){ 
                                if(sort) {SY=sort.int(y,decreasing=decreasing,index.return=TRUE);levels=x[SY$ix];y=y[SY$ix]
                                          if(!is.null(intbar)) intbar=intbar[SY$ix]
                                         } 
                                if(!is.null(levels)) x=factor(x,levels=levels) else x=factor(x)
                                L=length(x);levels=levels(x)
                                if(!is.null(data)) data[,xval]=factor(data[,xval],levels=levels)
                                  } 
    if(!is.null(data))
    {
      XLIM=range(BREAKS)
      BREAKS="Sturges";
      #cat("BREAKS:",BREAKS,"XLIM:",XLIM,"xval:",xval,"\n")
      if(is.factor(data[,xval])) mhistogram(data[,xval],xaxt="n",yaxt="n",xlab="",ylab="",main="",col=datacol,labels=TRUE)
      else hist(data[,xval],xlim=XLIM,main="",xaxt="n",yaxt="n",xlab="",ylab="",col=datacol,breaks=BREAKS)
      #axis(1) #
      axis(4)
      par(new=TRUE)#,mar=c(0,0,0,0))
    }

    if(is.factor(x)) { 
                                DIF=0.5;XLIM=c(1-DIF,L+DIF)
                                if(!is.null(intbar)) YLIM=c(min(y-intbar),max(y+intbar)) else YLIM=range(y)
                                plot(1,y[1],xlim=XLIM,ylim=YLIM,type="n",main=main,xlab=xlab,ylab=ylab,xaxt="n",yaxt=yaxt,panel.first=grid(Grid,Grid))
                                x=seq(1,L,length.out=L);
                                DIF2=0.01
                                if(showlevels[1]) text(x+DIF/2,y+DIF2,levels,cex=cex)
                                segments(x-DIF,y,x+DIF,y,lwd=2)
                                axis(1,x,levels,cex.axis=cex)
                              }
    else plot(x,y,xlim=XLIM,ylim=YLIM,type="b",yaxt=yaxt,main=main,xlab=xlab,ylab=ylab,lwd=2,pch=pch,panel.first=grid(Grid,Grid))
    if(yaxt=="n") axis(2,unique(y),ylevels,cex.axis=cex)
    if(!is.null(intbar)) { 
    			  DISP=0.18*min(DIF)
                          plotH(x,y,intbar,DISP=DISP,col=col,horiz=FALSE) 
                         }
  }
  else  # else length xval
  {
   if(length(lty)==1) lty=1:length(xval)
   if(length(pch)==1) pch=rep(pch,length(xval))
   if(length(col)==1) col=rep(col,length(xval))
   YMAX=-Inf;YMIN=Inf
   for(i in xval){ y=I$sresponses[[i]]$y;
                   if(is.factor(y)){ylevels=levels(y);y=as.numeric(y);yaxt="n"} else {ylevels=NULL;yaxt="s";}
                   if(class(y)=="matrix"||class(y)=="data.frame") y=y[,TC]
                   YMIN=min(YMIN,y);YMAX=max(YMAX,y)
                 }
   #plot(1,1,xlim=c(0,1),ylim=c(YMIN,YMAX),main=main,type="n",yaxt=yaxt,xlab=paste(xlab," (scaled)",sep=""),ylab=ylab,panel.first=grid(Grid,Grid))
   plot(1,1,xlim=c(0,1),ylim=c(YMIN,YMAX),main=main,type="n",yaxt=yaxt,xlab=paste(xlab," (scaled)",sep=""),xaxt="n",ylab=ylab,panel.first=grid(Grid,Grid))
   if(yaxt=="n") axis(2,unique(y),ylevels)
   yi=1
   for(i in xval)
   {
     x=I$sresponses[[i]]$x;y=I$sresponses[[i]]$y;
     if(is.factor(y)){ylevels=levels(y);y=as.numeric(y);yaxt="n"} else {ylevels=NULL;yaxt="s";}
     if(class(y)=="matrix"||class(y)=="data.frame") y=y[,TC]
     n=I$sresponses[[i]]$n;LY=length(y)
     if(class(x)=="character" || is.factor(x)) 
                               { 
                                if(sort) {SY=sort.int(y,decreasing=decreasing,index.return=TRUE);levels=x[SY$ix];y=y[SY$ix]
                                          if(!is.null(intbar)) intbar=intbar[SY$ix]
                                         } 
                                #if(!is.null(levels)) x=factor(x,levels=levels) else x=factor(x)
                                #L=length(x);levels=levels(x)
                                #if(!is.null(data)) data[,xval]=factor(data[,xval],levels=levels)
                                DIF=1/LY
                                x=seq(0,1-DIF,length.out=LY);
                                DIF2=0.01
                                #levels[3]="dec"
                                if(showlevels[1]) text(x+DIF/2,y+DIF2,levels,cex=cex)
                                segments(x,y,x+DIF,y,lwd=2,lty=lty[yi])
                               }
     else { x=seq(0,1,length.out=LY)
            lines(x,y,type="b",lwd=2,lty=lty[yi],pch=pch[yi])
          }
     yi=yi+1
   }
   if(!is.null(leg)) { if(is.list(leg)) { cleg=leg$pos; leg=leg$leg; } 
                       else cleg="topright"
                       if(length(cleg)==2) legend(cleg[1],cleg[2],leg,lty=lty,col=col,lwd=2,cex=cex) else legend(cleg,leg,lty=lty,col=col,lwd=2,cex=cex)
                     }
  } # endif length xval
 } # endif "VEC"
 else if(graph=="VECB")
 {
   if(is.null(col)) col="black"
   L=I$sresponses[[xval]]$l
   NC=ncol(I$sresponses[[xval]]$y)
   if(I$nclasses==0) SEQ=1:L else SEQ=seq(TC,L*NC,NC)
   B=I$sresponses[[xval]]$yy[,SEQ]
   xlabels=I$sresponses[[xval]]$x
   if(is.numeric(xlabels)) xlabels=round(xlabels,digits=digits[1])
   rmboxplot(B,MEAN=TRUE,LINE=TRUE,MIN=FALSE,MAX=FALSE,ALL=FALSE,BOXPLOT=TRUE,col=col,xlabels=xlabels,main=main,cex=cex,xlab=xlab,Grid=Grid)

 }
 else if(graph=="VECC" || graph=="VEC3" || graph=="VECC2")
 {
  # encoding here!

 if(is.null(col)) col="grayrange"
  #I=I5;xval=1;xlab="";ylab="";digits=c(1,1);zlab="Y";main="";main2="";TC=1;levels=NULL;decreasing=FALSE;sort=TRUE;sort2=FALSE;graph="VECC";levels2=NULL;levels2=c("D","E","C","B","A"))
  #I=I7;xval=4;xlab="";ylab="";digits=c(1,1);zlab="Y";main="";main2="";TC=1;levels=NULL;sort=TRUE;graph="VECC";levels2=NULL;levels=NULL;decreasing=FALSE
  #screen=list(z = -10,x=-60)
 if(is.null(I$xy)) # normal
  {
   if(length(xlab)==1 && xlab=="") xlab=I$sresponses[[xval]]$n[1]
   if(length(ylab)==1 && ylab=="") ylab=I$sresponses[[xval]]$n[2]
   x=I$sresponses[[xval]]$x
   y=I$sresponses[[xval]]$y

  if(is.factor(y)){ylevels=levels(y);y=as.numeric(y)} else {ylevels=NULL}
  if(class(y)=="matrix"||class(y)=="data.frame") y=y[,TC]
  n=I$sresponses[[xval]]$n
  L=I$sresponses[[xval]]$l[1]
  LY=length(y)
  L2=LY/I$sresponses[[xval]]$l
  NR=L;NC=L2;TRANSPOSE=FALSE;
  #cat("L:",L,"LY:",LY,"L2:",L2,"NC:",NC,"NR:",NR,"\n")
  if(x[1,1]==x[2,1]) { 
                       I=sort.int(as.numeric(x[,2]),index.return=TRUE)$ix
                       x=x[I,];y=y[I]
                     }
  # problems with "unif", change this code below: problem with "quantile"? check later...
  LAB1=unique(x[,1])
  #if(is.factor(LAB1)) LAB1=1:L
  LAB2=unique(x[,2])
#print(digits)
  # end of change this code below -----
  if(!is.factor(LAB1) && !is.null(digits))LAB1=round(LAB1,digits=digits[1])
  else if(!is.null(levels)) { 
                              TRANSPOSE=TRUE;NR=L2;NC=L;
                              levels(LAB1)=levels; 
                              I=dforder(x,1,levels);
                              x=x[I,]; y=y[I]; 
                            }
  else if(sort)
  { 
    TRANSPOSE=TRUE;NR=L2;NC=L;
    I=dforder(x,1,y=y,decreasing=decreasing);
    levels=I$levels;I=I$I;levels(LAB1)=levels;
    x=x[I,]; y=y[I]; 
  }
  if(!is.factor(LAB2) && !is.null(digits))LAB2=round(LAB2,digits=digits[2])
  else if(!is.null(levels2)) { levels(LAB2)=levels2;
                               I=dforder(x,2,levels2);
                               x=x[I,]; y=y[I]; 
                             } 
  else if(sort2)
  { 
    I=dforder(x,2,y=y,decreasing=decreasing);
    levels2=I$levels;I=I$I;levels(LAB2)=levels2;
    x=x[I,]; y=y[I]; 
  }
#print(LAB2)
#YY<<-y
  A=matrix(y,nrow=NR,ncol=NC)
  if(TRANSPOSE) A=t(A) # check this better...

  A=fenlarge(A,LAB1,LAB2)
#AB<<-A
 }
 else{
      TRANSPOSE=FALSE
#cat(" >> ELSE << \n")
      if(length(xlab)==1 && xlab=="") xlab=I$xlab
      if(length(xlab)==1 && ylab=="") ylab=I$ylab
      n=xlab;LAB1=I$x1;LAB2=I$x2;L=I$L;NLY=I$NLY;
      ylevels=NULL;
      #if(!is.factor(LAB1) && !is.null(digits))LAB1=round(LAB1,digits=digits[1])
      #if(!is.factor(LAB2) && !is.null(digits))LAB2=round(LAB2,digits=digits[2])
      if(is.numeric(LAB1) && !is.null(digits))LAB1=round(LAB1,digits=digits[1])
      if(is.numeric(LAB2) && !is.null(digits))LAB2=round(LAB2,digits=digits[2])
      if(showlevels[1]){xlab=paste(xlab,": ",round(I$value[TC],digits=2),sep="") }
      if(showlevels[2]){ylab=paste(ylab,": ",round(I$value2[TC],digits=2),sep="")}

      ini=(TC-1)*L[2]+1;end=ini+L[2]-1;
#print(ini:end)
      NR=L[1];NC=L[2]
      y=I$Y[,ini:end]
#III<<-I
if(FALSE){
      if(sort)
      { 
       TRANSPOSE=TRUE;I=dforder(x,1,y=y,decreasing=decreasing);
       LAB1=LAB1[I]; 
       #x=x[I,]; 
       y=y[I]; 
      }
      else if(sort2)
      { 
       I=dforder(x,2,y=y,decreasing=decreasing);
       LAB2=LAB2[I]
       #x=x[I,]; 
       y=y[I]; 
      }
}

      A=matrix(y,nrow=NR,ncol=NC)
#AB<<-A
      if(TRANSPOSE) A=t(A) # check this better...
      A=fenlarge(A,LAB1,LAB2)
#cat("TC:",TC,"L2:",L[2],"ini:",ini,"end:",end,"\n")
      #x=seq(0,1,length.out=L[1]);ax=x;
      #y=seq(0,1,length.out=L[2]);ay=y;
      #A=list(z=I$Y[,ini:end],x=x,y=y,ax=ax,ay=ay)
#AC<<-A
     }

#AA<<-A
  #A=t(matrix(y,nrow=L,ncol=L2))
  if(graph=="VEC3")
  {
   if(showlevels[1]==2){ scales=list(arrows=FALSE,cex=0.4,col="black",font=3,tck=1)} 
   else {
   scales=list()
   if(showlevels[1]){ if(is.factor(LAB1)) xlab=paste(xlab,": ",paste(LAB1,collapse=","),sep="")
                      else xlab=paste(xlab,": ",LAB1[1],",...,",LAB1[length(LAB1)],sep="")
                    }
   if(showlevels[2]){ if(is.factor(LAB2)) ylab=paste(ylab,": ",paste(LAB2,collapse=","),sep="")
                      else ylab=paste(ylab,": ",LAB2[1],",...,",LAB2[length(LAB2)],sep="")
                    }
   if(showlevels[3]) zlab=paste(zlab,": ",min(y),",...,",max(y),sep="")
        }
   if(length(col)==1 && col=="grayrange") col=gray(70:1/70)
#cat("YL:",ylevels,"\n")
#AAA<<-A
#COL<<-col
#YLEVELS<<-ylevels
#A=AAA;col=COL;ylevels=YLEVELS;main="";xlab="x1";zlab="y";ylab="x2";scales=list();screen=list(z=40,x=-60);zoom=1;
#   if(is.null(ylevels)) 
   if(length(zlab)==1 && zlab=="") zlab=list(label="y",cex=cex)
wireframe(A$z,row.values=A$x,scales=scales,column.values=A$y,xlab=xlab,ylab=ylab,zlab=zlab,screen=screen,main=main,drape=TRUE,col.regions=col,zoom=zoom)
#   else wireframe(A$z,row.values=A$x,column.values=A$y,scales=scales,xlab=xlab,ylab=ylab,zlab=zlab,screen=screen,main=main,drape=TRUE,zoom=zoom,
#                  col.regions=col,key.axes=axis(4,1:length(ylevels),ylevels))

  }
  else if(graph=="VECC")
  {
   #GC=gray
   #RY=range(y,finite=TRUE)
   #LUY=length(unique(y))
#   GC=gray.colors(LUY,start=0.9,end=0,gamma=1)
   #GC=gray.colors(1000,start=0.9,end=0,gamma=1)
   flevels = pretty(range(A$z,finite=TRUE),20);LFL=length(flevels)
   if(length(col)==1 && col=="grayrange") col=gray(LFL:1/LFL)
#cat("YL:",ylevels,"\n")
#print(A$ax)
#print(LAB1)
#print(A$ay)
#print(LAB2)
#AY<<-A$ay;LAB2<<-LAB2
   if(is.null(ylevels)) filled.contour(A$x,A$y,A$z,col=col,plot.title=title(main=main,xlab=xlab,ylab=ylab),key.title=title(main=main2),plot.axes = {axis(1,at=A$ax,labels=LAB1);axis(2, at=A$ay,labels=LAB2);})
   else filled.contour(A$x,A$y,A$z,col=col,plot.title=title(main=main,xlab=xlab,ylab=ylab),key.axes=axis(4,1:length(ylevels),ylevels),key.title=title(main=main2),plot.axes = {axis(1,at=A$ax,labels=LAB1);axis(2, at=A$ay,labels=LAB2);})
  }
 }
#} # special else
}

# ---- internal functions, do not use these:
dforder=function(d,col,levels=NULL,y=NULL,decreasing=FALSE) 
{ I=NULL; 
  if(!is.null(levels)) {for(i in 1:length(levels)) I=c(I,which(d[,col]==levels[i])); return(I)}
  else if(!is.null(y)) {
                        levels=levels(d[1,col]);L=length(levels);ymean=rep(NA,L)
                        for(i in 1:length(levels)) ymean[i]=mean(y[which(d[,col]==levels[i])])
                        levels=levels[sort.int(ymean,decreasing=decreasing,index.return=TRUE)$ix]
                        return(list(I=dforder(d,col,levels=levels,y=NULL),levels=levels)) 
                       } 
  else return(1:NROW(d))
}

enlarge=function(A,X=0,Y=0,LAB1=NULL,LAB2=NULL)
{
 if(X>0)
 {
  NR=NROW(A);NC=NCOL(A);k=1
  for(i in 1:NR)
  { New=matrix(data=NA,nrow=X,ncol=NC)
    if(i<NR) A=rbind(A[1:k,],New,A[(k+1):NROW(A),]) else A=rbind(A[1:k,],New)
    k=k+X+1
  }
 }
 if(Y>0)
 {
  NR=NROW(A);k=1
  for(i in 1:NC)
  { New=matrix(data=NA,nrow=NR,ncol=Y)
    if(i<NC) A=cbind(A[,1:k],New,A[,(k+1):NCOL(A)]) else A=cbind(A[,1:k],New)
    k=k+Y+1
  }
 }
 return (A)
}

fenlarge=function(A,LAB1=NULL,LAB2=NULL,DIF=0.001)
{
 if(is.factor(LAB1)) LX=levels(LAB1) else if(is.character(LAB1)) LX=LAB1 else LX=""
 if(length(LX)>1 || LX!="")
 {
  NL=length(LX);
  sx=sort(c(seq(0,1,length.out=(NL+1))-DIF,seq(0,1,length.out=(NL+1))+DIF))
  I=which(sx>DIF & sx<(1-DIF));x=c(0,sx[I],1);AX=1/length(x);ax=seq(AX,1,by=2*AX) 
  NR=NROW(A);NC=NCOL(A);k=1
  for(i in 1:NR)
  { A=rbind(A[1:k,],A[k:NROW(A),]) 
    k=k+1+1
  }
 }
 else {x=seq(0,1,length.out=NROW(A));ax=x;}

 if(is.factor(LAB2)) LY=levels(LAB2) else if(is.character(LAB2)) LY=LAB2 else LY=""
 if(length(LY)>1 || LY!="")
 {
  NL=length(LY);
  sx=sort(c(seq(0,1,length.out=(NL+1))-DIF,seq(0,1,length.out=(NL+1))+DIF))
  I=which(sx>DIF & sx<(1-DIF));y=c(0,sx[I],1);AY=1/length(y);ay=seq(AY,1,by=2*AY) 
  NR=NROW(A);NC=NCOL(A);k=1
  for(i in 1:NC)
  { A=cbind(A[,1:k],A[,k:NCOL(A)]) 
    k=k+1+1
  }
 }
 else {y=seq(0,1,length.out=NCOL(A));ay=y;}
 return (list(x=x,y=y,z=A,ax=ax,ay=ay))
}

# --- ts plots ---:
tsplot=function(file,ts,test=12,xlab="Time",ylab="Values",main="",disp=2,disp2=0,left=0.4)
{
 pdf(file=paste(file,".pdf",sep=""),paper="special",width=5,height=5)
 #c(bottom, left, top, right)
 par(mar=c(2.0,2.0,1.8,left))
 plot(ts,type="l",lwd=2,xlab=xlab,ylab=ylab,main=main)
 L=length(ts)
 Y=min(ts) 
 DISP=0.95
 abline(v=(L-test))
 text(round((L-test)/2),disp2+round(Y*DISP),"training set")
 text(round(mean(c(L,L-test)))+disp,disp2+round(Y*DISP),"test set")
 cat(round(mean(c(L,L-test)))+disp,":",disp2+round(Y*DISP),"\n")
 dev.off()
}

tsacf=function(file="",ts,test=12,xlab="Lags",ylab="Autocorrelations",main="",disp=0.1,left=0.4,lag.max,K=12)
{
 if(file!="") pdf(file=paste(file,".pdf",sep=""),paper="special",width=5,height=5)
 #c(bottom, left, top, right)
 #if(file!="") par(mar=c(2.0,2.0,1.8,left))
 if(file!="") par(mar=c(2.0,1.8,0.1,0.1))
 L=length(ts)
 LTR=L-test
 M=acf(ts[1:LTR],lag.max=lag.max,plot=FALSE)

 plot(M,type="h",lwd=2,xlab=xlab,ylab=ylab,main=main)
 text(x=round(lag.max/2),y=0.9,cex=1.25,labels=bquote(bold(.(main)))) 
 
 text(x=K,y=(M$acf[K]+disp),cex=1,paste("K=",K,sep=""))
 if(file!="") dev.off()
}

# experimental
forplot=function(file="",ts,test=12,PRED,xlab="lead time",ylab="Values",disp=2,disp2=0,left=0.4,names=NULL,leg=-1,MIN=0.9,MAX=1.1,start=1,main="")
{
 if(file!="") pdf(file=paste(file,".pdf",sep=""),paper="special",width=5,height=5)
 #c(bottom, left, top, right)
 if(file!="") par(mar=c(2.0,2.0,1.6,left))
 L=length(ts)
 LI=L-test+start

 TS=ts[LI:L]
 ymin=min(TS,PRED)*MIN
 ymax=max(TS,PRED)*MAX
 #cat("min:",ymin,"max:",ymax,"\n")
 iend=start+test
 if(start>1) { plot(ts[LI:L],type="l",lwd=2,ylim=c(ymin,ymax),xlab=xlab,ylab=ylab,main=main,xaxt="n",panel.first = grid(10,10))
               #AT=axis(1)
               AT=round(seq(start,iend,length=5))
               cat("AT:",AT,"\n")
               axis(1,(AT-start+1),AT)
             }
 else plot(ts[LI:L],type="l",lwd=2,ylim=c(ymin,ymax),xlab=xlab,ylab=ylab,main=main,panel.first=grid(10,10))
 N=NCOL(PRED); LP=NROW(PRED)
 for(i in 1:N)
 {
  P=PRED[(start:LP),i]
  lines(P,type="l",lwd=2,lty=(i+1)) 
 } 
 if(!is.null(names)) 
 {
   if(leg[1]==-1) legend("topright",legend=names,lwd=2,lty=1:(N+1))
   else legend(leg,legend=names,lwd=2,lty=1:(N+1))
 }
 if(file!="") dev.off()
}

# experimental
rmboxplot=function(x,MEAN=TRUE,LINE=FALSE,MIN=FALSE,MAX=FALSE,ALL=FALSE,BOXPLOT=TRUE,col="black",xlabels="",main="",cex=1.0,sort=FALSE,xlab="",ylab="",Grid=0)
{
 NC=NCOL(x)
 if(NC>1){ 
           Means=vector(length=NC);for(i in 1:NC) Means[i]=mean(x[,i]);
           Min=vector(length=NC);for(i in 1:NC) Min[i]=min(x[,i]);
           Max=vector(length=NC);for(i in 1:NC) Max[i]=max(x[,i]);
           n=NC
         } 
 else {Means=mean(x);
       Min=min(x)
       MAX=max(x)
       n=length(x);}

 decreasing=FALSE
 if(sort=="decreasing"){decreasing=TRUE;sort=TRUE} else if(sort=="increasing"){sort=TRUE}
 if(sort){ SI=sort.int(Means,decreasing=decreasing,index.return=TRUE)
           xlabels=xlabels[SI$ix]
           Means=Means[SI$ix]
           Min=Min[SI$ix]
           Max=Max[SI$ix]
           x=x[,SI$ix]
         }
 DIM=0.3
 PCEX=2.0
 #plot(1:n,seq(min(Min),max(Max),length.out=n),ylab="",xlim=c(1-DIM,n+DIM),ylim=c(-0.02,1.025),type="n",xaxt="n")
 plot(1:n,seq(min(Min),max(Max),length.out=n),xlab=xlab,ylab=ylab,xlim=c(1-DIM,n+DIM),type="n",xaxt="n",yaxt="n",main=main,panel.first=grid(Grid,Grid))
#cat("main:",main,"\n")
 #if(main!="") text(mean(c(1,n)),1.025,main,cex=1.3)
 if(length(xlabels)==1 && xlabels=="") xlabels=1:n
 axis(1,1:n,xlabels,cex.axis=cex)
 axis(2,cex.axis=cex)
 #else plot(1:n,seq(min(Min),max(Max),length.out=n),type="n")
 if(ALL){   
         if(NC>1){ NR=NROW(x); for(i in 1:NR) lines(1:n,x[i,],col="gray") }
         else    { NR=length(x); for(i in 1:NR) lines(1:n,x,col="gray")   }
        }
 if(BOXPLOT) boxplot(x,add=TRUE,yaxt="n",xaxt="n",range=0)
 if(MEAN) points(1:n,Means,pch=18,cex=PCEX,col=col)
 if(LINE>1) LWD=LINE else LWD=1
 if(LINE) lines(1:n,Means,col=col,lwd=LWD)
 if(MIN)  lines(1:n,Min)
 if(MAX)  lines(1:n,Max)
}

# matrix of pair importances, x1, x2, minimum thresholds
# min - length 1 -> percentage x ...
cmatrixplot=function(m1,m2=NULL,threshold=c(0.1,0.1),L=7,cex=1.0)
{
 if(is.null(m2)) {m2=m1$m2;m1=m1$m1;}
 #m1=bm$m1;m2=bm$m2;threshold=c(0.2,0.29);L=L=datalevels(d1[,AT],L=7);cex=1.0
 if(length(threshold)==1) { threshold=rep(max(m1,m2)*threshold,2) }
 #m=m/max(m)
 NR=nrow(m1)
 m=matrix(0,nrow=NR,ncol=NR)
 for(i in 1:NR)
 for(j in 1:NR)
 {
   if(m1[i,j]>threshold[1] && m2[i,j]>threshold[2]) m[i,j]=m1[i,j]+m2[i,j]
 }
 #print(m)
 #c(bottom, left, top, right)
 par(mar=c(2.0,2.0,0.4,6))
# C2D=color2D.matplot(m,extremes=c("white","black"))
 C2D=color2D.matplot(m,extremes=c("white","black"),axes=FALSE)
# A=axTicks(1)
#A=1:11
 #axis(1,cex.axis=cex)
# cat("USR:",C2D$usr[2],"\n")
# AT=A
# if(A[length(A)]!=C2D$usr[2]) AT=AT+(C2D$usr[2]-A[length(A)])
# cat("ATx:",AT,"\n")
# cat("ATy:",AT[length(AT):1],"\n")
 DIFF=0.5
#AAT<<-AT
 axis(1,at=(1:NR)-DIFF,labels=(1:NR),cex.axis=cex)
 axis(2,at=(1:NR)-DIFF,labels=(NR:1),cex.axis=cex)
 col.labels=seq(min(m),max(m),length.out=L)
 col.labels=round(col.labels,digits=2)
 testcol=gray(seq(1,0,length.out=L))
 Step=NR*0.1
 #color.legend(NR+0.5*Step,(NC/2),NR+1.5*Step,NC,col.labels,testcol,gradient="y",align="rb")
 color.legend(NR+0.5*Step,0,NR+1.25*Step,NR,col.labels,testcol,gradient="y",align="rb",cex=cex)
 return(m)
}
