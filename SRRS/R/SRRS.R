SRRS <-
function() 
{ main=function(X,Y,nME,gamma,int,criterion,Xh,Yh)
  { nME=as.integer(nME)
    if(nME<=0|nME>ncol(X)) nME=ncol(X)
    XME=X[,1:nME]
    if(!Xh) colnames(XME)=paste("X",c(1:nME),sep="")
    X2fi=col2fi=c()
    if(int==T)
    { decompose=c(); 
      counter=nME+1
      for(i in 1:(nME-1))
      { for(j in (i+1):nME)
        { X2fi=cbind(X2fi,X[,i]*X[,j])
          col2fi=cbind(col2fi,c(i,j)) }}
      colnames(X2fi)=paste(names(XME)[col2fi[1,]],names(XME)[col2fi[2,]],sep=":")
      X=cbind(XME,X2fi)
      composeInt=rbind(c((nME+1):(nME+ncol(col2fi))),col2fi) }
    else { X=XME }
    pass=T; ono=0; br=0
    sf=include=rep(0,ncol(X))
    include[1:nME]=1
    Y0=Y
    rho=cor(Y,XME)
    sftemp=sftemp0=which(abs(rho)==max(abs(rho[1:nME])))[1]
    model=lm(data.matrix(Y)~X[,c(sftemp,which(sf!=0))])
    beta1=beta=summary(model)$coef[2,1]
    if(gamma<=0) gamma=round(abs(beta1)/10,4)
    assign("ggGamma",gamma,envir=tt$env)
    while(pass)
    { Y=Y-X[,sftemp]*beta
      if(sftemp0==sftemp) br=br+1
      sftemp0=sftemp
      if(sf[sftemp]==0)
      { ono=ono+1
        sf[sftemp]=1 }
      if(int==T & sftemp<=nME)
      { for(i in (nME+1):ncol(X))
        { include[i]=ifelse(col2fi[1,i-nME]==sftemp | col2fi[2,i-nME]==sftemp,1,include[i]) }}
      rho=cor(Y,X)*include
      sftemp=which(abs(rho)==max(abs(rho)))[1]
      model=lm(data.matrix(Y)~data.matrix(X[,c(sftemp,which(sf!=0))]))
      beta=summary(model)$coef[2,1]
      pass=abs(summary(model)$coef[2,1])>=gamma & ono<(nrow(X)-2) }
    sftlist=which(sf!=0)
    sftlistME=sftlist[sftlist<=nME]
    listlength=min(length(sftlist),ceiling(nrow(X)/3))
    listMElen=min(length(sftlistME),ceiling(nrow(X)/3))
    modelchoice=c()
    for(cn in 1:listMElen)
    { modelcols=combinations(length(sftlistME),cn,v=sftlistME)
      modelchoice=rbind(modelchoice, cbind(modelcols, matrix(0,nrow(modelcols),listMElen-ncol(modelcols)))) }
    if(int==T)
    { sftlistInt=sftlist[sftlist>nME]
      decompose=c()
      for(dc in 1:length(sftlistInt)) decompose=cbind(decompose,composeInt[,composeInt[1,]==sftlistInt[dc]])
      Intchoice=c()
      Intloopn=min(listlength-1,length(sftlistInt))
      for(cn in 1:Intloopn)
      { Intcols=combinations(length(sftlistInt),cn,v=sftlistInt)
        Intchoice=rbind(Intchoice,cbind(Intcols,matrix(0,nrow(Intcols),Intloopn-ncol(Intcols)))) }
      newmodelchoice=c()
      for(mm in 1:nrow(modelchoice))
      { locME=modelchoice[mm,which(modelchoice[mm,]!=0)]
        addflag=rep(0,nrow(Intchoice))
        for(ii in 1:nrow(Intchoice))
        { Intcolumn=Intchoice[ii,Intchoice[ii,]!=0]
          for(jj in 1:length(Intcolumn))
          { locInt=which(decompose[1,]==Intcolumn[jj])
            totprod=1
            for(pp in 1:length(locME)) totprod=totprod*prod(abs(locME[pp]-as.vector(decompose[2:3,locInt])))
            addflag[ii]=addflag[ii]+totprod }}
        addInt=rbind(rep(0,ncol(Intchoice)),Intchoice[addflag==0,])
        if(nrow(addInt)==2) combinem=c(modelchoice[mm,],addInt[-1,])
        if(nrow(addInt)!=2) combinem=cbind(t(t(rep(1,nrow(addInt)-1)))%*%t(modelchoice[mm,]),addInt[-1,])
        newmodelchoice=rbind(newmodelchoice,combinem) }
      verifylen=newmodelchoice
      verifylen[verifylen>0]=1
      vnum=apply(verifylen,1,sum)
      newmodelchoice=newmodelchoice[vnum<=listlength,]
      modelchoice=rbind(cbind(modelchoice,matrix(0,nrow(modelchoice),ncol(newmodelchoice)-ncol(modelchoice))),newmodelchoice) }
    critm=rep(0,nrow(modelchoice))
    for(cc in 1:length(critm))
    { mmobj=lm(data.matrix(Y0)~data.matrix(X[,modelchoice[cc,modelchoice[cc,]!=0]]))
      nnobj=length(mmobj$resid)
      kkobj=length(mmobj$coeff)-1
      RSSobj=sum(resid(mmobj)^2)
      if(criterion=="mAIC") critm[cc]=nnobj*log(RSSobj/nnobj)+2*kkobj^2
      if(criterion=="AIC") critm[cc]=nnobj*log(RSSobj/nnobj)+2*kkobj }
    Result=c()
    for(llen in 1:5)
    { tarloc=which(critm==min(critm))
      mincritm=critm[tarloc]
      Result=rbind(Result,c(modelchoice[tarloc,],mincritm))
      modelchoice=modelchoice[-tarloc,]
      critm=critm[-tarloc] }
    objvalue=Result[,ncol(Result)]
    modelcolno=Result[,-ncol(Result)]
    Modelname=c()
    for(mn in 1:5)
    { temp=colnames(X)[modelcolno[mn,modelcolno[mn,]!=0]]
      Modelname[mn]=temp[1]
      if(length(temp)>1)
      { for(wlen in 2:length(temp)) Modelname[mn]=paste(Modelname[mn],temp[wlen],sep="+") }}
    Resultmatrix=cbind(c(1:5),Modelname,round(objvalue,4))
    Resultmatrix=rbind(c("Rank","Model",criterion),Resultmatrix)
    colnames(Resultmatrix)=c("","","")
    return(Resultmatrix) }
  doSRRS <- function() 
  { assign("nMEkk",as.integer(tclvalue(nME)),envir=tt$env)
    assign("Gammakk",as.double(tclvalue(Gamma)),envir=tt$env)
    cbValue1kk=as.character(tclvalue(cbValue1))
    assign("cbValue1kkk",as.logical(as.integer(cbValue1kk)),envir=tt$env)
    cbhDMvalkk=as.character(tclvalue(cbhDMval))
    assign("cbhDMvalkkk",as.logical(as.integer(cbhDMvalkk)),envir=tt$env)
    cbhRVvalkk=as.character(tclvalue(cbhRVval))
    assign("cbhRVvalkkk",as.logical(as.integer(cbhRVvalkk)),envir=tt$env)
    criterionsChoice <- criterions[as.numeric(tclvalue(tcl(comboBox,"getvalue")))+1]
    if(criterionsChoice=="mAIC") assign("criterions",criterionsChoice,envir=tt$env)
    if(criterionsChoice=="AIC") assign("criterions",criterionsChoice,envir=tt$env)
    gXh=gYh=F
    gDM=get("DMF",envir=tt$env)
    if(get("cbhDMvalkkk",envir=tt$env)) 
    { gDM=get("DMT",envir=tt$env)
      gXh=as.logical(1,envir=tt$env) }
    gRV=get("RVF",envir=tt$env)
    if(get("cbhRVvalkkk",envir=tt$env)) 
    { gRV=get("RVT",envir=tt$env)
      gYh=as.logical(1,envir=tt$env) }
    gnME=get("nMEkk",envir=tt$env)
    gGamma=get("Gammakk",envir=tt$env)
    gInt=get("cbValue1kkk",envir=tt$env)
    gCriterion=get("criterions",envir=tt$env)
    Rematr=main(gDM, gRV, gnME, gGamma, gInt, gCriterion, gXh, gYh)
    Result=tktoplevel()
    tkwm.title(Result,"Result")
    tclRequire("Tktable")
    tclarray=tclArray()
    for(i in (0:5))
    { for(j in (0:2))
      { tclarray[[i,j]]=Rematr[i+1,j+1] }}
    tclarray[[6,0]]="gamma"
    tclarray[[6,1]]=get("ggGamma",envir=tt$env)
    table1=tkwidget(Result,"table",variable=tclarray,rows=7,cols=ncol
(Rematr),selectmode="extended",colwidth=25,background="white")
    tkgrid(table1)
    OnOK=function() { tkdestroy(Result) }
    OK.but=tkbutton(Result,text=" OK ",command=OnOK)
    tkgrid(OK.but) }
  getDMfile <- function() 
  { name=tclvalue(tkgetOpenFile(filetypes="{{All files} *}"))
    if(name=="") return;
    zzF=read.table(name,header=F)
    assign("DMF",zzF,envir=tt$env)
    zzT=read.table(name,header=T)
    assign("DMT",zzT,envir=tt$env) }
  getRVfile=function() 
  { name=tclvalue(tkgetOpenFile(filetypes="{{All files} *}"))
    if(name=="") return;
    zzF=data.matrix(read.table(name,header=F))
    assign("RVF",zzF,envir=tt$env)
    zzT=data.matrix(read.table(name,header=T))
    assign("RVT",zzT,envir=tt$env) }
  PressedUpdate=function()
  { browseURL("http://www.stat.sinica.edu.tw/fredphoa") }
  PressedAbout=function()
  { tkmessageBox(message="SRRS version 0.1 (2013/08) is written, maintained and copyrighted by Frederick Kin Hing Phoa and Shu-
Ching Lin, Institute of Statistical Science, Academia Sinica, Taiwan ROC. We welcome any noncommer- cial uses of this program for 
your own researches. For free resource, the authors guarantee neither the correctness of functions in this program nor 
responsibility for the results of analyses. Please do not redistribute the package in any form without the official permission of 
Frederick Kin Hing Phoa. For redistributions or any commercial uses, please contact Frederick Kin Hing Phoa via email: 
fredphoa@stat.sinica.edu.tw.",icon="info") }
  PressedQ1=function()
  { tkmessageBox(message="Load a text file that contains a design matrix. If headers exist in the first row of the file, click 
the indicator box.",icon="question") }
  PressedQ2=function()
  { tkmessageBox(message="Load a text file that contains a column of response. If response name exists in the first entry, click 
the indicator box.",icon="question") }
  PressedQ3=function()
  { tkmessageBox(message="Enter the number of main effects in the design matrix. All columns in the design matrix are considered 
as main effects if 0 is entered.",icon="question") }
  PressedQ4=function()
  { tkmessageBox(message="Enter a threshold for the termination of factor screening. An automated value is set if 0 is 
entered",icon="question") }
  PressedQ5=function()
  { tkmessageBox(message="If two-factors interacton effects are considered with heredity principle, click the indicator 
box.",icon="question") }
  PressedQ6=function()
  { tkmessageBox(message="Choose the model selection criterion for model searching, either mAIC or AIC.",icon="question") }
  tt=tktoplevel()
  tkwm.title(tt,"The SRRS Method")
  Update.but=tkbutton(tt,text="Update",command=PressedUpdate)
  About.but=tkbutton(tt,text="About",command=PressedAbout)
  tkgrid(tklabel(tt,text="SRRS version 0.1               "),tklabel(tt,text=""),tklabel(tt,text=""),Update.but,About.but,tklabel
(tt,text=""))
  tkgrid(tklabel(tt,text="Data Files:"),sticky="w")
  cbhDM=tkcheckbutton(tt); cbhDMval=tclVar("0")
  tkconfigure(cbhDM,variable=cbhDMval)
  button.DMfile=tkbutton(tt,text="Open File",command=getDMfile)
  Q1.but=tkbutton(tt,text="?",command=PressedQ1)
  tkgrid(tklabel(tt,text="   Design Matrix:"),tklabel(tt,text=""),button.DMfile,tklabel(tt,text="   
Header?"),cbhDM,Q1.but,tklabel(tt,text="")) 
  cbhRV=tkcheckbutton(tt); cbhRVval=tclVar("0")
  tkconfigure(cbhRV,variable=cbhRVval)
  button.RVfile=tkbutton(tt,text="Open File",command=getRVfile)
  Q2.but=tkbutton(tt,text="?",command=PressedQ2)
  tkgrid(tklabel(tt,text="   Response Vector :"),tklabel(tt,text=""),button.RVfile,tklabel(tt,text="   
Header?"),cbhRV,Q2.but,tklabel(tt,text=""))
  tkgrid(tklabel(tt,text="    "))
  tkgrid(tklabel(tt,text="Optional Parameters:"),sticky="w")
  nME=tclVar("0"); Gamma=tclVar("0")
  entry.nME=tkentry(tt,width="5",textvariable=nME,validate="key")
  Q3.but=tkbutton(tt,text="?",command=PressedQ3)
  tkgrid(tklabel(tt,text="   Number of Main Effects:"),tklabel(tt,text=""),tklabel(tt,text=""),entry.nME,tklabel
(tt,text=""),Q3.but,tklabel(tt,text=""))
  entry.gamma=tkentry(tt,width="5",textvariable=Gamma,validate="key")
  Q4.but=tkbutton(tt,text="?",command=PressedQ4)
  tkgrid(tklabel(tt,text="   Specify Gamma :"),tklabel(tt,text=""),tklabel(tt,text=""),entry.gamma,tklabel
(tt,text=""),Q4.but,tklabel(tt,text=""))
  tkgrid(tklabel(tt,text="    "))
  tkgrid(tklabel(tt,text="Selection Conditions:"),sticky="w")
  cb1=tkcheckbutton(tt); cbValue1=tclVar("0")
  tkconfigure(cb1,variable=cbValue1)
  Q5.but=tkbutton(tt,text="?",command=PressedQ5)
  tkgrid(tklabel(tt,text="   2-Factor Interactions?"),tklabel(tt,text=""),tklabel(tt,text=""),cb1,tklabel
(tt,text=""),Q5.but,tklabel(tt,text=""))
  tclRequire("BWidget")
  criterions=c("mAIC","AIC")
  comboBox=tkwidget(tt,"ComboBox",editable=FALSE,values=criterions,width=8,justify="center")
  Q6.but=tkbutton(tt,text="?",command=PressedQ6)
  tkgrid(tklabel(tt,text="   Model Selction Criterion:"),tklabel(tt,text=""),tklabel(tt,text=""),comboBox,tklabel
(tt,text=""),Q6.but,tklabel(tt,text=""))
  tkgrid(tklabel(tt,text="    "))
  OK.but=tkbutton(tt,text="RUN",command=doSRRS)
  tkgrid(tklabel(tt,text=""),tklabel(tt,text=""),OK.but,tklabel(tt,text=""),tklabel(tt,text=""),tklabel(tt,text=""))
}
