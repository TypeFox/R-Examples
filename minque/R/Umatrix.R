
#library(MASS)
#library(lattice)
#library(Matrix)
#library(agridat)
#formula=mod
#mixed.data=function(Y,formula,data=list(),...)
#formula=Brate~Geno|Geno*POS+REP
#data=brate
##lm.data(cbind(y1,y2)~1|G*B)
#formula=y~X1+X2
#data=dat
#mixed_formula(formula)
#model_data(mod1,yv)

###A function to generate a list of information matrices #######
GenerateU=function(mat){
  names=colnames(mat)
  r=length(names)
  U=list(r)
  for(i in 1:r){
    U[[i]]=IndMatrix(mat[,i])
    colnames(U[[i]])=paste(names[i],"(",colnames(U[[i]]),")",sep="")
  }
  names(U)=names
  class(U)="Umatrix"
  return(U)
}

generateu=function(mat){
  names=colnames(mat)
  r=length(names)
  U=list(r)
  for(i in 1:r){
    U[[i]]=IndMatrix(mat[,i])
    colnames(U[[i]])=paste(names[i],"(",colnames(U[[i]]),")",sep="")
  }
  names(U)=names
  class(U)="Umatrix"
  return(U)
}


########################################################
## A very nice function to generate a indicator matrix##
IndMatrix <- function(v)
{  
   nr=length(v)
   if(nr>=2000) stop("Sample size is too large, please contact the author <qgtools@gmail.com> for additional assistance")
   #lv=sort(unique(v))
   lv=unique(v)
   nc=length(lv)
   m=matrix(0,nr,nc)
   for(i in 1:nr){
      j=which(v[i]==lv)
      m[i,j]=1
   }
   colnames(m)=lv
   return(m)
  
   #cl <- as.factor(cl)
   #x <- matrix(0, n, length(levels(cl)) )
   #x[(1:n) + n*(unclass(cl)-1)] <- 1
   #dimnames(x) <- list(names(cl), levels(cl))
   #x
}


rowsum= function(X){sum(X, na.rm=TRUE)}
lm.data=function(formula,data=list(),...){
  a=mixed_formula(formula)
  model=a$model
  yv=a$yv
  mod1=a$mod1
  if(model==2)mod2=a$mod2
  if(is.null(data))a1=model_data(mod1,yv)
  else a1=model_data(mod1,yv,data)
  df1=a1$df
  y=a1$y
  data1=a1$data
  name1=a1$name
  
  if(model==2){
     a2=model_data(mod2,yv,data)
     df2=a2$df
     data2=a2$data
     name2=a2$name
  }
  if(model==1)df=df1  
  if(model==2){
     e=1
     for(i in 1:length(name1)){
       #i=1
       id=which(name2==name1[i])
       if(length(id)>0){
         if(e==1)du_id=id
         else du_id=c(du_id,id)
         e=e+1
       }
     }
     if(e>1)name2=name2[-du_id]
     dfname1=colnames(df1)
     dfname2=colnames(df2)
     e=1
     for(i in 1:length(dfname1)){
       id=which(dfname2==dfname1[i])
       if(length(id)>0){
         if(e==1)id0=id
         else id0=c(id0,id)
         e=e+1
       }
     }
     df=data.frame(df1,df2[,-id0])
     colnames(df)=c(dfname1,dfname2[-id0])
     head(df)
  }
  
  for(i in 1:ncol(df)){
     if(i==1)cls=class(df[,i])
     else cls=c(cls,class(df[,i]))
  }
  name0=colnames(df)
  if(length(name1)==1){
    dat=matrix(1,nrow(df),1)
    #colnames(dat)="mu"
    cls1="numeric"
  }
  if(length(name1)>1){
    d1=fixed_data(name0,name1,cls,df)
    dat1=d1$dat
    cls1=c("numeric",d1$cls)
    dat=data.frame(1,dat1)
    uname1=d1$uname
  }
  if(model==2){
     d2=random_data(name0,name2,cls,df)
     dat2=d2$dat
     cls2=d2$cls
     uname2=d2$uname
     dat=data.frame(dat,dat2)
  }
   
  head(dat)
  #if(length(name1)>0)
  #dat=cbind(1,dat)
  dat=data.frame(dat)
  if(model==1){
     if(length(name1)==1)colnames(dat)="mu"
     else colnames(dat)=c("mu",uname1)
  }
  else if(model==2){
     if(length(name1)==1)colnames(dat)=c("mu",uname2)
     else if(length(name1)>1)colnames(dat)=c("mu",uname1,uname2)
  }
  
  for(i in 2:ncol(df)){
      if(cls[i-1]=="factor")dat[,i]=factor(dat[,i])
  }
  if(model==1)cls0=cls1
  else cls0=c(cls1,cls2)
  
  fid=length(cls1)
  id=which(cls0=="numeric")
  if(fid==1){
    X1=matrix(1,nrow(dat),1)
    colnames(X1)=colnames(dat)[1]
  }
  else{
    X0=dat[,1:fid]
    if(length(id)==1){
      X1=matrix(1,nrow(dat),1)
      colnames(X1)=colnames(dat)[1]
    }
    else X1=X0[,id]
  }
  X2=NULL
  if(fid>ncol(X1)){
     xmat=X0[,-id]
     if(ncol(X0)==(length(id)+1)){
        xmat=data.frame(xmat)
        colnames(xmat)=colnames(X0)[-id]
     }
     X2=GenerateU(xmat)
  }
 
  X=list()
  X4=matrix(0,nrow(X1),ncol(X1))
  for(i in 1:ncol(X1))X4[,i]=as.numeric(X1[,i])
  colnames(X4)=colnames(X1)
  X[[1]]=X4
  xk=length(X2)
  if(xk>0){
    for(i in 1:xk)X[[i+1]]=X2[[i]]
  }
  class(X)="Umatrix"
  names(X)=c("B0",names(X2))
  if(model==2){
     rmat=dat[,-(1:fid)]
     if(ncol(dat)==(fid+1)){
        rmat=data.frame(rmat)
        colnames(rmat)=colnames(dat)[fid+1]
     }
     U=GenerateU(rmat)
     class(U)="Umatrix"
  }
  #colnames(U[[1]])
  else U=NULL
  TraitNames=colnames(y)
  result=list(Y=y,dat=dat,X=X,U=U,class=cls0,TraitNames=TraitNames)
  return(result)
}


mixed_formula=function(formula,...){
  s1=toString(formula)
  s1=strsplit(s1,",")[[1]]
  mc=length(s1)  
  
  s2=strsplit(s1[mc],"\\+")[[1]]
  id=grep("\\|",s2)
  if(length(id)==0)model=1
  else if(length(id)>0) model=2
  
  ##to get two model equations
  if(model==1)mod1=s1[mc]
  if(model==2){
     fullmod=strsplit(s1[mc],"\\|")[[1]]
     mod1=fullmod[1]
     mod2=fullmod[2]
  }
  
  ##to get names of Y used for model statement
  s=toString(formula)
  s=gsub("\\~,","",s)
  s=gsub(" ","",s)
  id1=grep("),",s)
  if(length(id1)==0){
    s2=strsplit(s,",")[[1]]
    yv=s2[1]
  }
  if(length(id1)>0){
    s2=strsplit(s,"),")[[1]]
    #length(s2)
    s2[1]=paste(s2[1],")",sep="")
    yv=s2[1]
  }
  if(model==1)res=list(mod1=mod1,yv=yv,model=model)
  else res=list(mod1=mod1,mod2=mod2,yv=yv,model=model)
  return(res)
}
#mod=mod1
model_data=function(mod,yv,data=NULL,...){
  #mod1=gsub(" ","",mod)
  #yv="Y"
  names=c(yv,mod)
  fm1=paste(names, collapse="~")
  fm1=as.formula(fm1)
  if(is.null(data))data=NULL
  if(is.null(data))df1=model.frame(fm1)
  else df1=model.frame(fm1,data)
  
  dat1=model.matrix(fm1,data)
  y=df1[,1]
  y=data.frame(y)
  if(ncol(y)==1)colnames(y)=yv
  name1=colnames(dat1)
  
  a=list(df=df1,data=dat1,y=y,name=name1)
  return(a)
}

random_data=function(name0,name2,cls,df,...){
  M=matrix(0,length(name0),length(name2))
  for(i in 1:length(name0)){
     id=grep(name0[i],name2)
     M[i,id]=1
  }
  for(i in 1:length(name2)){
     id=which(M[,i]==1)
     str=paste(name0[id],collapse=":")
     if(i==1)vname=str
     else vname=c(vname,str)
  }
  uname2=unique(vname)
  M2=matrix(0,length(name0),length(uname2))
  for(i in 1:length(uname2)){
     id=which(vname==uname2[i])
     if(length(id)>0)M2[,i]=M[,id[1]]
     
  }
  colnames(M2)=uname2
  cls2=NULL
  for(i in 1:length(uname2)){
     #i=5
     id=which(M2[,i]==1)
     a=CFactor(df,id)
     cls2=c(cls2,"factor")
     if(i==1)dat=a
     else dat=data.frame(dat,a)
  }
  if(length(uname2)==1){dat=data.frame(dat);colnames(dat)=uname2}
  head(dat)
  a=list(dat=dat,cls=cls2,uname=uname2)
  return(a)
}

fixed_data=function(name0,name1,cls,df,...){
  M=matrix(0,length(name0),length(name1))
  for(i in 1:length(name0)){
     id=grep(name0[i],name1)
     M[i,id]=1
  }
  id=which(name1=="(Intercept)")
  if(length(id)>0){
     name1=name1[-id]
     M=M[,-id]
     if(length(name1)==1){
        M=data.frame(M)
        colnames(M)=name1
     }
     #class(M)
     
  }
  if(length(name1)==0){
    X1=rep(1,nrow(df))
    X1=data.frame(X1)
    cls1="numeric"
    dat=X1
  }
  else if(length(name1)>0){
    for(i in 1:length(name1)){
       id=which(M[,i]==1)
       str=paste(name0[id],collapse=":")
       if(i==1)vname=str
       else vname=c(vname,str)
    }
    uname1=unique(vname)
    M1=matrix(0,length(name0),length(uname1))
    for(i in 1:length(uname1)){
       id=which(vname==uname1[i])
       if(length(id)>0)M1[,i]=M[,id[1]]
     
    }
    colnames(M1)=uname1
    cls1=NULL
    for(i in 1:length(uname1)){
       #i=5
       id=which(M1[,i]==1)
       #if(length(id)==1)a=df[,id]
       #else if(length(id)>1){
       clst=cls[id]
       idn=which(clst=="numeric")
       idi=which(clst=="integer")
       idt=unique(c(idn,idi))
       if(length(idt)==length(id)){
          a=CNumeric(df,id)
          cls1=c(cls1,"numeric")
       }
       else {
          a=CFactor(df,id)
          cls1=c(cls1,"factor")
       }
       if(i==1)dat=a
       else dat=data.frame(dat,a)
    }
    if(length(uname1)==1){dat=data.frame(dat);colnames(dat)=uname1}
    head(dat)
  }
  a=list(dat=dat,cls=cls1,uname=uname1)
  return(a)
}


lm.data2=function(formula,data=list(),...){
  #formula=Brate~Geno*POS|Geno/REP

  s1=toString(formula)
  s1=strsplit(s1,",")[[1]]
  mc=length(s1)     ##model components to determine if a responsible variable included
  #if(mc==3)yv=s1[2]
  
  #mo=strsplit(s1[3],",")[[1]][3]
  s2=strsplit(s1[mc],"\\+")[[1]]
  id=grep("\\|",s2)
  if(length(id)==0){
     model=1
     #id=length(s2)
  }
  if(length(id)>0) model=2
  #mod=gsub("\\|","\\+",s1[mc])
  #mod=gsub(" ","",mod)
  #fo=paste("~",mod,sep="")
  #as.formula(fo)
  #s3=strsplit(mod,"\\+")[[1]]
  #s3=gsub(" ","",s3)
  #if(s3[1]=="1")flag=1
  #else flag=0
  #if(flag==0)id=id+1
  #fo=as.formula(fo)

  ##to get two model equations
  if(model==1)mod1=s1[mc]
  if(model==2){
     fullmod=strsplit(s1[mc],"\\|")[[1]]
     mod1=fullmod[1]
     mod2=fullmod[2]
  }
  
  ##to get names of Y used for model statement
  s=toString(formula)
  s=gsub("\\~,","",s)
  s=gsub(" ","",s)
  id1=grep("),",s)
  if(length(id1)==0){
    s2=strsplit(s,",")[[1]]
    yv=s2[1]
  }
  if(length(id1)>0){
    s2=strsplit(s,"),")[[1]]
    #length(s2)
    s2[1]=paste(s2[1],")",sep="")
    yv=s2[1]
  }
  ## generate model.frame and model.matrix for fixed model and random model##
  mod1=gsub(" ","",mod1)
  names=c(yv,mod1)
  fm1=paste(names, collapse="~")
  fm1=as.formula(fm1)
  if(is.null(data))data=NULL
  #if(is.null(data)){
  #   df1=model.frame(fm1)
  #   data1=model.matrix(fm11)
  #}
  #else{
     df1=model.frame(fm1,data)
     data1=model.matrix(fm1,data)
  #}
  #class(df1[,2])
  ##colnames(df1)
  #TraitNames=colnames(df1)[1]
  y=df1[,1]
  y=data.frame(y)
  if(ncol(y)==1)colnames(y)=yv
  name1=colnames(data1)
  colnames(df1)
  
  
  if(model==2){
     names=c(yv,mod2)
     fm2=paste(names, collapse="~")
     fm2=as.formula(fm2)
     #if(is.null(data)){
     #   df2=model.frame(fm2)
     #   data2=model.matrix(fm2)
     #}
     #else{
        df2=model.frame(fm2,data)
        data2=model.matrix(fm2,data)
     #}
     name2=colnames(data2)
  }
  if(model==2){
    e=1
    for(i in 1:length(name1)){
      #i=1
      id=which(name2==name1[i])
      if(length(id)>0){
         if(e==1)du_id=id
         else du_id=c(du_id,id)
         e=e+1
      }
    }
    if(e>1)name2=name2[-du_id]
    dfname1=colnames(df1)
    dfname2=colnames(df2)
    e=1
    for(i in 1:length(dfname1)){
      id=which(dfname2==dfname1[i])
      if(length(id)>0){
         if(e==1)id0=id
         else id0=c(id0,id)
         e=e+1
      }
    }
    df=data.frame(df1,df2[,-id0])
    colnames(df)=c(dfname1,dfname2[-id0])
  }
  else if(model==1){
    df=df1
  }  
  for(i in 1:ncol(df)){
     if(i==1)cls=class(df[,i])
     else cls=c(cls,class(df[,i]))
  }
  name0=colnames(df)
  ##Get data colums for fixed effects##
  M=matrix(0,length(name0),length(name1))
  for(i in 1:length(name0)){
     id=grep(name0[i],name1)
     M[i,id]=1
  }
  id=which(name1=="(Intercept)")
  if(length(id)>0){
     name1=name1[-id]
     M=M[,-id]
     if(length(name1)==1){
        M=data.frame(M)
        colnames(M)=name1
     }
     #class(M)
     
  }
  if(length(name1)==0){
    X1=rep(1,nrow(df))
    X1=data.frame(X1)
    cls1="numeric"
    dat=X1
  }
  else if(length(name1)>0){
    for(i in 1:length(name1)){
       id=which(M[,i]==1)
       str=paste(name0[id],collapse=":")
       if(i==1)vname=str
       else vname=c(vname,str)
    }
    uname1=unique(vname)
    M1=matrix(0,length(name0),length(uname1))
    for(i in 1:length(uname1)){
       id=which(vname==uname1[i])
       if(length(id)>0)M1[,i]=M[,id[1]]
     
    }
    colnames(M1)=uname1
    cls1=NULL
    for(i in 1:length(uname1)){
       #i=5
       id=which(M1[,i]==1)
       #if(length(id)==1)a=df[,id]
       #else if(length(id)>1){
       clst=cls[id]
       idn=which(clst=="numeric")
       idi=which(clst=="integer")
       idt=unique(c(idn,idi))
       if(length(idt)==length(id)){
          a=CNumeric(df,id)
          cls1=c(cls1,"numeric")
       }
       else {
          a=CFactor(df,id)
          cls1=c(cls1,"factor")
       }
       if(i==1)dat=a
       else dat=data.frame(dat,a)
    }
    if(length(uname1)==1){dat=data.frame(dat);colnames(dat)=uname1}
  }
  
  ##colnames(dat);head(dat)
  ##Get data columns for random effects ####
  if(model==2){
     M=matrix(0,length(name0),length(name2))
     for(i in 1:length(name0)){
        id=grep(name0[i],name2)
        M[i,id]=1
     }
     for(i in 1:length(name2)){
        id=which(M[,i]==1)
        str=paste(name0[id],collapse=":")
        if(i==1)vname=str
        else vname=c(vname,str)
     }
     uname2=unique(vname)
     M2=matrix(0,length(name0),length(uname2))
     for(i in 1:length(uname2)){
        id=which(vname==uname2[i])
        if(length(id)>0)M2[,i]=M[,id[1]]
     }
     colnames(M2)=uname2
     cls2=NULL
     for(i in 1:length(uname2)){
        id=which(M2[,i]==1)
        a=CFactor(df,id)
        #if(i==1)dat=a
        dat=data.frame(dat,a)
        cls2=c(cls2,"factor")
     }
  }
  #colnames(dat)=c(uname1,uname2)
  #names(df)
  #ncol(dat)
  #cls1
  
  
  head(dat)
  if(length(name1)>0)dat=data.frame(1,dat)
  dat=data.frame(dat)
  if(model==1)colnames(dat)=c("mu",uname1)
  if(length(name1)==0)colnames(dat)=c("mu",uname2)
  else if(model==2&&length(name1)>1)colnames(dat)=c("mu",uname1,uname2)
  for(i in 1:length(df[,-1])){
      if(cls[i]=="factor")dat[,i+1]=factor(dat[,i+1])
  }
  if(model==1)cls0=c("numeric",cls1)
  else cls0=c("numeric",cls1,cls2)
  
 
  if(length(name1)==0)fid=1
  else fid=length(cls1)+1
  id=which(cls0=="numeric")
 
  if(fid==1){
    #X1=dat[,1]
    X1=matrix(1,nrow(dat),1)
    colnames(X1)=colnames(dat)[1]
  }
  else{
    X0=dat[,1:fid]
    if(length(id)==1){
      X1=matrix(1,nrow(dat),1)
      colnames(X1)=colnames(dat)[1]
    }
    else X1=X0[,id]
  }
  X2=NULL
  if(fid>ncol(X1)){
     xmat=X0[,-id]
     if(ncol(X0)==(length(id)+1)){
        xmat=data.frame(xmat)
        colnames(xmat)=colnames(X0)[-id]
     }
     X2=GenerateU(xmat)
  }
 
  X=list()
  X[[1]]=X1
  xk=length(X2)
  if(xk>0){
    for(i in 1:xk)X[[i+1]]=X2[[i]]
  }
  class(X)="Umatrix"
  if(model==2){
     rmat=dat[,-(1:fid)]
     if(ncol(dat)==(fid+1)){
        rmat=data.frame(rmat)
        colnames(rmat)=colnames(dat)[fid+1]
     }
     U=GenerateU(rmat)
  }
  #colnames(U[[1]])
  else U=NULL
  TraitNames=colnames(y)
  result=list(Y=y,dat=dat,X=X,U=U,class=cls0,TraitNames=TraitNames)
  return(result)
}



lm.data0=function(formula=formula,data=list(),...){
   inf=lm.conversion0(formula)
   fm=inf$formula
   #fm=gsub("~","",fm)
   yv=inf$yv
   fid=inf$fid
   mc=inf$mc
   model=inf$mod
   #fm=gsub("1+","\\",fm)
   #fm=gsub("~+","\\",fm)
   fm1=paste(yv,fm,sep="")
   #fm1=gsub("\1+","",fm1)
   
   fm1=as.formula(fm1)
   if(is.null(data))dxf=model.frame(fm1)  ##ncol(dxf)
   else  dxf=model.frame(fm1,data)
   #fm1=paste(yv,fm,sep="")
   #fm1=as.formula(fm1)
   #dxf=model.frame(fm1)
   Y=model.response(dxf)
   #dim(Y)
   #traitnum=ncol(Y)
   
   #ncol(dx)
   #class(dx[,3])
   d2=model.matrix(fm1,dxf)
   #fm=as.formula(fm)
   #d2=model.matrix(fm,dx)
   dx=dxf[,-1]
   
   if(model==1){
      mu=1
      Xmatrix=d2
      dat=data.frame(mu,dx)

   }
   else if(model==2){
     if(mc==2)dx1=dxf
     else if(mc>2)dx1=dxf[,-1]
   
     for(i in 1:ncol(dx1)){
       if(i==1)cls=class(dx1[,i])
       else cls=c(cls,class(dx1[,i]))
     }
     if(mc>2){
       y=dxf[,1]
       name0=colnames(dxf)[-1]
       name1=colnames(d2)
     }
     else{
       name0=colnames(dxf)
       name1=colnames(d2)
     }
     M=matrix(0,length(name0),length(name1))
     for(i in 1:length(name0)){
       id=grep(name0[i],name1)
       M[i,id]=1
     }
     id=which(name1=="(Intercept)")
     if(length(id)>0){
       M=M[,-id]
       name1=name1[-id]
     }
     for(i in 1:length(name1)){
       id=which(M[,i]==1)
       str=paste(name0[id],collapse=":")
       if(i==1)vname=str
       else vname=c(vname,str)
     }
     uname=unique(vname)
     M1=matrix(0,length(name0),length(uname))

     for(i in 1:length(uname)){
       id=which(vname==uname[i])
       if(length(id)>0){
         M1[,i]=M[,id[1]]
       }
    }
    colnames(M1)=uname

    for(i in 1:length(uname)){
       id=which(M1[,i]==1)
       a=CFactor(dx1,id)
       if(i==1)dat=a
       else dat=data.frame(dat,a)
    }
    #if(mc==2)dat=data.frame(1,dat)
    #else 
    dat=data.frame(1,dat)
    dat=data.frame(dat)
    for(i in 1:length(dx1)){
      if(cls[i]=="factor")dat[,i+1]=factor(dat[,i+1])
    }
    #if(mc==2)
    colnames(dat)=c("mu",uname)
    #else colnames(dat)=c("mu",uname,names(dx)[1])
    #res=list(ndata=dat,fid=inf[[2]]+1)
    #fid=inf[[2]]
    #fid=fid+1
    names=colnames(dat)
    xnames=names[1:fid]
    xnames=xnames[-1]
    if(length(xnames)==0){
      Xmatrix=matrix(1,nrow(dat),1)
      Xmatrix=data.frame(Xmatrix)
      colnames(Xmatrix)="Intercept"
    }
    else if(length(xnames)>0){
      #if(mc>2)yname=names[length(names)]
      fmla <- as.formula(paste("~", paste(xnames, collapse= "+")))
      #model.frame(fmla,data=a$ndata)
      #Xmatrix=model.matrix(fmla,data=dat)
      if(is.null(data))Xmatrix=model.matrix(fmla,data)
      else Xmatrix=model.matrix(fmla,data)
    }
    if(mc==2){
      dr=dat[,((fid+1):length(names))]
      dr=data.frame(dr)
      colnames(dr)=names[(fid+1):length(names)]
      Umatrix=GenerateU(dr)
     
    }
    else if(mc>2){
      dr=dat[,((fid+1):length(names))]
      dr=data.frame(dr)
      colnames(dr)=names[(fid+1):length(names)]
      Umatrix=GenerateU(dr)
    }
  } ##end of loop else if(model==2)
  
  if(model==1)res=list(ndata=dat,Y=Y,X=Xmatrix,fid=(fid+1),model=model)
  else if(model==2)res=list(ndata=dat,Y=Y,X=Xmatrix,U=Umatrix,fid=fid,model=model)

  ##res=list(ndata=dat,X=Xmatrix,U=Umatrix,fid=fid)
  return(res)
}
#formula=mod1

lm.conversion0=function(formula=formula,...){
  s1=toString(formula)
  s1=strsplit(s1,",")[[1]]
  mc=length(s1)     ##model components to determine if a responsible variable included
  #if(mc==3)yv=s1[2]
  
  #mo=strsplit(s1[3],",")[[1]][3]
  s2=strsplit(s1[mc],"\\+")[[1]]
  id=grep("\\|",s2)
  if(length(id)==0){
     model=1
     id=length(s2)
  }
  else model=2
  mod=gsub("\\|","\\+",s1[mc])
  #if(length(mod)==1)
  mod=gsub(" ","",mod)
  #mod=gsub("1+","",mod)
  #if(mc==3)fo=paste(yv,"~",mod,sep="")
  #else 
  fo=paste("~",mod,sep="")
  as.formula(fo)
  s3=strsplit(mod,"\\+")[[1]]
  s3=gsub(" ","",s3)
  if(s3[1]=="1")flag=1
  else flag=0
  if(flag==0)id=id+1
  #fo=as.formula(fo)

  #fullmod=strsplit(s1[mc],"\\|")
  
  
  
  s=toString(formula)
  s=gsub("\\~,","",s)
  s=gsub(" ","",s)
  id1=grep("),",s)
  if(length(id1)==0){
    s2=strsplit(s,",")[[1]]
    yv=s2[1]
  }
  else if(length(id1)>0){
    s2=strsplit(s,"),")[[1]]
    #length(s2)
    s2[1]=paste(s2[1],")",sep="")
    yv=s2[1]
  }
  f=list(yv=yv,formula=fo,fid=id,mc=mc,mod=model)
  return(f)
}



mixed <- function(x, ...) UseMethod("mixed")
#a=mixed(y~1|G*B)
#formula=y~1|G*B
mixed.formula=function(formula,data=list(),...){
   ff=formula
   s1=toString(ff)
   #yv=s1[2]
   s1=gsub("\\|","+",s1)
   s1=gsub("\\~,","",s1)
   s1=gsub("\\,","~",s1)
   fo=as.formula(s1)
   if(is.null(data)==TRUE)d1=model.frame(fo)
   else if(is.null(data)==FALSE)d1=data
   #mf=model.frame(formula=ff,data=d1)
   #y=mf[,1]
   s1=toString(ff)
   s2=strsplit(s1,",")[[1]]
   s3=s2[length(s2)]
   yv=s2[2]

   flag=grep("\\|",s3)
   if(length(flag)==0)flag=0
   else if(length(flag)>0)flag=1
   f1=strsplit(s3,"\\|")[[1]]
   if(flag==0){
      mod=1
      ff1=paste(yv,"~", f1[1],sep="")
      ff1=as.formula(ff1)
      dx=model.frame(ff1,d1)
      y=model.response(dx)
      fdata=model.matrix(ff1,dx)
   }
   else if(flag==1){
     mod=2
     #f1=strsplit(s3,"\\|")[[1]]
     f11=f1[1]
     ff1=paste(yv,"~", f11,sep="")
     ff1=as.formula(ff1)
     dx=model.frame(ff1,d1)
     y=model.response(dx)
     fdata=model.matrix(ff1,dx)
     ff2=paste("~", f1[2],sep="")
     rdata=model.rdata(ff2,d1)
     U=GenerateU(rdata)
   }
   #y=as.numeric(y)
   if(mod==1)dat=list(Y=y,X=fdata,U=NULL,mod=mod)
   else if(mod>1)dat=list(Y=y,X=fdata,U=U,mod=mod)
   #length(y)
   #au=rep(1,length(U)+1)
   #res=MINQUE(au,U,fdata,y)
   #y=data.frame(y)
   #res=RandJackReal(U,y,JacNum=10,JacRep=10,X)
   return(dat)
}




model.rdata=function(formula,data=list(),...){

   d2=model.frame(formula=formula,data=data)
   #y=model.response(d2)
   dnames=names(d2)
   dn=d2
   s1=toString(formula)
   s2=strsplit(s1,",")
   s3=s2[[1]][length(s2[[1]])]
   s3=gsub("~","",s3)
   s4=strsplit(s3," \\+")[[1]]
   
   s4=gsub(" ","",s4)

   ID=matrix(0,length(dnames),length(s4))
   mk=length(s4)
   for(i in 1:length(dnames)){
      v=numeric(mk)
      id=grep(dnames[i],s4)
      v[id]=1
      ID[i,]=v
   }
   vnames=NULL
   d3=NULL
   for(i in 1:mk){
  
     id=which(ID[,i]==1)
     vn=paste(dnames[id], collapse= ":")
     if(i==1){
        vnames=vn
        d3=CFactor(dn,id)
     }
  
     else if(i>1){
        vnames=c(vnames,vn)
        d3=data.frame(d3,CFactor(dn,id))
    }
  }
  dat=data.frame(d3)
  colnames(dat)=vnames
  return(dat)
}
 


reg.formula=function(x,y){
   fmla <- as.formula(paste("y ~ ", paste(x, collapse= "+")))
   return(fmla)
} 


print.Umatrix=function(U,...){
   mk=length(U)
   cat("These are very complicated matrices:\n")
   cat("We only provide summarized information:\n")
   for(i in 1:mk){
      cat("Dimensions for matrix", i,":", "row=", nrow(U[[i]]), "col=", ncol(U[[i]]), "\n")
   }  
}
print.Matrix=function(X,...){
   row=nrow(X)
   if(row>=100){
      cat("This is very large data set:\n")
      cat("We only provide summarized information.\n")
      cat("The sample size is ", nrow(X), " and the number of variables is ", ncol(X), "\n")
   }
}   

########################################################
## A very nice function to generate a indicator matrix##
IndMatrix1 <- function(v)
{  
   if(class(v)=="factor")v=as.character(v)
   nr=length(v)
   lv=sort(unique(v))
   nc=length(lv)
   m=matrix(0,nr,nc)
   for(i in 1:nc){
      id=which(v==lv[i])
      m[id,i]=1
   }
   colnames(m)=lv
   return(m)
   
   #cl <- as.factor(cl)
   #x <- matrix(0, n, length(levels(cl)) )
   #x[(1:n) + n*(unclass(cl)-1)] <- 1
   #dimnames(x) <- list(names(cl), levels(cl))
   #x
}

######################################################

#### Implementation of design matrices for GYL model.

#######################################################################################
#######################################################################################
#######################################################################################

GYLMatrix=function(Year,Loc,Geno,Block){
   Y=Year
   L=Loc
   G=Geno
   B=Block
   y0=length(unique(Y))
   b0=length(unique(B))
   if(y0==1){
     if(b0==1)mat=data.frame(L,G)
     else if(b0>1){
        B=paste(G,B,sep="*")
        GL=paste(G,L,sep="*")
        mat=data.frame(L,G,GL,B)
     }
   }
   else if(y0>1){
      if(b0==1){
         GY=paste(G,Y,sep="*")
         GL=paste(G,L,sep="*")
         YL=paste(Y,L,sep="*")
         mat=data.frame(Y,L,G,YL,GY,GL)
      }
      if(b0>1){
         GY=paste(G,Y,sep="*")
         GL=paste(G,L,sep="*")
         YL=paste(Y,L,sep="*")
         GYL=paste(G,Y,L,sep="*")
         B=paste(Y,L,B,sep="*")
         mat=data.frame(Y,L,G,YL,GY,GL,GYL,B)
      }
   }
   U=GenerateU(mat)
   return(U)
}

#### Implementation of design matrices for 4-way AD model.

#######################################################################################
#######################################################################################
#######################################################################################

AD4Matrix=function(Env,F1,M1,F2,M2,Gen,Block,BlockID,CC=NULL,...){
   if(is.null(CC))CC=numeric(length(Env))
   Year=Env
   TP=sort(unique(c(F1,M1,F2,M2)))
   Ty=sort(unique(Year))
   Tb=unique(Block)
   
   n=length(F1)

   y0=length(Ty)
   p0=length(TP)
   b0=length(Tb)
   if(BlockID==1&&b0==1)BlockID=0

   ##########################################
   ## Additive effects#######################
   UA=matrix(0,nrow=n,ncol=p0)
   dim(UA)
   colnames(UA)=TP
   for(i in 1:n){
       if(Gen[i]==0){
           j=which(colnames(UA)==F1[i])
           UA[i,j]=2
       }
       else if(Gen[i]>0){
           if(CC[i]==0){
              j=which(colnames(UA)==F1[i])
              UA[i,j]=UA[i,j]+0.5
              j=which(colnames(UA)==M1[i])
              UA[i,j]=UA[i,j]+0.5
              j=which(colnames(UA)==F2[i])
              UA[i,j]=UA[i,j]+0.5
              j=which(colnames(UA)==M2[i])
              UA[i,j]=UA[i,j]+0.5
          }
          else{
             j=which(colnames(UA)==F1[i])
             UA[i,j]=UA[i,j]+0.25
             j=which(colnames(UA)==M1[i])
             UA[i,j]=UA[i,j]+0.25
             j=which(colnames(UA)==F2[i])
             UA[i,j]=UA[i,j]+0.5
             j=which(colnames(UA)==M2[i])
             UA[i,j]=UA[i,j]+1.0
         
         }

       }
   }
   colnames(UA)=paste("A(",TP,")",sep="")
   #head(UA)
   UA=DropColumns(UA)
   
   ##################################################
   ## Dominance effects  ############################
   TD=NULL

   TD=paste(TP,TP,sep="*")

   for(i in 1:(p0-1)){
       str=paste(TP[i],TP[(i+1):p0],sep="*")
       TD=c(TD,str)
   }

   UD=matrix(0,nrow=n,ncol=length(TD))
   
   colnames(UD)=TD

   for(i in 1:n){
      if(CC[i]==0){
          if(Gen[i]==0){
             j= which(colnames(UD)==paste(F1[i], F1[i], sep="*"))
             UD[i,j]=1
          }
      
          else if(Gen[i]==1){
            j= which(colnames(UD)==paste(F1[i], F2[i], sep="*"))
            UD[i,j]=UD[i,j]+0.25
          
            j= which(colnames(UD)==paste(F1[i], M2[i], sep="*"))
            UD[i,j]=UD[i,j]+0.25
          
            j= which(colnames(UD)==paste(M1[i], F2[i], sep="*"))
            UD[i,j]=UD[i,j]+0.25
          
            j= which(colnames(UD)==paste(M1[i], M2[i], sep="*"))
            UD[i,j]=UD[i,j]+0.25
         }

         else if(Gen[i]==2){
            j= which(colnames(UD)==paste(F1[i], F1[i], sep="*"))
            UD[i,j]=UD[i,j]+0.25/2
          
            j= which(colnames(UD)==paste(M1[i], M1[i], sep="*"))
            UD[i,j]=UD[i,j]+0.25/2
          
            j= which(colnames(UD)==paste(F2[i], F2[i], sep="*"))
            UD[i,j]=UD[i,j]+0.25/2
          
            j= which(colnames(UD)==paste(M2[i], M2[i], sep="*"))
            UD[i,j]=UD[i,j]+0.25/2

            j= which(colnames(UD)==paste(F1[i], F2[i], sep="*"))
            UD[i,j]=UD[i,j]+0.25/2
          
            j= which(colnames(UD)==paste(F1[i], M2[i], sep="*"))
            UD[i,j]=UD[i,j]+0.25/2
          
            j= which(colnames(UD)==paste(M1[i], F2[i], sep="*"))
            UD[i,j]=UD[i,j]+0.25/2
          
            j= which(colnames(UD)==paste(M1[i], M2[i], sep="*"))
            UD[i,j]=UD[i,j]+0.25/2
        }
      
        else if(Gen[i]==3){
           j= which(colnames(UD)==paste(F1[i], F1[i], sep="*"))
           UD[i,j]=UD[i,j]+3/16
          
           j= which(colnames(UD)==paste(M1[i], M1[i], sep="*"))
           UD[i,j]=UD[i,j]+3/16
          
           j= which(colnames(UD)==paste(F2[i], F2[i], sep="*"))
           UD[i,j]=UD[i,j]+3/16
          
           j= which(colnames(UD)==paste(M2[i], M2[i], sep="*"))
           UD[i,j]=UD[i,j]+3/16

           j= which(colnames(UD)==paste(F1[i], F2[i], sep="*"))
           UD[i,j]=UD[i,j]+0.25/4
          
           j= which(colnames(UD)==paste(F1[i], M2[i], sep="*"))
           UD[i,j]=UD[i,j]+0.25/4
          
           j= which(colnames(UD)==paste(M1[i], F2[i], sep="*"))
           UD[i,j]=UD[i,j]+0.25/4
          
           j= which(colnames(UD)==paste(M1[i], M2[i], sep="*"))
           UD[i,j]=UD[i,j]+0.25/4 
       }
     }
     else{
          if(Gen[i]==0){
             j= which(colnames(UD)==paste(F1[i], F1[i], sep="*"))
             UD[i,j]=1
          }
          else if(Gen[i]==1){
            j= which(colnames(UD)==paste(F1[i], M2[i], sep="*"))
            UD[i,j]=UD[i,j]+0.25
          
            j= which(colnames(UD)==paste(M1[i], M2[i], sep="*"))
            UD[i,j]=UD[i,j]+0.25
          
            j= which(colnames(UD)==paste(F2[i], M2[i], sep="*"))
            UD[i,j]=UD[i,j]+0.50
          
            #j= which(colnames(UD)==paste(M1[i], M2[i], sep="*"))
            #UD[i,j]=UD[i,j]+0.25
         }

         else if(Gen[i]==2){
            j= which(colnames(UD)==paste(F1[i], F1[i], sep="*"))
            UD[i,j]=UD[i,j]+1/16
          
            j= which(colnames(UD)==paste(M1[i], M1[i], sep="*"))
            UD[i,j]=UD[i,j]+1/16
          
            j= which(colnames(UD)==paste(F2[i], F2[i], sep="*"))
            UD[i,j]=UD[i,j]+1/8
          
            j= which(colnames(UD)==paste(M2[i], M2[i], sep="*"))
            UD[i,j]=UD[i,j]+1/4

            j= which(colnames(UD)==paste(F1[i], M2[i], sep="*"))
            UD[i,j]=UD[i,j]+1/8
          
            j= which(colnames(UD)==paste(M1[i], M2[i], sep="*"))
            UD[i,j]=UD[i,j]+1/8
          
            j= which(colnames(UD)==paste(F2[i], M2[i], sep="*"))
            UD[i,j]=UD[i,j]+0.25
            
        }
      
        else if(Gen[i]==3){
           j= which(colnames(UD)==paste(F1[i], F1[i], sep="*"))
            UD[i,j]=UD[i,j]+1/16+1/32
          
            j= which(colnames(UD)==paste(M1[i], M1[i], sep="*"))
            UD[i,j]=UD[i,j]+1/16+1/32
          
            j= which(colnames(UD)==paste(F2[i], F2[i], sep="*"))
            UD[i,j]=UD[i,j]+1/8+1/16
          
            j= which(colnames(UD)==paste(M2[i], M2[i], sep="*"))
            UD[i,j]=UD[i,j]+1/4+1/8

            j= which(colnames(UD)==paste(F1[i], M2[i], sep="*"))
            UD[i,j]=UD[i,j]+1/16
          
            j= which(colnames(UD)==paste(M1[i], M2[i], sep="*"))
            UD[i,j]=UD[i,j]+1/16
          
            j= which(colnames(UD)==paste(F2[i], M2[i], sep="*"))
            UD[i,j]=UD[i,j]+1/8

       }
     }
     
     

   }
   colnames(UD)=paste("D(",TD,")",sep="")
   UD= DropColumns(UD)
   

   ##############################################################
   ### Environmental effects:  random effect#####################
   if(y0>0){
       ### Environmental effects ###################################
       UY=matrix(0,nrow=n,ncol=length(Ty))
       colnames(UY)= Ty
       for(i in 1:n){
	     j=which(colnames(UY)==Year[i])
	     UY[i,j] = 1
       }
       colnames(UY)=paste("Y(",Ty,")",sep="")
       #UY
       #dim(UY)
       UY= DropColumns(UY)
   

       
       ##Additive*environment interaction effects
       
       ## comments: AE interaction
       TAE=NULL
       for(i in 1:y0){
           str=paste(Ty[i],TP[1:p0],sep="*")
           TAE=c(TAE,str)
       }
     
       UAE=matrix(0,nrow=n,ncol=length(TAE))
       #UAE
       colnames(UAE)=TAE
       #head(UAE)
       
       for(i in 1:n){
           if(CC[i]==0){
	        j= which(colnames(UAE)==paste(Year[i], F1[i], sep="*"))
	        UAE[i, j] = UAE[i,j]+0.5

              j= which(colnames(UAE)==paste(Year[i], M1[i], sep="*"))
	        UAE[i, j] = UAE[i,j]+0.5
           
              j= which(colnames(UAE)==paste(Year[i], F2[i], sep="*"))
	        UAE[i, j] = UAE[i,j]+0.5

              j= which(colnames(UAE)==paste(Year[i], M2[i], sep="*"))
	        UAE[i, j] = UAE[i,j]+0.5
           }
           else{
	        j= which(colnames(UAE)==paste(Year[i], F1[i], sep="*"))
	        UAE[i, j] = UAE[i,j]+0.25

              j= which(colnames(UAE)==paste(Year[i], M1[i], sep="*"))
	        UAE[i, j] = UAE[i,j]+0.25
           
              j= which(colnames(UAE)==paste(Year[i], F2[i], sep="*"))
	        UAE[i, j] = UAE[i,j]+0.5

              j= which(colnames(UAE)==paste(Year[i], M2[i], sep="*"))
	        UAE[i, j] = UAE[i,j]+1.0
           }

       }
       colnames(UAE)=paste("AE(",TAE,")",sep="")
       #dim(UAE)

       UAE= DropColumns(UAE)
       #dim(UAE)

       ###########################################################
       ## DE interaction effects #################################
       
       TDE=NULL

       for(i in 1:y0){
           for(j in 1:length(TD)){
              TDE=c(TDE,paste(Ty[i],TD[j],sep="*"))
           }
       }

       UDE=matrix(0,nrow=n,ncol=length(TDE))
       colnames(UDE)=TDE
          
       for(i in 1:n){
           if(Gen[i]==0){
               j= which(colnames(UDE)==paste(Year[i],F1[i], F1[i], sep="*"))
               UDE[i,j]=1

           }
           if(CC[i]==0){
              if(Gen[i]==1){
                 j= which(colnames(UDE)==paste(Year[i],F1[i], F2[i], sep="*"))
                 UDE[i,j]=UDE[i,j]+0.25
          
                 j= which(colnames(UDE)==paste(Year[i],F1[i], M2[i], sep="*"))
                 UDE[i,j]=UDE[i,j]+0.25
          
                 j= which(colnames(UDE)==paste(Year[i],M1[i], F2[i], sep="*"))
                 UDE[i,j]=UDE[i,j]+0.25
          
                 j= which(colnames(UDE)==paste(Year[i],M2[i], M2[i], sep="*"))
                 UDE[i,j]=UDE[i,j]+0.25
              }
             
              else if(Gen[i]==2){
                 j= which(colnames(UDE)==paste(Year[i],F1[i], F1[i], sep="*"))
                 UDE[i,j]=UDE[i,j]+0.25/2
          
                 j= which(colnames(UDE)==paste(Year[i],M1[i], M1[i], sep="*"))
                 UDE[i,j]=UDE[i,j]+0.25/2
          
                 j= which(colnames(UDE)==paste(Year[i],F2[i], F2[i], sep="*"))
                 UDE[i,j]=UDE[i,j]+0.25/2
          
                 j= which(colnames(UDE)==paste(Year[i],M2[i], M2[i], sep="*"))
                 UDE[i,j]=UDE[i,j]+0.25/2

                 j= which(colnames(UDE)==paste(Year[i],F1[i], F2[i], sep="*"))
                 UDE[i,j]=UDE[i,j]+0.25/2
          
                 j= which(colnames(UDE)==paste(Year[i],F1[i], M2[i], sep="*"))
                 UDE[i,j]=UDE[i,j]+0.25/2
          
                 j= which(colnames(UDE)==paste(Year[i],M1[i], F2[i], sep="*"))
                 UDE[i,j]=UDE[i,j]+0.25/2
          
                 j= which(colnames(UDE)==paste(Year[i],M1[i], M2[i], sep="*"))
                 UDE[i,j]=UDE[i,j]+0.25/2

             }

             else if(Gen[i]==3){

                j= which(colnames(UDE)==paste(Year[i],F1[i], F1[i], sep="*"))
                UDE[i,j]=UDE[i,j]+3/16
          
                j= which(colnames(UDE)==paste(Year[i],M1[i], M1[i], sep="*"))
                UDE[i,j]=UDE[i,j]+3/16

          
                j= which(colnames(UDE)==paste(Year[i],F2[i], F2[i], sep="*"))
                UDE[i,j]=UDE[i,j]+3/16

          
                j= which(colnames(UDE)==paste(Year[i],M2[i], M2[i], sep="*"))
                UDE[i,j]=UDE[i,j]+3/16


                j= which(colnames(UDE)==paste(Year[i],F1[i], F2[i], sep="*"))
                UDE[i,j]=UDE[i,j]+1/16
          
                j= which(colnames(UDE)==paste(Year[i],F1[i], M2[i], sep="*"))
                UDE[i,j]=UDE[i,j]+1/16
          
                j= which(colnames(UDE)==paste(Year[i],M1[i], F2[i], sep="*"))
                UDE[i,j]=UDE[i,j]+1/16
          
                j= which(colnames(UDE)==paste(Year[i],M1[i], M2[i], sep="*"))
                UDE[i,j]=UDE[i,j]+1/16
             }  
          }
          else{
              if(Gen[i]==1){
                 j= which(colnames(UDE)==paste(Year[i],F1[i], F2[i], sep="*"))
                 UDE[i,j]=UDE[i,j]+0.25
          
                 j= which(colnames(UDE)==paste(Year[i],M1[i], M2[i], sep="*"))
                 UDE[i,j]=UDE[i,j]+0.25
          
                 j= which(colnames(UDE)==paste(Year[i],F2[i], M2[i], sep="*"))
                 UDE[i,j]=UDE[i,j]+0.50
                          
              }
             
              else if(Gen[i]==2){
                 j= which(colnames(UDE)==paste(Year[i],F1[i], F1[i], sep="*"))
                 UDE[i,j]=UDE[i,j]+1/16
          
                 j= which(colnames(UDE)==paste(Year[i],M1[i], M1[i], sep="*"))
                 UDE[i,j]=UDE[i,j]+1/16
          
                 j= which(colnames(UDE)==paste(Year[i],F2[i], F2[i], sep="*"))
                 UDE[i,j]=UDE[i,j]+1/8
          
                 j= which(colnames(UDE)==paste(Year[i],M2[i], M2[i], sep="*"))
                 UDE[i,j]=UDE[i,j]+1/4

                 j= which(colnames(UDE)==paste(Year[i],F1[i], M2[i], sep="*"))
                 UDE[i,j]=UDE[i,j]+1/8
          
                 j= which(colnames(UDE)==paste(Year[i],M1[i], M2[i], sep="*"))
                 UDE[i,j]=UDE[i,j]+1/8
          
                 j= which(colnames(UDE)==paste(Year[i],F2[i], M2[i], sep="*"))
                 UDE[i,j]=UDE[i,j]+1/4
                          

             }

             else if(Gen[i]==3){

                 j= which(colnames(UDE)==paste(Year[i],F1[i], F1[i], sep="*"))
                 UDE[i,j]=UDE[i,j]+1/16+1/32
          
                 j= which(colnames(UDE)==paste(Year[i],M1[i], M1[i], sep="*"))
                 UDE[i,j]=UDE[i,j]+1/16+1/32
          
                 j= which(colnames(UDE)==paste(Year[i],F2[i], F2[i], sep="*"))
                 UDE[i,j]=UDE[i,j]+1/8+1/16
          
                 j= which(colnames(UDE)==paste(Year[i],M2[i], M2[i], sep="*"))
                 UDE[i,j]=UDE[i,j]+1/4+1/8

                 j= which(colnames(UDE)==paste(Year[i],F1[i], M2[i], sep="*"))
                 UDE[i,j]=UDE[i,j]+1/16
          
                 j= which(colnames(UDE)==paste(Year[i],M1[i], M2[i], sep="*"))
                 UDE[i,j]=UDE[i,j]+1/16
          
                 j= which(colnames(UDE)==paste(Year[i],F2[i], M2[i], sep="*"))
                 UDE[i,j]=UDE[i,j]+1/8
   }  
          }
       }
       colnames(UDE)=paste("DE(",TDE,")",sep="")
       UDE= DropColumns(UDE)
   }

   if(BlockID==1){
       if(y0==1&&b0>1){
           UB=matrix(0,nrow=n,ncol=b0)
           colnames(UB)=Tb
           for(i in 1:n){
                j=which(colnames(UB)==Block[i])
                UB[i,j]=1
           }
           colnames(UB)=paste("B(",Tb,")",sep="")
           #dim(UB)

           UB= DropColumns(UB)
           #dim(UB)
       }
       else if(y0>1&&b0>1){
           TB=NULL

           for(i in 1:y0){
               str=paste(Ty[i],Tb[1:b0],sep="*")
               TB=c(TB,str)
  
           }

           UB=matrix(0,nrow=n,ncol=length(TB))
           UB
           colnames(UB)=TB
         
           for(i in 1:n){
	
	          j= which(colnames(UB)==paste(Year[i],Block[i],sep="*"))
	
	          UB[i, j] = 1
           }
           colnames(UB)=paste("B(",TB,")",sep="")
           
           UB= DropColumns(UB)
           
        }
    }
    if(y0==1&&BlockID==0){
        U=list(UA,UD)
        names(U)=c("A","D")
    }
    else if(y0==1&&BlockID==1){
        U=list(UA,UD,UB)
        names(U)=c("A","D","B")
    }
    else if(y0>1&&BlockID==0){
        U=list(UY,UA,UD,UAE,UDE)
        names(U)=c("Y","A","D","AE","DE")
    }
    else if(y0>1&&BlockID==1){
        U=list(UY,UA,UD,UAE,UDE,UB)
        names(U)=c("Y","A","D","AE","DE","B")
    }
    class(U)="Umatrix"
    return(U)

}


#### Implementation of design matrices for 4-way ADC model.

#######################################################################################
#######################################################################################
#######################################################################################

ADC4Matrix=function(Env,F1,M1,F2,M2,Gen,Block,BlockID){
   Year=Env
   TP=sort(unique(c(F1,M1,F2,M2)))
   Ty=sort(unique(Year))
   Tb=sort(unique(Block))
   
   n=length(F1)

   y0=length(Ty)
   p0=length(TP)
   b0=length(Tb)
   if(BlockID==1&&b0==1)BlockID=0
   
   ##########################################
   ## Cytoplasmic effects#######################
   UC=matrix(0,nrow=n,ncol=p0)
   colnames(UC)=TP
   for(i in 1:n){
       j=which(colnames(UC)==F1[i])
       UC[i,j]=1
   }
   
   colnames(UC)=paste("C",TP,sep="")
   UC=DropColumns(UC)

   ##########################################
   ## Additive effects#######################
   UA=matrix(0,nrow=n,ncol=p0)
   dim(UA)
   colnames(UA)=TP
   for(i in 1:n){
       if(Gen[i]==0){
           j=which(colnames(UA)==F1[i])
           UA[i,j]=2
       }
       else if(Gen[i]>0){
           j=which(colnames(UA)==F1[i])
           UA[i,j]=UA[i,j]+0.5
           j=which(colnames(UA)==M1[i])
           UA[i,j]=UA[i,j]+0.5
           j=which(colnames(UA)==F2[i])
           UA[i,j]=UA[i,j]+0.5
           j=which(colnames(UA)==M2[i])
           UA[i,j]=UA[i,j]+0.5
       }
   }
   colnames(UA)=paste("A(",TP,")",sep="")
   UA=DropColumns(UA)
   
   ##################################################
   ## Dominance effects  ############################
   TD=NULL

   TD=paste(TP,TP,sep="*")

   for(i in 1:(p0-1)){
       str=paste(TP[i],TP[(i+1):p0],sep="*")
       TD=c(TD,str)
   }

   UD=matrix(0,nrow=n,ncol=length(TD))
   
   colnames(UD)=TD

   for(i in 1:n){
      if(Gen[i]==0){
          j= which(colnames(UD)==paste(F1[i], F1[i], sep="*"))
          UD[i,j]=1
      }
      
      else if(Gen[i]==1){
          j= which(colnames(UD)==paste(F1[i], F2[i], sep="*"))
          UD[i,j]=UD[i,j]+0.25
          
          j= which(colnames(UD)==paste(F1[i], M2[i], sep="*"))
          UD[i,j]=UD[i,j]+0.25
          
          j= which(colnames(UD)==paste(M1[i], F2[i], sep="*"))
          UD[i,j]=UD[i,j]+0.25
          
          j= which(colnames(UD)==paste(M1[i], M2[i], sep="*"))
          UD[i,j]=UD[i,j]+0.25
      }

      else if(Gen[i]==2){
          j= which(colnames(UD)==paste(F1[i], F1[i], sep="*"))
          UD[i,j]=UD[i,j]+0.25/2
          
          j= which(colnames(UD)==paste(M1[i], M1[i], sep="*"))
          UD[i,j]=UD[i,j]+0.25/2
          
          j= which(colnames(UD)==paste(F2[i], F2[i], sep="*"))
          UD[i,j]=UD[i,j]+0.25/2
          
          j= which(colnames(UD)==paste(M2[i], M2[i], sep="*"))
          UD[i,j]=UD[i,j]+0.25/2

          j= which(colnames(UD)==paste(F1[i], F2[i], sep="*"))
          UD[i,j]=UD[i,j]+0.25/2
          
          j= which(colnames(UD)==paste(F1[i], M2[i], sep="*"))
          UD[i,j]=UD[i,j]+0.25/2
          
          j= which(colnames(UD)==paste(M1[i], F2[i], sep="*"))
          UD[i,j]=UD[i,j]+0.25/2
          
          j= which(colnames(UD)==paste(M1[i], M2[i], sep="*"))
          UD[i,j]=UD[i,j]+0.25/2
      }
      
      else if(Gen[i]==3){
          j= which(colnames(UD)==paste(F1[i], F1[i], sep="*"))
          UD[i,j]=UD[i,j]+3/16
          
          j= which(colnames(UD)==paste(M1[i], M1[i], sep="*"))
          UD[i,j]=UD[i,j]+3/16
          
          j= which(colnames(UD)==paste(F2[i], F2[i], sep="*"))
          UD[i,j]=UD[i,j]+3/16
          
          j= which(colnames(UD)==paste(M2[i], M2[i], sep="*"))
          UD[i,j]=UD[i,j]+3/16

          j= which(colnames(UD)==paste(F1[i], F2[i], sep="*"))
          UD[i,j]=UD[i,j]+0.25/4
          
          j= which(colnames(UD)==paste(F1[i], M2[i], sep="*"))
          UD[i,j]=UD[i,j]+0.25/4
          
          j= which(colnames(UD)==paste(M1[i], F2[i], sep="*"))
          UD[i,j]=UD[i,j]+0.25/4
          
          j= which(colnames(UD)==paste(M1[i], M2[i], sep="*"))
          UD[i,j]=UD[i,j]+0.25/4 
      }
   }
   colnames(UD)=paste("D(",TD,")",sep="")
   UD= DropColumns(UD)
   

   ##############################################################
   ### Environmental effects:  random effect#####################
   if(y0>0){
       ### Environmental effects ###################################
       UY=IndMatrix(Env)
       colnames(UY)=paste("E(",colnames(UY),")",sep="")
       UY=matrix(0,nrow=n,ncol=length(Ty))
       
       ##Cytoplasm*environment interaction effects
       UCE=UCE2Matrix(Env,F1)
       
       #dim(UAE)
       
       ##Additive*environment interaction effects
       
       ## comments: AE interaction
       TAE=NULL
       for(i in 1:y0){
           str=paste(Ty[i],TP[1:p0],sep="*")
           TAE=c(TAE,str)
       }
     
       UAE=matrix(0,nrow=n,ncol=length(TAE))
       #UAE
       colnames(UAE)=TAE
       #head(UAE)
       for(i in 1:n){
	     j= which(colnames(UAE)==paste(Year[i], F1[i], sep="*"))
	     UAE[i, j] = UAE[i,j]+0.5

           j= which(colnames(UAE)==paste(Year[i], M1[i], sep="*"))
	     UAE[i, j] = UAE[i,j]+0.5
           
           j= which(colnames(UAE)==paste(Year[i], F2[i], sep="*"))
	     UAE[i, j] = UAE[i,j]+0.5

           j= which(colnames(UAE)==paste(Year[i], M2[i], sep="*"))
	     UAE[i, j] = UAE[i,j]+0.5
       }
       colnames(UAE)=paste("AE(",TAE,")",sep="")
       #dim(UAE)

       UAE= DropColumns(UAE)
       #dim(UAE)

       ###########################################################
       ## DE interaction effects #################################
       
       TDE=NULL

       for(i in 1:y0){
           for(j in 1:length(TD)){
              TDE=c(TDE,paste(Ty[i],TD[j],sep="*"))
           }
       }

       UDE=matrix(0,nrow=n,ncol=length(TDE))
       colnames(UDE)=TDE
          
       for(i in 1:n){
           if(Gen[i]==0){
               j= which(colnames(UDE)==paste(Year[i],F1[i], F1[i], sep="*"))
               UDE[i,j]=1

           }

           else if(Gen[i]==1){
              j= which(colnames(UDE)==paste(Year[i],F1[i], F2[i], sep="*"))
              UDE[i,j]=UDE[i,j]+0.25
          
              j= which(colnames(UDE)==paste(Year[i],F1[i], M2[i], sep="*"))
              UDE[i,j]=UDE[i,j]+0.25
          
              j= which(colnames(UDE)==paste(Year[i],M1[i], F2[i], sep="*"))
              UDE[i,j]=UDE[i,j]+0.25
          
              j= which(colnames(UDE)==paste(Year[i],M2[i], M2[i], sep="*"))
              UDE[i,j]=UDE[i,j]+0.25
           }
             
           else if(Gen[i]==2){
              j= which(colnames(UDE)==paste(Year[i],F1[i], F1[i], sep="*"))
              UDE[i,j]=UDE[i,j]+0.25/2
          
              j= which(colnames(UDE)==paste(Year[i],M1[i], M1[i], sep="*"))
              UDE[i,j]=UDE[i,j]+0.25/2
          
              j= which(colnames(UDE)==paste(Year[i],F2[i], F2[i], sep="*"))
              UDE[i,j]=UDE[i,j]+0.25/2
          
              j= which(colnames(UDE)==paste(Year[i],M2[i], M2[i], sep="*"))
              UDE[i,j]=UDE[i,j]+0.25/2

              j= which(colnames(UDE)==paste(Year[i],F1[i], F2[i], sep="*"))
              UDE[i,j]=UDE[i,j]+0.25/2
          
              j= which(colnames(UDE)==paste(Year[i],F1[i], M2[i], sep="*"))
              UDE[i,j]=UDE[i,j]+0.25/2
          
              j= which(colnames(UDE)==paste(Year[i],M1[i], F2[i], sep="*"))
              UDE[i,j]=UDE[i,j]+0.25/2
          
              j= which(colnames(UDE)==paste(Year[i],M1[i], M2[i], sep="*"))
              UDE[i,j]=UDE[i,j]+0.25/2

           }

           else if(Gen[i]==3){

              j= which(colnames(UDE)==paste(Year[i],F1[i], F1[i], sep="*"))
              UDE[i,j]=UDE[i,j]+3/16
          
              j= which(colnames(UDE)==paste(Year[i],M1[i], M1[i], sep="*"))
              UDE[i,j]=UDE[i,j]+3/16

          
              j= which(colnames(UDE)==paste(Year[i],F2[i], F2[i], sep="*"))
              UDE[i,j]=UDE[i,j]+3/16

          
              j= which(colnames(UDE)==paste(Year[i],M2[i], M2[i], sep="*"))
              UDE[i,j]=UDE[i,j]+3/16


              j= which(colnames(UDE)==paste(Year[i],F1[i], F2[i], sep="*"))
              UDE[i,j]=UDE[i,j]+1/16
          
              j= which(colnames(UDE)==paste(Year[i],F1[i], M2[i], sep="*"))
              UDE[i,j]=UDE[i,j]+1/16
          
              j= which(colnames(UDE)==paste(Year[i],M1[i], F2[i], sep="*"))
              UDE[i,j]=UDE[i,j]+1/16
          
              j= which(colnames(UDE)==paste(Year[i],M1[i], M2[i], sep="*"))
              UDE[i,j]=UDE[i,j]+1/16
           }  
       }
       colnames(UDE)=paste("DE(",TDE,")",sep="")
       UDE= DropColumns(UDE)
   }

   if(BlockID==1){
       if(y0==1)UB=IndMatrix(Block)
       else if(y0>1){
           B=paste(Env,Block,sep="*")
           UB=IndMatrix(Block)
       }
       colnames(UB)=paste("B(",colnames(UB),")",sep="")
   }
   
   if(y0==1&&BlockID==0){
        U=list(UC,UA,UD)
        names(U)=c("C","A","D")
   }
   else if(y0==1&&BlockID==1){
        U=list(UC,UA,UD,UB)
        names(U)=c("C","A","D","B")
   }
   else if(y0>1&&BlockID==0){
        U=list(UY,UC,UA,UD,UCE,UAE,UDE)
        names(U)=c("E","C","A","D","CE","AE","DE")
   }
   else if(y0>1&&BlockID==1){
        U=list(UY,UC,UA,UD,UC,UAE,UDE,UB)
        names(U)=c("Y","C","A","D","CE","AE","DE","B")
   }
   class(U)="Umatrix"
   return(U)

}




#### Implementation of design matrices for 4-way Marker AD model.

#######################################################################################
#######################################################################################
#######################################################################################

MarkerAD4Matrix=function(Env,MP,F1,M1,F2,M2,Gen,Block,BlockID=NULL){
   Year=Env
   TP=sort(unique(c(F1,M1,F2,M2)))
   Ty=sort(unique(Year))
   Tb=sort(unique(Block))
   Tm=sort(unique(MP))
      
   n=length(F1)

   y0=length(Ty)
   p0=length(TP)
   b0=length(Tb)
   if(is.null(BlockID))BlockID=0
   if(BlockID==1&&b0==1)BlockID=0
   
   ##########################################
   ## Marker Additive effects#######################
   UAm=matrix(0,nrow=n,ncol=length(Tm))
   #dim(UA)
   colnames(UAm)=Tm
   for(i in 1:n){
       j=which(colnames(UAm)==MP[F1[i]])
       UAm[i,j]=UAm[i,j]+0.5
       j=which(colnames(UAm)==MP[M1[i]])
       UAm[i,j]=UAm[i,j]+0.5
       j=which(colnames(UAm)==MP[F2[i]])
       UAm[i,j]=UAm[i,j]+0.5
       j=which(colnames(UAm)==MP[M2[i]])
       UAm[i,j]=UAm[i,j]+0.5
   }
   colnames(UAm)=paste("Am(",Tm,")",sep="")
   UAm=DropColumns(UAm)


  ##################################################
   ## Dominance effects  ############################
   TDm=NULL

   TDm=paste(Tm,Tm,sep="*")

   for(i in 1:(length(Tm)-1)){
       str=paste(Tm[i],Tm[(i+1):length(Tm)],sep="*")
       TDm=c(TDm,str)
   }

   UDm=matrix(0,nrow=n,ncol=length(TDm))
   
   colnames(UDm)=TDm

   for(i in 1:n){
      if(Gen[i]==0){
          j= which(colnames(UDm)==paste(MP[F1[i]], MP[F1[i]], sep="*"))
          UDm[i,j]=1
      }
      
      else if(Gen[i]==1){
          j= which(colnames(UDm)==paste(MP[F1[i]], MP[F2[i]], sep="*"))
          UDm[i,j]=UDm[i,j]+0.25
          
          j= which(colnames(UDm)==paste(MP[F1[i]], MP[M2[i]], sep="*"))
          UDm[i,j]=UDm[i,j]+0.25
          
          j= which(colnames(UDm)==paste(MP[M1[i]], MP[F2[i]], sep="*"))
          UDm[i,j]=UDm[i,j]+0.25
          
          j= which(colnames(UDm)==paste(MP[M1[i]], MP[M2[i]], sep="*"))
          UDm[i,j]=UDm[i,j]+0.25
      }

      else if(Gen[i]==2){
          j= which(colnames(UDm)==paste(MP[F1[i]], MP[F1[i]], sep="*"))
          UDm[i,j]=UDm[i,j]+0.25/2
          
          j= which(colnames(UDm)==paste(MP[M1[i]], MP[M1[i]], sep="*"))
          UDm[i,j]=UDm[i,j]+0.25/2
          
          j= which(colnames(UDm)==paste(MP[F2[i]], MP[F2[i]], sep="*"))
          UDm[i,j]=UDm[i,j]+0.25/2
          
          j= which(colnames(UDm)==paste(MP[M2[i]], MP[M2[i]], sep="*"))
          UDm[i,j]=UDm[i,j]+0.25/2

          j= which(colnames(UDm)==paste(MP[F1[i]], MP[F2[i]], sep="*"))
          UDm[i,j]=UDm[i,j]+0.25/2
          
          j= which(colnames(UDm)==paste(MP[F1[i]], MP[M2[i]], sep="*"))
          UDm[i,j]=UDm[i,j]+0.25/2
          
          j= which(colnames(UDm)==paste(MP[M1[i]], MP[F2[i]], sep="*"))
          UDm[i,j]=UDm[i,j]+0.25/2
          
          j= which(colnames(UDm)==paste(MP[M1[i]], MP[M2[i]], sep="*"))
          UDm[i,j]=UDm[i,j]+0.25/2
      }
      
      else if(Gen[i]==3){
          j= which(colnames(UDm)==paste(MP[F1[i]], MP[F1[i]], sep="*"))
          UDm[i,j]=UDm[i,j]+3/16
          
          j= which(colnames(UDm)==paste(MP[M1[i]], MP[M1[i]], sep="*"))
          UDm[i,j]=UDm[i,j]+3/16
          
          j= which(colnames(UDm)==paste(MP[F2[i]], MP[F2[i]], sep="*"))
          UDm[i,j]=UDm[i,j]+3/16
          
          j= which(colnames(UDm)==paste(MP[M2[i]], MP[M2[i]], sep="*"))
          UDm[i,j]=UDm[i,j]+3/16

          j= which(colnames(UDm)==paste(MP[F1[i]], MP[F2[i]], sep="*"))
          UDm[i,j]=UDm[i,j]+1/16
          
          j= which(colnames(UDm)==paste(MP[F1[i]], MP[M2[i]], sep="*"))
          UDm[i,j]=UDm[i,j]+1/16
          
          j= which(colnames(UDm)==paste(MP[M1[i]], MP[F2[i]], sep="*"))
          UDm[i,j]=UDm[i,j]+1/16
          
          j= which(colnames(UDm)==paste(MP[M1[i]], MP[M2[i]], sep="*"))
          UDm[i,j]=UDm[i,j]+1/16
      }
   }
   colnames(UDm)=paste("Dm(",TDm,")",sep="")
   UDm= DropColumns(UDm)

   
   #P=c(1:10)
   #Y=c(1,2)
   #TYP=NULL
   #for(i in 1:2){
   #  TYP=c(TYP,paste(Y[i],P,sep="*"))
   #}
   
   ##########################################
   ## Additive effects#######################
   UA=matrix(0,nrow=n,ncol=p0)
   dim(UA)
   colnames(UA)=TP
   for(i in 1:n){
       if(Gen[i]==0){
           j=which(colnames(UA)==F1[i])
           UA[i,j]=2
       }
       else if(Gen[i]>0){
           j=which(colnames(UA)==F1[i])
           UA[i,j]=UA[i,j]+0.5
           j=which(colnames(UA)==M1[i])
           UA[i,j]=UA[i,j]+0.5
           j=which(colnames(UA)==F2[i])
           UA[i,j]=UA[i,j]+0.5
           j=which(colnames(UA)==M2[i])
           UA[i,j]=UA[i,j]+0.5
       }
   }
   colnames(UA)=paste("A(",TP,")",sep="")
   UA=DropColumns(UA)
   
   ##################################################
   ## Dominance effects  ############################
   TD=NULL

   TD=paste(TP,TP,sep="*")

   for(i in 1:(p0-1)){
       str=paste(TP[i],TP[(i+1):p0],sep="*")
       TD=c(TD,str)
   }

   UD=matrix(0,nrow=n,ncol=length(TD))
   
   colnames(UD)=TD

   for(i in 1:n){
      if(Gen[i]==0){
          j= which(colnames(UD)==paste(F1[i], F1[i], sep="*"))
          UD[i,j]=1
      }
      
      else if(Gen[i]==1){
          j= which(colnames(UD)==paste(F1[i], F2[i], sep="*"))
          UD[i,j]=UD[i,j]+0.25
          
          j= which(colnames(UD)==paste(F1[i], M2[i], sep="*"))
          UD[i,j]=UD[i,j]+0.25
          
          j= which(colnames(UD)==paste(M1[i], F2[i], sep="*"))
          UD[i,j]=UD[i,j]+0.25
          
          j= which(colnames(UD)==paste(M1[i], M2[i], sep="*"))
          UD[i,j]=UD[i,j]+0.25
      }

      else if(Gen[i]==2){
          j= which(colnames(UD)==paste(F1[i], F1[i], sep="*"))
          UD[i,j]=UD[i,j]+0.25/2
          
          j= which(colnames(UD)==paste(M1[i], M1[i], sep="*"))
          UD[i,j]=UD[i,j]+0.25/2
          
          j= which(colnames(UD)==paste(F2[i], F2[i], sep="*"))
          UD[i,j]=UD[i,j]+0.25/2
          
          j= which(colnames(UD)==paste(M2[i], M2[i], sep="*"))
          UD[i,j]=UD[i,j]+0.25/2

          j= which(colnames(UD)==paste(F1[i], F2[i], sep="*"))
          UD[i,j]=UD[i,j]+0.25/2
          
          j= which(colnames(UD)==paste(F1[i], M2[i], sep="*"))
          UD[i,j]=UD[i,j]+0.25/2
          
          j= which(colnames(UD)==paste(M1[i], F2[i], sep="*"))
          UD[i,j]=UD[i,j]+0.25/2
          
          j= which(colnames(UD)==paste(M1[i], M2[i], sep="*"))
          UD[i,j]=UD[i,j]+0.25/2
      }
      
      else if(Gen[i]==3){
          j= which(colnames(UD)==paste(F1[i], F1[i], sep="*"))
          UD[i,j]=UD[i,j]+3/16
          
          j= which(colnames(UD)==paste(M1[i], M1[i], sep="*"))
          UD[i,j]=UD[i,j]+3/16
          
          j= which(colnames(UD)==paste(F2[i], F2[i], sep="*"))
          UD[i,j]=UD[i,j]+3/16
          
          j= which(colnames(UD)==paste(M2[i], M2[i], sep="*"))
          UD[i,j]=UD[i,j]+3/16

          j= which(colnames(UD)==paste(F1[i], F2[i], sep="*"))
          UD[i,j]=UD[i,j]+1/16
          
          j= which(colnames(UD)==paste(F1[i], M2[i], sep="*"))
          UD[i,j]=UD[i,j]+1/16
          
          j= which(colnames(UD)==paste(M1[i], F2[i], sep="*"))
          UD[i,j]=UD[i,j]+1/16
          
          j= which(colnames(UD)==paste(M1[i], M2[i], sep="*"))
          UD[i,j]=UD[i,j]+1/16 
      }
   }
   colnames(UD)=paste("D(",TD,")",sep="")
   UD= DropColumns(UD)
   

   ##############################################################
   ### Environmental effects:  random effect#####################
   if(y0>0){
       ### Environmental effects ###################################
       UY=matrix(0,nrow=n,ncol=length(Ty))
       colnames(UY)= Ty
       for(i in 1:n){
	     j=which(colnames(UY)==Year[i])
	     UY[i,j] = 1
       }
       colnames(UY)=paste("E(",Ty,")",sep="")
       #UY
       #dim(UY)
       UY= DropColumns(UY)
   
       ##Marker additive*environment interaction effects
       
       ## comments: AmE interaction
       TAmE=NULL
       for(i in 1:y0){
           str=paste(Ty[i],Tm,sep="*")
           TAmE=c(TAmE,str)
       }
     
       UAmE=matrix(0,nrow=n,ncol=length(TAmE))
       #UAE
       colnames(UAmE)=TAmE
       #head(UAE)
       for(i in 1:n){
	     j= which(colnames(UAmE)==paste(Year[i], MP[F1[i]], sep="*"))
	     UAmE[i, j] = UAmE[i,j]+0.5

           j= which(colnames(UAmE)==paste(Year[i], MP[M1[i]], sep="*"))
	     UAmE[i, j] = UAmE[i,j]+0.5
           
           j= which(colnames(UAmE)==paste(Year[i], MP[F2[i]], sep="*"))
	     UAmE[i, j] = UAmE[i,j]+0.5

           j= which(colnames(UAmE)==paste(Year[i], MP[M2[i]], sep="*"))
	     UAmE[i, j] = UAmE[i,j]+0.5
       }
       colnames(UmAE)=paste("AmE(",TAmE,")",sep="")
       #dim(UAE)

       UAmE= DropColumns(UAmE)
       #dim(UAE)

       ###########################################################
       ## DE interaction effects #################################
       
       TDmE=NULL

       for(i in 1:y0){
           for(j in 1:length(TDm)){
              TDmE=c(TDmE,paste(Ty[i],TDm[j],sep="*"))
           }
       }

       UDmE=matrix(0,nrow=n,ncol=length(TDmE))
       colnames(UDmE)=TDmE
          
       for(i in 1:n){
           if(Gen[i]==0){
               j= which(colnames(UDmE)==paste(Year[i],MP[F1[i]], MP[F1[i]], sep="*"))
               UDmE[i,j]=1

           }

           else if(Gen[i]==1){
              j= which(colnames(UDmE)==paste(Year[i],MP[F1[i]], MP[F2[i]], sep="*"))
              UDmE[i,j]=UDmE[i,j]+0.25
          
              j= which(colnames(UDmE)==paste(Year[i],MP[F1[i]], MP[M2[i]], sep="*"))
              UDmE[i,j]=UDmE[i,j]+0.25
          
              j= which(colnames(UDmE)==paste(Year[i],MP[M1[i]], MP[F2[i]], sep="*"))
              UDmE[i,j]=UDmE[i,j]+0.25
          
              j= which(colnames(UDmE)==paste(Year[i],MP[M2[i]], MP[M2[i]], sep="*"))
              UDmE[i,j]=UDmE[i,j]+0.25
           }
             
           else if(Gen[i]==2){
              j= which(colnames(UDmE)==paste(Year[i],MP[F1[i]], MP[F1[i]], sep="*"))
              UDmE[i,j]=UDmE[i,j]+0.25/2
          
              j= which(colnames(UDmE)==paste(Year[i],MP[M1[i]], MP[M1[i]], sep="*"))
              UDmE[i,j]=UDmE[i,j]+0.25/2
          
              j= which(colnames(UDmE)==paste(Year[i],MP[F2[i]], MP[F2[i]], sep="*"))
              UDmE[i,j]=UDmE[i,j]+0.25/2
          
              j= which(colnames(UDmE)==paste(Year[i],MP[M2[i]], MP[M2[i]], sep="*"))
              UDmE[i,j]=UDmE[i,j]+0.25/2

              j= which(colnames(UDmE)==paste(Year[i],MP[F1[i]], MP[F2[i]], sep="*"))
              UDmE[i,j]=UDmE[i,j]+0.25/2
          
              j= which(colnames(UDmE)==paste(Year[i],MP[F1[i]], MP[M2[i]], sep="*"))
              UDE[i,j]=UDE[i,j]+0.25/2
          
              j= which(colnames(UDmE)==paste(Year[i],MP[M1[i]], MP[F2[i]], sep="*"))
              UDmE[i,j]=UDmE[i,j]+0.25/2
          
              j= which(colnames(UDmE)==paste(Year[i],M1[i], M2[i], sep="*"))
              UDmE[i,j]=UDmE[i,j]+0.25/2

           }

           else if(Gen[i]==3){

              j= which(colnames(UDmE)==paste(Year[i],MP[F1[i]], MP[F1[i]], sep="*"))
              UDmE[i,j]=UDmE[i,j]+3/16
          
              j= which(colnames(UDmE)==paste(Year[i],MP[M1[i]], MP[M1[i]], sep="*"))
              UDmE[i,j]=UDmE[i,j]+3/16
          
              j= which(colnames(UDmE)==paste(Year[i],MP[F2[i]], MP[F2[i]], sep="*"))
              UDmE[i,j]=UDmE[i,j]+3/16
          
              j= which(colnames(UDmE)==paste(Year[i],MP[M2[i]], MP[M2[i]], sep="*"))
              UDmE[i,j]=UDmE[i,j]+3/16

              j= which(colnames(UDmE)==paste(Year[i],MP[F1[i]], MP[F2[i]], sep="*"))
              UDmE[i,j]=UDmE[i,j]+1/16
          
              j= which(colnames(UDmE)==paste(Year[i],MP[F1[i]], MP[M2[i]], sep="*"))
              UDmE[i,j]=UDmE[i,j]+1/16
          
              j= which(colnames(UDmE)==paste(Year[i],MP[M1[i]], MP[F2[i]], sep="*"))
              UDmE[i,j]=UDmE[i,j]+1/16
          
              j= which(colnames(UDmE)==paste(Year[i],M1[i], M2[i], sep="*"))
              UDmE[i,j]=UDmE[i,j]+1/16


           }  
       }
       colnames(UDmE)=paste("DmE(",TDmE,")",sep="")
       UDmE= DropColumns(UDmE)


       
       ##Additive*environment interaction effects
       
       ## comments: AE interaction
       TAE=NULL
       for(i in 1:y0){
           str=paste(Ty[i],TP[1:p0],sep="*")
           TAE=c(TAE,str)
       }
     
       UAE=matrix(0,nrow=n,ncol=length(TAE))
       #UAE
       colnames(UAE)=TAE
       #head(UAE)
       for(i in 1:n){
	     j= which(colnames(UAE)==paste(Year[i], F1[i], sep="*"))
	     UAE[i, j] = UAE[i,j]+0.5

           j= which(colnames(UAE)==paste(Year[i], M1[i], sep="*"))
	     UAE[i, j] = UAE[i,j]+0.5
           
           j= which(colnames(UAE)==paste(Year[i], F2[i], sep="*"))
	     UAE[i, j] = UAE[i,j]+0.5

           j= which(colnames(UAE)==paste(Year[i], M2[i], sep="*"))
	     UAE[i, j] = UAE[i,j]+0.5
       }
       colnames(UAE)=paste("AE(",TAE,")",sep="")
       #dim(UAE)

       UAE= DropColumns(UAE)
       #dim(UAE)

       ###########################################################
       ## DE interaction effects #################################
       
       TDE=NULL

       for(i in 1:y0){
           for(j in 1:length(TD)){
              TDE=c(TDE,paste(Ty[i],TD[j],sep="*"))
           }
       }

       UDE=matrix(0,nrow=n,ncol=length(TDE))
       colnames(UDE)=TDE
          
       for(i in 1:n){
           if(Gen[i]==0){
               j= which(colnames(UDE)==paste(Year[i],F1[i], F1[i], sep="*"))
               UDE[i,j]=1

           }

           else if(Gen[i]==1){
              j= which(colnames(UDE)==paste(Year[i],F1[i], F2[i], sep="*"))
              UDE[i,j]=UDE[i,j]+0.25
          
              j= which(colnames(UDE)==paste(Year[i],F1[i], M2[i], sep="*"))
              UDE[i,j]=UDE[i,j]+0.25
          
              j= which(colnames(UDE)==paste(Year[i],M1[i], F2[i], sep="*"))
              UDE[i,j]=UDE[i,j]+0.25
          
              j= which(colnames(UDE)==paste(Year[i],M2[i], M2[i], sep="*"))
              UDE[i,j]=UDE[i,j]+0.25
           }
             
           else if(Gen[i]==2){
              j= which(colnames(UDE)==paste(Year[i],F1[i], F1[i], sep="*"))
              UDE[i,j]=UDE[i,j]+0.25/2
          
              j= which(colnames(UDE)==paste(Year[i],M1[i], M1[i], sep="*"))
              UDE[i,j]=UDE[i,j]+0.25/2
          
              j= which(colnames(UDE)==paste(Year[i],F2[i], F2[i], sep="*"))
              UDE[i,j]=UDE[i,j]+0.25/2
          
              j= which(colnames(UDE)==paste(Year[i],M2[i], M2[i], sep="*"))
              UDE[i,j]=UDE[i,j]+0.25/2

              j= which(colnames(UDE)==paste(Year[i],F1[i], F2[i], sep="*"))
              UDE[i,j]=UDE[i,j]+0.25/2
          
              j= which(colnames(UDE)==paste(Year[i],F1[i], M2[i], sep="*"))
              UDE[i,j]=UDE[i,j]+0.25/2
          
              j= which(colnames(UDE)==paste(Year[i],M1[i], F2[i], sep="*"))
              UDE[i,j]=UDE[i,j]+0.25/2
          
              j= which(colnames(UDE)==paste(Year[i],M1[i], M2[i], sep="*"))
              UDE[i,j]=UDE[i,j]+0.25/2

           }

           else if(Gen[i]==3){

              j= which(colnames(UDE)==paste(Year[i],F1[i], F1[i], sep="*"))
              UDE[i,j]=UDE[i,j]+3/16
          
              j= which(colnames(UDE)==paste(Year[i],M1[i], M1[i], sep="*"))
              UDE[i,j]=UDE[i,j]+3/16

          
              j= which(colnames(UDE)==paste(Year[i],F2[i], F2[i], sep="*"))
              UDE[i,j]=UDE[i,j]+3/16

          
              j= which(colnames(UDE)==paste(Year[i],M2[i], M2[i], sep="*"))
              UDE[i,j]=UDE[i,j]+3/16


              j= which(colnames(UDE)==paste(Year[i],F1[i], F2[i], sep="*"))
              UDE[i,j]=UDE[i,j]+1/16
          
              j= which(colnames(UDE)==paste(Year[i],F1[i], M2[i], sep="*"))
              UDE[i,j]=UDE[i,j]+1/16
          
              j= which(colnames(UDE)==paste(Year[i],M1[i], F2[i], sep="*"))
              UDE[i,j]=UDE[i,j]+1/16
          
              j= which(colnames(UDE)==paste(Year[i],M1[i], M2[i], sep="*"))
              UDE[i,j]=UDE[i,j]+1/16
           }  
       }
       colnames(UDE)=paste("DE(",TDE,")",sep="")
       UDE= DropColumns(UDE)
   }

   if(BlockID==1){
       if(y0==1&&b0>1){
           UB=matrix(0,nrow=n,ncol=b0)
           colnames(UB)=Tb
           for(i in 1:n){
                j=which(colnames(UB)==Block[i])
                UB[i,j]=1
           }
           colnames(UB)=paste("B(",Tb,")",sep="")
           #dim(UB)

           UB= DropColumns(UB)
           #dim(UB)
       }
       else if(y0>1&&b0>1){
           TB=NULL

           for(i in 1:y0){
               str=paste(Ty[i],Tb[1:b0],sep="*")
               TB=c(TB,str)
  
           }

           UB=matrix(0,nrow=n,ncol=length(TB))
           #UB
           colnames(UB)=TB
         
           for(i in 1:n){
	
	          j= which(colnames(UB)==paste(Year[i],Block[i],sep="*"))
	
	          UB[i, j] = 1
           }
           colnames(UB)=paste("B(",TB,")",sep="")
           
           UB= DropColumns(UB)
           
        }
    }
    if(y0==1&&BlockID==0){
        U=list(UAm,UDm,UA,UD)
        names(U)=c("Am","Dm","A","D")
    }
    else if(y0==1&&BlockID==1){
        U=list(UAm,UDm,UA,UD,UB)
        names(U)=c("Am","Dm","A","D","B")
    }
    else if(y0>1&&BlockID==0){
       U=list(UY,UAm,UDm,UA,UD,UAmE,UDmE,UAE,UDE)
       names(U)=c("E","Am","Dm","A","D","AmE","DmE","AE","DE")
    }

    else if(y0>1&&BlockID==1){
      U=list(UY,UAm,UDm,UA,UD,UAmE,UDmE,UAE,UDE,UB)
      names(U)=c("E","Am","Dm","A","D","AmE","DmE","AE","DE","B")
    }
    class(U)="Umatrix"
    return(U)
}
             
       

#### Implementation of design matrices for 2-way ADM model.

#######################################################################################
#######################################################################################
#######################################################################################
ADM2Matrix=function(Env,F1,M1,Gen,Block,BlockID=NULL){
   Year=Env
   TP=sort(unique(c(F1,M1)))
   Ty=sort(unique(Year))
   Tb=sort(unique(Block))
   
   n=length(F1)

   y0=length(Ty)
   p0=length(TP)
   b0=length(Tb)
   if(is.null(Block))BlockID=0
   if(BlockID==1&&b0==1)BlockID=0
   GenM=Gen
   id=which(GenM>0)
   GenM[id]=GenM[id]-1
   M2=M1
   id=which(GenM==0)
   M2[id]=F1[id]
   UAm=UA2Matrix(F1,M2)
   colnames(UAm)=gsub("A","Am",colnames(UAm))
   UA=UA2Matrix(F1,M1)
   UDm=UD2Matrix(F1,M2,GenM)
   colnames(UDm)=gsub("D","Dm",colnames(UDm))
   UD=UD2Matrix(F1,M1,Gen)
   if(y0==1){
      if(BlockID==1)UB=IndMatrix(Block)
      colnames(UB)=paste("B(",colnames(UB),")",sep="")

   }
   if(y0>1){
     UY=IndMatrix(Env)
     colnames(UY)=paste("E(",colnames(UY),")",sep="")
     UAmE=UAE2Matrix(Env,F1,M2)
     colnames(UAmE)=gsub("A","Am",colnames(UAmE))
     UDmE=UDE2Matrix(Env,F1,M2,GenM)
     colnames(UDmE)=gsub("D","Dm",colnames(UDmE))
     UAE=UAE2Matrix(Env,F1,M1)
     UDE=UDE2Matrix(Env,F1,M1,Gen)
     if(BlockID==1){
       B=paste(Env,Block,sep="*")
       UB=IndMatrix(B)
       colnames(UB)=paste("B(",colnames(UB),")",sep="")
     }
   }

   
   if(y0==1&&BlockID==0){
       U=list(UAm,UDm,UA,UD)
       names(U)=c("Am","Dm","A","D")
   }
   else if(y0==1&&BlockID==1){
       U=list(UAm,UDm,UA,UD,UB)
       names(U)=c("Am","Dm","A","D","B")
   }

   else if(y0>1&&BlockID==0){
       U=list(UY,UAm,UDm,UA,UD,UAmE,UDmE,UAE,UDE)
       names(U)=c("E","Am","Dm","A","D","AmE","DmE","AE","DE")
   }

   else if(y0>1&&BlockID==1){
       U=list(UY,UAm,UDm,UA,UD,UAmE,UDmE,UAE,UDE,UB)
       names(U)=c("E","Am","Dm","A","D","AmE","DmE","AE","DE","B")
   }
   class(U)="Umatrix"
    
   return(U)

}

#### Implementation of design matrices for 2-way AD model.

#######################################################################################
#######################################################################################
#######################################################################################
#Env=dat$Loc
#F1=dat$Female
#M1=dat$Male
#Block=dat$Blk
#BlockID=1
#Gen=dat$Gen

AD2Matrix=function(Env,F1,M1,Gen,Block=NULL,BlockID=NULL){
   Year=Env
   F1=as.vector(F1)
   M1=as.vector(M1)
   TP=sort(unique(c(F1,M1)))
   Ty=sort(unique(Year))
   if(is.null(Block))Block=1
   Tb=sort(unique(Block))
   
   n=length(F1)

   y0=length(Ty)
   p0=length(TP)
   b0=length(Tb)
   if(is.null(BlockID))BlockID=0
   if(BlockID==1&&b0==1)BlockID=0

   UA=UA2Matrix(F1,M1)
   UD=UD2Matrix(F1,M1,Gen)
   #UAA=UAA2Matrix(F1,M1)
   if(y0==1){
      if(BlockID==1){
         UB=IndMatrix(Block)
         colnames(UB)=paste("B(",colnames(UB),")",sep="")
      }

   }
   if(y0>1){
     UY=IndMatrix(Env)
     colnames(UY)=paste("E(",colnames(UY),")",sep="")
     UAE=UAE2Matrix(Env,F1,M1)
     UDE=UDE2Matrix(Env,F1,M1,Gen)
     #UAAE=UAAE2Matrix(Env,F1,M1)
     if(BlockID==1){
       B=paste(Env,Block,sep="*")
       UB=IndMatrix(B)
       colnames(UB)=paste("B(",colnames(UB),")",sep="")
     }
   }

    if(y0==1&&BlockID==0){
       U=list(UA,UD)
       names(U)=c("A","D")
    }
    else if(y0==1&&BlockID==1){
       U=list(UA,UD,UB)
       names(U)=c("A","D","B")
    }

    else if(y0>1&&BlockID==0){
       U=list(UY,UA,UD,UAE,UDE)
       names(U)=c("E","A","D","AE","DE")
    }
    else if(y0>1&&BlockID==1){
       U=list(UY,UA,UD,UAE,UDE,UB)
       names(U)=c("E","A","D","AE","DE","B")
    }
    class(U)="Umatrix"
    return(U)

}

#### Implementation of design matrices for 2-way RCAD model.

#######################################################################################
#######################################################################################
#######################################################################################

RCAD2Matrix=function(Env,R,C,F1,M1,Gen,Block,BlockID){
   Year=Env
   TP=sort(unique(c(F1,M1)))
   Ty=sort(unique(Year))
   Tb=sort(unique(Block))
   Tr=sort(unique(R))
   Tc=sort(unique(C))
   
   n=length(F1)

   y0=length(Ty)
   p0=length(TP)
   b0=length(Tb)
   c0=length(Tc)
   r0=length(Tr)

   if(BlockID==1&&b0==1)BlockID=0

   ##########################################
   ## Additive effects#######################
   UA=matrix(0,nrow=n,ncol=p0)
   dim(UA)
   colnames(UA)=TP
   for(i in 1:n){
       if(Gen[i]==0){
           j=which(colnames(UA)==F1[i])
           UA[i,j]=2
       }
       else if(Gen[i]>0){
           j=which(colnames(UA)==F1[i])
           UA[i,j]=UA[i,j]+1
           j=which(colnames(UA)==M1[i])
           UA[i,j]=UA[i,j]+1
     
       }
   }
   colnames(UA)=paste("A(",TP,")",sep="")
   UA=DropColumns(UA)
   
   ##################################################
   ## Dominance effects  ############################
   TD=NULL

   TD=paste(TP,TP,sep="*")

   for(i in 1:(p0-1)){
       str=paste(TP[i],TP[(i+1):p0],sep="*")
       TD=c(TD,str)
   }

   UD=matrix(0,nrow=n,ncol=length(TD))
   
   colnames(UD)=TD

   for(i in 1:n){
      if(Gen[i]==0){
          j= which(colnames(UD)==paste(F1[i], M1[i], sep="*"))
          UD[i,j]=1
      }
      
      else if(Gen[i]==1){
          j= which(colnames(UD)==paste(F1[i], M1[i], sep="*"))
          UD[i,j]=UD[i,j]+1
     }

      else if(Gen[i]==2){
          j= which(colnames(UD)==paste(F1[i], F1[i], sep="*"))
          UD[i,j]=UD[i,j]+0.25
          
          j= which(colnames(UD)==paste(M1[i], M1[i], sep="*"))
          UD[i,j]=UD[i,j]+0.25
          
          j= which(colnames(UD)==paste(F1[i], M1[i], sep="*"))
          UD[i,j]=UD[i,j]+0.5
          
      }
      
      else if(Gen[i]==3){
          j= which(colnames(UD)==paste(F1[i], F1[i], sep="*"))
          UD[i,j]=UD[i,j]+3/8
          
          j= which(colnames(UD)==paste(M1[i], M1[i], sep="*"))
          UD[i,j]=UD[i,j]+3/8
          
          j= which(colnames(UD)==paste(F1[i], M1[i], sep="*"))
          UD[i,j]=UD[i,j]+0.25
          
      }
   }
   colnames(UD)=paste("D(",TD,")",sep="")
   UD= DropColumns(UD)
   

   ##############################################################
   ### Environmental effects:  random effect#####################
   if(y0>0){
       ### Environmental effects ###################################
       UY=matrix(0,nrow=n,ncol=length(Ty))
       colnames(UY)= Ty
       for(i in 1:n){
	     j=which(colnames(UY)==Year[i])
	     UY[i,j] = 1
       }
       colnames(UY)=paste("E(",Ty,")",sep="")
       #UY
       #dim(UY)
       UY= DropColumns(UY)
   
       
       ##Additive*environment interaction effects
       
       ## comments: AE interaction
       TAE=NULL
       for(i in 1:y0){
           str=paste(Ty[i],TP[1:p0],sep="*")
           TAE=c(TAE,str)
       }
     
       UAE=matrix(0,nrow=n,ncol=length(TAE))
       UAE
       colnames(UAE)=TAE
       head(UAE)
       for(i in 1:n){
	     j= which(colnames(UAE)==paste(Year[i], F1[i], sep="*"))
	     UAE[i, j] = UAE[i,j]+1

           j= which(colnames(UAE)==paste(Year[i], M1[i], sep="*"))
	     UAE[i, j] = UAE[i,j]+1
           
       }
       colnames(UAE)=paste("AE(",TAE,")",sep="")
       #dim(UAE)

       UAE= DropColumns(UAE)
       #dim(UAE)

       ###########################################################
       ## DE interaction effects #################################
       
       TDE=NULL

       for(i in 1:y0){
           for(j in 1:length(TD)){
              TDE=c(TDE,paste(Ty[i],TD[j],sep="*"))
           }
       }

       UDE=matrix(0,nrow=n,ncol=length(TDE))
       colnames(UDE)=TDE
          
       for(i in 1:n){
           if(Gen[i]==0){
               j= which(colnames(UDE)==paste(Year[i],F1[i], F1[i], sep="*"))
               UDE[i,j]=1

           }

           else if(Gen[i]==1){
              j= which(colnames(UDE)==paste(Year[i],F1[i], M1[i], sep="*"))
              UDE[i,j]=UDE[i,j]+1
           }
             
           else if(Gen[i]==2){
              j= which(colnames(UDE)==paste(Year[i],F1[i], F1[i], sep="*"))
              UDE[i,j]=UDE[i,j]+0.25
          
              j= which(colnames(UDE)==paste(Year[i],M1[i], M1[i], sep="*"))
              UDE[i,j]=UDE[i,j]+0.25
          
              j= which(colnames(UDE)==paste(Year[i],F1[i], M1[i], sep="*"))
              UDE[i,j]=UDE[i,j]+0.5
           }

           else if(Gen[i]==3){

              j= which(colnames(UDE)==paste(Year[i],F1[i], F1[i], sep="*"))
              UDE[i,j]=UDE[i,j]+3/8
          
              j= which(colnames(UDE)==paste(Year[i],M1[i], M1[i], sep="*"))
              UDE[i,j]=UDE[i,j]+3/8

          
              j= which(colnames(UDE)==paste(Year[i],F1[i], M1[i], sep="*"))
              UDE[i,j]=UDE[i,j]+1/4
           }  
       }
       colnames(UDE)=paste("DE(",TDE,")",sep="")
       UDE= DropColumns(UDE)
   }

   if(BlockID==1){
       if(y0==1&&b0>1){
           UB=matrix(0,nrow=n,ncol=b0)
           colnames(UB)=Tb
           for(i in 1:n){
                j=which(colnames(UB)==Block[i])
                UB[i,j]=1
           }
           colnames(UB)=paste("B(",Tb,")",sep="")
           #dim(UB)

           UB= DropColumns(UB)
           #dim(UB)
       }
       else if(y0>1&&b0>1){
           TB=NULL

           for(i in 1:y0){
               str=paste(Ty[i],Tb[1:b0],sep="*")
               TB=c(TB,str)
  
           }

           UB=matrix(0,nrow=n,ncol=length(TB))
           #UB
           colnames(UB)=TB
         
           for(i in 1:n){
	
	          j= which(colnames(UB)==paste(Year[i],Block[i],sep="*"))
	
	          UB[i, j] = 1
           }
           colnames(UB)=paste("B(",TB,")",sep="")
           
           UB= DropColumns(UB)
           
        }
   }

    
    if(y0==1){
       UR=matrix(0,nrow=n,ncol=r0)
       colnames(UR)=Tr
       UC=matrix(0,nrow=n,ncol=c0)
       colnames(UC)=Tc
       for(i in 1:n){
          j=which(colnames(UR)==R[i])
          UR[i,j]=1
          j=which(colnames(UC)==C[i])
          UC[i,j]=1
       }
       colnames(UR)=paste("Row(",Tr,")",sep="")
       colnames(UC)=paste("Col(",Tc,")",sep="")
           #dim(UB)

       UC=DropColumns(UC)
       UR=DropColumns(UR)

           #dim(UB)
    }

    else if(y0>1){
      TR=NULL
      TC=NULL
       
      for(i in 1:y0){
         TR=c(TR,paste(Ty[i],Tr,sep="*"))
         TC=c(TC,paste(Ty[i],Tc,sep="*"))
      }
      
      UR=matrix(0,nrow=n,ncol=length(TR))
      UC=matrix(0,nrow=n,ncol=length(TC))
      colnames(UR)=TR
      colnames(UC)=TC
      
      for(i in 1:n){
          j=which(colnames(UR)==paste(Year[i],R[i],sep="*"))
          UR[i,j]=1
          j=which(colnames(UC)==paste(Year[i],C[i],sep="*"))
          UC[i,j]=1
       }
       
       colnames(UR)=paste("Row(",TR,")",sep="")
       colnames(UC)=paste("Col(",TR,")",sep="")
       
       UC=DropColumns(UC)
       UR=DropColumns(UR)
    }

    if(y0==1&&BlockID==0){
        U=list(UR,UC,UA,UD)
        names(U)=c("R","C","A","D")
    }
    else if(y0==1&&BlockID==1){
        U=list(UR,UC,UA,UD,UB)
        names(U)=c("R","C","A","D","B")
    }
    else if(y0>1&&BlockID==0){
        U=list(UY,UR,UC,UA,UD,UAE,UDE)
        names(U)=c("E","R","C","A","D","AE","DE")
    }

    else if(y0>1&&BlockID==1){
        U=list(UY,UR,UC,UA,UD,UAE,UDE,UB)
        names(U)=c("E","R","C","A","D","AE","DE","B")
    }
    class(U)="Umatrix"
    return(U)

}
   

#### Implementation of design matrices for 2-way ADC model.

#######################################################################################
#######################################################################################
#######################################################################################

ADC2Matrix=function(Env,F1,M1,Gen,Block,BlockID){
   
   Year=Env
   TP=sort(unique(c(F1,M1)))
   Ty=sort(unique(Year))
   Tb=sort(unique(Block))
   
   n=length(F1)

   y0=length(Ty)
   p0=length(TP)
   b0=length(Tb)
   if(BlockID==1&&b0==1)BlockID=0
   UC=UC2Matrix(F1)
   UA=UA2Matrix(F1,M1)
   UD=UD2Matrix(F1,M1,Gen)
   #UAA=UAA2Matrix(F1,M1)
   if(y0==1){
      if(BlockID==1)UB=IndMatrix(Block)
      colnames(UB)=paste("B(",colnames(UB),")",sep="")

   }
   if(y0>1){
     UY=IndMatrix(Env)
     colnames(UY)=paste("E(",colnames(UY),")",sep="")
     UCE=UCE2Matrix(Env,F1)
     UAE=UAE2Matrix(Env,F1,M1)
     UDE=UDE2Matrix(Env,F1,M1,Gen)
     #UAAE=UAAE2Matrix(Env,F1,M1)
     if(BlockID==1){
       B=paste(Env,Block,sep="*")
       UB=IndMatrix(B)
       colnames(UB)=paste("B(",colnames(UB),")",sep="")
     }
   }
   if(y0==1&&BlockID==0){
       U=list(UC,UA,UD)
       names(U)=c("C","A","D")
   }
   else if(y0==1&&BlockID==1){
       U=list(UC,UA,UD,UB)
       names(U)=c("C","A","D","B")
   }
   else if(y0>1&&BlockID==0){
       U=list(UY,UC,UA,UD,UCE,UAE,UDE)
       names(U)=c("E","C","A","D","CE","AE","DE")
   }
   else if(y0>1&&BlockID==1){
       U=list(UY,UC,UA,UD,UCE,UAE,UDE,UB)
       names(U)=c("E","C","A","D","CE","AE","DE","B")
   }
   class(U)="Umatrix"      
   return(U)

}


#### Implementation of design matrices for 2-way ADAA model.

#######################################################################################
#######################################################################################
#######################################################################################

ADAA2Matrix=function(Env,F1,M1,Gen,Block,BlockID){
   Year=Env
   #F=as.vector(F1)
   TP=c(F1,M1)
   F1=as.vector(F1)
   M1=as.vector(M1)

   TP=unique(c(F1,M1))
   Ty=sort(unique(Year))
   Tb=sort(unique(Block))
   
   n=length(F1)

   y0=length(Ty)
   p0=length(TP)
   b0=length(Tb)
   if(BlockID==1&&b0==1)BlockID=0
   UA=UA2Matrix(F1,M1)
   UD=UD2Matrix(F1,M1,Gen)
   UAA=UAA2Matrix(F1,M1)
   if(y0==1){
      if(BlockID==1)UB=IndMatrix(Block)
      colnames(UB)=paste("B(",colnames(UB),")",sep="")

   }
   if(y0>1){
     UY=IndMatrix(Env)
     colnames(UY)=paste("E(",colnames(UY),")",sep="")
     UAE=UAE2Matrix(Env,F1,M1)
     UDE=UDE2Matrix(Env,F1,M1,Gen)
     UAAE=UAAE2Matrix(Env,F1,M1)
     if(BlockID==1){
       B=paste(Env,Block,sep="*")
       UB=IndMatrix(B)
       colnames(UB)=paste("B(",colnames(UB),")",sep="")
     }
   }
   if(y0==1&&BlockID==0){
       U=list(UA,UD,UAA)
       names(U)=c("A","D","AA")
   }
   else if(y0==1&&BlockID==1){
       U=list(UA,UD,UAA,UB)
       names(U)=c("A","D","AA","B")
   }
   else if(y0>1&&BlockID==0){
       U=list(UY,UA,UD,UAA,UAE,UDE,UAAE)
       names(U)=c("E","A","D","AA","AE","DE","AAE")
   }
   else if(y0>1&&BlockID==1){
       U=list(UY,UA,UD,UAA,UAE,UDE,UAAE,UB)
       names(U)=c("E","A","D","AA","AE","DE","AAE","B")
   }
   class(U)="Umatrix"
   return(U)
}



#### Implementation of design matrices for 2-way Single Marker AD model.

#######################################################################################
#######################################################################################
#######################################################################################

SMAD2Matrix=function(Env,MP,F1,M1,Gen,Block,BlockID){
   Year=Env
   TP=sort(unique(c(F1,M1)))
   Ty=sort(unique(Year))
   Tm=sort(unique(MP))   ## number of a marker alleles
   Tb=sort(unique(Block))
   
   n=length(F1)

   y0=length(Ty)
   p0=length(TP)
   b0=length(Tb)
   if(BlockID==1&&b0==1)BlockID=0


  ##########################################
   ## Marker additive effects#######################
   UAm=matrix(0,nrow=n,ncol=length(Tm))
   #dim(UA)
   colnames(UAm)=Tm
   for(i in 1:n){
       j=which(colnames(UAm)==MP[F1[i]])
       UAm[i,j]=UAm[i,j]+1
       
       j=which(colnames(UAm)==MP[M1[i]])
       UAm[i,j]=UAm[i,j]+1
   }
   
   colnames(UAm)=paste("Am(",Tm,")",sep="")
   UAm=DropColumns(UAm)
   
   ##################################################
   ## Marker dominance effects  ############################
   TDM=NULL

   TDm=paste(Tm,Tm,sep="*")

   for(i in 1:(length(Tm)-1)){
       str=paste(TP[i],TP[(i+1):length(Tm)],sep="*")
       TDm=c(TDm,str)
   }

   UDm=matrix(0,nrow=n,ncol=length(TDm))
   
   colnames(UDm)=TDm

   for(i in 1:n){
      if(Gen[i]==0||MP[F1[i]]==MP[M1[i]]){
          j= which(colnames(UDm)==paste(F1[i], M1[i], sep="*"))
          UDm[i,j]=1
      }
      else if(MP[F1[i]]!=MP[M1[i]]){
        if(Gen[i]==1){
            j= which(colnames(UDm)==paste(MP[F1[i]], MP[M1[i]], sep="*"))
            UDm[i,j]=UDm[i,j]+1
        }

        else if(Gen[i]==2){
            j= which(colnames(UDm)==paste(MP[F1[i]], MP[F1[i]], sep="*"))
            UDm[i,j]=UDm[i,j]+0.25
          
            j= which(colnames(UDm)==paste(MP[M1[i]], MP[M1[i]], sep="*"))
            UDm[i,j]=UDm[i,j]+0.25
          
            j= which(colnames(UDm)==paste(MP[F1[i]], MP[M1[i]], sep="*"))
            UDm[i,j]=UDm[i,j]+0.5
          
        }
      
        else if(Gen[i]==3){
            j= which(colnames(UDm)==paste(MP[F1[i]], MP[F1[i]], sep="*"))
            UDm[i,j]=UDm[i,j]+3/8
          
            j= which(colnames(UDm)==paste(MP[M1[i]], MP[M1[i]], sep="*"))
            UDm[i,j]=UDm[i,j]+3/8
          
            j= which(colnames(UDm)==paste(MP[F1[i]], MP[M1[i]], sep="*"))
            UDm[i,j]=UDm[i,j]+1/4
        }
      }
   }
   colnames(UDm)=paste("Dm(",TDm,")",sep="")
   UDm= DropColumns(UDm)


   ##########################################
   ## Additive effects#######################
   UA=matrix(0,nrow=n,ncol=p0)
   dim(UA)
   colnames(UA)=TP
   for(i in 1:n){
       if(Gen[i]==0){
           j=which(colnames(UA)==F1[i])
           UA[i,j]=2
       }
       else if(Gen[i]>0){
           j=which(colnames(UA)==F1[i])
           UA[i,j]=UA[i,j]+1
           j=which(colnames(UA)==M1[i])
           UA[i,j]=UA[i,j]+1
     
       }
   }
   colnames(UA)=paste("A(",TP,")",sep="")
   UA=DropColumns(UA)
   
   ##################################################
   ## Dominance effects  ############################
   TD=NULL

   TD=paste(TP,TP,sep="*")

   for(i in 1:(p0-1)){
       str=paste(TP[i],TP[(i+1):p0],sep="*")
       TD=c(TD,str)
   }

   UD=matrix(0,nrow=n,ncol=length(TD))
   
   colnames(UD)=TD

   for(i in 1:n){
      if(Gen[i]==0){
          j= which(colnames(UD)==paste(F1[i], M1[i], sep="*"))
          UD[i,j]=1
      }
      
      else if(Gen[i]==1){
          j= which(colnames(UD)==paste(F1[i], M1[i], sep="*"))
          UD[i,j]=UD[i,j]+1
     }

      else if(Gen[i]==2){
          j= which(colnames(UD)==paste(F1[i], F1[i], sep="*"))
          UD[i,j]=UD[i,j]+0.25
          
          j= which(colnames(UD)==paste(M1[i], M1[i], sep="*"))
          UD[i,j]=UD[i,j]+0.25
          
          j= which(colnames(UD)==paste(F1[i], M1[i], sep="*"))
          UD[i,j]=UD[i,j]+0.5
          
      }
      
      else if(Gen[i]==3){
          j= which(colnames(UD)==paste(F1[i], F1[i], sep="*"))
          UD[i,j]=UD[i,j]+3/8
          
          j= which(colnames(UD)==paste(M1[i], M1[i], sep="*"))
          UD[i,j]=UD[i,j]+3/8
          
          j= which(colnames(UD)==paste(F1[i], M1[i], sep="*"))
          UD[i,j]=UD[i,j]+0.25
          
      }
   }
   colnames(UD)=paste("D(",TD,")",sep="")
   UD= DropColumns(UD)
   

   ##############################################################
   ### Environmental effects:  random effect#####################
   if(y0>0){
       ### Environmental effects ###################################
       UY=matrix(0,nrow=n,ncol=length(Ty))
       colnames(UY)= Ty
       for(i in 1:n){
	     j=which(colnames(UY)==Year[i])
	     UY[i,j] = 1
       }
       colnames(UY)=paste("E(",Ty,")",sep="")
       #UY
       #dim(UY)
       UY= DropColumns(UY)
   

       
       ##Marker Additive*environment interaction effects
       
       ## comments: AmE interaction
       TAmE=NULL
       for(i in 1:y0){
           str=paste(Ty[i],Tm,sep="*")
           TAmE=c(TAmE,str)
       }
     
       UAmE=matrix(0,nrow=n,ncol=length(TAmE))
       #UAE
       colnames(UAmE)=TAmE
       #head(UAE)
       for(i in 1:n){
	     j= which(colnames(UAmE)==paste(Year[i], MP[F1[i]], sep="*"))
	     UAmE[i, j] = UAmE[i,j]+1

           j= which(colnames(UAmE)==paste(Year[i], MP[M1[i]], sep="*"))
	     UAmE[i, j] = UAmE[i,j]+1
           
       }
       colnames(UAmE)=paste("AmE(",TAmE,")",sep="")
       #dim(UAE)

       UAmE= DropColumns(UAmE)
       #dim(UAE)

       ###########################################################
       ## Marker DE interaction effects #################################
       
       TDmE=NULL

       for(i in 1:y0){
           for(j in 1:length(TDm)){
              TDmE=c(TDmE,paste(Ty[i],TDm[j],sep="*"))
           }
       }

       UDmE=matrix(0,nrow=n,ncol=length(TDmE))
       colnames(UDmE)=TDmE
          
       for(i in 1:n){
           if(Gen[i]==0){
               j= which(colnames(UDmE)==paste(Year[i],MP[F1[i]], MP[F1[i]], sep="*"))
               UDmE[i,j]=1
           }

           else if(Gen[i]==1){
              j= which(colnames(UDmE)==paste(Year[i],MP[F1[i]], MP[M1[i]], sep="*"))
              UDmE[i,j]=1
           }
             
           else if(Gen[i]==2){
              j= which(colnames(UDmE)==paste(Year[i],MP[F1[i]], MP[F1[i]], sep="*"))
              UDmE[i,j]=UDmE[i,j]+0.25
          
              j= which(colnames(UDmE)==paste(Year[i],MP[M1[i]], MP[M1[i]], sep="*"))
              UDmE[i,j]=UDmE[i,j]+0.25
          
              j= which(colnames(UDmE)==paste(Year[i],MP[F1[i]], MP[M1[i]], sep="*"))
              UDmE[i,j]=UDmE[i,j]+0.5
           }

           else if(Gen[i]==3){

              j= which(colnames(UDmE)==paste(Year[i],MP[F1[i]], MP[F1[i]], sep="*"))
              UDmE[i,j]=UDmE[i,j]+3/8
          
              j= which(colnames(UDmE)==paste(Year[i],MP[M1[i]], MP[M1[i]], sep="*"))
              UDmE[i,j]=UDmE[i,j]+3/8

          
              j= which(colnames(UDmE)==paste(Year[i],MP[F1[i]], MP[M1[i]], sep="*"))
              UDmE[i,j]=UDmE[i,j]+1/4
           }  
       }
       colnames(UDmE)=paste("DmE(",TDmE,")",sep="")
       UDmE= DropColumns(UDmE)

       ##Additive*environment interaction effects
       
       ## comments: AE interaction
       TAE=NULL
       for(i in 1:y0){
           str=paste(Ty[i],TP[1:p0],sep="*")
           TAE=c(TAE,str)
       }
     
       UAE=matrix(0,nrow=n,ncol=length(TAE))
       #UAE
       colnames(UAE)=TAE
       #head(UAE)
       for(i in 1:n){
	     j= which(colnames(UAE)==paste(Year[i], F1[i], sep="*"))
	     UAE[i, j] = UAE[i,j]+1

           j= which(colnames(UAE)==paste(Year[i], M1[i], sep="*"))
	     UAE[i, j] = UAE[i,j]+1
           
       }
       colnames(UAE)=paste("AE(",TAE,")",sep="")
       #dim(UAE)

       UAE= DropColumns(UAE)
       #dim(UAE)

       ###########################################################
       ## DE interaction effects #################################
       
       TDE=NULL

       for(i in 1:y0){
           for(j in 1:length(TD)){
              TDE=c(TDE,paste(Ty[i],TD[j],sep="*"))
           }
       }

       UDE=matrix(0,nrow=n,ncol=length(TDE))
       colnames(UDE)=TDE
          
       for(i in 1:n){
           if(Gen[i]==0){
               j= which(colnames(UDE)==paste(Year[i],F1[i], F1[i], sep="*"))
               UDE[i,j]=1

           }

           else if(Gen[i]==1){
              j= which(colnames(UDE)==paste(Year[i],F1[i], M1[i], sep="*"))
              UDE[i,j]=1
           }
             
           else if(Gen[i]==2){
              j= which(colnames(UDE)==paste(Year[i],F1[i], F1[i], sep="*"))
              UDE[i,j]=UDE[i,j]+0.25
          
              j= which(colnames(UDE)==paste(Year[i],M1[i], M1[i], sep="*"))
              UDE[i,j]=UDE[i,j]+0.25
          
              j= which(colnames(UDE)==paste(Year[i],F1[i], M1[i], sep="*"))
              UDE[i,j]=UDE[i,j]+0.5
           }

           else if(Gen[i]==3){

              j= which(colnames(UDE)==paste(Year[i],F1[i], F1[i], sep="*"))
              UDE[i,j]=UDE[i,j]+3/8
          
              j= which(colnames(UDE)==paste(Year[i],M1[i], M1[i], sep="*"))
              UDE[i,j]=UDE[i,j]+3/8

          
              j= which(colnames(UDE)==paste(Year[i],F1[i], M1[i], sep="*"))
              UDE[i,j]=UDE[i,j]+0.25
           }  
       }
       colnames(UDE)=paste("DE(",TDE,")",sep="")
       UDE= DropColumns(UDE)
   }

   if(BlockID==1){
       if(y0==1&&b0>1){
           UB=matrix(0,nrow=n,ncol=b0)
           colnames(UB)=Tb
           for(i in 1:n){
                j=which(colnames(UB)==Block[i])
                UB[i,j]=1
           }
           colnames(UB)=paste("B(",Tb,")",sep="")
           #dim(UB)

           UB= DropColumns(UB)
           #dim(UB)
       }
       else if(y0>1&&b0>1){
           TB=NULL

           for(i in 1:y0){
               str=paste(Ty[i],Tb[1:b0],sep="*")
               TB=c(TB,str)
  
           }

           UB=matrix(0,nrow=n,ncol=length(TB))
           #UB
           colnames(UB)=TB
         
           for(i in 1:n){
	
	          j= which(colnames(UB)==paste(Year[i],Block[i],sep="*"))
	
	          UB[i, j] = 1
           }
           colnames(UB)=paste("B(",TB,")",sep="")
           
           UB= DropColumns(UB)
           
        }
    }
    if(y0==1&&BlockID==0)U=list(UA,UD)
    else if(y0==1&&BlockID==1)U=list(UA,UD,UB)
    else if(y0>1&&BlockID==0)U=list(UY,UA,UD,UAE,UDE)
    else if(y0>1&&BlockID==1)U=list(UY,UA,UD,UAE,UDE,UB)
    
    return(U)

}


#### Implementation of design matrices for GGE model.


GGE=function(Env,Geno,Block, BlockID){
    y0=length(unique(Env))
    if(y0>1){
       UE=IndMatrix(Env)
       colnames(UE)=paste("E(",colnames(UE),")",sep="")
       ge=paste(Env,Geno,sep="*")
       UGE=IndMatrix(ge)
       colnames(UGE)=paste("GE(",colnames(UGE),")",sep="")
       
    }
    UG=IndMatrix(Geno)
    colnames(UG)=paste("G(",colnames(UG),")",sep="")

    if(BlockID==1){
       be=paste(Env,Block,sep="*")
       UB=IndMatrix(be)
       colnames(UB)=paste("B(",colnames(UB),")",sep="")

    }
    if(y0==1){
        if(BlockID==1){
            U=list(UG,UB)
            names(U)=c("G","B")
        }
        else if(BlockID!=1){
            U=list(UG)
            names(U)="G"
        }
    }
    else if(y0>1){
        if(BlockID==1){
           U=list(UE,UG,UGE,UB)
           names(U)=c("E","G","GE","B")
        }
        else if(BlockID!=1){
           U=list(UE,UG,UGE)
           names(U)=c("E","G","GE")
        }
    }
    class(U)="Umatrix"
    return(U)
}

#######################################################################################
#######################################################################################
#######################################################################################

GGEMatrix=function(Env,Geno,Block,BlockID=NULL){
   Year=Env
   Tg=unique(Geno)
   Ty=unique(Year)
   Tb=unique(Block)
   
   n=length(Geno)

   y0=length(Ty)
   g0=length(Tg)
   b0=length(Tb)
   if(is.null(BlockID))BlockID=0
   if(BlockID==1&&b0==1)BlockID=0
   
   ##########################################
   ## Genotypic effects#######################
   UG=matrix(0,nrow=n,ncol=length(Tg))
   
   colnames(UG)=Tg
   for(i in 1:n){
       j=which(colnames(UG)==Geno[i])
       UG[i,j]=1
   }
   colnames(UG)=paste("Geno(",Tg,")",sep="")
   UG=DropColumns(UG)


   ##########################################
   if(y0>1){
      ## Environmental effects#######################
      UE=matrix(0,nrow=n,ncol=length(Ty))
   
      colnames(UE)=Ty
      for(i in 1:n){
         j=which(colnames(UE)==Year[i])
         UE[i,j]=1
      }
      colnames(UE)=paste("Env(",Ty,")",sep="")
      UE=DropColumns(UE)

   
     ##################################################
     ## Genotype*environment effects  ############################
     TGE=NULL

     for(i in 1:y0){
        str=paste(Ty[i],Tg,sep="*")
        TGE=c(TGE,str)
     }

     UGE=matrix(0,nrow=n,ncol=length(TGE))
   
     colnames(UGE)=TGE

     for(i in 1:n){
         j= which(colnames(UGE)==paste(Year[i], Geno[i], sep="*"))
         UGE[i,j]=1
     }
     colnames(UGE)=paste("GE(",TGE,")",sep="")
     UGE=DropColumns(UGE)
   }
   

  
   if(BlockID==1){
       if(y0==1&&b0>1){
           UB=matrix(0,nrow=n,ncol=b0)
           colnames(UB)=Tb
           for(i in 1:n){
                j=which(colnames(UB)==Block[i])
                UB[i,j]=1
           }
           colnames(UB)=paste("B(",Tb,")",sep="")
           #dim(UB)

           UB= DropColumns(UB)
           #dim(UB)
       }
       else if(y0>1&&b0>1){
           TB=NULL

           for(i in 1:y0){
               str=paste(Ty[i],Tb[1:b0],sep="*")
               TB=c(TB,str)
  
           }

           UB=matrix(0,nrow=n,ncol=length(TB))
           UB
           colnames(UB)=TB
         
           for(i in 1:n){
	         j= which(colnames(UB)==paste(Year[i],Block[i],sep="*"))
	         UB[i, j] = 1
           }
           colnames(UB)=paste("B(",TB,")",sep="")
           UB= DropColumns(UB)
        }
    }

    if(BlockID==1){
       be=paste(Env,Block,sep="*")
       UB=IndMatrix(be)
       colnames(UB)=paste("B(",colnames(UB),")",sep="")

    }
    if(y0==1){
        if(BlockID==1){
            U=list(UG,UB)
            names(U)=c("G","B")
        }
        else if(BlockID!=1){
            U=list(UG)
            names(U)="G"
        }
    }
    else if(y0>1){
        if(BlockID==1){
           U=list(UE,UG,UGE,UB)
           names(U)=c("E","G","GE","B")
        }
        else if(BlockID!=1){
           U=list(UE,UG,UGE)
           names(U)=c("E","G","GE")
        }
    }
    class(U)="Umatrix"
    return(U)

}

#P=c(as.vector(unique(F1)),as.vector(unique(M1)))
#unique(P)

UA2Matrix=function(F1,M1){
   TP=unique(c(F1,M1))
   n=length(F1)
   p0=length(TP)
   TA=TP
   UA=matrix(0,nrow=n,ncol=length(TA))
   colnames(UA)=TA
   #head(UAE)
   for(i in 1:n){
	 j= which(colnames(UA)==F1[i])
	 UA[i, j] = UA[i,j]+1

       j= which(colnames(UA)==M1[i])
	 UA[i, j] = UA[i,j]+1
           
   }
   colnames(UA)=paste("A(",TA,")",sep="")
   #dim(UAE)
   UA= DropColumns(UA)
   return(UA)
}

UD2Matrix=function(F1,M1,Gen){
   TP=sort(unique(c(F1,M1)))
   n=length(F1)
   F=M=numeric(n)
   for(i in 1:n){
     if(F1[i]<=M1[i]){
        F[i]=F1[i]
        M[i]=M1[i]
     }
     else{
        F[i]=M1[i]
        M[i]=F1[i]
     }
   }
   F1=F
   M1=M
   p0=length(TP)
   #TD=NULL
   TD=paste(TP,TP,sep="*")
   str=unique(paste(F1,M1,sep="*"))
   #for(i in 1:(p0-1)){
   #    str=paste(TP[i],TP[(i+1):p0],sep="*")
   #    TD=c(TD,str)
   #}
   TD=c(TD,str)
   TD=unique(TD)
   UD=matrix(0,nrow=n,ncol=length(TD))
   colnames(UD)=TD

   for(i in 1:n){
      if(Gen[i]==0){
          j= which(colnames(UD)==paste(F1[i], M1[i], sep="*"))
          UD[i,j]=1
      }
      
      else if(Gen[i]==1){
          if(F1[i]<M1[i])j= which(colnames(UD)==paste(F1[i], M1[i], sep="*"))
          else j= which(colnames(UD)==paste(M1[i], F1[i], sep="*"))
          UD[i,j]=UD[i,j]+1

      }

      else if(Gen[i]==2){
          j= which(colnames(UD)==paste(F1[i], F1[i], sep="*"))
          UD[i,j]=UD[i,j]+0.25
          
          j= which(colnames(UD)==paste(M1[i], M1[i], sep="*"))
          UD[i,j]=UD[i,j]+0.25
          
          if(F1[i]<M1[i])j= which(colnames(UD)==paste(F1[i], M1[i], sep="*"))
          else j= which(colnames(UD)==paste(M1[i], F1[i], sep="*"))
          #j= which(colnames(UD)==paste(F1[i], M1[i], sep="*"))
          #UD[i,j]=UD[i,j]+0.5
          #j= which(colnames(UD)==paste(M1[i], F1[i], sep="*"))
          UD[i,j]=UD[i,j]+0.5

          
      }
      
      else if(Gen[i]==3){
          j= which(colnames(UD)==paste(F1[i], F1[i], sep="*"))
          UD[i,j]=UD[i,j]+3/8
          
          j= which(colnames(UD)==paste(M1[i], M1[i], sep="*"))
          UD[i,j]=UD[i,j]+3/8
          
          if(F1[i]<M1[i])j= which(colnames(UD)==paste(F1[i], M1[i], sep="*"))
          else j= which(colnames(UD)==paste(M1[i], F1[i], sep="*"))
          #j= which(colnames(UD)==paste(F1[i], M1[i], sep="*"))
          #UD[i,j]=UD[i,j]+0.25
          #j= which(colnames(UD)==paste(M1[i], F1[i], sep="*"))
          UD[i,j]=UD[i,j]+0.25

          
      }
   }
   colnames(UD)=paste("D(",TD,")",sep="")
   UD= DropColumns(UD)
   return(UD)
}

##data.frame(F1,M1)
UAA2Matrix=function(F1,M1){
   TP=sort(unique(c(F1,M1)))
   n=length(F1)
   p0=length(TP)
   F=M=numeric(n)
   for(i in 1:n){
     if(F1[i]<=M1[i]){
        F[i]=F1[i]
        M[i]=M1[i]
     }
     else{
        F[i]=M1[i]
        M[i]=F1[i]
     }
   }
   F1=F
   M1=M
   p0=length(TP)
   #TD=NULL
   TD=paste(TP,TP,sep="*")
   str=unique(paste(F1,M1,sep="*"))
   #TD=paste(TP,TP,sep="*")
   #for(i in 1:(p0-1)){
   #    str=paste(TP[i],TP[(i+1):p0],sep="*")
   #    TD=c(TD,str)
   #}
   TD=c(TD,str)
   TD=unique(TD)
   UAA=matrix(0,nrow=n,ncol=length(TD))
   colnames(UAA)=TD

   for(i in 1:n){
      #i=11
      j= which(colnames(UAA)==paste(F1[i], F1[i], sep="*"))
      UAA[i,j]=1
      j= which(colnames(UAA)==paste(M1[i], M1[i], sep="*"))
      UAA[i,j]=UAA[i,j]+1
      if(F1[i]<M1[i])j= which(colnames(UAA)==paste(F1[i], M1[i], sep="*"))
      else j= which(colnames(UAA)==paste(M1[i], F1[i], sep="*"))
      #j= which(colnames(UAA)==paste(F1[i], M1[i], sep="*"))
      UAA[i,j]=UAA[i,j]+2
   }
   #sum(UAA[,11])
   colnames(UAA)=paste("AA(",TD,")",sep="")
   UAA= DropColumns(UAA)
   return(UAA)
}



UAE2Matrix=function(Env,F1,M1){
   Year=Env
   TP=sort(unique(c(F1,M1)))
   Ty=sort(unique(Year))
   n=length(F1)
   y0=length(Ty)
   p0=length(TP)
   TAE=NULL
   for(i in 1:y0){
       str=paste(Ty[i],TP[1:p0],sep="*")
       TAE=c(TAE,str)
   }
   UAE=matrix(0,nrow=n,ncol=length(TAE))
   #UAE
   colnames(UAE)=TAE
   #head(UAE)
   for(i in 1:n){
	 j= which(colnames(UAE)==paste(Year[i], F1[i], sep="*"))
	 UAE[i, j] = UAE[i,j]+1

       j= which(colnames(UAE)==paste(Year[i], M1[i], sep="*"))
	 UAE[i, j] = UAE[i,j]+1
           
   }
   colnames(UAE)=paste("AE(",TAE,")",sep="")
   #dim(UAE)
   UAE= DropColumns(UAE)
   return(UAE)
}


UDE2Matrix=function(Env,F1,M1,Gen){
   Year=Env
   TP=sort(unique(c(F1,M1)))
   Ty=sort(unique(Year))
   n=length(F1)
   y0=length(Ty)
   p0=length(TP)
   F=M=numeric(n)
   for(i in 1:n){
     if(F1[i]<=M1[i]){
        F[i]=F1[i]
        M[i]=M1[i]
     }
     else{
        F[i]=M1[i]
        M[i]=F1[i]
     }
   }
   F1=F
   M1=M
   p0=length(TP)
   TD=paste(TP,TP,sep="*")
   str=unique(paste(F1,M1,sep="*"))
   #TD=paste(TP,TP,sep="*")
   #for(i in 1:(p0-1)){
   #    str=paste(TP[i],TP[(i+1):p0],sep="*")
   #    TD=c(TD,str)
   #}
   TD=c(TD,str)
   TD=unique(TD)
     
   #TD=paste(TP,TP,sep="*")
   #for(i in 1:(p0-1)){
   #    str=paste(TP[i],TP[(i+1):p0],sep="*")
   #    TD=c(TD,str)
   #}
   
   TDE=NULL
   for(i in 1:y0){
      for(j in 1:length(TD)){
         TDE=c(TDE,paste(Ty[i],TD[j],sep="*"))
      }
   }

   UDE=matrix(0,nrow=n,ncol=length(TDE))
   colnames(UDE)=TDE
          
   for(i in 1:n){
        if(Gen[i]==0){
           j= which(colnames(UDE)==paste(Year[i],F1[i], F1[i], sep="*"))
           UDE[i,j]=1

        }

        else if(Gen[i]==1){
           if(F1[i]<M1[i])j= which(colnames(UDE)==paste(Year[i],F1[i], M1[i], sep="*"))
           else j= which(colnames(UDE)==paste(Year[i],M1[i], F1[i], sep="*"))
           #j= which(colnames(UDE)==paste(Year[i],F1[i], M1[i], sep="*"))
           UDE[i,j]=UDE[i,j]+1
        }
             
        else if(Gen[i]==2){
           j= which(colnames(UDE)==paste(Year[i],F1[i], F1[i], sep="*"))
           UDE[i,j]=UDE[i,j]+0.25
          
           j= which(colnames(UDE)==paste(Year[i],M1[i], M1[i], sep="*"))
           UDE[i,j]=UDE[i,j]+0.25
          
           if(F1[i]<M1[i])j= which(colnames(UDE)==paste(Year[i],F1[i], M1[i], sep="*"))
           else j= which(colnames(UDE)==paste(Year[i],M1[i], F1[i], sep="*"))

           UDE[i,j]=UDE[i,j]+0.5
        }

        else if(Gen[i]==3){

           j= which(colnames(UDE)==paste(Year[i],F1[i], F1[i], sep="*"))
           UDE[i,j]=UDE[i,j]+3/8
          
           j= which(colnames(UDE)==paste(Year[i],M1[i], M1[i], sep="*"))
           UDE[i,j]=UDE[i,j]+3/8

          
           if(F1[i]<M1[i])j= which(colnames(UDE)==paste(Year[i],F1[i], M1[i], sep="*"))
           else j= which(colnames(UDE)==paste(Year[i],M1[i], F1[i], sep="*"))

           UDE[i,j]=UDE[i,j]+1/4
        }  
    }
    colnames(UDE)=paste("DE(",TDE,")",sep="")
    UDE= DropColumns(UDE)
    return(UDE)
}


UAAE2Matrix=function(Env,F1,M1){
   Year=Env
   TP=sort(unique(c(F1,M1)))
   Ty=sort(unique(Year))
   n=length(F1)
   y0=length(Ty)
   p0=length(TP)
   F=M=numeric(n)
   for(i in 1:n){
     if(F1[i]<=M1[i]){
        F[i]=F1[i]
        M[i]=M1[i]
     }
     else{
        F[i]=M1[i]
        M[i]=F1[i]
     }
   }
   F1=F
   M1=M
   p0=length(TP)
   #TD=NULL
   TAA=paste(TP,TP,sep="*")
   str=unique(paste(F1,M1,sep="*"))
   #TD=paste(TP,TP,sep="*")
   #for(i in 1:(p0-1)){
   #    str=paste(TP[i],TP[(i+1):p0],sep="*")
   #    TD=c(TD,str)
   #}
   TAA=c(TAA,str)
   TAA=unique(TAA)

      
   #TAA=paste(TP,TP,sep="*")
   #for(i in 1:(p0-1)){
   #    str=paste(TP[i],TP[(i+1):p0],sep="*")
   #    TAA=c(TAA,str)
   #}
   
   TAAE=NULL
   for(i in 1:y0){
      for(j in 1:length(TAA)){
         TAAE=c(TAAE,paste(Ty[i],TAA[j],sep="*"))
      }
   }

   UAAE=matrix(0,nrow=n,ncol=length(TAAE))
   colnames(UAAE)=TAAE
          
   for(i in 1:n){
        j= which(colnames(UAAE)==paste(Year[i],F1[i], F1[i], sep="*"))
        UAAE[i,j]=UAAE[i,j]+1
          
        j= which(colnames(UAAE)==paste(Year[i],M1[i], M1[i], sep="*"))
        UAAE[i,j]=UAAE[i,j]+1
        if(F1[i]<M1[i])j= which(colnames(UAAE)==paste(Year[i],F1[i], M1[i], sep="*"))
        else j= which(colnames(UAAE)==paste(Year[i],M1[i], F1[i], sep="*"))
  
        #j= which(colnames(UAAE)==paste(Year[i],F1[i], M1[i], sep="*"))
        UAAE[i,j]=UAAE[i,j]+2
   }
    colnames(UAAE)=paste("AAE(",TAAE,")",sep="")
    UAAE= DropColumns(UAAE)
    return(UAAE)
}

UC2Matrix=function(F1){
   UC=IndMatrix(F1)
   colnames(UC)=paste("C(",colnames(UC),")",sep="")
   return(UC)
}

UCE2Matrix=function(Env, F1){
   UCE=IndMatrix(paste(Env,F1,sep="*"))
   colnames(UCE)=paste("CE(",colnames(UCE),")",sep="")
   return(UCE)
}

#Env=dat$Loc
#Female=dat$Female
#Male=dat$Male
#Gen=dat$Gen
#Block=dat$Blk
#BlockID=1

AD2Matrix2=function(Env,Female,Male,Gen,Block=NULL,BlockID=NULL){
   
   #Year=Env
   #TP=unique(c(F1,M1))
   F1=Female
   M1=Male
   Ty=sort(unique(Env))
   
   
   #n=length(F1)

   y0=length(Ty)
   #p0=length(TP)
   b0=length(unique(Block))
   if(is.null(Block)){
     BlockID=0
   }else{
      Tb=unique(Block)
   }
   if(is.null(BlockID))BlockID=0
   if(BlockID==1&&b0==1)BlockID=0
   UA=UA2Matrix(F1,M1)
   UD=UD2Matrix(F1,M1,Gen)
   if(y0==1){
      if(BlockID==1){
         UB=IndMatrix(Block)
         colnames(UB)=paste("B(",colnames(UB),")",sep="")
      }
   }
   else if(y0>1){
      UE=IndMatrix(Env)
      colnames(UE)=paste("E(",colnames(UE),")",sep="")
      UAE=UAE2Matrix(Env,F1,M1)
      UDE=UDE2Matrix(Env,F1,M1,Gen)
      if(BlockID==1){
         YB=paste(Block,"<", Env, ">",sep="")
         UB=IndMatrix(YB)
         colnames(UB)=paste("B(",colnames(UB),")",sep="")
      }
   }
   
   if(y0==1&&BlockID==0){
      U=list(UA,UD)
      names(U)=c("A","D")
   }
   else if(y0==1&&BlockID==1){
      U=list(UA,UD,UB)
      names(U)=c("A","D","B")
   }   
   else if(y0>1&&BlockID==0){
      U=list(UE,UA,UD,UAE,UDE)
      names(U)=c("E","A","D","AE","DE")
   }
   else if(y0>1&&BlockID==1){
      U=list(UE,UA,UD,UAE,UDE,UB)
      names(U)=c("E","A","D","AE","DE","B")
   }
   class(U)="Umatrix"
    
   return(U)
}

 
   
AD2MatrixIBD=function(Env,Female,Male,Gen,Block,IBlock){
   #Env=dat$ENV
   #Female=dat$Female
   #Male=dat$Male
   #Gen=dat$Gen
   #Block=dat$Block
   #IBlock=dat$ibloc

   Year=Env
   Ty0=sort(unique(Env))
   Tb0=sort(unique(Block))
   
   y0=length(Ty0)
   b0=length(Tb0)
   
   UA=UA2Matrix(Female,Male)
   UD=UD2Matrix(Female,Male,Gen)

   if(y0==1)YB=Block
   else if(y0>1)YB=paste(Year,Block,sep="*")
   UB=IndMatrix(YB)
   colnames(UB)=paste("B(",colnames(UB),")",sep="")
   #head(UB)

   UAB=UAE2Matrix(YB,Female,Male)
   UDB=UDE2Matrix(YB,Female,Male,Gen)

   BI=paste(YB,IBlock,sep="*")
   UBI=IndMatrix(BI)
   colnames(UBI)=paste("BI(",colnames(UBI),")",sep="")

   if(y0>1){
      UE=IndMatrix(Env)
      colnames(UE)=paste("E(",colnames(UE),")",sep="")
      UAE=UAE2Matrix(Env,Female,Male)
      UDE=UDE2Matrix(Env,Female,Male,Gen)
   }
   head(UBI)
   
   if(y0==1){
      U=list(UA,UD,UB,UBI)
      names(U)=c("A","D","B","BI")
   }
   else if(y0>1){
      U=list(UE,UA,UD,UAE,UDE,UB,UBI)
      names(U)=c("E","A","D","AE","DE","B","BI")
   }
   class(U)="Umatrix"
  
   return(U)
}   
   


GGEMatrixIBD=function(Env,Geno,Block,IBlock){
   Year=Env
   Ty0=sort(unique(Env))
   Tb0=sort(unique(Block))
   
   y0=length(Ty0)
   b0=length(Tb0)
   
   UG=IndMatrix(Geno)
   colnames(UG)=paste("G(",colnames(UG),")",sep="")
   if(y0==1)YB=Block
   else if(y0>1)YB=paste(Year,Block,sep="*")
   UB=IndMatrix(YB)
   colnames(UB)=paste("B(",colnames(UB),")",sep="")
   
   GB=paste(YB,Geno,sep="*")
   UGB=IndMatrix(GB)
   colnames(UGB)=paste("GB(",colnames(UGB),")",sep="")

   BI=paste(YB,IBlock,sep="*")
   UBI=IndMatrix(BI)
   colnames(UBI)=paste("BI(",colnames(UBI),")",sep="")

   if(y0>1){
      UE=IndMatrix(Env)
      colnames(UE)=paste("E(",colnames(UE),")",sep="")
      GE=paste(Env,Geno,sep="*")
      UGE=IndMatrix(GE)
      colnames(UGE)=paste("GE(",colnames(UGE),")",sep="")


      
   }
   
   if(y0==1){
      U=list(UG,UB,UBI)
      names(U)=c("G","B","BI")
   }
   else if(y0>1){
      U=list(UE,UG,UGE,UB,UBI)
      names(U)=c("E","G","UGE","B","BI")
   }
   class(U)="Umatrix"
      
   return(U)
}



      
#M=X
#v=c(1,2)      

CFactor=function(M,v){
   
   n=nrow(M)
   #m1=M[,v]
   #Factor=NULL
   #NAME=NULL
   r=length(v)
   ok=complete.cases(M[,v])
   for(i in 1:r){
     if(i==1){
        NAME=colnames(M)[v[i]]
        Factor=M[,v[i]]
     }
     else if(i>1){
        NAME=paste(NAME,colnames(M)[v[i]],sep=":")
        Factor=paste(Factor,M[,v[i]],sep=":")
     }
   }
   #F=matrix(0,nrow=n,ncol=1)
   #F[,1]=Factor
   #colnames(F)=NAME
   n=length(ok)
   id=which(ok==F)
   Factor[id]=NA
   return(Factor)
}


Edata=function(A,Geno){
   n13=length(A)
   a=numeric(length(Geno))
   for(i in 1:n13){
      j=which(Geno==i)
      a[j]=A[i]
   }
   return(a)
}

PedSplit=function(S){
    r=length(S)
    P=matrix(0,r,4)
    Gen=numeric(r)
    for (i in 1:r){
       #i=1
       s=S[i]
       s1=strsplit(s,"//")[[1]]
       if(length(s1)==2){
          cross=4
          a1=s1[1]
          a2=s1[2]
       }
       else if(length(s1)==1){
          a1=s1
          a2=s1
       }
       p=strsplit(a1,"/")[[1]]
       if(length(p)==1){
          P[i,1]=a1
          P[i,2]=a1
       }
       else if(length(p)==2){
          P[i,1]=p[1]
          P[i,2]=p[2]
       }
       if(a1==a2){
          P[i,3]=P[i,1]
          P[i,4]=P[i,2]
       }
       else if(a1!=a2){
          p=strsplit(a2,"/")[[1]]
          if(length(p)==1){
             P[i,3]=a2
             P[i,4]=a2
          }
          else if(length(p)==2){
             P[i,3]=p[1]
             P[i,4]=p[2]
          }
       }
       t=unique(P[i,])
       if(length(t)==1)Gen[i]=0
       else Gen[i]=2
    }
    P=data.frame(P,Gen)
    return(P)
}

DropColumns = function(A){
	
	for(i in colnames(A)){
	   x=sum(ifelse(A[,which(colnames(A)==i)]==0, 0, 1))
	
	   if(x==0){
		A= A[, -which(colnames(A)==i)]
	   }

      }
	return(A)
	
}


CNumeric=function(M,v){
   
   n=nrow(M)
   r=length(v)
   y=numeric(n)
   for(i in 1:r){
     if(i==1){
        NAME=colnames(M)[v[i]]
        a=M[,v[i]]
     }
     else if(i>1){
        NAME=paste(NAME,colnames(M)[v[i]],sep=":")
        a=a*M[,v[i]]
     }
   }
   y=a
   return(y)
}
#X=gdata$X
BigX=function(X,...){
   
   X0=NULL
   xk=length(X)
   xnames=NULL
   #i=2
   for(i in 1:xk)X0=cbind(X0,as.matrix(X[[i]]))
   colnames(X0)   
   return(X0)
}

##Y=Y1
##Pedigree=Ped
#Y=y
qg.model=function(Y,Ped,Model,Cross,X=NULL,CC=NULL,...){
   if(is.null(X))return(genmod.data(Y,Ped,Model,Cross))
   else return(genmod.data(Y,Ped,Model,Cross,X))
}

genmod.data=function(Y,Ped,Model,Cross,X=NULL,CC=NULL,...){
    
   nc=ncol(Ped)
   #for(i in 1:nc)Ped[,i]=as.vector(Ped[,i])

   if(class(Y)=="numeric"){
      Y=data.frame(Y)
      colnames(Y)=paste(deparse(substitute(Y)))
   }
   n=nrow(Y)
   if(is.null(X)){
     #X=rep(1,n)
     X=matrix(1,n,1)
     colnames(X)="mu"
   }
   if(Cross=="Two"){
      nc=ncol(Ped)
      Env=Ped[,1]
      F1=Ped[,2]
      M1=Ped[,3]
      Gen=Ped[,4]
      if(nc==5){
         BlockID=1
         Block=Ped[,5]
      }
      else if(nc==4){
          BlockID=0
          Block=NULL
      }
      if(Model=="ADAA"){
          U=ADAA2Matrix(Env,F1,M1,Gen,Block,BlockID)
          mk=length(U)
          if(mk==3)VC=c(2,1,4,1)
          else if(mk==4)VC=c(2,1,4,1,1)
          else if(mk==7)VC=c(1,2,1,4,2,1,4,1)
          else if(mk==8)VC=c(1,2,1,4,2,1,4,1,1)
             
      }
      else if(Model=="AD"){
          U=AD2Matrix(Env,F1,M1,Gen,Block,BlockID)
          mk=length(U)
          if(mk==2)VC=c(2,1,1)
          else if(mk==3)VC=c(2,1,1,1)
          else if(mk==5)VC=c(1,2,1,2,1,1)
          else if(mk==6)VC=c(1,2,1,2,1,1,1)

      }
      else if(Model=="ADC"){
          U=ADC2Matrix(Env,F1,M1,Gen,Block,BlockID)
          mk=length(U)
          if(mk==3)VC=c(1,2,1,1)
          else if(mk==4)VC=c(1,2,1,1,1)
          else if(mk==7)VC=c(1,1,2,1,1,2,1,1)
          else if(mk==8)VC=c(1,1,2,1,1,2,1,1,1)

      }
      else if(Model=="ADM"){
          U=ADM2Matrix(Env,F1,M1,Gen,Block,BlockID)
          mk=length(U)
          if(mk==4)VC=c(2,1,2,1,1)
          else if(mk==5)VC=c(2,1,2,1,1,1)
          else if(mk==9)VC=c(1,2,1,2,1,2,1,2,1,1)
          else if(mk==10)VC=c(1,2,1,2,1,2,1,2,1,1,1)


      }
      #result=list(Y=Y,X=X,U=U,VC=VC,Model=Model)
      #result$call=match.call()
      #class(result)="genmod.data"
      #return(result)
   }
   if(Cross=="Four"){
      nc=ncol(Ped)
      Env=Ped[,1]
      F1=Ped[,2]
      M1=Ped[,3]
      F2=Ped[,4]
      M2=Ped[,5]
      Gen=Ped[,6]
      if(nc==7){
         BlockID=1
         Block=Ped[,7]
      }
      else if(nc==6){
         BlockID=0
         Block=NULL
      }
      if(Model=="ADAA"){
          #U=ADAA2Matrix(Env,F1,M1,Gen,Block,BlockID)
          #mk=length(U)
          #if(mk==3)VC=c(2,1,4,1)
          #else if(mk==4)VC=c(2,1,4,1,1)
          #else if(mk==7)VC=c(1,2,1,4,2,1,4,1)
          #else if(mk==8)VC=c(1,2,1,4,2,1,4,1,1)
          cat("No such a model available\n")
          return(0)
             
      }
      else if(Model=="AD"){
          if(is.null(CC))U=AD4Matrix(Env,F1,M1,F2,M2,Gen,Block,BlockID)
          else U=AD4Matrix(Env,F1,M1,F2,M2,Gen,Block,BlockID,CC)
          mk=length(U)
          if(mk==2)VC=c(2,1,1)
          else if(mk==3)VC=c(2,1,1,1)
          else if(mk==5)VC=c(1,2,1,2,1,1)
          else if(mk==6)VC=c(1,2,1,2,1,1,1)
      }
      else if(Model=="ADC"){ ##done
          U=ADC4Matrix(Env,F1,M1,F2,M2,Gen,Block,BlockID)
          mk=length(U)
          if(mk==3)VC=c(1,2,1,1)
          else if(mk==4)VC=c(1,2,1,1,1)
          else if(mk==7)VC=c(1,1,2,1,1,2,1,1)
          else if(mk==8)VC=c(1,1,2,1,1,2,1,1,1)

      }
      else if(Model=="ADM"){
          cat("No such a model available\n")
          return(0)
          #U=ADM2Matrix(Env,F1,M1,Gen,Block,BlockID)
          #mk=length(U)
          #if(mk==4)VC=c(2,1,2,1,1)
          #else if(mk==5)VC=c(2,1,2,1,1,1)
          #else if(mk==9)VC=c(1,2,1,2,1,2,1,2,1,1)
          #else if(mk==10)VC=c(1,2,1,2,1,2,1,2,1,1,1)
      }
   }
   X0=list()
   X0[[1]]=X
   result=list(Y=Y,X=X0,U=U,VC=VC,C=NULL)
   #result$call=match.call()
   #class(result)="genmod.data"
   return(result)
}

combid <- function(n, r, i) {

  # http://msdn.microsoft.com/en-us/library/aa289166(VS.71).aspx
  # http://en.wikipedia.org/wiki/Combinadic

  if(i < 1 | i > choose(n,r)) stop("'i' must be 0 < i <= n!/(n-r)!")

  largestV <- function(n, r, i) {
    #v <- n-1
    v <- n                                  # Adjusted for one-based indexing
    #while(choose(v,r) > i) v <- v-1
    while(choose(v,r) >= i) v <- v-1        # Adjusted for one-based indexing
    return(v)
  }

  res <- rep(NA,r)
  for(j in 1:r) {
    res[j] <- largestV(n,r,i)
    i <- i-choose(res[j],r)
    n <- res[j]
    r <- r-1
  }
  res <- res + 1
  return(res)
}

Comb=function(n,r,id1,id2){
   nrow=id2-id1+1
   A=matrix(0,nrow,r)
   for(i in id1:id2){
     a=sort(combid(n,r,i))
     A[(i-id1+1),]=a
   }
   return(A)
}



LargeInv=function(A){
   sys=1
   m=ncol(A)
   if(m<=20)return(ginv(A))
   else if(m<=50)return(subinv20(A))
   else if(m<=100)return(subinv50(A))
   else if(m<=200)return(subinv100(A))
   else if(m<=500)return(subinv200(A))
   else if(m<=1000)return(subinv500(A))
   else if(m<=2000)return(subinv1000(A))
   else if(m<=4000)return(subinv2000(A))
   else if(m<=8000)return(subinv4000(A))
   else if(m<=16000)return(subinv8000(A))
   else if(m<=32000)return(subinv16000(A))
   else if(m<=64000)return(subinv32000(A))
   else if(m<=128000)return(subinv64000(A))

   for(i in 1:(m-1)){
       for(j in (i+1):m){
          if(A[i,j]!=A[j,i]){
             sys=0
             break
          }
       }
       if(sys==0)break
   }
   p=as.integer(m/2)
   s=m-p
   P=A[1:p,1:p]
   R=A[(p+1):m,1:p]
   Q=A[1:p,(p+1):m]
   S=A[(p+1):m,(p+1):m]
    
   if(p<=20){
      S10=ginv(S)
      P1=ginv(P-Q%*%S10%*%R)
   }
   else if(p>20&&p<=50){
     S10=subinv20(S)
     P1=subinv20(P-Q%*%S10%*%R)
   }
   else if(p>50&&p<=100){
     S10=subinv50(S)
     P1=subinv50(P-Q%*%S10%*%R)
   }
   else if(p>100&&p<=200){
     S10=subinv100(S)
     P1=subinv100(P-Q%*%S10%*%R)
   }
   else if(p>200&&p<=500){
     S10=subinv200(S)
     P1=subinv200(P-Q%*%S10%*%R)
   }

  else if(p>500&&p<=1000){
     S10=subinv500(S)
     P1=subinv500(P-Q%*%S10%*%R)
  }
  else if(p>1000&&p<=2000){
     S10=subinv1000(S)
     P1=subinv1000(P-Q%*%S10%*%R)
  }
  else if(p>2000&&p<=4000){
     S10=subinv2000(S)
     P1=subinv2000(P-Q%*%S10%*%R)
   }
  else if(p>4000&&p<=8000){
     S10=subinv4000(S)
     P1=subinv4000(P-Q%*%S10%*%R)
  }
  else if (p>8000&&p<=16000){
      S10=subinv8000(S)
      P1=subinv8000(P-Q%*%S10%*%R)
   }
  else if (p>16000&&p<=32000){
      S10=subinv16000(S)
      P1=subinv16000(P-Q%*%S10%*%R)
   }
   else if (p>32000&&p<=64000){
      S10=subinv32000(S)
      P1=subinv32000(P-Q%*%S10%*%R)
   }
   else if (p>64000&&p<=128000){
      S10=subinv64000(S)
      P1=subinv64000(P-Q%*%S10%*%R)
   }
   else{
      S10=subinv128000(S)
      P1=subinv128000(P-Q%*%S10%*%R)
   }
   Q1=-P1%*%Q%*%S10
   if(sys==0)R1=-(S10%*%R)%*%P1
   S1=S10+(S10%*%R)%*%P1%*%Q%*%S10
   A1=cbind(P1,Q1)
   if(sys==0)A1=rbind(cbind(P1,Q1),cbind(R1,S1))
   else A1=rbind(cbind(P1,Q1),cbind(t(Q1),S1))
   return(A1)

}

subinv128000=function(A){
   sys=1
   m=ncol(A)
   if(m<=20)return(ginv(A))
   else if(m<=50)return(subinv20(A))
   else if(m<=100)return(subinv50(A))
   else if(m<=200)return(subinv100(A))
   else if(m<=500)return(subinv200(A))
   else if(m<=1000)return(subinv500(A))
   else if(m<=2000)return(subinv1000(A))
   else if(m<=4000)return(subinv2000(A))
   else if(m<=8000)return(subinv4000(A))
   else if(m<=16000)return(subinv8000(A))
   else if(m<=32000)return(subinv16000(A))
   else if(m<=64000)return(subinv32000(A))
   else if(m<=128000)return(subinv64000(A))

   for(i in 1:(m-1)){
       for(j in (i+1):m){
          if(A[i,j]!=A[j,i]){
             sys=0
             break
          }
       }
       if(sys==0)break
   }

   p=as.integer(m/2)
   s=m-p
   P=A[1:p,1:p]
   R=A[(p+1):m,1:p]
   Q=A[1:p,(p+1):m]
   S=A[(p+1):m,(p+1):m]
    
   if(p<=20){
      S10=ginv(S)
      P1=ginv(P-Q%*%S10%*%R)
   }
   else if(p>20&&p<=50){
      S10=subinv20(S)
      P1=subinv20(P-Q%*%S10%*%R)
   }
     
   else if(p>50&&p<=100){
     S10=subinv50(S)
     P1=subinv50(P-Q%*%S10%*%R)
   }
   else if(p>100&&p<=200){
     S10=subinv100(S)
     P1=subinv100(P-Q%*%S10%*%R)
   }
   else if(p>200&&p<=500){
     S10=subinv200(S)
     P1=subinv200(P-Q%*%S10%*%R)
   }
   else if(p>500&&p<=1000){
     S10=subinv500(S)
     P1=subinv500(P-Q%*%S10%*%R)
   }
   else if(p>1000&&p<=2000){
     S10=subinv1000(S)
     P1=subinv1000(P-Q%*%S10%*%R)
   }
   else if(p>2000&&p<=4000){
     S10=subinv2000(S)
     P1=subinv2000(P-Q%*%S10%*%R)
   }
   else if(p>4000&&p<=8000){
     S10=subinv4000(S)
     P1=subinv4000(P-Q%*%S10%*%R)
   }
   else if(p>8000&&p<=16000){
     S10=subinv8000(S)
     P1=subinv8000(P-Q%*%S10%*%R)
   }
   else if(p>16000&&p<=32000){
     S10=subinv16000(S)
     P1=subinv16000(P-Q%*%S10%*%R)
   }
   else if(p>32000&&p<=64000){
     S10=subinv32000(S)
     P1=subinv32000(P-Q%*%S10%*%R)
   }
   else if(p>64000){
     S10=subinv64000(S)
     P1=subinv64000(P-Q%*%S10%*%R)
   }

   Q1=-P1%*%Q%*%S10
   if(sys==0)R1=-(S10%*%R)%*%P1
   S1=S10+(S10%*%R)%*%P1%*%Q%*%S10
   A1=cbind(P1,Q1)
   if(sys==0)A1=rbind(cbind(P1,Q1),cbind(R1,S1))
   else A1=rbind(cbind(P1,Q1),cbind(t(Q1),S1))
   return(A1)
}



subinv64000=function(A){
   sys=1
   m=ncol(A)
   if(m<=20)return(ginv(A))
   else if(m<=50)return(subinv20(A))
   else if(m<=100)return(subinv50(A))
   else if(m<=200)return(subinv100(A))
   else if(m<=500)return(subinv200(A))
   else if(m<=1000)return(subinv500(A))
   else if(m<=2000)return(subinv1000(A))
   else if(m<=4000)return(subinv2000(A))
   else if(m<=8000)return(subinv4000(A))
   else if(m<=16000)return(subinv8000(A))
   else if(m<=32000)return(subinv16000(A))
   else if(m<=64000)return(subinv32000(A))


   
   for(i in 1:(m-1)){
       for(j in (i+1):m){
          if(A[i,j]!=A[j,i]){
             sys=0
             break
          }
       }
       if(sys==0)break
   }

   p=as.integer(m/2)
   s=m-p
   P=A[1:p,1:p]
   R=A[(p+1):m,1:p]
   Q=A[1:p,(p+1):m]
   S=A[(p+1):m,(p+1):m]
    
   if(p<=20){
      S10=ginv(S)
      P1=ginv(P-Q%*%S10%*%R)
   }
   else if(p>20&&p<=50){
      S10=subinv20(S)
      P1=subinv20(P-Q%*%S10%*%R)
   }
     
   else if(p>50&&p<=100){
     S10=subinv50(S)
     P1=subinv50(P-Q%*%S10%*%R)
   }
   else if(p>100&&p<=200){
     S10=subinv100(S)
     P1=subinv100(P-Q%*%S10%*%R)
   }
   else if(p>200&&p<=500){
     S10=subinv200(S)
     P1=subinv200(P-Q%*%S10%*%R)
   }
   else if(p>500&&p<=1000){
     S10=subinv500(S)
     P1=subinv500(P-Q%*%S10%*%R)
   }
   else if(p>1000&&p<=2000){
     S10=subinv1000(S)
     P1=subinv1000(P-Q%*%S10%*%R)
   }
   else if(p>2000&&p<=4000){
     S10=subinv2000(S)
     P1=subinv2000(P-Q%*%S10%*%R)
   }
   else if(p>4000&&p<=8000){
     S10=subinv4000(S)
     P1=subinv4000(P-Q%*%S10%*%R)
   }
   else if(p>8000&&p<=16000){
     S10=subinv8000(S)
     P1=subinv8000(P-Q%*%S10%*%R)
   }
   else if(p>16000&&p<=32000){
     S10=subinv16000(S)
     P1=subinv16000(P-Q%*%S10%*%R)
   }
   else if(p>32000){
     S10=subinv32000(S)
     P1=subinv32000(P-Q%*%S10%*%R)
   }


   Q1=-P1%*%Q%*%S10
   if(sys==0)R1=-(S10%*%R)%*%P1
   S1=S10+(S10%*%R)%*%P1%*%Q%*%S10
   A1=cbind(P1,Q1)
   if(sys==0)A1=rbind(cbind(P1,Q1),cbind(R1,S1))
   else A1=rbind(cbind(P1,Q1),cbind(t(Q1),S1))
   return(A1)
}




subinv32000=function(A){
   sys=1
   m=ncol(A)
   if(m<=20)return(ginv(A))
   else if(m<=50)return(subinv20(A))
   else if(m<=100)return(subinv50(A))
   else if(m<=200)return(subinv100(A))
   else if(m<=500)return(subinv200(A))
   else if(m<=1000)return(subinv500(A))
   else if(m<=2000)return(subinv1000(A))
   else if(m<=4000)return(subinv2000(A))
   else if(m<=8000)return(subinv4000(A))
   else if(m<=16000)return(subinv8000(A))
   else if(m<=32000)return(subinv16000(A))


   
   for(i in 1:(m-1)){
       for(j in (i+1):m){
          if(A[i,j]!=A[j,i]){
             sys=0
             break
          }
       }
       if(sys==0)break
   }

   p=as.integer(m/2)
   s=m-p
   P=A[1:p,1:p]
   R=A[(p+1):m,1:p]
   Q=A[1:p,(p+1):m]
   S=A[(p+1):m,(p+1):m]
    
   if(p<=20){
      S10=ginv(S)
      P1=ginv(P-Q%*%S10%*%R)
   }
   else if(p>20&&p<=50){
      S10=subinv20(S)
      P1=subinv20(P-Q%*%S10%*%R)
   }
     
   else if(p>50&&p<=100){
     S10=subinv50(S)
     P1=subinv50(P-Q%*%S10%*%R)
   }
   else if(p>100&&p<=200){
     S10=subinv100(S)
     P1=subinv100(P-Q%*%S10%*%R)
   }
   else if(p>200&&p<=500){
     S10=subinv200(S)
     P1=subinv200(P-Q%*%S10%*%R)
   }
   else if(p>500&&p<=1000){
     S10=subinv500(S)
     P1=subinv500(P-Q%*%S10%*%R)
   }
   else if(p>1000&&p<=2000){
     S10=subinv1000(S)
     P1=subinv1000(P-Q%*%S10%*%R)
   }
   else if(p>2000&&p<=4000){
     S10=subinv2000(S)
     P1=subinv2000(P-Q%*%S10%*%R)
   }
   else if(p>4000&&p<=8000){
     S10=subinv4000(S)
     P1=subinv4000(P-Q%*%S10%*%R)
   }
   else if(p>8000&&p<=16000){
     S10=subinv8000(S)
     P1=subinv8000(P-Q%*%S10%*%R)
   }
   else if(p>16000){
     S10=subinv16000(S)
     P1=subinv16000(P-Q%*%S10%*%R)
   }

   Q1=-P1%*%Q%*%S10
   if(sys==0)R1=-(S10%*%R)%*%P1
   S1=S10+(S10%*%R)%*%P1%*%Q%*%S10
   A1=cbind(P1,Q1)
   if(sys==0)A1=rbind(cbind(P1,Q1),cbind(R1,S1))
   else A1=rbind(cbind(P1,Q1),cbind(t(Q1),S1))
   return(A1)
}




subinv16000=function(A){
   sys=1
   m=ncol(A)
   if(m<=20)return(ginv(A))
   else if(m<=50)return(subinv20(A))
   else if(m<=100)return(subinv50(A))
   else if(m<=200)return(subinv100(A))
   else if(m<=500)return(subinv200(A))
   else if(m<=1000)return(subinv500(A))
   else if(m<=2000)return(subinv1000(A))
   else if(m<=4000)return(subinv2000(A))
   else if(m<=8000)return(subinv4000(A))
   else if(m<=16000)return(subinv8000(A))

   for(i in 1:(m-1)){
       for(j in (i+1):m){
          if(A[i,j]!=A[j,i]){
             sys=0
             break
          }
       }
       if(sys==0)break
   }

   p=as.integer(m/2)
   s=m-p
   P=A[1:p,1:p]
   R=A[(p+1):m,1:p]
   Q=A[1:p,(p+1):m]
   S=A[(p+1):m,(p+1):m]
    
   if(p<=20){
      S10=ginv(S)
      P1=ginv(P-Q%*%S10%*%R)
   }
   else if(p>20&&p<=50){
      S10=subinv20(S)
      P1=subinv20(P-Q%*%S10%*%R)
   }
     
   else if(p>50&&p<=100){
     S10=subinv50(S)
     P1=subinv50(P-Q%*%S10%*%R)
   }
   else if(p>100&&p<=200){
     S10=subinv100(S)
     P1=subinv100(P-Q%*%S10%*%R)
   }
   else if(p>200&&p<=500){
     S10=subinv200(S)
     P1=subinv200(P-Q%*%S10%*%R)
   }
   else if(p>500&&p<=1000){
     S10=subinv500(S)
     P1=subinv500(P-Q%*%S10%*%R)
   }
   else if(p>1000&&p<=2000){
     S10=subinv1000(S)
     P1=subinv1000(P-Q%*%S10%*%R)
   }
   else if(p>2000&&p<=4000){
     S10=subinv2000(S)
     P1=subinv2000(P-Q%*%S10%*%R)
   }
   else if(p>4000&&p<=8000){
     S10=subinv4000(S)
     P1=subinv4000(P-Q%*%S10%*%R)
   }
   else if(p>8000){
     S10=subinv8000(S)
     P1=subinv8000(P-Q%*%S10%*%R)
   }

   Q1=-P1%*%Q%*%S10
   if(sys==0)R1=-(S10%*%R)%*%P1
   S1=S10+(S10%*%R)%*%P1%*%Q%*%S10
   A1=cbind(P1,Q1)
   if(sys==0)A1=rbind(cbind(P1,Q1),cbind(R1,S1))
   else A1=rbind(cbind(P1,Q1),cbind(t(Q1),S1))
   return(A1)
}

subinv8000=function(A){
   sys=1
   m=ncol(A)
   if(m<=20)return(ginv(A))
   else if(m<=50)return(subinv20(A))
   else if(m<=100)return(subinv50(A))
   else if(m<=200)return(subinv100(A))
   else if(m<=500)return(subinv200(A))
   else if(m<=1000)return(subinv500(A))
   else if(m<=2000)return(subinv1000(A))
   else if(m<=4000)return(subinv2000(A))
   else if(m<=8000)return(subinv8000(A))
   
   for(i in 1:(m-1)){
       for(j in (i+1):m){
          if(A[i,j]!=A[j,i]){
             sys=0
             break
          }
       }
       if(sys==0)break
   }

   p=as.integer(m/2)
   s=m-p
   P=A[1:p,1:p]
   R=A[(p+1):m,1:p]
   Q=A[1:p,(p+1):m]
   S=A[(p+1):m,(p+1):m]
    
   if(p<=20){
      S10=ginv(S)
      P1=ginv(P-Q%*%S10%*%R)
   }
   else if(p>20&&p<=50){
      S10=subinv20(S)
      P1=subinv20(P-Q%*%S10%*%R)
   }
     
   else if(p>50&&p<=100){
     S10=subinv50(S)
     P1=subinv50(P-Q%*%S10%*%R)
   }
   else if(p>100&&p<=200){
     S10=subinv100(S)
     P1=subinv100(P-Q%*%S10%*%R)
   }
   else if(p>200&&p<=500){
     S10=subinv200(S)
     P1=subinv200(P-Q%*%S10%*%R)
   }
   else if(p>500&&p<=1000){
     S10=subinv500(S)
     P1=subinv500(P-Q%*%S10%*%R)
   }
   else if(p>1000&&p<=2000){
     S10=subinv1000(S)
     P1=subinv1000(P-Q%*%S10%*%R)
   }
   else if(p>2000&&p<=4000){
     S10=subinv2000(S)
     P1=subinv2000(P-Q%*%S10%*%R)
   }
   else if(p>4000){
     S10=subinv4000(S)
     P1=subinv4000(P-Q%*%S10%*%R)
   }

   Q1=-P1%*%Q%*%S10
   if(sys==0)R1=-(S10%*%R)%*%P1
   S1=S10+(S10%*%R)%*%P1%*%Q%*%S10
   A1=cbind(P1,Q1)
   if(sys==0)A1=rbind(cbind(P1,Q1),cbind(R1,S1))
   else A1=rbind(cbind(P1,Q1),cbind(t(Q1),S1))
   return(A1)
}

subinv4000=function(A){
   sys=1
   m=ncol(A)
   if(m<=20)return(ginv(A))
   else if(m<=50)return(subinv20(A))
   else if(m<=100)return(subinv50(A))
   else if(m<=200)return(subinv100(A))
   else if(m<=500)return(subinv200(A))
   else if(m<=1000)return(subinv500(A))
   else if(m<=2000)return(subinv1000(A))
   else if(m<=4000)return(subinv2000(A))


   for(i in 1:(m-1)){
       for(j in (i+1):m){
          if(A[i,j]!=A[j,i]){
             sys=0
             break
          }
       }
       if(sys==0)break
   }

   p=as.integer(m/2)
   s=m-p
   P=A[1:p,1:p]
   R=A[(p+1):m,1:p]
   Q=A[1:p,(p+1):m]
   S=A[(p+1):m,(p+1):m]
    
   if(p<=20){
      S10=ginv(S)
      P1=ginv(P-Q%*%S10%*%R)
   }
   else if(p>20&&p<=50){
      S10=subinv20(S)
      P1=subinv20(P-Q%*%S10%*%R)
   }
     
   else if(p>50&&p<=100){
     S10=subinv50(S)
     P1=subinv50(P-Q%*%S10%*%R)
   }
   else if(p>100&&p<=200){
     S10=subinv100(S)
     P1=subinv100(P-Q%*%S10%*%R)
   }
   else if(p>200&&p<=500){
     S10=subinv200(S)
     P1=subinv200(P-Q%*%S10%*%R)
   }
   else if(p>500&&p<=1000){
     S10=subinv500(S)
     P1=subinv500(P-Q%*%S10%*%R)
   }
   else if(p>1000&&p<=2000){
     S10=subinv1000(S)
     P1=subinv1000(P-Q%*%S10%*%R)
   }
   else if(p>2000){
     S10=subinv2000(S)
     P1=subinv2000(P-Q%*%S10%*%R)
   }


   Q1=-P1%*%Q%*%S10
   if(sys==0)R1=-(S10%*%R)%*%P1
   S1=S10+(S10%*%R)%*%P1%*%Q%*%S10
   A1=cbind(P1,Q1)
   if(sys==0)A1=rbind(cbind(P1,Q1),cbind(R1,S1))
   else A1=rbind(cbind(P1,Q1),cbind(t(Q1),S1))
   return(A1)
}




subinv2000=function(A){
   sys=1
   m=ncol(A)
   if(m<=20)return(ginv(A))
   else if(m<=50)return(subinv20(A))
   else if(m<=100)return(subinv50(A))
   else if(m<=200)return(subinv100(A))
   else if(m<=500)return(subinv200(A))
   else if(m<=1000)return(subinv500(A))
   else if(m<=2000)return(subinv1000(A))

   for(i in 1:(m-1)){
       for(j in (i+1):m){
          if(A[i,j]!=A[j,i]){
             sys=0
             break
          }
       }
       if(sys==0)break
   }

   p=as.integer(m/2)
   s=m-p
   P=A[1:p,1:p]
   R=A[(p+1):m,1:p]
   Q=A[1:p,(p+1):m]
   S=A[(p+1):m,(p+1):m]
    
   if(p<=20){
      S10=ginv(S)
      P1=ginv(P-Q%*%S10%*%R)
   }
   else if(p>20&&p<=50){
      S10=subinv20(S)
      P1=subinv20(P-Q%*%S10%*%R)
   }
     
   else if(p>50&&p<=100){
     S10=subinv50(S)
     P1=subinv50(P-Q%*%S10%*%R)
   }
   else if(p>100&&p<=200){
     S10=subinv100(S)
     P1=subinv100(P-Q%*%S10%*%R)
   }
   else if(p>200&&p<=500){
     S10=subinv200(S)
     P1=subinv200(P-Q%*%S10%*%R)
   }
   else if(p>500&&p<=1000){
     S10=subinv500(S)
     P1=subinv500(P-Q%*%S10%*%R)
   }
   else if(p>1000){
     S10=subinv1000(S)
     P1=subinv1000(P-Q%*%S10%*%R)
   }

   Q1=-P1%*%Q%*%S10
   if(sys==0)R1=-(S10%*%R)%*%P1
   S1=S10+(S10%*%R)%*%P1%*%Q%*%S10
   A1=cbind(P1,Q1)
   if(sys==0)A1=rbind(cbind(P1,Q1),cbind(R1,S1))
   else A1=rbind(cbind(P1,Q1),cbind(t(Q1),S1))
   return(A1)
}

subinv1000=function(A){
   sys=1
   m=ncol(A)
   if(m<=20)return(ginv(A))
   else if(m<=50)return(subinv20(A))
   else if(m<=100)return(subinv50(A))
   else if(m<=200)return(subinv100(A))
   else if(m<=500)return(subinv200(A))
   else if(m<=1000)return(subinv500(A))
   for(i in 1:(m-1)){
       for(j in (i+1):m){
          if(A[i,j]!=A[j,i]){
             sys=0
             break
          }
       }
       if(sys==0)break
   }

   p=as.integer(m/2)
   s=m-p
   P=A[1:p,1:p]
   R=A[(p+1):m,1:p]
   Q=A[1:p,(p+1):m]
   S=A[(p+1):m,(p+1):m]
    
   if(p<=20){
      S10=ginv(S)
      P1=ginv(P-Q%*%S10%*%R)
   }
   else if(p>20&&p<=50){
      S10=subinv20(S)
      P1=subinv20(P-Q%*%S10%*%R)
   }
     
   else if(p>50&&p<=100){
     S10=subinv50(S)
     P1=subinv50(P-Q%*%S10%*%R)
   }
   else if(p>100&&p<=200){
     S10=subinv100(S)
     P1=subinv100(P-Q%*%S10%*%R)
   }
   else if(p>200&&p<=500){
     S10=subinv200(S)
     P1=subinv200(P-Q%*%S10%*%R)
   }
   else if(p>500){
     S10=subinv500(S)
     P1=subinv500(P-Q%*%S10%*%R)
   }
 
   Q1=-P1%*%Q%*%S10
   if(sys==0)R1=-(S10%*%R)%*%P1
   S1=S10+(S10%*%R)%*%P1%*%Q%*%S10
   A1=cbind(P1,Q1)
   if(sys==0)A1=rbind(cbind(P1,Q1),cbind(R1,S1))
   else A1=rbind(cbind(P1,Q1),cbind(t(Q1),S1))
   return(A1)
}



subinv500=function(A){
   sys=1
   m=ncol(A)
   if(m<=20)return(ginv(A))
   else if(m<=50)return(subinv20(A))
   else if(m<=100)return(subinv50(A))
   else if(m<=200)return(subinv100(A))
   else if(m<=500)return(subinv200(A))
   for(i in 1:(m-1)){
       for(j in (i+1):m){
          if(A[i,j]!=A[j,i]){
             sys=0
             break
          }
       }
       if(sys==0)break
   }
   p=as.integer(m/2)
   s=m-p
   P=A[1:p,1:p]
   R=A[(p+1):m,1:p]
   Q=A[1:p,(p+1):m]
   S=A[(p+1):m,(p+1):m]
    
   if(p<=20){
      S10=ginv(S)
      P1=ginv(P-Q%*%S10%*%R)
   }
   else if(p>20&&p<=50){
      S10=subinv20(S)
      P1=subinv20(P-Q%*%S10%*%R)
   }
     
   else if(p>50&&p<=100){
     S10=subinv50(S)
     P1=subinv50(P-Q%*%S10%*%R)
   }
   else if(p>100&&p<=200){
     S10=subinv100(S)
     P1=subinv100(P-Q%*%S10%*%R)
   }
   else if(p>200){
     S10=subinv200(S)
     P1=subinv200(P-Q%*%S10%*%R)
   }

 
   Q1=-P1%*%Q%*%S10
   if(sys==0)R1=-(S10%*%R)%*%P1
   S1=S10+(S10%*%R)%*%P1%*%Q%*%S10
   A1=cbind(P1,Q1)
   if(sys==0)A1=rbind(cbind(P1,Q1),cbind(R1,S1))
   else A1=rbind(cbind(P1,Q1),cbind(t(Q1),S1))
   return(A1)
}

subinv200=function(A){
   sys=1
   m=ncol(A)
   if(m<=20)return(ginv(A))
   else if(m<=50)return(subinv20(A))
   else if(m<=100)return(subinv50(A))
   else if(m<=200)return(subinv100(A))
   for(i in 1:(m-1)){
       for(j in (i+1):m){
          if(A[i,j]!=A[j,i]){
             sys=0
             break
          }
       }
       if(sys==0)break
   }

   p=as.integer(m/2)
   s=m-p
   P=A[1:p,1:p]
   R=A[(p+1):m,1:p]
   Q=A[1:p,(p+1):m]
   S=A[(p+1):m,(p+1):m]
    
   if(p<=20){
      S10=ginv(S)
      P1=ginv(P-Q%*%S10%*%R)
   }
   else if(p>20&&p<=50){
      S10=subinv20(S)
      P1=subinv20(P-Q%*%S10%*%R)
   }
     
   else if(p>50&&p<=100){
     S10=subinv50(S)
     P1=subinv50(P-Q%*%S10%*%R)
   }
   else if(p>100){
     S10=subinv100(S)
     P1=subinv100(P-Q%*%S10%*%R)
   }
 
   Q1=-P1%*%Q%*%S10
   if(sys==0)R1=-(S10%*%R)%*%P1
   S1=S10+(S10%*%R)%*%P1%*%Q%*%S10
   A1=cbind(P1,Q1)
   if(sys==0)A1=rbind(cbind(P1,Q1),cbind(R1,S1))
   else A1=rbind(cbind(P1,Q1),cbind(t(Q1),S1))
   return(A1)

}

#A=vi
subinv100=function(A){
   sys=1
   m=ncol(A)
   if(m<=20)return(ginv(A))
   else if(m<=50)return(subinv20(A))
   else if(m<=100)return(subinv50(A))
   for(i in 1:(m-1)){
       for(j in (i+1):m){
          if(A[i,j]!=A[j,i]){
             sys=0
             break
          }
       }
       if(sys==0)break
   }

   p=as.integer(m/2)
   s=m-p
   P=A[1:p,1:p]
   R=A[(p+1):m,1:p]
   Q=A[1:p,(p+1):m]
   S=A[(p+1):m,(p+1):m]
    
   if(p<=20){
      S10=ginv(S)
      P1=ginv(P-Q%*%S10%*%R)
   }
   else if(p>20&&p<=50){
      S10=subinv20(S)
      P1=subinv20(P-Q%*%S10%*%R)
   }
     
   else{
     S10=subinv50(S)
     P1=subinv50(P-Q%*%S10%*%R)
   }
   Q1=-P1%*%Q%*%S10
   if(sys==0)R1=-(S10%*%R)%*%P1
   S1=S10+(S10%*%R)%*%P1%*%Q%*%S10
   A1=cbind(P1,Q1)
   if(sys==0)A1=rbind(cbind(P1,Q1),cbind(R1,S1))
   else A1=rbind(cbind(P1,Q1),cbind(t(Q1),S1))
   return(A1)

}


subinv50=function(A){
   sys=1
   m=ncol(A)
   if(m<=20)return(ginv(A))
   else if(m<=50)return(subinv20(A))
   for(i in 1:(m-1)){
       for(j in (i+1):m){
          if(A[i,j]!=A[j,i]){
             sys=0
             break
          }
       }
       if(sys==0)break
   }

   p=as.integer(m/2)
   s=m-p
   P=A[1:p,1:p]
   R=A[(p+1):m,1:p]
   Q=A[1:p,(p+1):m]
   S=A[(p+1):m,(p+1):m]
    
   if(p<=4){
      S10=ginv(S)
      P1=ginv(P-Q%*%S10%*%R)
   }
   else{
     S10=subinv20(S)
     P1=subinv20(P-Q%*%S10%*%R)
   }
  
  
   Q1=-P1%*%Q%*%S10
   if(sys==0)R1=-(S10%*%R)%*%P1
   S1=S10+(S10%*%R)%*%P1%*%Q%*%S10
   A1=cbind(P1,Q1)
   if(sys==0)A1=rbind(cbind(P1,Q1),cbind(R1,S1))
   else A1=rbind(cbind(P1,Q1),cbind(t(Q1),S1))
   return(A1)

}

subinv20=function(A){
    sys=1
    m=ncol(A)
    if(m<=20)return(ginv(A))
    for(i in 1:(m-1)){
       for(j in (i+1):m){
          if(A[i,j]!=A[j,i]){
             sys=0
             break
          }
       }
       if(sys==0)break
    }

    p=as.integer(m/2)
    s=m-p
    P=A[1:p,1:p]
    R=A[(p+1):m,1:p]
    Q=A[1:p,(p+1):m]
    S=A[(p+1):m,(p+1):m]
    
    S10=ginv(S)
    P1=ginv(P-Q%*%S10%*%R)
    Q1=-P1%*%Q%*%S10
    if(sys==0)R1=-(S10%*%R)%*%P1
    S1=S10+(S10%*%R)%*%P1%*%Q%*%S10
    A1=cbind(P1,Q1)
    if(sys==0)A1=rbind(cbind(P1,Q1),cbind(R1,S1))
    else A1=rbind(cbind(P1,Q1),cbind(t(Q1),S1))
    return(A1)
}

PedRecode=function(Ped,p,...){
   T1=NULL
   for(i in 2:(p+1))T1=c(T1,as.vector(Ped[,i]))
   T1=crecode(T1)$nv
   TM=matrix(T1,ncol=p,byrow=F)
   Ped[,2:(p+1)]=TM
   if(p==2)colnames(Ped)[2:3]=c("F","M")
   else colnames(Ped)[2:5]=c("F1","M1","F2","M2")
   return(Ped)
}
