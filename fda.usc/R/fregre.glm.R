fregre.glm=function(formula,family = gaussian(),data,basis.x=NULL,
basis.b=NULL,CV=FALSE,...){
 tf <- terms.formula(formula)
 terms <- attr(tf, "term.labels")
 nt <- length(terms)
 if (attr(tf, "response") > 0) {
        response <- as.character(attr(tf, "variables")[2])
        pf <- rf <- paste(response, "~", sep = "")
    } else pf <- rf <- "~"
 vtab<-rownames(attr(tf,"factors"))
 vnf=intersect(terms,names(data$df))
# vnf2=intersect(vtab[-1],names(data$df)[-1])
 vfunc2=setdiff(terms,vnf)
 vint=setdiff(terms,vtab)
 vfunc=setdiff(vfunc2,vint)
# vnf=c(vnf2,vint)
 off<-attr(tf,"offset")
 name.coef=nam=par.fregre=beta.l=list()
 kterms=1
#print(response)
 if (length(vnf)>0) {
 XX=data[["df"]][,c(response,vnf)] #data.frame el 1er elemento de la lista
 for ( i in 1:length(vnf)){
#     print(paste("Non functional covariate:",vnf[i]))
     if (kterms > 1)   pf <- paste(pf, "+", vnf[i], sep = "")
     else pf <- paste(pf, vnf[i], sep = "")
     kterms <- kterms + 1
     }
if   (attr(tf,"intercept")==0) {
     pf<- paste(pf,-1,sep="")
     }
}
else {
 XX=data$df[response]
 names(XX)=response
}
#print(paste("Non functional covariate:",vnf))
#print(paste("Functional covariate:",vfunc))
 bsp1<-bsp2<-TRUE
 lenfunc<-length(vfunc)>0
 if (lenfunc) {
  mean.list=vs.list=JJ=list()
  for (i in 1:length(vfunc)) {
    if (class(data[[vfunc[i]]])[1]=="fdata"){
      tt<-data[[vfunc[i]]][["argvals"]]
      rtt<-data[[vfunc[i]]][["rangeval"]]
      fdat<-data[[vfunc[i]]];      dat<-data[[vfunc[i]]]$data
      if (is.null(basis.x[[vfunc[i]]]))  basis.x[[vfunc[i]]]<-create.fdata.basis(fdat,l=1:7)
      else   if (basis.x[[vfunc[i]]]$type=="pc" | basis.x[[vfunc[i]]]$type=="pls") bsp1=FALSE
      if (is.null(basis.b[[vfunc[i]]])& bsp1)  basis.b[[vfunc[i]]]<-create.fdata.basis(fdat)
      else           if (class(basis.x[[vfunc[i]]])=="fdata.comp" | basis.x[[vfunc[i]]]$type=="pls") bsp2=FALSE
      if (bsp1 & bsp2) {
          if (is.null(rownames(dat)))    rownames(fdat$data)<-1:nrow(dat)
          fdnames=list("time"=tt,"reps"=rownames(fdat[["data"]]),"values"="values")
          xcc<-fdata.cen(data[[vfunc[i]]])
          mean.list[[vfunc[i]]]=xcc[[2]]
          if (!is.null( basis.x[[vfunc[i]]]$dropind)) {
              int<-setdiff(1:basis.x[[vfunc[i]]]$nbasis,basis.x[[vfunc[i]]]$dropind)
              basis.x[[vfunc[i]]]$nbasis<-length(int)
              basis.x[[vfunc[i]]]$dropind<-NULL
              basis.x[[vfunc[i]]]$names<-basis.x[[vfunc[i]]]$names[int]
              }
          if (!is.null( basis.b[[vfunc[i]]]$dropind)) {
              int<-setdiff(1:basis.b[[vfunc[i]]]$nbasis,basis.b[[vfunc[i]]]$dropind)
              basis.b[[vfunc[i]]]$nbasis<-length(int)
              basis.b[[vfunc[i]]]$dropind<-NULL
              basis.b[[vfunc[i]]]$names<-basis.b[[vfunc[i]]]$names[int]
              }
        x.fd = Data2fd(argvals = tt, y = t(xcc[[1]]$data),basisobj = basis.x[[vfunc[i]]],fdnames=fdnames)
          r=x.fd[[2]][[3]]
          J=inprod(basis.x[[vfunc[i]]],basis.b[[vfunc[i]]])
          Z =t(x.fd$coefs) %*% J
          colnames(J)=colnames(Z) = name.coef[[vfunc[i]]]=paste(vfunc[i],".",basis.b[[vfunc[i]]]$names,sep="")
          XX = cbind(XX,Z)
          for ( j in 1:length(colnames(Z))){
           if (kterms >= 1)  pf <- paste(pf, "+", colnames(Z)[j], sep = "")
           else pf <- paste(pf, colnames(Z)[j], sep = "")
           kterms <- kterms + 1
           }
        JJ[[vfunc[i]]]<-J
    }
      else {   
        l<-basis.x[[vfunc[i]]]$l
        lenl<-length(l)
        vs <- t(basis.x[[vfunc[i]]]$basis$data)
        Z<-basis.x[[vfunc[i]]]$x[,l,drop=FALSE]     
        response = "y"
        colnames(Z) = name.coef[[vfunc[i]]]=paste(vfunc[i], ".",colnames(Z),sep ="")      
#       colnames(Z) = name.coef[[vfunc[i]]]=paste(vfunc[i],".",rownames(basis.x[[vfunc[i]]]$basis$data),sep ="")
        name.coef[[vfunc[i]]]<-colnames(Z)
        XX = cbind(XX,Z)
        vs.list[[vfunc[i]]]=basis.x[[vfunc[i]]]$basis
        mean.list[[vfunc[i]]]=basis.x[[vfunc[i]]]$mean
        for ( j in 1:length(colnames(Z))){
            if (kterms >= 1)  pf <- paste(pf, "+", name.coef[[vfunc[i]]][j], sep = "")
           else pf <- paste(pf, name.coef[[vfunc[i]]][j], sep = "")
           kterms <- kterms + 1
           }       
          }
    }
 else {
 if(class(data[[vfunc[i]]])[1]=="fd"){
      fdat<-data[[vfunc[i]]]
      if (is.null(basis.x[[vfunc[i]]]))  basis.x[[vfunc[i]]]<-fdat$basis
      else   if (class(basis.x[[vfunc[i]]])=="pca.fd") bsp1=FALSE
      if (is.null(basis.b[[vfunc[i]]])& bsp1)
         basis.b[[vfunc[i]]]<-create.fdata.basis(fdat,
         l=1:max(5,floor(basis.x[[vfunc[i]]]$nbasis/5)),type.basis=basis.x[[vfunc[i]]]$type,
         rangeval=fdat$basis$rangeval)
      else           if (class(basis.x[[vfunc[i]]])=="pca.fd") bsp2=FALSE
      if (bsp1 & bsp2) {
          r=fdat[[2]][[3]]
          if (!is.null( basis.x[[vfunc[i]]]$dropind)) {
              int<-setdiff(1:basis.x[[vfunc[i]]]$nbasis,basis.x[[vfunc[i]]]$dropind)
              basis.x[[vfunc[i]]]$nbasis<-length(int)
              basis.x[[vfunc[i]]]$dropind<-NULL
              basis.x[[vfunc[i]]]$names<-basis.x[[vfunc[i]]]$names[int]
              }
          if (!is.null( basis.b[[vfunc[i]]]$dropind)) {
              int<-setdiff(1:basis.b[[vfunc[i]]]$nbasis,basis.b[[vfunc[i]]]$dropind)
              basis.b[[vfunc[i]]]$nbasis<-length(int)
              basis.b[[vfunc[i]]]$dropind<-NULL
              basis.b[[vfunc[i]]]$names<-basis.b[[vfunc[i]]]$names[int]
              }
          J=inprod(basis.x[[vfunc[i]]],basis.b[[vfunc[i]]])
          mean.list[[vfunc[i]]]<-mean.fd(x.fd)
          x.fd<-center.fd(x.fd)
          Z =t(x.fd$coefs) %*% J
          colnames(J)=colnames(Z) = name.coef[[vfunc[i]]]=paste(vfunc[i],".",basis.b[[vfunc[i]]]$names,sep="")
          XX = cbind(XX,Z)
          for ( j in 1:length(colnames(Z))){
           if (kterms >= 1)  pf <- paste(pf, "+", colnames(Z)[j], sep = "")
           else pf <- paste(pf, colnames(Z)[j], sep = "")
           kterms <- kterms + 1
           }
        JJ[[vfunc[i]]]<-J
}
      else {
        l<-ncol(basis.x[[vfunc[i]]]$scores)
        vs <- basis.x[[vfunc[i]]]$harmonics$coefs
        Z<-basis.x[[vfunc[i]]]$scores
        response = "y"
        colnames(Z) = name.coef[[vfunc[i]]]=paste(vfunc[i], ".",colnames(basis.x[[vfunc[i]]]$harmonics$coefs),sep ="")
        XX = cbind(XX,Z)
        vs.list[[vfunc[i]]]=vs
        mean.list[[vfunc[i]]]=basis.x[[vfunc[i]]]$meanfd
        for ( j in 1:length(colnames(Z))){
            if (kterms >= 1)  pf <- paste(pf, "+", name.coef[[vfunc[i]]][j], sep = "")
           else pf <- paste(pf, name.coef[[vfunc[i]]][j], sep = "")
           kterms <- kterms + 1
           }
         }
    }
   else stop("Please, enter functional covariate")
   }
  }  }
####################################
    if (!is.data.frame(XX)) XX=data.frame(XX)
    par.fregre$formula=pf
    par.fregre$data=XX
    ndatos<-nrow(XX)
    yp<-rep(NA,ndatos)
    if (CV) {
     for (k in 1:ndatos) {
      z=glm(formula=pf,data=XX[-k,],family=family,x=TRUE,y=TRUE,...)
      for (i in 1:length(vfunc)) {

      if (bsp1) beta.est=fdata(fd(z[[1]][name.coef[[vfunc[i]]]],basis.b[[vfunc[i]]]),tt)
      else{
       if (class(data[[vfunc[i]]])[1]=="fdata"){
            beta.est<-z$coefficients[name.coef[[vfunc[i]]]]*vs.list[[vfunc[i]]]
#            beta.est$data<-apply(beta.est$data,2,sum)
            beta.est$data<-colSums(beta.est$data)
            beta.est$names$main<-"beta.est"
            beta.est$data <- matrix(as.numeric(beta.est$data),nrow=1)
            }
       else {
            beta.est<-z$coefficients[name.coef[[vfunc[i]]]]*t(vs.list[[vfunc[i]]])
#            beta.est<-apply(beta.est,2,sum)
            beta.est<-colSums(beta.est)
            beta.est<-fdata(fd(beta.est,basis.x[[vfunc[i]]]$harmonics$basis),tt)
       }
      }
      if (k==1) beta.l[[vfunc[i]]]<-beta.est
      else beta.l[[vfunc[i]]]<-c(beta.l[[vfunc[i]]],beta.est)
      }
      yp[k]=predict.glm(object=z,newdata =XX[k,],type="response",x=TRUE,y=TRUE,...)
     }
        z=glm(formula=pf,data=XX,family=family,x=TRUE,y=TRUE,...)
        dev.resids<-family$dev.resids(XX[[response]],yp,1) #sum(dev.resids=z$deviance)
        z$CV<-list("y.pred"=yp,"beta.l"=beta.l,"dev.resids"=dev.resids)
    }
    else {    
      if (length(vfunc)==0 & length(vnf)==0)      {
        pf<-as.formula(paste(pf,1,sep=""))
        z=glm(formula=pf,data=XX,family=family,x=TRUE,y=TRUE,...)
        class(z)<-c(class(z),"fregre.glm")
        z$formula.ini=pf
        z$XX=XX
        z$data<-data
        return(z )
       }
      else    z=glm(formula=pf,data=XX,family=family,x=TRUE,y=TRUE,...)
      
     }

    z$call<-z$call[1:2]
if (lenfunc) {    
    for (i in 1:length(vfunc)) {
      if (bsp1) beta.l[[vfunc[i]]]=fd(z[[1]][name.coef[[vfunc[i]]]],basis.b[[vfunc[i]]])
       else{
       if (class(data[[vfunc[i]]])[1]=="fdata"){
            beta.est<-z$coefficients[name.coef[[vfunc[i]]]]*vs.list[[vfunc[i]]]
#            beta.est$data<-apply(beta.est$data,2,sum)
            beta.est$data<-colSums(beta.est$data)
            beta.est$names$main<-"beta.est"
            beta.est$data <- matrix(as.numeric(beta.est$data),nrow=1)
            if  (basis.x[[vfunc[i]]]$type=="pls") {
             if (basis.x[[vfunc[i]]]$norm)  {
              sd.X <- sqrt(apply(data[[vfunc[i]]]$data, 2, var))
              beta.est$data<-  beta.est$data/sd.X
             }      
            }  
            beta.l[[vfunc[i]]]<-beta.est
            }
       else {
            beta.est<-z$coefficients[name.coef[[vfunc[i]]]]*t(vs.list[[vfunc[i]]])
#            beta.est<-apply(beta.est,2,sum)
            beta.est<-colSums(beta.est)
            beta.l[[vfunc[i]]]<-fd(beta.est,basis.x[[vfunc[i]]]$harmonics$basis)
      }
      }
    }
# z$sr2<-sum(z$residuals^2)/z$df.residual
 z$beta.l=beta.l
 z$basis.x=basis.x
 z$basis.b=basis.b
 z$mean=mean.list
 z$formula=pf
 z$JJ<-JJ
 z$vs.list=vs.list
 }
 z$formula.ini=formula
 z$data=z$data
 z$XX=XX
 class(z)<-c(class(z),"fregre.glm")
 z
}

  fregre.lm<-function(formula,data,basis.x=NULL,basis.b=NULL,
rn,lambda,weights=rep(1,n),...){
 tf <- terms.formula(formula)
 terms <- attr(tf, "term.labels")
 nt <- length(terms)
 if (attr(tf, "response") > 0) {
        response <- as.character(attr(tf, "variables")[2])
        pf <- rf <- paste(response, "~", sep = "")
    } else pf <- rf <- "~"
 vtab<-rownames(attr(tf,"factors"))                                                    
 vnf=intersect(terms,names(data$df))
# vnf2=intersect(vtab[-1],names(data$df)[-1])
 vfunc2=setdiff(terms,vnf)
 vint=setdiff(terms,vtab)
 vfunc=setdiff(vfunc2,vint)
# vnf=c(vnf2,vint)
 off<-attr(tf,"offset")
 name.coef=nam=par.fregre=beta.l=list()
 kterms=1
 n<-length(data[["df"]][,response])
 XX=data.frame(data[["df"]][,c(response)],weights)
 namxx=names(XX)=c(response,"weights")
 lenvnf<-length(vnf)
 df<-0
 if (lenvnf>0) { 
   df<-lenvnf +df         
#print(paste("Functional covariate:",vnf))
# XX=data.frame(XX,data[["df"]][,c(vnf)])    # si es facotr poner k-1 variables dummies
# names(XX)=c(namxx,vnf)
  for ( i in 1:lenvnf){
#     print(paste("Non functional covariate:",vnf[i]))
     if (kterms > 1)   pf <- paste(pf, "+", vnf[i], sep = "")
#     else pf <- paste(pf, terms[i], sep = "")
     else pf <- paste(pf, vnf[i], sep = "")
     kterms <- kterms + 1
     }
if   (attr(tf,"intercept")==0) {
#     print("No intecept")
     pf<- paste(pf,-1,sep="")
     df<-0
     }
else df<-1     
#  contrasts = NULL
    mf<-as.data.frame(model.matrix(formula(pf),data$df))
#    vnf2<-names(mf)[-respo]
    pf <- rf <- paste(response, "~", sep = "")
    for ( i in 1:length(vnf))  pf<-paste(pf, "+", vnf[i], sep = "")
    XX <- data.frame(XX,mf)
}
else {
    XX <- data.frame(XX,model.matrix(formula(paste(pf, "1")),data$df))
    names(XX)[3]<-"(Intercept)"
    df<-1
}
#print(paste("Functional covariate:",vfunc))
if (missing(rn))    {    rn0=FALSE;                    rn=list()}
else rn0<-TRUE
if (missing(lambda))    {    lambda0=FALSE;                    lambda=list()}
else lambda0<-TRUE
mat<-rep(0,len=ncol(XX)-2)
imat2<-ncol(XX)-2
mat2<-diag(0,nrow=imat2)
lenvfunc<-length(vfunc)
hay.pls<-FALSE
if (lenvfunc>0) {
 mean.list=vs.list=JJ=list()
 bsp1<-bsp2<-TRUE
 for (i in 1:lenvfunc) {
	if(class(data[[vfunc[i]]])[1]=="fdata"){
      tt<-data[[vfunc[i]]][["argvals"]]
      rtt<-data[[vfunc[i]]][["rangeval"]]
      np<-length(tt)
      fdat<-data[[vfunc[i]]];      dat<-data[[vfunc[i]]]$data
      if (is.null(basis.x[[vfunc[i]]]))  basis.x[[vfunc[i]]]<-create.fdata.basis(fdat,l=1:7)
      else   if (basis.x[[vfunc[i]]]$type=="pc" | basis.x[[vfunc[i]]]$type=="pls") bsp1=FALSE
      if (is.null(basis.b[[vfunc[i]]])& bsp1)  basis.b[[vfunc[i]]]<-create.fdata.basis(fdat)
      else           if (basis.x[[vfunc[i]]]$type=="pc" | basis.x[[vfunc[i]]]$type=="pls") bsp2=FALSE
      if (bsp1 & bsp2) {
          if (is.null(rownames(dat)))    rownames(fdat$data)<-1:nrow(dat)
          fdnames=list("time"=tt,"reps"=rownames(fdat[["data"]]),"values"="values")
          xcc<-fdata.cen(data[[vfunc[i]]])
          mean.list[[vfunc[i]]]=xcc[[2]]
          if (!is.null( basis.x[[vfunc[i]]]$dropind)) {
              int<-setdiff(1:basis.x[[vfunc[i]]]$nbasis,basis.x[[vfunc[i]]]$dropind)
              basis.x[[vfunc[i]]]$nbasis<-length(int)
              basis.x[[vfunc[i]]]$dropind<-NULL
              basis.x[[vfunc[i]]]$names<-basis.x[[vfunc[i]]]$names[int]
              }
#          else int<-1:basis.x[[vfunc[i]]]$nbasis
#                  nam[[vfunc[i]]]<-  basis.x[[vfunc[i]]]$names
          if (!is.null( basis.b[[vfunc[i]]]$dropind)) {
              int<-setdiff(1:basis.b[[vfunc[i]]]$nbasis,basis.b[[vfunc[i]]]$dropind)
              basis.b[[vfunc[i]]]$nbasis<-length(int)
              basis.b[[vfunc[i]]]$dropind<-NULL
              basis.b[[vfunc[i]]]$names<-basis.b[[vfunc[i]]]$names[int]
              }
    	    x.fd = Data2fd(argvals = tt, y = t(xcc[[1]]$data),basisobj = basis.x[[vfunc[i]]],fdnames=fdnames)
          r=x.fd[[2]][[3]]
#          Theta = getbasismatrix(tt,basis.b[[vfunc[i]]])
#          Psi = getbasismatrix(tt,basis.x[[vfunc[i]]])
#          J = t(Psi) %*% Theta
          J=inprod(basis.x[[vfunc[i]]],basis.b[[vfunc[i]]])
          Z =t(x.fd$coefs) %*% J
          colnames(J)=colnames(Z) = name.coef[[vfunc[i]]]=paste(vfunc[i],".",basis.b[[vfunc[i]]]$names,sep="")
          XX = data.frame(XX,Z)
          for ( j in 1:length(colnames(Z))){
           if (kterms >= 1)  pf <- paste(pf, "+", colnames(Z)[j], sep = "")
           else pf <- paste(pf, colnames(Z)[j], sep = "")
           kterms <- kterms + 1
           }
        	JJ[[vfunc[i]]]<-J
           nbasisb<- basis.b[[vfunc[i]]]$nbasis
         R=diag(0,ncol= nbasisb,nrow=nbasisb)
         R<-eval.penalty(basis.b[[vfunc[i]]],Lfdobj = int2Lfd(2))
        mat2<-diag(0,nrow=imat2+nbasisb)
        if (!is.null(lambda[[vfunc[i]]]))  {
        mat2<-diag(0,nrow=imat2+nbasisb)
        mat2[(imat2+1):(imat2+nbasisb),(imat2+1):(imat2+nbasisb)]<-lambda[[vfunc[i]]]*R
        imat2<-imat2+nbasisb
        }     
				}
      else {
        basis<-basis.x[[vfunc[i]]]
        l<-basis$l
        vs <- t(basis$basis$data)           
        basis$x<-basis$x[,l,drop=FALSE]
        Z<-basis$x
        response = "y"
        colnames(Z) = name.coef[[vfunc[i]]]=paste(vfunc[i], ".",rownames(basis$basis$data),sep ="")
        XX =data.frame(XX,Z)
        vs.list[[vfunc[i]]]=basis$basis
        mean.list[[vfunc[i]]]=basis$mean
        for ( j in 1:length(colnames(Z))){
            if (kterms >= 1)  pf <- paste(pf, "+", name.coef[[vfunc[i]]][j], sep = "")
           else pf <- paste(pf, name.coef[[vfunc[i]]][j], sep = "")
           kterms <- kterms + 1
           }
        if (!is.null(rn[[vfunc[i]]]))  {
           mat<-c(mat,rep(rn[[vfunc[i]]],len=length(l)))
           }     
        else          mat<-c(mat,rep(0,len=length(l)))
        mat2<-diag(0,nrow=imat2+length(l))
        if (!is.null(lambda[[vfunc[i]]]))  {
        mat2<-diag(0,nrow=imat2+length(l))
        R<-P.penalty(1:length(l))
        mat2[(imat2+1):(imat2+length(l)),(imat2+1):(imat2+length(l))]<-R
        imat2<-imat2+length(l)
        }                    
   }
  }
 	else {
 		if(class(data[[vfunc[i]]])[1]=="fd"){
      fdat<-data[[vfunc[i]]]
      if (is.null(basis.x[[vfunc[i]]]))  basis.x[[vfunc[i]]]<-fdat$basis
      else   if (class(basis.x[[vfunc[i]]])=="pca.fd") bsp1=FALSE
      if (is.null(basis.b[[vfunc[i]]])& bsp1)
         basis.b[[vfunc[i]]]<-create.fdata.basis(fdat,
         l=1:max(5,floor(basis.x[[vfunc[i]]]$nbasis/5)),type.basis=basis.x[[vfunc[i]]]$type,
         rangeval=fdat$basis$rangeval)
      else           if (class(basis.x[[vfunc[i]]])=="pca.fd") bsp2=FALSE
      if (bsp1 & bsp2) {
          r=fdat[[2]][[3]]
#          tt = seq(r[1], r[2], len = length(fdat[[3]]$time))
          if (!is.null( basis.x[[vfunc[i]]]$dropind)) {
              int<-setdiff(1:basis.x[[vfunc[i]]]$nbasis,basis.x[[vfunc[i]]]$dropind)
              basis.x[[vfunc[i]]]$nbasis<-length(int)
              basis.x[[vfunc[i]]]$dropind<-NULL
              basis.x[[vfunc[i]]]$names<-basis.x[[vfunc[i]]]$names[int]
              }
          if (!is.null( basis.b[[vfunc[i]]]$dropind)) {
              int<-setdiff(1:basis.b[[vfunc[i]]]$nbasis,basis.b[[vfunc[i]]]$dropind)
              basis.b[[vfunc[i]]]$nbasis<-length(int)
              basis.b[[vfunc[i]]]$dropind<-NULL
              basis.b[[vfunc[i]]]$names<-basis.b[[vfunc[i]]]$names[int]
              }
          J=inprod(basis.x[[vfunc[i]]],basis.b[[vfunc[i]]])
          mean.list[[vfunc[i]]]<-mean.fd(x.fd)
          x.fd<-center.fd(x.fd)
          Z =t(x.fd$coefs) %*% J
          colnames(J)=colnames(Z) = name.coef[[vfunc[i]]]=paste(vfunc[i],".",basis.b[[vfunc[i]]]$names,sep="")
          XX = data.frame(XX,Z)
          for ( j in 1:length(colnames(Z))){
           if (kterms >= 1)  pf <- paste(pf, "+", colnames(Z)[j], sep = "")
           else pf <- paste(pf, colnames(Z)[j], sep = "")
           kterms <- kterms + 1
           }
        	JJ[[vfunc[i]]]<-J
				}
      else {
        l<-ncol(basis.x[[vfunc[i]]]$scores)
        vs <- basis.x[[vfunc[i]]]$harmonics$coefs
        Z<-basis.x[[vfunc[i]]]$scores
        response = "y"
        colnames(Z) = name.coef[[vfunc[i]]]=paste(vfunc[i], ".",colnames(basis.x[[vfunc[i]]]$harmonics$coefs),sep ="")
        XX = data.frame(XX,Z)
        vs.list[[vfunc[i]]]=vs
        mean.list[[vfunc[i]]]=basis.x[[vfunc[i]]]$meanfd
        for ( j in 1:length(colnames(Z))){
            if (kterms >= 1)  pf <- paste(pf, "+", name.coef[[vfunc[i]]][j], sep = "")
           else pf <- paste(pf, name.coef[[vfunc[i]]][j], sep = "")
           kterms <- kterms + 1
           }
         }
    }
   else stop("Please, enter functional covariate")
   }
  }  
if    (basis.x[[vfunc[i]]]$type=="pls") {
   hay.pls<-TRUE
   rn0<-TRUE # SI ES Pls ASI TIENE EN CUENTA LOS DOF QUE EN EL LM() NO LO TENDRIA
   lenl<-ncol(basis.x[[vfunc[i]]]$x)
   df<-df+basis.x[[vfunc[i]]]$df[lenl]
   }
}       
  
if (!is.data.frame(XX)) XX=data.frame(XX)
    par.fregre$formula=pf
    par.fregre$data=XX
    y<-XX[,1]    
    scores<-as.matrix(XX[,-(1:2)])     
#    scores<-as.matrix(XX)     
    W<-diag(weights)  
    if (!rn0 & !lambda0) {
      if (lenvfunc==0 & length(vnf)==0)      {
      z=lm(formula=paste(pf,1,sep=""),data=XX,x=TRUE,...)   
       class(z)<-c("lm","fregre.lm")
      return(z )
      }
      else       z=lm(formula=pf,data=XX,x=TRUE,...)
      e<-z$residuals
#      S<-solve(t(scores)%*%W%*%scores)          
       S<-diag(coef(summary(z))[,2])
      class(z)<-c(class(z),"fregre.lm")
      }      
    else {
       if (lambda0) { S<-solve(t(scores)%*%W%*%scores+mat2)           }             
       if (rn0) {
             mat<-diag(mat)
             S<-solve(t(scores)%*%W%*%scores+mat)       
       }      
       if (rn0 & lambda0) warning("Only ridge penalization is done by rn argument (lambda argument is ignored)")               
#      ddd<-t(scores)%*%W  
#      S<-solve(t(scores)%*%W%*%scores+mat+mat2)       #incluir pesos solve(W)
      Cinv<-S%*%t(scores)%*%W                    #incluir pesos W repetri proceso hasta que no cambie la prediccion   
      ycen = y - mean(y)
      coefs<-Cinv%*%XX[,1]
      z<-list()      
      z$fitted.values<-yp<-drop(scores%*%coefs)
      e<-z$residuals<-XX[,1]- z$fitted.values      
################################################################################
################################################################################
################################################################################           
      H<-scores%*%Cinv
      if  (!hay.pls) df<-traza(H)  
      coefs<-drop(coefs)
#      cnames<-names(XX[,-(1:2)])
#      names(coefs)<-c("Intercept",cnames)  
      z$coefficients<-coefs
      z$mean.list<-mean.list
      z$df.residual<-n-df
      z$H<-H
      z$r2 <- 1 - sum(z$residuals^2)/sum(ycen^2)       
      if  (class(basis.x[[vfunc[1]]])=="basisfd") {
        z$call[[1]] = "fregre.basis"
        }
       else {
        if  (basis.x[[vfunc[1]]]$type=="pc")  z$call[[1]] = "fregre.pc"
        if  (basis.x[[vfunc[1]]]$type=="pls")  z$call[[1]] = "fregre.pls"        
        }             
    class(z)<-c("fregre.fd","fregre.lm")
    rdf<-n-df
    sr2 <- sum(e^2)/ rdf
    r2 <- 1 - sum(e^2)/sum(ycen^2)
    r2.adj<- 1 - (1 - r2) * ((n -    1)/ rdf)
    GCV <- sum(e^2)/(n - df)^2
        z$residuals <- drop(e)
        z$fitted.values <- yp
        z$y <- y
        z$rank <- df
        z$df.residual <-  rdf
        Z=cbind(rep(1,len=n),Z)
        colnames(Z)[1] = "(Intercept)"
        std.error = sqrt(diag(S) *sr2)
        t.value = coefs/std.error
        p.value = 2 * pt(abs(t.value), n - df, lower.tail = FALSE)
        coefficients <- cbind(coefs, std.error, t.value, p.value)
        colnames(coefficients) <- c("Estimate", "Std. Error",
            "t value", "Pr(>|t|)")
        z$coefs<-coefficients    
        class(z) <- "lm"
}       
#    z$call<-z$call[1:2]
for (i in 1:length(vfunc)) {
 if (bsp1) beta.l[[vfunc[i]]]=fd(z[[1]][name.coef[[vfunc[i]]]],basis.b[[vfunc[i]]])
 else{
	if(class(data[[vfunc[i]]])[1]=="fdata"){
#     beta.est<-z$coefficients[name.coef[[vfunc[i]]]]*vs.list[[vfunc[i]]]
     beta.est<-z$coefficients[name.coef[[vfunc[i]]]]*vs.list[[vfunc[i]]]
     beta.est$data<-colSums(beta.est$data)
     beta.est$names$main<-"beta.est"
     beta.est$data <- matrix(as.numeric(beta.est$data),nrow=1)
     beta.est$names$main<-"beta.est"
     beta.est$data <- matrix(as.numeric(beta.est$data),nrow=1)
           if  (basis.x[[vfunc[i]]]$type=="pls") {
             if (basis.x[[vfunc[i]]]$norm)  {
              sd.X <- sqrt(apply(data[[vfunc[i]]]$data, 2, var))
              beta.est$data<-  beta.est$data/sd.X
             }      
            }  
 
     beta.l[[vfunc[i]]]<-beta.est     
     }
 else {
     beta.est<-z$coefficients[name.coef[[vfunc[i]]]]*t(vs.list[[vfunc[i]]])
#     beta.est<-apply(beta.est,2,sum)
     beta.est<-colSums(beta.est)
     beta.l[[vfunc[i]]]<-fd(beta.est,basis.x[[vfunc[i]]]$harmonics$basis)
      }
}
}
 z$sr2<-sum(e^2)/z$df.residual
 z$Vp=z$sr2*S
 z$beta.l=beta.l
 z$formula=pf
 z$mean=mean.list
 z$formula.ini=formula
 z$basis.x=basis.x
 z$basis.b=basis.b
 z$JJ<-JJ
 z$data=z$data
 z$XX=XX
 z$data<-data
 z$fdataobj<-data[[vfunc[1]]]
 z$rn<-rn0
 z$lambda<-lambda0
 z$vs.list=vs.list   
## z$ar<-aa
 class(z)<-c("lm","fregre.lm")
 z
}     
 


