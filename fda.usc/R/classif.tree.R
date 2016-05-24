classif.tree=function(formula,data,basis.x=NULL,basis.b=NULL,CV=FALSE,...){
 C<-match.call()
 a<-list()
 mf <- match.call(expand.dots = FALSE)
 m <- match(c("formula","data","basis.x","basis.b","CV"), names(mf), 0L)
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
newdata<-data
group<-newy<-y<-data$df[[response]]    
#data$df[[response]]<-as.numeric(y)
 if (length(vnf)>0) {
 XX=data[[1]][,c(response,vnf)] #data.frame el 1er elemento de la lista
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
 XX=data.frame(y)
 names(XX)=response
}
if (!is.factor(y)) y<-as.factor(y)
prob2<-prob<-ngroup<-nlevels(y)
ny<-levels(y)
n<-length(y)
 mean.list=vs.list=JJ=list()
 #print(paste("Functional covariate:",vfunc))
if (length(vfunc)>0) {  
 bsp1<-bsp2<-TRUE
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
        vs <- t(basis.x[[vfunc[i]]]$basis$data)
#        Z<-basis.x[[vfunc[i]]]$x
        Z<-basis.x[[vfunc[i]]]$x[,l,drop=FALSE]  
        response = "y"
#        colnames(Z) = name.coef[[vfunc[i]]]=paste(vfunc[i], ".",rownames(basis.x[[vfunc[i]]]$basis$data),sep ="")       
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
#   else stop("Please, enter functional covariate")
   }
  }  }
  
    if (!is.data.frame(XX)) XX=data.frame(XX)
    par.fregre$formula=pf
    par.fregre$data=XX
    ndatos<-nrow(XX)
    yp<-rep(NA,ndatos)
    if (CV) {
     for (k in 1:ndatos) {
      z<-rpart(formula=pf,data=XX[-k,],...)
      yp[k]<-predict(object=z, newdata = XX[k,], type = "class")
     }
      z<-rpart(formula=pf,data=XX,...)
      yp<-predict(object=z, newdata = XX, type = "class")
      z$CV<-list("y.pred"=yp)
    }
    else {
      z<-rpart(formula=pf,data=XX,...)
      yp<-predict(object=z, type = "class")
     }
z$call<-z$call[1:2]
prob.group<-NULL
#y<-XX[,names(XX)[1]]
yp<-factor(yp,levels=ny)
tab<-table(yp,group)
prob2<-prob<-ngroup<-nlevels(y)
for (i in 1:ngroup) {     prob[i]=tab[i,i]/sum(tab[,i])     }
max.prob=sum(yp==y)/n
names(prob)<-ny
output<-list(formula=formula,fit=z,C=C,prob.group=prob.group,prob.classification=prob,
max.prob=max.prob,group.est=yp,group=y,formula.ini=formula,formula=pf,m=m,mean=mean.list,
basis.x=basis.x,basis.b=basis.b,JJ=JJ,data=z$data,XX=XX,vs.list=vs.list)
 class(output)="classif"
output
}
