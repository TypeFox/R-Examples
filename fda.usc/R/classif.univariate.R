######################
#######################
classif.glm2boost=function(group,fdataobj,family=binomial(),basis.x=NULL,
basis.b=NULL,CV=FALSE,...){
C<-match.call()
a<-list()
mf <- match.call(expand.dots = FALSE)
m <- match(c("group","fdataobj","family","basis.x","basis.b","CV"), names(mf), 0L)
   if (is.fdata(fdataobj))   {
       dataf<-data.frame("y"=group)
       ldata<-list("df"=dataf,"X"=fdataobj)
       formula<-formula(y~X)       
       }
   else {
        if (is.matrix(fdataobj)) nms<-colnames(fdataobj)
        else nms<-names(fdataobj)        
#quita [ o $ del names
#for (i in 1:length(nms)) nms[i]<-unlist(strsplit(nms[i], "[$]"))[[1]]
#for (i in 1:length(nms)) nms[i]<-unlist(strsplit(nms[i], "[[]"))[[1]]
        dataf<-data.frame("y"=group,fdataobj)
        ldata<-list("df"=dataf)
        aaa<-paste(nms,collapse="+")
        formula<-formula(paste("y~",aaa,sep=""))      
     }
   newy<-y<-ldata$df$y
   if (!is.factor(y)) y<-as.factor(y)
   n<-length(y)
   newdata<-ldata
   ngroup<-nlevels(y)
   ny<-levels(y)
   prob<-rep(NA,ngroup)
       if (!is.null(basis.x)) basis.x=list("X"=basis.x)
       if (!is.null(basis.b)) basis.b=list("X"=basis.b)  
   #nlevels(y) cambiar todos
if (ngroup==2) {
#      ny<-as.numeric(names(table(y)))
      newy<-ifelse(y==ny[1],0,1)
      newdata$df$y<-newy
      a[[1]]<-fregre.glm(formula,data=newdata,family=binomial,basis.x=basis.x,
              basis.b=basis.b,CV=CV)                
      yest<-ifelse(a[[1]]$fitted.values<.5,ny[1],ny[2])
      tab<-table(yest,y)
      prob[1]=tab[1,1]/sum(tab[,1])
      dtab<-dim(tab)
      if (dtab[1]==dtab[2])    {
           prob[2]=tab[2,2]/sum(tab[,2])
           names(prob)<-ny
         }
      else prob[2]<-0
      prob.group<-a$fitted.values
      yest<-factor(yest,levels=ny)

   }
else {
#   ny<-as.numeric(names(table(y)))
   prob.group<-array(NA,dim=c(n,ngroup))
   colnames(prob.group)<-ny
   for (i in 1:ngroup) {
              newy<-ifelse(y==ny[i],0,1)
              newdata$df$y<-newy
              a[[i]]<-fregre.glm(formula,data=newdata,family=family,basis.x=basis.x,
              basis.b=basis.b,CV=CV)
              prob.group[,i]<-a[[i]]$fitted.values
              }
   yest<-ny[apply(prob.group,1,which.min)]
   yest<-factor(yest,levels=ny)
   tab<-table(yest,y)
   for (ii in 1:ngroup) {prob[ii]=tab[ii,ii]/sum(tab[,ii])}
   names(prob)<-ny
}
max.prob=sum(diag(tab))/sum(tab)
output<-list(fdataobj=fdataobj,group=y,group.est=as.factor(yest),
prob.classification=prob,prob.group=prob.group,C=C,m=m,max.prob=max.prob,fit=a
,formula=formula,data=newdata)
class(output)="classif"
return(output)
}
#############################
classif.gsam2boost=function(group,fdataobj,family=binomial,weights=NULL,
basis.x=NULL,basis.b=NULL,par.gsam=NULL,CV=FALSE,...){
C<-match.call()
a<-list()
mf <- match.call(expand.dots = FALSE)
m <- match(c( "group","fdataobj","family","basis.x","basis.b","par.gsam","CV"), names(mf), 0L)
numg=nlevels(as.factor(group)) 

  if (is.fdata(fdataobj))   {
     gsam<-C[[3]];
     if (is.null(par.gsam$func)&is.null(par.gsam$k))   {gsam<-paste("s(X)",sep="")}
     else  {
        if (!is.null(par.gsam$func)&!is.null(par.gsam$k))
            gsam<-paste(par.gsam$func,"(X,k=",par.gsam$k,")",sep="")
        else if (is.null(par.gsam$func)&!is.null(par.gsam$k))
            gsam<-paste("s(X,k=",par.gsam$k,")",sep="")
            else         if (!is.null(par.gsam$func)&is.null(par.gsam$k))
            gsam<-paste(par.gsam$func,"(X)",sep="")
           }
#      pf2<-as.formula(paste(C[[2]], "~", gsam,sep = ""))
      pf2<-as.formula(paste("y~", gsam,sep = ""))
      dataf<-data.frame("y"=group)
      ldata<-list("df"=dataf,"X"=fdataobj)
      X<-C[[3]]       
      }
   else {   
        if (is.null(par.gsam$func)) par.gsam$func<-"s"
        if (is.null(par.gsam$k)) par.gsam$k<--1        
        dataf<-data.frame("y"=group,fdataobj)
        ldata<-list("df"=dataf)
        if (is.matrix(fdataobj)) nms<-colnames(fdataobj)
        else nms<-names(fdataobj)               
        gsam<-paste("+",par.gsam$func,"(",nms,",k=",par.gsam$k,")",sep="",collapse="")
#      pf2<-as.formula(paste(C[[2]], "~", gsam,sep = ""))
        pf2<-as.formula(paste("y~", gsam,sep = ""))
              
#        if (is.null(par.gsam$formula)) {
#          aaa<-paste("s(",nms,")",collapse="+")
#          pf2<-formula(paste("y~",aaa,sep=""))
#          }
#        else pf2<-par.gsam$formula  
        X<-C[[3]]              
     }   
newy<-y<-ldata$df$y
if (!is.factor(y)) y<-as.factor(y)
n<-length(y);
newdata<-ldata

   ngroup<-nlevels(y)
   ny<-levels(y)
prob<-ngroup<-length(table(y))
if (!is.null(basis.x)) basis.x=list("X"=basis.x)
if (!is.null(basis.b)) basis.b=list("X"=basis.b)
formula<-pf2
if (ngroup==2) {
      ny<-as.numeric(names(table(y)))
      newy<-ifelse(y==ny[1],0,1)
      newdata$df$y<-newy
#      formula<-formula(paste("y~",gsam),sep="")
      a[[1]]<-fregre.gsam(formula,data=newdata,family=family,weights=weights,basis.x=basis.x,
      basis.b=basis.b,CV=CV,...)
      yest<-ifelse(a[[1]]$fitted.values<.5,ny[1],ny[2])
      yest<-factor(yest,levels=ny)
      tab<-table(yest,y)
      prob[1]=tab[1,1]/sum(tab[,1])
      dtab<-dim(tab)
      if (dtab[1]==dtab[2])    {
           prob[2]=tab[2,2]/sum(tab[,2])
           names(prob)<-ny
         }
      else prob[2]<-0
      prob.group<-a$fitted.values
      #devolver a mayores y estimada
   }
else {
   ny<-as.numeric(names(table(y)))
   prob.group<-array(NA,dim=c(n,ngroup))
   colnames(prob.group)<-ny
   for (i in 1:ngroup) {
              newy<-ifelse(y==ny[i],0,1)
              newdata$df$y<-newy
#              formula<-formula(paste("y~",gsam),sep="")
              a[[i]]<-fregre.gsam(formula,data=newdata,family=family,
              weights=weights,basis.x=basis.x,basis.b=basis.b,CV=CV,...)
              prob.group[,i]<-a[[i]]$fitted.values
            }
   yest<-ny[apply(prob.group,1,which.min)]######no sera which.max
   yest<-factor(yest,levels=ny)
   tab<-table(yest,y)
   for (i in 1:ngroup) {     prob[i]=tab[i,i]/sum(tab[,i])     }
   names(prob)<-ny
}
max.prob=sum(diag(tab))/sum(tab)
output<-list(fdataobj=fdataobj,group=y,group.est=as.factor(yest),
prob.classification=prob,prob.group=prob.group,C=C,m=m,max.prob=max.prob,fit=a,
formula=formula,data=newdata)
class(output)="classif"
return(output)
}

#######################
classif.tree2boost=function(group,fdataobj,basis.x=NULL,basis.b=NULL,...){   
if (!is.factor(group)) group<-as.factor(group)
   C<-match.call()
   mf <- match.call(expand.dots = FALSE)
   m <- match(c("group","fdataobj","basis.x"), names(mf), 0L)

#  if (is.null(basis.x))  {
#   basis.x2<-create.pc.basis(fdataobj,l=1:5)
#   X<-basis.x<- basis.x2$x
#   }
#else {
#    basis.x2<-basis.x
#    X<-switch(basis.x$type,
#    "pc"= basis.x$x,
#    "raw"=basis.x$basis$data,
#    "bspline"=t(Data2fd(argvals = argvals(fdataobj), y = t(fdataobj$data),basisobj = basis.x)$coef)
#    )
#}
   if (!is.null(basis.x))  basis.x<-list("X"=basis.x)
   if (!is.null(basis.b))  basis.b<-list("X"=basis.b)     
   dataf<-data.frame("y"=group)
   ldata<-list("df"=dataf,"X"=fdataobj)
#   dataf<-list("y"=group,"X"=X)
   formula<-formula(y~X)     
#   basis.x<-list("X"=basis.x2)
   fit<-classif.tree(formula,data=ldata,basis.x=basis.x,basis.b=basis.b,...)       
  output<-list(fit=fit,formula=formula,fdataobj=dataf,basis.x=basis.x,
  basis.b=basis.b,group=group,C=C,max.prob=fit$max.prob,"prob.group"=fit$prob.group,
  "prob.classification"=fit$prob.classification,"group.est"=fit$group.est)
class(output)="classif"
return(output)
}



################################################################################
classif.gkam2boost=function(group,fdataobj,family = binomial(),weights=rep(1,n),
par.metric = NULL,par.np=NULL, offset=NULL,
control = list(maxit = 100,epsilon = 0.001, trace = FALSE,inverse="solve"),...)  {
formula<-formula(y~X)
C<-match.call()
a<-list()
mf <- match.call(expand.dots = FALSE)
m <- match(c("group","fdataobj","family","weights","par.metric","par.np","offset",
"control"), names(mf),0L)
   dataf<-data.frame("y"=group)
   ldata<-list("df"=dataf,"X"=fdataobj)
   newy<-y<-ldata$df$y
   if (!is.factor(y)) y<-as.factor(y)
   n<-length(y)
   newdata<-ldata
   ngroup<-nlevels(y)
   prob<-rep(NA,ngroup)
   ny<-levels(y)
#newdata<-data
if (!is.null(par.np)) {
   par.np=list("X"=par.np)
   if (is.null(par.np[["X"]]$Ker)) par.np[["X"]]$Ker=AKer.norm
   if (is.null(par.np[["X"]]$type.S)) par.np[["X"]]$type.S="S.NW"
   }
else          par.np =list("X"=list(Ker=AKer.norm,type.S="S.NW"))

if (ngroup==2) {
      newy<-ifelse(y==ny[1],0,1)
      newdata$df$y<-newy
      a[[1]]<-fregre.gkam(formula,data=newdata,family=family,weights=weights,
      par.metric=par.metric,par.np=par.np,offset=offset,control=control,...)
      yest<-ifelse(a[[1]]$fitted.values<.5,ny[1],ny[2])
            tab<-table(yest,y)
      prob[1]=tab[1,1]/sum(tab[,1])
      dtab<-dim(tab)
      if (dtab[1]==dtab[2])    {
           prob[2]=tab[2,2]/sum(tab[,2])
           names(prob)<-ny
         }
      else prob[2]<-0
      prob.group<-a$fitted.values
      yest<-factor(yest,levels=ny)
   }
else {
#   ny<-levels(y)
   prob.group<-array(NA,dim=c(n,ngroup))
   colnames(prob.group)<-ny
   for (i in 1:ngroup) {
              newy<-ifelse(y==ny[i],0,1)
              newdata$df$y<-newy
              a[[i]]<-fregre.gkam(formula,data=newdata,family=family,weights=weights,
              par.metric=par.metric,par.np=par.np,offset=offset,control=control,...)
              prob.group[,i]<-a[[i]]$fitted.values
              }
   yest<-ny[apply(prob.group,1,which.min)]#no sera which.max
   yest<-factor(yest,levels=ny)
   tab<-table(yest,y)
   for (ii in 1:ngroup) {
       prob[ii]=tab[ii,ii]/sum(tab[,ii])
       }
  names(prob)<-ny
}
max.prob=sum(diag(tab))/sum(tab)
output<-list(fdataobj=fdataobj,group=y,group.est=yest,
prob.classification=prob,prob.group=prob.group,C=C,m=m,max.prob=max.prob,fit=a
,formula=formula,data=newdata)
class(output)="classif"
return(output)
}

#######################
#######################
classif.glm.cv=function(formula,data,w=NULL,family = binomial(),
basis.x=NULL,basis.b=NULL,CV=FALSE,...){
C<-match.call()
a<-list()
mf <- match.call(expand.dots = FALSE)
m <- match(c( "formula","data","w","family","basis.x","basis.b","CV"), names(mf), 0L)
 tf <- terms.formula(formula)
 terms <- attr(tf, "term.labels")
 nt <- length(terms)
 if (attr(tf, "response") > 0) {
        response <- as.character(attr(tf, "variables")[2])
        pf <- rf <- paste(response, "~", sep = "")
    } else pf <- rf <- "~"
newy<-y<-data$df[[response]]
if (!is.factor(y)) y<-as.factor(y)
n<-length(y)
ny<-levels(y)
newdata<-data
prob2<-prob<-ngroup<-length(table(y))
if (ngroup==2) {
      ny<-as.numeric(names(table(y)))
      newy<-ifelse(y==ny[1],0,1)
      newdata$df$y<-newy
      a[[1]]<-fregre.glm(formula,data,family=family,basis.x=basis.x,basis.b=basis.b,CV=CV,...)

      if (CV) prediction<-a[[1]]$pred.cv
      else prediction<-a[[1]]$fitted.values
      yest<-ifelse(prediction<.5,0,1)
      yest<-factor(yest,levels=ny)
      tab<-table(yest,y)
   prob[1]=tab[1,1]/sum(tab[,1])
      dtab<-dim(tab)
      if (dtab[1]==dtab[2])    {
           prob[2]=tab[2,2]/sum(tab[,2])
           names(prob)<-ny
         }
      else prob[2]<-0
      prob.group<-prediction

      #devolver a mayores y estimada!
   }
else {
   ny<-as.numeric(names(table(y)))
   prob.group<-array(NA,dim=c(n,ngroup))
   colnames(prob.group)<-ny
   for (i in 1:ngroup) {
              newy<-ifelse(y==ny[i],0,1)
              newdata$df$y<-newy
              a[[i]]<-fregre.glm(formula,data,family=family,basis.x=basis.x,basis.b=basis.b,CV=CV,...)
      if (CV) prediction<-a[[1]]$pred.cv
      else prediction<-a[[1]]$fitted.values
            }
   yest<-ny[apply(prob.group,1,which.min)]
   yest<-factor(yest,levels=ny)
   tab<-table(yest,y)
   for (i in 1:ngroup) { prob[i]=tab[i,i]/sum(tab[,i])     }
           names(prob)<-ny
}
lev<-levels(y)
max.prob=sum(diag(tab))/sum(tab)
output<-list(formula=formula,data=data,group=y,group.est=factor(yest,levels=lev),
prob.classification=prob,prob.group=prob.group,C=C,fit=a,m=m,max.prob=max.prob)
class(output)="classif"
return(output)
}