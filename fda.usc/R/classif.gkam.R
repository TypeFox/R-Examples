#######################
#######################
classif.gkam=function(formula,family = binomial(),data,weights = rep(1, nobs),
 par.metric = NULL,par.np=NULL, offset=NULL,
 control = list(maxit = 100,epsilon = 0.001, trace = FALSE, inverse="solve"),...){
C<-match.call()
a<-list()
mf <- match.call(expand.dots = FALSE)
m <- match(c("formula","family","data","weigths","par.metric","par.np","offset","control"), names(mf),0L)
 tf <- terms.formula(formula)
 terms <- attr(tf, "term.labels")
 nt <- length(terms)
 if (attr(tf, "response") > 0) {
        response <- as.character(attr(tf, "variables")[2])
        pf <- rf <- paste(response, "~", sep = "")
    } else pf <- rf <- "~"
newy<-y<-data$df[[response]]
    nobs <- if (is.matrix(y))   nrow(y)
            else length(y)
if (!is.factor(y)) y<-as.factor(y)
n<-length(y)
newdata<-data
ny<-levels(y)
prob<-ngroup<-nlevels(y)
prob.group<-array(NA,dim=c(n,ngroup))
colnames(prob.group)<-ny
if (ngroup==2) {
      newy<-ifelse(y==ny[1],0,1)
      newdata$df[[response]]<-newy
      a[[1]]<-fregre.gkam(formula,family=family,data=newdata,weights=weights,par.metric=par.metric,par.np=par.np,offset=offset,control=control)
      yest<-ifelse(a[[1]]$fitted.values<.5,ny[1],ny[2])
      tab<-table(yest,y)
      prob[1]=tab[1,1]/sum(tab[,1])
      dtab<-dim(tab)
      if (dtab[1]==dtab[2])    {
           prob[2]=tab[2,2]/sum(tab[,2])
           names(prob)<-ny
         }
      else prob[2]<-0
      prob.group[,1]<-1-a[[1]]$fitted.values
      prob.group[,2]<-a[[1]]$fitted.values   
   }
else {
   for (i in 1:ngroup) {
               newy<-ifelse(y==ny[i],0,1)
              newdata$df[[response]]<-newy
              a[[i]]<-fregre.gkam(formula,family=family,data=newdata,weigths=weights,
              par.metric=par.metric,par.np=par.np,offset=offset,control=control,...)
              prob.group[,i]<-1-a[[i]]$fitted.values
            }
   yest<-ny[apply(prob.group,1,which.max)]
   yest<-factor(yest,levels=ny)
   tab<-table(yest,y)
   prob.group<-prob.group/apply(prob.group,1,sum)   
    for (i in 1:ngroup) {     prob[i]=tab[i,i]/sum(tab[,i])     }
               names(prob)<-ny
}
max.prob=sum(diag(tab))/sum(tab)
output<-list(formula=formula,data=data,group=y,group.est=yest,
prob.classification=prob,prob.group=prob.group,C=C,m=m,max.prob=max.prob,fit=a)
class(output)="classif"
return(output)
}
