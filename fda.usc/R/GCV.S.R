GCV.S=function(y,S,criteria="GCV",W=NULL,trim=0,draw=FALSE,metric=metric.lp,...){
    isfdata<-is.fdata(y)
    tab=list("GCV","AIC","FPE","Shibata","Rice")
    type.i=pmatch(criteria,tab)
    n=ncol(S);l=1:n
    if (isfdata)  {
         nn<-nrow(y)
         if (is.null(W)) W<-diag(nn)
         y2=t(y$data)
         y.est=t(S%*%y2)
         y.est<-fdata(y.est,y$argvals, y$rangeval, y$names)
         e <- y - y.est
#         e$data<-sqrt(W)%*%(e$data)   
         ee <- drop(norm.fdata(e,metric=metric,...)[,1]^2)
         if (trim>0) {
            e.trunc=quantile(ee,probs=(1-trim),na.rm=TRUE,type=4)
            ind<-ee<=e.trunc
            if (draw)  plot(y,col=(2-ind))
            l<-which(abs(ee)<=e.trunc)
            res = mean(ee[ind],na.rm=TRUE)
            }
        else  res = mean(ee, na.rm = TRUE)
         }
   else {
        if (is.null(W)) W<-diag(n)
        if (is.matrix(y)&&(ncol(y)==1) ){y2<-y;draw<-FALSE}
        else if (is.vector(y)){y2<-y;draw<-FALSE}
        else stop("y is not a fdata,  vector or matrix")
     y.est=S%*%y2
     e=y2-y.est
     if (trim>0) {
             ee = t(e)
             e.trunc=quantile(abs(ee),probs=(1-trim),na.rm=TRUE,type=4)
             l<-which(abs(ee)<=e.trunc)
             res=traza(t(e[l])%*%W[l,l]%*%e[l])
             }
     else    res=traza(t(e)%*%W%*%e)
    }
    d<-diag(S)[l] 
    df<-sum(d)
    if (is.na(type.i))   {
                   if (mean(d,na.rm=TRUE)>0.5) vv=Inf
                   else   vv= 1/(1-2*mean(d,na.rm=TRUE))        }
    else {
        vv<-switch(type.i,
                   "1"=if (type.i==1)  vv=(1-mean(d,na.rm=TRUE))^(-2),
                   "2"=if (type.i==2)  vv= exp(2*mean(d)),
                   "3"=if (type.i==3)  vv=(1+mean(d,na.rm=TRUE))/(1-mean(d,na.rm=TRUE)),
                   "4"=if (type.i==4)  vv=(1+2*mean(d,na.rm=TRUE)),
                   "5"=if (type.i==5)  {
                           if (mean(d,na.rm=TRUE)>0.5) vv=Inf
                           else   vv= 1/(1-2*mean(d,na.rm=TRUE))   }
                  )
         }
out<-res*vv/n
 attr(out, "df") <- df           
return(out) }
  
