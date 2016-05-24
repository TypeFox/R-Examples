vbfr<-function(age=NULL,L=NULL,agephi=NULL,agepsi=NULL,
               graph=TRUE,gestimate=TRUE,
               Lphiparms=c(NA,NA,NA),Lchiparms=c(NA,NA,NA),
               Lpsiparms=c(NA,NA,NA),control=list(maxiter=10000)){
     if(is.null(age)) stop ("age vector does not exist")
     if(is.null(L)) stop ("L vector does not exist")
     if(is.null(agephi)) stop ("agephi does not exist")
     if(is.null(agepsi)) stop ("agepsi does not exist")
     if(agephi>agepsi) stop("agephi should be smaller than agepsi")

#get data 
     x<-as.data.frame(cbind(age,L)) 
     x<-x[!is.na(x$age) & !is.na(x$L),]  

   if(gestimate==TRUE){
       Lphi.init <-  with(subset(x, round(age,0)==agephi), mean(L))
       Lchi.init <-  with(subset(x, round(age,0)==round((agephi+agepsi)/2,0)), mean(L))
       Lpsi.init <-  with(subset(x, round(age,0)==agepsi), mean(L))
	  Lphiparms<-c(Lphi.init,0.5 * Lphi.init,1.5*Lphi.init)
	  Lchiparms<-c(Lchi.init,0.5 * Lchi.init,1.5*Lchi.init)
	  Lpsiparms<-c(Lpsi.init,0.5 * Lpsi.init,1.5*Lpsi.init)
       }
   if(gestimate==FALSE){
      if(any(is.na(Lphiparms))) stop ("Missing values in Lphiparms")
      if(any(is.na(Lchiparms))) stop ("Missing values in Lchiparms")
      if(any(is.na(Lpsiparms))) stop ("Missing values in Lpsiparms")
      if(length(Lphiparms)>3|length(Lchiparms)>3|length(Lpsiparms)>3) stop("There should be only 3 values in Lphiparms,Lchiparms or Lpsiparms")
     if(length(Lphiparms)<3|length(Lchiparms)<3|length(Lpsiparms)<3) stop("There should be 3 values in Lphiparms,Lchiparms or Lpsiparms") 
      if(Lphiparms[2]>Lphiparms[1]|Lphiparms[3]<Lphiparms[1]) stop("Lower or upper value < or > start value of Lphiparms")
      if(Lchiparms[2]>Lchiparms[1]|Lchiparms[3]<Lchiparms[1]) stop("Lower or upper value < or > start value of Lchiparms")
      if(Lpsiparms[2]>Lpsiparms[1]|Lpsiparms[3]<Lpsiparms[1]) stop("Lower or upper value < or > start value of Lpsiparms")
   }
parms<-c(Lphiparms[1],Lchiparms[1],Lpsiparms[1])
 if(gestimate==FALSE){
 	mod1<-try(nls(L~lphi+(lpsi-lphi)*(1-((lpsi-lchi)/(lchi-lphi))^(2*(x$age-agephi)/(agepsi-agephi)))/
           (1-((lpsi-lchi)/(lchi-lphi))^2), data=x,start=list(lphi=parms[1],lchi=parms[2],lpsi=parms[3]),
           control=control,algorithm="port",lower=c(Lphiparms[2],Lchiparms[2],Lpsiparms[2]),
           upper=c(Lphiparms[3],Lchiparms[3],Lpsiparms[3])))
 }
if(gestimate==TRUE){
 	mod1<-try(nls(L~lphi+(lpsi-lphi)*(1-((lpsi-lchi)/(lchi-lphi))^(2*(x$age-agephi)/(agepsi-agephi)))/
           (1-((lpsi-lchi)/(lchi-lphi))^2), data=x,start=list(lphi=parms[1],lchi=parms[2],lpsi=parms[3]),
           control=control))
 }
 if(class(mod1)!="try-error"){
     if(graph==TRUE){
       par(mfrow=c(1,2))
       plot(x$L~x$age,xlab="Age",ylab="Length")
       preds<-as.data.frame(unique(cbind(predict(mod1),x$age)))
       preds<-preds[order(preds[,2]),]
       lines(preds[,1]~preds[,2],col="red")
       plot(summary(mod1)$residuals~x$age,xlab="Age",ylab="Residuals")
       abline(h=0,col="red")
     }
      return(mod1)
  }
  if(class(mod1)=="try-error"){
    return("Fit Failed.")
   }
}#end function
