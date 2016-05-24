lmp<-function(formula,data,p) UseMethod("lmp")

lmp.default<-function(formula,data=list(),p=NULL){
 ll<-lm(formula=formula, data=data)
 r <- residuals(ll)
 s <- sqrt(deviance(ll)/df.residual(ll))
 hii <- lm.influence(ll)$hat
 rs <- r/(s * sqrt(1 - hii))
 cl<-match.call()
 cl[[1]]<-as.name("lmp")
 b<-vector()
 bb<-vector()
 nvar<-ll$rank
 for (i in 1:nvar) {bb[i]<-ll$coef[[i]]}
 fit<-vector()
 N<-nrow(ll$model)
 y<-ll$model[[1]]
 M<-as.matrix(ll$model)
 M[,1]<-1     
 res<-ll$residual
  f<-function(b){
  sum(abs(y-M%*%b)^pp)
  }
 if (is.null(p)){
  knp<-FALSE
  pp<-estimatep(res,mean(res),2)
  op<-optim(bb,f,method="BFGS")
  for (i in 1:nvar){b[i]<-op$par[i]}
  fit<-c(M%*%b) 
  res<-c(y-fit)
 p<-estimatep(res,mean(res),2)
 i1<-0
 iter<-0
 while (abs(pp-p)>0.0001 || abs(bb-b)>0.0001){
 pp<-p
 bb<-b
 op<-optim(bb,f,method="BFGS")
 for (i in 1:nvar){b[i]<-op$par[i]}
 fit<-c(M%*%b) 
 res<-c(y-fit)
 p<-estimatep(res,mean(res),2)
 i1<-i1+1
 if (i1==100) {iter<-1; break}
        }
        }
 else {
 knp<-TRUE
 pp<-p
 op<-optim(bb,f)
 for (i in 1:nvar){b[i]<-op$par[i]}
 fit<-c(M%*%b) 
 res<-c(y-fit)
 }
 for (i in 1:nvar){ll$coefficients[[i]]<-b[i]}
 names(res)<-as.character(1:N)
 ll$residuals<-res
 ll$call<-cl
 ll$knp<-knp
 ll$p<-p 
 names(fit)<-as.character(1:N)
 ll$fitted.values<-fit 
 if(knp==FALSE) ll$iter<-iter
 ll$rs<-rs
 class(ll)<-c("lmp","lm")
 ll
}

