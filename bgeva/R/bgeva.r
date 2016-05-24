bgeva <- function(formula.eq, data=list(), tau=-0.25, Hes=TRUE, gIM="a",  
                             iterlimSP=50, pr.tol=1e-6, 
                             gamma=1, aut.sp=TRUE, fp=FALSE, start.v=NULL, start.vo=1, rinit=1, rmax=100, 
                             fterm=sqrt(.Machine$double.eps), mterm=sqrt(.Machine$double.eps), 
        		     control=list(maxit=50,tol=1e-6,step.half=25,rank.tol=sqrt(.Machine$double.eps))){

                          
  conv.sp <- iter.sp <- iter.inner <- bs.mgfit <- wor.c <- qu.mag <- NULL
             

  # logit, probit, cauchit, log and cloglog  
                                                
  gam.fit <- gam(formula.eq, binomial(link="logit"), gamma=gamma, data=data)
  
  X     <- model.matrix(gam.fit); X.d2 <- dim(X)[2]
  l.sp <- length(gam.fit$smooth) 
  y <- gam.fit$y; sp <- c(gam.fit$sp)

  if(is.null(start.v)){ 

    if(start.vo==1) start.v <- c(coef(gam.fit))
    if(start.vo==2) start.v <- c(log(-log(mean(y))),rep(0,X.d2-1))
    if(start.vo==3) start.v <- c(log(-log(mean(y))),coef(gam.fit)[-1])  

  }

  names(start.v) <- names(coef(gam.fit))


  if(l.sp!=0 && fp==FALSE) qu.mag <- S.m(gam.fit)

  if(gIM=="a") func.opt <- bgeva.gIMa else func.opt <- bgeva.gIMb 

             
    fit  <- trust(func.opt, start.v, rinit=rinit, rmax=rmax, y=y, X=X,
                  X.d2=X.d2, tau=tau,sp=sp, qu.mag=qu.mag, gam.fit=gam.fit, fp=fp, l.sp=l.sp, Hes=Hes, blather=TRUE,
                  fterm=fterm, mterm=mterm, iterlim = 1e+4) 

    iter.if <- fit$iterations 
     
    if(aut.sp==TRUE && l.sp!=0 && fp==FALSE){

          stoprule.SP <- 1; conv.sp <- TRUE; iter.inner <- iter.sp <- 0

	  while( stoprule.SP > pr.tol ){ 
       
             coefo <- fit$argument
             spo <- sp
  
		 wor.c <- try(working.comp(fit,X[fit$good,],X.d2,fit$n)); if(class(wor.c)=="try-error") break
             
                	bs.mgfit <- try(magic(y=wor.c$rW.Z,X=wor.c$rW.X,sp=sp,S=qu.mag$Ss,
                        	           off=qu.mag$off,rank=qu.mag$rank,
                                	   gcv=FALSE,gamma=gamma,control=control))
                	if(class(bs.mgfit)=="try-error") {conv.sp <- FALSE; break} 
                	sp <- bs.mgfit$sp; iter.sp <- iter.sp + 1 
                               
             o.ests <- c(fit$argument)
            
             fit <- try(trust(func.opt, fit$argument, rinit=rinit, rmax=rmax, y=y, X=X,  
                              X.d2=X.d2, tau=tau, sp=sp, qu.mag=qu.mag, gam.fit=gam.fit, fp=fp, l.sp=l.sp, Hes=Hes, blather=TRUE,  
                              iterlim=1e+4, fterm=fterm, mterm=mterm),silent=TRUE)
             iter.inner <- iter.inner + fit$iterations

              if(class(fit)=="try-error"){ 
               fit  <- trust(func.opt, coefo, rinit=rinit, rmax=rmax, y=y, X=X, 
                             X.d2=X.d2, tau=tau, sp=spo, qu.mag=qu.mag, fp=fp, l.sp=l.sp, Hes=Hes, blather=TRUE, 
                             iterlim=1e+4, fterm=fterm, mterm=mterm)
               conv.sp <- FALSE; break
                                          } 
              
             n.ests <- c(fit$argument)

             if(iter.sp>iterlimSP){
              conv.sp <- FALSE; break
                                }
             stoprule.SP <- max(abs(o.ests-n.ests))    
           
          }
      
    }
    
    n <- fit$n 
    He <- fit$hessian
    logL <- fit$l
    He.eig <- eigen(He)
    k.e <- sum(as.numeric(He.eig$val<sqrt(.Machine$double.eps)))
    
     if(k.e!=0){
      ind.e <- (length(He.eig$val)-(k.e-1)):length(He.eig$val)
      min.e <- min(He.eig$val[1:(ind.e[1]-1)])
      for(i in 1:k.e) He.eig$val[ind.e[i]] <- min.e/10^i  
      Vb <- He.eig$vec%*%tcrossprod(diag(1/He.eig$val),He.eig$vec)      
     }else{
      Vb <- He.eig$vec%*%tcrossprod(diag(1/He.eig$val),He.eig$vec) 
    }
        
                                      
if (l.sp!= 0 && fp == FALSE) {
        HeSh <- He - fit$S.h
        F <- Vb %*% HeSh
    }
    else{HeSh <- He; F <- diag(rep(1,dim(Vb)[1]))}
 
  t.edf <- sum(diag(F))

L <- list(fit=fit, coefficients=fit$argument, gam.fit=gam.fit, sp=sp, fp=fp,
          iter.if=iter.if, conv.sp=conv.sp, iter.sp=iter.sp, iter.inner=iter.inner,  
          tau=tau, n=n, X=X, Xr=fit$Xr, good=fit$good, X.d2=X.d2, logL = logL, 
          l.sp=l.sp, He=He, HeSh=HeSh, Vb=Vb, F=F, 
          t.edf=t.edf, bs.mgfit=bs.mgfit, wor.c=wor.c, eta=fit$eta)

class(L) <- "bgeva"

L

}
