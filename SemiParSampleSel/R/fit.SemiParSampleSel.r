fit.SemiParSampleSel <- function(funcop, start.v, dat, VC, qu.mag, sp, iterlimsp, pr.tolsp, rinit, rmax, univariate) {


  parsc <- rep(VC$parscale, length(start.v) )  

  fit <- trust(funcop, start.v, dat=dat, VC=VC, qu.mag=qu.mag, sp=sp, rinit=rinit, rmax=rmax, blather=TRUE, parscale=parsc, iterlim=1e+4)

  iter.if <- fit$iterations  

  conv.sp <- iter.sp <- iter.inner <- bs.mgfit <- wor.c <- magpp <- NULL
  
  if(univariate == TRUE) l.sp1 <- 0 else l.sp1 <- VC$l.sp1 
  
  
  if(VC$fp==FALSE && (l.sp1!=0 || VC$l.sp2!=0 || VC$l.sp3!=0 || VC$l.sp4!=0 || VC$l.sp5!=0)){
                    
      stoprule.SP <- 1; conv.sp <- TRUE; iter.inner <- iter.sp <- 0                                  
                                      

      while( stoprule.SP > pr.tolsp ){   
      
                 spo <- sp 

		 wor.c <- working.comp(fit) 
             
                 bs.mgfit <- try(magic(y=wor.c$Z, X=wor.c$X, sp=sp, S=qu.mag$Ss, off=qu.mag$off,
                                       rank=qu.mag$rank, gcv=FALSE, gamma=VC$infl.fac))
                 if(class(bs.mgfit)=="try-error") { conv.sp <- FALSE; break } 
                 
                 sp <- bs.mgfit$sp; iter.sp <- iter.sp + 1; names(sp) <- names(spo) 
                 
                 o.ests <- c(fit$argument)  

                 fit  <- trust(funcop, o.ests, dat=dat, VC=VC, qu.mag=qu.mag, parscale=parsc,
                                   sp=sp, rinit=rinit, rmax=rmax, blather=TRUE, iterlim=1e+4)

                 iter.inner <- iter.inner + fit$iterations  
                              
             if(iter.sp > iterlimsp){conv.sp <- FALSE; break }

             stoprule.SP <- max(abs(o.ests - c(fit$argument))) 
           
        }

magpp <- magic.post.proc(wor.c$X, bs.mgfit)

  }else{
    
    wor.c <- working.comp(fit) 
    
    bs.mgfit <- magic(wor.c$Z, wor.c$X, numeric(0), list(), numeric(0))    
    magpp    <- magic.post.proc(wor.c$X, bs.mgfit)
    
    
    
    }

                  list(fit = fit, 
                       iter.if = iter.if, 
                       conv.sp = conv.sp, 
                       iter.sp = iter.sp, 
                       iter.inner = iter.inner, 
                       bs.mgfit = bs.mgfit, 
                       wor.c = wor.c, 
                       sp = sp, magpp = magpp)

}




