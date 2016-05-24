SemiParBIVProbit.fit  <- function(func.opt, start.v, 
                                   rinit, rmax, iterlim, iterlimsp, tolsp,
                                   respvec, VC,
                                   sp = NULL, qu.mag = NULL){ 

  parsc <- rep(VC$parscale, length(start.v) )  
  
  sc <- TRUE

  fit  <- try( trust(func.opt, start.v, rinit = rinit, rmax = rmax, parscale = parsc,
                     respvec = respvec, VC = VC, sp = sp, qu.mag = qu.mag, blather = TRUE, 
                     iterlim = iterlim), silent = sc)   

  if(class(fit) == "try-error" || is.null(fit$l)){
  
    fit  <- try( trust(func.opt, start.v, rinit = rinit, rmax = rmax, parscale = parsc,
                       respvec = respvec, VC = VC, sp = sp, qu.mag = qu.mag, blather = TRUE, 
                       iterlim = iterlim/4), silent = sc)  
                     
        if(class(fit) == "try-error"|| is.null(fit$l)){
        
            fit  <- try( trust(func.opt, start.v, rinit = rinit, rmax = rmax, parscale = parsc,
                               respvec = respvec, VC = VC, sp = sp, qu.mag = qu.mag, blather = TRUE, 
                               iterlim = iterlim/10), silent = sc)   

if((class(fit) == "try-error" || is.null(fit$l)) && VC$Cont == "NO"  && VC$gamlssfit == FALSE && VC$ccss == "no"  && VC$triv == FALSE) stop("It is not possible to fit the model. Try re-fitting the model and setting gamlssfit = TRUE.\n Also, read the WARNINGS section in ?SemiParBIVProbit.")
if((class(fit) == "try-error" || is.null(fit$l)) && VC$Cont == "NO"  && VC$gamlssfit == TRUE  && VC$ccss == "no"  && VC$triv == FALSE) stop("It is not possible to fit the model. Read the WARNINGS section in ?SemiParBIVProbit.")         
if((class(fit) == "try-error" || is.null(fit$l)) && VC$Cont == "YES" && VC$gamlssfit == FALSE && VC$ccss == "no"  && VC$triv == FALSE) stop("It is not possible to fit the model. Try re-fitting the model and setting gamlssfit = TRUE.\n Also, read the WARNINGS section in ?copulaReg.")
if((class(fit) == "try-error" || is.null(fit$l)) && VC$Cont == "YES" && VC$gamlssfit == TRUE  && VC$ccss == "no"  && VC$triv == FALSE) stop("It is not possible to fit the model. Read the WARNINGS section in ?copulaReg.")
if((class(fit) == "try-error" || is.null(fit$l)) && VC$Cont == "NO"  && VC$gamlssfit == FALSE && VC$ccss == "yes" && VC$triv == FALSE) stop("It is not possible to fit the model. Try re-fitting the model and setting gamlssfit = TRUE.\n Also, read the WARNINGS section in ?copulaSampleSel.")
if((class(fit) == "try-error" || is.null(fit$l)) && VC$Cont == "NO"  && VC$gamlssfit == TRUE  && VC$ccss == "yes" && VC$triv == FALSE) stop("It is not possible to fit the model. Read the WARNINGS section in ?copulaSampleSel.")             
if((class(fit) == "try-error" || is.null(fit$l)) && VC$triv == TRUE)                                                                   stop("It is not possible to fit the model. Read the WARNINGS section in ?SemiParTRIVProbit.")             
             
                                     }
                     
  }


  iter.if <- fit$iterations  

  conv.sp <- iter.sp <- iter.inner <- bs.mgfit <- wor.c <- magpp <- NULL
  
  #####################################################################
  
  l.sp1 <- VC$l.sp1 
  l.sp2 <- VC$l.sp2 
  l.sp3 <- VC$l.sp3 
  l.sp4 <- VC$l.sp4 
  l.sp5 <- VC$l.sp5 
  l.sp6 <- VC$l.sp6 
  l.sp7 <- VC$l.sp7 
  
  if(VC$Cont == "NO" && VC$margins[2] %in% VC$m2 && respvec$univ == 1  && VC$triv == FALSE) l.sp1 <- l.sp4 <- l.sp5 <- l.sp6 <- l.sp7 <- 0
  if(VC$Cont == "NO" && VC$margins[2] %in% VC$m3 && respvec$univ == 1  && VC$triv == FALSE) l.sp1 <- l.sp5 <- l.sp6 <- l.sp7 <- 0  
  
  if(VC$Cont == "YES" && VC$margins[1] %in% VC$m2 && respvec$univ == 2 && VC$triv == FALSE) l.sp2 <- l.sp4 <- l.sp5 <- l.sp6 <- l.sp7 <- 0
  if(VC$Cont == "YES" && VC$margins[1] %in% VC$m3 && respvec$univ == 2 && VC$triv == FALSE) l.sp2 <- l.sp4 <- l.sp6 <- l.sp7 <- 0

  if(VC$Cont == "YES" && VC$margins[2] %in% VC$m2 && respvec$univ == 3 && VC$triv == FALSE) l.sp1 <- l.sp3 <- l.sp5 <- l.sp6 <- l.sp7 <- 0
  if(VC$Cont == "YES" && VC$margins[2] %in% VC$m3 && respvec$univ == 3 && VC$triv == FALSE) l.sp1 <- l.sp3 <- l.sp5 <- l.sp7 <- 0

  #####################################################################
  
  
    if((VC$fp==FALSE && (l.sp1!=0 || l.sp2!=0 || l.sp3!=0 || l.sp4!=0 || l.sp5!=0 || l.sp6!=0 || l.sp7!=0)) ){

       stoprule.SP <- 1; conv.sp <- TRUE; iter.inner <- iter.sp <- 0  

	  while( stoprule.SP > tolsp ){ 

             fito   <- fit$l
             o.ests <- c(fit$argument) 
             spo    <- sp 
             
             wor.c <- working.comp(fit) 
             
              qu.magM <- qu.mag 
              if(VC$triv == TRUE) { if(!(VC$penCor %in% c("unpen"))) qu.magM <- fit$qu.mag.a} 
  
                	bs.mgfit <- try(magic(y = wor.c$Z,  
                	                      X = wor.c$X,
                	                      sp= sp, S = qu.magM$Ss,
                        	              off = qu.mag$off, rank = qu.mag$rank,
                                	      gcv = FALSE,
                                	      gamma = VC$infl.fac), silent = sc)
                		if(class(bs.mgfit)=="try-error") {conv.sp <- FALSE; break} 
                		
                	sp <- bs.mgfit$sp; iter.sp <- iter.sp + 1; names(sp) <- names(spo) 
                        

             fit <- try( trust(func.opt, o.ests, rinit=rinit, rmax = rmax,  parscale = parsc,  
                          respvec = respvec, VC = VC, 
                          sp = sp, qu.mag = qu.mag, 
                          blather = TRUE, iterlim = iterlim), silent = sc) 
                          
                          
                          
                          if(class(fit) == "try-error" || is.null(fit$l)){conv.sp <- FALSE
                          
                                                         fit <- try( trust(func.opt, o.ests, rinit=rinit, rmax = rmax,  parscale = parsc,  
			                                              respvec = respvec, VC = VC, 
			                                              sp = spo, qu.mag = qu.mag, 
                                                                      blather = TRUE, iterlim = iterlim), silent = sc)  
                                                                      
if((class(fit) == "try-error" || is.null(fit$l)) && VC$Cont == "NO"  && VC$gamlssfit == FALSE && VC$ccss == "no"  && VC$triv == FALSE) stop("It is not possible to fit the model. Try re-fitting the model and setting gamlssfit = TRUE.\n Also, read the WARNINGS section in ?SemiParBIVProbit.")
if((class(fit) == "try-error" || is.null(fit$l)) && VC$Cont == "NO"  && VC$gamlssfit == TRUE  && VC$ccss == "no"  && VC$triv == FALSE) stop("It is not possible to fit the model. Read the WARNINGS section in ?SemiParBIVProbit.")         
if((class(fit) == "try-error" || is.null(fit$l)) && VC$Cont == "YES" && VC$gamlssfit == FALSE && VC$ccss == "no"  && VC$triv == FALSE) stop("It is not possible to fit the model. Try re-fitting the model and setting gamlssfit = TRUE.\n Also, read the WARNINGS section in ?copulaReg.")
if((class(fit) == "try-error" || is.null(fit$l)) && VC$Cont == "YES" && VC$gamlssfit == TRUE  && VC$ccss == "no"  && VC$triv == FALSE) stop("It is not possible to fit the model. Read the WARNINGS section in ?copulaReg.")
if((class(fit) == "try-error" || is.null(fit$l)) && VC$Cont == "NO"  && VC$gamlssfit == FALSE && VC$ccss == "yes" && VC$triv == FALSE) stop("It is not possible to fit the model. Try re-fitting the model and setting gamlssfit = TRUE.\n Also, read the WARNINGS section in ?copulaSampleSel.")
if((class(fit) == "try-error" || is.null(fit$l)) && VC$Cont == "NO"  && VC$gamlssfit == TRUE  && VC$ccss == "yes" && VC$triv == FALSE) stop("It is not possible to fit the model. Read the WARNINGS section in ?copulaSampleSel.")             
if((class(fit) == "try-error" || is.null(fit$l)) && VC$triv == TRUE)                                                                   stop("It is not possible to fit the model. Read the WARNINGS section in ?SemiParTRIVProbit.")             
     
    
                                                        # break
                          } 
                          
                                                             
             iter.inner <- iter.inner + fit$iterations   	                                    
              
             if(iter.sp >= iterlimsp){conv.sp <- FALSE; break }

             stoprule.SP <- abs(fit$l - fito)/(0.1 + abs(fit$l))  # max(abs(o.ests - c(fit$argument)))
                  
           
          } # end smoothing fitting loop
          

       if(VC$gc.l == TRUE) gc()   
       
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




