#######################################################################
##
## Function: anchors.chopit()
## Author  : Jonathan Wand <wand@stanford.edu>
## Created : 2007-08-20
##
## Extracted and refined from former chopit()  
##
#######################################################################
anchors.chopit.fit <- function(data, parm, count, options) {
  
  fitted    <- options$fitted
  verbose   <- options$verbose
  debug     <- options$debug
  print.oob <- options$print.oob
  offset    <- 0
  do.gr     <- options$analytic
  get.LL <- FALSE  
  
  ###################################################################
  ##
  ## A.2 MODEL/LIKELIHOOD
  ##
  ###################################################################
  if (options$debug > 0) {
    cat("CHECKPOINT: specify model log-likelihood\n")
  }
  tidx <- function(i.self,i.cat,tmp.n.cat=NULL) {
    if (!is.null(tmp.n.cat))
      n.cat <- tmp.n.cat
    return( (n.cat-1)*(i.self-1)+i.cat )
  }
  
  chopit.llog <- function(b,do.gr=FALSE,verbose=FALSE) {
    if (debug > 1) cat("chopit.llog",do.gr,verbose,"\n")
    
    ## JW: what if n.self == 0?
    ## JW: I have sub'ed in count$n.tau.set for n.self in gamma = ... 
    if (!do.gr) {
    } else {
      GR <- list(gamma1  = rep(0, count$nvars.gamma1),
                 gamma   = rep(0, count$nvars.gamma * (count$n.cat-2)*count$n.tau.set ),
                 se.re   = NULL,
                 se.self = NULL,
                 se.vign = NULL,
                 theta   = NULL,
                 beta    = NULL)

      if (options$vign.var == "homo") 
        GR$se.vign <- 0
      
    }
    
    #### A.2.a EXTRACTION of parameters
    parm$pvec <- b
    parm <- unpackv(parm)
    #print(parm$start)
          
    ## SET-UP...
    LL.vign <- LL.self <- NULL
    vign.prob <- list()
    self.prob <- list()

    ## linear vs. exponentiated parameterization
    if (!options$linear) {
      sigma.vign <- exp(parm$start$sigma.vign)
      sigma.self <- exp(parm$start$sigma.self)
      sigma.re   <- exp(parm$start$sigma.re)
    } else {
      sigma.vign <- parm$start$sigma.vign
      sigma.self <- parm$start$sigma.self
      sigma.re   <- parm$start$sigma.re
    }
    beta <- parm$start$beta
    gamma<- parm$start$gamma
    gamma1<- parm$start$gamma1
    for (iij in 1:count$n.vign.set) {
      zti <- ifelse( count$n.vign.set > 1, iij, "")
      txt1 <- paste("theta",zti," <- parm$start$theta",zti,sep="")
      eval(parse(text=txt1))
    }

     if (options$debug > 3  ) {
       cat("TEST: beta value in chopit.llog():\n")
       print("b")        ;print(b)
       print("beta")     ;print(beta)
       print("gamma")    ;print(gamma)
       print("gamma1")   ;print(gamma1)
#       print("theta")   ;print(theta)
       print("sigma.vign");print(sigma.vign)
       print("sigma.self");print(sigma.self)
       print("sigma.re")  ;print(sigma.re)  
     }
    
    ## 
    ## SANITY CHECKS!
    ## 
    if (any(sigma.vign <= 0)) {
      if (print.oob) {
        cat("\nILLEGAL vignette variance value(s):",sigma.vign,"\n")
        cat("Possible reasons for this error:\n")
        cat("1. you have given invalid (zero or negative) starting values\n")
        cat("2. your model is severely misspecified\n")
        cat("3. the optimizer is temporarily exploring (and rejecting) an illegal part of the parameter space\n\n")
      }
      if (options$optimizer!="genoud") return(NA) else return(1e16)
    }
    if (any(sigma.self <= 0)) {
      if (print.oob) {
        cat("\nILLEGAL self-question variance value(s):",sigma.self,"\n")
        cat("Possible reasons for this error:\n")
        cat("1. you have given invalid (zero or negative) starting values\n")
        cat("2. your model is severely misspecified\n")
        cat("3. the optimizer is temporarily exploring (and rejecting) an illegal part of the parameter space\n\n")
      }
      if (options$optimizer!="genoud") return(NA) else return(1e16)
    }
    if (any(sigma.re <= 0)) {
      if (print.oob) {
        cat("\nILLEGAL random effect variance value(s):",sigma.re,"\n")
        cat("Possible reasons for this error:\n")
        cat("1. you have given invalid (zero or negative) starting values\n")
        cat("2. your model is severely misspecified \n")
        cat("3. the optimizer is temporarily exploring (and rejecting) an illegal part of the parameter space\n")
        cat("4. you do not need to estimate a random effect\n\n")
      }
      if (options$optimizer!="genoud") return(NA) else return(1e16)
    }
    
    ##
    ## End sanity check
    ## 

    if (options$debug > 1) 
      cat("anchors.chopit.fit::chopit.llog : end sanity checks\n")

    ## JW: GAMMA and TAU --- replaced with function...
#    print("going into calc...\n")
    Vtaus <- anchors.chopit.tau(data$v0v,
                                gamma,
                                data$v0v1,
                                gamma1,
                                offset,
                                count$nobs.vign,
                                count$n.cat,
                                count$n.tau.set,
                                count$nvars.gamma,
                                count$tau.start.idx,
                                options$linear,
                                options$debug,
                                do.gr,
                                options$verbose)
#     cat("Vtaus\n"); print(Vtaus)

    ## Error trap
    if (is.null(Vtaus) || !all(is.finite(Vtaus)) ) {
      if (options$debug > 1) {
        cat("chopit.llog: infinite value trap 1\n")
        if (do.gr) {
          cat("N deriv",sum(unlist(parm$estimated)) ,"\n")
          print(Vtaus)
          print(b)
        }
        
      }
      if (!do.gr) {
        if (options$optimizer!="genoud")  return(1e16)
        return(Inf)
      } else {
        if (options$optimizer!="genoud")  return(rep(1e16,sum(unlist(parm$estimated))))
        return( rep(Inf,sum(unlist(parm$estimated))) )
      }
    }
    
    if (count$n.self > 0) {
      Staus <- anchors.chopit.tau(data$v0s,gamma,data$v0s1,gamma1,offset,
                        count$nobs.self,count$n.cat,count$n.tau.set,
                        count$nvars.gamma, count$tau.start.idx,
                        options$linear,options$debug,
                                  do.gr,options$verbose)
#       cat("Staus\n"); print(Staus)
  
  
      ## Error trap
      if (is.null(Staus) || is.null(Vtaus)) {
        if (debug > 0) cat("chopit.llog: infinite value trap 2\n")
        if (options$optimizer!="genoud")  return(1e16)
        return(Inf)
      }
      ## More sanity checks
      if (!all(is.finite(Staus)) || is.null(Staus) ||
          !all(is.finite(Vtaus)) || is.null(Vtaus)
          ) {
        if (print.oob) {
          cat("\nIllegal tau values.\n")
        }
        if (debug > 0) cat("chopit.llog: infinite value trap 3\n")
        if (options$optimizer!="genoud") return(NA) else return(1e16)
      }
    }

    
    if (options$rprof)
      Rprof("chopit.llog",append=TRUE,interval=0.005)
    
    #### A.2.b VIGNETTE COMPONENT OF LIKELIHOOD
    ## NEW: sigma.vign (vec)          -- in parm 
    ##      thetaN (vec)          -- in parm
    ##      univ   (list of list) -- DONE
    ##      n.vign (vec)          -- DONE
    ##      vign.self (vec)       -- in opt
    llik1 <- 0
    ## Loop unnecessary with 1 set of vigns
    if (count$n.vign.set > 0) {
      for (idx.vign.set in 1:count$n.vign.set) {
        vign.prob[[idx.vign.set]] <- list()

        zti <- ifelse( count$n.vign.set > 1, idx.vign.set, "")
        
        tmp.theta <- eval(parse(text=paste("theta",zti,sep="")))
        ## which of the taus will be used for this particular set of vign?
        voffset <- ( parm$vign.map[idx.vign.set] -1 )*count$nobs.vign*(count$n.cat-1)
        tmp.taus <- Vtaus[ voffset + c(1:(count$nobs.vign*(count$n.cat-1)))]

#        cat("V OFFSET",voffset,"tmp.taus",
#            range(voffset + c(1:(count$nobs.vign*(count$n.cat-1)))),"\n")
        ## extract vignette responses
#        print(tmp.taus)
        
        ## IF NO vign.set LOOP, then just replace zz with "correct" single version
        #zz <- eval(parse(text=paste("z",idx.vign.set,sep="")))
        zz <- data$z0
        
        ## LOOP over vignettes in set
        for (idx.vign in 1:count$n.vign) {

          if (options$debug > 1)
            cat("Doing vignette",idx.vign,"in set",idx.vign.set,"\n")

          if (options$vign.cut == "hetero") {
            voffset <- ( idx.vign - 1 )*count$nobs.vign*(count$n.cat-1)
            tmp.taus <- Vtaus[ voffset + c(1:(count$nobs.vign*(count$n.cat-1)))]
          }
          
          if (options$vign.var=="homo") {
            tmp.vign.se <- sigma.vign
          } else if (idx.vign.set == 1) {
            tmp.vign.se <- sigma.vign[idx.vign]
          } else {
            ## JW: is this right index??
            tmp.vign.se <- sigma.vign[sum( count$n.vign[1:(idx.vign.set-1)] ) + idx.vign]
          }
          
          if (!do.gr) {
            ## JW: do a fitted == TRUE version!          
#             if (options$debug > 0) {
#               cat("ll.oprobit1\n") 
# #              print(zz[,idx.vign])
# #             print(tmp.taus)
# #              print(rep( tmp.theta[idx.vign], count$nobs.vign))
#               print(idx.vign)
#               print( tmp.taus[1:10])
#               print(count$nobs.vign)
#               print(tmp.theta[idx.vign])
#               print(tmp.vign.se)
#               print(n.cat)
#             }
            
            val <- ll.oprobit(zz[,idx.vign],
                              rep( tmp.theta[idx.vign], count$nobs.vign),
                              tmp.vign.se,
                              tmp.taus,
                              count$n.cat,
                              options$debug)

            llik1 <- llik1 + sum(val)

            ## keep track of -LL
            if (options$debug > 1)
              cat("VIGNETTE LLIK1:",
                  idx.vign.set,
                  idx.vign,
                  val,"\n")
            if (get.LL) {
              LL.vign <- c(LL.vign,-sum(val))
            }
            
          } else {
            gval <- gr.oprobit.vign.C(zz[,idx.vign],
                                    rep( tmp.theta[idx.vign], count$nobs.vign),
                                    tmp.vign.se,
                                    tmp.taus,
                                    data$v0v,
                                    data$v0v1,
                                    count$n.cat,
                                    options$debug)

            if (debug > 1) {
              cat("Gradients of vign\n")
              print(gval)
            }
            
            GR$theta  <- c(GR$theta, gval$theta)

            GR$gamma1 <- GR$gamma1+ gval$gamma1
            GR$gamma  <- GR$gamma + gval$gamma
            
#             ## JW??
#             ii <- (idx.vign.set-1) * (n.cat-2)*nvars.gamma + c(1:((n.cat-1)*nvars.gamma))
# 
# #            ## JW
# #            print("GAMMA")
# #            print(ii)
# #            print(gval$gamma)
# #            print(GR$gamma)
#             
#             ## JW??
#             GR$gamma[ii]  <- GR$gamma[ii] + gval$gamma
            
            if (options$vign.var=="homo") {
              GR$se.vign<- GR$se.vign + gval$sigma
            } else {
              GR$se.vign<- c(GR$se.vign, gval$sigma)
            }
          }
        }
      }
    }

#    ## JW
#    if (do.gr) {
#      print("VIGN")
#      print(GR)
#    }
  
#### A.2.c SELF-REPORT COMPONENT OF LIKELIHOOD WITH RANDOM EFFECT
    llik2 <- 0
    if (count$n.self > 0) {
      xb <- data$x0 %*% beta
      if (parm$estimated$sigma.re == FALSE || count$n.self == 1) { ## NO random effect
        
        for (idx.self in 1:count$n.self) { ## loop over self questions w/o RE

          ## extract
          voffset <- ( idx.self -1 )*count$nobs.self*(count$n.cat-1)
          tmp.taus <- Staus[ voffset + c(1:(count$nobs.self*(count$n.cat-1)))]
          tmp.sigma.self <- as.numeric(sigma.self)[1]  #[idx.self]

          #cat("S VOFFSET",voffset,"tmp.taus",range(voffset + c(1:(count$nobs.self*(count$n.cat-1)))),"\n")

          
#          ## and select
#          if (parm$estimated$sigma.re == FALSE) {
#            tmp.sigma <- tmp.sigma.self
#          } else if (count$n.self == 1) {
#            ## single self-response and random effect -- computational trick
#            ## NOT possible to currently use given normalization of model
#            tmp.sigma <- sqrt(tmp.sigma.self^2+tmp.sigma.re^2)
#          }

          
          if (!do.gr) {
            if (options$debug > 1) { cat("ll.oprobit2\n") }
            val <- ll.oprobit(data$y0[,idx.self],
                              xb,
                              tmp.sigma.self,
                              tmp.taus,
                              count$n.cat,
                              options$debug)

            llik2 <- llik2 + sum(val)

            if (options$debug > 1)
              cat("SELF LLIK2 (no RE):",idx.self,val,"\n")

          }
          
        }
        if (do.gr) {
          gval <- gr.oprobit.self.C(data$y0,
                                    xb,
                                    as.numeric(sigma.self)[1],
                                    Staus, ## JW: why is this not tmp.taus
                                    data$v0s,
                                    data$v0s1,
                                    data$x0,
                                    count$n.cat,
                                    count$n.self,
                                    options$debug)
          
          GR$beta   <-   gval$beta
          GR$se.self<-   gval$sigma
          GR$gamma  <-   GR$gamma  + gval$gamma
          GR$gamma1 <-   GR$gamma1 + gval$gamma1
        }
          
      } else { ## WITH random effect
        if (options$int.meth=="gh") {## estimated with gaussian hermite quadrature

          tmp.sigma.self <- as.numeric(sigma.self)[1]
          tmp.sigma.re   <- sigma.re
          
          llik2 <- ll.self.gh(data$y0, xb,
                            tmp.sigma.self,tmp.sigma.re,
                            Staus, gh, count$n.cat, options$debug>0)
          
          if (!fitted) {
            if (options$debug > 1)
              cat("SELF LLIK2 (GH):",llik2,"\n")
          } else {
            self.prob[[1]] <- NULL
          }
          
        } else {  # estimated with adaptive integration

          stop("Only Gaussian hermite quadrature currently operational\n")
          
#          tmp.sigma.self <- as.real(sigma.self)[1]
#          tmp.sigma.re   <- sigma.re
#            
#          llik2 <-ll.self.ai(data$y0,data$x0%*%beta,
#                           tmp.sigma.self,tmp.sigma.re,
#                           tau2,count$n.cat, options,options$debug>0)
#
#          if (!fitted) {
#            if (options$debug > 1)
#              cat("SELF LLIK2 (AI):",llik2,"\n")
#          } else {
#            self.prob[[1]] <- NULL
#          }
        } 
      }
    }
    
    if (options$rprof)
      Rprof(NULL)

    
    
    if (do.gr) {
      ## only keep things which are estimated!
      gf <- function(x,y) {
        if (sum(y) > 0 ) {
          x   <- x[y]
        } else {
          x  <- NULL
        }
        return(x)
      }

      ## zero out non-estimated values!
      GR$gamma  <- gf( GR$gamma , parm$estimated$gamma )
      GR$gamma1 <- gf( GR$gamma1, parm$estimated$gamma1 )
       GR$beta  <- gf( GR$beta , parm$estimated$beta  )
      GR$theta  <- gf( GR$theta, Estimated.theta.all  )
      GR$se.vign<- gf( GR$se.vign, parm$estimated$sigma.vign )
      GR$se.self<- gf( GR$se.self, parm$estimated$sigma.self )
      
      GR.out <- c(GR$gamma1,
                  GR$gamma,
                  GR$se.self,
                  GR$se.vign,
                  GR$theta,
                  GR$beta )
      if (options$debug>1) { 
        cat("GR:\n")
        print(-GR.out)
      }

      if (debug > 1) {
        cat("GR out\n")
        print(GR)
        print(GR.out)
      }
      
      return(-GR.out)
    }

    if (get.LL) {
      LL.self <- -llik2
      if (options$debug > 0) 
        cat("SELF LLIK2 get.ll:",LL.self,"\n")
    }
    
    ## just getting probabilities? end here
    if (get.LL)
      return( list( LL.self=LL.self, LL.vign=LL.vign) )
    if (fitted)
      return( list(self=self.prob,vign=vign.prob))
    
    #### A.2.d OVERALL -Log-Likelihood 
    llik <- -(llik1+llik2)

    if (options$debug > 0) {
      cat("vign -LL :",-llik1,"avg",-llik1/count$nobs.vign,"\n")
      cat("self -LL :",-llik2,"avg",-llik2/count$nobs.self,"\n")
      cat("-LL:",llik,"\n")
    }

    if (is.finite(llik)) {
      return( llik )
    } else {
      if (options$debug > 0) cat("Bad -LL:",llik,"\n")
      return(count$nobs.penalty*200)
    }
    
  } ## END of loglikelihood

  ## GRADIENT
  chopit.gr <- function(b) {
    rval <- chopit.llog(b,do.gr=TRUE,verbose=verbose)
    return(rval)
  }
  
  ## and given (potentially) overwritten theta estimated flags:
  Estimated.theta.all <- NULL
  for (i in 1:count$n.vign.set) {
    zti <- ifelse( count$n.vign.set > 1, i, "")
    txt <- paste("Estimated.theta.all <- c(Estimated.theta.all, parm$estimated$theta",zti,")",sep="")
    eval(parse(text=txt))
  }
  
  
  #############################################################
  ## 
  ## A.6 OPTIMIZING THE FUNCTION
  ## 
  #############################################################
  if (options$debug > 0) 
    cat("anchors.chopit.fit: optimize the likelihood function\n\n")

#  if (normalize != "self")
#    do.oprobit <- FALSE
  
  
  ## testing, testing...
  if (verbose) {
    cat("\nStarting parameters:\n")
    print(parm$start)
    tmp <- as.data.frame(as.matrix(parm$pvec))
    rownames(tmp) <- parm$nvec;
    print(tmp,digits=16)
    cat("\n\n")

#    if (options$analytic) {
#      cat("\nStarting gradients:\n")
#      tmp.GR <- chopit.llog(parm$pvec,do.gr=TRUE);
#      print(as.matrix(tmp.GR),digits=16)
#    }
    
  }

  if (verbose) 
    cat("\n\nStarting chopit() estimation...\n")


  if (!is.null(options$int.gh))
    gh   <- as.matrix(options$int.gh)
  else
    gh <- as.matrix( ghweights(options$int.ghorder) )
  
  if (options$debug > 5) {
    ## we may need GH for the test of the LL
    ## so use the first version... which is where we are going to start anyway

    out <- chopit.llog(parm$pvec);


  rv <- list(data   = data,
             parm   = parm,
             count  = count,
             options= options,
             optim  = out,
             hess   = NULL,
             LL.vign   = NULL,
             LL.self   = NULL,
#             prob      = est.prob,
             gr        = NULL,
             time      = NULL
             )
    
    class(rv) <- "anchors.chopit.fit"
    return(invisible(rv))
  }


  if (options$analytic && options$linear && !options$random ) {
    cat("\nNOTE: analytical gradients are being employed\n\n")
    tmp.chopit.gr <- chopit.gr
  } else {
    cat("\nNOTE: numerical gradients are being employed\n\n")
    tmp.chopit.gr <- NULL
  }

  if (options$optimizer == "genoud") {
    if (verbose) 
      cat("anchors.chopit.fit: entering GENOUD optimization\n")
    
    # loaded automatically by Dependencies: 
	# require(rgenoud) || stop("rgenoud library required to invoke optimizer='genoud' option")

    domain <- options$domain
    
    d.gamma1 <- rep( -domain, sum(c(parm$estimated$gamma1, parm$estimated$gamma)) )
    d.gamma2 <- rep(  domain, sum(c(parm$estimated$gamma1, parm$estimated$gamma)) )
#    print(d.gamma1)
#    print(    parm$start$gamma > 0 )
#    d.gamma1[parm$start$gamma >0 ] <- -domain
    d.gamma  <- as.numeric(rbind(d.gamma1,d.gamma2))

    lb.se <- ifelse(options$linear, 0.001, -5)
    
    d.re     <- rep( c(lb.se, domain), length(parm$estimated$sigma.re))
    d.self   <- rep( c(lb.se, domain), length(parm$estimated$sigma.self))
    d.vign   <- rep( c(lb.se, domain), length(parm$estimated$sigma.vign))
    d.theta  <- rep( c(-domain   , domain), length(parm$estimated$theta)) ## theta2..
    d.beta   <- rep( c(-50  ,50), length(parm$estimated$beta  ))
    
    Domains <- matrix(c(d.gamma,
                        d.re,
                        d.self,
                        d.vign,
                        d.theta,
                        d.beta
                        ), ncol=2,byrow=TRUE)

    if (options$debug > 0) {
      cat("Genoud Domains\n")
      print(Domains)
      cat("Parameter estimated on/off flags\n")
      print( as.matrix(unlist(parm$estimated)))
    }

    
    
    Domains <- Domains[ unlist(parm$estimated) , ]
    

    if (NROW(Domains) != length(parm$pvec)) {
      cat("ERROR! Genoud domain not of correct dimension:", dim(Domains),": length of pvec",length(parm$pvec),"\n")
      print(parm$pvec)
    }
    
    print.oob <- FALSE
    stime <- system.time(out <- 
           genoud(chopit.llog,
                  nvars           = length(parm$pvec),
                  pop.size        = options$pop.size,
                  starting.values = parm$pvec,
                  gr              = tmp.chopit.gr, 
                  wait.generations= options$wait.generations,
                  MemoryMatrix    = options$MemoryMatrix,
                  Domains         = Domains,
                  BFGS            = options$BFGS,
                  hessian         = options$hess,
                  max.generations = options$max.generations,
                  print.level     = options$print.level,
                  ))
    print.oob <- options$print.oob
                  
                  
  } else {
          
    print.oob <- FALSE
    stime <- system.time(out <- optim(par = parm$pvec,
                                      fn  = chopit.llog,
                                      gr  = tmp.chopit.gr, 
                                      method =options$optim.method,
                                      control=list(trace=options$trace,
                                                   maxit=options$maxit,
                                                   reltol=options$reltol),
                                      hessian=options$hess
                                      ))
    print.oob <- options$print.oob
  }
  
  final.gr <- NULL
  if (options$analytic) {
    print.oob <- TRUE
    final.gr <-  chopit.gr(out$par)
    if (verbose)  {
      cat("FINAL GRADIENTS\n")
      print(final.gr)
    }
  }
  
  if (options$debug > 0)  {
    cat("\nCompleted fitting:\n")
    cat("optim() estimation time:",stime,"\n")
  }

  ## Save final parameters into parm list
  parm$pvec <- out$par
  parm <- packv(unpackv(parm))

  info.se <- NULL
  ## now get SE from information matrix:
  chopit.hess <- NULL
  if (options$hess==TRUE) {
    chopit.hess <- out$hessian
  }

  if (options$debug > 0)
    cat("anchors.chopit.fit: go get final LL")

  get.LL <- TRUE
  out.LL <- chopit.llog(parm$pvec)
  get.LL <- FALSE  


  rv <- list(data   = data,
             parm   = parm,
             count  = count,
             options= options,
             optim  = out,
             hess   = chopit.hess,
             LL.vign   = out.LL$LL.vign,
             LL.self   = out.LL$LL.self,
#             prob      = est.prob,
             gr        = final.gr,
             time      = stime
             )

  class(rv) <- "anchors.chopit.fit"
  
  return(invisible(rv))

}

