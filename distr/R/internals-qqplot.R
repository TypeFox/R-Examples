################################################################
# QQ - Plot helper functions in package distr
################################################################


## helpers
.inGaps <- function(x,gapm){
  if(is.null(gapm)) return(rep(FALSE,length(x)))
  fct <- function(x,m){ m[,2]>=x & m[,1]<=x}
  sapply(x, function(y) length(which(fct(y,gapm)))>0)
}

.isReplicated <- function(x, tol=.Machine$double.eps){
  tx <- table(x)
  rx <- as.numeric(names(tx[tx>1]))
  sapply(x, function(y) any(abs(y-rx)<tol))
}

.NotInSupport <- function(x,D){
  if(length(x)==0) return(logical(0))
  nInSupp <- which(x < q(D)(0))
  nInSupp <- unique(sort(c(nInSupp,which(x > q(D)(1)))))

  nInSuppo <-
      if("support" %in% names(getSlots(class(D))))
         which( ! x %in% support(D)) else numeric(0)
  if("gaps" %in% names(getSlots(class(D)))){
         InGap <- which( .inGaps(x,gaps(D)))
         if("support" %in% names(getSlots(class(D))))
            nInSupp <- unique(sort(c(nInSupp, intersect(InGap,nInSuppo))))
         else
            nInSupp <- unique(sort(c(nInSupp, InGap)))
  }else{
         nInSupp <- unique(sort(c(nInSupp, nInSuppo)))
  }
  return((1:length(x)) %in% nInSupp)
}


.SingleDiscrete <- function(x,D){
  ## produces a logical vector of
  ##     0  : discrete mass point
  ##     1  : within continuous support
  ##     2  : left gap point
  ##     3  : right gap point
  ##     4  : not in support
  lx <- x * 0

  lx[.NotInSupport(x,D)] <- 4

  idx.0 <- ((x>q(D)(1)) | (x<q(D)(0)))
  iG <- rep(FALSE,length(x))

  if(is(D, "DiscreteDistribution")){
     return(lx)
  }
  if("gaps" %in% names(getSlots(class(D)))){
     if(!is.null(gaps(D))){
        lx[apply(sapply(gaps(D)[,1], function(u) .isEqual(u,x)),1,any)] <- 2
        lx[apply(sapply(gaps(D)[,2], function(u) .isEqual(u,x)),1,any)] <- 3
        iG <- .inGaps(x,gaps(D))
        lx[!idx.0 & !iG] <- 1
     }else{
        lx[!idx.0 & !iG] <- 1
     }
  }
  if("support" %in% names(getSlots(class(D)))){
     idx <- x %in% support(D)
     if("acPart" %in% names(getSlots(class(D))))
         idx.0 <- ((x>q.ac(D)(1)) | (x<q.ac(D)(0)))
     lx[idx & (idx.0|iG)] <- 0
  }

  return(lx)
}


.makeLenAndOrder <- function(x,ord){
   n <- length(ord)
   x <- rep(x, length.out=n)
   x[ord]
}


.pk2 <- if(getRversion()<"2.16.0") function(p0, n){
                 .C("pkolmogorov2x", p = as.double(p0),
                     as.integer(n), PACKAGE = "stats")$p
        }else function(p0,n){
                 .Call("pKolmogorov2x", p0, n) #, PACKAGE = "stats")
        }
.pks2 <- if(getRversion()<"2.16.0") function(x, tol){
                 .C("pkstwo", as.integer(1),
                    p = as.double(x), as.double(tol), PACKAGE = "stats")$p
        }else function(x, tol){
                 .Call("pKS2", p = x, tol) #, PACKAGE = "stats")
        }


.q2kolmogorov <- function(alpha,n,exact=(n<100), silent0 = TRUE){ ## Kolmogorovstat
 if(is.numeric(alpha)) alpha <- as.vector(alpha)
 else stop("Level alpha must be numeric.")
 if(any(is.na(alpha))) stop("Level alpha must not contain missings.")
 if(exact){
   fct <- function(p0){
   ### from ks.test from package stats:
      .pk2(p0,n) -alpha
   }
   i <- 0
   oK <- FALSE
   del <- 0.01
   while(!oK && i < 20){
       i <- i + 1
       res <- try(uniroot(fct,lower=del,upper=3*(1-del)/sqrt(n))$root*sqrt(n), silent=silent0)
       del <- del / 10
       if(!is(res, "try-error")) oK <- TRUE
   }
 }else{
   fct <- function(p0){
   ### from ks.test from package stats:
        1 - .pks2(p0,1e-09)-alpha  }
   i <- 0
   oK <- FALSE
   del <- 0.01
   while(!oK && i < 20){
       i <- i + 1
       res <- try(uniroot(fct,lower=del,upper=3*(1-del))$root, silent=silent0)
       del <- del / 10
       if(!is(res, "try-error")){
          oK <- TRUE
       }else{
          if(!silent0){
              cat("Problem in uniroot in .q2kolmogorov:\nSearch bounds were")
              print(c(lower=del,upper=sqrt(n)*(1-del)))
          }
       }
   }
 }
 return(res)
}

.BinomCI.in <- function(t,p.bi,x.i, del.i=0,D.i,n.i,alpha.i){
   p.bi.u <- p(D.i)(x.i+(t+del.i)/sqrt(n.i))
   p.bi.l <- p.l(D.i)(x.i-(t-del.i)/sqrt(n.i))
   d.r <- if(n.i*p.bi>floor(n.i*p.bi)) 0 else
        dbinom(x = n.i*p.bi, size = n.i, prob = pmax(p.bi.u,1))
   p.r <- pbinom(q = n.i*p.bi, size = n.i, prob = pmin(p.bi.u,1),lower.tail=FALSE)+d.r
   p.l <- pbinom(q = n.i*p.bi, size = n.i, prob = pmax(p.bi.l,0),lower.tail=FALSE)
   r <- p.r -p.l - alpha.i
#   print(c(r=r,p=p.bi,x=x.i,p.u=p.bi.u,p.l=p.bi.l,r.r=p.r,r.l=p.l,t=t,np=n*p.bi))
   r
  }


.BinomCI <- function(x,p.b,D,n,alpha,silent0 = TRUE){
  if(length(x)==0) return(NA)
  res <- sapply(1:length(x),
           function(i){
              .nm.i <- sqrt(n)*max(x[i],n-x[i]) + 1
              .fct.i <- function(t){
                         .BinomCI.in(t, p.bi = p.b[i], x.i = x[i], del.i = 0,
                                     D.i = D, n.i = n, alpha.i = alpha)
              }
              ii <- 0
              oK <- FALSE
              del <- 0.01
              while(ii < 20 && !oK){
                 ii <- ii + 1
                 res0 <- try(uniroot(.fct.i, tol = 1e-9, lower = del,
                                upper=.nm.i*(1-del))$root,
                             silent=silent0)
                 del <- del / 10
                 if(!is(res0, "try-error")){
                     oK <- TRUE
                 }else{
                     if(!silent0){
                         cat("Problem in uniroot in .BinomCI:\nSearch bounds were")
                         print(c(lower=del,upper=.nm.i*(1-del)))
                         cat("Further settings:\n")
                         print(c(i=i,p.bi = p.b[i], x.i = x[i], del=del))
                     }
                     res0 <- NA
                 }
              }
              return(res0)
           }
         )
  return(cbind(left=-res, right=res))
}

.BinomCI.nosym <- function(x,p.b,D,n,alpha,silent0 = TRUE){
  if(length(x)==0) return(NA)
  res0 <- sapply(1:length(x),
            function(i){
              .nm.i <- max(x[i],n-x[i])*sqrt(n) + 1
              .fct.i.d <- function(t,del.0){
                         .BinomCI.in(t, p.bi = p.b[i], x.i = x[i], del.i = del.0,
                                     D.i = D, n.i = n, alpha.i = alpha)
              }
              get.t <- function(del.o){
                 ii <- 0
                 oK <- FALSE
                 del.l <- 0.01
                 while(ii < 20 && !oK){
                    ii <- ii + 1
                    res0 <- try(uniroot(.fct.i.d, tol = 1e-9, lower = del.l,
                                        upper = .nm.i*(1-del.l),
                                        del.0 = del.o)$root, silent = silent0)
                    del.l <- del.l / 10
                    if(!is(res0, "try-error")){
                       oK <- TRUE
                    }else{
                       if(!silent0){
                           cat("Problem in uniroot in .BinomCI.nosym:\nSearch bounds were")
                           print(c(lower=del.l,upper=.nm.i*(1-del.l)))
                           cat("Further settings:\n")
                           print(c(i=i,p.bi = p.b[i], x.i = x[i], del.i=del.o,
                                   del.l=del.l))
                       }
                       res0 <- NA
                    }
                 }
                 return(res0)
              }
    res <- try(optimize(get.t, lower=-.nm.i,  upper = .nm.i, tol = 1e-9),
               silent = silent0)
    if(is(res,"try-error")){
       if(!silent0){
           cat("Problem in optimize in .BinomCI.nosym:\nSearch bounds were")
           print(c(lower=-.nm.i,upper=.nm.i))
           cat("Further settings:\n")
           print(c(i=i,p.bi = p.b[i], x.i = x[i]))
       }
       return(c(NA,NA))
    }
    t.o <- res$objective
    del <- res$minimum
    c(left=-t.o+del, right=t.o+del)
    })
   return(t(res0))
}


.q2pw <- function(x,p.b,D,n,alpha,exact=(n<100),nosym=FALSE, silent0=TRUE){
 if(exact){
    fct <- if(nosym) .BinomCI.nosym else .BinomCI
    ro <- fct(x,p.b,D,n,alpha,silent0)
    return(ro)
 }
 pq <- log(p.b)+log(1-p.b)
 if(is(D,"AbscontDistribution")){
    dp <- d(D)(x,log=TRUE)
    dsupp.p <- dsupp.m <- 1
 }else{ ## have E and sd available ?
    if(!.distrExInstalled) stop("")
    supp.ind <- sapply(x, function(y)
                 which.min(abs(y-support(D))))
    nx <- length(support(D))
    supp.ind.p <- pmax(supp.ind + 1 ,nx)
    supp.ind.m <- pmax(supp.ind - 1 ,1)
    dsupp.p <- support(D)[supp.ind.p] - support(D)[supp.ind]
    dsupp.m <- support(D)[supp.ind] - support(D)[supp.ind.m]
    m <- distrEx::E(D)
    s <- sqrt(distrEx::E(D, fun = function(x) (x-m)^2))
    dp <- log(pnorm((x+dsupp.p/2-m)/s) - pnorm((x-dsupp.m/2-m)/s))
 }
 ro <- exp(pq/2-dp)*(dsupp.p+dsupp.m)/2*qnorm((1+alpha)/2)
 return(cbind(left=-ro,right=ro))
}




.confqq <- function(x,D, datax = FALSE, withConf.pw  = TRUE,
                    withConf.sim = TRUE, alpha,
                    col.pCI, lty.pCI, lwd.pCI, pch.pCI, cex.pCI,
                    col.sCI, lty.sCI, lwd.sCI, pch.sCI, cex.sCI,
                    n,exact.sCI=(n<100),exact.pCI=(n<100), nosym.pCI = FALSE,
                    with.legend = TRUE, legend.bg = "white",
                    legend.pos = "topleft", legend.cex = 0.8,
                    legend.pref = "", legend.postf = "", 
                    legend.alpha = alpha, qqb0=NULL, transf0=NULL, debug = FALSE){

   x <- sort(unique(x))
   if("gaps" %in% names(getSlots(class(D))))
       {if(!is.null(gaps(D)))
            x <- sort(unique(c(x, gaps(D))))
       }
   SI <- .SingleDiscrete(x,D)
#   print(SI)
   SI.in <- SI<4
   SIi <- SI[SI.in]
   SI.c <- SIi>0
   x.in <- x[SI.in]
   x.c <- x.in[SI.c]
   x.d <- x.in[!SI.c]        
   tx.c <- if(is.null(transf0)) x.c else transf0(x.c)
   tx.d <- if(is.null(transf0)) x.d else transf0(x.d)

   qqb <- if(is.null(qqb0)) qqbounds(x,D,alpha,n,withConf.pw, withConf.sim,
                   exact.sCI,exact.pCI,nosym.pCI, debug) else qqb0
                   
   qqb$crit <- if(is.null(qqb0)) qqb$crit[SI.in,]

   if(qqb$err["pw"]){
      if(sum(SI.c)>0){
         if(datax){
            lines(tx.c, qqb$crit[SI.c,"pw.right"],
               col=col.pCI,lty=lty.pCI,lwd=lwd.pCI)
            lines(tx.c, qqb$crit[SI.c,"pw.left"],
               col=col.pCI,lty=lty.pCI,lwd=lwd.pCI)
         }else{
            lines(qqb$crit[SI.c,"pw.right"], tx.c,
               col=col.pCI,lty=lty.pCI,lwd=lwd.pCI)
            lines(qqb$crit[SI.c,"pw.left"], tx.c,
               col=col.pCI,lty=lty.pCI,lwd=lwd.pCI)
         }
      }
      if(sum(!SI.c)>0){
         if(datax){
            points(tx.d, qqb$crit[!SI.c,"pw.right"],
               col=col.pCI, pch=pch.pCI, cex = cex.pCI)
            points(tx.d, qqb$crit[!SI.c,"pw.left"],
               col=col.pCI, pch=pch.pCI, cex = cex.pCI)
         }else{
            points(qqb$crit[!SI.c,"pw.right"], tx.d,
               col=col.pCI, pch=pch.pCI, cex = cex.pCI)
            points(qqb$crit[!SI.c,"pw.left"], tx.d,
               col=col.pCI, pch=pch.pCI, cex = cex.pCI)
         }
      }
   }
   if(qqb$err["sim"]){
      if(sum(SI.c)>0){
         if(datax){
            lines(tx.c, qqb$crit[SI.c,"sim.right"],
               col=col.sCI,lty=lty.sCI,lwd=lwd.sCI)
            lines(tx.c, qqb$crit[SI.c,"sim.left"],
               col=col.sCI,lty=lty.sCI,lwd=lwd.sCI)
         }else{
            lines(qqb$crit[SI.c,"sim.right"], tx.c,
               col=col.sCI,lty=lty.sCI,lwd=lwd.sCI)
            lines(qqb$crit[SI.c,"sim.left"], tx.c,
               col=col.sCI,lty=lty.sCI,lwd=lwd.sCI)
         }
      }
      if(sum(!SI.c)>0){
         if(datax){
            points(tx.d, qqb$crit[!SI.c,"sim.right"],
                col=col.sCI, pch=pch.sCI, cex = cex.sCI)
            points(tx.d, qqb$crit[!SI.c,"sim.left"],
                col=col.sCI, pch=pch.sCI, cex = cex.sCI)
         }else{
            points(qqb$crit[!SI.c,"sim.right"], tx.d,
                col=col.sCI, pch=pch.sCI, cex = cex.sCI)
            points(qqb$crit[!SI.c,"sim.left"], tx.d,
                col=col.sCI, pch=pch.sCI, cex = cex.sCI)
         }
      }
   }
   if(with.legend){
      if( qqb$err["pw"] ||  qqb$err["sim"] ){
         expression1 <- substitute(
            legpf~nosym0~"pointw."~ex.p~alpha==alpha0~"%- conf. interval"~legpof,
            list(legpf = legend.pref, legpof = legend.postf,
                 ex.p = if(exact.pCI) "exact" else "asympt.",
                 alpha0 = round(legend.alpha*100,2),
                 nosym0 = if(nosym.pCI&&exact.pCI) "shortest asymm." else "symm"))
         expression2 <- substitute(
            legpf~"simult."~ex.s~alpha==alpha0~"%- conf. interval"~legpof,
            list(legpf = legend.pref, legpof = legend.postf,
                 ex.s = if(exact.sCI) "exact" else "asympt.",
                 alpha0 = round(legend.alpha*100,2)))

         lcl <- list()
         if(!qqb$err["sim"]){
            expression3 <- expression1
            lcl$pch <- if(sum(!SI.c)>0) pch.pCI else NULL
            lcl$lty <- if(sum(SI.c)>0)  lty.pCI else NULL
            lcl$col <- col.pCI
            lcl$lwd <- if(sum(SI.c)>0)  2 else NULL
         }
         if(!qqb$err["pw"]){
            expression3 <- expression2
            lcl$pch <- if(sum(!SI.c)>0) pch.sCI else NULL
            lcl$lty <- if(sum(SI.c)>0)  lty.sCI else NULL
            lcl$col <- col.sCI
            lcl$lwd <- if(sum(SI.c)>0)  2 else NULL
         }
         if( qqb$err["pw"] && qqb$err["sim"]){
            expression3 <- eval(substitute(expression(expression1, expression2)))
            lcl$pch <- if(sum(!SI.c)>0) c(pch.pCI, pch.sCI) else NULL
            lcl$lty <- if(sum(SI.c)>0)  c(lty.pCI, lty.sCI) else NULL
            lcl$col <- c(col.pCI,col.sCI)
            lcl$lwd <- if(sum(SI.c)>0)  2 else NULL
         }
         do.call(legend, c(list(legend.pos, legend = expression3, bg = legend.bg,
                                merge = FALSE, cex = legend.cex), lcl))
      }
   }
  return(invisible(qqb))
}

.deleteItemsMCL <- function(mcl){
    mcl$n <- NULL
    mcl$col.IdL <- mcl$alpha.CI <- mcl$lty.IdL <-  NULL
    mcl$col.NotInSupport <- mcl$check.NotInSupport <- NULL
    mcl$exact.sCI <- mcl$exact.pCI <- NULL
    mcl$withConf <- mcl$withConf.sim <- mcl$withConf.pw <- NULL
    mcl$withIdLine <- mcl$distance <- NULL
    mcl$col.pCI <- mcl$lty.pCI <- mcl$col.sCI <- mcl$lty.sCI <- NULL
    mcl$lwd.IdL <- mcl$lwd.pCI <- mcl$lwd.sCI <- NULL
    mcl$withLab <- mcl$lab.pts <- mcl$which.lbs <- NULL
    mcl$which.Order <- mcl$order.traf  <- NULL
    mcl$col.pch <- mcl$cex.pch  <- mcl$jit.fac <- NULL
    mcl$col.lbl <- mcl$cex.lbl  <- mcl$adj.lbl <- NULL
    mcl$exp.cex2.pch <- mcl$exp.cex2.lbl <- NULL
    mcl$exp.fadcol.pch <- mcl$exp.fadcol.lbl <- NULL
    mcl$nosym.pCI <- mcl$n.CI <- mcl$n.adj <- NULL
    mcl$legend.cex <- mcl$with.legend <- mcl$legend.bg <- NULL
    mcl$legend.pos <- mcl$legend.pref <- mcl$legend.postf <- NULL
    mcl$legend.alpha <- NULL
mcl}


