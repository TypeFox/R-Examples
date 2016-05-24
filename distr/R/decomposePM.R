#######################################################################

setMethod("decomposePM", "AbscontDistribution",
      function(object){
         ep <- getdistrOption("TruncQuantile")
         p0 <- as.vector(object@p(0))
         p01 <- if(.inArgs("lower.tail",object@p))
                object@p(0, lower.tail = FALSE) else 1-p0
         neg <- if(p0>ep){
                qnew.n <- function(p, lower.tail = TRUE,
                               log.p = FALSE){}
                body(qnew.n) <- substitute({
                              if (log.p) p <- exp(p)
                              if(!lower.tail) p <- 1-p
                              if (any((p < -.Machine$double.eps)|
                                      (p > 1+.Machine$double.eps)))
                                  warning(gettextf(
                                  "q method of %s produced NaN's ", objN))
                              i01 <- (-.Machine$double.eps<=p)&
                                     (p<=1+.Machine$double.eps)
                              p0x <- p[i01] ## only values in [0,1] are used
                              q0  <- p*0
                              q0[!i01] <- NaN
                              q0[ i01] <- pmin(object@q(p0x*p0),0)
                              return(as.numeric(q0))},
                            list(objN = quote(.getObjName(1))))
                AbscontDistribution(
                    p = function(q, lower.tail = TRUE, log.p = FALSE){
                               px <- object@p(q, lower.tail = lower.tail,
                                          log.p = FALSE)
                               if(lower.tail)
                                  px <- pmin(p0,px)/p0
                               else
                                  px <- pmax(0,px-p01)/p0
                               if(log.p) px <- log(px)
                               return(px)},
                    d = function(x, log = FALSE){
                               dx <- object@d(x, log = log)
                               d0 <- (x<=0)
                               if(log) dx+log(d0)-log(p0)
                                  else dx*d0/p0},
                    q = qnew.n,
                    r = function(n) qnew.n(runif(n))
                    )} else -abs(Norm())
         pos <- if(p01>ep){
                qnew <- function(p, lower.tail = TRUE,
                               log.p = FALSE){}
                body(qnew) <- substitute({
                              if (log.p) p <- exp(p)
                              if(!lower.tail) p <- 1-p
                              if (any((p < -.Machine$double.eps)|
                                      (p > 1+.Machine$double.eps)))
                                  warning(gettextf(
                                  "q method of %s produced NaN's ", objN))
                              i01 <- (-.Machine$double.eps<=p)&
                                     (p<=1+.Machine$double.eps)
                              p0x <- p[i01] ## only values in [0,1] are used
                              q0  <- p*0
                              q0[!i01] <- NaN
                              q0[ i01] <- pmax(object@q(pmin(p0+p0x*p01,1)),0)
                              return(as.numeric(q0))},
                            list(objN = quote(.getObjName(1))))
                AbscontDistribution(
                    p = function(q, lower.tail = TRUE, log.p = FALSE){
                               px <- object@p(q, lower.tail = lower.tail,
                                          log.p = FALSE)
                               if(lower.tail)
                                  px <- pmax(0,px-p0)/p01
                               else
                                  px <- pmin(p01,px)/p01
                               if(log.p) px <- log(px)
                               return(px)},
                    d = function(x, log = FALSE){
                               dx <- object@d(x, log = log)
                               d0 <- (x>=0)
                               if(log) dx+log(d0)-log(p01)
                                  else dx*d0/p01},
                    q = qnew,
                    r = function(n) qnew(runif(n))
                    )} else abs(Norm())
         list(
         "neg"=list(D = neg, w = p0*(p0>ep)),
         "pos"=list(D = pos, w = p01*(p01>ep))
          )})

setMethod("decomposePM", "DiscreteDistribution",
      function(object){
         supp <- support(object)

         suppp <- supp[supp>0]
         suppm <- supp[supp<0]

         pp <- max(1-p(object)(0),0)
         pm <- max(p(object)(0)-d(object)(0),0)

         p0 <- max(d(object)(0),0)

         Dm <- if(length(suppm))
               DiscreteDistribution(supp = suppm, prob = d(object)(suppm)/pm)
         else  Dirac(-1)
         Dp <- if(length(suppp))
              DiscreteDistribution(supp = suppp, prob = d(object)(suppp)/pp)
         else Dirac(1)
         D0 <- Dirac(0)
         list("neg" = list( D = Dm, w = pm),
              "0"   = list( D = D0, w = p0),
              "pos" = list( D = Dp, w = pp))
          })

setMethod("decomposePM", "UnivarLebDecDistribution", function(object){
         p.a <- acPart(object)
         p.d <- discretePart(object)

         w.a <- acWeight(object)
         w.d <- discreteWeight(object)

         p.a.DC <- if(w.a) decomposePM(p.a) else
                       list(pos = list(D =  Exp(), w = 1/2),
                            neg = list(D = -Exp(), w = 1/2))
         p.d.DC <- if(w.d) decomposePM(p.d) else
                       list(pos = list(D =  Dirac(1),  w = 0),
                            neg = list(D =  Dirac(-1), w = 0),
                            "0" = list(D = Dirac(0), w = 1))


         w.p0 <- p.d.DC$pos$w * w.d + p.a.DC$pos$w * w.a
         w.m0 <- p.d.DC$neg$w * w.d + p.a.DC$neg$w * w.a
         w.p <- w.p0 + (w.p0==0)
         w.m <- w.m0 + (w.m0==0)

         pos <- UnivarLebDecDistribution(
                 discretePart = p.d.DC$pos$D,
                 acPart =  p.a.DC$pos$D,
                 discreteWeight = max(min(p.d.DC$pos$w*w.d/w.p,1),0))
         neg <- UnivarLebDecDistribution(
                 discretePart = p.d.DC$neg$D,
                 acPart =  p.a.DC$neg$D,
                 discreteWeight = max(min(p.d.DC$neg$w*w.d/w.m,1),0))
         list( "pos" = list( D = pos, w = w.p0),
               "neg" = list( D = neg, w = w.m0),
               "0" = list( D = Dirac(0), w = 1-w.m0-w.p0))

})

