setMethod("*", c("AcDcLcDistribution","AcDcLcDistribution"),
function(e1,e2){

 e1 <- .ULC.cast(e1)
 e2 <- .ULC.cast(e2)

#         if( is(e1,"AffLinUnivarLebDecDistribution"))
#             e1 <- as(e1, "UnivarLebDecDistribution")
#         if( is(e2,"AffLinUnivarLebDecDistribution"))
#             e2 <- as(e2, "UnivarLebDecDistribution")
#
#         if( is(e1,"AbscontDistribution"))
#             e1 <- as(as(e1,"AbscontDistribution"), "UnivarLebDecDistribution")
#         if( is(e2,"AbscontDistribution"))
#             e2 <- as(as(e2,"AbscontDistribution"), "UnivarLebDecDistribution")
#         if(is(e1,"DiscreteDistribution"))
#             e1 <- as(as(e1,"DiscreteDistribution"), "UnivarLebDecDistribution")
#         if(is(e2,"DiscreteDistribution"))
#             e2 <- as(as(e2,"DiscreteDistribution"), "UnivarLebDecDistribution")

         ep <- getdistrOption("TruncQuantile")

         e1DC <- decomposePM(e1)
         e2DC <- decomposePM(e2)
         w12pp <- e1DC$pos$w*e2DC$pos$w
         w12mm <- e1DC$neg$w*e2DC$neg$w
         w12pm <- e1DC$pos$w*e2DC$neg$w
         w12mp <- e1DC$neg$w*e2DC$pos$w
         mixCoeff <- c(w12pp,w12mm,w12pm,w12mp)
         mixCoeff <- c(mixCoeff,1-sum(mixCoeff))

         e12pp <- if(w12pp>ep)
                       as(exp(log(e1DC$pos$D)+log(e2DC$pos$D)),
                          "UnivarLebDecDistribution")
                  else as(Dirac(1), "UnivarLebDecDistribution")

         e12mm <- if(w12mm>ep)
                       as(exp(log(-e1DC$neg$D)+log(-e2DC$neg$D)),
                          "UnivarLebDecDistribution")
                  else as(Dirac(1), "UnivarLebDecDistribution")
         e12pm <- if(w12pm>ep)
                       as(-exp(log(e1DC$pos$D)+log(-e2DC$neg$D)),
                          "UnivarLebDecDistribution")
                  else as(Dirac(-1), "UnivarLebDecDistribution")
         if(identical(e1,e2)) e12mp <- e12pm
         else e12mp <- if(w12mp>ep)
                            as( -exp(log(-e1DC$neg$D)+log(e2DC$pos$D)),
                                 "UnivarLebDecDistribution")
                       else as(Dirac(-1), "UnivarLebDecDistribution")

         e12pm <- .del0dmixfun(e12pm)
         e12mp <- .del0dmixfun(e12mp)

         obj <- flat.LCD(mixCoeff = mixCoeff,
                         e12pp, e12mm, e12pm, e12mp,
                         as(Dirac(0),"UnivarLebDecDistribution"))

         if(getdistrOption("simplifyD"))
            obj <- simplifyD(obj)

         rnew <- function(n, ...){}
         body(rnew) <- substitute({ g1(n, ...) * g2(n, ...) },
                                    list(g1 = e1@r, g2 = e2@r))
         obj@r <- rnew

          if(is(e1@Symmetry,"SphericalSymmetry")&&
             is(e2@Symmetry,"SphericalSymmetry"))
             if(.isEqual(SymmCenter(e1@Symmetry),0) && 
                .isEqual(SymmCenter(e2@Symmetry),0))
                 obj@Symmetry <-  SphericalSymmetry(0)   

         return(obj)
         })

setMethod("/", c("numeric",
                 "AcDcLcDistribution"),
function(e1,e2){
         e2s <- as.character(deparse(match.call(
                call = sys.call(sys.parent(1)))$e2))

 e2 <- .ULC.cast(e2)

#         if( is(e2,"AffLinUnivarLebDecDistribution"))
#             e2 <- as(e2, "UnivarLebDecDistribution")

#         if( is(e2,"AbscontDistribution"))
#             e2 <- as(as(e2, "AbscontDistribution"),
#                      "UnivarLebDecDistribution")

#         if( is(e2,"DiscreteDistribution"))
#             e2 <- as(as(e2, "DiscreteDistribution"),
#                      "UnivarLebDecDistribution")

         if (discreteWeight(e2)>getdistrOption("TruncQuantile"))
            if (d.discrete(e2)(0)>getdistrOption("TruncQuantile"))
                stop(gettextf("1 / %s is not well-defined with positive probability ", e2s))

         e2DC <- decomposePM(e2)
         w2p <- e2DC$pos$w
         w2m <- e2DC$neg$w

         e2p <- as(exp(-log(e2DC$pos$D)), "UnivarLebDecDistribution")
         e2m <- as(-exp(-log(-e2DC$neg$D)), "UnivarLebDecDistribution")

         e2D <- flat.LCD(mixCoeff = c(w2p, w2m), e2p, e2m)

         if(getdistrOption("simplifyD"))
            e2D <- simplifyD(e2D)

         obj <- e1*e2D

         rnew <- function(n, ...){}
         body(rnew) <- substitute({ g1 / g2(n, ...) },
                                    list(g1 = e1, g2 = e2@r))
         obj@r <- rnew

          if(is(e2@Symmetry,"SphericalSymmetry"))
             if(.isEqual(SymmCenter(e2@Symmetry),0))
                  obj@Symmetry <-  SphericalSymmetry(0)   

         return(obj)
         })

setMethod("/", c("AcDcLcDistribution",
                 "AcDcLcDistribution"),
function(e1,e2){
         e2s <- as.character(deparse(match.call(
                call = sys.call(sys.parent(1)))$e2))
#         if( is(e2,"AbscontDistribution"))
#             e2 <- as(as(e2,"AbscontDistribution"), "UnivarLebDecDistribution")

#         if( is(e2,"DiscreteDistribution"))
#             e2 <- as(as(e2,"DiscreteDistribution"), "UnivarLebDecDistribution")

 e2 <- .ULC.cast(e2)

         if (discreteWeight(e2)>getdistrOption("TruncQuantile"))
         if (d.discrete(e2)(0)>getdistrOption("TruncQuantile"))
            stop(gettextf("1 / %s is not well-defined with positive probability ", e2s))

         obj <- e1 * (1/e2)

         rnew <- function(n, ...){}
         body(rnew) <- substitute({ g1(n, ...) / g2(n, ...) },
                                    list(g1 = e1@r, g2 = e2@r))
         obj@r <- rnew

          if(is(e1@Symmetry,"SphericalSymmetry")&&
             is(e2@Symmetry,"SphericalSymmetry"))
             if(.isEqual(SymmCenter(e1@Symmetry),0) && 
                .isEqual(SymmCenter(e2@Symmetry),0))
                 obj@Symmetry <-  SphericalSymmetry(0)   

         return(obj)
         })

setMethod("^", c("AcDcLcDistribution","Integer"),
function(e1,e2){

           if(.isEqual(e2,0)) return(Dirac(1))
           if(.isEqual(e2,1)) return(e1)

           ep <- getdistrOption("TruncQuantile")
           d00 <- discretePart(e1)@d(0)
           d0 <- discreteWeight(e1)*d00
           if(d0 > ep){
                if(.isEqual(d00,1)){
                   e1 <- acPart(e1)
                }else{
                   su <- support(discretePart(e1))
                   pr <- d(discretePart(e1))(su)
                   acW <- acWeight(e1)
                   discreteP <- DiscreteDistribution(
                                supp = su[su!=0],
                                prob = pr[su!=0]/(1-d00))
                   e1 <- UnivarLebDecDistribution(acPart = acPart(e1),
                         discretePart = discreteP, acWeight = acW)
                }
               }

           e1DC <- decomposePM(e1)
           mixCoeff <- c(e1DC$pos$w,e1DC$neg$w)
           mixCoeff <- mixCoeff/sum(mixCoeff)
           e1p <- if(mixCoeff[1]>ep)
                    as(exp(e2*log(e1DC$pos$D)),"UnivarLebDecDistribution")
                 else as(Dirac(1), "UnivarLebDecDistribution")
           e1m <- if(mixCoeff[2]>ep)
                    as((-1)^e2*exp(e2*log(-e1DC$neg$D)),"UnivarLebDecDistribution")
                 else as(Dirac((-1)^e2), "UnivarLebDecDistribution")
           erg <- flat.LCD(mixCoeff = mixCoeff, e1p, e1m)

#
           if(d0 > ep){
                if(.isEqual(d00,1)){
                    erg <- UnivarLebDecDistribution(acPart = acPart(erg),
                           discretePart = Dirac(0), acWeight = acW)
                }else{
                    su <- support(discretePart(erg))
                    su0 <- c(su,0)
                    o <- order(su0)
                    pr <- c(d(discretePart(erg))(su) * (1-d00), d00)
                    suo <- su0[o]
                    pro <- pr[o]
                    discreteP <- DiscreteDistribution(supp = suo, prob = pro)
                    erg <- UnivarLebDecDistribution(acPart = acPart(erg),
                           discretePart = discreteP, acWeight = acW)
                }
             }
           if(getdistrOption("simplifyD"))
                erg <- simplifyD(erg)

           rnew <- function(n, ...){}
           body(rnew) <- substitute({ g1(n, ...)^g2 },
                                    list(g1 = e1@r, g2 = e2))
           erg@r <- rnew
           return(erg)
           })

setMethod("^", c("AcDcLcDistribution","numeric"),
function(e1,e2){
  if (is(try(mc <- match.call(call = sys.call(sys.parent(1))),
         silent=TRUE), "try-error"))
      {e1s <- "e1"; e2s <- "e2"}
  else {e1s <- as.character(deparse(mc$e1))
        e2s <- as.character(deparse(mc$e2))}

  if (length(e2)>1) stop("length of operator must be 1")
  if (isTRUE(all.equal(e2,1))) return(e1)
  if (isTRUE(all.equal(e2,0))) return(Dirac(1))


 e1 <- .ULC.cast(e1)
#  if( is(e1,"AbscontDistribution") || is(e1,"DiscreteDistribution") ||
#      is(e1,"AffLinUnivarLebDecDistribution"))
#      e1 <- as(e1, "UnivarLebDecDistribution")

  if (e2<0) return((1/e1)^(-e2))

  if (.isNatural(e2, tol = 1e-10))
      return(get("^")(e1 = e1, e2 = as(e2,"Integer")))

  ep <- getdistrOption("TruncQuantile")
  d00 <- discretePart(e1)@d(0)
  d0 <- discreteWeight(e1)*d00
  p0 <- p(e1)(0)

  if ((p0 > ep && e2 < 0) || (p0 > d0+ep))
      stop(gettextf("%s^%s is not well-defined with positive probability ",
                    e1s, e2s))

  ### special treatment if e2>=0 and d.discrete(e1)>0
  if(d0 > ep){
     if(.isEqual(d00,1)){
        e1 <- acPart(e1)
     }else{
        su <- support(discretePart(e1))
        pr <- d(discretePart(e1))(su)
        acW <- acWeight(e1)#/(1-d0)
        discreteP <- DiscreteDistribution(
                     supp = su[su!=0],
                     prob = pr[su!=0]/(1-d00))
        e1 <- UnivarLebDecDistribution(acPart = acPart(e1),
              discretePart = discreteP, acWeight = acW)
     }
   }

   erg <- exp( e2 * log(e1))

   ### special treatment if e2>=0 and d.discrete(e1)>0
   if(d0 > ep){
      if(.isEqual(d00,1)){
          erg <- UnivarLebDecDistribution(acPart = acPart(erg),
                 discretePart = Dirac(0), acWeight = acW)
      }else{
#      d1 <- d0 # 1-(1-d0)^e2
          acW <- acWeight(erg) #* (1-d1)
          su <- support(discretePart(erg))
          su0 <- c(su,0)
          o <- order(su0)
          pr <- c(d(discretePart(erg))(su) * (1-d00), d00)
          suo <- su0[o]
          pro <- pr[o] #/(1-acW)
          discreteP <- DiscreteDistribution(supp = suo, prob = pro)
          erg <- UnivarLebDecDistribution(acPart = acPart(erg),
                      discretePart = discreteP, acWeight = acW)
     }
   }

  if(getdistrOption("simplifyD"))
            erg <- simplifyD(erg)

  rnew <- function(n, ...){}
  body(rnew) <- substitute({ g1(n, ...)^g2 },
                            list(g1 = e1@r, g2 = e2))
  erg@r <- rnew

  return(erg)
  }
)

setMethod("^", c("AcDcLcDistribution","Dirac"),
function(e1,e2)e1^location(e2))

setMethod("^", c("AcDcLcDistribution","AcDcLcDistribution"),
function(e1,e2){
 ### check if there are problems
  if (is((e1s <- as.character(deparse(match.call(
                call = sys.call(sys.parent(1)))$e1))), "try-error"))
      e1s <- "e1"
  if (is((e2s <- as.character(deparse(match.call(
                call = sys.call(sys.parent(1)))$e2))), "try-error"))
      e2s <- "e2"

# if( is(e1,"AffLinUnivarLebDecDistribution"))
#     e1 <- as(e1, "UnivarLebDecDistribution")
# if( is(e2,"AffLinUnivarLebDecDistribution"))
#     e2 <- as(e2, "UnivarLebDecDistribution")

# if( is(e1,"AbscontDistribution"))
#     e1 <- as(as(e1,"AbscontDistribution"), "UnivarLebDecDistribution")
# if( is(e2,"AbscontDistribution"))
#     e2 <- as(as(e2,"AbscontDistribution"), "UnivarLebDecDistribution")


# if( is(e1,"DiscreteDistribution"))
#     e1 <- as(as(e1,"DiscreteDistribution"), "UnivarLebDecDistribution")
# if( is(e2,"DiscreteDistribution"))
#     e2 <- as(as(e2,"DiscreteDistribution"), "UnivarLebDecDistribution")

 e1 <- .ULC.cast(e1)
 e2 <- .ULC.cast(e2)


 ep <- getdistrOption("TruncQuantile")

 if(p(e2)(0)-discreteWeight(e2)*d.discrete(e2)(0)>ep)
    { ## must be able to work with negative exponents

         if (d.discrete(e1)(0)*discreteWeight(e1) > ep)
            stop(gettextf("%s^%s is not well-defined with positive probability ",
                 e1s, e2s))

         if ((discreteWeight(e2)>1-ep) && all(.isInteger(support(e2)))){
              Dlist <- lapply(support(e2), function(x)
                  as(do.call("^",list(e1=e1,e2=x)), "UnivarLebDecDistribution"))
              erg <- as(simplifyD( do.call(flat.LCD,
                        c(Dlist, alist(mixCoeff = d.discrete(e2)(support(e2)))))),
                        "UnivarLebDecDistribution")
              if(getdistrOption("simplifyD")) erg <- simplifyD(erg)
              return(erg)
              }

         if (p(e1)(0) > ep)
            stop(gettextf("%s^%s is not well-defined with positive probability ",
                 e1s, e2s))
    }

 if(p(e1)(0)>ep)
    { ## works only for purely natural e2


         if ((discreteWeight(e2)>1-ep) && all(.isInteger(support(e2)))){
              Dlist <- lapply(support(e2), function(x)
                  as(do.call("^",list(e1=e1,e2=x)), "UnivarLebDecDistribution"))
              erg <- as(simplifyD( do.call(flat.LCD,
                        c(Dlist, alist(mixCoeff = d.discrete(e2)(support(e2)))))),
                        "UnivarLebDecDistribution")
              if(getdistrOption("simplifyD")) erg <- simplifyD(erg)
              return(erg)
              }

         stop(gettextf("%s^%s is not well-defined with positive probability ",
                 e1s, e2s))
    }

 le1 <- log(e1)
 le <- le1 * e2
 erg <- exp(le)
 if(getdistrOption("simplifyD")) erg <- simplifyD(erg)

 rnew <- function(n, ...){}
 body(rnew) <- substitute({ g1(n, ...)^g2(n, ...) },
                            list(g1 = e1@r, g2 = e2@r))
 erg@r <- rnew

 return(erg)
})


setMethod("^", c("numeric","AcDcLcDistribution"),
function(e1,e2){
 ### check if there are problems
  if (is((e1s <- as.character(deparse(match.call(
                call = sys.call(sys.parent(1)))$e1))), "try-error"))
      e1s <- "e1"
  if (is((e2s <- as.character(deparse(match.call(
                call = sys.call(sys.parent(1)))$e2))), "try-error"))
      e2s <- "e2"

 e2 <- .ULC.cast(e2)
 #e2 <- .if( is(e2,"AffLinUnivarLebDecDistribution"))
 #    e2 <- as(e2, "UnivarLebDecDistribution")
 #if( is(e2,"AbscontDistribution"))
 #    e2 <- as(as(e2,"AbscontDistribution"), "UnivarLebDecDistribution")
 #if( is(e2,"DiscreteDistribution"))
 #    e2 <- as(as(e2,"DiscreteDistribution"), "UnivarLebDecDistribution")

 ep <- getdistrOption("TruncQuantile")
 if(p(e2)(0)-discreteWeight(e2)*d.discrete(e2)(0)>ep)
    { ## must be able to work with negative exponents

         if (abs(e1) < ep)
             stop(gettextf("%s^%s is not well-defined with positive probability ",
                 e1s, e2s))

         if ((discreteWeight(e2)>1-ep) && all(.isInteger(support(e2)))){
              erg <- DiscreteDistribution(e1^support(e2),
                                          d.discrete(e2)(support(e2)))
              if(!getdistrOption("simplifyD"))
                  erg <- as(erg,"UnivarLebDecDistribution")
              return(erg)
             }

         if (e1 < -ep)
            stop(gettextf("%s^%s is not well-defined with positive probability ",
                 e1s, e2s))
    }
 if(e1< -ep)
    { ## works only for purely natural e2

         if ((discreteWeight(e2)>1-ep) && all(.isInteger(support(e2)))){
              erg <- DiscreteDistribution(e1^support(e2),
                                          d.discrete(e2)(support(e2)))
             if(!getdistrOption("simplifyD"))
                 erg <- as(erg,"UnivarLebDecDistribution")
             return(erg)
             }

         stop(gettextf("%s^%s is not well-defined with positive probability ",
                 e1s, e2s))
    }
  le1 <- log(e1)
  le <- le1 * e2
  erg <- exp(le)
  if(getdistrOption("simplifyD")) erg <- simplifyD(erg)

  rnew <- function(n, ...){}
  body(rnew) <- substitute({ g1^g2(n, ...) },
                            list(g1 = e1, g2 = e2@r))
  erg@r <- rnew

  return(erg)
})

setMethod("+", signature(e1="AcDcLcDistribution", e2="AcDcLcDistribution"),
           function(e1,e2)(.ULC.cast(e1)+(-.ULC.cast(e2))))
setMethod("-", signature(e1="AcDcLcDistribution", e2="AcDcLcDistribution"),
           function(e1,e2)(.ULC.cast(e1)+(-.ULC.cast(e2))))


setMethod("sign", "AcDcLcDistribution",
            function(x){ 
            if(is(x,"AbscontDistribution")) d0 <-0
            else if(is(x,"DiscreteDistribution")) d0 <- d(x)(0)
            else d0 <- d.discrete(as(x,UnivarLebDecDistribution))(0)
            pm <- p(x)(-getdistrOption("TruncQuantile"))
            pp <- p(x)(getdistrOption("TruncQuantile"), lower=FALSE)
            Symmetry <- NoSymmetry()
            if(.isEqual(pm,pp))
               Symmetry <- SphericalSymmetry(0)
            DiscreteDistribution(supp = c(-1,0,1), prob = c(pm, d0, pp))
            })

setMethod("sqrt", "AcDcLcDistribution",
            function(x) x^0.5)

setMethod("Math", "AcDcLcDistribution",
          function(x) callGeneric(.ULC.cast(x)))
