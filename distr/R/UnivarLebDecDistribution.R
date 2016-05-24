UnivarLebDecDistribution <- function(acPart, discretePart, acWeight, discreteWeight,
                                     r = NULL, e = NULL, n = NULL, y = NULL){

    if(missing(discretePart) && missing(acPart) && missing(r) )
       stop("at least one of 'discretePart' resp. 'acPart'  resp. 'r' must be given")
    if(missing(discretePart) && missing(acPart))
       {if(missing(e)) e <- getdistrOption("RtoDPQ.e")
        if(missing(n)) n <- getdistrOption("DefaultNrGridPoints")
        return(RtoDPQ.LC(r = r, e = e, n = n, y = y))}

    if(missing(discretePart)) {discreteWeight <- 0; acWeight <- 1;
                               discretePart <- new("Dirac")}
    if(missing(acPart)) {discreteWeight <- 1; acWeight <- 0;
                        acPart <- new("Norm")}

    if(!missing(discretePart)&&!missing(acPart)){
       if(missing(discreteWeight) && missing(acWeight))
             stop("at least one of 'discreteWeight' resp. 'acWeight' must be given")
       if(missing(discreteWeight)) discreteWeight <- 1-acWeight
       if(missing(acWeight)) acWeight <- 1-discreteWeight
       if(discreteWeight <0 || acWeight<0 || acWeight+discreteWeight>1)
             stop("no proper weights given")
    }
    if(discreteWeight > 1 - getdistrOption("TruncQuantile"))
       {return(
           new("UnivarLebDecDistribution", p = discretePart@p,
                r = discretePart@r, d = NULL, q = discretePart@q,
                mixCoeff = c(acWeight = 0, discreteWeight = 1),
                mixDistr = new("UnivarDistrList", list(acPart = acPart,
                                discretePart = discretePart)),
                support = support(discretePart),
                gaps = gaps(acPart)                
              )
             )}
    if(discreteWeight < getdistrOption("TruncQuantile"))
       return(
           new("UnivarLebDecDistribution", p = acPart@p,
                r = acPart@r, d = NULL, q = acPart@q,
                mixCoeff = c(acWeight = 1, discreteWeight = 0),
                mixDistr = new("UnivarDistrList", list(acPart = acPart,
                                discretePart = discretePart)),
                support = support(discretePart),
                gaps = gaps(acPart)
              )
             )

    mixDistr <- new("UnivarDistrList", list(acPart = acPart,
                     discretePart = discretePart))
    mixCoeff <- c(acWeight = acWeight, discreteWeight = discreteWeight)
    rnew <- function(n)
             {U <- rbinom(n, size = 1, prob = acWeight)
              AC <- acPart@r(n); DISCRETE <- discretePart@r(n)
              U*AC+(1-U)*DISCRETE}

    pnew <- function(q, lower.tail = TRUE, log.p = FALSE )
                 {p <- cbind(acPart@p(q, lower.tail = lower.tail),
                         discretePart@p(q, lower.tail = lower.tail))
                  p <- as.vector(p%*%mixCoeff)
                  if(log.p) p <- log(p)
                  return(p)
                 }

    supp <- discretePart@support
    gaps <- .mergegaps(acPart@gaps,discretePart@support)
    #mixDistr[[1]]@gaps <- gaps

    qL1 <- min(getLow(acPart), getLow(discretePart))
    qU1 <- max(getUp(acPart), getUp(discretePart))
    n <- getdistrOption("DefaultNrGridPoints")
    h <- (qU1-qL1)/n
    ep <- getdistrOption("DistrResolution")
    xseq <- unique(sort(c(seq(from = qL1, to = qU1, by = h),gaps-ep,gaps,
                   gaps+ep,supp-ep,supp, supp+ep )))
    px.l <- pnew(xseq, lower.tail = TRUE)
    px.u <- pnew(xseq, lower.tail = FALSE)

    qL2 <- min(acPart@q(0),discretePart@q(0))
    qU2 <- max(acPart@q(1),discretePart@q(1))

    qnew <- .makeQNew(xseq, px.l, px.u, FALSE, qL2, qU2)
    .withSim <- acPart@.withSim | discretePart@.withSim
    .withArith <- acPart@.withArith | discretePart@.withArith
    .lowerExact <- .lowerExact(acPart) & .lowerExact(discretePart)

    Symmetry <- NoSymmetry()
    if(is(acPart@Symmetry,"SphericalSymmetry")&& 
       is(discretePart@Symmetry,"SphericalSymmetry"))
       if(.isEqual(SymmCenter(Symmetry(acPart)),
                   SymmCenter(Symmetry(discretePart))))
       Symmetry <- SphericalSymmetry(SymmCenter(Symmetry(acPart)))   



    new("UnivarLebDecDistribution", p = pnew, r = rnew, d = NULL, q = qnew,
         mixCoeff = mixCoeff, mixDistr = mixDistr, .withSim = .withSim,
         .withArith = .withArith, .lowerExact = .lowerExact, support = supp,
         gaps = gaps, Symmetry = Symmetry)
}

############################## Accessor / Replacement functions

setMethod("discretePart", "UnivarLebDecDistribution",
           function(object) object@mixDistr[[2]])
setReplaceMethod("discretePart", "UnivarLebDecDistribution",
          function(object,value){
          obj <- UnivarLebDecDistribution(acPart = acPart(object),
                 discretePart = value,
                 acWeight = acWeight(object),
                 discreteWeight = discreteWeight(object))
          obj})

setMethod("discreteWeight", "UnivarLebDecDistribution",
           function(object) object@mixCoeff[2])
setReplaceMethod("discreteWeight", "UnivarLebDecDistribution",
          function(object,value){
          obj <- UnivarLebDecDistribution(acPart = acPart(object),
                 discretePart = discretePart(object),
                 acWeight = 1-value,
                 discreteWeight = value)
          obj})

setMethod("acPart", "UnivarLebDecDistribution",
           function(object) object@mixDistr[[1]])
setReplaceMethod("acPart", "UnivarLebDecDistribution",
          function(object,value){
          obj <- UnivarLebDecDistribution(acPart = value,
                 discretePart = discretePart(object),
                 acWeight = acWeight(object),
                 discreteWeight = discreteWeight(object))
          obj})

setMethod("acWeight", "UnivarLebDecDistribution",
           function(object) object@mixCoeff[1])
setReplaceMethod("acWeight", "UnivarLebDecDistribution",
          function(object,value){
          obj <- UnivarLebDecDistribution(acPart = acPart(object),
                 discretePart = discretePart(object),
                 acWeight = value,
                 discreteWeight = 1-value)
          obj})

#setMethod("support", "UnivarLebDecDistribution",
#           function(object) object@mixDistr[[2]]@support)

#setMethod("gaps", "UnivarLebDecDistribution",
#           function(object) object@mixDistr[[1]]@gaps)


setMethod("p.discrete", "UnivarLebDecDistribution",
           function(object, CondOrAbs="cond"){ 
                  CondOrAbs0 <- pmatch(CondOrAbs,c("cond","abs"),nomatch=1)
                  pd <- object@mixDistr[[2]]@p
                  if(CondOrAbs0==1)
                       return(pd)
                  else {wd <- discreteWeight(object)
                        return(function(q, lower.tail = TRUE, log.p = FALSE ){
                               wd * pd(q, lower.tail = lower.tail, log.p = log.p)
                        })
                  }   
           })
setMethod("d.discrete", "UnivarLebDecDistribution",
           function(object, CondOrAbs="cond"){ 
                  CondOrAbs0 <- pmatch(CondOrAbs,c("cond","abs"),nomatch=1)
                  dd <- object@mixDistr[[2]]@d
                  if(CondOrAbs0==1)
                       return(dd)
                  else {wd <- discreteWeight(object)
                        return(function(x, log = FALSE ){
                               wd * dd(x, log = log)
                        })
                  }   
           })
setMethod("q.discrete", "UnivarLebDecDistribution",
           function(object) object@mixDistr[[2]]@q)
setMethod("r.discrete", "UnivarLebDecDistribution",
           function(object) object@mixDistr[[2]]@r)
setMethod("p.ac", "UnivarLebDecDistribution",
           function(object, CondOrAbs="cond"){ 
                  CondOrAbs0 <- pmatch(CondOrAbs,c("cond","abs"),nomatch=1)
                  pd <- object@mixDistr[[1]]@p
                  if(CondOrAbs0==1)
                       return(pd)
                  else {wa <- acWeight(object)
                        return(function(q, lower.tail = TRUE, log.p = FALSE ){
                               wa * pd(q, lower.tail = lower.tail, log.p = log.p)
                        })
                  }   
           })
setMethod("d.ac", "UnivarLebDecDistribution",
           function(object, CondOrAbs="cond"){ 
                  CondOrAbs0 <- pmatch(CondOrAbs,c("cond","abs"),nomatch=1)
                  dd <- object@mixDistr[[1]]@d
                  if(CondOrAbs0==1)
                       return(dd)
                  else {wa <- acWeight(object)
                        return(function(x, log = FALSE ){
                               wa * dd(x, log = log)
                        })
                  }   
           })
setMethod("q.ac", "UnivarLebDecDistribution",
           function(object) object@mixDistr[[1]]@q)
setMethod("r.ac", "UnivarLebDecDistribution",
           function(object) object@mixDistr[[1]]@r)



setMethod("p.l", "UnivarLebDecDistribution", function(object){
           ep <- getdistrOption("TruncQuantile")
           w.d <- discreteWeight(object)
           w.a <- acWeight(object)   
           if(w.d<ep) return(p(object))
           mixCoeff <- c(w.a,w.d)
           p.a <- p(acPart(object))
           p.d <- p.l(discretePart(object))
           return(function(q, lower.tail = TRUE, log.p = FALSE){
                  p <- cbind(p.a(q, lower.tail = lower.tail),
                             p.d(q, lower.tail = lower.tail))
                  p <- as.vector(p%*%mixCoeff)
                  if(log.p) p <- log(p)
                  return(p)
                 })
       })

### right continuous quantile function

setMethod("q.r", "UnivarLebDecDistribution", function(object){
    ep <- getdistrOption("TruncQuantile")
    if(discreteWeight(object)<ep) return(q(object))
    supp <- support(object)
    gaps <- gaps(object)
    aP <- acPart(object)
    dP <- discretePart(object)
    pl <- p.l(object)
    qL1 <- min(getLow(aP), getLow(dP))
    qU1 <- max(getUp(aP), getUp(dP))
    n <- getdistrOption("DefaultNrGridPoints")
    h <- (qU1-qL1)/n
    xseq <- unique(sort(c(seq(from = qL1, to = qU1, by = h),
                   gaps-ep,gaps,gaps+ep,
                   supp-ep,supp, supp+ep )))
    px.l <- pl(q=xseq, lower.tail = TRUE)
    px.u <- pl(q=xseq, lower.tail = FALSE)

    qL2 <- min(aP@q(0),dP@q(0))
    qU2 <- max(aP@q(1),dP@q(1))

    return( .makeQNew(xseq, px.l, px.u, FALSE, qL2, qU2))
})






setMethod("prob", "UnivarLebDecDistribution", 
function(object) {pr0 <- prob(as(object@mixDistr[[2]],"DiscreteDistribution"))
                  d <- discreteWeight(object)
                  return(rbind("cond"=pr0,"abs"=d*pr0))})
                  


############################## setAs relations

setAs("AbscontDistribution", "UnivarLebDecDistribution",
       function(from) UnivarLebDecDistribution(acPart = from))

setAs("DiscreteDistribution", "UnivarLebDecDistribution",
       function(from) UnivarLebDecDistribution(discretePart = from))



############################## Arithmetics

setMethod("+", c("UnivarLebDecDistribution","UnivarLebDecDistribution"),
function(e1,e2){
         ep <- getdistrOption("TruncQuantile")

         simple <- 0

         if(discreteWeight(e1)>1-ep)
            {simple <- 1
             e10 <- discretePart(e1)
             if(is(e10,"Dirac"))
               return(e2+location(as(e10,"Dirac")))}
         if(discreteWeight(e1) < ep)
            {simple <- 1
             e10 <- acPart(e1)}

         if(discreteWeight(e2)>1-ep)
            {simple <- simple + 1
             e20 <- discretePart(e2)
             if(is(e20,"Dirac"))
               return(e1+location(as(e20,"Dirac")))}
         if(discreteWeight(e2) < ep)
            {simple <- simple + 1
             e20 <- acPart(e2)}

         if(simple == 2)
            return(as(e10+e20,"UnivarLebDecDistribution"))

         p1.a <- acPart(e1)
         p2.a <- acPart(e2)
         p1.d <- discretePart(e1)
         p2.d <- discretePart(e2)
         w1.a <- acWeight(e1)
         w2.a <- acWeight(e2)
         w1.d <- discreteWeight(e1)
         w2.d <- discreteWeight(e2)


         f1 <- if(w1.a*w2.a>ep)  p1.a + p2.a else Norm()
         f2 <- if(w1.a*w2.d>ep)  p1.a + p2.d else Norm()
         if (identical(e1,e2)) f3 <- f2
         else {f3 <- if(w1.d*w2.a>ep)  p1.d + p2.a else Norm()}
         f.d <- if(w1.d*w2.d>ep) p1.d + p2.d else Dirac(0)


         w.d <- w2.d * w1.d
         w.a <- 1 - w.d

         rnew <- function(n) {
                 U <- sample(1:3, size = n, replace = TRUE, prob =
                             c(w1.a*w2.a, w1.a*w2.d, w1.d*w2.a))
                 R1 <- r(f1)(n)
                 R2 <- r(f2)(n)
                 R3 <- r(f3)(n)
                 (U==1)*R1 + (U==2)*R2 + (U==3)*R3
                 }
         pnew <- function(q, lower.tail = TRUE, log.p = FALSE)
                     {p1 <- w1.a * w2.a *
                         p(f1)(q, lower.tail = lower.tail, log.p = log.p)
                      p2 <- w1.a * w2.d *
                         p(f2)(q, lower.tail = lower.tail, log.p = log.p)
                      p3 <- w1.d * w2.a *
                         p(f3)(q, lower.tail = lower.tail, log.p = log.p)
                      (p1 + p2 + p3)/w.a
                      }
         dnew <- function(x, log = FALSE)
                     {(w1.a * w2.a * d(f1)(x, log = log)+
                       w1.a * w2.d * d(f2)(x, log = log)+
                       w1.d * w2.a * d(f3)(x, log = log)) / w.a }

         qL2 <- min(f1@q(0),f2@q(0),f3@q(0))
         qU2 <- max(f1@q(1),f2@q(1),f3@q(1))

         qL <-  min(getLow(f1), getLow(f2), getLow(f3))
         qU <-  max(getUp(f1), getUp(f2), getUp(f3))
         n <- getdistrOption("DefaultNrGridPoints")
         h <- (qU-qL)/n
         xseq <- seq(from = qL, to = qU, by = h)
         px.l <- pnew(xseq, lower.tail = TRUE)
         px.u <- pnew(xseq, lower.tail = FALSE)

         qnew <- .makeQNew(xseq, px.l, px.u, FALSE, qL2, qU2)

         f.a <- AbscontDistribution(
                    r = rnew, p = pnew, q = qnew, d = dnew, .withSim = FALSE,
                    .withArith = TRUE)

         obj <- UnivarLebDecDistribution(discretePart = f.d, acPart = f.a,
                                  discreteWeight = w.d)
         if(getdistrOption("simplifyD"))
            obj <- simplifyD(obj)

         return(obj)
         })

setMethod("*", c("UnivarLebDecDistribution","numeric"),
          function(e1, e2) {
          if (length(e2)>1) stop("length of operator must be 1")

          if (isTRUE(all.equal(e2,1))) return(e1)
          if (isTRUE(all.equal(e2,0)))
               return(new("Dirac", location = 0))

          Distr <- UnivarLebDecDistribution(
                     discretePart = discretePart(e1)*e2,
                     acPart = acPart(e1)*e2,
                     discreteWeight = discreteWeight(e1),
                     acWeight = acWeight(e1))

          object <- new("AffLinUnivarLebDecDistribution",
                    r = Distr@r, d = Distr@d, p = Distr@p,
                    q = Distr@q, X0 = e1, mixDistr = Distr@mixDistr,
                    mixCoeff = Distr@mixCoeff,
                    a = e2, b = 0, .withSim  = e1@.withSim,
                    .logExact = .logExact(e1), .lowerExact = .lowerExact(e1),
                    .withArith = TRUE, support= support(Distr),
                    gaps = gaps(Distr))
          object})

setMethod("+", c("UnivarLebDecDistribution","numeric"),
          function(e1, e2) {
          if (length(e2)>1) stop("length of operator must be 1")
          if (isTRUE(all.equal(e2,0))) return(e1)

          Distr <- UnivarLebDecDistribution(
                     discretePart = discretePart(e1)+e2,
                     acPart = acPart(e1)+e2,
                     discreteWeight = discreteWeight(e1),
                     acWeight = acWeight(e1))

          object <- new("AffLinUnivarLebDecDistribution",
                    r = Distr@r, d = Distr@d, p = Distr@p,
                    q = Distr@q, X0 = e1, mixDistr = Distr@mixDistr,
                    mixCoeff = Distr@mixCoeff,
                    a = 1, b = e2, .withSim  = e1@.withSim,
                    .withArith = TRUE, support= support(Distr),
                    .logExact = .logExact(e1), .lowerExact = .lowerExact(e1),
                    gaps = gaps(Distr))
          object})

#setMethod("*", c("numeric","UnivarLebDecDistribution"),
#          function(e1, e2) e2*e1)
#setMethod("+", c("numeric","UnivarLebDecDistribution"),
#          function(e1, e2) e2+e1)


## Group Math for absolutly continuous distributions
setMethod("Math", "UnivarLebDecDistribution",
          function(x){
            rnew <- function(n, ...){}
            body(rnew) <- substitute({ f(g(n, ...)) },
                              list(f = as.name(.Generic), g = x@r))
            object <- UnivarLebDecDistribution(r = rnew)

            if(getdistrOption("simplifyD"))
               object <- simplifyD(object)
            object
          })

## ??? setGroupGeneric("exFct", def=function(object) standardGeneric("exFct"), group = list(exp, abs,...))

## exact: abs for absolutly continuous distributions
setMethod("abs", "UnivarLebDecDistribution",
          function(x)UnivarLebDecDistribution(
                     discretePart = abs(discretePart(x)),
                     acPart = abs(acPart(x)),
                     discreteWeight = discreteWeight(x),
                     acWeight = acWeight(x)))

## exact: exp for absolutly continuous distributions
setMethod("exp", "UnivarLebDecDistribution",
          function(x){
           if(acWeight(x)>1-getdistrOption("TruncQuantile"))
              {x <- force(acPart(x))
               return(as(exp(x),"UnivarLebDecDistribution"))}
           if(acWeight(x) < getdistrOption("TruncQuantile"))
              {x <- force(discretePart(x))
               return(as(exp(x),"UnivarLebDecDistribution"))}
           return(UnivarLebDecDistribution(
                     discretePart = exp(discretePart(x)),
                     acPart = exp(acPart(x)),
                     discreteWeight = discreteWeight(x),
                     acWeight = acWeight(x)))})


setMethod("log", "UnivarLebDecDistribution",
          function(x, base = exp(1)){
           xs <- as.character(deparse(match.call(
                 call = sys.call(sys.parent(1)))$x))
           ep <- getdistrOption("TruncQuantile")
           basl <- log(base)
           if(p(x)(0)>ep) 
                stop(gettextf("log(%s) is not well-defined with positive probability ", xs))
           
           if(acWeight(x)>1-getdistrOption("TruncQuantile")) 
              {x <- force(acPart(x))
               return(as(log(x, base = base),"UnivarLebDecDistribution"))}    
           if(acWeight(x) < getdistrOption("TruncQuantile")) 
              {x <- force(discretePart(x))
               return(as(log(x, base = base),"UnivarLebDecDistribution"))}    
           return(UnivarLebDecDistribution(
                     discretePart = log(discretePart(x), base = base),
                     acPart = log(acPart(x), base = base),
                     discreteWeight = discreteWeight(x),
                     acWeight = acWeight(x)))})

setMethod("log10", "UnivarLebDecDistribution",
          function(x) log(x = x, base = 10))

setMethod("sign", "UnivarLebDecDistribution",
          function(x){ 
          d0 <- d.discrete(as(x,UnivarLebDecDistribution))(0)
          DiscreteDistribution(supp=c(-1,0,1), 
              prob=c(p(x)(-getdistrOption("TruncQuantile")),
                     d0,
                     p(x)(getdistrOption("TruncQuantile"), lower=FALSE)))                     
          })

setMethod("sqrt", "UnivarLebDecDistribution",
            function(x) x^0.5)
#######################

