"cdfkmu" <-
function(x, para, paracheck=TRUE, getmed=TRUE, qualo=NA, quahi=NA,
         marcumQ=TRUE, marcumQmethod=c("chisq", "delta", "integral")) {

   marcumQmethod <- match.arg(marcumQmethod)
   if(paracheck == TRUE & ! are.parkmu.valid(para)) return()

   names(para$para) <- NULL
   K  <- para$para[1]
   MU <- para$para[2]

   MEDIAN <- NA
   if(getmed) MEDIAN <- quakmu(0.5, para, getmed=FALSE, qualo=qualo, quahi=quahi,
                                 marcumQ=marcumQ, marcumQmethod=marcumQmethod)
   was.getmed.true <- getmed
   getmed <- FALSE

   diracdelta <- 0
   fixedM <- FALSE
   if(! is.finite(K)) { M <- MU; fixedM <- TRUE }
   if(fixedM) {
      tmpA <- 4*M/exp(2*M)
      tmpB <- sqrt(2*M*pi)/exp(M)
      toI <- 4*M*0
      B1 <- besselI(toI, nu=1)
      if(! is.finite(B1)) B1 <- 0
      B2 <- besselI(M, nu=1/2)
      if(! is.finite(B2)) B2 <- 0
      diracdelta <- tmpA*B1 + (1 - tmpB*B2)
      names(diracdelta) <- "Dirac Delta x=0"
      # This resetting of the Dirac Delta functions part of the parameter list is made
      # so that the conditional tests numerical equivalence is effectively bypassed.
      # This feature is provided so that should a user create their own parameter list manually
      # and not through the function vec2par that it will be assumed that the Dirac computed here is fine.
      if(length(para$diracdelta) == 0 | is.na(para$diracdelta)) para$diracdelta <- diracdelta
      if(diracdelta != para$diracdelta) {
          warning("Dirac delta (x=0) computed herein does not match that embedded in the parameter object, ",
                  "going to use the freshly computed one")
          warning("Dirac delta = ", diracdelta,"   and embedded = ", para$diracdelta)
      }
      #message("Note: The Dirac Delta function for (x=0) for this parameterized ",
      #        "Kappa-Mu distribution provides ",round(diracdelta, digits=6)," of total probability.\n")
   }

    LARGE <- sqrt(.Machine$double.xmax)
    SMALL <-  .Machine$double.eps
    "marcumq.integral" <- function(a, b, nu=NULL) {
       if(is.null(nu)) {
          warning("nu is NULL for Marcum Q function, returning NA")
          return(NA)
       }
       if(a == 0) a <- SMALL
       "afunc" <- function(t) {
           B <- vector(mode="numeric", length=length(t))
           B <- sapply(1:length(B), function(i) {
                  toI <- a * t[i]
                  b <- NULL
                  try(b <- besselI(toI, nu=nu-1,
                                   expon.scaled=TRUE)/exp(-toI),
                                   silent=TRUE)
                  if(is.null(b) | is.nan(b) | ! is.finite(b)) return(LARGE)
                  return(b) })
           z <- t^nu * exp(-(t^2 + a^2)/2) * B
           return( z )
       }
       int1 <- NULL
       try( int1 <- integrate(afunc, lower=b, upper=Inf) )
       if(is.null(int1)) return(NA)
       return( int1$value/a^(nu-1) )
    }

    # Shi, Q., Karasawa, Y., 2012, An intuitive methodology
    # for efficient evaluation of the Nuttall Q-function and
    # performance analysis of energy detection in fading
    # channels: IEEE Wireless Communications Letters,
    # v. 1, no. 2, pp. 109--112.
    "marcumq.bydelta" <- function(a, b, nu=NULL) {
       if(is.null(nu)) {
          warning("nu is NULL for Marcum Q function, returning NA")
          return(NA)
       }
       delta <- nu %% as.integer(nu); nuint <- as.integer(nu)
       if(is.nan(delta)) { delta <- nu; nuint <- 0 }
       if(nu >  0) { beg <- 0;     end <- nuint - 1; sign <-  1 }
       if(nu <= 0) { beg <- nuint; end <-       - 1; sign <- -1 }
       if(nuint == 0) { sign <- 0 }
       if(a == 0) a <- SMALL
       Qdelta <- marcumq.integral(a, b, nu=delta)
       tmp <- ifelse(sign == 0, 0,
                sum(sapply(beg:end, function(i) { (b/a)^(i+delta) *
                        besselI(a*b, nu=i+delta) })))
       Qnuint <- sign * exp(-(a^2+b^2)/2) * tmp
       if(is.nan(Qnuint)) Qnuint <- 0
       Qnu <- Qdelta + Qnuint
       #message("nuint=", nuint, "  delta=",delta,"  beg=",beg,"  end=",end,"   sign=",sign,"\n")
       #message("Qdelta=",Qdelta,"  Qnuint=",Qnuint,"\n")
       if(is.nan(Qnu)) Qnu <- 0
       return(Qnu)
    }

    if(marcumQmethod == "chisq") {
       marcumq <- function(a, b, nu=1) {
                 pchisq(b^2, df=2*nu, ncp=a^2, lower.tail=FALSE) }
    } else if(marcumQmethod == "delta") {
       marcumq <- marcumq.bydelta
    } else {
       marcumq <- marcumq.integral
    }

    f <- sapply(1:length(x), function(i) {
                   xi <- x[i]
                   if(xi <  0)                     return(NA)
                   if(xi == 0)                     return(diracdelta)
                   if(! is.na(qualo) & xi < qualo) return(diracdelta)
                   if(! is.na(quahi) & xi > quahi) return(1)
                   if(! is.finite(xi))             return(1)
                   if(is.na(xi))                   return(NA)
                   if(marcumQ & K != 0 & is.finite(K)) {
                      Q <- marcumq(sqrt(2*K*MU), sqrt(2*(1+K)*MU)*xi, nu=MU)
                      #message(" marcumQ=",Q)
                      return(1 - Q)
                   }
                   int1 <- NULL
                   try( int1 <- integrate(pdfkmu, SMALL, xi,
                                para=para, paracheck=FALSE), silent=TRUE )
                   ifelse(is.null(int1), return(NA),
                          return(int1$value + diracdelta)) })
    names(f) <- NULL
    f[! is.finite(f)] <- NA
    if(was.getmed.true) {
       #print(MEDIAN); print(f)
       f[! is.na(MEDIAN) & x < MEDIAN & f > 0.5] <- 0
       f[! is.na(MEDIAN) & x < MEDIAN & f < 0  ] <- 0
    }
    return(f)
}


