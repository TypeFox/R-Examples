"EMPIRcop" <-
function(u,v, para=NULL, ctype=c("weibull", "hazen", "1/n", "bernstein"),
                                           bernprogress=FALSE, ...) {
  bernstein <- FALSE
  ctype <- match.arg(ctype)
  uobs <- vobs <- NA
  if(exists("para", para)) {
     uobs <- para$para[,1]
     vobs <- para$para[,2]
     if(exists("bernstein",    para)) bernstein    <- para$bernstein
     if(exists("bernprogress", para)) bernprogress <- para$bernprogress
     para <- para$para # now reset the para to ONLY have the uobs and vobs
     # because if bernstein is triggered the a secondary call to this function
     # is made and bernstein must be turned off.
  }
  if(length(names(para)) != 2) {
     warning("a data.frame having only two columns is required in the para ",
             "argument or para$para")
     return(NULL)
  }

  if(ctype == "bernstein") bernstein <- TRUE

  uobs <- para[,1]
  vobs <- para[,2]
  n <- length(uobs)

  nu <- length(u); nv <- length(v)
  if(nu > 1 & nv > 1 & nu != nv) {
     warning("length u = ",nu, " and length v = ",nv)
     warning("longer object length is not a multiple of shorter object length, ",
             "no recycling in EMPIRcop()")
     return(NA)
  }
  # The extra hassle of vectorization made here is to handle situations
  # in which nested integrals are used where uneven vectors can be passed.
  if(nu == 1) {
     u <- rep(u, nv)
  } else if(nv == 1) {
     v <- rep(v, nu)
  }
  nu <- length(u) # make sure to RESET IT!

  if(bernstein) { # will potentially burn SERIOUS CPU time
     if(bernprogress) {
         message("Bernstein triggered, index of each {u,v} shown: ",
                 appendLF=FALSE)
     }
     ber <- sapply(1:nu, function(k) {
                   if(bernprogress) message(k,"-", appendLF=FALSE);
                   tmpA <- sapply(1:n, function(i) {
                        A <- choose(n, i) * u[k]^i * (1-u[k])^(n-i)
                        tmpB <- sapply(1:n, function(j) {
                                  B <- choose(n, j) * v[k]^j *(1-v[k])^(n-j)
                         return(EMPIRcop(i/n, j/n, para=para, ctype="1/n")*A*B) })
                                  return(sum(tmpB)) })
                        return(sum(tmpA)) })
     if(bernprogress) message("\n", appendLF=FALSE)
     return(ber)
  } else {
     R <- rank(uobs);     S <- rank(vobs)
     if(ctype == "hazen") {
        R <- (R - 0.5)/n; S <- (S - 0.5)/n
     } else if(ctype == "weibull") {
        R <- R/(n+1);     S <- S/(n+1)
     } else {
        R <- R/n;         S <- S/n
     }
     sapply(1:nu, function(i) sum(as.numeric(R <= u[i] & S <= v[i]))/n )
  }
}

