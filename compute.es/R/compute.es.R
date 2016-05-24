##=========================== compute.es ==============================##

# to do:
# 1. add dependent (paired) group options to each
# 2. help function for easy recall (something with args())
# 3. create gui for this (and web application!)


# ## UPDATE:  SINGLE UNIFYING FUNCTION--IN PROCESS
# es <- function(...){
#   mf <- match.call(expand.dots = TRUE)
#       args <- match(c("d", "n.1", "n.2", 
#                       "r", "var.r", "n",  
#                       "m.1", "m.2", "sd.1", "sd.2", "n.1", "n.2",  
#                       "s.pooled",  
#                       "p", "tail",  
#                       "t",  
#                       "lor", "var.lor",  
#                       "chi.sq",  
#                       "f", 
#                       "B", "D", "n.0", 
#                       "p1", "p2", "n.ab", "n.cd",   
#                       "R", "q",  
#                       "m.1.adj", "m.2.adj", "sd.adj",
#                       "level", "dig", "id","data"), names(mf), 
#                     0)
#     mf <- mf[c(1, args)]
#     mf$drop.unused.levels <- TRUE
# #   if(!is.null(mf$d)) {    
# #     mf.d <- mf[[match("d", names(mf))]]
# #     es <- eval(mf.d, data, enclos = sys.frame(sys.parent()))
# #     mf.n.1 <- mf[[match("n.1", names(mf))]]
# #     n.1 <- eval(mf.n.1, data, enclos = sys.frame(sys.parent()))
# #     mf.n.2 <- mf[[match("n.2", names(mf))]]
# #     n.2 <- eval(mf.n.2, data, enclos = sys.frame(sys.parent()))
# #     mf.id <- mf[[match("id", names(mf))]]
# #     id <- eval(mf.id, data, enclos = sys.frame(sys.parent()))
# #     
# #   }
# #   if(!is.null(mf$r)) {
# #     mf.r <- mf[[match("r", names(mf))]]
# #     es <- eval(mf.r, data, enclos = sys.frame(sys.parent()))
# #     mf.var.r <- mf[[match("var.r", names(mf))]]
# #     var.r <- eval(mf.var.r, data, enclos = sys.frame(sys.parent()))
# #     mf.n <- mf[[match("n", names(mf))]]
# #     n <- eval(mf.n, data, enclos = sys.frame(sys.parent()))
# #     mf.id <- mf[[match("id", names(mf))]]
# #     id <- eval(mf.id, data, enclos = sys.frame(sys.parent()))
# #    
# #  }
#   
#   out <- mf #data.frame(mf) #data.frame(id, es)
#   return(out)
#   #out <- data.frame(d,n.1,n.2) #mf
#   #return(out)
#   #dat <- list(...)             
#   #return(dat)
#   
# }
# 
# 
# es <- function(d, n.1, n.2, # DES(d, n.1, n.2)
#    r, var.r = NULL, n,  # RES(r, var.r = NULL, n)
#    m.1, m.2, sd.1, sd.2,  # MES(m.1, m.2, sd.1, sd.2, n.1, n.2)
#    s.pooled,  # MES2(m.1, m.2, s.pooled, n.1, n.2,)
#    p, tail = "two",  # PES(p, n.1, n.2, tail = "two")
#    t,  # TES(t, n.1, n.2)
#    lor, var.lor,  # LORES(lor, var.lor, n.1, n.2)
#    chi.sq,  # CHIES(chi.sq, n)
#    f, # FES(f, n.1, n.2)
#    B, D, n.0, # FAILES(B, D, n.1, n.0)
#    p1, p2, n.ab, n.cd,   # PROPES(p1, p2, n.ab, n.cd)
#    R, q,   # A.FES(f, n.1, n.2, R, q)
#    m.1.adj, m.2.adj, sd.adj,  # A.MES(m.1.adj, m.2.adj, sd.adj); A.MES2(m.1.adj, m.2.adj, s.pooled, n.1, n.2, R, q, ...)
#      # A.PES(p, n.1, n.2, R, q, )
#      # A.TES(t, n.1, n.2, R, q, )
#    level = 95, cer = 0.2, dig = 2, id = NULL, data = NULL, ...) {
#       #if (!is.null(data)) {
#         mf <- match.call()
#         temp <- as.list(mf)
#         
#         args <- match(c("d", "n.1", "n.2", 
#                     "r", "var.r", "n",  
#                     "m.1", "m.2", "sd.1", "sd.2", "n.1", "n.2",  
#                     "s.pooled",  
#                     "p", "tail",  
#                     "t",  
#                     "lor", "var.lor",  
#                     "chi.sq",  
#                     "f", 
#                     "B", "D", "n.0", 
#                     "p1", "p2", "n.ab", "n.cd",   
#                     "R", "q",  
#                     "m.1.adj", "m.2.adj", "sd.adj",
#                     "level", "dig", "id","data"), names(mf), 
#                   0)
#         mf <- mf[c(1, args)]
#         mf$drop.unused.levels <- TRUE
#         mf.d <- mf[[match("d", names(mf))]]
#         d <- eval(mf.d, data, enclos = sys.frame(sys.parent()))
#         mf.n.1 <- mf[[match("n.1", names(mf))]]
#         n.1 <- eval(mf.n.1, data, enclos = sys.frame(sys.parent()))
#         mf.n.2 <- mf[[match("n.2", names(mf))]]
#         n.2 <- eval(mf.n.2, data, enclos = sys.frame(sys.parent()))
#         mf.r <- mf[[match("r", names(mf))]]
#         r <- eval(mf.r, data, enclos = sys.frame(sys.parent()))
#         mf.var.r <- mf[[match("var.r", names(mf))]]
#         var.r <- eval(mf.var.r, data, enclos = sys.frame(sys.parent()))
#         mf.n <- mf[[match("n", names(mf))]]
#         n <- eval(mf.n, data, enclos = sys.frame(sys.parent()))
#         mf.id <- mf[[match("id", names(mf))]]
#         id <- eval(mf.id, data, enclos = sys.frame(sys.parent()))
#      # }
#      #out <-  temp  #mf #data.frame(mf) #data.frame(id, es)
#      #return(out)
# #   if(!missing(d) && !missing(n.1) && !missing(n.2)) { 
# #     des(d, n.1, n.2, ...)
# #   }
# #   if(!missing(r) && !missing(n))  { 
# #     res(r, n, ...)
# #     }
# #     else{
# #     stop("Missing argument(s) to calculate effect sizes")
# #     }
# #   #if(!missing(r) && !missing(n)) { res(r, n, ... )}
# #   #}
# 
#   if(!is.null(temp$d) && !is.null(temp$n.1) && !is.null(temp$n.2)) { 
#     des(d, n.1, n.2, ...)
#   }
#   if(!is.null(temp$r) && !is.null(temp$n))  { 
#     res(r, n, ...)
#     }
#     else{
#     stop("Missing argument(s) to calculate effect sizes")
#     }
#   #if(!missing(r) && !missing(n)) { res(r, n, ... )}
#   #}
#   
# }
# 
# 
# es(d=.3, n.1=30, n.2=30)
# es(d=d, n.1=nT, n.2=nC, id=id, data=dat)
# es(r=.3, n=30)
# es(r=r, n=n, id=id, data=dar)
# dar <- data.frame(id=1:30,
#                   r=rnorm(30, 0.5, 0.1),
#                   n=round(rnorm(30, 30, 5), 0))


# CONVENIENCE FUNCTIONS

ztor <- function(z) (exp(2 * z) - 1)/(1 + exp(2 * z))

# d to es

# ## in progress
# des2 <- function(d, n.1, n.2=NULL, type = "indep", r = NULL) {
#    mf <- match.call()
#   if(type == "indep") {
#     var.d<-(n.1+n.2)/(n.1*n.2)+ (d^2)/(2*(n.1+n.2))
#     df<- (n.1+n.2)-2
#     j<-1-(3/(4*df-1))
#     g<-j*d
#     var.g<-j^2*var.d
#     a <- ((n.1 + n.2)^2)/(n.1*n.2)  # will correct for inbalanced n, if applicable
#     n <- n.1 + n.2
#     r <- d/sqrt((d^2) + a) # to compute r from d
#     var.r <- (a^2*var.d)/(d^2 + a)^3
#     lor <- d*(pi/sqrt(3))
#     var.lor <- var.d*(pi^2/3)
#     z <-  0.5*log((1 + r)/(1-r))  
#     var.z <- 1/(n-3) 
#     #z.score <- r/sqrt(var.r)
#   } 
#   if(type == "paired") {
#     var.d<- ((1/n.1) + ((d^2)/(2*n.1)))*2*(1-r) # where r = cor of pre and post test
#     df<- (n.1)-1
#     j<-1-(3/(4*df-1))
#     g<-j*d
#     var.g<-j^2*var.d
#     #a <- ((n.1 )^2)/(n.1*n.2)  # will correct for inbalanced n, if applicable
#     n <- n.1 + n.2
#     r <- d/sqrt((d^2) + a) # to compute r from d
#     var.r <- (a^2*var.d)/(d^2 + a)^3
#     lor <- d*(pi/sqrt(3))
#     var.lor <- var.d*(pi^2/3)
#     z <-  0.5*log((1 + r)/(1-r))  
#     var.z <- 1/(n-3) 
#     #z.score <- r/sqrt(var.r)
#   } 
#   out<-list(MeanDifference= c(d = d,var.d = var.d, g= g, var.g = var.g), 
#            Correlation= c(r= r, var.r = var.r),
#            Log_Odds = c(log_odds = lor, var.log_odds = var.lor),
#            Fishers_z = c(z = z, var.z = var.z),
#            TotalSample = c(n= n))
#   return(out)
# }



des <- function(d, n.1, n.2, level=95, cer=.2,  dig=2, verbose=TRUE, id=NULL, data=NULL) {
  if (!is.null(data)) {
    mf <- match.call()
    args <- match(c("d", "n.1", "n.2", "level", "dig", "id","data"), names(mf), 
                  0)
    mf <- mf[c(1, args)]
    mf$drop.unused.levels <- TRUE
    mf.d <- mf[[match("d", names(mf))]]
    d <- eval(mf.d, data, enclos = sys.frame(sys.parent()))
    mf.n.1 <- mf[[match("n.1", names(mf))]]
    n.1 <- eval(mf.n.1, data, enclos = sys.frame(sys.parent()))
    mf.n.2 <- mf[[match("n.2", names(mf))]]
    n.2 <- eval(mf.n.2, data, enclos = sys.frame(sys.parent()))
    mf.id <- mf[[match("id", names(mf))]]
    id <- eval(mf.id, data, enclos = sys.frame(sys.parent()))
  }
  var.d<-(n.1+n.2)/(n.1*n.2)+ (d^2)/(2*(n.1+n.2))  # POOLED SD
  df<- (n.1+n.2)-2
  j<-1-(3/(4*df-1))
  g<-j*d
  var.g<-j^2*var.d
  
  a <- ((n.1 + n.2)^2)/(n.1*n.2)  # will correct for inbalanced n, if applicable
  n <- n.1 + n.2
  r <- d/sqrt((d^2) + a) # to compute r from d
  var.r <- (a^2*var.d)/(d^2 + a)^3
  
  lor <- d*(pi/sqrt(3))
  var.lor <- var.d*(pi^2/3)
  
  z <-  0.5*log((1 + r)/(1-r))  
  var.z <- 1/(n-3)
  
  # UPDATE: ADDED FOLLOWING TO EACH FUNCTION (FOR P-VAL, CI, U3, CLES, CLIFF'S DELTA, NNT)
  alpha <- (100 - level)/100
  crit <- qt(alpha/2, df, lower.tail = FALSE)
  #crit <- qnorm(alpha)
  zval.d <- d/sqrt(var.d)
  pval.d <- 2 * pt(abs(zval.d), df , lower.tail = FALSE)
  lower.d <- d - crit * sqrt(var.d)
  upper.d <- d + crit * sqrt(var.d)
  U3.d <- pnorm(d)*100
  cl.d<- (pnorm((d)/sqrt(2)))*100
  cliffs.d <- 2 * pnorm(d/sqrt(2)) - 1
  zval.g <- g/sqrt(var.g)
  pval.g <- 2 * pt(abs(zval.g), df , lower.tail = FALSE)
  lower.g <- g - crit * sqrt(var.g)
  upper.g <- g + crit * sqrt(var.g)
  U3.g <- pnorm(g)*100
  cl.g <- (pnorm((g)/sqrt(2)))*100
  
  zval.z <- z/sqrt(var.z)
  pval.z <- 2 * pt(abs(zval.z), df , lower.tail = FALSE)
  lower.z <- z - crit * sqrt(var.z)
  upper.z <- z + crit * sqrt(var.z)
  
  pval.r <- pval.z
  lower.r <- ztor(lower.z)
  upper.r <- ztor(upper.z)
  
  zval.lor <- lor/sqrt(var.lor)
  pval.lor <- 2 * pt(abs(zval.lor), df , lower.tail = FALSE)
  lower.lor <- lor - crit * sqrt(var.lor)
  upper.lor <- lor + crit * sqrt(var.lor)
  
  nnt <- 1/(pnorm(d- qnorm(1-cer))-cer)
  
  if (!is.null(data)) {
#     cat("\n")
#     message("    EFFECT SIZE CALCULATION (FOR VECTOR INPUT)")
#     cat("\n")
    
    out <- round(data.frame(id=id, N.total=n, n.1=n.1, n.2=n.2, d=d, var.d=var.d, l.d=lower.d, u.d=upper.d, U3.d=U3.d,
                            cl.d=cl.d, cliffs.d=cliffs.d, pval.d=pval.d, g=g, var.g=var.g,l.g=lower.g, 
                            u.g=upper.g, U3.g=U3.g, cl.g=cl.g, pval.g=pval.g, r=r, var.r=var.r, l.r=lower.r,
                            u.r=upper.r, pval.r=pval.r, fisher.z=z, var.z=var.z, l.z=lower.z, u.z=upper.z,
                            OR=exp(lor), l.or=exp(lower.lor), u.or=exp(upper.lor), pval.or=pval.lor, 
                            lOR=lor, l.lor=lower.lor, u.lor=upper.lor, pval.lor=pval.lor,
                            NNT=nnt), dig)
    return(out)
    
  }
  else{
    if (verbose) {
#     cat("\n")
#     message("    EFFECT SIZE CALCULATION (FOR SINGLE INPUT)")
#     cat("\n")
    
    cat(#"\n", "\n", "   EFFECT SIZE CALCULATION (FOR SINGLE INPUT)", "\n","\n",
      "Mean Differences ES:","\n", "\n", 
      "d [",level,"%CI] =", round(d, dig),
      "[",round(lower.d, dig),",",round(upper.d,dig),"]","\n",
      "\ var(d) =", round(var.d,dig), "\n",
      "\ p-value(d) =", round(pval.d,dig), "\n",
      "\ U3(d) =", round(U3.d, dig),"%", "\n", 
      "\ CLES(d) =", round(cl.d, dig),"%","\n",
      "\ Cliff's Delta =", round(cliffs.d, dig),"\n", "\n",
      "g [",level,"%CI] =", round(g, dig),
      "[",round(lower.g, dig),",",round(upper.g,dig),"]","\n",
      "\ var(g) =", round(var.g,dig), "\n",
      "\ p-value(g) =", round(pval.g,dig), "\n",
      "\ U3(g) =", round(U3.g, dig),"%", "\n", 
      "\ CLES(g) =", round(cl.g, dig),"%","\n", "\n",
      
      "Correlation ES:","\n", "\n", 
      "r [",level,"%CI] =", round(r, dig),
      "[",round(lower.r, dig),",",round(upper.r,dig),"]","\n",
      "\ var(r) =", round(var.r,dig), "\n",
      "\ p-value(r) =", round(pval.r,dig), "\n","\n",
      "z [",level,"%CI] =", round(z, dig),
      "[",round(lower.z, dig),",",round(upper.z,dig),"]","\n", 
      "\ var(z) =", round(var.z,dig), "\n",
      "\ p-value(z) =", round(pval.z,dig), "\n", "\n",
      
      "Odds Ratio ES:","\n", "\n", 
      "OR [",level,"%CI] =", round(exp(lor), dig),
      "[",round(exp(lower.lor), dig),",",round(exp(upper.lor),dig),"]","\n", 
      "\ p-value(OR) =", round(pval.lor,dig), "\n", "\n",
      "Log OR [",level,"%CI] =", round(lor, dig),
      "[",round(lower.lor, dig),",",round(upper.lor,dig),"]","\n",
      "\ var(lOR) =", round(var.lor,dig), "\n",
      "\ p-value(Log OR) =", round(pval.lor,dig), "\n", "\n", 
      
      "Other:","\n", "\n", 
      "NNT =", round(nnt, dig),"\n", 
      "Total N =", n
    )
    }
    out <- round(data.frame(N.total = n, n.1 = n.1, 
                            n.2 = n.2, d = d, var.d = var.d, l.d = lower.d, u.d = upper.d, 
                            U3.d = U3.d, cl.d = cl.d, cliffs.d = cliffs.d, pval.d = pval.d, 
                            g = g, var.g = var.g, l.g = lower.g, u.g = upper.g, 
                            U3.g = U3.g, cl.g = cl.g, pval.g = pval.g, r = r, 
                            var.r = var.r, l.r = lower.r, u.r = upper.r, pval.r = pval.r, 
                            fisher.z = z, var.z = var.z, l.z = lower.z, u.z = upper.z, 
                            OR = exp(lor), l.or = exp(lower.lor), u.or = exp(upper.lor), 
                            pval.or = pval.lor, lOR = lor, l.lor = lower.lor, 
                            u.lor = upper.lor, pval.lor = pval.lor, NNT = nnt), dig)
    invisible(out)
  }
}


# Formulas for computing effect sizes (d, r, log odds ratio)
# in designs with independent groups.

# Computing effect sizes d and g, independent groups
# (0) Study reported: 
# m.1 (post-test mean of treatment), m.2 (post-test mean of comparison),
# sd.1 (treatment standard deviation at post-test), sd.2 (comparison 
# standard deviation at post-test), n.1 (treatment), n.2 (comparison/control).

# EXAMPLE
# md <- data.frame(m.1=rnorm(30,30),m.2=rnorm(30,30),sd.1=rnorm(30,30),sd.2=rnorm(30,30),
#                  n.1=rnorm(30,30), n.2=rnorm(30,30), id=1:30)
# mes(m.1,m.2,sd.1,sd.2,n.1, n.2,level=95, cer=.2, dig=2, verbose=TRUE, id=id, data=md)

mes <- function(m.1,m.2,sd.1,sd.2,n.1, n.2,level=95, cer=.2, dig=2, verbose=TRUE, id=NULL, data=NULL) {
  if (!is.null(data)) {
    mf <- match.call()
    args <- match(c("m.1","m.2","sd.1","sd.2", "n.1", "n.2", "level", "dig", "id","data"), names(mf), 
                  0)
    mf <- mf[c(1, args)]
    mf$drop.unused.levels <- TRUE
    mf.m.1 <- mf[[match("m.1", names(mf))]]
    m.1 <- eval(mf.m.1, data, enclos = sys.frame(sys.parent()))
    mf.m.2 <- mf[[match("m.2", names(mf))]]
    m.2 <- eval(mf.m.2, data, enclos = sys.frame(sys.parent()))
    mf.sd.1 <- mf[[match("sd.1", names(mf))]]
    sd.1 <- eval(mf.sd.1, data, enclos = sys.frame(sys.parent()))
    mf.sd.2 <- mf[[match("sd.2", names(mf))]]
    sd.2 <- eval(mf.sd.2, data, enclos = sys.frame(sys.parent()))    
    mf.n.1 <- mf[[match("n.1", names(mf))]]
    n.1 <- eval(mf.n.1, data, enclos = sys.frame(sys.parent()))
    mf.n.2 <- mf[[match("n.2", names(mf))]]
    n.2 <- eval(mf.n.2, data, enclos = sys.frame(sys.parent()))
    mf.id <- mf[[match("id", names(mf))]]
    id <- eval(mf.id, data, enclos = sys.frame(sys.parent()))
  }
  s.within<-sqrt(((n.1-1)*sd.1^2+(n.2-1)*sd.2^2)/(n.1+n.2-2))
  d<-(m.1-m.2)/s.within
  var.d<-(n.1+n.2)/(n.1*n.2)+ (d^2)/(2*(n.1+n.2))
  df<- (n.1+n.2)-2
  j<-1-(3/(4*df-1))
  g<-j*d
  var.g<-j^2*var.d
  a <- ((n.1 + n.2)^2)/(n.1*n.2)  # will correct for inbalanced n, if applicable
  r <- d/sqrt((d^2) + a)
  var.r <- (a^2*var.d)/(d^2 + a)^3
  lor <- pi*d/sqrt(3)
  var.lor <- pi^2*var.d/3
  n <- n.1 + n.2
  z <-  0.5*log((1 + r)/(1-r))  #computing r to z' for each study
  var.z <- 1/(n-3) 
  # UPDATE: ADDED FOLLOWING TO EACH FUNCTION (FOR P-VAL, CI, U3, CLES, CLIFF'S DELTA, NNT)
  alpha <- (100 - level)/100
  crit <- qt(alpha/2, df, lower.tail = FALSE)
  #crit <- qnorm(alpha)
  zval.d <- d/sqrt(var.d)
  pval.d <- 2 * pt(abs(zval.d), df , lower.tail = FALSE)
  lower.d <- d - crit * sqrt(var.d)
  upper.d <- d + crit * sqrt(var.d)
  U3.d <- pnorm(d)*100
  cl.d<- (pnorm((d)/sqrt(2)))*100
  cliffs.d <- 2 * pnorm(d/sqrt(2)) - 1
  zval.g <- g/sqrt(var.g)
  pval.g <- 2 * pt(abs(zval.g), df , lower.tail = FALSE)
  lower.g <- g - crit * sqrt(var.g)
  upper.g <- g + crit * sqrt(var.g)
  U3.g <- pnorm(g)*100
  cl.g <- (pnorm((g)/sqrt(2)))*100
  
  zval.z <- z/sqrt(var.z)
  pval.z <- 2 * pt(abs(zval.z), df , lower.tail = FALSE)
  lower.z <- z - crit * sqrt(var.z)
  upper.z <- z + crit * sqrt(var.z)
  
  pval.r <- pval.z
  lower.r <- ztor(lower.z)
  upper.r <- ztor(upper.z)
  
  zval.lor <- lor/sqrt(var.lor)
  pval.lor <- 2 * pt(abs(zval.lor), df , lower.tail = FALSE)
  lower.lor <- lor - crit * sqrt(var.lor)
  upper.lor <- lor + crit * sqrt(var.lor)
  
  nnt <- 1/(pnorm(d- qnorm(1-cer))-cer)
  
  if (!is.null(data)) {
#     cat("\n")
#     message("    EFFECT SIZE CALCULATION (FOR VECTOR INPUT)")
#     cat("\n")
    
    out <- round(data.frame(id=id, N.total=n, n.1=n.1, n.2=n.2, d=d, var.d=var.d, l.d=lower.d, u.d=upper.d, U3.d=U3.d,
                            cl.d=cl.d, cliffs.d=cliffs.d, pval.d=pval.d, g=g, var.g=var.g,l.g=lower.g, 
                            u.g=upper.g, U3.g=U3.g, cl.g=cl.g, pval.g=pval.g, r=r, var.r=var.r, l.r=lower.r,
                            u.r=upper.r, pval.r=pval.r, fisher.z=z, var.z=var.z, l.z=lower.z, u.z=upper.z,
                            OR=exp(lor), l.or=exp(lower.lor), u.or=exp(upper.lor), pval.or=pval.lor,
                            lOR=lor, l.lor=lower.lor, u.lor=upper.lor, pval.lor=pval.lor,                          lOR=lor, l.lor=lower.lor, u.lor=upper.lor, pval.lor=pval.lor,
                            NNT=nnt), dig)
    return(out)
    
  }
  else{
#     cat("\n")
#     message("    EFFECT SIZE CALCULATION (FOR SINGLE INPUT)")
#     cat("\n")
    if (verbose) {
    cat(#"\n", "\n", "   EFFECT SIZE CALCULATION (FOR SINGLE INPUT)", "\n","\n",
      "Mean Differences ES:","\n", "\n", 
      "d [",level,"%CI] =", round(d, dig),
      "[",round(lower.d, dig),",",round(upper.d,dig),"]","\n",
      "\ var(d) =", round(var.d,dig), "\n",
      "\ p-value(d) =", round(pval.d,dig), "\n",
      "\ U3(d) =", round(U3.d, dig),"%", "\n", 
      "\ CLES(d) =", round(cl.d, dig),"%","\n",
      "\ Cliff's Delta =", round(cliffs.d, dig),"\n", "\n",
      "g [",level,"%CI] =", round(g, dig),
      "[",round(lower.g, dig),",",round(upper.g,dig),"]","\n",
      "\ var(g) =", round(var.g,dig), "\n",
      "\ p-value(g) =", round(pval.g,dig), "\n",
      "\ U3(g) =", round(U3.g, dig),"%", "\n", 
      "\ CLES(g) =", round(cl.g, dig),"%","\n", "\n",
      
      "Correlation ES:","\n", "\n", 
      "r [",level,"%CI] =", round(r, dig),
      "[",round(lower.r, dig),",",round(upper.r,dig),"]","\n",
      "\ var(r) =", round(var.r,dig), "\n",
      "\ p-value(r) =", round(pval.r,dig), "\n","\n",
      "z [",level,"%CI] =", round(z, dig),
      "[",round(lower.z, dig),",",round(upper.z,dig),"]","\n", 
      "\ var(z) =", round(var.z,dig), "\n",
      "\ p-value(z) =", round(pval.z,dig), "\n", "\n",
      
      "Odds Ratio ES:","\n", "\n", 
      "OR [",level,"%CI] =", round(exp(lor), dig),
      "[",round(exp(lower.lor), dig),",",round(exp(upper.lor),dig),"]","\n", 
      "\ p-value(OR) =", round(pval.lor,dig), "\n", "\n",
      "Log OR [",level,"%CI] =", round(lor, dig),
      "[",round(lower.lor, dig),",",round(upper.lor,dig),"]","\n",
      "\ var(lOR) =", round(var.lor,dig), "\n",
      "\ p-value(Log OR) =", round(pval.lor,dig), "\n", "\n", 
      
      "Other:","\n", "\n", 
      "NNT =", round(nnt, dig),"\n", 
      "Total N =", n
    )
    }
    out <- round(data.frame(N.total = n, n.1 = n.1, 
                            n.2 = n.2, d = d, var.d = var.d, l.d = lower.d, u.d = upper.d, 
                            U3.d = U3.d, cl.d = cl.d, cliffs.d = cliffs.d, pval.d = pval.d, 
                            g = g, var.g = var.g, l.g = lower.g, u.g = upper.g, 
                            U3.g = U3.g, cl.g = cl.g, pval.g = pval.g, r = r, 
                            var.r = var.r, l.r = lower.r, u.r = upper.r, pval.r = pval.r, 
                            fisher.z = z, var.z = var.z, l.z = lower.z, u.z = upper.z, 
                            OR = exp(lor), l.or = exp(lower.lor), u.or = exp(upper.lor), 
                            pval.or = pval.lor, lOR = lor, l.lor = lower.lor, 
                            u.lor = upper.lor, pval.lor = pval.lor, NNT = nnt), dig)
    invisible(out)
  }
}

# (1) Study reported: 
# m.1 (post-test mean of treatment), m.2 (post-test mean of comparison),
# s.pooled (pooled standard deviation), n.1 (treatment), 
# n.2 (comparison/control).

mes2 <- function(m.1,m.2,s.pooled,n.1, n.2,level=95, cer=.2, dig=2, verbose=TRUE, id=NULL, data=NULL) {
  if (!is.null(data)) {
    mf <- match.call()
    args <- match(c("m.1","m.2","s.pooled", "n.1", "n.2", "level", "dig", "id","data"), names(mf), 
                  0)
    mf <- mf[c(1, args)]
    mf$drop.unused.levels <- TRUE
    mf.m.1 <- mf[[match("m.1", names(mf))]]
    m.1 <- eval(mf.m.1, data, enclos = sys.frame(sys.parent()))
    mf.m.2 <- mf[[match("m.2", names(mf))]]
    m.2 <- eval(mf.m.2, data, enclos = sys.frame(sys.parent()))
    mf.s.pooled <- mf[[match("s.pooled", names(mf))]]
    s.pooled <- eval(mf.s.pooled, data, enclos = sys.frame(sys.parent()))
    mf.n.1 <- mf[[match("n.1", names(mf))]]
    n.1 <- eval(mf.n.1, data, enclos = sys.frame(sys.parent()))
    mf.n.2 <- mf[[match("n.2", names(mf))]]
    n.2 <- eval(mf.n.2, data, enclos = sys.frame(sys.parent()))
    mf.id <- mf[[match("id", names(mf))]]
    id <- eval(mf.id, data, enclos = sys.frame(sys.parent()))
  }
  d<-(m.1-m.2)/s.pooled
  var.d<-(n.1+n.2)/(n.1*n.2)+ (d^2)/(2*(n.1+n.2))
  df<- (n.1+n.2)-2
  j<-1-(3/(4*df-1))
  g<-j*d
  var.g<-j^2*var.d
  a <- ((n.1 + n.2)^2)/(n.1*n.2)  # will correct for inbalanced n, if applicable
  r <- d/sqrt((d^2) + a)
  var.r <- (a^2*var.d)/(d^2 + a)^3
  lor <- pi*d/sqrt(3)
  var.lor <- pi^2*var.d/3
  n <- n.1 + n.2
  z <-  0.5*log((1 + r)/(1-r))  
  var.z <- 1/(n-3) 
  # UPDATE: ADDED FOLLOWING TO EACH FUNCTION (FOR P-VAL, CI, U3, CLES, CLIFF'S DELTA, NNT)
  alpha <- (100 - level)/100
  crit <- qt(alpha/2, df, lower.tail = FALSE)
  #crit <- qnorm(alpha)
  zval.d <- d/sqrt(var.d)
  pval.d <- 2 * pt(abs(zval.d), df , lower.tail = FALSE)
  lower.d <- d - crit * sqrt(var.d)
  upper.d <- d + crit * sqrt(var.d)
  U3.d <- pnorm(d)*100
  cl.d<- (pnorm((d)/sqrt(2)))*100
  cliffs.d <- 2 * pnorm(d/sqrt(2)) - 1
  zval.g <- g/sqrt(var.g)
  pval.g <- 2 * pt(abs(zval.g), df , lower.tail = FALSE)
  lower.g <- g - crit * sqrt(var.g)
  upper.g <- g + crit * sqrt(var.g)
  U3.g <- pnorm(g)*100
  cl.g <- (pnorm((g)/sqrt(2)))*100
  
  zval.z <- z/sqrt(var.z)
  pval.z <- 2 * pt(abs(zval.z), df , lower.tail = FALSE)
  lower.z <- z - crit * sqrt(var.z)
  upper.z <- z + crit * sqrt(var.z)
  
  pval.r <- pval.z
  lower.r <- ztor(lower.z)
  upper.r <- ztor(upper.z)
  
  zval.lor <- lor/sqrt(var.lor)
  pval.lor <- 2 * pt(abs(zval.lor), df , lower.tail = FALSE)
  lower.lor <- lor - crit * sqrt(var.lor)
  upper.lor <- lor + crit * sqrt(var.lor)
  
  nnt <- 1/(pnorm(d- qnorm(1-cer))-cer)
  
  if (!is.null(data)) {
#     cat("\n")
#     message("    EFFECT SIZE CALCULATION (FOR VECTOR INPUT)")
#     cat("\n")
    
    out <- round(data.frame(id=id, N.total=n, n.1=n.1, n.2=n.2, d=d, var.d=var.d, l.d=lower.d, u.d=upper.d, U3.d=U3.d,
                            cl.d=cl.d, cliffs.d=cliffs.d, pval.d=pval.d, g=g, var.g=var.g,l.g=lower.g, 
                            u.g=upper.g, U3.g=U3.g, cl.g=cl.g, pval.g=pval.g, r=r, var.r=var.r, l.r=lower.r,
                            u.r=upper.r, pval.r=pval.r, fisher.z=z, var.z=var.z, l.z=lower.z, u.z=upper.z,
                            OR=exp(lor), l.or=exp(lower.lor), u.or=exp(upper.lor), pval.or=pval.lor, 
                            lOR=lor, l.lor=lower.lor, u.lor=upper.lor, pval.lor=pval.lor,
                            NNT=nnt), dig)
    return(out)
    
  }
  else{
#     cat("\n")
#     message("    EFFECT SIZE CALCULATION (FOR SINGLE INPUT)")
#     cat("\n")
    if (verbose) {
    cat(#"\n", "\n", "   EFFECT SIZE CALCULATION (FOR SINGLE INPUT)", "\n","\n",
      "Mean Differences ES:","\n", "\n", 
      "d [",level,"%CI] =", round(d, dig),
      "[",round(lower.d, dig),",",round(upper.d,dig),"]","\n",
      "\ var(d) =", round(var.d,dig), "\n",
      "\ p-value(d) =", round(pval.d,dig), "\n",
      "\ U3(d) =", round(U3.d, dig),"%", "\n", 
      "\ CLES(d) =", round(cl.d, dig),"%","\n",
      "\ Cliff's Delta =", round(cliffs.d, dig),"\n", "\n",
      "g [",level,"%CI] =", round(g, dig),
      "[",round(lower.g, dig),",",round(upper.g,dig),"]","\n",
      "\ var(g) =", round(var.g,dig), "\n",
      "\ p-value(g) =", round(pval.g,dig), "\n",
      "\ U3(g) =", round(U3.g, dig),"%", "\n", 
      "\ CLES(g) =", round(cl.g, dig),"%","\n", "\n",
      
      "Correlation ES:","\n", "\n", 
      "r [",level,"%CI] =", round(r, dig),
      "[",round(lower.r, dig),",",round(upper.r,dig),"]","\n",
      "\ var(r) =", round(var.r,dig), "\n",
      "\ p-value(r) =", round(pval.r,dig), "\n","\n",
      "z [",level,"%CI] =", round(z, dig),
      "[",round(lower.z, dig),",",round(upper.z,dig),"]","\n", 
      "\ var(z) =", round(var.z,dig), "\n",
      "\ p-value(z) =", round(pval.z,dig), "\n", "\n",
      
      "Odds Ratio ES:","\n", "\n", 
      "OR [",level,"%CI] =", round(exp(lor), dig),
      "[",round(exp(lower.lor), dig),",",round(exp(upper.lor),dig),"]","\n", 
      "\ p-value(OR) =", round(pval.lor,dig), "\n", "\n",
      "Log OR [",level,"%CI] =", round(lor, dig),
      "[",round(lower.lor, dig),",",round(upper.lor,dig),"]","\n",
      "\ var(lOR) =", round(var.lor,dig), "\n",
      "\ p-value(Log OR) =", round(pval.lor,dig), "\n", "\n", 
      
      "Other:","\n", "\n", 
      "NNT =", round(nnt, dig),"\n", 
      "Total N =", n
    )
    }
    out <- round(data.frame(N.total = n, n.1 = n.1, 
                            n.2 = n.2, d = d, var.d = var.d, l.d = lower.d, u.d = upper.d, 
                            U3.d = U3.d, cl.d = cl.d, cliffs.d = cliffs.d, pval.d = pval.d, 
                            g = g, var.g = var.g, l.g = lower.g, u.g = upper.g, 
                            U3.g = U3.g, cl.g = cl.g, pval.g = pval.g, r = r, 
                            var.r = var.r, l.r = lower.r, u.r = upper.r, pval.r = pval.r, 
                            fisher.z = z, var.z = var.z, l.z = lower.z, u.z = upper.z, 
                            OR = exp(lor), l.or = exp(lower.lor), u.or = exp(upper.lor), 
                            pval.or = pval.lor, lOR = lor, l.lor = lower.lor, 
                            u.lor = upper.lor, pval.lor = pval.lor, NNT = nnt), dig)
    invisible(out)
  }
}

# (2) Study reported: 
# t (t-test value of treatment v comparison), n.1 (treatment),
# n.2 (comparison/control).

tes <- function(t, n.1, n.2,level=95, cer=.2, dig=2, verbose=TRUE, id=NULL, data=NULL) {
  if (!is.null(data)) {
    mf <- match.call()
    args <- match(c("t", "n.1", "n.2", "level", "dig", "id","data"), names(mf), 
                  0)
    mf <- mf[c(1, args)]
    mf$drop.unused.levels <- TRUE
    mf.t <- mf[[match("t", names(mf))]]
    t <- eval(mf.t, data, enclos = sys.frame(sys.parent()))
    mf.n.1 <- mf[[match("n.1", names(mf))]]
    n.1 <- eval(mf.n.1, data, enclos = sys.frame(sys.parent()))
    mf.n.2 <- mf[[match("n.2", names(mf))]]
    n.2 <- eval(mf.n.2, data, enclos = sys.frame(sys.parent()))
    mf.id <- mf[[match("id", names(mf))]]
    id <- eval(mf.id, data, enclos = sys.frame(sys.parent()))
  } #If only total n reported just split .5
  d<-t*sqrt((n.1+n.2)/(n.1*n.2))
  var.d<-(n.1+n.2)/(n.1*n.2)+ (d^2)/(2*(n.1+n.2))
  df<- (n.1+n.2)-2
  j<-1-(3/(4*df-1))
  g<-j*d
  var.g<-j^2*var.d
  n <- n.1 + n.2
  r <- sqrt((t^2)/(t^2 + n-2))
  var.r <- ((1-r^2)^2)/(n-1)
  lor <- pi*d/sqrt(3)
  var.lor <- pi^2*var.d/3
  z <-  0.5*log((1 + r)/(1-r))  
  var.z <- 1/(n-3) 
  # UPDATE: ADDED FOLLOWING TO EACH FUNCTION (FOR P-VAL, CI, U3, CLES, CLIFF'S DELTA, NNT)
  alpha <- (100 - level)/100
  crit <- qt(alpha/2, df, lower.tail = FALSE)
  #crit <- qnorm(alpha)
  zval.d <- d/sqrt(var.d)
  pval.d <- 2 * pt(abs(zval.d), df , lower.tail = FALSE)
  lower.d <- d - crit * sqrt(var.d)
  upper.d <- d + crit * sqrt(var.d)
  U3.d <- pnorm(d)*100
  cl.d<- (pnorm((d)/sqrt(2)))*100
  cliffs.d <- 2 * pnorm(d/sqrt(2)) - 1
  zval.g <- g/sqrt(var.g)
  pval.g <- 2 * pt(abs(zval.g), df , lower.tail = FALSE)
  lower.g <- g - crit * sqrt(var.g)
  upper.g <- g + crit * sqrt(var.g)
  U3.g <- pnorm(g)*100
  cl.g <- (pnorm((g)/sqrt(2)))*100
  
  zval.z <- z/sqrt(var.z)
  pval.z <- 2 * pt(abs(zval.z), df , lower.tail = FALSE)
  lower.z <- z - crit * sqrt(var.z)
  upper.z <- z + crit * sqrt(var.z)
  
  pval.r <- pval.z
  lower.r <- ztor(lower.z)
  upper.r <- ztor(upper.z)
  
  zval.lor <- lor/sqrt(var.lor)
  pval.lor <- 2 * pt(abs(zval.lor), df , lower.tail = FALSE)
  lower.lor <- lor - crit * sqrt(var.lor)
  upper.lor <- lor + crit * sqrt(var.lor)
  
  nnt <- 1/(pnorm(d- qnorm(1-cer))-cer)
  
  if (!is.null(data)) {
#     cat("\n")
#     message("    EFFECT SIZE CALCULATION (FOR VECTOR INPUT)")
#     cat("\n")
    
    out <- round(data.frame(id=id, N.total=n, n.1=n.1, n.2=n.2, d=d, var.d=var.d, l.d=lower.d, u.d=upper.d, U3.d=U3.d,
                            cl.d=cl.d, cliffs.d=cliffs.d, pval.d=pval.d, g=g, var.g=var.g,l.g=lower.g, 
                            u.g=upper.g, U3.g=U3.g, cl.g=cl.g, pval.g=pval.g, r=r, var.r=var.r, l.r=lower.r,
                            u.r=upper.r, pval.r=pval.r, fisher.z=z, var.z=var.z, l.z=lower.z, u.z=upper.z,
                            OR=exp(lor), l.or=exp(lower.lor), u.or=exp(upper.lor), pval.or=pval.lor, 
                            lOR=lor, l.lor=lower.lor, u.lor=upper.lor, pval.lor=pval.lor,
                            NNT=nnt), dig)
    return(out)
    
  }
  else{
#     cat("\n")
#     message("    EFFECT SIZE CALCULATION (FOR SINGLE INPUT)")
#     cat("\n")
    if (verbose) {
    cat(#"\n", "\n", "   EFFECT SIZE CALCULATION (FOR SINGLE INPUT)", "\n","\n",
      "Mean Differences ES:","\n", "\n", 
      "d [",level,"%CI] =", round(d, dig),
      "[",round(lower.d, dig),",",round(upper.d,dig),"]","\n",
      "\ var(d) =", round(var.d,dig), "\n",
      "\ p-value(d) =", round(pval.d,dig), "\n",
      "\ U3(d) =", round(U3.d, dig),"%", "\n", 
      "\ CLES(d) =", round(cl.d, dig),"%","\n",
      "\ Cliff's Delta =", round(cliffs.d, dig),"\n", "\n",
      "g [",level,"%CI] =", round(g, dig),
      "[",round(lower.g, dig),",",round(upper.g,dig),"]","\n",
      "\ var(g) =", round(var.g,dig), "\n",
      "\ p-value(g) =", round(pval.g,dig), "\n",
      "\ U3(g) =", round(U3.g, dig),"%", "\n", 
      "\ CLES(g) =", round(cl.g, dig),"%","\n", "\n",
      
      "Correlation ES:","\n", "\n", 
      "r [",level,"%CI] =", round(r, dig),
      "[",round(lower.r, dig),",",round(upper.r,dig),"]","\n",
      "\ var(r) =", round(var.r,dig), "\n",
      "\ p-value(r) =", round(pval.r,dig), "\n","\n",
      "z [",level,"%CI] =", round(z, dig),
      "[",round(lower.z, dig),",",round(upper.z,dig),"]","\n",
      "\ var(z) =", round(var.z,dig), "\n",
      "\ p-value(z) =", round(pval.z,dig), "\n", "\n",
      
      "Odds Ratio ES:","\n", "\n", 
      "OR [",level,"%CI] =", round(exp(lor), dig),
      "[",round(exp(lower.lor), dig),",",round(exp(upper.lor),dig),"]","\n", 
      "\ p-value(OR) =", round(pval.lor,dig), "\n", "\n",
      "Log OR [",level,"%CI] =", round(lor, dig),
      "[",round(lower.lor, dig),",",round(upper.lor,dig),"]","\n",
      "\ var(lOR) =", round(var.lor,dig), "\n",
      "\ p-value(Log OR) =", round(pval.lor,dig), "\n", "\n", 
      
      "Other:","\n", "\n", 
      "NNT =", round(nnt, dig),"\n", 
      "Total N =", n
    )
    }
    out <- round(data.frame(N.total = n, n.1 = n.1, 
                            n.2 = n.2, d = d, var.d = var.d, l.d = lower.d, u.d = upper.d, 
                            U3.d = U3.d, cl.d = cl.d, cliffs.d = cliffs.d, pval.d = pval.d, 
                            g = g, var.g = var.g, l.g = lower.g, u.g = upper.g, 
                            U3.g = U3.g, cl.g = cl.g, pval.g = pval.g, r = r, 
                            var.r = var.r, l.r = lower.r, u.r = upper.r, pval.r = pval.r, 
                            fisher.z = z, var.z = var.z, l.z = lower.z, u.z = upper.z, 
                            OR = exp(lor), l.or = exp(lower.lor), u.or = exp(upper.lor), 
                            pval.or = pval.lor, lOR = lor, l.lor = lower.lor, 
                            u.lor = upper.lor, pval.lor = pval.lor, NNT = nnt), dig)
    invisible(out)
  }
}

# (3) Study reported: 
# f (F-test value of treatment v comparison), n.1 (treatment),
# n.2 (comparison/control).


fes <- function(f,n.1, n.2,level=95, cer=.2, dig=2, verbose=TRUE, id=NULL, data=NULL) {
  if (!is.null(data)) {
    mf <- match.call()
    args <- match(c("f", "n.1", "n.2", "level", "dig", "id","data"), names(mf), 
                  0)
    mf <- mf[c(1, args)]
    mf$drop.unused.levels <- TRUE
    mf.f <- mf[[match("f", names(mf))]]
    f <- eval(mf.f, data, enclos = sys.frame(sys.parent()))
    mf.n.1 <- mf[[match("n.1", names(mf))]]
    n.1 <- eval(mf.n.1, data, enclos = sys.frame(sys.parent()))
    mf.n.2 <- mf[[match("n.2", names(mf))]]
    n.2 <- eval(mf.n.2, data, enclos = sys.frame(sys.parent()))
    mf.id <- mf[[match("id", names(mf))]]
    id <- eval(mf.id, data, enclos = sys.frame(sys.parent()))
  }
  d<-sqrt(f*(n.1+n.2)/(n.1*n.2))
  var.d<-(n.1+n.2)/(n.1*n.2)+ (d^2)/(2*(n.1+n.2))
  df<- (n.1+n.2)-2
  j<-1-(3/(4*df-1))
  g<-j*d
  var.g<-j^2*var.d
  a <- ((n.1 + n.2)^2)/(n.1*n.2)  # will correct for inbalanced n, if applicable
  t <- sqrt(f)
  n <- n.1 + n.2
  #r<- sqrt((t^2)/(t^2 + n-2)) # to compute r from f with 1 df
  r <- d/sqrt((d^2) + a) # to compute r from d
  var.r <- (a^2*var.d)/(d^2 + a)^3
  lor <- pi*d/sqrt(3)
  var.lor <- pi^2*var.d/3
  z <-  0.5*log((1 + r)/(1-r))  
  var.z <- 1/(n-3) 
  # UPDATE: ADDED FOLLOWING TO EACH FUNCTION (FOR P-VAL, CI, U3, CLES, CLIFF'S DELTA, NNT)
  alpha <- (100 - level)/100
  crit <- qt(alpha/2, df, lower.tail = FALSE)
  #crit <- qnorm(alpha)
  zval.d <- d/sqrt(var.d)
  pval.d <- 2 * pt(abs(zval.d), df , lower.tail = FALSE)
  lower.d <- d - crit * sqrt(var.d)
  upper.d <- d + crit * sqrt(var.d)
  U3.d <- pnorm(d)*100
  cl.d<- (pnorm((d)/sqrt(2)))*100
  cliffs.d <- 2 * pnorm(d/sqrt(2)) - 1
  zval.g <- g/sqrt(var.g)
  pval.g <- 2 * pt(abs(zval.g), df , lower.tail = FALSE)
  lower.g <- g - crit * sqrt(var.g)
  upper.g <- g + crit * sqrt(var.g)
  U3.g <- pnorm(g)*100
  cl.g <- (pnorm((g)/sqrt(2)))*100
  
  zval.z <- z/sqrt(var.z)
  pval.z <- 2 * pt(abs(zval.z), df , lower.tail = FALSE)
  lower.z <- z - crit * sqrt(var.z)
  upper.z <- z + crit * sqrt(var.z)
  
  pval.r <- pval.z
  lower.r <- ztor(lower.z)
  upper.r <- ztor(upper.z)
  
  zval.lor <- lor/sqrt(var.lor)
  pval.lor <- 2 * pt(abs(zval.lor), df , lower.tail = FALSE)
  lower.lor <- lor - crit * sqrt(var.lor)
  upper.lor <- lor + crit * sqrt(var.lor)
  
  nnt <- 1/(pnorm(d- qnorm(1-cer))-cer)
  
  if (!is.null(data)) {
#     cat("\n")
#     message("    EFFECT SIZE CALCULATION (FOR VECTOR INPUT)")
#     cat("\n")
    
    out <- round(data.frame(id=id, N.total=n, n.1=n.1, n.2=n.2, d=d, var.d=var.d, l.d=lower.d, u.d=upper.d, U3.d=U3.d,
                            cl.d=cl.d, cliffs.d=cliffs.d, pval.d=pval.d, g=g, var.g=var.g,l.g=lower.g, 
                            u.g=upper.g, U3.g=U3.g, cl.g=cl.g, pval.g=pval.g, r=r, var.r=var.r, l.r=lower.r,
                            u.r=upper.r, pval.r=pval.r, fisher.z=z, var.z=var.z, l.z=lower.z, u.z=upper.z,
                            OR=exp(lor), l.or=exp(lower.lor), u.or=exp(upper.lor), pval.or=pval.lor, 
                            lOR=lor, l.lor=lower.lor, u.lor=upper.lor, pval.lor=pval.lor,
                            NNT=nnt), dig)
    return(out)
    
  }
  else{
#     cat("\n")
#     message("    EFFECT SIZE CALCULATION (FOR SINGLE INPUT)")
#     cat("\n")
    if (verbose) {
    cat(#"\n", "\n", "   EFFECT SIZE CALCULATION (FOR SINGLE INPUT)", "\n","\n",
      "Mean Differences ES:","\n", "\n", 
      "d [",level,"%CI] =", round(d, dig),
      "[",round(lower.d, dig),",",round(upper.d,dig),"]","\n",
      "\ var(d) =", round(var.d,dig), "\n",
      "\ p-value(d) =", round(pval.d,dig), "\n",
      "\ U3(d) =", round(U3.d, dig),"%", "\n", 
      "\ CLES(d) =", round(cl.d, dig),"%","\n",
      "\ Cliff's Delta =", round(cliffs.d, dig),"\n", "\n",
      "g [",level,"%CI] =", round(g, dig),
      "[",round(lower.g, dig),",",round(upper.g,dig),"]","\n",
      "\ var(g) =", round(var.g,dig), "\n",
      "\ p-value(g) =", round(pval.g,dig), "\n",
      "\ U3(g) =", round(U3.g, dig),"%", "\n", 
      "\ CLES(g) =", round(cl.g, dig),"%","\n", "\n",
      
      "Correlation ES:","\n", "\n", 
      "r [",level,"%CI] =", round(r, dig),
      "[",round(lower.r, dig),",",round(upper.r,dig),"]","\n",
      "\ var(r) =", round(var.r,dig), "\n",
      "\ p-value(r) =", round(pval.r,dig), "\n","\n",
      "z [",level,"%CI] =", round(z, dig),
      "[",round(lower.z, dig),",",round(upper.z,dig),"]","\n", 
      "\ var(z) =", round(var.z,dig), "\n",
      "\ p-value(z) =", round(pval.z,dig), "\n", "\n",
      
      "Odds Ratio ES:","\n", "\n", 
      "OR [",level,"%CI] =", round(exp(lor), dig),
      "[",round(exp(lower.lor), dig),",",round(exp(upper.lor),dig),"]","\n", 
      "\ p-value(OR) =", round(pval.lor,dig), "\n", "\n",
      "Log OR [",level,"%CI] =", round(lor, dig),
      "[",round(lower.lor, dig),",",round(upper.lor,dig),"]","\n",
      "\ var(lOR) =", round(var.lor,dig), "\n",
      "\ p-value(Log OR) =", round(pval.lor,dig), "\n", "\n", 
      
      "Other:","\n", "\n", 
      "NNT =", round(nnt, dig),"\n", 
      "Total N =", n
    )
    }
    out <- round(data.frame(N.total = n, n.1 = n.1, 
                            n.2 = n.2, d = d, var.d = var.d, l.d = lower.d, u.d = upper.d, 
                            U3.d = U3.d, cl.d = cl.d, cliffs.d = cliffs.d, pval.d = pval.d, 
                            g = g, var.g = var.g, l.g = lower.g, u.g = upper.g, 
                            U3.g = U3.g, cl.g = cl.g, pval.g = pval.g, r = r, 
                            var.r = var.r, l.r = lower.r, u.r = upper.r, pval.r = pval.r, 
                            fisher.z = z, var.z = var.z, l.z = lower.z, u.z = upper.z, 
                            OR = exp(lor), l.or = exp(lower.lor), u.or = exp(upper.lor), 
                            pval.or = pval.lor, lOR = lor, l.lor = lower.lor, 
                            u.lor = upper.lor, pval.lor = pval.lor, NNT = nnt), dig)
    invisible(out)
  }
}

# (4) Study reported: 
# p-value,  n.1 (treatment), n.2 (comparison/control), tail (one or two tailed?).

pes <- function(p, n.1, n.2, tail = "two",level=95, cer=.2, dig=2, verbose=TRUE, id=NULL, data=NULL) {
  if (!is.null(data)) {
    mf <- match.call()
    args <- match(c("p", "n.1", "n.2", "level", "dig", "id","data"), names(mf), 
                  0)
    mf <- mf[c(1, args)]
    mf$drop.unused.levels <- TRUE
    mf.p <- mf[[match("p", names(mf))]]
    p <- eval(mf.p, data, enclos = sys.frame(sys.parent()))
    mf.n.1 <- mf[[match("n.1", names(mf))]]
    n.1 <- eval(mf.n.1, data, enclos = sys.frame(sys.parent()))
    mf.n.2 <- mf[[match("n.2", names(mf))]]
    n.2 <- eval(mf.n.2, data, enclos = sys.frame(sys.parent()))
    mf.id <- mf[[match("id", names(mf))]]
    id <- eval(mf.id, data, enclos = sys.frame(sys.parent()))
  }
  n <- n.1 + n.2
  df<- (n.1+n.2)-2
  j<-1-(3/(4*df-1))
  a <- ((n.1 + n.2)^2)/(n.1*n.2) 
  if(tail == "one") {
    pxtwo<-p*2
    TINV<-qt((1-pxtwo/2),df)
    d<-TINV*sqrt((n.1+n.2)/(n.1*n.2))
    var.d<-(n.1+n.2)/(n.1*n.2)+ (d^2)/(2*(n.1+n.2))
    g<-j*d
    var.g<-j^2*var.d
    r <- d/sqrt((d^2) + a) # to compute r from d
    var.r <- (a^2*var.d)/(d^2 + a)^3
    lor <- pi*d/sqrt(3)
    var.lor <- pi^2*var.d/3
    z <-  0.5*log((1 + r)/(1-r))  
    var.z <- 1/(n-3) 
    
  }
  if(tail == "two") {
    TINV<-qt((1-p/2),df)
    d<-TINV*sqrt((n.1+n.2)/(n.1*n.2))
    var.d<-(n.1+n.2)/(n.1*n.2)+ (d^2)/(2*(n.1+n.2))
    g<-j*d
    var.g<-j^2*var.d
    r <- d/sqrt((d^2) + a) # to compute r from d
    var.r <- (a^2*var.d)/(d^2 + a)^3
    lor <- pi*d/sqrt(3)
    var.lor <- pi^2*var.d/3
    z <-  0.5*log((1 + r)/(1-r))  
    var.z <- 1/(n-3) 
  }
  # UPDATE: ADDED FOLLOWING TO EACH FUNCTION (FOR P-VAL, CI, U3, CLES, CLIFF'S DELTA, NNT)
  alpha <- (100 - level)/100
  crit <- qt(alpha/2, df, lower.tail = FALSE)
  #crit <- qnorm(alpha)
  zval.d <- d/sqrt(var.d)
  pval.d <- 2 * pt(abs(zval.d), df , lower.tail = FALSE)
  lower.d <- d - crit * sqrt(var.d)
  upper.d <- d + crit * sqrt(var.d)
  U3.d <- pnorm(d)*100
  cl.d<- (pnorm((d)/sqrt(2)))*100
  cliffs.d <- 2 * pnorm(d/sqrt(2)) - 1
  zval.g <- g/sqrt(var.g)
  pval.g <- 2 * pt(abs(zval.g), df , lower.tail = FALSE)
  lower.g <- g - crit * sqrt(var.g)
  upper.g <- g + crit * sqrt(var.g)
  U3.g <- pnorm(g)*100
  cl.g <- (pnorm((g)/sqrt(2)))*100
  
  zval.z <- z/sqrt(var.z)
  pval.z <- 2 * pt(abs(zval.z), df , lower.tail = FALSE)
  lower.z <- z - crit * sqrt(var.z)
  upper.z <- z + crit * sqrt(var.z)
  
  pval.r <- pval.z
  lower.r <- ztor(lower.z)
  upper.r <- ztor(upper.z)
  
  zval.lor <- lor/sqrt(var.lor)
  pval.lor <- 2 * pt(abs(zval.lor), df , lower.tail = FALSE)
  lower.lor <- lor - crit * sqrt(var.lor)
  upper.lor <- lor + crit * sqrt(var.lor)
  
  nnt <- 1/(pnorm(d- qnorm(1-cer))-cer)
  
  if (!is.null(data)) {
#     cat("\n")
#     message("    EFFECT SIZE CALCULATION (FOR VECTOR INPUT)")
#     cat("\n")
    
    out <- round(data.frame(id=id, N.total=n, n.1=n.1, n.2=n.2, d=d, var.d=var.d, l.d=lower.d, u.d=upper.d, U3.d=U3.d,
                            cl.d=cl.d, cliffs.d=cliffs.d, pval.d=pval.d, g=g, var.g=var.g,l.g=lower.g, 
                            u.g=upper.g, U3.g=U3.g, cl.g=cl.g, pval.g=pval.g, r=r, var.r=var.r, l.r=lower.r,
                            u.r=upper.r, pval.r=pval.r, fisher.z=z, var.z=var.z, l.z=lower.z, u.z=upper.z,
                            OR=exp(lor), l.or=exp(lower.lor), u.or=exp(upper.lor), pval.or=pval.lor, 
                            lOR=lor, l.lor=lower.lor, u.lor=upper.lor, pval.lor=pval.lor,
                            NNT=nnt), dig)
    return(out)
    
  }
  else{
#     cat("\n")
#     message("    EFFECT SIZE CALCULATION (FOR SINGLE INPUT)")
#     cat("\n")
    if (verbose) {
    cat(#"\n", "\n", "   EFFECT SIZE CALCULATION (FOR SINGLE INPUT)", "\n","\n",
      "Mean Differences ES:","\n", "\n", 
      "d [",level,"%CI] =", round(d, dig),
      "[",round(lower.d, dig),",",round(upper.d,dig),"]","\n",
      "\ var(d) =", round(var.d,dig), "\n",
      "\ p-value(d) =", round(pval.d,dig), "\n",
      "\ U3(d) =", round(U3.d, dig),"%", "\n", 
      "\ CLES(d) =", round(cl.d, dig),"%","\n",
      "\ Cliff's Delta =", round(cliffs.d, dig),"\n", "\n",
      "g [",level,"%CI] =", round(g, dig),
      "[",round(lower.g, dig),",",round(upper.g,dig),"]","\n",
      "\ var(g) =", round(var.g,dig), "\n",
      "\ p-value(g) =", round(pval.g,dig), "\n",
      "\ U3(g) =", round(U3.g, dig),"%", "\n", 
      "\ CLES(g) =", round(cl.g, dig),"%","\n", "\n",
      
      "Correlation ES:","\n", "\n", 
      "r [",level,"%CI] =", round(r, dig),
      "[",round(lower.r, dig),",",round(upper.r,dig),"]","\n",
      "\ var(r) =", round(var.r,dig), "\n",
      "\ p-value(r) =", round(pval.r,dig), "\n","\n",
      "z [",level,"%CI] =", round(z, dig),
      "[",round(lower.z, dig),",",round(upper.z,dig),"]","\n", 
      "\ var(z) =", round(var.z,dig), "\n",
      "\ p-value(z) =", round(pval.z,dig), "\n", "\n",
      
      "Odds Ratio ES:","\n", "\n", 
      "OR [",level,"%CI] =", round(exp(lor), dig),
      "[",round(exp(lower.lor), dig),",",round(exp(upper.lor),dig),"]","\n", 
      "\ p-value(OR) =", round(pval.lor,dig), "\n", "\n",
      "Log OR [",level,"%CI] =", round(lor, dig),
      "[",round(lower.lor, dig),",",round(upper.lor,dig),"]","\n",
      "\ var(lOR) =", round(var.lor,dig), "\n",
      "\ p-value(Log OR) =", round(pval.lor,dig), "\n", "\n", 
      
      "Other:","\n", "\n", 
      "NNT =", round(nnt, dig),"\n", 
      "Total N =", n
    )
    }
    out <- round(data.frame(N.total = n, n.1 = n.1, 
                            n.2 = n.2, d = d, var.d = var.d, l.d = lower.d, u.d = upper.d, 
                            U3.d = U3.d, cl.d = cl.d, cliffs.d = cliffs.d, pval.d = pval.d, 
                            g = g, var.g = var.g, l.g = lower.g, u.g = upper.g, 
                            U3.g = U3.g, cl.g = cl.g, pval.g = pval.g, r = r, 
                            var.r = var.r, l.r = lower.r, u.r = upper.r, pval.r = pval.r, 
                            fisher.z = z, var.z = var.z, l.z = lower.z, u.z = upper.z, 
                            OR = exp(lor), l.or = exp(lower.lor), u.or = exp(upper.lor), 
                            pval.or = pval.lor, lOR = lor, l.lor = lower.lor, 
                            u.lor = upper.lor, pval.lor = pval.lor, NNT = nnt), dig)
    invisible(out)
  }
}



#  Pearson r to effect size.
res <- function(r, var.r = NULL, n ,level=95, cer=.2, dig=2, verbose=TRUE, id=NULL, data=NULL) {
  if (!is.null(data)) {
    mf <- match.call()
    args <- match(c("r", "var.r", "n", "level", "dig", "id","data"), names(mf), 
                  0)
    mf <- mf[c(1, args)]
    mf$drop.unused.levels <- TRUE
    mf.r <- mf[[match("r", names(mf))]]
    r <- eval(mf.r, data, enclos = sys.frame(sys.parent()))
    mf.var.r <- mf[[match("var.r", names(mf))]]
    var.r <- eval(mf.var.r, data, enclos = sys.frame(sys.parent()))
    mf.n <- mf[[match("n", names(mf))]]
    n <- eval(mf.n, data, enclos = sys.frame(sys.parent()))
    mf.id <- mf[[match("id", names(mf))]]
    id <- eval(mf.id, data, enclos = sys.frame(sys.parent()))
  } # If var.r not reported use n
  if(is.null(var.r)){
    r <- r
    var.r <-(1-r^2)^2/(n-1) 
    #d<-2*r*sqrt((n-1)/(n*(1-r^2)))*abs(r)/r 
    d <- (2*r)/(sqrt(1-r^2))
    var.d <- 4*var.r/(1-r^2)^3
    lor <- pi*d/sqrt(3)
    var.lor <- pi^2*var.d/3
    z <-  0.5*log((1 + r)/(1-r))  
    var.z <- 1/(n-3) 
  }
  r <- r
  var.r <- var.r
  d <- 2*r/sqrt(1-r^2)
  var.d <- 4*var.r/(1-r^2)^3 
  lor <- pi*d/sqrt(3)
  var.lor <- pi^2*var.d/3
  z <-  0.5*log((1 + r)/(1-r))  
  var.z <- 1/(n-3) 
  df <- n-1
  # UPDATE: ADDED FOLLOWING TO EACH FUNCTION (FOR P-VAL, CI, U3, CLES, CLIFF'S DELTA, NNT)
  alpha <- (100 - level)/100
  crit <- qt(alpha/2, df, lower.tail = FALSE)
  #crit <- qnorm(alpha)
  zval.d <- d/sqrt(var.d)
  pval.d <- 2 * pt(abs(zval.d), df , lower.tail = FALSE)
  lower.d <- d - crit * sqrt(var.d)
  upper.d <- d + crit * sqrt(var.d)
  U3.d <- pnorm(d)*100
  cl.d<- (pnorm((d)/sqrt(2)))*100
   cliffs.d <- 2 * pnorm(d/sqrt(2)) - 1
#   zval.g <- g/sqrt(var.g)
#   pval.g <- 2 * pt(abs(zval.g), df , lower.tail = FALSE)
#   lower.g <- g - crit * sqrt(var.g)
#   upper.g <- g + crit * sqrt(var.g)
#   U3.g <- pnorm(g)*100
#   cl.g <- (pnorm((g)/sqrt(2)))*100
  
  zval.z <- z/sqrt(var.z)
  pval.z <- 2 * pt(abs(zval.z), df , lower.tail = FALSE)
  lower.z <- z - crit * sqrt(var.z)
  upper.z <- z + crit * sqrt(var.z)
  
  pval.r <- pval.z
  lower.r <- ztor(lower.z)
  upper.r <- ztor(upper.z)
  
  zval.lor <- lor/sqrt(var.lor)
  pval.lor <- 2 * pt(abs(zval.lor), df , lower.tail = FALSE)
  lower.lor <- lor - crit * sqrt(var.lor)
  upper.lor <- lor + crit * sqrt(var.lor)
  
  nnt <- 1/(pnorm(d- qnorm(1-cer))-cer)
  
  if (!is.null(data)) {
#     cat("\n")
#     message("    EFFECT SIZE CALCULATION (FOR VECTOR INPUT)")
#     cat("\n")
    
    out <- round(data.frame(id=id, N.total=n,  d=d, var.d=var.d, l.d=lower.d, u.d=upper.d, U3.d=U3.d,
                            cl.d=cl.d, cliffs.d=cliffs.d, pval.d=pval.d, 
                            #g=g, var.g=var.g,l.g=lower.g, u.g=upper.g, U3.g=U3.g, cl.g=cl.g, pval.g=pval.g,
                            r=r, var.r=var.r, l.r=lower.r,
                            u.r=upper.r, pval.r=pval.r, fisher.z=z, var.z=var.z, l.z=lower.z, u.z=upper.z,
                            OR=exp(lor), l.or=exp(lower.lor), u.or=exp(upper.lor), pval.or=pval.lor, 
                            lOR=lor, l.lor=lower.lor, u.lor=upper.lor, pval.lor=pval.lor,
                            NNT=nnt), dig)
    return(out)    
  }
  else{
#   cat("\n")
#   message("    EFFECT SIZE CALCULATION (FOR SINGLE INPUT)")
#   cat("\n")
    if (verbose) {
  cat(#"\n", "\n", "   EFFECT SIZE CALCULATION (FOR SINGLE INPUT)", "\n","\n",
    "Mean Differences ES:","\n", "\n", 
    "d [",level,"%CI] =", round(d, dig),
    "[",round(lower.d, dig),",",round(upper.d,dig),"]","\n",
    "\ var(d) =", round(var.d,dig), "\n",
    "\ p-value(d) =", round(pval.d,dig), "\n",
    "\ U3(d) =", round(U3.d, dig),"%", "\n", 
    "\ CLES(d) =", round(cl.d, dig),"%","\n",
    "\ Cliff's Delta =", round(cliffs.d, dig),"\n", "\n",
#     "g [",level,"%CI] =", round(g, dig),
#     "[",round(lower.g, dig),",",round(upper.g,dig),"]","\n",
#     "\ var(g) =", round(var.g,dig), "\n",
#     "\ p-value(g) =", round(pval.g,dig), "\n",
#     "\ U3(g) =", round(U3.g, dig),"%", "\n", 
#     "\ CLES(g) =", round(cl.g, dig),"%","\n", "\n",
    
    "Correlation ES:","\n", "\n", 
    "r [",level,"%CI] =", round(r, dig),
    "[",round(lower.r, dig),",",round(upper.r,dig),"]","\n",
    "\ var(r) =", round(var.r,dig), "\n",
    "\ p-value(r) =", round(pval.r,dig), "\n","\n",
    "z [",level,"%CI] =", round(z, dig),
    "[",round(lower.z, dig),",",round(upper.z,dig),"]","\n", 
    "\ var(z) =", round(var.z,dig), "\n",
    "\ p-value(z) =", round(pval.z,dig), "\n", "\n",
    
    "Odds Ratio ES:","\n", "\n", 
    "OR [",level,"%CI] =", round(exp(lor), dig),
    "[",round(exp(lower.lor), dig),",",round(exp(upper.lor),dig),"]","\n", 
    "\ p-value(OR) =", round(pval.lor,dig), "\n", "\n",
    "Log OR [",level,"%CI] =", round(lor, dig),
    "[",round(lower.lor, dig),",",round(upper.lor,dig),"]","\n",
    "\ var(lOR) =", round(var.lor,dig), "\n",
    "\ p-value(Log OR) =", round(pval.lor,dig), "\n", "\n", 
    
    "Other:","\n", "\n", 
    "NNT =", round(nnt, dig),"\n", 
    "Total N =", n
    )
    }
    out <- round(data.frame(N.total = n, d = d, var.d = var.d, l.d = lower.d, u.d = upper.d, 
                        U3.d = U3.d, cl.d = cl.d, cliffs.d = cliffs.d, pval.d = pval.d, 
                        #g = g, var.g = var.g, l.g = lower.g, u.g = upper.g, 
                        #U3.g = U3.g, cl.g = cl.g, pval.g = pval.g, 
                        r = r, 
                        var.r = var.r, l.r = lower.r, u.r = upper.r, pval.r = pval.r, 
                        fisher.z = z, var.z = var.z, l.z = lower.z, u.z = upper.z, 
                        OR = exp(lor), l.or = exp(lower.lor), u.or = exp(upper.lor), 
                        pval.or = pval.lor, lOR = lor, l.lor = lower.lor, 
                        u.lor = upper.lor, pval.lor = pval.lor, NNT = nnt), dig)
  invisible(out)
  }
}

# Formulas for computing effect sizes in designs with independent groups
# using ANCOVA. 
# Study reported: 
# m.1.adj (adjusted mean of treatment from ANCOVA),
# m.2.adj (adjusted mean of comparison/control from ANCOVA),
# s.adj (adjusted standard deviation), n.1 (treatment),
# n.2 (comparison/control), R (covariate outcome correlation or multiple
# correlation), q (number of covariates).

a.mes <- function(m.1.adj,m.2.adj,sd.adj,n.1, n.2, R, q,level=95, cer=.2, dig=2, verbose=TRUE, id=NULL, data=NULL) {
  if (!is.null(data)) {
    mf <- match.call()
    args <- match(c("m.1.adj","m.2.adj","sd.adj", "n.1", "n.2", "R", "q",
            "level", "dig", "id","data"), names(mf), 
            0)
    mf <- mf[c(1, args)]
    mf$drop.unused.levels <- TRUE
    mf.m.1 <- mf[[match("m.1.adj", names(mf))]]
    m.1.adj <- eval(mf.m.1, data, enclos = sys.frame(sys.parent()))
    mf.m.2 <- mf[[match("m.2.adj", names(mf))]]
    m.2.adj <- eval(mf.m.2, data, enclos = sys.frame(sys.parent()))
    mf.sd.adj <- mf[[match("sd.adj", names(mf))]]
    sd.adj <- eval(mf.sd.adj, data, enclos = sys.frame(sys.parent()))
    mf.n.1 <- mf[[match("n.1", names(mf))]]
    n.1 <- eval(mf.n.1, data, enclos = sys.frame(sys.parent()))
    mf.n.2 <- mf[[match("n.2", names(mf))]]
    n.2 <- eval(mf.n.2, data, enclos = sys.frame(sys.parent()))
    mf.R <- mf[[match("R", names(mf))]]
    R <- eval(mf.R, data, enclos = sys.frame(sys.parent()))
    mf.q <- mf[[match("q", names(mf))]]
    q <- eval(mf.q, data, enclos = sys.frame(sys.parent()))
    mf.id <- mf[[match("id", names(mf))]]
    id <- eval(mf.id, data, enclos = sys.frame(sys.parent()))
  }
   s.within<-sd.adj/sqrt(1-R^2)
   d<-(m.1.adj-m.2.adj)/s.within
   var.d<-((n.1+n.2)*(1-R^2))/(n.1*n.2)+ (d^2)/(2*(n.1+n.2))
   df<- (n.1+n.2)-2 - q
   j<-1-(3/(4*df-1))
   g<-j*d
   var.g<-j^2*var.d
   a <- ((n.1 + n.2)^2)/(n.1*n.2)  # will correct for inbalanced n, if applicable
   r <- d/sqrt((d^2) + a)
   var.r <- (a^2*var.d)/(d^2 + a)^3
   lor <- pi*d/sqrt(3)
   var.lor <- pi^2*var.d/3
   n <- n.1 + n.2
   z <-  0.5*log((1 + r)/(1-r))  #computing r to z' for each study
   var.z <- 1/(n-3) 
    # UPDATE: ADDED FOLLOWING TO EACH FUNCTION (FOR P-VAL, CI, U3, CLES, CLIFF'S DELTA, NNT)
  alpha <- (100 - level)/100
  crit <- qt(alpha/2, df, lower.tail = FALSE)
                    #crit <- qnorm(alpha)
                    zval.d <- d/sqrt(var.d)
                    pval.d <- 2 * pt(abs(zval.d), df , lower.tail = FALSE)
                    lower.d <- d - crit * sqrt(var.d)
                    upper.d <- d + crit * sqrt(var.d)
                    U3.d <- pnorm(d)*100
                    cl.d<- (pnorm((d)/sqrt(2)))*100
                    cliffs.d <- 2 * pnorm(d/sqrt(2)) - 1
                    zval.g <- g/sqrt(var.g)
                    pval.g <- 2 * pt(abs(zval.g), df , lower.tail = FALSE)
                    lower.g <- g - crit * sqrt(var.g)
                    upper.g <- g + crit * sqrt(var.g)
                    U3.g <- pnorm(g)*100
                    cl.g <- (pnorm((g)/sqrt(2)))*100
                    
                    zval.z <- z/sqrt(var.z)
                    pval.z <- 2 * pt(abs(zval.z), df , lower.tail = FALSE)
                    lower.z <- z - crit * sqrt(var.z)
                    upper.z <- z + crit * sqrt(var.z)
                    
                    pval.r <- pval.z
                    lower.r <- ztor(lower.z)
                    upper.r <- ztor(upper.z)
                    
                    zval.lor <- lor/sqrt(var.lor)
                    pval.lor <- 2 * pt(abs(zval.lor), df , lower.tail = FALSE)
                    lower.lor <- lor - crit * sqrt(var.lor)
                    upper.lor <- lor + crit * sqrt(var.lor)
                    
                    nnt <- 1/(pnorm(d- qnorm(1-cer))-cer)
                    
                    if (!is.null(data)) {
#                     cat("\n")
#                     message("    EFFECT SIZE CALCULATION (FOR VECTOR INPUT)")
#                     cat("\n")
                    
                    out <- round(data.frame(id=id, N.total=n, n.1=n.1, n.2=n.2, d=d, var.d=var.d, l.d=lower.d, u.d=upper.d, U3.d=U3.d,
                    cl.d=cl.d, cliffs.d=cliffs.d, pval.d=pval.d, g=g, var.g=var.g,l.g=lower.g, 
                    u.g=upper.g, U3.g=U3.g, cl.d=cl.d, pval.g=pval.g, r=r, var.r=var.r, l.r=lower.r,
                    u.r=upper.r, pval.r=pval.r, fisher.z=z, var.z=var.z, l.z=lower.z, u.z=upper.z,
                    OR=exp(lor), l.or=exp(lower.lor), u.or=exp(upper.lor), pval.or=pval.lor, 
                    lOR=lor, l.lor=lower.lor, u.lor=upper.lor, pval.lor=pval.lor,
                    NNT=nnt), dig)
                    return(out)
                    
      }
      else{
#         cat("\n")
#         message("    EFFECT SIZE CALCULATION (FOR SINGLE INPUT)")
#         cat("\n")
        if (verbose) {
        cat(#"\n", "\n", "   EFFECT SIZE CALCULATION (FOR SINGLE INPUT)", "\n","\n",
          "Mean Differences ES:","\n", "\n", 
          "d [",level,"%CI] =", round(d, dig),
          "[",round(lower.d, dig),",",round(upper.d,dig),"]","\n",
          "\ var(d) =", round(var.d,dig), "\n",
          "\ p-value(d) =", round(pval.d,dig), "\n",
          "\ U3(d) =", round(U3.d, dig),"%", "\n", 
          "\ CLES(d) =", round(cl.d, dig),"%","\n",
          "\ Cliff's Delta =", round(cliffs.d, dig),"\n", "\n",
          "g [",level,"%CI] =", round(g, dig),
          "[",round(lower.g, dig),",",round(upper.g,dig),"]","\n",
          "\ var(g) =", round(var.g,dig), "\n",
          "\ p-value(g) =", round(pval.g,dig), "\n",
          "\ U3(g) =", round(U3.g, dig),"%", "\n", 
          "\ CLES(g) =", round(cl.g, dig),"%","\n", "\n",
          
          "Correlation ES:","\n", "\n", 
          "r [",level,"%CI] =", round(r, dig),
          "[",round(lower.r, dig),",",round(upper.r,dig),"]","\n",
          "\ var(r) =", round(var.r,dig), "\n",
          "\ p-value(r) =", round(pval.r,dig), "\n","\n",
          "z [",level,"%CI] =", round(z, dig),
          "[",round(lower.z, dig),",",round(upper.z,dig),"]","\n",
          "\ var(z) =", round(var.z,dig), "\n",
          "\ p-value(z) =", round(pval.z,dig), "\n", "\n",
          
          "Odds Ratio ES:","\n", "\n", 
          "OR [",level,"%CI] =", round(exp(lor), dig),
          "[",round(exp(lower.lor), dig),",",round(exp(upper.lor),dig),"]","\n", 
          "\ p-value(OR) =", round(pval.lor,dig), "\n", "\n",
          "Log OR [",level,"%CI] =", round(lor, dig),
          "[",round(lower.lor, dig),",",round(upper.lor,dig),"]","\n",
          "\ var(lOR) =", round(var.lor,dig), "\n",
          "\ p-value(Log OR) =", round(pval.lor,dig), "\n", "\n", 
          
          "Other:","\n", "\n", 
          "NNT =", round(nnt, dig),"\n", 
          "Total N =", n
        )
        }
        out <- round(data.frame(N.total = n, n.1 = n.1, 
                                n.2 = n.2, d = d, var.d = var.d, l.d = lower.d, u.d = upper.d, 
                                U3.d = U3.d, cl.d = cl.d, cliffs.d = cliffs.d, pval.d = pval.d, 
                                g = g, var.g = var.g, l.g = lower.g, u.g = upper.g, 
                                U3.g = U3.g, cl.g = cl.g, pval.g = pval.g, r = r, 
                                var.r = var.r, l.r = lower.r, u.r = upper.r, pval.r = pval.r, 
                                fisher.z = z, var.z = var.z, l.z = lower.z, u.z = upper.z, 
                                OR = exp(lor), l.or = exp(lower.lor), u.or = exp(upper.lor), 
                                pval.or = pval.lor, lOR = lor, l.lor = lower.lor, 
                                u.lor = upper.lor, pval.lor = pval.lor, NNT = nnt), dig)
        invisible(out)
      }
}


# Study reported: 
# m.1.adj (adjusted mean of treatment from ANCOVA),
# m.2.adj (adjusted mean of comparison/control from ANCOVA),
# s.pooled (pooled standard deviation), n.1 (treatment), 
# n.2 (comparison/control), R (covariate outcome correlation or multiple
# correlation), q (number of covariates).

a.mes2 <- function(m.1.adj, m.2.adj, s.pooled, n.1, n.2, R, q,level=95, cer=.2, dig=2, verbose=TRUE, id=NULL, data=NULL) {
  if (!is.null(data)) {
    mf <- match.call()
    args <- match(c("m.1.adj","m.2.adj","s.pooled", "n.1", "n.2", "R", "q",
                    "level", "dig", "id","data"), names(mf), 
                  0)
    mf <- mf[c(1, args)]
    mf$drop.unused.levels <- TRUE
    mf.m.1 <- mf[[match("m.1.adj", names(mf))]]
    m.1.adj <- eval(mf.m.1, data, enclos = sys.frame(sys.parent()))
    mf.m.2 <- mf[[match("m.2.adj", names(mf))]]
    m.2.adj <- eval(mf.m.2, data, enclos = sys.frame(sys.parent()))
    mf.s.pooled <- mf[[match("s.pooled", names(mf))]]
    s.pooled <- eval(mf.s.pooled, data, enclos = sys.frame(sys.parent()))
    mf.n.1 <- mf[[match("n.1", names(mf))]]
    n.1 <- eval(mf.n.1, data, enclos = sys.frame(sys.parent()))
    mf.n.2 <- mf[[match("n.2", names(mf))]]
    n.2 <- eval(mf.n.2, data, enclos = sys.frame(sys.parent()))
    mf.R <- mf[[match("R", names(mf))]]
    R <- eval(mf.R, data, enclos = sys.frame(sys.parent()))
    mf.q <- mf[[match("q", names(mf))]]
    q <- eval(mf.q, data, enclos = sys.frame(sys.parent()))
    mf.id <- mf[[match("id", names(mf))]]
    id <- eval(mf.id, data, enclos = sys.frame(sys.parent()))
  }
  d<-(m.1.adj-m.2.adj)/s.pooled
  var.d<-((n.1+n.2)*(1-R^2))/(n.1*n.2)+ (d^2)/(2*(n.1+n.2))
  df<- (n.1+n.2)-2 - q
  j<-1-(3/(4*df-1))
  g<-j*d
  var.g<-j^2*var.d
  a <- ((n.1 + n.2)^2)/(n.1*n.2)  # will correct for inbalanced n, if applicable
  r <- d/sqrt((d^2) + a)
  var.r <- (a^2*var.d)/(d^2 + a)^3
  lor <- pi*d/sqrt(3)
  var.lor <- pi^2*var.d/3
  n <- n.1 + n.2
  z <-  0.5*log((1 + r)/(1-r))  #computing r to z' for each study
  var.z <- 1/(n-3) 
  # UPDATE: ADDED FOLLOWING TO EACH FUNCTION (FOR P-VAL, CI, U3, CLES, CLIFF'S DELTA, NNT)
  alpha <- (100 - level)/100
  crit <- qt(alpha/2, df, lower.tail = FALSE)
  #crit <- qnorm(alpha)
  zval.d <- d/sqrt(var.d)
  pval.d <- 2 * pt(abs(zval.d), df , lower.tail = FALSE)
  lower.d <- d - crit * sqrt(var.d)
  upper.d <- d + crit * sqrt(var.d)
  U3.d <- pnorm(d)*100
  cl.d<- (pnorm((d)/sqrt(2)))*100
  cliffs.d <- 2 * pnorm(d/sqrt(2)) - 1
  zval.g <- g/sqrt(var.g)
  pval.g <- 2 * pt(abs(zval.g), df , lower.tail = FALSE)
  lower.g <- g - crit * sqrt(var.g)
  upper.g <- g + crit * sqrt(var.g)
  U3.g <- pnorm(g)*100
  cl.g <- (pnorm((g)/sqrt(2)))*100
  
  zval.z <- z/sqrt(var.z)
  pval.z <- 2 * pt(abs(zval.z), df , lower.tail = FALSE)
  lower.z <- z - crit * sqrt(var.z)
  upper.z <- z + crit * sqrt(var.z)
  
  pval.r <- pval.z
  lower.r <- ztor(lower.z)
  upper.r <- ztor(upper.z)
  
  zval.lor <- lor/sqrt(var.lor)
  pval.lor <- 2 * pt(abs(zval.lor), df , lower.tail = FALSE)
  lower.lor <- lor - crit * sqrt(var.lor)
  upper.lor <- lor + crit * sqrt(var.lor)
  
  nnt <- 1/(pnorm(d- qnorm(1-cer))-cer)
  
  if (!is.null(data)) {
#     cat("\n")
#     message("    EFFECT SIZE CALCULATION (FOR VECTOR INPUT)")
#     cat("\n")
    
    out <- round(data.frame(id=id, N.total=n, n.1=n.1, n.2=n.2, d=d, var.d=var.d, l.d=lower.d, u.d=upper.d, U3.d=U3.d,
                            cl.d=cl.d, cliffs.d=cliffs.d, pval.d=pval.d, g=g, var.g=var.g,l.g=lower.g, 
                            u.g=upper.g, U3.g=U3.g, cl.g=cl.g, pval.g=pval.g, r=r, var.r=var.r, l.r=lower.r,
                            u.r=upper.r, pval.r=pval.r, fisher.z=z, var.z=var.z, l.z=lower.z, u.z=upper.z,
                            OR=exp(lor), l.or=exp(lower.lor), u.or=exp(upper.lor), pval.or=pval.lor, 
                            lOR=lor, l.lor=lower.lor, u.lor=upper.lor, pval.lor=pval.lor,
                            NNT=nnt), dig)
    return(out)
    
  }
  else{
#     cat("\n")
#     message("    EFFECT SIZE CALCULATION (FOR SINGLE INPUT)")
#     cat("\n")
    if (verbose) {
    cat(#"\n", "\n", "   EFFECT SIZE CALCULATION (FOR SINGLE INPUT)", "\n","\n",
      "Mean Differences ES:","\n", "\n", 
      "d [",level,"%CI] =", round(d, dig),
      "[",round(lower.d, dig),",",round(upper.d,dig),"]","\n",
      "\ var(d) =", round(var.d,dig), "\n",
      "\ p-value(d) =", round(pval.d,dig), "\n",
      "\ U3(d) =", round(U3.d, dig),"%", "\n", 
      "\ CLES(d) =", round(cl.d, dig),"%","\n",
      "\ Cliff's Delta =", round(cliffs.d, dig),"\n", "\n",
      "g [",level,"%CI] =", round(g, dig),
      "[",round(lower.g, dig),",",round(upper.g,dig),"]","\n",
      "\ var(g) =", round(var.g,dig), "\n",
      "\ p-value(g) =", round(pval.g,dig), "\n",
      "\ U3(g) =", round(U3.g, dig),"%", "\n", 
      "\ CLES(g) =", round(cl.g, dig),"%","\n", "\n",
      
      "Correlation ES:","\n", "\n", 
      "r [",level,"%CI] =", round(r, dig),
      "[",round(lower.r, dig),",",round(upper.r,dig),"]","\n",
      "\ var(r) =", round(var.r,dig), "\n",
      "\ p-value(r) =", round(pval.r,dig), "\n","\n",
      "z [",level,"%CI] =", round(z, dig),
      "[",round(lower.z, dig),",",round(upper.z,dig),"]","\n", 
      "\ var(z) =", round(var.z,dig), "\n",
      "\ p-value(z) =", round(pval.z,dig), "\n", "\n",
      
      "Odds Ratio ES:","\n", "\n", 
      "OR [",level,"%CI] =", round(exp(lor), dig),
      "[",round(exp(lower.lor), dig),",",round(exp(upper.lor),dig),"]","\n", 
      "\ p-value(OR) =", round(pval.lor,dig), "\n", "\n",
      "Log OR [",level,"%CI] =", round(lor, dig),
      "[",round(lower.lor, dig),",",round(upper.lor,dig),"]","\n",
      "\ var(lOR) =", round(var.lor,dig), "\n",
      "\ p-value(Log OR) =", round(pval.lor,dig), "\n", "\n", 
      
      "Other:","\n", "\n", 
      "NNT =", round(nnt, dig),"\n", 
      "Total N =", n
    )
    }
    out <- round(data.frame(N.total = n, n.1 = n.1, 
                            n.2 = n.2, d = d, var.d = var.d, l.d = lower.d, u.d = upper.d, 
                            U3.d = U3.d, cl.d = cl.d, cliffs.d = cliffs.d, pval.d = pval.d, 
                            g = g, var.g = var.g, l.g = lower.g, u.g = upper.g, 
                            U3.g = U3.g, cl.g = cl.g, pval.g = pval.g, r = r, 
                            var.r = var.r, l.r = lower.r, u.r = upper.r, pval.r = pval.r, 
                            fisher.z = z, var.z = var.z, l.z = lower.z, u.z = upper.z, 
                            OR = exp(lor), l.or = exp(lower.lor), u.or = exp(upper.lor), 
                            pval.or = pval.lor, lOR = lor, l.lor = lower.lor, 
                            u.lor = upper.lor, pval.lor = pval.lor, NNT = nnt), dig)
    invisible(out)
  }
}


# Study reported: 
# t (t-test value from ANCOVA), n.1 (treatment),
# n.2 (comparison/control), R (covariate outcome correlation or multiple
# correlation), q (number of covariates).

a.tes <- function(t, n.1, n.2, R, q,level=95, cer=.2, dig=2, verbose=TRUE, id=NULL, data=NULL) {
  if (!is.null(data)) {
    mf <- match.call()
    args <- match(c("t", "n.1", "n.2", "R", "q",
                    "level", "dig", "id","data"), names(mf), 
                  0)
    mf <- mf[c(1, args)]
    mf$drop.unused.levels <- TRUE
    mf.t <- mf[[match("t", names(mf))]]
    t <- eval(mf.t, data, enclos = sys.frame(sys.parent()))
    mf.n.1 <- mf[[match("n.1", names(mf))]]
    n.1 <- eval(mf.n.1, data, enclos = sys.frame(sys.parent()))
    mf.n.2 <- mf[[match("n.2", names(mf))]]
    n.2 <- eval(mf.n.2, data, enclos = sys.frame(sys.parent()))
    mf.R <- mf[[match("R", names(mf))]]
    R <- eval(mf.R, data, enclos = sys.frame(sys.parent()))
    mf.q <- mf[[match("q", names(mf))]]
    q <- eval(mf.q, data, enclos = sys.frame(sys.parent()))
    mf.id <- mf[[match("id", names(mf))]]
    id <- eval(mf.id, data, enclos = sys.frame(sys.parent()))
  }
  d<-t*sqrt((n.1+n.2)/(n.1*n.2))*sqrt(1-R^2)
  var.d<-((n.1+n.2)*(1-R^2))/(n.1*n.2)+ (d^2)/(2*(n.1+n.2))
  df<- (n.1+n.2)-2 - q
  j<-1-(3/(4*df-1))
  g<-j*d
  var.g<-j^2*var.d
  a <- ((n.1 + n.2)^2)/(n.1*n.2)  # will correct for inbalanced n, if applicable
  r <- d/sqrt((d^2) + a)
  var.r <- (a^2*var.d)/(d^2 + a)^3
  lor <- pi*d/sqrt(3)
  var.lor <- pi^2*var.d/3
  n <- n.1 + n.2
  z <-  0.5*log((1 + r)/(1-r))  #computing r to z' for each study
  var.z <- 1/(n-3) 
  # UPDATE: ADDED FOLLOWING TO EACH FUNCTION (FOR P-VAL, CI, U3, CLES, CLIFF'S DELTA, NNT)
  alpha <- (100 - level)/100
  crit <- qt(alpha/2, df, lower.tail = FALSE)
  #crit <- qnorm(alpha)
  zval.d <- d/sqrt(var.d)
  pval.d <- 2 * pt(abs(zval.d), df , lower.tail = FALSE)
  lower.d <- d - crit * sqrt(var.d)
  upper.d <- d + crit * sqrt(var.d)
  U3.d <- pnorm(d)*100
  cl.d<- (pnorm((d)/sqrt(2)))*100
  cliffs.d <- 2 * pnorm(d/sqrt(2)) - 1
  zval.g <- g/sqrt(var.g)
  pval.g <- 2 * pt(abs(zval.g), df , lower.tail = FALSE)
  lower.g <- g - crit * sqrt(var.g)
  upper.g <- g + crit * sqrt(var.g)
  U3.g <- pnorm(g)*100
  cl.g <- (pnorm((g)/sqrt(2)))*100
  
  zval.z <- z/sqrt(var.z)
  pval.z <- 2 * pt(abs(zval.z), df , lower.tail = FALSE)
  lower.z <- z - crit * sqrt(var.z)
  upper.z <- z + crit * sqrt(var.z)
  
  pval.r <- pval.z
  lower.r <- ztor(lower.z)
  upper.r <- ztor(upper.z)
  
  zval.lor <- lor/sqrt(var.lor)
  pval.lor <- 2 * pt(abs(zval.lor), df , lower.tail = FALSE)
  lower.lor <- lor - crit * sqrt(var.lor)
  upper.lor <- lor + crit * sqrt(var.lor)
  
  nnt <- 1/(pnorm(d- qnorm(1-cer))-cer)
  
  if (!is.null(data)) {
#     cat("\n")
#     message("    EFFECT SIZE CALCULATION (FOR VECTOR INPUT)")
#     cat("\n")
    
    out <- round(data.frame(id=id, N.total=n, n.1=n.1, n.2=n.2, d=d, var.d=var.d, l.d=lower.d, u.d=upper.d, U3.d=U3.d,
                            cl.d=cl.d, cliffs.d=cliffs.d, pval.d=pval.d, g=g, var.g=var.g,l.g=lower.g, 
                            u.g=upper.g, U3.g=U3.g, cl.g=cl.g, pval.g=pval.g, r=r, var.r=var.r, l.r=lower.r,
                            u.r=upper.r, pval.r=pval.r, fisher.z=z, var.z=var.z, l.z=lower.z, u.z=upper.z,
                            OR=exp(lor), l.or=exp(lower.lor), u.or=exp(upper.lor), pval.or=pval.lor, 
                            lOR=lor, l.lor=lower.lor, u.lor=upper.lor, pval.lor=pval.lor,
                            NNT=nnt), dig)
    return(out)
    
  }
  else{
#     cat("\n")
#     message("    EFFECT SIZE CALCULATION (FOR SINGLE INPUT)")
#     cat("\n")
    if (verbose) {
    cat(#"\n", "\n", "   EFFECT SIZE CALCULATION (FOR SINGLE INPUT)", "\n","\n",
      "Mean Differences ES:","\n", "\n", 
      "d [",level,"%CI] =", round(d, dig),
      "[",round(lower.d, dig),",",round(upper.d,dig),"]","\n",
      "\ var(d) =", round(var.d,dig), "\n",
      "\ p-value(d) =", round(pval.d,dig), "\n",
      "\ U3(d) =", round(U3.d, dig),"%", "\n", 
      "\ CLES(d) =", round(cl.d, dig),"%","\n",
      "\ Cliff's Delta =", round(cliffs.d, dig),"\n", "\n",
      "g [",level,"%CI] =", round(g, dig),
      "[",round(lower.g, dig),",",round(upper.g,dig),"]","\n",
      "\ var(g) =", round(var.g,dig), "\n",
      "\ p-value(g) =", round(pval.g,dig), "\n",
      "\ U3(g) =", round(U3.g, dig),"%", "\n", 
      "\ CLES(g) =", round(cl.g, dig),"%","\n", "\n",
      
      "Correlation ES:","\n", "\n", 
      "r [",level,"%CI] =", round(r, dig),
      "[",round(lower.r, dig),",",round(upper.r,dig),"]","\n",
      "\ var(r) =", round(var.r,dig), "\n",
      "\ p-value(r) =", round(pval.r,dig), "\n","\n",
      "z [",level,"%CI] =", round(z, dig),
      "[",round(lower.z, dig),",",round(upper.z,dig),"]","\n", 
      "\ var(z) =", round(var.z,dig), "\n",
      "\ p-value(z) =", round(pval.z,dig), "\n", "\n",
      
      "Odds Ratio ES:","\n", "\n", 
      "OR [",level,"%CI] =", round(exp(lor), dig),
      "[",round(exp(lower.lor), dig),",",round(exp(upper.lor),dig),"]","\n", 
      "\ p-value(OR) =", round(pval.lor,dig), "\n", "\n",
      "Log OR [",level,"%CI] =", round(lor, dig),
      "[",round(lower.lor, dig),",",round(upper.lor,dig),"]","\n",
      "\ var(lOR) =", round(var.lor,dig), "\n",
      "\ p-value(Log OR) =", round(pval.lor,dig), "\n", "\n", 
      
      "Other:","\n", "\n", 
      "NNT =", round(nnt, dig),"\n", 
      "Total N =", n
    )
    }
    out <- round(data.frame(N.total = n, n.1 = n.1, 
                            n.2 = n.2, d = d, var.d = var.d, l.d = lower.d, u.d = upper.d, 
                            U3.d = U3.d, cl.d = cl.d, cliffs.d = cliffs.d, pval.d = pval.d, 
                            g = g, var.g = var.g, l.g = lower.g, u.g = upper.g, 
                            U3.g = U3.g, cl.g = cl.g, pval.g = pval.g, r = r, 
                            var.r = var.r, l.r = lower.r, u.r = upper.r, pval.r = pval.r, 
                            fisher.z = z, var.z = var.z, l.z = lower.z, u.z = upper.z, 
                            OR = exp(lor), l.or = exp(lower.lor), u.or = exp(upper.lor), 
                            pval.or = pval.lor, lOR = lor, l.lor = lower.lor, 
                            u.lor = upper.lor, pval.lor = pval.lor, NNT = nnt), dig)
    invisible(out)
  }
}

# Study reported: 
# f (F-test value from ANCOVA) with independent groups, n.1 (treatment),
# n.2 (comparison/control),R (covariate outcome correlation or multiple
# correlation), q (number of covariates).

a.fes<-function(f,n.1, n.2, R, q,level=95, cer=.2, dig=2, verbose=TRUE, id=NULL, data=NULL) {
  if (!is.null(data)) {
    mf <- match.call()
    args <- match(c("f", "n.1", "n.2", "R", "q",
                    "level", "dig", "id","data"), names(mf), 
                  0)
    mf <- mf[c(1, args)]
    mf$drop.unused.levels <- TRUE
    mf.f <- mf[[match("f", names(mf))]]
    f <- eval(mf.f, data, enclos = sys.frame(sys.parent()))
    mf.n.1 <- mf[[match("n.1", names(mf))]]
    n.1 <- eval(mf.n.1, data, enclos = sys.frame(sys.parent()))
    mf.n.2 <- mf[[match("n.2", names(mf))]]
    n.2 <- eval(mf.n.2, data, enclos = sys.frame(sys.parent()))
    mf.R <- mf[[match("R", names(mf))]]
    R <- eval(mf.R, data, enclos = sys.frame(sys.parent()))
    mf.q <- mf[[match("q", names(mf))]]
    q <- eval(mf.q, data, enclos = sys.frame(sys.parent()))
    mf.id <- mf[[match("id", names(mf))]]
    id <- eval(mf.id, data, enclos = sys.frame(sys.parent()))
  }
  d<-sqrt(f*(n.1+n.2)/(n.1*n.2))*sqrt(1-R^2)
  var.d<-((n.1+n.2)*(1-R^2))/(n.1*n.2)+ (d^2)/(2*(n.1+n.2))
  df<- (n.1+n.2)-2 - q
  j<-1-(3/(4*df-1))
  g<-j*d
  var.g<-j^2*var.d
  a <- ((n.1 + n.2)^2)/(n.1*n.2)  # will correct for inbalanced n, if applicable
  r <- d/sqrt((d^2) + a)
  var.r <- (a^2*var.d)/(d^2 + a)^3
  lor <- pi*d/sqrt(3)
  var.lor <- pi^2*var.d/3
  n <- n.1 + n.2
  z <-  0.5*log((1 + r)/(1-r))  #computing r to z' for each study
  var.z <- 1/(n-3) 
  # UPDATE: ADDED FOLLOWING TO EACH FUNCTION (FOR P-VAL, CI, U3, CLES, CLIFF'S DELTA, NNT)
  alpha <- (100 - level)/100
  crit <- qt(alpha/2, df, lower.tail = FALSE)
  #crit <- qnorm(alpha)
  zval.d <- d/sqrt(var.d)
  pval.d <- 2 * pt(abs(zval.d), df , lower.tail = FALSE)
  lower.d <- d - crit * sqrt(var.d)
  upper.d <- d + crit * sqrt(var.d)
  U3.d <- pnorm(d)*100
  cl.d<- (pnorm((d)/sqrt(2)))*100
  cliffs.d <- 2 * pnorm(d/sqrt(2)) - 1
  zval.g <- g/sqrt(var.g)
  pval.g <- 2 * pt(abs(zval.g), df , lower.tail = FALSE)
  lower.g <- g - crit * sqrt(var.g)
  upper.g <- g + crit * sqrt(var.g)
  U3.g <- pnorm(g)*100
  cl.g <- (pnorm((g)/sqrt(2)))*100
  
  zval.z <- z/sqrt(var.z)
  pval.z <- 2 * pt(abs(zval.z), df , lower.tail = FALSE)
  lower.z <- z - crit * sqrt(var.z)
  upper.z <- z + crit * sqrt(var.z)
  
  pval.r <- pval.z
  lower.r <- ztor(lower.z)
  upper.r <- ztor(upper.z)
  
  zval.lor <- lor/sqrt(var.lor)
  pval.lor <- 2 * pt(abs(zval.lor), df , lower.tail = FALSE)
  lower.lor <- lor - crit * sqrt(var.lor)
  upper.lor <- lor + crit * sqrt(var.lor)
  
  nnt <- 1/(pnorm(d- qnorm(1-cer))-cer)
  
  if (!is.null(data)) {
#     cat("\n")
#     message("    EFFECT SIZE CALCULATION (FOR VECTOR INPUT)")
#     cat("\n")
    
    out <- round(data.frame(id=id, N.total=n, n.1=n.1, n.2=n.2, d=d, var.d=var.d, l.d=lower.d, u.d=upper.d, U3.d=U3.d,
                            cl.d=cl.d, cliffs.d=cliffs.d, pval.d=pval.d, g=g, var.g=var.g,l.g=lower.g, 
                            u.g=upper.g, U3.g=U3.g, cl.g=cl.g, pval.g=pval.g, r=r, var.r=var.r, l.r=lower.r,
                            u.r=upper.r, pval.r=pval.r, fisher.z=z, var.z=var.z, l.z=lower.z, u.z=upper.z,
                            OR=exp(lor), l.or=exp(lower.lor), u.or=exp(upper.lor), pval.or=pval.lor, 
                            lOR=lor, l.lor=lower.lor, u.lor=upper.lor, pval.lor=pval.lor,
                            NNT=nnt), dig)
    return(out)
    
  }
  else{
#     cat("\n")
#     message("    EFFECT SIZE CALCULATION (FOR SINGLE INPUT)")
#     cat("\n")
    if (verbose) {
    cat(#"\n", "\n", "   EFFECT SIZE CALCULATION (FOR SINGLE INPUT)", "\n","\n",
      "Mean Differences ES:","\n", "\n", 
      "d [",level,"%CI] =", round(d, dig),
      "[",round(lower.d, dig),",",round(upper.d,dig),"]","\n",
      "\ var(d) =", round(var.d,dig), "\n",
      "\ p-value(d) =", round(pval.d,dig), "\n",
      "\ U3(d) =", round(U3.d, dig),"%", "\n", 
      "\ CLES(d) =", round(cl.d, dig),"%","\n",
      "\ Cliff's Delta =", round(cliffs.d, dig),"\n", "\n",
      "g [",level,"%CI] =", round(g, dig),
      "[",round(lower.g, dig),",",round(upper.g,dig),"]","\n",
      "\ var(g) =", round(var.g,dig), "\n",
      "\ p-value(g) =", round(pval.g,dig), "\n",
      "\ U3(g) =", round(U3.g, dig),"%", "\n", 
      "\ CLES(g) =", round(cl.g, dig),"%","\n", "\n",
      
      "Correlation ES:","\n", "\n", 
      "r [",level,"%CI] =", round(r, dig),
      "[",round(lower.r, dig),",",round(upper.r,dig),"]","\n",
      "\ var(r) =", round(var.r,dig), "\n",
      "\ p-value(r) =", round(pval.r,dig), "\n","\n",
      "z [",level,"%CI] =", round(z, dig),
      "[",round(lower.z, dig),",",round(upper.z,dig),"]","\n", 
      "\ var(z) =", round(var.z,dig), "\n",
      "\ p-value(z) =", round(pval.z,dig), "\n", "\n",
      
      "Odds Ratio ES:","\n", "\n", 
      "OR [",level,"%CI] =", round(exp(lor), dig),
      "[",round(exp(lower.lor), dig),",",round(exp(upper.lor),dig),"]","\n", 
      "\ p-value(OR) =", round(pval.lor,dig), "\n", "\n",
      "Log OR [",level,"%CI] =", round(lor, dig),
      "[",round(lower.lor, dig),",",round(upper.lor,dig),"]","\n",
      "\ var(lOR) =", round(var.lor,dig), "\n",
      "\ p-value(Log OR) =", round(pval.lor,dig), "\n", "\n", 
      
      "Other:","\n", "\n", 
      "NNT =", round(nnt, dig),"\n", 
      "Total N =", n
    )
    }
    out <- round(data.frame(N.total = n, n.1 = n.1, 
                            n.2 = n.2, d = d, var.d = var.d, l.d = lower.d, u.d = upper.d, 
                            U3.d = U3.d, cl.d = cl.d, cliffs.d = cliffs.d, pval.d = pval.d, 
                            g = g, var.g = var.g, l.g = lower.g, u.g = upper.g, 
                            U3.g = U3.g, cl.g = cl.g, pval.g = pval.g, r = r, 
                            var.r = var.r, l.r = lower.r, u.r = upper.r, pval.r = pval.r, 
                            fisher.z = z, var.z = var.z, l.z = lower.z, u.z = upper.z, 
                            OR = exp(lor), l.or = exp(lower.lor), u.or = exp(upper.lor), 
                            pval.or = pval.lor, lOR = lor, l.lor = lower.lor, 
                            u.lor = upper.lor, pval.lor = pval.lor, NNT = nnt), dig)
    invisible(out)
  }
}

# Study reported: 
# p-value (for ONE-tailed test, from ANCOVA), n.1 (treatment), 
# n.2 (comparison/control), R (covariate outcome correlation or multiple
# correlation), q (number of covariates).

a.pes <- function(p, n.1, n.2, R, q, tail = "two",level=95, cer=.2, dig=2, verbose=TRUE, id=NULL, data=NULL) {
  if (!is.null(data)) {
    mf <- match.call()
    args <- match(c("p", "n.1", "n.2", "R", "q",
                    "level", "dig", "id","data"), names(mf), 
                  0)
    mf <- mf[c(1, args)]
    mf$drop.unused.levels <- TRUE
    mf.p <- mf[[match("p", names(mf))]]
    p<- eval(mf.p, data, enclos = sys.frame(sys.parent()))
    mf.n.1 <- mf[[match("n.1", names(mf))]]
    n.1 <- eval(mf.n.1, data, enclos = sys.frame(sys.parent()))
    mf.n.2 <- mf[[match("n.2", names(mf))]]
    n.2 <- eval(mf.n.2, data, enclos = sys.frame(sys.parent()))
    mf.R <- mf[[match("R", names(mf))]]
    R <- eval(mf.R, data, enclos = sys.frame(sys.parent()))
    mf.q <- mf[[match("q", names(mf))]]
    q <- eval(mf.q, data, enclos = sys.frame(sys.parent()))
    mf.id <- mf[[match("id", names(mf))]]
    id <- eval(mf.id, data, enclos = sys.frame(sys.parent()))
  }
  n <- n.1 + n.2
  df<- (n.1+n.2)-2 - q
  j<-1-(3/(4*df-1))
  a <- ((n.1 + n.2)^2)/(n.1*n.2) 
  if(tail == "one") {
    pxtwo<-p*2
    TINV<-qt((1-pxtwo/2),df)
    d<-TINV*sqrt((n.1+n.2)/(n.1*n.2))*sqrt(1-R^2)
    var.d<-((n.1+n.2)*(1-R^2))/(n.1*n.2)+ (d^2)/(2*(n.1+n.2))
    g<-j*d
    var.g<-j^2*var.d
    r <- d/sqrt((d^2) + a) # to compute r from d
    var.r <- (a^2*var.d)/(d^2 + a)^3
    lor <- pi*d/sqrt(3)
    var.lor <- pi^2*var.d/3
    z <-  0.5*log((1 + r)/(1-r))  
    var.z <- 1/(n-3) 
  }
  if(tail == "two") {
    TINV<-qt((1-p/2),df)
    d<-TINV*sqrt((n.1+n.2)/(n.1*n.2))*sqrt(1-R^2)
    var.d<-((n.1+n.2)*(1-R^2))/(n.1*n.2)+ (d^2)/(2*(n.1+n.2))
    g<-j*d
    var.g<-j^2*var.d
    r <- d/sqrt((d^2) + a) # to compute r from d
    var.r <- (a^2*var.d)/(d^2 + a)^3
    lor <- pi*d/sqrt(3)
    var.lor <- pi^2*var.d/3
    z <-  0.5*log((1 + r)/(1-r))  
    var.z <- 1/(n-3) 
  }
  # UPDATE: ADDED FOLLOWING TO EACH FUNCTION (FOR P-VAL, CI, U3, CLES, CLIFF'S DELTA, NNT)
  alpha <- (100 - level)/100
  crit <- qt(alpha/2, df, lower.tail = FALSE)
  #crit <- qnorm(alpha)
  zval.d <- d/sqrt(var.d)
  pval.d <- 2 * pt(abs(zval.d), df , lower.tail = FALSE)
  lower.d <- d - crit * sqrt(var.d)
  upper.d <- d + crit * sqrt(var.d)
  U3.d <- pnorm(d)*100
  cl.d<- (pnorm((d)/sqrt(2)))*100
  cliffs.d <- 2 * pnorm(d/sqrt(2)) - 1
  zval.g <- g/sqrt(var.g)
  pval.g <- 2 * pt(abs(zval.g), df , lower.tail = FALSE)
  lower.g <- g - crit * sqrt(var.g)
  upper.g <- g + crit * sqrt(var.g)
  U3.g <- pnorm(g)*100
  cl.g <- (pnorm((g)/sqrt(2)))*100
  
  zval.z <- z/sqrt(var.z)
  pval.z <- 2 * pt(abs(zval.z), df , lower.tail = FALSE)
  lower.z <- z - crit * sqrt(var.z)
  upper.z <- z + crit * sqrt(var.z)
  
  pval.r <- pval.z
  lower.r <- ztor(lower.z)
  upper.r <- ztor(upper.z)
  
  zval.lor <- lor/sqrt(var.lor)
  pval.lor <- 2 * pt(abs(zval.lor), df , lower.tail = FALSE)
  lower.lor <- lor - crit * sqrt(var.lor)
  upper.lor <- lor + crit * sqrt(var.lor)
  
  nnt <- 1/(pnorm(d- qnorm(1-cer))-cer)
  
  if (!is.null(data)) {
#     cat("\n")
#     message("    EFFECT SIZE CALCULATION (FOR VECTOR INPUT)")
#     cat("\n")
    
    out <- round(data.frame(id=id, N.total=n, n.1=n.1, n.2=n.2, d=d, var.d=var.d, l.d=lower.d, u.d=upper.d, U3.d=U3.d,
                            cl.d=cl.d, cliffs.d=cliffs.d, pval.d=pval.d, g=g, var.g=var.g,l.g=lower.g, 
                            u.g=upper.g, U3.g=U3.g, cl.g=cl.g, pval.g=pval.g, r=r, var.r=var.r, l.r=lower.r,
                            u.r=upper.r, pval.r=pval.r, fisher.z=z, var.z=var.z, l.z=lower.z, u.z=upper.z,
                            OR=exp(lor), l.or=exp(lower.lor), u.or=exp(upper.lor), pval.or=pval.lor, 
                            lOR=lor, l.lor=lower.lor, u.lor=upper.lor, pval.lor=pval.lor,
                            NNT=nnt), dig)
    return(out)
    
  }
  else{
#     cat("\n")
#     message("    EFFECT SIZE CALCULATION (FOR SINGLE INPUT)")
#     cat("\n")
    if (verbose) {
    cat(#"\n", "\n", "   EFFECT SIZE CALCULATION (FOR SINGLE INPUT)", "\n","\n",
      "Mean Differences ES:","\n", "\n", 
      "d [",level,"%CI] =", round(d, dig),
      "[",round(lower.d, dig),",",round(upper.d,dig),"]","\n",
      "\ var(d) =", round(var.d,dig), "\n",
      "\ p-value(d) =", round(pval.d,dig), "\n",
      "\ U3(d) =", round(U3.d, dig),"%", "\n", 
      "\ CLES(d) =", round(cl.d, dig),"%","\n",
      "\ Cliff's Delta =", round(cliffs.d, dig),"\n", "\n",
      "g [",level,"%CI] =", round(g, dig),
      "[",round(lower.g, dig),",",round(upper.g,dig),"]","\n",
      "\ var(g) =", round(var.g,dig), "\n",
      "\ p-value(g) =", round(pval.g,dig), "\n",
      "\ U3(g) =", round(U3.g, dig),"%", "\n", 
      "\ CLES(g) =", round(cl.g, dig),"%","\n", "\n",
      
      "Correlation ES:","\n", "\n", 
      "r [",level,"%CI] =", round(r, dig),
      "[",round(lower.r, dig),",",round(upper.r,dig),"]","\n",
      "\ var(r) =", round(var.r,dig), "\n",
      "\ p-value(r) =", round(pval.r,dig), "\n","\n",
      "z [",level,"%CI] =", round(z, dig),
      "[",round(lower.z, dig),",",round(upper.z,dig),"]","\n", 
      "\ var(z) =", round(var.z,dig), "\n",
      "\ p-value(z) =", round(pval.z,dig), "\n", "\n",
      
      "Odds Ratio ES:","\n", "\n", 
      "OR [",level,"%CI] =", round(exp(lor), dig),
      "[",round(exp(lower.lor), dig),",",round(exp(upper.lor),dig),"]","\n", 
      "\ p-value(OR) =", round(pval.lor,dig), "\n", "\n",
      "Log OR [",level,"%CI] =", round(lor, dig),
      "[",round(lower.lor, dig),",",round(upper.lor,dig),"]","\n",
      "\ var(lOR) =", round(var.lor,dig), "\n",
      "\ p-value(Log OR) =", round(pval.lor,dig), "\n", "\n", 
      
      "Other:","\n", "\n", 
      "NNT =", round(nnt, dig),"\n", 
      "Total N =", n
    )
    }
    out <- round(data.frame(N.total = n, n.1 = n.1, 
                            n.2 = n.2, d = d, var.d = var.d, l.d = lower.d, u.d = upper.d, 
                            U3.d = U3.d, cl.d = cl.d, cliffs.d = cliffs.d, pval.d = pval.d, 
                            g = g, var.g = var.g, l.g = lower.g, u.g = upper.g, 
                            U3.g = U3.g, cl.g = cl.g, pval.g = pval.g, r = r, 
                            var.r = var.r, l.r = lower.r, u.r = upper.r, pval.r = pval.r, 
                            fisher.z = z, var.z = var.z, l.z = lower.z, u.z = upper.z, 
                            OR = exp(lor), l.or = exp(lower.lor), u.or = exp(upper.lor), 
                            pval.or = pval.lor, lOR = lor, l.lor = lower.lor, 
                            u.lor = upper.lor, pval.lor = pval.lor, NNT = nnt), dig)
    invisible(out)
  }
}

# computing es from log odds ratio

lores <- function(lor, var.lor, n.1, n.2,level=95, cer=.2, dig=2, verbose=TRUE, id=NULL, data=NULL) {
  if (!is.null(data)) {
    mf <- match.call()
    args <- match(c("lor", "var.lor", "n.1", "n.2", "level", "dig", "id","data"), names(mf), 
                  0)
    mf <- mf[c(1, args)]
    mf$drop.unused.levels <- TRUE
    mf.lor <- mf[[match("lor", names(mf))]]
    lor <- eval(mf.lor, data, enclos = sys.frame(sys.parent()))
    mf.var.lor <- mf[[match("var.lor", names(mf))]]
    var.lor <- eval(mf.var.lor, data, enclos = sys.frame(sys.parent()))
    mf.n.1 <- mf[[match("n.1", names(mf))]]
    n.1 <- eval(mf.n.1, data, enclos = sys.frame(sys.parent()))
    mf.n.2 <- mf[[match("n.2", names(mf))]]
    n.2 <- eval(mf.n.2, data, enclos = sys.frame(sys.parent()))
    mf.id <- mf[[match("id", names(mf))]]
    id <- eval(mf.id, data, enclos = sys.frame(sys.parent()))
  } # n.1 =  tmt grp
  d <- lor*sqrt(3)/pi
  var.d <- 3*var.lor/pi^2
  df<- (n.1+n.2)-2 
  j<-1-(3/(4*df-1))
  g<-j*d
  var.g<-j^2*var.d
  a <- ((n.1 + n.2)^2)/(n.1*n.2)  # will correct for inbalanced n, if applicable
  r <- d/sqrt((d^2) + a)
  var.r <- (a^2*var.d)/(d^2 + a)^3
  lor <- pi*d/sqrt(3)
  var.lor <- pi^2*var.d/3
  n <- n.1 + n.2
  z <-  0.5*log((1 + r)/(1-r))  #computing r to z' for each study
  var.z <- 1/(n-3) 
  # UPDATE: ADDED FOLLOWING TO EACH FUNCTION (FOR P-VAL, CI, U3, CLES, CLIFF'S DELTA, NNT)
  alpha <- (100 - level)/100
  crit <- qt(alpha/2, df, lower.tail = FALSE)
  #crit <- qnorm(alpha)
  zval.d <- d/sqrt(var.d)
  pval.d <- 2 * pt(abs(zval.d), df , lower.tail = FALSE)
  lower.d <- d - crit * sqrt(var.d)
  upper.d <- d + crit * sqrt(var.d)
  U3.d <- pnorm(d)*100
  cl.d<- (pnorm((d)/sqrt(2)))*100
  cliffs.d <- 2 * pnorm(d/sqrt(2)) - 1
  zval.g <- g/sqrt(var.g)
  pval.g <- 2 * pt(abs(zval.g), df , lower.tail = FALSE)
  lower.g <- g - crit * sqrt(var.g)
  upper.g <- g + crit * sqrt(var.g)
  U3.g <- pnorm(g)*100
  cl.g <- (pnorm((g)/sqrt(2)))*100
  
  zval.z <- z/sqrt(var.z)
  pval.z <- 2 * pt(abs(zval.z), df , lower.tail = FALSE)
  lower.z <- z - crit * sqrt(var.z)
  upper.z <- z + crit * sqrt(var.z)
  
  pval.r <- pval.z
  lower.r <- ztor(lower.z)
  upper.r <- ztor(upper.z)
  
  zval.lor <- lor/sqrt(var.lor)
  pval.lor <- 2 * pt(abs(zval.lor), df , lower.tail = FALSE)
  lower.lor <- lor - crit * sqrt(var.lor)
  upper.lor <- lor + crit * sqrt(var.lor)
  
  nnt <- 1/(pnorm(d- qnorm(1-cer))-cer)
  
  if (!is.null(data)) {
#     cat("\n")
#     message("    EFFECT SIZE CALCULATION (FOR VECTOR INPUT)")
#     cat("\n")
    
    out <- round(data.frame(id=id, N.total=n, n.1=n.1, n.2=n.2, d=d, var.d=var.d, l.d=lower.d, u.d=upper.d, U3.d=U3.d,
                            cl.d=cl.d, cliffs.d=cliffs.d, pval.d=pval.d, g=g, var.g=var.g,l.g=lower.g, 
                            u.g=upper.g, U3.g=U3.g, cl.g=cl.g, pval.g=pval.g, r=r, var.r=var.r, l.r=lower.r,
                            u.r=upper.r, pval.r=pval.r, fisher.z=z, var.z=var.z, l.z=lower.z, u.z=upper.z,
                            OR=exp(lor), l.or=exp(lower.lor), u.or=exp(upper.lor), pval.or=pval.lor, 
                            lOR=lor, l.lor=lower.lor, u.lor=upper.lor, pval.lor=pval.lor,
                            NNT=nnt), dig)
    return(out)
    
  }
  else{
#     cat("\n")
#     message("    EFFECT SIZE CALCULATION (FOR SINGLE INPUT)")
#     cat("\n")
    if (verbose) {
    cat(#"\n", "\n", "   EFFECT SIZE CALCULATION (FOR SINGLE INPUT)", "\n","\n",
      "Mean Differences ES:","\n", "\n", 
      "d [",level,"%CI] =", round(d, dig),
      "[",round(lower.d, dig),",",round(upper.d,dig),"]","\n",
      "\ var(d) =", round(var.d,dig), "\n",
      "\ p-value(d) =", round(pval.d,dig), "\n",
      "\ U3(d) =", round(U3.d, dig),"%", "\n", 
      "\ CLES(d) =", round(cl.d, dig),"%","\n",
      "\ Cliff's Delta =", round(cliffs.d, dig),"\n", "\n",
      "g [",level,"%CI] =", round(g, dig),
      "[",round(lower.g, dig),",",round(upper.g,dig),"]","\n",
      "\ var(g) =", round(var.g,dig), "\n",
      "\ p-value(g) =", round(pval.g,dig), "\n",
      "\ U3(g) =", round(U3.g, dig),"%", "\n", 
      "\ CLES(g) =", round(cl.g, dig),"%","\n", "\n",
      
      "Correlation ES:","\n", "\n", 
      "r [",level,"%CI] =", round(r, dig),
      "[",round(lower.r, dig),",",round(upper.r,dig),"]","\n",
      "\ var(r) =", round(var.r,dig), "\n",
      "\ p-value(r) =", round(pval.r,dig), "\n","\n",
      "z [",level,"%CI] =", round(z, dig),
      "[",round(lower.z, dig),",",round(upper.z,dig),"]","\n", 
      "\ var(z) =", round(var.z,dig), "\n",
      "\ p-value(z) =", round(pval.z,dig), "\n", "\n",
      
      "Odds Ratio ES:","\n", "\n", 
      "OR [",level,"%CI] =", round(exp(lor), dig),
      "[",round(exp(lower.lor), dig),",",round(exp(upper.lor),dig),"]","\n", 
      "\ p-value(OR) =", round(pval.lor,dig), "\n", "\n",
      "Log OR [",level,"%CI] =", round(lor, dig),
      "[",round(lower.lor, dig),",",round(upper.lor,dig),"]","\n",
      "\ var(lOR) =", round(var.lor,dig), "\n",
      "\ p-value(Log OR) =", round(pval.lor,dig), "\n", "\n", 
      
      "Other:","\n", "\n", 
      "NNT =", round(nnt, dig),"\n", 
      "Total N =", n
    )
    }
    out <- round(data.frame(N.total = n, n.1 = n.1, 
                            n.2 = n.2, d = d, var.d = var.d, l.d = lower.d, u.d = upper.d, 
                            U3.d = U3.d, cl.d = cl.d, cliffs.d = cliffs.d, pval.d = pval.d, 
                            g = g, var.g = var.g, l.g = lower.g, u.g = upper.g, 
                            U3.g = U3.g, cl.g = cl.g, pval.g = pval.g, r = r, 
                            var.r = var.r, l.r = lower.r, u.r = upper.r, pval.r = pval.r, 
                            fisher.z = z, var.z = var.z, l.z = lower.z, u.z = upper.z, 
                            OR = exp(lor), l.or = exp(lower.lor), u.or = exp(upper.lor), 
                            pval.or = pval.lor, lOR = lor, l.lor = lower.lor, 
                            u.lor = upper.lor, pval.lor = pval.lor, NNT = nnt), dig)
    invisible(out)
  }
} 

# compute or from proportions

propes <- function(p1, p2, n.ab, n.cd,level=95, cer=.2, dig=2, verbose=TRUE, id=NULL, data=NULL) {
  if (!is.null(data)) {
    mf <- match.call()
    args <- match(c("p1", "p2", "n.ab", "n.cd", "level", "dig", "id","data"), names(mf), 
                  0)
    mf <- mf[c(1, args)]
    mf$drop.unused.levels <- TRUE
    mf.p1 <- mf[[match("p1", names(mf))]]
    p1 <- eval(mf.p1, data, enclos = sys.frame(sys.parent()))
    mf.p2 <- mf[[match("p2", names(mf))]]
    p2 <- eval(mf.p2, data, enclos = sys.frame(sys.parent()))
    mf.n.ab <- mf[[match("n.ab", names(mf))]]
    n.ab <- eval(mf.n.ab, data, enclos = sys.frame(sys.parent()))
    mf.n.cd <- mf[[match("n.cd", names(mf))]]
    n.cd <- eval(mf.n.cd, data, enclos = sys.frame(sys.parent()))
    mf.id <- mf[[match("id", names(mf))]]
    id <- eval(mf.id, data, enclos = sys.frame(sys.parent()))
  }
  or <-(p1*(1-p2))/(p2*(1-p1))
  lor <- log(or)
  var.lor <- 1/(n.ab*p1*(1-p1))+1/(n.cd*p2*(1-p2))
  d <- lor*sqrt(3)/pi
  var.d <- 3*var.lor/pi^2
  df<- (n.ab+n.cd)-2 
  j<-1-(3/(4*df-1))
  g<-j*d
  var.g<-j^2*var.d
  a <- ((n.ab + n.cd)^2)/(n.ab*n.cd)  # not sure if this is appropriate*
  r <- d/sqrt((d^2) + a)
  var.r <- (a^2*var.d)/(d^2 + a)^3
  lor <- pi*d/sqrt(3)
  var.lor <- pi^2*var.d/3
  n <- n.ab+n.cd
  z <-  0.5*log((1 + r)/(1-r))  #computing r to z' for each study
  var.z <- 1/(n-3) 
  # UPDATE: ADDED FOLLOWING TO EACH FUNCTION (FOR P-VAL, CI, U3, CLES, CLIFF'S DELTA, NNT)
  alpha <- (100 - level)/100
  crit <- qt(alpha/2, df, lower.tail = FALSE)
  #crit <- qnorm(alpha)
  zval.d <- d/sqrt(var.d)
  pval.d <- 2 * pt(abs(zval.d), df , lower.tail = FALSE)
  lower.d <- d - crit * sqrt(var.d)
  upper.d <- d + crit * sqrt(var.d)
  U3.d <- pnorm(d)*100
  cl.d<- (pnorm((d)/sqrt(2)))*100
  cliffs.d <- 2 * pnorm(d/sqrt(2)) - 1
  zval.g <- g/sqrt(var.g)
  pval.g <- 2 * pt(abs(zval.g), df , lower.tail = FALSE)
  lower.g <- g - crit * sqrt(var.g)
  upper.g <- g + crit * sqrt(var.g)
  U3.g <- pnorm(g)*100
  cl.g <- (pnorm((g)/sqrt(2)))*100
  
  zval.z <- z/sqrt(var.z)
  pval.z <- 2 * pt(abs(zval.z), df , lower.tail = FALSE)
  lower.z <- z - crit * sqrt(var.z)
  upper.z <- z + crit * sqrt(var.z)
  
  pval.r <- pval.z
  lower.r <- ztor(lower.z)
  upper.r <- ztor(upper.z)
  
  zval.lor <- lor/sqrt(var.lor)
  pval.lor <- 2 * pt(abs(zval.lor), df , lower.tail = FALSE)
  lower.lor <- lor - crit * sqrt(var.lor)
  upper.lor <- lor + crit * sqrt(var.lor)
  
  nnt <- 1/(pnorm(d- qnorm(1-cer))-cer)
  
  if (!is.null(data)) {
#     cat("\n")
#     message("    EFFECT SIZE CALCULATION (FOR VECTOR INPUT)")
#     cat("\n")
    
    out <- round(data.frame(id=id, N.total=n, d=d, var.d=var.d, l.d=lower.d, u.d=upper.d, U3.d=U3.d,
                            cl.d=cl.d, cliffs.d=cliffs.d, pval.d=pval.d, g=g, var.g=var.g,l.g=lower.g, 
                            u.g=upper.g, U3.g=U3.g, cl.g=cl.g, pval.g=pval.g, r=r, var.r=var.r, l.r=lower.r,
                            u.r=upper.r, pval.r=pval.r, fisher.z=z, var.z=var.z, l.z=lower.z, u.z=upper.z,
                            OR=exp(lor), l.or=exp(lower.lor), u.or=exp(upper.lor), pval.or=pval.lor, 
                            lOR=lor, l.lor=lower.lor, u.lor=upper.lor, pval.lor=pval.lor,
                            NNT=nnt), dig)
    return(out)
    
  }
  else{
#     cat("\n")
#     message("    EFFECT SIZE CALCULATION (FOR SINGLE INPUT)")
#     cat("\n")
    if (verbose) {
    cat(#"\n", "\n", "   EFFECT SIZE CALCULATION (FOR SINGLE INPUT)", "\n","\n",
      "Mean Differences ES:","\n", "\n", 
      "d [",level,"%CI] =", round(d, dig),
      "[",round(lower.d, dig),",",round(upper.d,dig),"]","\n",
      "\ var(d) =", round(var.d,dig), "\n",
      "\ p-value(d) =", round(pval.d,dig), "\n",
      "\ U3(d) =", round(U3.d, dig),"%", "\n", 
      "\ CLES(d) =", round(cl.d, dig),"%","\n",
      "\ Cliff's Delta =", round(cliffs.d, dig),"\n", "\n",
      "g [",level,"%CI] =", round(g, dig),
      "[",round(lower.g, dig),",",round(upper.g,dig),"]","\n",
      "\ var(g) =", round(var.g,dig), "\n",
      "\ p-value(g) =", round(pval.g,dig), "\n",
      "\ U3(g) =", round(U3.g, dig),"%", "\n", 
      "\ CLES(g) =", round(cl.g, dig),"%","\n", "\n",
      
      "Correlation ES:","\n", "\n", 
      "r [",level,"%CI] =", round(r, dig),
      "[",round(lower.r, dig),",",round(upper.r,dig),"]","\n",
      "\ var(r) =", round(var.r,dig), "\n",
      "\ p-value(r) =", round(pval.r,dig), "\n","\n",
      "z [",level,"%CI] =", round(z, dig),
      "[",round(lower.z, dig),",",round(upper.z,dig),"]","\n", 
      "\ var(z) =", round(var.z,dig), "\n",
      "\ p-value(z) =", round(pval.z,dig), "\n", "\n",
      
      "Odds Ratio ES:","\n", "\n", 
      "OR [",level,"%CI] =", round(exp(lor), dig),
      "[",round(exp(lower.lor), dig),",",round(exp(upper.lor),dig),"]","\n", 
      "\ p-value(OR) =", round(pval.lor,dig), "\n", "\n",
      "Log OR [",level,"%CI] =", round(lor, dig),
      "[",round(lower.lor, dig),",",round(upper.lor,dig),"]","\n",
      "\ var(lOR) =", round(var.lor,dig), "\n",
      "\ p-value(Log OR) =", round(pval.lor,dig), "\n", "\n", 
      
      "Other:","\n", "\n", 
      "NNT =", round(nnt, dig),"\n", 
      "Total N =", n
    )
    }
    out <- round(data.frame(N.total = n, d = d, var.d = var.d, l.d = lower.d, u.d = upper.d, 
                            U3.d = U3.d, cl.d = cl.d, cliffs.d = cliffs.d, pval.d = pval.d, 
                            g = g, var.g = var.g, l.g = lower.g, u.g = upper.g, 
                            U3.g = U3.g, cl.g = cl.g, pval.g = pval.g, r = r, 
                            var.r = var.r, l.r = lower.r, u.r = upper.r, pval.r = pval.r, 
                            fisher.z = z, var.z = var.z, l.z = lower.z, u.z = upper.z, 
                            OR = exp(lor), l.or = exp(lower.lor), u.or = exp(upper.lor), 
                            pval.or = pval.lor, lOR = lor, l.lor = lower.lor, 
                            u.lor = upper.lor, pval.lor = pval.lor, NNT = nnt), dig)
    invisible(out)
  }
}


# Odds Ratio to es: if have info for 'failure' in both conditions 
# (B = # tmt failure; D = # non-tmt failure) and the sample size
# for each group (n.1 & n.0 respectively):

failes <- function(B, D, n.1, n.0,level=95, cer=.2, dig=2, verbose=TRUE, id=NULL, data=NULL) {
  if (!is.null(data)) {
    mf <- match.call()
    args <- match(c("B", "D", "n.1", "n.0", "level", "dig", "id","data"), names(mf), 
                  0)
    mf <- mf[c(1, args)]
    mf$drop.unused.levels <- TRUE
    mf.B <- mf[[match("B", names(mf))]]
    B <- eval(mf.B, data, enclos = sys.frame(sys.parent()))
    mf.D <- mf[[match("D", names(mf))]]
    D <- eval(mf.D, data, enclos = sys.frame(sys.parent()))
    mf.n.1 <- mf[[match("n.1", names(mf))]]
    n.1 <- eval(mf.n.1, data, enclos = sys.frame(sys.parent()))
    mf.n.0 <- mf[[match("n.0", names(mf))]]
    n.0 <- eval(mf.n.0, data, enclos = sys.frame(sys.parent()))
    mf.id <- mf[[match("id", names(mf))]]
    id <- eval(mf.id, data, enclos = sys.frame(sys.parent()))
  }
  A <- n.1 - B  # tmt success
  B <- B        # tmt failure
  C <- n.0 - D  # non-tmt success
  D <- D        # non-tmt failure
  p1 <- A/n.1   # proportion 1 
  p2 <- C/n.0   # proportion 2
  n.ab <-  A+B  # n of A+B
  n.cd <-  C+D  # n of C+D        
  or <- (p1 * (1 - p2))/(p2 * (1 - p1))  # odds ratio
  lor <- log(or)  # log odds ratio
  var.lor <-  1/A + 1/B + 1/C + 1/D  # variance of log odds ratio
  #var.lor <- 1/(n.ab*p1*(1-p1))+1/(n.cd*p2*(1-p2))
  d <- lor * sqrt(3)/pi  # conversion to d
  var.d <- 3 * var.lor/pi^2  # variance of d
  df<- (n.1 + n.0)-2 
  j<-1-(3/(4*df-1))
  g<-j*d
  var.g<-j^2*var.d
  a <- ((n.1 + n.0)^2)/(n.1*n.0)  # not sure if this is appropriate*
  r <- d/sqrt((d^2) + a)
  var.r <- (a^2*var.d)/(d^2 + a)^3
  lor <- pi*d/sqrt(3)
  var.lor <- pi^2*var.d/3
  n <- n.1 + n.0
  z <-  0.5*log((1 + r)/(1-r))  #computing r to z' for each study
  var.z <- 1/(n-3) 
  # UPDATE: ADDED FOLLOWING TO EACH FUNCTION (FOR P-VAL, CI, U3, CLES, CLIFF'S DELTA, NNT)
  alpha <- (100 - level)/100
  crit <- qt(alpha/2, df, lower.tail = FALSE)
  #crit <- qnorm(alpha)
  zval.d <- d/sqrt(var.d)
  pval.d <- 2 * pt(abs(zval.d), df , lower.tail = FALSE)
  lower.d <- d - crit * sqrt(var.d)
  upper.d <- d + crit * sqrt(var.d)
  U3.d <- pnorm(d)*100
  cl.d<- (pnorm((d)/sqrt(2)))*100
  cliffs.d <- 2 * pnorm(d/sqrt(2)) - 1
  zval.g <- g/sqrt(var.g)
  pval.g <- 2 * pt(abs(zval.g), df , lower.tail = FALSE)
  lower.g <- g - crit * sqrt(var.g)
  upper.g <- g + crit * sqrt(var.g)
  U3.g <- pnorm(g)*100
  cl.g <- (pnorm((g)/sqrt(2)))*100
  
  zval.z <- z/sqrt(var.z)
  pval.z <- 2 * pt(abs(zval.z), df , lower.tail = FALSE)
  lower.z <- z - crit * sqrt(var.z)
  upper.z <- z + crit * sqrt(var.z)
  
  pval.r <- pval.z
  lower.r <- ztor(lower.z)
  upper.r <- ztor(upper.z)
  
  zval.lor <- lor/sqrt(var.lor)
  pval.lor <- 2 * pt(abs(zval.lor), df , lower.tail = FALSE)
  lower.lor <- lor - crit * sqrt(var.lor)
  upper.lor <- lor + crit * sqrt(var.lor)
  
  nnt <- 1/(pnorm(d- qnorm(1-cer))-cer)
  
  if (!is.null(data)) {
#     cat("\n")
#     message("    EFFECT SIZE CALCULATION (FOR VECTOR INPUT)")
#     cat("\n")
    
    out <- round(data.frame(id=id, N.total=n, n.1=n.1, n.2=n.0, d=d, var.d=var.d, l.d=lower.d, u.d=upper.d, U3.d=U3.d,
                            cl.d=cl.d, cliffs.d=cliffs.d, pval.d=pval.d, g=g, var.g=var.g,l.g=lower.g, 
                            u.g=upper.g, U3.g=U3.g, cl.g=cl.g, pval.g=pval.g, r=r, var.r=var.r, l.r=lower.r,
                            u.r=upper.r, pval.r=pval.r, fisher.z=z, var.z=var.z, l.z=lower.z, u.z=upper.z,
                            OR=exp(lor), l.or=exp(lower.lor), u.or=exp(upper.lor), pval.or=pval.lor, 
                            lOR=lor, l.lor=lower.lor, u.lor=upper.lor, pval.lor=pval.lor,
                            NNT=nnt), dig)
    return(out)
    
  }
  else{
#     cat("\n")
#     message("    EFFECT SIZE CALCULATION (FOR SINGLE INPUT)")
#     cat("\n")
    if (verbose) {
    cat(#"\n", "\n", "   EFFECT SIZE CALCULATION (FOR SINGLE INPUT)", "\n","\n",
      "Mean Differences ES:","\n", "\n", 
      "d [",level,"%CI] =", round(d, dig),
      "[",round(lower.d, dig),",",round(upper.d,dig),"]","\n",
      "\ var(d) =", round(var.d,dig), "\n",
      "\ p-value(d) =", round(pval.d,dig), "\n",
      "\ U3(d) =", round(U3.d, dig),"%", "\n", 
      "\ CLES(d) =", round(cl.d, dig),"%","\n",
      "\ Cliff's Delta =", round(cliffs.d, dig),"\n", "\n",
      "g [",level,"%CI] =", round(g, dig),
      "[",round(lower.g, dig),",",round(upper.g,dig),"]","\n",
      "\ var(g) =", round(var.g,dig), "\n",
      "\ p-value(g) =", round(pval.g,dig), "\n",
      "\ U3(g) =", round(U3.g, dig),"%", "\n", 
      "\ CLES(g) =", round(cl.g, dig),"%","\n", "\n",
      
      "Correlation ES:","\n", "\n", 
      "r [",level,"%CI] =", round(r, dig),
      "[",round(lower.r, dig),",",round(upper.r,dig),"]","\n",
      "\ var(r) =", round(var.r,dig), "\n",
      "\ p-value(r) =", round(pval.r,dig), "\n","\n",
      "z [",level,"%CI] =", round(z, dig),
      "[",round(lower.z, dig),",",round(upper.z,dig),"]","\n", 
      "\ var(z) =", round(var.z,dig), "\n",
      "\ p-value(z) =", round(pval.z,dig), "\n", "\n",
      
      "Odds Ratio ES:","\n", "\n", 
      "OR [",level,"%CI] =", round(exp(lor), dig),
      "[",round(exp(lower.lor), dig),",",round(exp(upper.lor),dig),"]","\n", 
      "\ p-value(OR) =", round(pval.lor,dig), "\n", "\n",
      "Log OR [",level,"%CI] =", round(lor, dig),
      "[",round(lower.lor, dig),",",round(upper.lor,dig),"]","\n",
      "\ var(lOR) =", round(var.lor,dig), "\n",
      "\ p-value(Log OR) =", round(pval.lor,dig), "\n", "\n", 
      
      "Other:","\n", "\n", 
      "NNT =", round(nnt, dig),"\n", 
      "Total N =", n
    )
    }
    out <- round(data.frame(N.total = n, n.1=n.1, n.2=n.0, 
                            d = d, var.d = var.d, l.d = lower.d, u.d = upper.d, 
                            U3.d = U3.d, cl.d = cl.d, cliffs.d = cliffs.d, pval.d = pval.d, 
                            g = g, var.g = var.g, l.g = lower.g, u.g = upper.g, 
                            U3.g = U3.g, cl.g = cl.g, pval.g = pval.g, r = r, 
                            var.r = var.r, l.r = lower.r, u.r = upper.r, pval.r = pval.r, 
                            fisher.z = z, var.z = var.z, l.z = lower.z, u.z = upper.z, 
                            OR = exp(lor), l.or = exp(lower.lor), u.or = exp(upper.lor), 
                            pval.or = pval.lor, lOR = lor, l.lor = lower.lor, 
                            u.lor = upper.lor, pval.lor = pval.lor, NNT = nnt), dig)
    invisible(out)
  }
}

# Converting Chi-squared statistic with 1 df to es

chies <- function(chi.sq,  n,level=95, cer=.2, dig=2, verbose=TRUE, id=NULL, data=NULL) {
  if (!is.null(data)) {
    mf <- match.call()
    args <- match(c("chi.sq", "n",  "level", "dig", "id","data"), names(mf), 
                  0)
    mf <- mf[c(1, args)]
    mf$drop.unused.levels <- TRUE
    mf.chi.sq <- mf[[match("chi.sq", names(mf))]]
    chi.sq <- eval(mf.chi.sq, data, enclos = sys.frame(sys.parent()))
    mf.n <- mf[[match("n", names(mf))]]
    n <- eval(mf.n, data, enclos = sys.frame(sys.parent()))
    mf.id <- mf[[match("id", names(mf))]]
    id <- eval(mf.id, data, enclos = sys.frame(sys.parent()))
  }
  r <- sqrt(chi.sq/n)
  var.r <-(1-r^2)^2/(n-1) 
  d<-2*r*sqrt((n-1)/(n*(1-r^2)))*abs(r)/r 
  var.d <- 4*var.r/(1-r^2)^3
  df<- (n)-2 
  j<-1-(3/(4*df-1))
  g<-j*d
  var.g<-j^2*var.d
  lor <- pi*d/sqrt(3)
  var.lor <- pi^2*var.d/3
  z <-  0.5*log((1 + r)/(1-r))  
  var.z <- 1/(n-3) 
  # UPDATE: ADDED FOLLOWING TO EACH FUNCTION (FOR P-VAL, CI, U3, CLES, CLIFF'S DELTA, NNT)
  alpha <- (100 - level)/100
  crit <- qt(alpha/2, df, lower.tail = FALSE)
  #crit <- qnorm(alpha)
  zval.d <- d/sqrt(var.d)
  pval.d <- 2 * pt(abs(zval.d), df , lower.tail = FALSE)
  lower.d <- d - crit * sqrt(var.d)
  upper.d <- d + crit * sqrt(var.d)
  U3.d <- pnorm(d)*100
  cl.d<- (pnorm((d)/sqrt(2)))*100
  cliffs.d <- 2 * pnorm(d/sqrt(2)) - 1
  zval.g <- g/sqrt(var.g)
  pval.g <- 2 * pt(abs(zval.g), df , lower.tail = FALSE)
  lower.g <- g - crit * sqrt(var.g)
  upper.g <- g + crit * sqrt(var.g)
  U3.g <- pnorm(g)*100
  cl.g <- (pnorm((g)/sqrt(2)))*100
  
  zval.z <- z/sqrt(var.z)
  pval.z <- 2 * pt(abs(zval.z), df , lower.tail = FALSE)
  lower.z <- z - crit * sqrt(var.z)
  upper.z <- z + crit * sqrt(var.z)
  
  pval.r <- pval.z
  lower.r <- ztor(lower.z)
  upper.r <- ztor(upper.z)
  
  zval.lor <- lor/sqrt(var.lor)
  pval.lor <- 2 * pt(abs(zval.lor), df , lower.tail = FALSE)
  lower.lor <- lor - crit * sqrt(var.lor)
  upper.lor <- lor + crit * sqrt(var.lor)
  
  nnt <- 1/(pnorm(d- qnorm(1-cer))-cer)
  
  if (!is.null(data)) {
#     cat("\n")
#     message("    EFFECT SIZE CALCULATION (FOR VECTOR INPUT)")
#     cat("\n")
    
    out <- round(data.frame(id=id, N.total=n, d=d, var.d=var.d, l.d=lower.d, u.d=upper.d, U3.d=U3.d,
                            cl.d=cl.d, cliffs.d=cliffs.d, pval.d=pval.d, g=g, var.g=var.g,l.g=lower.g, 
                            u.g=upper.g, U3.g=U3.g, cl.g=cl.g, pval.g=pval.g, r=r, var.r=var.r, l.r=lower.r,
                            u.r=upper.r, pval.r=pval.r, fisher.z=z, var.z=var.z, l.z=lower.z, u.z=upper.z,
                            OR=exp(lor), l.or=exp(lower.lor), u.or=exp(upper.lor), pval.or=pval.lor, 
                            lOR=lor, l.lor=lower.lor, u.lor=upper.lor, pval.lor=pval.lor,
                            NNT=nnt), dig)
    return(out)
    
  }
  else{
#     cat("\n")
#     message("    EFFECT SIZE CALCULATION (FOR SINGLE INPUT)")
#     cat("\n")
    if (verbose) {
    cat(#"\n", "\n", "   EFFECT SIZE CALCULATION (FOR SINGLE INPUT)", "\n","\n",
      "Mean Differences ES:","\n", "\n", 
      "d [",level,"%CI] =", round(d, dig),
      "[",round(lower.d, dig),",",round(upper.d,dig),"]","\n",
      "\ var(d) =", round(var.d,dig), "\n",
      "\ p-value(d) =", round(pval.d,dig), "\n",
      "\ U3(d) =", round(U3.d, dig),"%", "\n", 
      "\ CLES(d) =", round(cl.d, dig),"%","\n",
      "\ Cliff's Delta =", round(cliffs.d, dig),"\n", "\n",
      "g [",level,"%CI] =", round(g, dig),
      "[",round(lower.g, dig),",",round(upper.g,dig),"]","\n",
      "\ var(g) =", round(var.g,dig), "\n",
      "\ p-value(g) =", round(pval.g,dig), "\n",
      "\ U3(g) =", round(U3.g, dig),"%", "\n", 
      "\ CLES(g) =", round(cl.g, dig),"%","\n", "\n",
      
      "Correlation ES:","\n", "\n", 
      "r [",level,"%CI] =", round(r, dig),
      "[",round(lower.r, dig),",",round(upper.r,dig),"]","\n",
      "\ var(r) =", round(var.r,dig), "\n",
      "\ p-value(r) =", round(pval.r,dig), "\n","\n",
      "z [",level,"%CI] =", round(z, dig),
      "[",round(lower.z, dig),",",round(upper.z,dig),"]","\n", 
      "\ var(z) =", round(var.z,dig), "\n",
      "\ p-value(z) =", round(pval.z,dig), "\n", "\n",
      
      "Odds Ratio ES:","\n", "\n", 
      "OR [",level,"%CI] =", round(exp(lor), dig),
      "[",round(exp(lower.lor), dig),",",round(exp(upper.lor),dig),"]","\n", 
      "\ p-value(OR) =", round(pval.lor,dig), "\n", "\n",
      "Log OR [",level,"%CI] =", round(lor, dig),
      "[",round(lower.lor, dig),",",round(upper.lor,dig),"]","\n",
      "\ var(lOR) =", round(var.lor,dig), "\n",
      "\ p-value(Log OR) =", round(pval.lor,dig), "\n", "\n", 
      
      "Other:","\n", "\n", 
      "NNT =", round(nnt, dig),"\n", 
      "Total N =", n
    )
    }
    out <- round(data.frame(N.total = n, d = d, var.d = var.d, l.d = lower.d, u.d = upper.d, 
                            U3.d = U3.d, cl.d = cl.d, cliffs.d = cliffs.d, pval.d = pval.d, 
                            g = g, var.g = var.g, l.g = lower.g, u.g = upper.g, 
                            U3.g = U3.g, cl.g = cl.g, pval.g = pval.g, r = r, 
                            var.r = var.r, l.r = lower.r, u.r = upper.r, pval.r = pval.r, 
                            fisher.z = z, var.z = var.z, l.z = lower.z, u.z = upper.z, 
                            OR = exp(lor), l.or = exp(lower.lor), u.or = exp(upper.lor), 
                            pval.or = pval.lor, lOR = lor, l.lor = lower.lor, 
                            u.lor = upper.lor, pval.lor = pval.lor, NNT = nnt), dig)
    invisible(out)
  }
}





