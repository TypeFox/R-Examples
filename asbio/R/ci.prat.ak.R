ci.prat.ak <- function(y1, n1 , pi2 = NULL, method = "ac", conf = .95, bonf = FALSE, bootCI.method = "perc", R = 1000, sigma.t = NULL, r = length(y1), gamma.hyper = 1, beta.hyper = 1){

rnm <- data.frame(CImethod = c("norm", "basic", "perc", "BCa", "student"),
row = 1:5, new = c("normal","basic","percentile", "BCa","studentized"))
rn  <- rnm[,2][rnm[,1] == bootCI.method]
name <- matrix(rnm[,3])[rn]

indices <- c("ac", "bayes", "boot", "wald", "fixed-log", "noether-fixed")
    method <- match.arg(method, indices)

   if (length(y1) != length(n1)) stop("y1 and n1 must have equal length")
   if(any(n1 - y1 < 0)) stop("n1 cannot be less than y1")
   
alpha <- 1 - conf
oconf <- conf
conf <- ifelse(bonf == FALSE, conf, 1 - alpha/r)
z.star <- qnorm(1 - (1 - conf)/2)

ci.prat.ak1 <- function(y1, n1, pi2 = NULL, method = "ac", conf = .95, bonf = FALSE, sigma.t = sigma.t, bootCI.method = "perc", R = 1000){

#--------------- Bootstrap ---------------#

if(method == "boot"){
  pih1 <- y1/n1
  rat <- pih1/pi2
  if((y1 == 0 & pi2 == 0)|(y1 == 0 & pi2 != 0)|(y1 != 0 & pi2 == 0)|(y1 == n1)){
      if(y1 == 0 & pi2 == 0) {CIL <- 0;  CIU <- Inf; rat = 0; varhat = NA}
      if(y1 == 0 & pi2 != 0){
        if(pi2 == 1)stop("pi2 = 1 but y1 = 0?")    
        if(pi2 != 1){
          CIL <- 0
          yn1 <- 0.5
          pihn1 <- yn1/n1
          nrat <- (yn1/n1)/pi2
        varhat <- (pihn1 * (1 - pihn1))/(n1 * pi2^2)
        CIU <- nrat + (z.star * sqrt(varhat))
        }
      }
      if(pi2 == 0 & y1 != 0)stop("pi2 = 0 but y1 != 0?") 
      if(pi2 != 0 & y1 == n1){
        if(pi2 == 1){
          varhat <- 0
          CIL <- 1
          CIU <- 1
        }
        if(pi2 != 1){
          yn1 <- n1 - 0.5
          pihn1 <- yn1/n1
          nrat <- pihn1/pi2
          varhat <- (pihn1 * (1 - pihn1))/(n1 * pi2^2)
          CIL <- nrat - (z.star * sqrt(varhat))
          CIU <- nrat + (z.star * sqrt(varhat))
         }
       } 
    } 
  else{
    num.data <- c(rep(1, y1), rep(0, (n1 - y1)))
    brat <- function(x){(sum(x)/n1)/pi2}
    b <- bootstrap(num.data, brat, R = R)
    c <- ci.boot(b, method = bootCI.method, sigma.t = sigma.t, conf = conf)$res[rn,]
    CIU <- c[2]
    CIL <- c[1]
    varhat <- b$res[4]
  }
  CI <- c(rat, CIL, CIU)
}

#--------------- Manly/Wald ---------------#

if(method == "wald"){
  pih1 <- y1/n1
  rat <- pih1/pi2
  if((y1 == 0 & pi2 == 0)|(y1 == 0 & pi2 != 0)|(y1 != 0 & pi2 == 0)|(y1 == n1)){
    if(y1 == 0 & pi2 == 0) {CIL <- 0;  CIU <- Inf; rat = 0; varhat = NA}
    if(y1 == 0 & pi2 != 0){
      if(pi2 == 1)stop("pi2 = 1 but y1 = 0?")    
      if(pi2 != 1){
        CIL <- 0
        yn1 <- 0.5
        pihn1 <- yn1/n1
        nrat <- (yn1/n1)/pi2
      varhat <- (pihn1 * (1 - pihn1))/(n1 * pi2^2)
      CIU <- nrat + (z.star * sqrt(varhat))
      }
    }
    if(pi2 == 0 & y1 != 0)stop("pi2 = 0 but y1 != 0?") 
    if(pi2 != 0 & y1 == n1){
      if(pi2 == 1){
        varhat <- 0
        CIL <- 1
        CIU <- 1
      }
      if(pi2 != 1){
        yn1 <- n1 - 0.5
        pihn1 <- yn1/n1
        nrat <- pihn1/pi2
        varhat <- (pihn1 * (1 - pihn1))/(n1 * pi2^2)
        CIL <- nrat - (z.star * sqrt(varhat))
        CIU <- nrat + (z.star * sqrt(varhat))
       }
     } 
  }
  else{
  varhat <- (pih1 * (1 - pih1))/(n1 * pi2^2)
CIL <- rat - (z.star*sqrt(varhat))
CIU <- rat + (z.star*sqrt(varhat))
}
CIL <- ifelse(CIL < 0, 0, CIL)
CI <- c(rat, CIL, CIU)
}

#-------------------noether-fixed----------------#

if(method == "noether-fixed"){
  pih1 <- y1/n1
  rat <- pih1/pi2
  if((y1 == 0 & pi2 == 0)|(y1 == 0 & pi2 != 0)|(y1 != 0 & pi2 == 0)|(y1 == n1)){
    if(y1 == 0 & pi2 == 0) {CIL <- 0;  CIU <- Inf; rat = 0; varhat = NA}
    if(y1 == 0 & pi2 != 0){
      if(pi2 == 1)stop("pi2 = 1 but y1 = 0?")    
      if(pi2 != 1){
        CIL <- 0
        yn1 <- 0.5
        pihn1 <- yn1/n1
        nrat <- (yn1/n1)/pi2
      varhat <- (1 - pihn1)/(n1 * pihn1)
      CIU <- ((pihn1/pi2)/(1 + (z.star^2/n1))) * (1 + ((z.star^2/(2*yn1)))) + z.star * sqrt(varhat + (z.star^2/(4 * yn1^2)))
      }
    }
    if(pi2 == 0 & y1 != 0)stop("pi2 = 0 but y1 != 0?") 
    if(pi2 != 0 & y1 == n1){
      if(pi2 == 1){
        varhat <- 0
        CIL <- 1
        CIU <- 1
      }
      if(pi2 != 1){
        yn1 <- n1 - 0.5
        pihn1 <- yn1/n1
        nrat <- pihn1/pi2
        varhat <- (1 - pihn1)/(n1 * pihn1)
        CIL <- ((pihn1/pi2)/(1 + (z.star^2/n1))) * ((1 + ((z.star^2/(2*yn1)))) - z.star * sqrt(varhat + (z.star^2/(4 * yn1^2))))
        CIU <- ((pihn1/pi2)/(1 + (z.star^2/n1))) * ((1 + ((z.star^2/(2*yn1)))) + z.star * sqrt(varhat + (z.star^2/(4 * yn1^2))))
      }
    } 
  }
  else{
    varhat <- (1 - pih1)/(n1 * pih1)
    CIL <- ((pih1/pi2)/(1 + (z.star^2/n1))) * ((1 + ((z.star^2/(2*y1)))) - z.star * sqrt(varhat + (z.star^2/(4 * y1^2))))
    CIU <- ((pih1/pi2)/(1 + (z.star^2/n1))) * ((1 + ((z.star^2/(2*y1)))) + z.star * sqrt(varhat + (z.star^2/(4 * y1^2))))
    
  }
CIL <- ifelse(CIL < 0, 0, CIL)
  CI <- c(rat, CIL, CIU)
}

#------------------fixed-log--------------------#

if(method == "fixed-log"){
  pih1 <- y1/n1
  rat <- pih1/pi2
  if((y1 == 0 & pi2 == 0)|(y1 == 0 & pi2 != 0)|(y1 != 0 & pi2 == 0)|(y1 == n1)){
    if(y1 == 0 & pi2 == 0) {CIL <- 0;  CIU <- Inf; rat = 0; varhat = NA}
    if(y1 == 0 & pi2 != 0){
      if(pi2 == 1)stop("pi2 = 1 but y1 = 0?")    
      if(pi2 != 1){
        CIL <- 0
        yn1 <- 0.5
        pihn1 <- yn1/n1
        nrat <- (yn1/n1)/pi2
      varhat <- (1 - pihn1)/(n1 * pihn1)
      CIU <- nrat * exp(z.star * sqrt(varhat))
      }
    }
    if(pi2 == 0 & y1 != 0)stop("pi2 = 0 but y1 != 0?") 
    if(pi2 != 0 & y1 == n1){
      if(pi2 == 1){
        varhat <- 0
        CIL <- 1
        CIU <- 1
      }
      if(pi2 != 1){
        yn1 <- n1 - 0.5
        pihn1 <- yn1/n1
        nrat <- pihn1/pi2
        varhat <- (1 - pihn1)/(n1 * pihn1)
        CIL <- nrat * exp(-z.star * sqrt(varhat))
      CIU <- nrat * exp(z.star * sqrt(varhat))
      }
    } 
  }
  else{
    varhat <- (1 - pih1)/(n1 * pih1)
  CIL <- rat * exp(-z.star * sqrt(varhat))
  CIU <- rat * exp(z.star * sqrt(varhat))
}
CIL <- ifelse(CIL < 0, 0, CIL)
CI <- c(rat, CIL, CIU)
}

#------------------Agresti-Coull--------------------#

if(method == "ac"){
  if(pi2 == 0 & y1 != 0)stop("pi2 = 0 but y1 != 0?")
  if(pi2 == 1 & y1 == 0)stop("pi2 = 1 but y1 = 0?")
  
  y1 <- y1 + 2; n1 <- n1 + 4
  
  pih1 <- y1/n1
  rat <- pih1/pi2
  varhat <- pih1 * (1 - pih1)/((n1) * pi2^2)
  CIL <- rat - (z.star*sqrt(varhat))
  CIU <- rat + (z.star*sqrt(varhat))
  CIL <- ifelse(CIL < 0, 0, CIL)
  CI <- c(rat, CIL, CIU)
}


#------------------- Bayes ------------------------#
 
if(method == "bayes"){
  alpha <- 1 - conf
  pih1 <- y1/n1 
  rat <- pih1/pi2
  a <- gamma.hyper + y1
  b <- beta.hyper + n1 - y1
  varhat <- a * b/((a + b)^2*(a + b + 1))
  CIL <- qbeta(alpha/2, a, b)/pi2 
  CIU <- qbeta(1-alpha/2, a, b)/pi2
  CI <- c(rat, CIL, CIU)
}

#--------------------------------------------------#

  res <- list(CI = CI, varhat = varhat)
  res
}
     
  CI <- matrix(ncol = 3, nrow = length(y1))
  vh <- rep(NA, length(y1))
  for (i in 1 : length(y1)) {
  temp <- ci.prat.ak1(y1 = y1[i], n1 = n1[i], pi2 = pi2[i], conf = conf, method = method, bonf = bonf, bootCI.method = bootCI.method, sigma.t=sigma.t, R = R)
  CI[i,] <- temp$CI
  vh[i] <- temp$varhat
  }

CI <- data.frame(CI)
if (length(y1) == 1)row.names(CI) <- ""

head <- paste(paste(as.character(oconf * 100), "%", sep = ""), c("Confidence interval for ratio of binomial proportions"))
head1 <- paste(paste(as.character(oconf * 100), "%", sep = ""), c("Credible interval for ratio of binomial proportions"))

   if (method == "bayes")  head <- paste(head1, "(method=Bayes-beta)", sep ="")
   if (method == "boot")  head <- paste(head, "(method=", name, " bootstrap)", sep ="")
   if (method == "fixed-log")  head <- paste(head, "(method=fixed-log)")
   if (method == "wald") head <- paste(head, "(method=Wald-adjusted)")
   if(method == "ac") head <- paste(head, "(method=AC-adjusted)")
   if(method == "noether-fixed") head <- paste(head, "(method=Noether-fixed)")
   if(bonf == TRUE) head <- paste(head, "\n Bonferroni simultaneous intervals, r = ", bquote(.(r)), 
"\n Marginal confidence = ", bquote(.(conf)), "\n", sep = "")                                                                                                                                              


ends <- c("Estimate", paste(as.character(c((1 - oconf)/2, 1 - ((1 - oconf)/2)) * 100), "%", sep = ""))
res <- list(varhat = vh, ci = CI, ends = ends, head = head)
class(res) <- "ci"
res
}
