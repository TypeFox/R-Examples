ci.prat <- function(y1, n1, y2, n2, conf = 0.95, method = "katz.log", bonf = FALSE, tol = .Machine$double.eps^0.25, R = 1000, r = length(y1)){

x <- y1; m <- n1; y <- y2; n <- n2 

indices <- c("adj.log","bailey","boot","katz.log","koopman","noether","sinh-1")
method <- match.arg(method, indices)


if(any(c(length(m),length(y),length(n))!= length(x))) stop("y1, n1, y2, and n2 vectors must have equal length")

alpha <- 1 - conf
oconf <- conf
conf <- ifelse(bonf == FALSE, conf, 1 - alpha/r)
z.star <- qnorm(1 - (1 - conf)/2)
x2 <- qchisq(conf, 1)
     
ci.prat1 <- function(x, m, y, n, conf = 0.95, method = "katz.log", bonf = FALSE){
if((x > m)|(y > n)) stop("Use correct parameterization for y1, y2, n1, and n2")

#-------------------------- Adj-log ------------------------------#

if(method == "adj.log"){
    if((x == m & y == n)){
       rat <- (x/m)/(y/n); x <- m - 0.5; y <- n - 0.5; nrat <- ((x+0.5)/(m+0.5))/((y+0.5)/(n+0.5)); varhat <- (1/(x+0.5)) - (1/(m+0.5)) + (1/(y+0.5)) - (1/(n+0.5))
       CIL <- nrat * exp(-1 * z.star * sqrt(varhat))
       CIU <- nrat * exp(z.star * sqrt(varhat))
       } else if(x == 0 & y == 0){CIL = 0; CIU = Inf; rat = 0; varhat <- (1/(x+0.5)) - (1/(m+0.5)) + (1/(y+0.5)) - (1/(n+0.5))
       }else{ 
     rat <- (x/m)/(y/n); nrat <- ((x+0.5)/(m+0.5))/((y+0.5)/(n+0.5)); varhat <- (1/(x+0.5)) - (1/(m+0.5)) + (1/(y+0.5)) - (1/(n+0.5))
     CIL <- nrat * exp(-1 * z.star * sqrt(varhat))
     CIU <- nrat * exp(z.star * sqrt(varhat))}
CI <- c(rat, CIL, CIU)
}
                               
#-------------------------------Bailey-----------------------------#

if(method == "bailey"){
   rat <- (x/m)/(y/n)
   varhat <- ifelse((x == m) & (y == n),(1/(m-0.5)) - (1/(m)) + (1/(n-0.5)) - (1/(n)),(1/(x)) - (1/(m)) + (1/(y)) - (1/(n)))

   p.hat1 <- x/m; p.hat2 <- y/n; 
   q.hat1 <- 1 - p.hat1; q.hat2 <- 1 - p.hat2 

   if(x == 0 | y == 0){
      xn <- ifelse(x == 0, 0.5, x)
      yn <- ifelse(y == 0, 0.5, y)
      nrat <- (xn/m)/(yn/n)
      p.hat1 <- xn/m; p.hat2 <- yn/n; 
      q.hat1 <- 1 - p.hat1; q.hat2 <- 1 - p.hat2
      if(xn == m | yn == n){
         xn <- ifelse(xn == m, m - 0.5, xn)
         yn <- ifelse(yn == n, n - 0.5, yn)
         nrat <- (xn/m)/(yn/n)
         p.hat1 <- xn/m; p.hat2 <- yn/n; 
         q.hat1 <- 1 - p.hat1; q.hat2 <- 1 - p.hat2
      }
   }

   if(x == 0 | y == 0){
      if(x == 0 & y == 0){
         rat <- Inf
         CIL <- 0
         CIU <- Inf
      }
      if(x == 0 & y != 0){
         CIL <- 0
         CIU <- nrat * ((1+ z.star * sqrt((q.hat1/xn) + (q.hat2/yn) - (z.star^2 * q.hat1 * q.hat2)/(9 * xn * yn))/3)/((1 - (z.star^2 * q.hat2)/(9 * yn))))^3
      }
      if(y == 0 & x != 0){
         CIU = Inf                                                                                                                               
         CIL <- nrat * ((1- z.star * sqrt((q.hat1/xn) + (q.hat2/yn) - (z.star^2 * q.hat1 * q.hat2)/(9 * xn * yn))/3)/((1 - (z.star^2 * q.hat2)/(9 * yn))))^3 
      } 
   }else if(x == m | y == n){
       xn <- ifelse(x == m, m - 0.5, x)
       yn <- ifelse(y == n, n - 0.5, y)
       nrat <- (xn/m)/(yn/n)
       p.hat1 <- xn/m; p.hat2 <- yn/n; 
       q.hat1 <- 1 - p.hat1; q.hat2 <- 1 - p.hat2
       CIL <- nrat * ((1- z.star * sqrt((q.hat1/xn) + (q.hat2/yn) - (z.star^2 * q.hat1 * q.hat2)/(9 * xn * yn))/3)/((1 - (z.star^2 * q.hat2)/(9 * yn))))^3 
       CIU <- nrat * ((1+ z.star * sqrt((q.hat1/xn) + (q.hat2/yn) - (z.star^2 * q.hat1 * q.hat2)/(9 * xn * yn))/3)/((1 - (z.star^2 * q.hat2)/(9 * yn))))^3 
   }else{
       CIL <- rat * ((1- z.star * sqrt((q.hat1/x) + (q.hat2/y) - (z.star^2 * q.hat1 * q.hat2)/(9 * x * y))/3)/((1 - (z.star^2 * q.hat2)/(9 * y))))^3 
       CIU <- rat * ((1+ z.star * sqrt((q.hat1/x) + (q.hat2/y) - (z.star^2 * q.hat1 * q.hat2)/(9 * x * y))/3)/((1 - (z.star^2 * q.hat2)/(9 * y))))^3 
   }
CI <- c(rat, CIL, CIU)
}

#-------------------------- Boot ------------------------------#

if(method == "boot"){
rat <- (x/m)/(y/n)
if((x == 0 & y == 0)|(x == 0 & y != 0)|(x != 0 & y == 0)){
    if(x == 0 & y == 0) {CIL <- 0;  CIU <- Inf; rat = 0; varhat = NA}
    if(x == 0 & y != 0) {CIL <- 0;  rat <- (x/m)/(y/n); x <- 0.5; nrat <- (x/m)/(y/n)
    varhat <- (1/x) - (1/m) + (1/y) - (1/n)
    CIU <- nrat * exp(z.star * sqrt(varhat))}
    if(x != 0 & y == 0) {CIU <- Inf;  rat <- (x/m)/(y/n); y <- 0.5; nrat <- (x/m)/(y/n)
    varhat <- (1/x) - (1/m) + (1/y) - (1/n)
    CIL <- nrat * exp(-1 * z.star * sqrt(varhat))}
    } else{
num.data <- c(rep(1, x), rep(0, m - x))
den.data <- c(rep(1, y), rep(0, n - y))
nd <- matrix(ncol = R, nrow = m)
dd <- matrix(ncol = R, nrow = n)
brat <- 1:R
	for(i in 1:R){
		nd[,i] <- sample(num.data, m, replace = TRUE)
		dd[,i] <- sample(den.data, n, replace = TRUE)
		brat[i] <- (sum(nd[,i])/m)/(sum(dd[,i])/n)
		}
alpha <- 1 - conf
CIU <- quantile(brat, 1 - alpha/2, na.rm = TRUE)
CIL <- quantile(brat, alpha/2, na.rm = TRUE)
varhat <- var(brat)
}
CI <- c(rat, CIL, CIU)
}	

#-------------------------- Katz-log ------------------------------#

if(method == "katz.log"){
    if((x == 0 & y == 0)|(x == 0 & y != 0)|(x != 0 & y == 0)|(x == m & y == n)){
    if(x == 0 & y == 0) {CIL <- 0;  CIU <- Inf; rat = 0; varhat = NA}
    if(x == 0 & y != 0) {CIL <- 0;  rat <- (x/m)/(y/n); x <- 0.5; nrat <- (x/m)/(y/n)
    varhat <- (1/x) - (1/m) + (1/y) - (1/n)
    CIU <- nrat * exp(z.star * sqrt(varhat))}
    if(x != 0 & y == 0) {CIU <- Inf;  rat <- (x/m)/(y/n); y <- 0.5; nrat <- (x/m)/(y/n)
    varhat <- (1/x) - (1/m) + (1/y) - (1/n)
    CIL <- nrat * exp(-1 * z.star * sqrt(varhat))}
    if(x == m & y == n) {
       rat <- (x/m)/(y/n); x <- m - 0.5; y <- n - 0.5; nrat <- (x/m)/(y/n); varhat <- (1/x) - (1/m) + (1/y) - (1/n); CIL <- nrat * exp(-1 * z.star * sqrt(varhat))
       x <- m - 0.5; y <- n - 0.5; nrat <- (x/m)/(y/n); varhat <- (1/x) - (1/m) + (1/y) - (1/n); CIU <- nrat * exp(z.star * sqrt(varhat))
       }
    } else
    {rat <- (x/m)/(y/n); varhat <- (1/x) - (1/m) + (1/y) - (1/n)
     CIL <- rat * exp(-1 * z.star * sqrt(varhat)) 
     CIU <- rat * exp(z.star * sqrt(varhat))}
CI <- c(rat, CIL, CIU)
}

#-------------------------- Koopman ------------------------------#

if(method == "koopman"){

    if(x == 0 & y == 0) {CIL <- 0;  CIU <- Inf; rat = 0; varhat = NA
    } else {

        a1 = n * (n * (n + m) * x + m * (n + x) * (z.star^2))
        a2 = -n * (n * m * (y + x) + 2 * (n + m) * y *
            x + m * (n + y + 2 * x) * (z.star^2))
        a3 = 2 * n * m * y * (y + x) + (n + m) * (y^2) *
            x + n * m * (y + x) * (z.star^2)
        a4 = -m * (y^2) * (y + x)
        b1 = a2/a1; b2 = a3/a1; b3 = a4/a1
        c1 = b2 - (b1^2)/3;  c2 = b3 - b1 * b2/3 + 2 * (b1^3)/27
        ceta = suppressWarnings(acos(sqrt(27) * c2/(2 * c1 * sqrt(-c1))))
        t1 <- suppressWarnings(-2 * sqrt(-c1/3) * cos(pi/3 - ceta/3))
        t2 <- suppressWarnings(-2 * sqrt(-c1/3) * cos(pi/3 + ceta/3))
        t3 <- suppressWarnings(2 * sqrt(-c1/3) * cos(ceta/3))
        p01 = t1 - b1/3; p02 = t2 - b1/3; p03 = t3 - b1/3
        p0sum = p01 + p02 + p03; p0up = min(p01, p02, p03); p0low = p0sum - p0up - max(p01, p02, p03)
        
        
        U <- function(a){
         p.hatf <- function(a){
        (a * (m + y) + x + n - ((a * (m + y) + x + n)^2 - 4 * a * (m + n) * (x + y))^0.5)/(2 * (m + n))
        }
        p.hat <- p.hatf(a)
        (((x - m * p.hat)^2)/(m * p.hat * (1 - p.hat)))*(1 + (m * (a - p.hat))/(n * (1 - p.hat))) - x2
        }
        
        rat <- (x/m)/(y/n); nrat <- (x/m)/(y/n); varhat <- (1/x) - (1/m) + (1/y) - (1/n) 
        if((x == 0) & (y != 0)) {nrat <- ((x + 0.5)/m)/(y/n); varhat <- (1/(x + 0.5)) - (1/m) + (1/y) - (1/n)} 
        if((y == 0) & (x != 0)) {nrat <- (x/m)/((y + 0.5)/n); varhat <- (1/x) - (1/m) + (1/(y + 0.5)) - (1/n)}
        if((y == n) & (x == m)) {nrat <- 1; varhat <- (1/(m - 0.5)) - (1/m) + 1/(n - 0.5) - (1/n)}
        
        La <- nrat * exp(-1 * z.star * sqrt(varhat)) * 1/4
        Lb <- nrat
        Ha <- nrat
        Hb <- nrat * exp(z.star * sqrt(varhat)) * 4
        
     #----------------------------------------------------------------------------#
     
        if((x != 0) & (y == 0)) {
        if(x == m){
        CIL = (1 - (m - x) * (1 - p0low)/(y + m - (n + m) * p0low))/p0low
        CIU <- Inf
        }
        else{         
        CIL <- uniroot(U, c(La, Lb), tol=tol)$root
        CIU <- Inf
        }
        }

        #------------------------------------------------------------#  

        if((x == 0) & (y != n)) {
        CIU <- uniroot(U, c(Ha, Hb), tol=tol)$root
        CIL <- 0
        }
        
        #------------------------------------------------------------#  

        if(((x == m)|(y == n)) & (y != 0)){

          
            if((x == m)&(y == n)){
                U.0 <- function(a){if(a <= 1) {m * (1 - a)/a - x2}
                else{(n * (a - 1)) - x2}
                }                                         
                CIL <- uniroot(U.0, c(La, rat), tol = tol)$root
                CIU <- uniroot(U.0, c(rat, Hb), tol = tol)$root
            }
            
            #------------------------------------------------------------#
            
            if((x == m) & (y != n)){ 
                           
                phat1 = x/m; phat2 = y/n
                phihat = phat2/phat1
                phiu = 1.1 * phihat
                r = 0
                while (r >= -z.star) {
                    a = (m + n) * phiu
                    b = -((x + n) * phiu + y + m)
                    c = x + y
                    p1hat = (-b - sqrt(b^2 - 4 * a * c))/(2 * a)
                    p2hat = p1hat * phiu
                    q2hat = 1 - p2hat
                    var = (m * n * p2hat)/(n * (phiu - p2hat) +
                      m * q2hat)
                    r = ((y - n * p2hat)/q2hat)/sqrt(var)
                    phiu1 = phiu
                    phiu = 1.0001 * phiu1
                }
                CIU = (1 - (m - x) * (1 - p0up)/(y + m - (n + m) * p0up))/p0up
                CIL = 1/phiu1
            }
           
            #------------------------------------------------------------#
              
            if((y == n) & (x != m)){
                p.hat2 = y/n; p.hat1 = x/m; phihat = p.hat1/p.hat2
                phil = 0.95 * phihat; r = 0
              if(x != 0){  
                while(r <= z.star) {
                    a = (n + m) * phil
                    b = -((y + m) * phil + x + n)
                    c = y + x
                    p1hat = (-b - sqrt(b^2 - 4 * a * c))/(2 * a)
                    p2hat = p1hat * phil
                    q2hat = 1 - p2hat
                    var = (n * m * p2hat)/(m * (phil - p2hat) +
                      n * q2hat)
                    r = ((x - m * p2hat)/q2hat)/sqrt(var)
                    CIL = phil
                    phil = CIL/1.0001
                    }
                 }

                phiu = 1.1 * phihat
                
                if(x == 0){CIL = 0; phiu <- ifelse(n < 100, 0.01, 0.001)}
                
                r = 0
                while(r >= -z.star) {
                    a = (n + m) * phiu
                    b = -((y + m) * phiu + x  + n)
                    c = y + x
                    p1hat = (-b - sqrt(b^2 - 4 * a * c))/(2 * a)
                    p2hat = p1hat * phiu
                    q2hat = 1 - p2hat
                    var = (n * m * p2hat)/(m * (phiu - p2hat) +
                      n * q2hat)
                    r = ((x  - m * p2hat)/q2hat)/sqrt(var)
                    phiu1 = phiu
                    phiu = 1.0001 * phiu1
                    }
                CIU <- phiu1
                }
          } else if((y != n) & (x != m) & (x != 0) & (y != 0)){       
            CIL <- uniroot(U, c(La, Lb), tol=tol)$root
            CIU <- uniroot(U, c(Ha, Hb), tol=tol)$root
          }
    }
CI <- c(rat, CIL, CIU)
}

#-------------------------- Noether ------------------------------#

if(method == "noether"){
    if((x == 0 & y == 0)|(x == 0 & y != 0)|(x != 0 & y == 0)|(x == m & y == n)){
    if(x == 0 & y == 0) {CIL <- 0;  CIU <- Inf; rat = 0; se.hat <- NA; varhat = NA}
    if(x == 0 & y != 0) {rat <- (x/m)/(y/n); CIL <- 0;  x <- 0.5
    nrat <- (x/m)/(y/n); se.hat <- nrat * sqrt((1/x) - (1/m) + (1/y) - (1/n))
    CIU <- nrat + z.star * se.hat}
    if(x != 0 & y == 0) {rat <- Inf; CIU <- Inf;  y <- 0.5
    nrat <- (x/m)/(y/n); se.hat <- nrat * sqrt((1/x) - (1/m) + (1/y) - (1/n))
    CIL <- nrat - z.star * se.hat}
    if(x == m & y == n) {
       rat <- (x/m)/(y/n); x <- m - 0.5; y <- n - 0.5; nrat <- (x/m)/(y/n); se.hat <- nrat * sqrt((1/x) - (1/m) + (1/y) - (1/n))
       CIU <- nrat + z.star * se.hat 
       CIL <- nrat - z.star * se.hat
       }
    } else
    {
    rat <- (x/m)/(y/n)
    se.hat <- rat * sqrt((1/x) - (1/m) + (1/y) - (1/n))
    CIL <- rat - z.star * se.hat 
    CIU <- rat + z.star * se.hat
    }
varhat <- ifelse(is.na(se.hat), NA, se.hat^2)
CI <- c(rat, max(0,CIL), CIU)
}

#------------------------- sinh-1 -----------------------------#

if(method == "sinh-1"){

if((x == 0 & y == 0)|(x == 0 & y != 0)|(x != 0 & y == 0)|(x == m & y == n)){
    if(x == 0 & y == 0) {CIL <- 0;  CIU <- Inf; rat = 0; varhat = NA}
    if(x == 0 & y != 0) {rat <- (x/m)/(y/n); CIL <- 0;  x <- z.star^2
    nrat <- (x/m)/(y/n); varhat <- 2 * asinh((z.star/2)*sqrt(1/x + 1/y - 1/m - 1/n))
    CIU <- exp(log(nrat) + varhat)}
    if(x != 0 & y == 0) {rat = Inf; CIU <- Inf;  y <- z.star^2
    nrat <- (x/m)/(y/n); varhat <- 2 * asinh((z.star/2)*sqrt(1/x + 1/y - 1/m - 1/n))
    CIL <- exp(log(nrat) - varhat)}
    if(x == m & y == n) {
       rat <- (x/m)/(y/n); x <- m - 0.5; y <- n - 0.5; nrat <- (x/m)/(y/n); varhat <- 2 * asinh((z.star/2)*sqrt(1/x + 1/y - 1/m - 1/n))
       CIL <- exp(log(nrat) - varhat) 
       CIU <- exp(log(nrat) + varhat) 
       }
    } else
    {rat <- (x/m)/(y/n); varhat <- 2 * asinh((z.star/2)*sqrt(1/x + 1/y - 1/m - 1/n))
    CIL <- exp(log(rat) - varhat) 
    CIU <- exp(log(rat) + varhat)
    } 
CI <- c(rat, CIL, CIU)
}

#------------------------Results ------------------------------#
                  
res <- list(CI = CI, varhat = varhat)
res
}

CI <- matrix(ncol = 3, nrow = length(y1))
vh <- rep(NA, length(y1))

for(i in 1 : length(y1)){
temp <- ci.prat1(x = x[i], m = m[i], y = y[i], n = n[i], conf = conf, method = method, bonf = bonf)
CI[i,] <- temp$CI
vh[i] <- temp$varhat
}

CI <- data.frame(CI)
if(length(y1) == 1) row.names(CI) <- ""
head <- paste(paste(as.character(oconf * 100),"%",sep=""), c("Confidence interval for ratio of binomial proportions"))
if(method == "adj.log")head <- paste(head,"(method=adj-log)")
if(method == "bailey")head <- paste(head,"(method=Bailey)") 
if(method == "boot")head <- paste(head,"(method=percentile bootstrap)") 
if(method == "katz.log")head <- paste(head,"(method=Katz-log)")
if(method == "koopman")head <- paste(head,"(method=Koopman)")
if(method == "noether")head <- paste(head,"(method=Noether)")
if(method == "sinh")head <- paste(head,"(method=sinh^-1)")

if(bonf == TRUE)head <- paste(head, "\n Bonferroni simultaneous intervals, r = ", bquote(.(r)), 
"\n Marginal confidence = ", bquote(.(conf)), "\n", sep = "") 

ends <- c("Estimate",paste(as.character(c((1-oconf)/2,1-((1-oconf)/2))*100),"%", sep=""))
res <- list(varhat = vh, ci = CI, ends = ends, head = head)
class(res) <- "ci"
res
}

