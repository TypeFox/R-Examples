resp.check <- function(y, margin = "N", 
                           main = "Histogram and Density of Response",
                           xlab = "Response", print.par = FALSE, plots = TRUE, loglik = FALSE, ...){

m2 <- c("N","GU","rGU","LO","LN","WEI","iG","GA","GAi","BE","FISK")
m3 <- c("DAGUM","SM")
nu <- NULL

if(!(margin %in% c(m2,m3)) ) stop("Error in margin value. It can be: N, GU, rGU, LO, LN, WEI, iG, GA, GAi, DAGUM, SM, BE, FISK.") 

if(margin %in% c("LN","WEI","iG","GA","GAi","DAGUM","SM","FISK") && min(y, na.rm = TRUE) <= 0) stop("The response must be positive.")
if(margin %in% c("BE") && (min(y, na.rm = TRUE) <= 0 || max(y, na.rm = TRUE) >= 1) ) stop("The response must be in the interval (0,1).")

y <- na.omit(y)

if(margin == "LN") y <- log(y)

margins <- c("probit", margin) # not important to chance probit here

VC <- list(X1 = matrix(1, nrow = length(y), ncol = 1), X1.d2 = 1,
           X2 = matrix(1, nrow = length(y), ncol = 1), X2.d2 = 1,
           X3 = matrix(1, nrow = length(y), ncol = 1), X3.d2 = 1,
           X4 = matrix(1, nrow = length(y), ncol = 1), X4.d2 = 1,
           X5 = matrix(1, nrow = length(y), ncol = 1), X5.d2 = 1,
           X6 = matrix(1, nrow = length(y), ncol = 1), X6.d2 = 1,
           X7 = matrix(1, nrow = length(y), ncol = 1), X7.d2 = 1,
           l.sp1 = 0, l.sp2 = 0, l.sp3 = 0, l.sp4 = 0, l.sp5 = 0, l.sp6 = 0, l.sp7 = 0, 
           weights = 1, 
           margins = margins, fp = TRUE,
           extra.regI = "t", Cont = "NO", ccss = "no", triv = FALSE)

respvec <- list(y2 = y, univ = 0)
           
if( margin %in% c("N","LN") )     st.v <- c( mean((y + mean(y))/2) ,           log( var(y) ) )  
if( margin %in% c("LO") )         st.v <- c( mean((y + mean(y))/2) ,           log(  3*var(y)/pi^2 ) )  
if( margin %in% c("iG") )         st.v <- c( log( mean((y + mean(y))/2) ) , log( var(y)/mean(y)^3)  )    
if( margin %in% c("GU") )         st.v <- c( mean(y) + 0.57722*sqrt(var(y)/1.64493) ,  log(6*var(y)/pi^2) )    
if( margin %in% c("rGU") )        st.v <- c( mean(y) - 0.57722*sqrt(var(y)/1.64493) ,  log(6*var(y)/pi^2) )   
if( margin %in% c("WEI") )        st.v <- c( log( mean( exp(log(y) + 0.5772/(1.283/sqrt(var(log(y))))) )  ) , log( ( 1.283/sqrt(var(log(y))) )^2 ) ) 
if( margin %in% c("GA") )         st.v <- c( log(mean((y + mean(y))/2)), log(var(y)/mean(y)^2)  ) # log( 1^2 )             
if( margin %in% c("GAi") )        st.v <- c( mean((y + mean(y))/2), log(var(y)/mean(y)^2)  ) # log( 1^2 )             
if( margin %in% c("DAGUM","SM") ) st.v <- c( log(mean((y + mean(y))/2)), log(sqrt(2)), log(1) )   # log(0.01) #  log(sqrt(2))       # 0.1    
if( margin %in% c("FISK") )       st.v <- c( log(mean((y + mean(y))/2)), log(sqrt(2)))    
if( margin %in% c("BE") )         st.v <- c( qlogis(mean((y + mean(y))/2)), qlogis( var(y)/( mean(y)*(1-mean(y)) )  )  )              



if( margin %in% m2 ) names(st.v) <- c("mu.star", "sigma2.star")
if( margin %in% m3 ) names(st.v) <- c("mu.star", "sigma2.star", "nu.star")

#if( margin %in% m2 ) names(st.v) <- c("mu", "sigma")
#if( margin %in% m3 ) names(st.v) <- c("mu", "sigma", "nu")


if(margin %in% m2) univfit <-  try(trust(bprobgHsContUniv, st.v, rinit = 1, rmax = 100, respvec = respvec, 
                                            VC = VC, sp = NULL, qu.mag = NULL, blather = TRUE), silent = TRUE)
                                            
if(margin %in% m3) univfit <-  try(trust(bprobgHsContUniv3, st.v, rinit = 1, rmax = 100, respvec = respvec, 
                                            VC = VC, sp = NULL, qu.mag = NULL, blather = TRUE), silent = TRUE)                                            
                 
if(class(univfit) == "try-error") stop("The parameters of the chosen distribution could not be estimated. Try a different distribution.")                  
    
    
if(plots == TRUE){ ##       
             

if(margin == "LN") y <- exp(y)

if(margin %in% m2)  pp <- distrHsAT(y, univfit$argument[1], esp.tr(univfit$argument[2], margin)$vrb, 1, margin)
if(margin %in% m3)  pp <- distrHsAT(y, univfit$argument[1], esp.tr(univfit$argument[2], margin)$vrb, exp(univfit$argument[3]), margin)

p <- pp$p2
d <- pp$pdf2

par(mfrow = c(1, 2))
hist(y, freq = FALSE, ylim=c(0, max(d, hist(y, plot = FALSE)$density) ),
     main=main,
     xlab=xlab, ...)

lines(sort(y),d[order(y)],lwd=2)


if(any(is.na(p)) == TRUE) stop("It is not possible to produce a QQ-plot.\nThe chosen distribution (unconditional on covariates)\nis not probably a good fit.")

qqnorm(qnorm(p))
abline(0, 1, col = "red")


if(print.par == TRUE) print( univfit$argument )  


                 } ##
                 
                 
         
#if(print.par == TRUE){         
#         
#mu    <- eta.tr(univfit$argument[1], margin)
#mupos <- c("LN","WEI","iG","GA","DAGUM","SM")
#mub   <- c("BE")
#if(margin %in% mupos) mu <- exp(mu)
#if(margin %in% mub)   mu <- plogis(mu)
#
#sigma <- esp.tr(univfit$argument[2], margin)$vrb
#
#if(margin %in% m3) nu <- esp.tr(univfit$argument[3], margin)$vrb                 
#                 
#}


                 
#if(print.par == TRUE && plots == TRUE) print( c(mu,sigma,nu) )                

if(plots == FALSE && print.par == TRUE && loglik == FALSE) return( univfit$argument )  



if(loglik == TRUE){ ##

if(margin == "LN") lk <- sum(log(d))
if(margin != "LN") lk <- -univfit$l

attr(lk, "nobs") <- length(y)
attr(lk, "df") <- length(st.v)
class(lk) <- "logLik"
lk

                  } ##


}

