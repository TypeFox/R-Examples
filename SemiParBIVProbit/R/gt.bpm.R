gt.bpm <- function(x) {


if(x$Cont == "YES") stop("This function is not suitable for bivariate models with continuous margins.")


dl.drho <- d.n1n2 <- p01 <- p00 <- p10 <- p11 <- eta1 <- eta2 <- 1
if(x$BivD!="N") stop("This test's implementation is currently valid for bivariate normal errors only.")
if(!is.null(x$X3) ) stop("This test is not designed for a varying correlation coefficient model.")
if( !(x$margins[2] == "probit") ) stop("This test is not designed for bivariate models with continuous response.")

eta1 <- (x$gam1$linear.predictors)
eta2 <- (x$VC$X2 %*% coef(x$gam2) )

inde <- rep(TRUE, length(eta1))

if(x$Model=="B")  p11 <- pbinorm(eta1, eta2, cov12 = 0)
if(x$Model=="BSS"){inde <- x$y1 > 0; p11 <- pbinorm(eta1[inde], eta2, cov12 = 0)}


if(x$Model=="B")   p10 <- pnorm(eta1) - p11
if(x$Model=="BSS") p10 <- pnorm(eta1[inde]) - p11

if(x$Model=="B"){ p01 <- pnorm(eta2) - p11
                  p00 <- 1 - p11 - p10 - p01
                 }

if(x$Model=="B")   d.n1n2 <- dbinorm(eta1, eta2, cov12 = 0)
if(x$Model=="BSS") d.n1n2 <- dbinorm(eta1[inde], eta2, cov12 = 0)

if(x$Model=="B")   dl.drho <- x$weights*d.n1n2*(x$respvec$y1.y2/p11 - x$respvec$y1.cy2/p10 - x$respvec$cy1.y2/p01 + x$respvec$cy1.cy2/p00) 
if(x$Model=="BSS") dl.drho <- x$weights[inde]*d.n1n2*(x$respvec$y1.y2/p11 - x$respvec$y1.cy2/p10)
if(x$Model=="BPO") dl.drho <- x$weights*d.n1n2*(x$y1/p11 - (1-x$y1)/(1-p11))  

G <- as.numeric(sum(dl.drho)*x$theta)

rm(dl.drho, d.n1n2, p01, p00, p10, p11, eta1, eta2)

return(pchisq(G, 1, lower.tail = FALSE))

}
