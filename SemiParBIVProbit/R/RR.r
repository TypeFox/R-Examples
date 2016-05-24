RR <- function(x, nm.end, E = TRUE, treat = TRUE, type = "bivariate", ind = NULL,
   n.sim = 100, prob.lev = 0.05, length.out = NULL, hd.plot = FALSE, rr.plot = FALSE, 
   main = "Histogram and Kernel Density of Simulated Risk Ratios", 
   xlab = "Simulated Risk Ratios", ...){

if(x$Cont == "YES") stop("This function is not suitable for bivariate models with continuous margins.")
if(x$Cont == "NO" && x$VC$ccss == "yes" ) stop("This function is not suitable for bivariate selection models with continuous margin.")



CIs <- est.AT <- NULL



etap.noi <- X.int <- X.noi <- eti1 <- eti0 <- etno <- indS <- bs <- ind.excl <- p.int1 <- p.int0 <- d.int1 <- d.int0 <- p.etn <- d.etn <- ass.p <- ass.pst <- C.11 <- C.10 <- sig2 <- peti1s <- peti0s <- sigma2.st <- sigma2s <- eti1s <- eti0s <- d0 <- d1 <- p.etns <- etnos <- etds <- ass.ps <- delta.AT <- 1
diffEf <- fy1.y2 <- est.ATso <- y2 <- CIF <- Pr <- Effects <- NULL

m2  <- x$VC$m2 
m3  <- x$VC$m3 
bin.link <- x$VC$bl  
end <- 0
epsilon <- 0.0000001 # 0.9999999 sqrt(.Machine$double.eps)
max.p   <- 0.9999999
est.ATb <- NA
indD <- list()

if(x$v1[1] %in% x$v2[-1]) {end <- 1; eq <- 2} 
if(x$v2[1] %in% x$v1[-1]) {end <- 2; eq <- 1}


if( !( type %in% c("naive","univariate","bivariate") ) ) stop("Error in parameter type value. It should be one of: naive, univariate or bivariate.")
if( !(x$margins[2] %in% bin.link) && eq == 2) stop("It does not make sense to calculate RR when the outcome is continuous.")
if(missing(nm.end)) stop("You must provide the name of the endogenous variable.")


if(x$Model=="BSS" || x$Model=="BPO" || x$Model=="BPO0" || end==0) stop("Calculation of this effect is valid for recursive models only.")
if(is.character(nm.end)==FALSE) stop("nm.end is not a character!")
if( !is.null(ind) && E == FALSE) stop("ind is not designed to be used when some observations are excluded from the RR's calculation.")
if( type == "naive" && E == FALSE) stop("It does not make sense to calculate the naive estimate from the treated only.")

if( !is.null(ind) ){ 

    if(is.logical(ind) == FALSE) stop("ind must be a logical variable.")
    if(length(table(ind))!=2 ) stop("ind must be a logical binary variable.")
    if( length(ind) != x$n ) stop("ind must have the same length as the number of observations used in fitting.")   

}


if( is.null(ind) ) ind <- 1:x$n


if(E == FALSE ) {

 if(!(x$margins[2] %in% bin.link)) ind <- 1:x$n  
 
 if(eq==1) X.int <- as.matrix(x$X1[ind,])
 if(eq==2) X.int <- as.matrix(x$X2[ind,]) 

    if(treat == TRUE)  ind <- as.logical(X.int[, nm.end]) 
    if(treat == FALSE) ind <- as.logical(X.int[, nm.end])!=TRUE
                                              
}




##############################################################################################################

if(type == "naive" && x$margins[2] %in% bin.link){ # it looks ok from comparing the naive and univariate options

if(eq==2){
y1 <- x$y1[ind] 
y2 <- x$y2[ind]
}

if(eq==1){
y1 <- x$y2[ind] 
y2 <- x$y1[ind]
}

tab2 <- table(y1, y2)                                  


pY1cT1 <- prop.table(tab2,1)[4] 
pY1cT0 <- prop.table(tab2,1)[3] 

est.AT <- (pY1cT1 / pY1cT0)

sv <- qnorm(prob.lev/2,lower.tail = FALSE) * sqrt( (1-pY1cT1)/(table(y1)[2]*pY1cT1) + (1-pY1cT0)/(table(y1)[1]*pY1cT0) )

CIs <- exp( c(log(est.AT) - sv, log(est.AT) + sv) )


est.ATb <- est.ATso <- NULL

}



##############################################################################################################

if(type == "naive" && !(x$margins[2] %in% bin.link)) stop("Please fit a bivariate model with intercept and endogenous variable only and then use RR with the univariate type option.")

##############################################################################################################





if(type != "naive" && x$margins[2] %in% bin.link) {


########################################################
# Set-up
########################################################


if(type == "bivariate"){

	indD[[1]] <- 1:x$X1.d2 
	indD[[2]] <- x$X1.d2+(1:x$X2.d2)
	
}


if(eq==1){ 
           X.int <- as.matrix(x$X1[ind,])
           if(type == "bivariate") ind.int <- indD[[1]]  
}

if(eq==2){ 
           X.int <- as.matrix(x$X2[ind,])
           if(type == "bivariate") ind.int <- indD[[2]]
}


if(type == "bivariate") coef.int <- as.numeric(coef(x)[ind.int])
	   


               
d0 <- d1 <- X.int
d0[,nm.end] <- 0
d1[,nm.end] <- 1



if(type == "bivariate"){

	eti1 <- d1%*%coef.int 
	eti0 <- d0%*%coef.int 
  
}

if(type == "univariate"){

if(eq==1) ngam <- x$gam1
if(eq==2) ngam <- x$gam2

        eti1 <- d1%*%coef(ngam) 
        eti0 <- d0%*%coef(ngam) 

}




#############################################################################
# RR
#############################################################################

p.int1 <- probm(eti1, x$margins[eq])$pr 
p.int0 <- probm(eti0, x$margins[eq])$pr  

est.AT <- mean(p.int1, na.rm = TRUE) / mean(p.int0, na.rm = TRUE) 
                        
#############################################################################
# CIs RR
#############################################################################


 if(type == "univariate") {bs <- rMVN(n.sim, mean = coef(ngam), sigma=ngam$Vp)
                           eti1s <- d1%*%t(bs)
                           eti0s <- d0%*%t(bs) 
                           }
 if(type == "bivariate")  {bs <- rMVN(n.sim, mean = coef(x), sigma=x$Vb)
                           eti1s <- d1%*%t(bs[,ind.int])
                           eti0s <- d0%*%t(bs[,ind.int]) 
                           } 


 peti1s <- probm(eti1s, x$margins[eq])$pr   
 peti0s <- probm(eti0s, x$margins[eq])$pr 

 est.ATb <- colMeans(peti1s, na.rm = TRUE) / colMeans(peti0s, na.rm = TRUE) 


CIs <- as.numeric(quantile(est.ATb, c(prob.lev/2,1-prob.lev/2), na.rm = TRUE))


  if(hd.plot == TRUE){
  
  hist(est.ATb, freq = FALSE, main=main, 
       xlab=xlab, 
       ylim=c(0,max(density(est.ATb)$y,hist(est.ATb, plot = FALSE)$density)), ...)
  lines(density(est.ATb))

                     }
                     
                     
 
}


















###############################################


if(type != "naive" && !(x$margins[2] %in% bin.link)) {


 
if(is.null(length.out)) length.out <- length( seq( min(ceiling(x$y2)) , max(floor(x$y2)) ) ) 
y2   <- round( seq( min(ceiling(x$y2)) , max(floor(x$y2)), length.out = length.out  ), 2 ) 
 
 ly2  <- length(y2)
 data <- x$dataset[ind,]
 
 if(type == "bivariate")  {
                           ind.int <- 1:x$X1.d2
                           bs <- rMVN(n.sim, mean = coef(x), sigma = x$Vb) 
                           coefe  <- x$coef[ind.int] 
                           coefes <- t(bs[, ind.int]) 
 
                          }
 
 if(type == "univariate") {bs <- rMVN(n.sim, mean = coef(x$gam1), sigma = x$gam1$Vp) 
                           coefe  <- x$gam1$coefficient 
                           coefes <- t(bs) 
                          }
 
 
 
 
 
 
sratio <- function(x1, x2) x1 / x2  
fy1.y2 <- fy1.y2S <- list()
diffE  <- NA 

diffES <- list()
diffEfSquant <- as.data.frame(matrix(NA, ly2 - 1, 2))


for(i in 1:ly2) {

data[, 2]   <- y2[i]
lpm    <- as.matrix( predict.gam(x$gam1, newdata = data, type = "lpmatrix") )
eta1   <- lpm%*%coefe
etins  <- lpm%*%coefes

# lpm%*%bs[5, ind.int]

fy1.y2[[i]]  <- mean( probm(eta1, x$margins[eq])$pr )
fy1.y2S[[i]] <- colMeans( probm(etins, x$margins[eq])$pr  )

}




for(i in 1:(ly2-1)) {

  diffE[i]          <- sratio(fy1.y2[[i+1]] , fy1.y2[[i]])
  diffES[[i]]       <- sratio(fy1.y2S[[i+1]], fy1.y2S[[i]])      
  diffEfSquant[i, ] <- quantile(diffES[[i]], probs = c(prob.lev/2,1-prob.lev/2), na.rm = TRUE) 
                            } 




Effects <- data.frame(Ratios = diffE, diffEfSquant)  
names(Effects)[2:3] <- names(quantile(c(1,1), probs = c(prob.lev/2,1-prob.lev/2)))
dimnames(Effects)[[1]] <- y2[2:ly2]


if(rr.plot == TRUE){

plot(y2[2:ly2], diffE, log = "y", ylab = "RR", xlab = "Treatment", pch = 16, ylim = c(min(diffEfSquant[,1]),max(diffEfSquant[,2])))
lines(y2[2:ly2], diffE, type = "l")
for (i in 1:(ly2-1)) lines( y = c(diffEfSquant[i,1], diffEfSquant[i,2]), x = c(y2[i+1],y2[i+1]))

}



}



####################



res <- c(CIs[1], est.AT, CIs[2])



rm(etap.noi, X.int, X.noi, eti1, eti0, etno, indS, bs, ind.excl, p.int1, p.int0, d.int1, d.int0,
   p.etn, d.etn, ass.p, ass.pst, C.11, C.10, sig2, peti1s, peti0s, sigma2.st, sigma2s, eti1s, eti0s, d0, d1,
   p.etns, etnos, etds, ass.ps)  
   
   
out <- list(res=res, prob.lev=prob.lev, sim.RR=est.ATb, mar2=x$margins[2], type = type,
            Ratios = Effects, treat = y2, eq = eq, bl = x$VC$bl) # Pr = Pr,
 
  
 
class(out) <- "RR"

out





}





















