############################################################################
#                              EXAMPLE 1:                                  # 
#            Simple example for a Gaussian density derivative              #
############################################################################        

x <- rnorm(100)
SD0 <- dkde(x,deriv.order=0)
SD1 <- dkde(x,deriv.order=1)
SD2 <- dkde(x,deriv.order=2)
SD3 <- dkde(x,deriv.order=3)
dev.new()
par(mfrow=c(2,2))
plot(SD0)
plot(SD1)
plot(SD2)
plot(SD3)

############################################################################
#                              EXAMPLE 2                                   # 
#                Trimodal Gaussian density derivative                      # 
#                Computing bandwidths with UCV methods                     #
############################################################################

data(trimodal)
h.ucv(trimodal,deriv.order=0,kernel="gaussian")
h.ucv(trimodal,deriv.order=1,kernel="gaussian")
h.ucv(trimodal,deriv.order=2,kernel="gaussian")
h.ucv(trimodal,deriv.order=3,kernel="gaussian")

############################################################################
#                             Example 3                                    #
#          Computing bandwidths with BCV methods different kernels         #       
#                        derivative order = 0                              #
############################################################################

data(outlier)
h.bcv(outlier,deriv.order=0,kernel="gaussian")
h.bcv(outlier,deriv.order=0,kernel="triweight")
h.bcv(outlier,deriv.order=0,kernel="tricube")
h.bcv(outlier,deriv.order=0,kernel="biweight")
h.bcv(outlier,deriv.order=0,kernel="cosine")

############################################################################
#                        derivative order = 1                              #
############################################################################

h.bcv(outlier,deriv.order=1,kernel="gaussian")
h.bcv(outlier,deriv.order=1,kernel="triweight")
h.bcv(outlier,deriv.order=1,kernel="tricube")
h.bcv(outlier,deriv.order=1,kernel="biweight")
h.bcv(outlier,deriv.order=1,kernel="cosine")

############################################################################
#                        derivative order = 2                              #
############################################################################

h.bcv(outlier,deriv.order=2,kernel="gaussian")
h.bcv(outlier,deriv.order=2,kernel="triweight")
h.bcv(outlier,deriv.order=2,kernel="tricube")
h.bcv(outlier,deriv.order=2,kernel="biweight")
h.bcv(outlier,deriv.order=2,kernel="cosine")

############################################################################
#                             Example 4                                    #
#                   Bimodal Gaussian density derivative                    #
############################################################################
fx  <- function(x) 0.5 * dnorm(x,-1.5,0.5) + 0.5 * dnorm(x,1.5,0.5)
fx1 <- function(x) 0.5 *(-4*x-6)* dnorm(x,-1.5,0.5) + 0.5 *(-4*x+6) * 
                   dnorm(x,1.5,0.5)
############################################################################
#                        derivative order = 0                              #
############################################################################

kernels <- eval(formals(dkde.default)$kernel)
dev.new()
plot(dkde(bimodal,h=0.3),sub=paste("Derivative order = 0",";",
     "Bandwidth =0.3 "),ylim=c(0,0.5), main = "Bimodal Gaussian Density")
for(i in 2:length(kernels))
   lines(dkde(bimodal, h = 0.3, kernel =  kernels[i]), col = i)
curve(fx,add=TRUE,lty=8)
legend("topright", legend = c(TRUE,kernels), col = c("black",seq(kernels)),
          lty = c(8,rep(1,length(kernels))),cex=0.7, inset = .015)
	   
############################################################################
#                        derivative order = 1                              #
############################################################################

kernels <- eval(formals(dkde.default)$kernel)[-3]
dev.new()
plot(dkde(bimodal,deriv.order=1,h=0.6),main = "Bimodal Gaussian Density Derivative",sub=paste
         ("Derivative order = 1",";","Bandwidth =0.6"),ylim=c(-0.6,0.6))
for(i in 2:length(kernels))
   lines(dkde(bimodal,deriv.order=1, h = 0.6, kernel =  kernels[i]), col = i)
curve(fx1,add=TRUE,lty=8)
legend("topright", legend = c(TRUE,kernels), col = c("black",seq(kernels)),
          lty = c(8,rep(1,length(kernels))),cex=0.7, inset = .015)

############################################################################
#                             Example 5                                    #
#                    Show the bandwidth selection                          #       
#               kernel = "gaussian" ; derivative order = 0                 #
############################################################################		  

############################################################################
#                 KDE of f (bimodal gaussian density)                      #
############################################################################

hbcv1  <- h.bcv(x=bimodal,whichbcv = 1,deriv.order = 0)$h
hbcv2  <- h.bcv(x=bimodal,whichbcv = 2,deriv.order = 0)$h
hucv   <- h.ucv(x=bimodal,deriv.order = 0)$h
htcv   <- h.tcv(x=bimodal,deriv.order = 0)$h
hccv   <- h.ccv(x=bimodal,deriv.order = 0)$h
hmcv   <- h.mcv(x=bimodal,deriv.order = 0)$h
h0 <- c(hbcv1,hbcv2,hucv,htcv,hccv,hmcv)
h0
dev.new()
plot(dkde(x=bimodal,deriv.order = 0,h=h0[1]),ylim=c(0,0.5),
      sub=paste("Kernel Gaussian",";","Derivative order = 0"),
	  main="Bimodal Gaussian density")
for(i in 1:length(h0)) lines(dkde(x=bimodal,deriv.order = 0,h=h0[i]), col = i)
curve(fx,lty=8,add=TRUE)
legend("topright",title="Bandwidth", c("True",expression(h[bcv1]),
       expression(h[bcv2]),expression(h[ucv]),expression(h[tcv]),
       expression(h[ccv]),expression(h[mcv])),
       lty=c(8,rep(1,length(h0))),col= c("black",seq(h0)),inset = .015)

############################################################################
#             KDDE of d/dx f (bimodal gaussian density)                    #
############################################################################

hbcv1  <- h.bcv(x=bimodal,whichbcv = 1,deriv.order = 1,upper=0.5)$h
hbcv2  <- h.bcv(x=bimodal,whichbcv = 2,deriv.order = 1,upper=0.5)$h
hucv   <- h.ucv(x=bimodal,deriv.order = 1)$h
htcv   <- h.tcv(x=bimodal,deriv.order = 1)$h
hccv   <- h.ccv(x=bimodal,deriv.order = 1)$h
hmcv   <- h.mcv(x=bimodal,deriv.order = 1,upper=0.5)$h
h1 <- c(hbcv1,hbcv2,hucv,htcv,hccv,hmcv)
h1
dev.new()
plot(dkde(x=bimodal,deriv.order = 1,h=h1[1]),ylim=c(-0.7,0.7),
          sub=paste("Kernel Gaussian",";","Derivative order = 1"),
          main="Bimodal Gaussian density derivative")
for(i in 1:length(h1)) lines(dkde(x=bimodal,deriv.order = 1,h=h1[i]), col = i)
curve(fx1,lty=8,add=TRUE)
legend("topright",title="Bandwidth", c("True",expression(h[bcv1]),
       expression(h[bcv2]),expression(h[ucv]),expression(h[tcv]),
       expression(h[ccv]),expression(h[mcv])),
       lty=c(8,rep(1,length(h1))),col= c("black",seq(h1)),inset = .015)		  
		  

############################################################################
#                             Example 6                                    #
#                   Bimodal Gaussian density derivative                    #
#                         CCV and MCV plot                                 #
#                        derivative order = 0                              #
############################################################################

data(bimodal)
dev.new()
plot(h.ccv(bimodal),main="CCV vs MCV",ylab="")
lines(h.mcv(bimodal),col="red")
legend("topright", c("CCV","MCV"),lty=c(1,1),col=c("black","red"), inset = .015)

############################################################################
#                        derivative order = 1                              #
############################################################################

dev.new()
plot(h.ccv(bimodal,deriv.order=1),main="CCV vs UCV",ylab="",ylim=c(-0.7,0.3),
     seq.bws=seq(0.05,1,length=50))
lines(h.ucv(bimodal,deriv.order=1),col="red")
legend("topright", c("CCV","UCV"),lty=c(1,1),col=c("black","red"), inset = .015)

############################################################################
#                        derivative order = 2                              #
############################################################################

dev.new()
plot(h.ccv(bimodal,deriv.order=2,upper=0.5),seq.bws=seq(0.1,0.6,length=50),
     main="CCV vs MCV",ylab="")
lines(h.ucv(bimodal,deriv.order=2),col="red")
legend("topright", c("CCV","UCV"),lty=c(1,1),col=c("black","red"), inset = .015)
