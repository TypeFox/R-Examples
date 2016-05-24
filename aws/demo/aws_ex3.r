require(aws)
if(exists("X11")) X11(,10,10)
fw0 <- function(){
         xy <- rbind(rep(0:255,256),rep(0:255,rep(256,256)))
         indw <- c(1:12,29:48,73:100,133:168,209:256)
         w0 <- matrix(rep(FALSE,256*256),ncol=256)
         w0[indw,] <- TRUE
         w0[,indw] <- !w0[,indw]
         w0 <- w0-.5
	 z <- (xy[1,]-129)^2+(xy[2,]-129)^2
         w0[z<=10000&z>=4900] <- 0
         w0[abs(xy[1,]-xy[2,])<=20&z<4900] <- 0
	 z <- (xy[1,]-225)^2+2*(xy[2,]-30)^2-(xy[1,]-225)*(xy[2,]-30)
         w0[z<=625] <- 0
         w0[z<=625&xy[2,]>27&xy[2,]<31] <- -.5
         w0[z<=625&xy[1,]>223&xy[1,]<227] <- .5
	 z <- ((xy[2,]-225)^2+2*(xy[1,]-30)^2)+(xy[2,]-225)*(xy[1,]-30)
         w0[z<=625] <- 0
         w0[z<=625&xy[1,]>27&xy[1,]<31] <- -.5
         w0[z<=625&xy[2,]>223&xy[2,]<227] <- .5
         w0[((xy[2,]-225)^2+(xy[1,]-225)^2)+1*(xy[2,]-225)*(xy[1,]-225)<=400] <- 0
         w0[((xy[2,]-30)^2+(xy[1,]-30)^2)<=256] <- 0
	 w0
	 }
w0 <- fw0()
image(w0,col=gray((0:255)/255))
title("Original image")
sigma <- readline("Standard deviation of noise:\n Press 'Enter' for sigma=0.25, otherwise provide value of sigma:")
if(is.na(as.numeric(sigma))) sigma <- 0.25 else sigma <- as.numeric(sigma)
if(sigma <= 0) sigma <- 0.1
y <- w0+rnorm(w0,0,sigma)
image(y,col=gray((0:255)/255))
title("Noisy image")
hmax <- readline("Maximal bandwidth:\n Press 'Enter' for hmax=10, otherwise provide value of hmax:")
if(is.na(as.numeric(hmax))) hmax <- 10 else hmax <- as.numeric(hmax)
if(hmax <= 1) hmax <- 10
risk <- readline("Report risks (N/Y):")
if(risk %in% c("y","Y")) u <-w0 else u <- NULL
cat("Run aws \n")
yhat <- aws(y,hmax=hmax,graph=TRUE,u=u)
readline("Press ENTER to show results")
oldpar <- par(mfrow=c(2,2),mar=c(1,1,2,.25),mgp=c(2,1,0))
image(w0,col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Artificial image"))
image(y,col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Noisy image (sigma=",signif(sigma,3),")"))
image(awsdata(yhat,"est"),col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Reconstruction  hmax=",signif(yhat@hmax,3),")"))
image(awsdata(yhat,"sd"),col=gray((0:255)/255),xaxt="n",yaxt="n")
title(paste("Standard deviation of estimates (min:",signif(min(awsdata(yhat,"sd")),3)," max:",signif(max(awsdata(yhat,"sd")),3),")"))
par(oldpar)
if(! readline("keep files and device (N/Y) :") %in% c("y","Y")){ 
rm(fw0,w0,y,sigma,hmax,yhat,u,risk)
dev.off()
}
