"plot.nDesign" <-
function(x,...)

{

args<-list(...)

op<-par(no.readonly = TRUE)

 if(x$alternative=="less")
  {alt.hyp <- paste("true proportion is less than ",x$p.hyp )
   ptrue <- paste(" assumed true proportion = ", x$p.hyp - x$delta)}

 if(x$alternative=="greater")
  {alt.hyp <- paste("true proportion is greater than", x$p.hyp )
   ptrue <- paste(" assumed true proportion = ", x$p.hyp + x$delta)}

 if(x$alternative=="two.sided")
  {alt.hyp <- paste("true proportion is not equal to ",x$p.hyp )
   ptrue <- paste(" assumed true proportion = ", x$p.hyp - x$delta," or ", x$p.hyp + x$delta)}

 stop.plot<-x$maxit
 powervec=x$powerit[1:stop.plot]
 biasvec=x$biasit[1:stop.plot]
 nvec=x$nit[1:stop.plot]

# set the parameters of plotting region 

# Plot of power iteration

pargs<-args

pargs$y <- powervec
pargs$x <- nvec

if(is.null(pargs$type))
 {pargs$type<-"l"}

if(is.null(pargs$lty))
 {pargs$lty<-1}

if(is.null(pargs$lwd))
 {pargs$lwd<-2}

if(is.null(pargs$ylim))
 {pargs$ylim <- c(0,1)}

if(is.null(pargs$xlab))
 {pargs$xlab <- "number of groups n"}

if(is.null(pargs$ylab))
 {pargs$ylab <- "power"}

if(is.null(pargs$main))
 {pargs$main <- paste(" Alternative: ", alt.hyp,",
   ",ptrue )}
if(is.null(pargs$cex.lab))
 { pargs$cex.lab<-1.2}

# layout(mat=matrix(1:2,ncol=1), heights=c(2,1))
 

par(mar=c(5,5,3,1), oma=c(0,0,0,0), mfrow=c(2,1))

do.call("plot", pargs)


 abline(h=x$powerout, lty=3)
 abline(v=x$nout, lty=2)



bargs<-args

bargs$y <- biasvec
bargs$x <- nvec

if(is.null(bargs$type))
 {bargs$type<-"l"}

if(is.null(bargs$lty))
 {bargs$lty<-1}

if(is.null(bargs$lwd))
 {bargs$lwd<-2}

if(is.null(bargs$ylim))
 {bargs$ylim <- c(0, 2*x$biasrest)}

if(is.null(bargs$xlab))
 {bargs$xlab <- "number of groups n"}

if(is.null(bargs$ylab))
 {bargs$ylab <- "bias of estimator"}

if(is.null(bargs$main))
 {bargs$main <- paste("bias of estimator if",ptrue)}

if(is.null(bargs$cex.lab))
 { bargs$cex.lab<-1.2}

# layout(mat=matrix(1:2,ncol=1), heights=c(2,1))
# par(mar=c(4,5,3,1), oma=c(0,0,0,0))

# par(mar=c(5,5,2,1))

do.call("plot", bargs)

 abline(h=x$bias.rest, lty=3)
 abline(v=x$nout, lty=2)
 abline(h=x$biasout, lty=2)

par(op)

}

