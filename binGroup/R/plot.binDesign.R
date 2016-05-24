"plot.binDesign" <-
function(x,...)

{

args<-list(...)

 if(x$alternative=="less")
  {alt.hyp <- paste("true proportion is less than ",x$p.hyp )
   ptrue <- paste(" assumed true proportion = ", x$p.hyp - x$delta)}

 if(x$alternative=="greater")
  {alt.hyp <- paste("true proportion is greater than", x$p.hyp )
   ptrue <- paste(" assumed true proportion = ", x$p.hyp + x$delta)}

 if(x$alternative=="two.sided")
  {alt.hyp <- paste("true proportion is not equal to ",x$p.hyp )
   ptrue <- paste(" assumed true proportion = ", x$p.hyp - x$delta," or ", x$p.hyp + x$delta)}

 powerout=round(x$powerout, digits=3)
 nout=x$nout
 nvec=x$nit[1:x$maxit]
 powervec=x$powerit[1:x$maxit]


args$y <- powervec
args$x <- nvec
if(is.null(args$type))
 {args$type<-"l"}

if(is.null(args$lty))
 {args$lty<-1}

if(is.null(args$lwd))
 {args$lwd<-2}

if(is.null(args$ylim))
 {args$ylim <- c(0,1)}

if(is.null(args$xlab))
 {args$xlab <- "sample size n"}

if(is.null(args$ylab))
 {args$ylab <- "power"}

if(is.null(args$main))
 {args$main <- paste(" Alternative: ", alt.hyp,",
   ",ptrue )}
if(is.null(args$cex.lab))
 { args$cex.lab<-1.2}

 par(mar=c(5,5,5,1))

 layout(matrix(c(1)))

do.call("plot", args)

 if(x$power.reached==TRUE)
   {abline(h=powerout, lty=3)
    abline(v=nout, lty=2)}

 if(x$power.reached==FALSE && powerout!=0)
   {abline(h=x$powerout, lty=3)
    abline(v=nout, lty=2)}

}

