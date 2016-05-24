"print.nhat" <- 
function( x, width=10, ... ){

if( !is.null(names(x$n.hat)) ){
	nms <- names(x$n.hat) 
} else {
	nms <- 1:length(x$n.hat)
}

if( x$nhat.v.meth == 1 ){
	mess <- "(Huggins variance estimator)"
} else if (x$nhat.v.meth == 2){
	mess <- "(Huggins variance estimator + higher terms)"
} else if (x$nhat.v.meth == 3){
	mess <- "(McDonald and Amststrup variance estimator)"
} else if (x$nhat.v.meth == 4){
    mess <- "(Model averaged Huggins variance)"
} else if (x$nhat.v.meth == 5){
    mess <- "(Model averaged Huggins variance + higher terms)"
} else if (x$nhat.v.meth == 6){
    mess <- "(Model averaged McDonald and Amstrup)"
} else {
	mess <- "(Unknown)"
}
	
ind <- !is.na( x$n.hat )

if( x$nhat.v.meth >= 4 ){
    cat("MODEL AVERAGES\n")
}
cat(paste( x$n.hat.conf*100, "% Confidence Intervals\n", sep="" ))
cat(paste( "Standard errors computed by method=", x$nhat.v.meth, " ", mess, "\n\n", sep=""))

wid <- width
z.lower.ci <- formatC(x$n.hat.lower, digits=1, format="f", flag=" ", width=wid)
z.upper.ci <- formatC(x$n.hat.upper, digits=1, format="f", flag=" ", width=wid)
n.hat <- formatC(x$n.hat,digits=1, format="f", flag=" ", width=wid)
se.n.hat <- formatC(x$se.n.hat, digits=1, format="f", flag=" ", width=wid)

#z.lower.ci[ind] <- round(x$n.hat.lower[ind])
#z.upper.ci[ind] <- round(x$n.hat.upper[ind])
#n.hat[ind] <- round(n.hat[ind], 1)
#se.n.hat[ind] <- round(se.n.hat[ind], 1)


if( is.null( se.n.hat ) ) se.n.hat<-rep(NA, length(x$n.hat)) 
if( is.null( z.lower.ci ) ) z.lower.ci<-rep(NA, length(x$n.hat))
if( is.null( z.upper.ci ) ) z.upper.ci<-rep(NA, length(x$n.hat))

dashes <- paste(rep("-",width), collapse="")

#n.hat <- formatC( n.hat, format="f" )

out <- paste( formatC( c(" ", "Occasion", dashes, nms), width=wid, format="s"), 
	formatC( c(" ", "N est", dashes, n.hat), width=wid, format="s" ),  
	formatC( c(" ", "Se(N)", dashes, se.n.hat ), width=wid, format="s" ),
	formatC( c("Norm", "Low", dashes, z.lower.ci), width=wid, format="s"), 
	formatC( c("Norm", "High", dashes, z.upper.ci), width=wid, format="s"), 
	sep="  ")

if( exists( "supsmu" ) ){
    smu.n <- supsmu( seq( along=nms[ind] ), x$n.hat[ind], bass=.5 )
	smu.n <- formatC( smu.n$y, digits=1, format="f", flag=" ", width=wid)
	out <- paste( out, 
		formatC( c("Smooth","N est", dashes, c(NA, smu.n) ), width=wid, format="s"),
		"\n", sep="  ")
} else {
	out <- paste( out, "\n")
}

cat(" ")
cat( out )

invisible()

}

