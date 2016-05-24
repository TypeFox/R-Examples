plotnPV <-
function(x, NPVpar=NULL, PPVpar=NULL, legpar=NULL, ...)
{

NSET <- x$NSET
propP <- x$propP

NPVn<-NULL
PPVn<-NULL

for(i in 1:NSET){
NPVn <- rbind(NPVn, rep(x$nlist[[i]][["NPV"]][["n"]], length.out=x$nsteps))
PPVn <- rbind(PPVn, rep(x$nlist[[i]][["PPV"]][["n"]], length.out=x$nsteps))
}

# default graphic parameters
if(is.null(NPVpar)){NPVPAR<-list(lty=1:NSET, lwd=rep(1,NSET), col=rep("black",NSET), pch=rep(1,NSET), type=rep("l", NSET))}
else{
NPVPAR<-list()
if(!is.list(NPVpar)){stop("If specified, NPVpar should be a list (with elements lts, lwd, col, pch)!")}
if(is.null(NPVpar$lty)){NPVPAR$lty<-rep(1, NSET)}else{NPVPAR$lty<-rep(NPVpar$lty, length.out=NSET)}
if(is.null(NPVpar$lwd)){NPVPAR$lwd<-rep(1, NSET)}else{NPVPAR$lwd<-rep(NPVpar$lwd, length.out=NSET)}
if(is.null(NPVpar$col)){NPVPAR$col<-rep("black", NSET)}else{NPVPAR$col<-rep(NPVpar$col, length.out=NSET)}
if(is.null(NPVpar$pch)){NPVPAR$pch<-rep(1, NSET)}else{NPVPAR$pch<-rep(NPVpar$pch, length.out=NSET)}
}

if(is.null(PPVpar)){PPVPAR<-list(lty=1:NSET, lwd=rep(1,NSET), col=rep("grey",NSET), pch=rep(1,NSET), type=rep("l", NSET))}
else{
PPVPAR<-list()
if(!is.list(PPVpar)){stop("If specified, PPVpar should be a list (with elements lts, lwd, col, pch, and/or type)!")}
if(is.null(PPVpar$lty)){PPVPAR$lty<-rep(1, NSET)}else{PPVPAR$lty<-rep(PPVpar$lty, length.out=NSET)}
if(is.null(PPVpar$lwd)){PPVPAR$lwd<-rep(1, NSET)}else{PPVPAR$lwd<-rep(PPVpar$lwd, length.out=NSET)}
if(is.null(PPVpar$col)){PPVPAR$col<-rep("black", NSET)}else{PPVPAR$col<-rep(PPVpar$col, length.out=NSET)}
if(is.null(PPVpar$pch)){PPVPAR$pch<-rep(1, NSET)}else{PPVPAR$pch<-rep(PPVpar$pch, length.out=NSET)}
}

PARGS<-list(...)

RN<-range(na.omit(NPVn), na.omit(PPVn))
RP<-range(na.omit(propP))

if(is.null(PARGS$ylim)){PARGS$ylim=RN}
if(is.null(PARGS$xlim)){PARGS$xlim=RP}
if(is.null(PARGS$ylab)){PARGS$ylab="Total sample size, (n0+n1)"}
if(is.null(PARGS$xlab)){PARGS$xlab="Proportion of true positives, n1/(n0+n1)"}
if(is.null(PARGS$type)){LTYPE<-"l"}else{LTYPE<-PARGS$type}

PARGS$type<-"n"
PARGS$x<-RP
PARGS$y<-RN

do.call(what="plot", args=PARGS)

for(i in 1:NSET)
{
lines(y=NPVn[i,], x=propP, lty=NPVPAR$lty[i], lwd=NPVPAR$lwd[i], col=NPVPAR$col[i], pch=NPVPAR$pch[i], type=LTYPE)
lines(y=PPVn[i,], x=propP, lty=PPVPAR$lty[i], lwd=PPVPAR$lwd[i], col=PPVPAR$col[i], pch=PPVPAR$pch[i], type=LTYPE)
}

# building the legend

if(is.null(legpar)){legpar<-list()
legpar$bty="n"
legpar$x<-PARGS$xlim[1]+0.1*abs(diff(PARGS$xlim)); legpar$y<-PARGS$ylim[2]
lnam<-rownames(x$outDAT)
legpar$legend<-c(paste("NPV", lnam), paste("PPV", lnam))
legpar$lty=c(NPVPAR$lty, PPVPAR$lty)
legpar$lwd=c(NPVPAR$lwd, PPVPAR$lwd)
legpar$col=c(NPVPAR$col, PPVPAR$col)
if(LTYPE %in%c("b", "p")){legpar$pch=c(NPVPAR$pch, PPVPAR$pch)}else{legpar$pch<-NULL}
do.call(what="legend", args=legpar)

}else{
if(!is.list(legpar) && (is.numeric(legpar) &  legpar==0)){
legpar$bty="n"
legpar$x<-PARGS$xlim[1]+0.1*abs(diff(PARGS$xlim)); legpar$y<-PARGS$ylim[2]
lnam<-rownames(x$outDAT)
legpar$legend<-c(paste("NPV", lnam), paste("PPV", lnam))
legpar$lty=c(NPVPAR$lty, PPVPAR$lty)
legpar$lwd=c(NPVPAR$lwd, PPVPAR$lwd)
legpar$col=c(NPVPAR$col, PPVPAR$col)
if(LTYPE %in%c("b", "p")){legpar$pch=c(NPVPAR$pch, PPVPAR$pch)}else{legpar$pch<-NULL}
return(legpar)
}else{
if(!is.list(legpar)){stop("If specified, legpar should be a list!")}
if(is.null(legpar$bty)){legpar$bty="n"}
if(is.null(legpar$x)){legpar$x<-PARGS$xlim[1]+0.1*abs(diff(PARGS$xlim))}
if(is.null(legpar$y)){legpar$y<-PARGS$ylim[2]}
lnam<-rownames(x$outDAT)
legpar$legend<-c(paste("NPV", lnam), paste("PPV", lnam))
legpar$lty=c(NPVPAR$lty, PPVPAR$lty)
legpar$lwd=c(NPVPAR$lwd, PPVPAR$lwd)
legpar$col=c(NPVPAR$col, PPVPAR$col)
if(LTYPE %in%c("b", "p")){legpar$pch=c(NPVPAR$pch, PPVPAR$pch)}else{legpar$pch<-NULL}
do.call(what="legend", args=legpar)
}}

}

