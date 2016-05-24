plotnPV2 <-
function(x, NPVlty=1, PPVlty=3, ...)
{

makeMAT<-function(n, byrow=FALSE)
{
NCOL=ceiling(sqrt(n))
NROW=ceiling(n/NCOL)
MAX<-NCOL*NROW
IND<-integer(length=MAX)
IND[1:n]<-1:n
if(MAX>n){IND[(n+1):MAX]<-(n+1)}
return(matrix(IND, ncol=NCOL, nrow=NROW, byrow=byrow))
}

PARGS<-list(...)

NPLOTS<-x$NSETS
XLIM<-range(x$propP)

ind<-1:NPLOTS

SETNAMES<-rownames(x$outDAT)

defplot<-function(x, NPVlty, PPVlty, ...)
{

PARGS<-list(...)

 if(is.null(PARGS$xlab))
  {PARGS$xlab<-"Ratio n1/(n0+n1)"}
 if(is.null(PARGS$ylab))
  {PARGS$ylab<-"Sample size"}
 if(is.null(PARGS$type))
  {PARGS$type<-"l"}
 if(is.null(PARGS$lty))
  {PARGS$lty<-NPVlty}
 if(is.null(PARGS$xlim))
  {PARGS$xlim<-XLIM}

print(PARGS)

ACT<-TRUE

if(all(is.na(x$NPV$n)) & all(is.na(x$PPV$n)))
 {ACT<-FALSE
  plot(x=0, y=0, xaxt="n", yaxt="n", bty="n", xlab="", ylab="", main=PARGS$main)
 }

if(ACT && !all(is.na(x$NPV$n)) & all(is.na(x$PPV$n)))
 {
  ACT<-FALSE
 if(is.null(PARGS$ylim))
  {PARGS$ylim<-range(c(x$NPV$n))}
  PARGS$y<-x$NPV$n
  PARGS$x<-x$NPV$propP
  do.call("plot", PARGS)
 }

if(ACT && all(is.na(x$NPV$n)) & !all(is.na(x$PPV$n)))
 {
  ACT<-FALSE
 if(is.null(PARGS$ylim))
  {PARGS$ylim<-range(c(x$PPV$n))}
  PARGS$y<-x$PPV$n
  PARGS$x<-x$PPV$propP
  do.call("plot", PARGS)
 }

if(ACT && !all(is.na(x$NPV$n)) & !all(is.na(x$PPV$n)))
 {
 if(is.null(PARGS$ylim))
  {PARGS$ylim<-range(c(x$NPV$n, x$PPV$n))}
  PARGS$y<-x$NPV$n
  PARGS$x<-x$NPV$propP
  do.call("plot", PARGS)
  lines(y=x$PPV$n, x=x$PPV$propP, lty=PPVlty)
 }

}

MAT<-makeMAT(NPLOTS)

old.par<-par(no.readonly = TRUE)
par(mar=c(4.1,4.1,3.1,1), oma=c(0,0,0,0))
layout(MAT)

for(i in ind)
{
PARGS <- list(...)
PARGS$NPVlty <- NPVlty
PARGS$PPVlty <- PPVlty
PARGS$x <- x$nlist[[i]]

if(is.null(PARGS$main))
 {PARGS$main<-SETNAMES[i]}

do.call("defplot", PARGS)
}

par(old.par)

}

