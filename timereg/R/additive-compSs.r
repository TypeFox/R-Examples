
additive.compSs<-function (formula = formula(data), data = sys.parent(), 
start.time=0,max.time=NULL,id=NULL,scale=FALSE,silent=0,omit=NULL) 
{ ## {{{
   call <- match.call()
    m <- match.call(expand.dots = FALSE)
    m$start.time <- m$max.time <-  m$id <- m$scale<- m$silent<- m$omit<-NULL
    special <- c("const")
    Terms <- if (missing(data)) terms(formula, special)
    else terms(formula, special, data = data)
    m$formula <- Terms
    m[[1]] <- as.name("model.frame")
    m <- eval(m, sys.parent())
    mt <- attr(m, "terms")
    intercept <- attr(mt, "intercept")
    Y <- model.extract(m, "response")
    if (!inherits(Y, "Surv")) stop("Response must be a survival object")

des<-read.design(m,Terms)
X<-des$X; Z<-des$Z; npar<-des$npar; px<-des$px; pz<-des$pz;
covnamesX<-des$covnamesX; covnamesZ<-des$covnamesZ
if (scale==TRUE) Z<-scale(Z); 
pxz <- px + pz; 

clusters=NULL;
survs<-read.surv(m,id,FALSE,clusters,start.time,max.time)
times<-survs$times;id<-id.call<-survs$id.cal;
clusters<-cluster.call<-survs$clusters; 
time2<-survs$stop; time<-survs$start
status<-survs$status; 

if (!is.null(omit)) {
   time<-time[-omit]
   time2<-time2[-omit]; status<-status[-omit]
   X<-X[-omit,]; Z<-Z[-omit,]; id<-id[omit]; 
   times<-time2
   if (is.null(max.time) == TRUE) maxtimes <- max(times)
   else maxtimes <- max.time
   times <- times[times <= maxtimes]
   times <- c(times, maxtimes)
   times <- unique(times)
}
   Nalltimes <- length(times)
   Ntimes <- 
   sum(status[(time2 > times[1]) & (time <= times[Nalltimes])]) + 1
   nx<-nrow(X); 

   fix<-0; 
   if ( (attr(m[, 1], "type") == "right" ) ) {
     fix<-1
     if (fix==1) ot<-order(time2,status==1); # order in time, status=1 first for ties
     if (fix==2) ot<-order(time2,status==0); # order in time, status=1 first for ties
     time2<-time2[ot]; status<-status[ot]; time<-time[ot]; 
     X<-as.matrix(X[ot,])
     Z<-as.matrix(Z[ot,])
     survs$stop<-time2; 
     id<-(1:nx)-1; 
  }

   deltaweight<-1; 
   intZHZ<-matrix(0,pz,pz); intZHdN<-matrix(0,pz,1); 


if (fix==0) {
semiout<-.C("compSs",
as.double(times),as.integer(Nalltimes),as.integer(Ntimes),
as.double(X),as.integer(nx),as.integer(px),
 as.double(Z),as.integer(nx),as.integer(pz),
 as.integer(survs$antpers),as.double(time),as.double(time2),
 as.integer(id),as.integer(status), as.integer(deltaweight),
 as.double(intZHZ),as.double(intZHdN),as.integer(silent)
,package="timereg")
}
if (fix==1) {
semiout<-.C("compSsrev",
as.double(times),as.integer(Nalltimes),as.integer(Ntimes),
as.double(X),as.integer(nx),as.integer(px),
as.double(Z),as.integer(nx),as.integer(pz),
as.integer(survs$antpers),as.double(time),as.double(time2),
as.integer(id),as.integer(status), as.integer(deltaweight),
as.double(intZHZ),as.double(intZHdN),as.integer(silent)
,PACKAGE="timereg")
} 
if (fix==2) {
semiout<-.C("compSsforward",
as.double(times),as.integer(Nalltimes),as.integer(Ntimes),
as.double(X),as.integer(nx),as.integer(px),
as.double(Z),as.integer(nx),as.integer(pz),
as.integer(survs$antpers),as.double(time),as.double(time2),
as.integer(id),as.integer(status), as.integer(deltaweight),
as.double(intZHZ),as.double(intZHdN),as.integer(silent)
,PACKAGE="timereg")
}

intZHZ=matrix(semiout[[16]],pz,pz); 
intZHdN=matrix(semiout[[17]],pz,1); 

ud<-list(intZHZ=intZHZ,intZHdN=intZHdN)
class(ud) <- "pls"
attr(ud, "Call") <- call
attr(ud, "Formula") <- formula
ud$call <- call

return(ud)
} ## }}}


