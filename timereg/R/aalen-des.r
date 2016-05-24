des.aalen<-function (formula = formula(data),
data = sys.parent(), start.time = 0, max.time = NULL, 
id=NULL, clusters=NULL, deltaweight=1,approx="dt")
{
  m <- match.call(expand.dots = FALSE)
  m$start.time <- m$max.time <-m$id <-
    m$clusters <- m$deltaweight<-m$approx<-NULL
  special <- c("const","cluster")
  Terms <- if (missing(data)) 
    terms(formula, special)
  else terms(formula, special, data = data)
  m$formula <- Terms
  m[[1]] <- as.name("model.frame")
  m <- eval(m, sys.parent())
  mt <- attr(m, "terms")
  intercept <- attr(mt, "intercept")
  Y <- model.extract(m, "response")
  if (!inherits(Y, "Surv")) 
    stop("Response must be a survival object")

  des<-read.design(m,Terms)
  X<-des$X; Z<-des$Z; npar<-des$npar; px<-des$px; pz<-des$pz;
  covnamesX<-des$covnamesX; covnamesZ<-des$covnamesZ

  if(is.null(clusters)) clusters <- des$clusters  
  
  pxz <- px + pz;

  if (approx=="death-times" & des$npar==FALSE) npar<-TRUE; 
  survs<-read.surv(m,id,npar,clusters,start.time,max.time)
  times<-survs$times;id<-id.call<-survs$id.cal;
  if (approx=="death-times" & des$npar==FALSE) npar<-FALSE; 

  id<-id.call<-survs$id; 
  clusters<-cluster.call<-survs$clusters; 
  time2 <- survs$stop;  status <- survs$status; 
  start = survs$start; stop = survs$stop;
  antpers = survs$antpers; antclust = survs$antclust;

  Ntimes<-length(times)
  NDtimes<-sum(status[(stop>times[1]) & (stop<=times[Ntimes])])+1;
  ng=nx=antpers; 
  if (pz==0) Z=0; 

  desret<-matrix(0,antpers*(Ntimes-1),(pz+4)); 

  semiout<-.C("aalendesL",
              as.double(times),as.integer(Ntimes),as.integer(NDtimes),
              as.double(X),as.integer(nx),as.integer(px),
              as.double(Z),as.integer(ng),as.integer(pz),
              as.integer(antpers),as.double(start),as.double(stop),
              as.integer(status),as.integer(id),as.integer(clusters),
              as.integer(antclust),as.integer(deltaweight),as.double(desret),
              PACKAGE="timereg") 
  desret<-matrix(semiout[[18]],antpers*(Ntimes-1),pz+4)
  ud<-apply(desret[,1:pz],1,sum)!=0
  desret<-desret[ud,]
  colnames(desret)<-c(covnamesZ,"Y","dtime","id","time")
  return(data.frame(desret))
}
