pe.sasieni<-function (formula = formula(data),data = sys.parent(), 
id=NULL,start.time=0,max.time=NULL,offsets=0,Nit=50,detail=0,n.sim=500)
{
  call <- match.call()
  m <- match.call(expand.dots = FALSE)
  m$id<-m$Nit<-m$detail<-m$start.time <- m$max.time <-
    m$offsets<-m$n.sim <- NULL 
  Terms <- if (missing(data)) 
    terms(formula ) else terms(formula, data = data)
  m$formula <- Terms
  m[[1]] <- as.name("model.frame")
  m <- eval(m, sys.parent())
  mt <- attr(m, "terms")
  intercept <- attr(mt, "intercept")
  Y <- model.extract(m, "response")
  if (!inherits(Y, "Surv")) 
    stop("Response must be a survival object")

  XZ<-model.matrix(Terms,m)[, drop = FALSE]
  cols<-attributes(XZ)$assign
  l.cols<-length(cols)

  X<-as.matrix(XZ[,1]); covnamesX <- dimnames(XZ)[[2]][1];
  dimnames(X)[[2]]<-covnamesX; 
  Z<-as.matrix(XZ[,-1]); covnamesZ <- dimnames(XZ)[[2]][-1]; 
  px <- ncol(X); pz <- ncol(Z); pxz <- px + pz

 if ( (nrow(X)!=nrow(data)) & (!is.null(id))) stop("Missing values in design matrix not allowed with id \n"); 
###  if (nrow(Z)!=nrow(data)) stop("Missing values in design matrix not allowed\n"); 

  if (attr(m[, 1], "type") == "right") {
    X <- data.matrix(X); Z <- data.matrix(Z); 
    time2 <- m[, 1][, "time"]; time <- rep(0, length(time2))
    status <- m[, 1][, "status"]
  }
  else if (attr(m[, 1], "type") == "counting") {
    X <- data.matrix(X); Z <- data.matrix(Z) ; 
    time <- m[, 1][, 1]; time2 <- m[, 1][, 2]; status <- m[, 1][, 3]
  }
  else { stop("only right-censored or counting processes data") }

  if (sum(duplicated(time2[status==1]))>0) {
    #cat("Non unique survival times: break ties ! \n")
    # cat("Break ties yourself\n");
    ties<-TRUE
    dtimes<-time2[status==1]
    index<-(1:length(time2))[status==1]
    ties<-duplicated(dtimes); nties<-sum(ties); index<-index[ties]
    dt<-diff(sort(time2)); dt<-min(dt[dt>0]); 
    time2[index]<-time2[index]+runif(nties,0,min(0.001,dt/2));
  } else ties<-FALSE; 

  times<-unique(time2); 
  times <- c(start.time, times[times>start.time]); times <- sort(times)
  if (is.null(max.time) == TRUE) maxtimes <- max(times) else maxtimes <- max.time
  times<-times[times<=maxtimes];
  times<-c(times,maxtimes)  
  times<-unique(times); Ntimes <- length(times); 
  ldata <- list(start = time, stop =time2)

  ntot <- ncol(XZ); px <- ncol(X)
  Nalltimes <- length(times);  
  Ntimes<-sum(status[(time2>times[1]) & (time2<=times[Nalltimes])])+1;

  X<-as.matrix(X); Z<-as.matrix(Z); 
  px <- as.integer(dim(X)[2]); nx <- as.integer(dim(X)[1]);
  pg <- as.integer(dim(Z)[2]); ng <- as.integer(dim(Z)[1]);

  if (length(offsets)==1) mof<-0 else mof<-1;
  mw<-0; weights<-rep(1,nx); 

  cum<-Vcum<-matrix(0,Ntimes,px+1); 
  Ut<-matrix(0,Nalltimes,pg+1); 
  gamma<-intZHdN<-rep(0,pg); 
  Vargam<-intZHZ<-matrix(0,pg,pg); 
  antpers<-length(unique(id))
  dUt<-matrix(0,Ntimes,pg*pg);
  if (n.sim >0) {testOBS<-rep(0,pg); test<-matrix(0,n.sim,pg);}
  else {testOBS<-0; test<-0;}
  rani<- -round(runif(1)*10000)

  #dyn.load("pes.so"); 

  semiout<-.C("pes",
              as.double(times),as.integer(Nalltimes),as.integer(Ntimes),
              as.double(X),as.integer(nx),as.integer(px),
              as.double(Z),as.integer(ng),as.integer(pg),
              as.integer(antpers),as.double(time),as.double(time2), 
              as.double(cum),as.double(Vcum),as.double(gamma),
              as.double(Vargam),as.integer(status),as.double(Ut),
              as.double(intZHZ),as.double(intZHdN),as.integer(mof),
              as.double(offsets),as.integer(mw),as.double(weights),
              as.integer(Nit),as.integer(detail),
              as.integer(rani),as.integer(n.sim),as.double(test),
              PACKAGE="timereg"); 

  cum <-matrix(semiout[[13]],Ntimes,px+1); 
  Vcum <-matrix(semiout[[14]],Ntimes,px+1); 
  gamma<-matrix(semiout[[15]],pg,1); Vargam<-matrix(semiout[[16]],pg,pg); 
  intZHZ<-matrix(semiout[[19]],pg,pg); intZHdN<-matrix(semiout[[20]],pg,1); 
  Ut<-matrix(semiout[[18]],Nalltimes,pg+1); 

  #dUt<-matrix(semiout[[27]],Ntimes,pg*pg); dUt.list<-list();
  #for (i in 1:Ntimes) dUt.list[[i]]<-matrix(dUt[i,],pg,pg);

  if (n.sim>0) {test<-matrix(semiout[[29]],n.sim,pg);
                testOBS<-apply(abs(Ut),2,max)[-1]; testval<-c(); 
                for (i in 1:pg) testval<-c(testval,pval(test[,i],testOBS[i]));
                pval.Prop<-testval; 
                names(pval.Prop) <- names(testOBS)<-covnamesZ
              } else {pval.Prop<-NULL;testOBS<-NULL;}

  ud<-list(cum=cum,var.cum=Vcum,gamma=gamma,var.gamma=Vargam,
           Ut=Ut,D2linv=intZHZ,score=intZHdN,test.Prop=testOBS,pval.Prop=pval.Prop); 

  colnames(ud$cum) <- colnames(ud$var.cum) <- c("time", covnamesX)
  ud$gamma<-as.matrix(ud$gamma);

  rownames(ud$gamma) <- c(covnamesZ)
  colnames(ud$gamma) <- "estimate"
  colnames(ud$var.gamma) <- c(covnamesZ)
  rownames(ud$var.gamma) <- c(covnamesZ)
  rownames(ud$score) <- c(covnamesZ)
  colnames(ud$D2linv) <- c(covnamesZ)

  attr(ud, "Call") <- call
  attr(ud, "Formula") <- formula
  attr(ud, "start") <- start.time
  attr(ud, "time2") <- time2
  class(ud) <- "pe.sasieni"
  ud$call<-call
  return(ud)
}

summary.pe.sasieni <- function (object,digits=3,...) {
  obj <- object; rm(object);
  if (!inherits(obj, 'pe.sasieni')) 
    stop ("Must be a Proportional Excess Survival Model object based")
  if (is.null(obj$gamma)==TRUE) stop(" No proportional  terms"); 
    
  cat("\nProportional Excess Survival Model \n\n")

  cat("Proportional terms:  \n"); 
  res <- cbind(obj$gamma,diag(obj$var.gamma)^.5) 
  z<-c((res[,1]/res[,2])); 
  pval<-1-pchisq(z^2,1)
  res<-as.matrix(cbind(res,z,pval)); 
  colnames(res) <- c("coef", "se(coef)","z","p") 
  prmatrix(signif(res, digits)); cat("   \n");  

  if (is.null(obj$pval.Prop)==TRUE)  ptest<-FALSE else ptest<-TRUE;
  if (ptest==TRUE) { cat("Test for Proportionality\n"); 
                     testP<-cbind(obj$test.Prop,obj$pval.Prop); testP<-as.matrix(testP);
                     colnames(testP) <- c("sup|  hat U(t) |","p-value H_0 ")
                     prmatrix(signif(testP,digits)); cat("\n"); }

  cat("  Call: "); dput(attr(obj, "Call")); cat("\n")
}


print.pe.sasieni <- function (x,digits=3,...) {
  obj <- x; rm(x);
  if (!inherits(obj, 'pe.sasieni')) 
    stop ("Must be a Proportional Excess Survival Model object based")
  if (is.null(obj$gamma)==TRUE) stop(" No proportional  terms");

  if (is.null(obj$gamma)==TRUE) cox<-FALSE else cox<-TRUE
      
  cat(" Proportional Excess Survival Model,\n   using Sasieni proportional excess risk \n\n")
  cat(" Nonparametric terms: "); 
  cat(colnames(obj$cum)[-1]);
  cat("   \n");  
  if (cox) {
    cat(" Proportional terms:  "); 
    cat(rownames(obj$gamma)); 
    cat("   \n");
  }
  cat("   \n");    

  cat("  Call: "); dput(attr(obj, "Call")); cat("\n")
}
