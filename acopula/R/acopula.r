# ------------------------------------------ function repository

# --- utility functions (move to utils.r afterwards)

#numerical (linear) approximation of partial derivatives of a real function
#fun can have arbitrary number of arguments, as well as derivatives can by of arbitrary order
#'order' is a sequence of the same length as is the number of arguments
#'difference' controls accuracy, 'area' allows for one-sided derivatives
#Example: nderive(function(x,y) x^2*y,c(5,11),c(2,0)) #22
#fun arguments can be sequence as well as vector 
nderive <- function(fun,point=c(0),order=c(1),difference=1e-4,area=c(0,1,-1)) {
  diffup <- difference*(1+sign(area[1]))/2; difflo <- difference*(1-sign(area[1]))/2
  ind <- t(expand.grid(lapply(order,function(x) seq(0,x,1))))
  arg <- point + (order-ind)*diffup - ind*difflo 
  if(length(formals(fun)) == length(point)) { #treatment of arguments in sequence e.g. fun(x,y) instead of fun(c(x,y))
    rownames(arg)<-NULL #due to unwanted passing of argument names in do.call
    func <- function(x) do.call(fun,as.list(x)) 
  } 
  else func <- fun
  sum(
    (-1)^apply(ind,2,sum) * 
      apply(choose(order,ind),2,prod) * 
      apply(arg,2,func)
    ) / (diffup+difflo)^sum(order)
}

#numerical integration of arbitrarily dimensional function
#arguments can form either vector or sequence
#possible to reduce integral with some bounds being equal
nintegrate <- function(fun, lower, upper, subdivisions=100, differences=(upper-lower)/subdivisions) {
  argseq <- mapply(seq.int,from=lower,to=upper,by=differences,SIMPLIFY=F,USE.NAMES=F) #make sequence
  #argseq <- mapply(function(x,y) unique(c(x,y)),argseq,upper,SIMPLIFY=F,USE.NAMES=F) #add upper limits to the sequence and remove if double
  differences <- rep(differences,length.out=length(argseq)) #ensure differences has the proper length
  indargseq <- lapply(argseq,seq_along) #indices of the sequence values
  nind <- sapply(indargseq,length) #number of grid nodes in each direction
  argsgrid <- as.matrix(expand.grid(argseq)); dimnames(argsgrid) <- NULL #grid nodes
  indargsgrid <- as.matrix(expand.grid(indargseq)); dimnames(indargsgrid) <- NULL #and their indices
  mult <- apply(indargsgrid,1,function(x) 2^sum(x>1 & x<nind)) #multiple of each node
  if(length(formals(fun)) > 1) func <- function(x) do.call(fun,as.list(x)) #treatment of arguments in sequence e.g. fun(x,y) instead of fun(c(x,y))
  else func <- fun
  fvalues <- apply(argsgrid,1,func) # compute function values in each node
  prod(differences[nind>1])*sum(fvalues*mult)/2^length(argseq[nind>1]) #final magic
}

#splits vector x to subvectors of specified lengths; create list or simplify to matrix if possible (identical lengths)
vpartition <- function(x, lengths, matrixify=TRUE) {
  n <- length(lengths)
  up <- cumsum(lengths)
  lo <- c(0,up[-n])
  sapply(mapply(function(i,j) c(x[1],x)[i:j], lo+1, up+1, SIMPLIFY = F),function(y) y[-1],simplify=matrixify)#treats zero lengths
}

#CDF of 4-parametric Pareto distribution (type IV), pars[1:3] > 0, location pars[4] in R
pPareto <- function(t,pars) ifelse(t>=pars[4], 1-(1+((t-pars[4])/pars[1])^(1/pars[3]))^(-pars[2]), 0) 
qPareto <- function(t,pars) pars[1]*((1-t)^(-1/pars[2])-1)^(pars[3])+pars[4]

# --- copula-related functions

#create function from 'base' and feed with 'data' (both being matrix or data frame)
pCopulaEmpirical <- function(data,base=data) {
  data <- rbind(data, deparse.level=0 ) #ensure data is matrix/data.frame
  fCe <- function(x) sum(apply(t(base)<=x,2,prod))
  apply(data,1,fCe)/nrow(base)
}

genLog <- function(...) {
  output <- list(
    parameters = NULL,
    pcopula = function(t, pars) prod(t),
    dcopula = function(t, pars) 1,
    rcopula = function(dim, pars) runif(dim),
    gen = function(t, pars) -log(t),
    gen.der = function(t, pars) -1/t,
    gen.der2 = function(t, pars) 1/t^2,
    gen.inv = function(t, pars) exp(-t),
    gen.inv.der = function(t, pars) -exp(-t),
    gen.inv.der2 = function(t, pars) exp(-t),
    lower = NULL,   
    upper = NULL,
    id="log"
  )
  output[names(list(...))] <- list(...) #allow user to modify(=replace) definition
  output
}

genGumbel <- function(...) {
  output <- list(
    parameters = c(4),
    pcopula = function(t, pars) exp(-sum((-log(t))^pars[1])^(1/pars[1])),
    gen = function(t, pars) (-log(t))^pars[1],
    gen.der = function(t, pars) -pars[1]*(-log(t))^(pars[1]-1)/t,
    gen.der2 = function(t, pars) pars[1]*(-log(t))^(pars[1]-2)*(pars[1]-1-log(t))/t^2,
    gen.inv = function(t, pars) exp(-t^(1/pars[1])),
    gen.inv.der = function(t, pars) -exp(-t^(1/pars[1]))*t^(1/pars[1]-1)/pars[1],
    gen.inv.der2 = function(t, pars) exp(-t^(1/pars[1]))*t^(1/pars[1]-2)*(pars[1]+t^(1/pars[1])-1)/pars[1]^2,
    kendall = list(coef=function(t) 1-1/t, icoef=function(t) 1/(1-t), bounds=c(0,1)),
    spearman = list(coef=function(t) pPareto(t,c(1.41917,2.14723,1,1)), icoef=function(t) qPareto(t,c(1.41917,2.14723,1,1)), bounds=c(0,1)),
    lower = 1,     #Pi, g(t)=-ln(t)
    upper = Inf,		#M
    id="Gumbel"
  )
  output[names(list(...))] <- list(...)
  output
}

#to fix: g(u)/(g(u)+g(v)) sends NaN  #edit: still true?
genClayton <- function(...) {
  output <- list(
    parameters = c(2),
    gen = function(t, pars) t^(-pars[1]) - 1,
    gen.der = function(t,pars) -pars[1]*t^(-pars[1] - 1),
    gen.der2 = function(t, pars) pars[1]*(1+pars[1])*t^(-2-pars[1]),
    gen.inv = function(t, pars) (1+t)^(-1/pars[1]),
    gen.inv.der = function(t, pars) -(1+t)^(-1-1/pars[1])/pars[1],
    gen.inv.der2 = function(t, pars) (1 + 1/pars[1])*(1 + t)^(-2 - 1/pars[1])/pars[1],
    kendall = list(coef=function(t) t/(t+2), icoef=function(t) 2*t/(1-t), bounds=c(0,1)),
    spearman = list(coef=function(t) pPareto(t,c(0.496534,1.95277,0.986814,0)), icoef=function(t) qPareto(t,c(0.496534,1.95277,0.986814,0)), bounds=c(0,1)),
    lower = 0.01, upper = Inf,
    id="Clayton"
  )
  output[names(list(...))] <- list(...)
  output
}

genFrank <- function(...) {
  output <- list(
    parameters = c(5),
    gen = function(t, pars) if(abs(pars[1]) < 0.00001) return(-log(t)) else -log( (exp(-pars[1]*t)-1)/(exp(-pars[1])-1) ),
    gen.der = function(t, pars) if(abs(pars[1]) < 0.00001) return(-1/t) else pars[1]*exp(-pars[1]*t)/(exp(-pars[1]*t)-1),
    gen.der2 = function(t, pars) if(abs(pars[1]) < 0.00001) return(1/t^2) else pars[1]^2*exp(pars[1]*t)/(1-exp(pars[1]*t))^2,
    gen.inv = function(t, pars) if(abs(pars[1]) < 0.00001) return(exp(-t)) else -log(1-exp(-t)+exp(-t-pars[1]))/pars[1],
    gen.inv.der = function(t, pars) if(abs(pars[1]) < 0.00001) return(-exp(-t)) else -exp(-t)/(-exp(-t)+1/(1-exp(-pars[1])))/pars[1],
    gen.inv.der2 = function(t, pars) if(abs(pars[1]) < 0.00001) return(exp(-t)) else (exp(pars[1]+t)*(-1+exp(pars[1])))/((1-exp(pars[1])+exp(pars[1]+t))^2*pars[1]),
    kendall = list(coef=function(t) pPareto(t,c(7.46846,1.25305,0.854724,0)), icoef=function(t) qPareto(t,c(7.46846,1.25305,0.854724,0)), bounds=c(0,1)), #only positive dep.
    spearman = list(coef=function(t) pPareto(t,c(13.2333,3.67989,0.857562,0)), icoef=function(t) qPareto(t,c(13.2333,3.67989,0.857562,0)), bounds=c(0,1)), #only positive dep.
    lower = -Inf, upper = Inf,
    id="Frank"
  )
  output[names(list(...))] <- list(...)
  output
}

genJoe <- function(...) {
  output <- list(
    parameters = c(3),
    gen = function(t, pars) -log(1-(1-t)^pars[1]),
    gen.der = function(t, pars) -pars[1]*(1-t)^(pars[1]-1)/(1-(1-t)^pars[1]),
    gen.der2 = function(t, pars) pars[1]*(pars[1]-1+(1-t)^pars[1])*(1-t)^(pars[1]-2)/(1-(1-t)^pars[1])^2,
    gen.inv = function(t, pars) 1-(1-exp(-t))^(1/pars[1]),
    gen.inv.der = function(t, pars) (1 - exp(-t))^(1/pars[1])/(pars[1]*(1 - exp(t))),
    gen.inv.der2 = function(t, pars) -(1-exp(-t))^(1/pars[1])*(1-pars[1]*exp(t))/(pars[1]*(1-exp(t)))^2,
    kendall = list(coef=function(t) pPareto(t,c(1.901055,1.015722,1.038055,1)), icoef=function(t) qPareto(t,c(1.901055,1.015722,1.038055,1)), bounds=c(0,1)), #only positive dep.
    spearman = list(coef=function(t) pPareto(t,c(2.42955,1.99806,1.02288,1)), icoef=function(t) qPareto(t,c(2.42955,1.99806,1.02288,1)), bounds=c(0,1)), #only positive dep.
    lower = 1,     #Pi, g(t)=-ln(t)
    upper = Inf,   #M
    id="Joe"
  )
  output[names(list(...))] <- list(...)
  output
}

genAMH <- function(...) {
  output <- list(
    parameters = c(0.1),
    gen = function(t, pars) log((1-(1-t)*pars[1])/t),
    gen.der = function(t, pars) (pars[1]-1)/(t*(1-(1-t)*pars[1])),
    gen.der2 = function(t, pars) (pars[1]-1)*(1+pars[1]*(2*t-1))/(t*(1-(1-t)*pars[1]))^2,
    gen.inv = function(t, pars) (1-pars[1])/(exp(t)-pars[1]),
    gen.inv.der = function(t, pars) -exp(t)*(1-pars[1])/(exp(t)-pars[1])^2,
    gen.inv.der2 = function(t, pars) exp(t)*(exp(t)+pars[1])*(1-pars[1])/(exp(t)-pars[1])^3,
    kendall = list(coef=function(t) (3*t-2)/(3*t)-(2*(1-t)^2*log(1-t))/(3*t^2), icoef=NULL, bounds=c(-0.1817258,1/3)), 
    spearman = list(coef=function(t) 0.641086*(-1 + 1.71131^t), icoef=function(t) log(t/0.641086+1,1.71131), bounds=c(-0.271065,0.478418)),
    lower = -1,   #par=0 => Pi   
    upper = 0.9999,
    id="AMH"
  )
  output[names(list(...))] <- list(...)
  output
}

generator <- function(name,...) {
  switch(name[1],
         "AMH"=genAMH(...),
         "Clayton"=genClayton(...),
         "Frank"=genFrank(...),
         "Gumbel"=genGumbel(...),
         "Joe"=genJoe(...),
         "log"=genLog(...),
         ({warning("Generator not recognized. Defaults to genLog().", immediate.=T);
           genLog(...)})
         )
}

#upper bound of all Pickands' dependence functions
dep1 <- function(...) {
  output <- list(
    parameters = NULL,
    dep = function(t,pars) 1,
    dep.der = function(t,pars) 0,
    dep.der2 = function(t,pars) 0,
    lower = NULL,
    upper = NULL,
    id="1"
  )
  output[names(list(...))] <- list(...)
  output
}

#lower bound of all Pickands' dependence functions
depMax <- function(power=10,...) {
  output <- list(
    parameters = NULL,
    dep = function(t,pars=NULL) (sum(t^power) + (1 - sum(t))^power)^(1/power),
    dep.der = function(t,pars=NULL) (t^(power-1)-(1-t)^(power-1))*(t^power + (1 - t)^power)^(1/power-1),
    dep.der2 = function(t,pars=NULL) (power-1)*((1-t)*t)^(power-2)*(t^power + (1 - t)^power)^(1/power-2),
    lower = NULL,
    upper = NULL,
    id="max"
  )
  output[names(list(...))] <- list(...)
  output
}

#depMax <- list(parameters = NULL,
#  dep = function(t,pars) max(t,1-sum(t)),
#  lower = NULL,
#  upper = NULL
#)

depGumbel <- function(...) {
  output <- list(
    parameters = c(2),
    pcopula = function(t, pars) exp(-sum((-log(t))^pars[1])^(1/pars)),
    dcopula = NULL,
    rcopula = NULL,
    dep = function(t,pars) (sum(t^pars[1]) + (1 - sum(t))^pars[1])^(1/pars[1]),
    dep.der = function(t,pars) (t^(pars[1]-1)-(1-t)^(pars[1]-1)) * (t^pars[1]+(1-t)^pars[1])^(1/pars[1]-1),
    dep.der2 = function(t,pars) (pars[1]-1)*(-(t-1)*t)^(pars[1]-2)*(t^pars[1] + (1 - t)^pars[1])^(1/pars[1]-2),
    kendall = list(coef=function(t) 1-1/t, icoef=function(t) 1/(1-t), bounds=c(0,1)),
    spearman = list(coef=function(t) pPareto(t,c(1.41917,2.14723,1,1)), icoef=function(t) qPareto(t,c(1.41917,2.14723,1,1)), bounds=c(0,1)),
    lower = c(1), #A=1, upper bound of all Pickands' dependence functions  
    upper = Inf,   #A=max(t1,t2,...,1-t1-t2-...), lower bound of all Pickands' dependence functions
    id="Gumbel"
  )
  output[names(list(...))] <- list(...)
  output
}

#so far only 2D
depGalambos <- function(...) {
  output <- list(
    parameters = c(0.5), 
    dep = function(t,pars) {
      #check for boundary arguments of A (prevents NAN)
      if(any(
        c(sapply(1:(length(t)+1),function(x) sum(c(t,1-sum(t))[-x]) > 1 )),
        pars[1] < 1e-5
      )
      ) return(1)
      1 - (sum(t^-pars[1]) + (1 - sum(t))^-pars[1])^(-1/pars[1])
    },
    dep.der = function(t,pars) if(pars[1] < 1e-5) return(0) else ((1-t)^(-1-pars[1])-t^(-1-pars[1]))/((1-t)^-pars[1]+t^-pars[1])^(1+1/pars[1]),
    dep.der2 = function(t,pars) {
      if(pars[1] < 1e-5) return(0) 
      if(abs(pars[1]-1) < 1e-5) return(2)
      ((1+pars[1])*(-(-1+t)*t)^(-2+pars[1])) * ((1-t)^-pars[1]+t^-pars[1])^(-1/pars[1])/((1-t)^pars[1]+t^pars[1])^2
    },
    kendall = list(coef=function(t) pPareto(t,c(0.579021,0.382682,0.48448,0)), icoef=function(t) qPareto(t,c(0.579021,0.382682,0.48448,0)), bounds=c(0,1)), 
    spearman = list(coef=function(t) pPareto(t,c(0.876443,0.484209,0.307277,0)), icoef=function(t) qPareto(t,c(0.876443,0.484209,0.307277,0)), bounds=c(0,1)),
    lower = c(0),    #A=1
    upper = c(10),		#A=max(t1,t2,...,1-t1-t2-...) 
    id="Galambos"
  )
  output[names(list(...))] <- list(...)
  output
}


#asymmetric logistic Picands' dependence function (Tawn, J. A. (1988). Bivariate extreme value theory: models and estimation)
#the first is the dependence parameter, all other are the asymmetry parameters (which when equal, results in symmetric logistic depfu)
depTawn <- function(dim=2,...) {
  parameters <- c(2,rep.int(0.5,dim))
  dep <- function(t,pars) {
    1 - pars[-1]%*%c(t,1-sum(t)) + sum((pars[-1]*c(t,1-sum(t)))^pars[1])^(1/pars[1])
  }  
  dep.der <- function(t,pars) {
    pars[3] - pars[2] + (1/pars[1])*((pars[2]*t)^pars[1]+(pars[3]*(1-t))^pars[1])^(1/pars[1]-1)*
      (pars[2]^pars[1]*pars[1]*t^(pars[1]-1) - pars[3]^pars[1]*pars[1]*(1-t)^(pars[1]-1))
  }
  dep.der2 <- function(t,pars) {
    (pars[2]*pars[3])^pars[1] * (pars[1]-1) * ((1-t)*t)^(pars[1]-2) * ((pars[2]*t)^pars[1]+(pars[3]*(1-t))^pars[1])^(1/pars[1]-2)
  }
  lower <- c(1,rep.int(0,dim))  #A=1, indep.
  upper <- c(Inf,rep.int(1,dim))	#A=max(t1,t2,...,1-t1-t2-...), perf.pos.dep
  output <- list(
    parameters=parameters,dep=dep,
    dep.der=if(dim==2) dep.der else NULL, dep.der2=if(dim==2) dep.der2 else NULL,
    lower=lower,upper=upper,
    id="Tawn"
  )
  output[names(list(...))] <- list(...)
  output
}


#2D only
depHuslerReiss <- function(...) {
  output <- list(
    parameters = c(0.5),
    dep = function(t,pars) {
      if(pars[1]==0) return(1)
      t*pnorm(1/pars[1] + pars[1]/2*log(t/(1-t))) + (1-t)*pnorm(1/pars[1] - pars[1]/2*log(t/(1-t)))
    },
    dep.der = function(t,pars) {
      if(pars[1]==0) return(0)
      pnorm(1/pars[1]+pars[1]/2*log(t/(1-t))) + pnorm(-1/pars[1]+pars[1]/2*log(t/(1-t))) - 1
    },
    dep.der2 = function(t,pars) {
      if(is.finite(t) && (t <= 0 || t >= 1)) return(0)
      exp(-(4+pars[1]^4*log(t/(1-t))^2)/(8*pars[1]^2)) * pars[1] * sqrt(t/(2*pi*(1-t))) / (2*(1-t)*t^2)
    },
    kendall = list(coef=function(t) pPareto(t,c(0.799473,0.25111,0.298499,0)), icoef=function(t) qPareto(t,c(0.799473,0.25111,0.298499,0)), bounds=c(0,1)), 
    spearman = list(coef=function(t) pPareto(t,c(0.877543,0.485643,0.307806,0)), icoef=function(t) qPareto(t,c(0.876443,0.484209,0.307277,0)), bounds=c(0,1)),
    lower = c(0),  #A=1, indep.
    upper = Inf,  	#A=max(t,1-t), perf.pos.dep
    id="Husler-Reiss"
  )
  output[names(list(...))] <- list(...)
  output
}

#general convex combination of dependence functions
#this function creates object (list) to be used in construction of archimax copula
#allows to include also parameters of constituting dependence functions (lined up before the combination parameters)
#the dimension of random vector need to be specified
#combination parameters are ordered variable-wise (i.e. first dim parameters are related to the first dependence function)
#with symmetry=TRUE the function depCC is called.
#requires: vpartition()
depGCC <- function(depfun=list(dep1(),depGumbel()),
                   dparameters=lapply(depfun,function(x) rep(list(NULL),max(1,length(x$parameters)))),
                   dim=2,symmetry=FALSE) {
  if(symmetry) return(depCC(depfun=depfun,dparameters=dparameters,dim=dim))
  dparameters <- lapply(dparameters,unlist)
  nd <- length(depfun)
  A <- lapply(depfun,function(x) x$dep) #extract dependence functions
  Ad <- lapply(depfun,function(x) x$dep.der) 
  Add <- lapply(depfun,function(x) x$dep.der2)
  np <- mapply(function(dep,par) {if(is.null(par)) length(dep$parameters) else 0}, depfun, dparameters) #count depfus' parameters that are to be estimated
  iD <- vpartition(1:sum(np),np); iC <- 1:(nd*dim)+sum(np) #indices of the depfus'(as list) and the combination parameters (as vector)
  dep <- function(t, pars) {
    pD <- mapply(function(static,idynamic) c(static,pars[idynamic]), dparameters, iD, SIMPLIFY=F,USE.NAMES=F) #extract depfus' pars
    pC <- combpars(pars); i0 <- pC[["i0"]]; pC <- pC[["par"]] #treat zero rows/cols, normalize (colsums=1), extract indices of depfus kept
    tpars <- matrix(c(t,1-sum(t)),byrow=T,nrow=nrow(pC[i0, ,drop=F]),ncol=dim) * pC[i0, ,drop=F] #t * pars (normalized and trimed)
    i00 <- as.logical(apply(tpars,1,sum)) # treat the zero rows occurence
    sum(mapply(
      function(PiDF,par,arg) sum(arg)*PiDF(arg[-dim]/sum(arg),par), 
      A[i0][i00], #remove those depfus and their parameters, for which the combination parameters (and also the multiplication with arguments) were zero (zero rows)
      pD[i0][i00], 
      lapply(apply(tpars[i00, ,drop=F],1,list),unlist) #remove zero rows, isolate the remaining rows
      )) 
  }
  dep.der <- function(t, pars) { 
    pD <- mapply(function(static,idynamic) c(static,pars[idynamic]), dparameters, iD, SIMPLIFY=F,USE.NAMES=F) #extract depfus' pars
    pC <- combpars(pars); i0 <- pC[["i0"]]; pC <- pC[["par"]][i0, ,drop=F] #treat zero rows/cols, normalize (colsums=1), extract indices of depfus kept
    tpars <- matrix(c(t,1-sum(t)),byrow=T,nrow=nrow(pC),ncol=dim)*pC
    i00 <- as.logical(apply(tpars,1,sum)) # treat the zero rows occurence
    sum(mapply(
      function(PiDF,PiDFd,par,arg,pc) {(pc[1]-pc[2])*PiDF(arg[1]/sum(arg),par) + pc[1]*pc[2]/sum(arg)*PiDFd(arg[1]/sum(arg),par)}, 
      A[i0][i00], #remove those depfus and their parameters, for which the combination parameters (and also the multiplication with arguments) were zero (zero rows)
      Ad[i0][i00],
      pD[i0][i00], 
      lapply(apply(tpars[i00, ,drop=F],1,list),unlist), #remove zero rows, isolate the remaining rows
      lapply(apply(pC[i00, ,drop=F],1,list),unlist)
    ))    
  }
  dep.der2 <- function(t, pars) { 
    pD <- mapply(function(static,idynamic) c(static,pars[idynamic]), dparameters, iD, SIMPLIFY=F,USE.NAMES=F) #extract depfus' pars
    pC <- combpars(pars); i0 <- pC[["i0"]]; pC <- pC[["par"]][i0, ,drop=F] #treat zero rows/cols, normalize (colsums=1), extract indices of depfus kept
    tpars <- matrix(c(t,1-sum(t)),byrow=T,nrow=nrow(pC),ncol=dim)*pC
    i00 <- as.logical(apply(tpars,1,sum)) # treat the zero rows occurence
    sum(mapply(
      function(PiDFdd,par,arg,pc) (pc[1]*pc[2])^2/sum(arg)^3*PiDFdd(arg[1]/sum(arg),par), 
      Add[i0][i00],#remove those depfus and their parameters, for which the combination parameters (and also the multiplication with arguments) were zero (zero rows)
      pD[i0][i00], 
      lapply(apply(tpars[i00, ,drop=F],1,list),unlist), #remove zero rows, isolate the remaining rows
      lapply(apply(pC[i00, ,drop=F],1,list),unlist)
    ))    
  }
  combpars <- function(pars) { #true combination parameters piled to matrix (#rows = #dep.functions) and indicators of nonzero rows(i0)/columns(j0)
    pC <- matrix(pars[iC],ncol=dim,nrow=nd,byrow=T) #extract combination parameters, pile them to matrix {p_ji} i=1...dim j=1...length(dep)
    if(all(pC==0)) pC <- pC + 1 # if all comb.parameters are 0, change them to 1  
    i0 <- as.logical(apply(pC,1,sum));  # find non-zero rows
    j0 <- as.logical(apply(pC,2,sum)); pC[i0,!j0] <- 1 # find and fill zero columns with a non-zero constant (e.g. 1) = treat zeros as equal weights
    pCn <- apply(pC,2,function(x) x/sum(x)); dim(pCn) <- dim(pC) #normalize the parameters so that column sums = 1, prevent reducing dimension
    list(par=pCn, i0=i0, j0=j0)
  }
  rescalepars <- function(pars,names=TRUE) {
    idep <- 0:(min(iC)-1)
    result <- c(dep=pars[idep],comb=t(combpars(pars)$par))
    if(!isTRUE(names)) names(result) <- NULL
    result
  }
  stpar <- function(lo,up) {
    lo[is.infinite(lo)] <- pmin(up-1,-11)[is.infinite(lo)]
    up[is.infinite(up)] <- pmax(lo+1, 11)[is.infinite(up)]
    (lo+up)/2
  }
  lower <- c(unlist(mapply(function(x,y) x$lower[0:y], depfun, np)),rep.int(0,nd*dim)) 
  upper <- c(unlist(mapply(function(x,y) x$upper[0:y], depfun, np)),rep.int(1,nd*dim))
  parameters <- stpar(lower,upper)
  list(
    parameters=parameters,dep=dep,
    dep.der=if(dim==2) dep.der else NULL, dep.der2=if(dim==2) dep.der2 else NULL,
    lower=lower,upper=upper,combpars=combpars,rescalepars=rescalepars,id="gcc"
  )  
}

#convex sum of dependence functions
#this function creates object (list) to be used in construction of archimax copula
#allows to include also parameters of constituting dependence functions (lined up before the combination parameters)
#the dimension of random vector need to be specified
#requires: vpartition()
depCC <- function(depfun=list(dep1(),depGumbel()),dparameters=lapply(depfun,function(x) rep(list(NULL),max(1,length(x$parameters)))),dim=2) {
  dparameters <- lapply(dparameters,unlist)
  nd <- length(depfun)
  A <- lapply(depfun,function(x) x$dep) #extract dependence functions
  Ad <- lapply(depfun,function(x) x$dep.der) 
  Add <- lapply(depfun,function(x) x$dep.der2)
  np <- mapply(function(dep,par) {if(is.null(par)) length(dep$parameters) else 0}, depfun, dparameters) #count depfus' parameters that are to be estimated
  iD <- vpartition(1:sum(np),np); iC <- 1:(nd)+sum(np) #indices of the depfus'(as list) and the combination parameters (as vector)
  dep <- function(t, pars) {
    pD <- mapply(function(static,idynamic) c(static,pars[idynamic]), dparameters, iD, SIMPLIFY=F,USE.NAMES=F) #extract depfus' pars
    pC <- combpars(pars); i0 <- pC[["i0"]]; pC <- pC[["par"]] #extract, treat zero vector, normalize (sum=1), get indices of nonzero values
    sum(mapply(function(PiDF,par) PiDF(t,par), A[i0], pD[i0]) * pC[i0])
  }
  dep.der <- function(t, pars) { 
    pD <- mapply(function(static,idynamic) c(static,pars[idynamic]), dparameters, iD, SIMPLIFY=F,USE.NAMES=F) #extract depfus' pars
    pC <- combpars(pars); i0 <- pC[["i0"]]; pC <- pC[["par"]] #extract, treat zero vector, normalize (sum=1), get indices of nonzero values
    sum(mapply(function(PiDFd,par) PiDFd(t,par), Ad[i0], pD[i0]) * pC[i0])
  }
  dep.der2 <- function(t, pars) { 
    pD <- mapply(function(static,idynamic) c(static,pars[idynamic]), dparameters, iD, SIMPLIFY=F,USE.NAMES=F) #extract depfus' pars
    pC <- combpars(pars); i0 <- pC[["i0"]]; pC <- pC[["par"]] #extract, treat zero vector, normalize (sum=1), get indices of nonzero values
    sum(mapply(function(PiDFdd,par) PiDFdd(t,par), Add[i0], pD[i0]) * pC[i0])
  }
  combpars <- function(pars) { #true combination parameters
    pC <- pars[iC] #extract combination parameters
    if(all(pC==0)) pC <- pC + 1 # if all comb.parameters are 0, change them to 1
    i0 <- (pC > 0) # find non-zero (and non-negative) values
    list(par=pC/sum(pC),i0=i0)
  }
  rescalepars <- function(pars,names=TRUE) {
    idep <- 0:(min(iC)-1)
    result <- c(dep=pars[idep],comb=t(combpars(pars)$par))
    if(!isTRUE(names)) names(result) <- NULL
    result
  }
  stpar <- function(lo,up) {
    lo[is.infinite(lo)] <- pmin(up-1,-11)[is.infinite(lo)]
    up[is.infinite(up)] <- pmax(lo+1, 11)[is.infinite(up)]
    (lo+up)/2
  }
  lower <- c(unlist(mapply(function(x,y) x$lower[0:y], depfun, np)),rep.int(0,nd)) 
  upper <- c(unlist(mapply(function(x,y) x$upper[0:y], depfun, np)),rep.int(1,nd))
  parameters <- stpar(lower,upper)
  list(
    parameters=parameters,dep=dep,
    dep.der=if(dim==2) dep.der else NULL, dep.der2=if(dim==2) dep.der2 else NULL,
    lower=lower,upper=upper,combpars=combpars,rescalepars=rescalepars,id="cc"
  )  
}

#return dependece function list or (unnamed) list of dependence functions list
depfun <- function(name,...) {
  ldep <- lapply(
    name,
    function(x) switch(x,
                       "1"=dep1(...),
                       "Galambos"=depGalambos(...),
                       "Gumbel"=depGumbel(...),
                       "Husler-Reiss"=depHuslerReiss(...),
                       "max"=depMax(...),
                       "Tawn"=depTawn(...),
                       "gcc"=depGCC(...),
                       "cc"=depCC(...),
                       ({warning("Dependence function not recognized. Defaults to dep1().", immediate.=T); 
                        dep1()})
                       )
    )
  if(length(ldep)==1) ldep[[1]] else ldep
}

#is there a bug when apply(somematrix,1,function(x) x/sum(x)) ??

#lapply/mapply are not able to return correct list of functions (contains only copies of the last function)
#e.g. mapply(function(a,b) c(function(c) a*b+c), list(1,2,3),list(10,100,1000))

#list of 5 dependence functions corresponding to all possible partitions of vector {1,2,3}, where 3 is number of dimensions
ldepPartition3D <- function(power=8) {
  dmax <- if(is.infinite(power)) {function(t) max(t)} else {function(t) sum(t^power)^(1/power)}
  mapply(
    function(x,y) list(parameters=NULL,dep=x,lower=NULL,upper=NULL,id=y), 
    list(
      function(t,pars) 1, #P={{1},{2},{3}}
      function(t,pars) dmax(c(t,1-sum(t))), #P={{1,2,3}}
      function(t,pars) dmax(c(t[-1],1-sum(t)))+t[1], #P={{1},{2,3}}
      function(t,pars) dmax(c(t[-2],1-sum(t)))+t[2], #P={{2},{1,3}}
      function(t,pars) dmax(t)+1-sum(t) #P={{1,2},{3}}
      ),
    list("1","max","max+1","max+2","max+3"),
    SIMPLIFY=FALSE
    )
}

copProduct <- function(...) {
  output <- list(
    parameters = NULL,
    pcopula = function(t, pars) prod(t),
    dcopula = function(t, pars) 1,
    rcopula = function(dim, pars) runif(dim),
    lower = NULL,   
    upper = NULL,
    id="product"
  )
  output[names(list(...))] <- list(...)
  output
}

copGumbel <- function(...) {
  output <- list(
    parameters = c(4),
    pcopula = function(t, pars) exp(-sum((-log(t))^pars[1])^(1/pars[1])),
    kendall = list(coef=function(t) 1-1/t, icoef=function(t) 1/(1-t), bounds=c(0,1)),
    spearman = list(coef=function(t) pPareto(t,c(1.41917,2.14723,1,1)), icoef=function(t) qPareto(t,c(1.41917,2.14723,1,1)), bounds=c(0,1)),
    lower = 1,     #Pi, g(t)=-ln(t)
    upper = Inf,  	#M
    id="Gumbel"
  )
  output[names(list(...))] <- list(...)
  output
}

copFGM <- function(...) {
  output<- list(
    parameters = c(0.5), #if 0 then product copula Pi
    pcopula = function(t, pars) prod(t)*(pars[1]*prod(1-t)+1),
    dcopula = function(t, pars) 1 + pars[1]*(prod(2*t-1)),
    kendall = list(coef=function(t) 2*t/9, icoef=function(t) 9*t/2, bounds=c(-2/9,2/9)),
    spearman = list(coef=function(t) t/3, icoef=function(t) 3*t, bounds=c(-1/3,1/3)),
    lower = -1,  #
    upper = 1,    #
    id="FGM"
  )
  output[names(list(...))] <- list(...)
  output
}

#only 2D
copPlackett <- function(...) {
  output <- list(
    parameters = c(2), #if 1 then product copula Pi
    pcopula = function(t, pars) {
      if(abs(pars[1]-1) < 0.00001) prod(t) 
      else {
        a <- 1+(pars[1]-1)*sum(t)
        (a-sqrt(a^2-4*prod(t)*pars[1]*(pars[1]-1)))/(2*(pars[1]-1))
      }
    },
    dcopula = function(t, pars) {
      if(abs(pars[1]-1) < 0.00001) 1 
      else pars[1]*(1+(sum(t)-2*prod(t))*(pars[1]-1))/((1+(pars[1]-1)*sum(t))^2-4*prod(t)*pars[1]*(pars[1]-1))^(3/2)
    },
    rcopula = function(n, pars) {
      u <- runif(n); v <- runif(n);
      if(abs(pars[1]-1) < 0.00001) return(cbind(u,v))
      a <- v*(1-v) 
      b <- sqrt(pars[1]*(pars[1]+4*a*u*(1-u)*(1-pars[1])^2))
      cbind(u,(2*a*(u*pars[1]^2+1-u)+pars[1]*(1-2*a)-(1-2*v)*b)/(2*pars[1]+2*a*(pars[1]-1)^2))
    },
    kendall = list(coef=function(t) pPareto(t,c(3.24135,0.538913,1.21742,1)), icoef=function(t) qPareto(t,c(3.24135,0.538913,1.21742,1)), bounds=c(0,1)),
    spearman = list(coef=function(t) (t^2-1-2*t*log(t))/(t-1)^2, icoef=NULL, bounds=c(-1,1)),
    lower = 1e-05,  #
    upper = Inf,    #
    id="Plackett"
  )
  output[names(list(...))] <- list(...)
  output
}

copNormal <- function(dim=2,...) {
  parvec2matrix <- function(parvec) {
    if(dim != (1+sqrt(8*length(parvec)+1))/2) stop("Dimension does not correspond to number of parameters") 
    m <- matrix(0,dim,dim)
    m[lower.tri(m)] <- parvec #fill the lower triangle (without diagonal) with parameters 
    m+t(m)+diag(1,dim) #mirror the lower triangle and add 1's
  }
  if(require(mvtnorm)) {
    pcopula <- function(t,pars) pmvnorm(lower=-Inf, upper=sapply(t,qnorm), corr=parvec2matrix(pars))[1]
    rcopula <- function(n,pars) pnorm(rmvnorm(n=n,sigma=parvec2matrix(pars)))
  }
  else {
    pcopula <- NULL
    rcopula <- function(n,pars) t(pnorm(t(chol(parvec2matrix(pars)))%*%replicate(n,rnorm(dim))))
  }
  npar <- factorial(dim)/factorial(dim-2)/2 #variacia bez opakovania / 2
  output <- list(
    parameters = rep.int(0.5,npar), #if 0 then product copula Pi
    pcopula = pcopula,
    dcopula = function(t, pars) {
      sig <- parvec2matrix(pars)
      v <- sapply(t, qnorm)
      exp(-0.5*t(v)%*%(solve(sig)-diag(1,length(t)))%*%v)/sqrt(det(sig))
    },
    rcopula = rcopula,
    kendall = list(coef=function(t) 2*asin(t)/pi, icoef=function(t) sin(t*pi/2), bounds=c(-1,1)),
    spearman = list(coef=function(t) 6*asin(t/2)/pi, icoef=function(t) 2*sin(t*pi/6), bounds=c(-1,1)),
    lower = rep.int(-0.999,npar),  #W (by limit)
    upper = rep.int(0.999,npar),  #M (by limit)
    id="normal"
  )
  output[names(list(...))] <- list(...)
  output
}

copula <- function(name,...) {
  switch(name[1],
         "FGM"=copFGM(...),
         "Gumbel"=copGumbel(...),
         "normal"=copNormal(...),
         "Plackett"=copPlackett(...),
         "product"=copProduct(...),
         ({warning("Copula not recognized. Defaults to genLog().", immediate.=T);
           copProduct(...)})
  )
}

#cumulative distribution function (or probability dis.fun., that's why "p")
pCopula <- function(data, generator=genGumbel(), depfun=dep1(), copula=NULL,
                    gpars=generator$parameters, dpars=depfun$parameters, pars=if(is.null(copula)) list(gpars,dpars) else copula$parameters, 
                    subdivisions=50,
                    quantile=NULL,probability=data[,quantile]) {
  data <- rbind(data, deparse.level=0 ) #ensure data is matrix/data.frame
  data[data>1] <- 1; data[data<0] <- 0 #crop data to [0,1]
  names(data) <- NULL
  dim <- ncol(data)
  eps <- .Machine$double.eps^0.5
  # --- cdf definition ---
  #Archimax copula
  if(is.null(copula)) { 
    if(!is.null(generator$pcopula) && depfun$id=="1") fun <- function(t) generator$pcopula(t,pars[[1]]) 
    else if(!is.null(depfun$pcopula) && generator$id=="log") fun <- function(t) depfun$pcopula(t,pars[[2]])
    else {
      g <- function(t) generator$gen(t, pars[[1]]) #generator
      gi <- function(t) generator$gen.inv(t, pars[[1]])
      A <- function(t) depfun$dep(t, pars[[2]]) #Pickands depedence function
      Aeq1 <- (abs(A(rep.int(1,dim-1)/dim)-1) < eps)
      fun <- function(t) {
        names(t) <- NULL
        #t <- pmin(pmax(t,0),1)
        #if( Aeq1 ) return(gi(sum(sapply(t,g)))) #reduce to archimedean
        if( abs(prod(t)-0) < eps ) return(0) #0 is annihilator
        if( abs(prod(t)-1) < eps ) return(1) #C(1,...,1)=1
        gen <- sapply(t,g)
        #if(any(!is.finite(gen))) browser()
        arg <- gen/(sum(gen)+eps) #prevent an accidental overflow caused by rounding at the last decimal place
        #if(!is.finite(A(arg[-dim]))) browser()
        if(any(arg > 1-eps, Aeq1)) cdf <- gi( sum(gen) )  #reduce to archimedean
        else cdf <- gi( sum(gen) * A(arg[-dim]) )
        #if(any(t > 1)) browser()
        if(!is.finite(cdf)) {
          #print(c(t=t,gpar=pars[[1]],dpar=pars[[2]],generator=gen,A=A(gen[-dim]/sum(gen)),cdf=cdf))
          return(0)
        }
        cdf
      }
    }
  }
  #arbitrary copula
  else { 
    pars <- unlist(pars) #ensure pars is a vector instead of a list    
    if(!is.null(copula$pcopula)) {
      fun <- function(t) {
        if( abs(prod(t)-0) < eps ) return(0) #0 is annihilator
        cdf <- copula$pcopula(t,pars)
        cdf
      }      
    }    
    else {
      pdf <- function(t) copula$dcopula(t,pars)
      if(!is.finite(pdf(lower <- rep.int(0,dim)))) lower <- 1e-4
      fun <- function(t) nintegrate(pdf,lower=lower,upper=t,subdivisions=subdivisions)    
    }
  }
  # --- evaluation ---
  # quantile (if asked for)
  if(!is.null(quantile)) { #to improve: treat explicitely if only copula density is available
    qua <- quantile[1]
    if(quantile > dim) stop("quantile index > copula dimension")
    data[,qua] <- rep(probability,length.out=nrow(data))
    if(any(data[,qua] > apply(data[,-qua,drop=F],1,min))) stop("probability > data")
    qcop <- function(v) { #v contains probability in place of the wanted quantile (qua-th element)
      fun1 <- function(t) {p <- v[qua]; v[qua] <- t; fun(v) - p}
      uniroot(fun1, interval=c(0,1))$root
    }
    apply(data,1,qcop)
  }
  #copula value (by default)
  else apply(data,1,fun)
}

# copula density, the function checks for closed form availability
dCopula <- function(data,generator=genGumbel(),depfun=dep1(),copula=NULL,
                    gpars=generator$parameters, dpars=depfun$parameters,
                    pars=if(is.null(copula)) list(gpars,dpars) else copula$parameters,
                    difference=1e-4,area=c(0),shrinkdiff=FALSE) {
  data <- rbind(data, deparse.level=0 ) #ensure data is matrix/data.frame
  # --- Archimax copula ---
  if(is.null(copula)) {
    if(!is.null(generator$dcopula) && depfun$id=="1") fun <- function(t) generator$dcopula(t,pars[[1]]) 
    else if(!is.null(depfun$dcopula) && generator$id=="log") fun <- function(t) depfun$dcopula(t,pars[[2]])
    else if(ncol(data)==2) {
      g <- function(t) generator$gen(t, pars[[1]])
      gd <- function(t) generator$gen.der(t, pars[[1]])
      gid <- function(t) generator$gen.inv.der(t, pars[[1]])
      gidd <- function(t) generator$gen.inv.der2(t, pars[[1]])
      A <- function(t) depfun$dep(t, pars[[2]])
      Ad <- function(t) depfun$dep.der(t, pars[[2]])
      Add <- function(t) depfun$dep.der2(t, pars[[2]])
      if( abs(A(0.5)-1) < 1e-5 )  fun <- function(t) gd(t[1])*gidd(g(t[1])+g(t[2]))*gd(t[2])
      else 
        fun <- function(t) {
          g1 <- g(t[1]); g2 <- g(t[2]); arg <- g1/(g1+g2)
          gd(t[1]) * (
            gidd((g1+g2)*A(arg))*(A(arg)-arg*Ad(arg))*(A(arg)+(1-arg)*Ad(arg)) - arg*(1-arg)*gid((g1+g2)*A(arg))*Add(arg)/(g1+g2)
            ) * gd(t[2])
        }
    }
    else {
      cdf <- function(t) pCopula(rbind(t),generator=generator,depfun=depfun,pars=pars,gpars=gpars,dpars=dpars)
      if((k=min(1-max(data),min(data))) <= difference/(abs(area)+1) && shrinkdiff) {
        difference=k*9/10 #treat leaking of numeric derivative over [0,1] 
        warning("dCopula: Difference in nderive decreased due to [0,1] overflow.",call.=FALSE)
      }
      fun <- function(t) nderive(cdf,point=t,order=rep.int(1,length(t)),difference=difference,area=area)
    }
  }
  # --- arbitrary copula ---
  else { 
    pars <- unlist(pars)
    if(!is.null(copula$dcopula)) fun <- function(t) copula$dcopula(t,pars) #take explicit formula for density if it exists
    else {
      cdf <- function(t) copula$pcopula(t,pars)
      if((min(1-max(data),min(data))) <= difference/(abs(area)+1)) {
        cdf <- function(t) copula$pcopula(pmin.int(pmax.int(t,0),1),pars) #collapse surroundings to boundary (instead of shrinking difference) #test and use above if good
      }
      fun <- function(t) nderive(cdf,point=t,order=rep.int(1,length(t)),difference=difference,area=area)
    }
  }
  # --- evaluation ---
  apply(data,1,fun)
}

cCopula <- function(data, conditional.on=c(1), generator=genGumbel(), depfun=dep1(), copula=NULL,
                    gpars=generator$parameters, dpars=depfun$parameters, pars=if(is.null(copula)) list(gpars,dpars) else copula$parameters,
                    difference=1e-4,area=c(0),
                    quantile=NULL,probability=data[,quantile]) {
  data <- rbind(data, deparse.level=0 ) #ensure data is matrix/data.frame
  dim <- ncol(data)
  if(max(conditional.on)>dim) stop("conditional.on > ncol(data)")
  con <- conditional.on
  eps <- .Machine$double.eps^0.5
  # --- Archimax copula ---
  if(is.null(copula)) {
    if(dim==2) {
      g <- function(t) generator$gen(t, pars[[1]])
      gd <- function(t) generator$gen.der(t, pars[[1]])
      gid <- function(t) generator$gen.inv.der(t, pars[[1]])
      A <- function(t) depfun$dep(t, pars[[2]])
      Ad <- function(t) depfun$dep.der(t, pars[[2]])
      ccop <- function(v) { #t is unknown quantile, v is the rest of arguments (with probability)
        if( min(abs(v)) < eps ) return(0) #0 is annihilator, prevent unboudedness of strict generator
        g1 <- g(v[1]); g2 <- g(v[2]); arg <- g1/(g1+g2)
        gid((g1+g2)*A(arg)) * gd((2-con)*v[1]+(con-1)*v[2]) * (A(arg)+Ad(arg)*(con-1-arg))          
      }
    }
    else {
      cdf <- function(t) pCopula(t,generator=generator,depfun=depfun,pars=pars,gpars=gpars,dpars=dpars)
      if((min(1-max(data),min(data))) <= difference/(abs(area)+1)) 
        cdf <- function(t) pCopula(pmin.int(pmax.int(t,0),1),generator=generator,depfun=depfun,pars=pars,gpars=gpars,dpars=dpars) #collapse surroundings to boundary (instead of shrinking difference) #test and use above if good
      ord <- rep.int(0,dim); ord[con] <- 1 #1 marks variable w.r.t. which the function will be differentiated
      dcop <- function(v) nderive(cdf,point=v,order=ord,difference=difference,area=area)
      ccop <- function(v) dcop(v)/local({v[setdiff(1:dim,con)] <- 1; dcop(v)})
    }
  }
  # --- arbitrary copula ---
  else { #to tidy up: if arbitrary copula will not contain dim=2 special case, join the 2<dimensional ccop-s with Archimax case
    pars <- unlist(pars)
    cdf <- function(t) copula$pcopula(t,pars)
    if((min(1-max(data),min(data))) <= difference/(abs(area)+1)) {
      cdf <- function(t) copula$pcopula(pmin.int(pmax.int(t,0),1),pars) #collapse surroundings to boundary (instead of shrinking difference) #test and use above if good
    }
    ord <- rep.int(0,dim); ord[con] <- 1 #1 marks variable w.r.t. which the function will be differentiated
    dcop <- function(v) nderive(cdf,point=v,order=ord,difference=difference,area=area)
    ccop <- function(v) dcop(v)/local({v[setdiff(1:dim,con)] <- 1; dcop(v)})
  }
  # --- evaluation ---
  if(!is.null(quantile)) {
    qua <- quantile
    if(qua %in% con) stop("quantile must NOT be an element of conditional.on")
    if(quantile > dim) stop("quantile index > copula dimension")
    data[,qua] <- rep(probability,length.out=nrow(data))
    qcop <- function(v) { #v contains probability in place of the wanted quantile (qua-th element)
      fun <- function(t) {p <- v[qua]; v[qua] <- t; ccop(v) - p}
      uniroot(fun, interval=c(0,1))$root
    }
    apply(data,1,qcop)
  }
  else apply(data,1,ccop)
}

#quantile of a copula (un/conditional) distribution function
qCopula <- function(data, quantile=1, probability=0.95, conditional.on=NULL, 
                    generator=genGumbel(), depfun=dep1(), copula=NULL, 
                    gpars=generator$parameters, dpars=depfun$parameters, pars=if(is.null(copula)) list(gpars,dpars) else copula$parameters,
                    difference=1e-4, area=c(0)) {
  data <- rbind(data, deparse.level=0 ) #ensure data is matrix/data.frame
  dim <- ncol(data) + 1
  qua <- quantile
  data <- cbind(data[,0:(qua-1)],probability,data[,if(qua==dim) 0 else qua:(dim-1)],deparse.level=0)
  if(is.null(conditional.on)) pCopula(data=data,quantile=quantile,generator=generator,depfun=depfun,copula=copula,pars=pars,gpars=gpars,dpars=dpars)
  else cCopula(data=data,conditional.on=conditional.on,quantile=quantile,generator=generator,depfun=depfun,copula=copula,pars=pars,gpars=gpars,dpars=dpars,difference=difference,area=area)
}

# - check the correctness of LS-nls, while LS-optim(nlminb) completely fails; 
# - allow to enter initial values for optimisation routines!!
# - do not use 'fnscale' in 'control' of 'optim'
# - to improve: print summary

eCopula <- function(data, generator=genGumbel(), depfun=dep1(), copula=NULL,
                    glimits=list(generator$lower,generator$upper), 
                    dlimits=list(depfun$lower,depfun$upper),
                    limits=list(copula$lower,copula$upper),
                    ggridparameters=if(!is.null(unlist(glimits))) do.call(function(...) mapply(c,...,length.out=pgrid,SIMPLIFY=FALSE), glimits) else NULL, 
                    dgridparameters=if(!is.null(unlist(dlimits))) do.call(function(...) mapply(c,...,length.out=pgrid,SIMPLIFY=FALSE), dlimits) else NULL, 
                    gridparameters=if(!is.null(unlist(limits))) do.call(function(...) mapply(c,...,length.out=pgrid,SIMPLIFY=FALSE), limits) else NULL,
                    technique=c("ML","LS","icorr"), procedure=c("optim","nlminb","nls","grid"), method="default", corrtype=c("kendall","spearman"), 
                    control=NULL, pgrid=10) {
  if(is.null(copula)) { 
    eCopulaArchimax(data=data,generator=generator,depfun=depfun,glimits=glimits,dlimits=dlimits,ggridparameters=ggridparameters,dgridparameters=dgridparameters,technique=technique,procedure=procedure,method=method,control=control,pgrid=pgrid,corrtype=corrtype)
  }
  else {
    eCopulaGeneric(data=data,copula=copula,limits=limits,gridparameters=gridparameters,technique=technique,procedure=procedure,method=method,control=control,pgrid=pgrid,corrtype=corrtype)
  }
}

## fitting archimax  
eCopulaArchimax <- function(data, generator, depfun=dep1(),
                            glimits=list(generator$lower,generator$upper), 
                            dlimits=list(depfun$lower,depfun$upper),
                            ggridparameters=if(!is.null(unlist(glimits))) do.call(function(...) mapply(c,...,length.out=pgrid,SIMPLIFY=FALSE), glimits) else NULL, 
                            dgridparameters=if(!is.null(unlist(dlimits))) do.call(function(...) mapply(c,...,length.out=pgrid,SIMPLIFY=FALSE), dlimits) else NULL, 
                            technique=c("ML","LS","icorr"), procedure=c("optim","nlminb","nls","grid"), method="default", corrtype=c("kendall","spearman"), 
                            control=NULL, pgrid=10) {
  #initial (local) variables
  npg <- length(generator$parameters); npd <- length(depfun$parameters)  #number of gen and dep parameters
  ig <- min(1,npg):npg; id <- {if(npd < 1) 0 else (1:npd)+npg}   #indices of generator and depfun parameters
  ind <- apply(data==1,1,prod) + apply(data==0,1,prod) # indices for removing unit and zero rows in data
  if(technique[1]=="LS") empC <- pCopulaEmpirical(data)[!ind]
  data <- data[!ind,]  
  splitad <- function(pars) list(gpars=pars[ig],dpars=pars[id])	#separate gen/dep parameters
  #temporarily simple procedure inserted  to provide basic functionality of "icorr" technique
  if(technique[1]=="icorr") { 
    if(npg+npd>1) stop("Technique icorr can currently estimate only 1-parameter copula families.")
    if(depfun$id=="1") copula <- generator
    else copula <- depfun
    if(is.null(corrlist <- copula[[corrtype[1]]])) stop("No relation of correlation coefficient to copula parameter defined.")
    coef <- cor(data,method=corrtype)[1,2]
    if((coef > corrlist$bounds[2]) || (coef < corrlist$bounds[1])) stop("Dependence measure out of bounds")
    if(is.null(corrlist$icoef)) {
      lim <- unlist(c(glimits,dlimits)); lim <- replace(lim,lim==Inf,1000); lim <- replace(lim,lim==-Inf,-1000)
      parestim <- uniroot(function(t) corrlist$coef(t)-coef,interval=lim)$root 
    }
    else parestim <- corrlist$icoef(coef)
    output <- list(parameters=splitad(parestim),approach=c(technique[1]),fvalue=NULL,procedure.output=NULL)
    class(output) <- "eCopulaArchimax"
    return(output)
  }
  #switch to fitting technique
  switch(technique[1],
         ML={
           fun <- function(pars) {
             pars <- comboe(pars)
             #truncate by "almost zero" to prevent negative values occured when numerical derivatives are used
             copuladensity <- pmax(dCopula(data, generator=generator, depfun=depfun, pars=list(pars[ig],pars[id])), exp(-50)) 
             sum(log(copuladensity))
           }
           fscale <- -1
         },
         LS={
           fun <- function(pars) {
             pars <- comboe(pars)
             sum((pCopula(data, generator=generator, depfun=depfun, pars=list(pars[ig],pars[id])) - empC)^2)
           }
           fscale <- +1
           funNLS <- function(u,pars) {
             pars <- comboe(pars)
             pCopula(u, generator=generator, depfun=depfun, pars=list(pars[ig],pars[id]))
           }
         }
  )
  #decision tree of grid estimation option 
  if(procedure[1]=="grid") {	
    #preparing parameters
    makeseq <- function(x) {	#make sequence if any argument for seq() is recognized  
      if(sum(match(names(x),c("from","to","by","length.out","along.with"),nomatch =0)) > 0) {
        x <- replace(x,x==Inf,100); x <- replace(x,x==-Inf,-100)
        return(unique(do.call(seq,as.list(x))))
      }
      else return(x)
    }
    ggridparameters <- lapply(ggridparameters,makeseq); dgridparameters <- lapply(dgridparameters,makeseq)
    ggridparameters <- mapply(function(x,y,z) x[x>=y & x<=z], ggridparameters, generator$lower, generator$upper, SIMPLIFY=FALSE)
    dgridparameters <- mapply(function(x,y,z) x[x>=y & x<=z], dgridparameters, depfun$lower, depfun$upper, SIMPLIFY=FALSE)
    parsgrid <- as.matrix(expand.grid(c(ggridparameters,dgridparameters), KEEP.OUT.ATTRS = FALSE))
    dimnames(parsgrid) <- NULL	#prevent apply from returning infinite values
    comboe <- identity #just for compatibility with (unconditional) optimization procedures
    #evaluation 
    ind <- order(fvalue <- fscale*apply(parsgrid,1,fun))
    res <- list(
      parestim=unlist(parsgrid[ind[1],], use.names = FALSE),
      fvalue=fscale*fvalue[ind[1]],
      complete=list(
        parsgrid=t(parsgrid),
        fvalue=fscale*fvalue,
        relfvalues=(max(fvalue)-fvalue)/(max(fvalue)-min(fvalue))
      )
    )
  }		
  else {
    #immediate exceptions
    if(technique[1]=="ML" && procedure[1]=="nls" ) stop("ML and nls cannot be combined")
    st <- c(generator$parameters, depfun$parameters)     #starting values (to improve: should be optional)
    lo <- c(glimits[[1]], dlimits[[1]]); up <- c(glimits[[2]], dlimits[[2]])
    if( npg+npd != length(lo) || npg+npd != length(up)) stop("differing number of parameters")
    ind <- st < lo;   st[ind] <- lo[ind]
    ipo <- (up - lo) <= 0; ipog <- ipo[ig]; ipod <- ipo[id]  #T/F index of parameters omited in estimation
    npe <- sum(!ipo)			#number of parameters that will be estimated
    ste <- st[!ipo]; loe <- lo[!ipo]; upe <- up[!ipo]   #starting and bounding values of estimated parameters
    comboe <- function(pars) {                               #combine omited and estimated parameters
      parst <- numeric(npg+npd); parst[ipo] <- lo[ipo]; parst[!ipo] <- pars; parst
    }
    #switch to fitting procedure
    switch(procedure[1],
           optim={
             res <- optim( par=st[!ipo], fn=fun, lower=lo[!ipo], upper=up[!ipo], 
                           method=ifelse(method=="default","L-BFGS-B",method), control=c(list(fnscale=fscale,factr=1e12),control) 
             )
             res <- list(parestim=res$par,fvalue=res$value,procedure.output=res)
           },
           nlminb={
             fun1 <- function(pars) fscale*fun(pars)        #check the "port" help to circumvent this line
             res <- nlminb( start=st[!ipo], objective=fun1, lower=lo[!ipo], upper=up[!ipo],control=c(list(),control));
             res <- list(parestim=res$par,fvalue=fscale*res$objective,procedure.output=res)
           },
           nls={
             nlsmethod <- ifelse(method=="default","port",method) 
             if(ncol(data)!=3) stop("nls currently available only for 3D") #tip for improvement: handle formula separately in advance
             emp <- as.data.frame(data); colnames(emp) <- c('u1','u2','u3'); emp <-cbind(empC,emp);
             res <- switch(npe,
                           nls(C ~ funNLS(c(u1, u2,u3), c(par1)),
                               data = data,
                               start = list(par1 = ste[1]),
                               lower = list(par1 = loe[1]),
                               upper = list(par1 = upe[1]),
                               algorithm = nlsmethod
                           ),
                           nls(C ~ funNLS(c(u1, u2,u3), c(par1,par2)),
                               data = emp,
                               start = list(par1 = ste[1],par2 = ste[2]),
                               lower = list(par1 = loe[1],par2 = loe[2]),
                               upper = list(par1 = upe[1],par2 = upe[2]),
                               algorithm = nlsmethod
                           ),
                           nls(C ~ funNLS(c(u1, u2,u3), c(par1,par2,par3)),
                               data = emp,
                               start = list(par1 = ste[1],par2 = ste[2],par3 = ste[3]),
                               lower = list(par1 = loe[1],par2 = loe[2],par3 = loe[3]),
                               upper = list(par1 = upe[1],par2 = upe[2],par3 = upe[3]),
                               algorithm = nlsmethod
                           ),
                           nls(C ~ funNLS(c(u1, u2,u3), c(par1,par2,par3,par4)),
                               data = emp,
                               start = list(par1 = ste[1],par2 = ste[2],par3 = ste[3],par4 = ste[4]),
                               lower = list(par1 = loe[1],par2 = loe[2],par3 = loe[3],par4 = loe[4]),
                               upper = list(par1 = upe[1],par2 = upe[2],par3 = upe[3],par4 = upe[4]),
                               algorithm = nlsmethod
                           ),
                           nls(C ~ funNLS(c(u1, u2,u3), c(par1,par2,par3,par4,par5)),
                               data = emp,
                               start = list(par1 = ste[1],par2 = ste[2],par3 = ste[3],par4 = ste[4],par5 = ste[5]),
                               lower = list(par1 = loe[1],par2 = loe[2],par3 = loe[3],par4 = loe[4],par5 = loe[5]),
                               upper = list(par1 = upe[1],par2 = upe[2],par3 = upe[3],par4 = upe[4],par5 = upe[5]),
                               algorithm = nlsmethod
                           )
             )
             res <- list(parestim=as.vector(coef(res)),fvalue=sum(residuals(res)^2),procedure.output=res)
           }
    )
  }
  parestim <- splitad(comboe(res$parestim))
  if("rescalepars" %in% names(depfun)) parestim$dpars <- depfun$rescalepars(parestim$dpars) # rescale to true parameters of depfun (if appliable)
  output <- list(parameters=parestim,approach=c(technique[1],procedure[1],method[1]),fvalue=res$fvalue,procedure.output=res)
  class(output) <- "eCopulaArchimax"
  output
}

## fitting generic copula
eCopulaGeneric <- function(data, copula=copGumbel(),
                           limits=list(copula$lower,copula$upper),
                           gridparameters=if(!is.null(unlist(limits))) do.call(function(...) mapply(c,...,length.out=pgrid,SIMPLIFY=FALSE), limits) else NULL,
                           technique=c("ML","LS","icorr"), procedure=c("optim","nlminb","grid"), method="default", corrtype=c("kendall","spearman"), 
                           control=NULL, pgrid=10) {
  #initial (local) variables
  np <- length(copula$parameters)  #number copula parameters
  ind <- apply(data==1,1,prod) + apply(data==0,1,prod) # indices for removing unit and zero rows in data
  if(technique[1]=="LS") empC <- pCopulaEmpirical(data)[!ind]
  data <- data[!ind,]
  #temporarily simple procedure inserted  to provide basic functionality of "icorr" technique
  if(technique[1]=="icorr") { 
    if(np>1) stop("Technique icorr can currently estimate only 1-parameter copula families.")
    if(is.null(corrlist <- copula[[corrtype[1]]])) stop("No relation of correlation coefficient to copula parameter defined.")
    coef <- cor(data,method=corrtype)[1,2]
    if((coef > corrlist$bounds[2]) || (coef < corrlist$bounds[1])) stop("Dependence measure out of bounds")
    if(is.null(corrlist$icoef)) {
      lim <- unlist(limits); lim <- replace(lim,lim==Inf,1000); lim <- replace(lim,lim==-Inf,-1000)
      parestim <- uniroot(function(t) corrlist$coef(t)-coef,interval=lim)$root 
    }
    else parestim <- corrlist$icoef(coef)
    output <- list(parameters=parestim,approach=c(technique[1]),fvalue=NULL,procedure.output=NULL)
    class(output) <- "eCopulaGeneric"
    return(output)
  }
  #switch to fitting technique
  switch(technique[1],
         ML={
           fun <- function(pars) {
             pars <- comboe(pars)
             #truncate by "almost zero" to prevent negative values occured when numerical derivatives are used
             copuladensity <- pmax(dCopula(data, copula=copula, pars=pars), exp(-50)) 
             sum(log(copuladensity))
           }
           fscale <- -1
           
         },
         LS={
           fun <- function(pars) {
             pars <- comboe(pars)
             sum((pCopula(data,  copula=copula, pars=pars) - empC)^2)
           }
           fscale <- +1
         }
  )
  #decision tree of grid estimation option 
  if(procedure[1]=="grid") {	
    #preparing parameters
    makeseq <- function(x) {	#make sequence if any argument for seq() is recognized  
      if(sum(match(names(x),c("from","to","by","length.out","along.with"),nomatch =0)) > 0) {
        x <- replace(x,x==Inf,100); x <- replace(x,x==-Inf,-100)
        return(unique(do.call(seq,as.list(x))))
      }
      else return(x)
    }
    gridparameters <- lapply(gridparameters,makeseq)
    gridparameters <- mapply(function(x,y,z) x[x>=y & x<=z], gridparameters, generator$lower, generator$upper, SIMPLIFY=FALSE)
    parsgrid <- as.matrix(expand.grid(gridparameters, KEEP.OUT.ATTRS = FALSE))
    dimnames(parsgrid) <- NULL	#prevent apply from returning infinite values
    comboe <- identity #just for compatibility with (unconditional) optimization procedures
    #evaluation 
    ind <- order(fvalue <- fscale*apply(parsgrid,1,fun))
    res <- list(
      parestim=unlist(parsgrid[ind[1],], use.names = FALSE),
      fvalue=fscale*fvalue[ind[1]],
      complete=list(
        parsgrid=t(parsgrid),
        fvalue=fscale*fvalue,
        relfvalues=(max(fvalue)-fvalue)/(max(fvalue)-min(fvalue))
      )
    )
  }		
  else {
    #immediate exceptions
    st <- copula$parameters     #starting values (to improve: should be optional)
    lo <- limits[[1]]; up <- limits[[2]]
    if( np != length(lo) || np != length(up)) stop("Differing number of parameters.")
    ind <- st < lo;   st[ind] <- lo[ind]
    ipo <- (up - lo) <= 0  #T/F index of parameters omited in estimation
    npe <- sum(!ipo)			#number of parameters that will be estimated
    ste <- st[!ipo]; loe <- lo[!ipo]; upe <- up[!ipo]   #starting and bounding values of estimated parameters
    comboe <- function(pars) {                          #combine omited and estimated parameters
      parst <- numeric(np); parst[ipo] <- lo[ipo]; parst[!ipo] <- pars; parst
    }
    #switch to fitting procedure
    switch(procedure[1],
           optim={
             method <- ifelse(method=="default","L-BFGS-B",method)
             if(method %in% c("Nelder-Mead", "BFGS", "CG", "SANN", "Brent")) { #for these optim methods limit the parameters range implicitely
               fun <- function(pars) {
                 parscut <- pmin.int(pmax.int(pars,loe),upe)
                 if(any(parscut != pars)) hndcp <- (1+sum(abs(pars-parscut)))^fscale else hndcp <- 1 #handicap out-of-bounds parameters
                 fun(parscut)*hndcp                 
               }
               res <- optim( par=ste, fn=fun, method=method, control=c(list(fnscale=fscale,factr=1e12),control) )
             }
             else { #if "L-BFGS-B" method is used, fun arguments will not surpass the limits during optimization
               res <- optim( par=ste, fn=fun, lower=loe, upper=upe, method=method, control=c(list(fnscale=fscale,factr=1e12),control) )
             }
             res <- list(parestim=res$par,fvalue=res$value,procedure.output=res) #compose output
           },
           nlminb={
             fun1 <- function(pars) fscale*fun(pars)        #consult "port" help to circumvent this line
             res <- nlminb( start=ste, objective=fun1, lower=loe, upper=upe,control=c(list(),control));
             res <- list(parestim=res$par,fvalue=fscale*res$objective,procedure.output=res)
           }
    )
  }
  parestim <- comboe(res$parestim)
  if("rescalepars" %in% names(copula)) parestim <- copula$rescalepars(parestim) # rescale to true parameters of copula (if appliable)
  output <- list(parameters=parestim,approach=c(technique[1],procedure[1],method[1]),fvalue=res$fvalue,procedure.output=res)
  class(output) <- "eCopulaGeneric"
  output
}


#methods for printing eCopula results list
print.eCopulaArchimax <- function(x,...) {
  cat("generator parameters: ", x$parameters$gpars, "\n")
  cat("   depfun parameters: ", x$parameters$dpars, "\n")
  cat("  ", x$approach[1], "function value: ", x$fvalue, "\n")
  cat("    convergence code: ", x$procedure.output$procedure.output$convergence, "\n")
}
print.eCopulaGeneric <- function(x,...) {
  cat("   copula parameters: ", x$parameters, "\n")
  cat("  ", x$approach[1], "function value: ", x$fvalue, "\n")
  cat("    convergence code: ", x$procedure.output$procedure.output$convergence, "\n")
}

#simulation of dim-dimensional random vector from archimax copula
#note the thinned runif range due to nderive 'difference' settings
rCopula <- function(n,dim=2,generator=genGumbel(),depfun=dep1(),copula=NULL,
                    gpars=generator$parameters, dpars=depfun$parameters, pars=if(is.null(copula)) list(gpars,dpars) else copula$parameters) {
  if(dim==2 && is.null(copula)) return( rCopulaArchimax2D(n,generator=generator,depfun=depfun,pars=pars) )
  if(!is.null(copula$rcopula)) return( copula$rcopula(n,pars) )
  cdf <- function(u) pCopula(u,generator=generator,depfun=depfun,copula=copula,pars=pars)
  mcop <- function(u) cdf(rbind(c(u,rep.int(1,dim-length(u)))))
  NDdiff <- 1e-4
  ccop <- function(v,u) nderive(fun=mcop,point=c(u,v),order=c(rep.int(1,length(u)),0),difference=NDdiff)/nderive(fun=mcop,point=u,order=rep.int(1,length(u)),difference=NDdiff)
  qcop <- function(p,u) uniroot(function(t) ccop(t,u) - p, interval=c(0,1))$root 
  p <- matrix(runif(dim*n,min=0+1.1*NDdiff,max=1-1.1*NDdiff),ncol=dim) #thin the runif range using NDdiff to avoid overflow in nderive
  result <- NULL
  for(i in 1:n) {
    x <- p[i,1]
    for(j in 2:dim) x <- c(x, qcop(p[i,j],x))
    result <- rbind(result,x,deparse.level=0)
  }
  result
}

# simulation of 2D copula
rCopulaArchimax2D <- function(n, generator=genLog(), depfun=dep1(),
                              gpars=generator$parameters, dpars=depfun$parameters, pars=list(gpars,dpars)) {
  g <- function(t) generator$gen(t, pars[[1]])
  gd <- function(t) generator$gen.der(t, pars[[1]])
  gi <- function(t) generator$gen.inv(t, pars[[1]])
  A <- function(t) depfun$dep(t, pars[[2]])
  Ad <- function(t) depfun$dep.der(t, pars[[2]])
  Add <- function(t) depfun$dep.der2(t, pars[[2]])
  K <- function(t) t - ifelse(t>0 && t<1, g(t)/gd(t), 0)
  Ki <- function(t) uniroot(f = function(x) K(x) - t, interval=c(0,1))$root 
  H <- function(t) t + ifelse(t>0 && t<1, t*(1-t)*Ad(t)/A(t), 0)
  Hi <- function(t) uniroot(f = function(x) H(x) - t, interval=c(0,1))$root 
  rCA <- function(n) {                       #simulation from Archimedean copula
    s <- runif(n); w <- runif(n)
    w <- sapply(w,Ki); gw <- sapply(w,g)
    cbind(u=sapply(s*gw,gi),v=sapply((1-s)*gw,gi))
  }
  Aeq1 <- isTRUE(all.equal(A(0.5),1))
  if( Aeq1 )  return( rCA(n) )
  p <- function(t) {
    A <- A(t); Ad <- Ad(t); Add <- Add(t); D <- Ad/A; Dd <- (Add*A - Ad^2)/A^2
    h <- 1 + (1-2*t)*D + t*(1-t)*Dd
    t*(1-t)*Add/(h*A)
  }
  z <- sapply(runif(n),Hi); u <- runif(n)
  ind <- (u <= sapply(z,p))
  w <- numeric(n); w[ind] <- runif(sum(ind)); w[!ind] <- sapply(runif(n-sum(ind)),Ki)
  gA <- sapply(w,g)/sapply(z,A)
  cbind(u=sapply(z*gA,gi), v=sapply((1-z)*gA,gi))
}


#Blanket goodness-of-fit test 
#require library(parallel) iff ncores > 1 (functionality not appliable on Windows OS)
#Example: gCopula(u123[,1:2],generator=genGumbel,dep=dep1,N=10,nc=1,m=400)
gCopula <- function(data, generator, depfun=dep1(), copula=NULL,
                    glimits=list(generator$lower,generator$upper), 
                    dlimits=list(depfun$lower,depfun$upper),
                    limits=list(copula$lower,copula$upper),
                    etechnique=c("ML","LS","icorr"), eprocedure=c("optim","nlminb","nls"), emethod="default", ecorrtype=c("kendall","spearman"), econtrol=NULL,
                    N=100, m=nrow(data), ncores=1) {  
  if(is.list(data)) return(gCopulaEmpirical(data=data,N=N,ncores=ncores))
  n <- nrow(data)
  dim <- ncol(data)
  Ein <- pCopulaEmpirical(data)
  fKn <- function(x,y) sum(y <= x)/length(y)
  fparest <- function(data) eCopula(data,generator=generator,depfun=depfun,copula=copula,glimits=glimits,dlimits=dlimits,limits=limits,
                                    technique=etechnique,procedure=eprocedure,method=emethod,corrtype=ecorrtype,control=econtrol)    
  fsimcop <- function(parameters) rCopula(m,dim=dim,generator=generator,depfun=depfun,copula=copula, pars=parameters) 
  estim <- fparest(data)
  #---2D Archimax---
  if(ncol(data)==2 && is.null(copula)) {  
    m <- n
    fK <- function(vec,pars) {
      tauA <- 0
      if(depfun$id!="1") {
        integrand <- Vectorize(function(t) t*(1-t)*depfun$dep.der2(t, pars[[2]])/depfun$dep(t, pars[[2]]), vectorize.args="t")
        tauA <- integrate(integrand,lower=0,upper=1)$value
      }
      sapply(vec, function(t) t - ifelse(t>0 && t<1, (1-tauA)*generator$gen(t, pars[[1]])/generator$gen.der(t, pars[[1]]), 0))
    }
    Sn <- sum((rank(Ein)/n - fK(Ein,estim$parameters))^2)
    pb <- txtProgressBar(min=0,max=N,style=3)
    loop <- function(k) {
      simcop <- fsimcop(estim$parameters)
      Ein <- pCopulaEmpirical(simcop)
      simcopRescaled <- apply(simcop,2,rank)/(n+1)
      parest <-  fparest(simcopRescaled)$parameters
      setTxtProgressBar(pb, k)
      #if(k%%10==0) print(k)
      sum((rank(Ein)/n - fK(Ein,parest))^2)
    }    
  }
  #---other copulas---
  else {  
    m <- max(m,n)
    simcop <- fsimcop(estim$parameters)
    Vi <- pCopulaEmpirical(simcop)
    Bm <- rank(Vi)/m
    Kn <- sapply(Vi,fKn,y=Ein)
    Sn <- (n/m)*sum((Kn-Bm)^2)
    pb <- txtProgressBar(min=0,max=N,style=3)
    loop <- function(k) {
      simcop <- fsimcop(estim$parameters)
      Ein <- pCopulaEmpirical(simcop)
      simcopRescaled <- apply(simcop,2,rank)/(m+1)
      parest <-  fparest(simcopRescaled)$parameters
      simcop <- fsimcop(parest)
      Vi <- pCopulaEmpirical(simcop)
      Bm <- rank(Vi)/m
      Kn <- sapply(Vi,fKn,y=Ein)
      setTxtProgressBar(pb, k)
      #if(k%%10==0) print(k)
      (n/m)*sum((Kn-Bm)^2)
    }
  }
  Snk <- if(ncores > 1) do.call(c, parallel::mclapply(1:N,loop,mc.cores=ncores)) else sapply(1:N,loop)
  close(pb)
  #cat("\n")  
  output <- list(
    statistic=Sn,
    q95=quantile(Snk,0.95,names=F),
    p.value=sum(Snk > Sn)/N,
    estimate=c(if(is.null(copula)) c(gpars=estim$parameters$gpars,dpars=estim$parameters$dpars) else pars=estim$parameters, fvalue=estim$fvalue),
    data.name=deparse(substitute(data)),
    method="Blanket GOF test based on Kendall's transform",
    fitting_method=as.vector(c(etechnique,eprocedure,emethod)),
    copula_id=if(is.null(copula)) c(generator=generator$id,depfun=depfun$id) else copula$id
    )
  class(output) <- "gCopula"
  output
}

gCopulaEmpirical <- function(data,N=100,ncores=1) {
  data1 <- data[[1]]; data2 <- data[[2]]
  dim <- ncol(data1); if(dim!=ncol(data2)) stop("Different number of columns.")
  n1 <- nrow(data1); n2 <- nrow(data2)
  #Cramer-von Mises test statistic
  Sn <- local({
    split1 <- split(data1,col(data1)); split2 <- split(data2,col(data2))
    s1 <- sum(Reduce("*",lapply(split1,function(a) outer(a,a,function(...) 1-pmax(...)))))
    s2 <- sum(Reduce("*",mapply(function(a,b) outer(a,b,function(...) 1-pmax(...)),split1,split2,SIMPLIFY=F)))
    s3 <- sum(Reduce("*",lapply(split2,function(a) outer(a,a,function(...) 1-pmax(...)))))
    1/(1/n1+1/n2)*(s1/n1^2-s2*2/n1/n2+s3/n2^2)
  })
  #empirical copulas and their approximate derivatives
  Cn <- function(x) sum(apply(t(data1) <= x,2,prod))/n1
  Dn <- function(x) sum(apply(t(data2) <= x,2,prod))/n2
  dCn <- function(x,i,h) {
    e <- numeric(dim); e[i] <- h
    (Cn(x+e)-Cn(x-e))/2/h
  } 
  dDn <- function(x,i,h) {
    e <- numeric(dim); e[i] <- h
    (Dn(x+e)-Dn(x-e))/2/h
  }
  #gaussian process and multiplier technique functions
  alpha <- function(x) sum(apply(t(data1) <= x,2,prod)*xi)/sqrt(n1)
  gamma <- function(x) sum(apply(t(data2) <= x,2,prod)*zeta)/sqrt(n2)
  beta <- function(x,i) sum((data1[,i]<=x)*xi)/sqrt(n1)
  delta <- function(x,i) sum((data2[,i]<=x)*zeta)/sqrt(n2)
  Ck <- function(x) alpha(x) - sum(sapply(1:dim,function(i) beta(x[i])*dCn(x,i,1/sqrt(n1))))
  Dk <- function(x) gamma(x) - sum(sapply(1:dim,function(i) delta(x[i])*dDn(x,i,1/sqrt(n2))))
  EPSk <- function(x) (sqrt(n2)*Ck(x)-sqrt(n1)*Dk(x)) #divisor sqrt(n1+n2) moved to integral
  #merged data
  data12 <- merge(data1,data2,all=T); n12 <- nrow(data12)
  #iteration
  pb <- txtProgressBar(min=0,max=N,style=3)
  xi <- zeta <- numeric(0)
  loop <- function(k) {
    xii <- rnorm(n1); xi <<- xii - mean(xii)
    zetaa <- rnorm(n2); zeta <<- zetaa - mean(zetaa)
    setTxtProgressBar(pb, k)
    sum(apply(data12,1,EPSk)^2)/(n1+n2)/n12    
  }
  Snk <- if(ncores > 1) do.call(c, parallel::mclapply(1:N,loop,mc.cores=ncores)) else sapply(1:N,loop)
  close(pb)
  #compose output
  output <- list(
    statistic=Sn,
    q95=quantile(Snk,0.95,names=F),
    p.value=sum(Snk > Sn)/N,
    estimate=NULL,
    data.name=if(is.null(names(data))) deparse(substitute(data)) else names(data),
    method="Test of equality between 2 empirical copulas (Remillard & Scaillet 2009)",
    fitting_method=NULL,
    copula_id=NULL
  )
  class(output) <- "gCopula"
  output  
}


#method for printing gCopula results list
print.gCopula <- function(x,...) {
  cat("\n\t\t", x$method, "\n\n")
  print.default(c(statistic=x$statistic,q95=x$q95,p.value=x$p.value))
  cat("-----------------------------\n")
  cat("data: ", x$data.name, "\n")
  cat("copula: ",x$copula_id,"\n")
  cat("estimates:\n")
  print.default(x$estimate)
  invisible(x)
}

#check d-increasingness, 0 as annihilator and 1 as neutral element for every parameters combination
#Example: isCopula(generator=genJoe,dep=depGalambos,dim=3)
isCopula <- function(generator=genLog(),depfun=dep1(),copula=NULL,
                     glimits=list(generator$lower,generator$upper), 
                     dlimits=list(depfun$lower,depfun$upper),
                     limits=list(copula$lower,copula$upper),  
                     ggridparameters=if(!is.null(unlist(glimits))) do.call(function(...) mapply(c,...,length.out=pgrid,SIMPLIFY=FALSE), glimits) else NULL,
                     dgridparameters=if(!is.null(unlist(dlimits))) do.call(function(...) mapply(c,...,length.out=pgrid,SIMPLIFY=FALSE), dlimits) else NULL, 
                     gridparameters=if(!is.null(unlist(limits))) do.call(function(...) mapply(c,...,length.out=pgrid,SIMPLIFY=FALSE), limits) else NULL,
                     dagrid=10, pgrid=10, dim=3, tolerance=1e-15) {
  #initial (local) variables
  npg <- length(generator$parameters); npd <- length(depfun$parameters)  #number of gen and dep parameters
  ig <- min(1,npg):npg; id <- {if(npd < 1) 0 else (1:npd)+npg}   #indices of generator and depfun parameters
  splitad <- function(pars) list(gpars=pars[ig],dpars=pars[id])  #func to separate gen/dep parameters (in the end)
  #preparing parameters
  makeseq <- function(x) {  #make sequence if any argument for seq() is recognized  
    if(sum(match(names(x),c("from","to","by","length.out","along.with"),nomatch =0)) > 0) {
      x <- replace(x,x==Inf,32); x <- replace(x,x==-Inf,-32) #replace infinite boundary
      return(unique(do.call(seq,as.list(x))))
    }
    else return(x)
  }
  gparameters <- lapply(ggridparameters,makeseq) #make sequence 
  dparameters <- lapply(dgridparameters,makeseq) 
  parameters <- lapply(gridparameters,makeseq)  
  gparameters <- mapply(function(x,y,z) x[x>=y & x<=z], gparameters, generator$lower, generator$upper, SIMPLIFY=FALSE) #ensure the parameter sequence is within permitted range
  dparameters <- mapply(function(x,y,z) x[x>=y & x<=z], dparameters, depfun$lower, depfun$upper, SIMPLIFY=FALSE)
  parameters <- mapply(function(x,y,z) x[x>=y & x<=z], parameters, copula$lower, copula$upper, SIMPLIFY=FALSE)
  #preparing data 
  axes <- mapply(seq.int,rep.int(0,dim),rep.int(1,dim),MoreArgs=list(length.out=dagrid)) #expand grid of all (n=m^d) points
  allpoints <- do.call(expand.grid,split(t(axes),1:dim))  
  #distinguish archimax and generic copula
  if(is.null(copula)) {
    parameters <- c(gpar=gparameters,dparameters=dparameters)
    cdfun <- function(pars) pCopula(allpoints, generator=generator, depfun=depfun, pars=list(pars[ig],pars[id]))
  }
  else {
    names(parameters) <- paste("cpar",1:length(parameters),sep="")
    cdfun <- function(pars) pCopula(allpoints, copula=copula, pars=pars)
  }
  #minimal dim-difference (C-volume) in hypercube data grid
  monot <- function(hcdata) {
    out <- numeric(dim)
    for(i in 1:dim) {
      sub1 <- subd <- rep.int(",",dim) #dim equals length(dim(hcdata))
      sub1[i] <- -1; subd[i] <- -dim(hcdata)[i]
      subcop <- function(sub) eval(parse(text = paste(c("hcdata[", sub, ", drop = FALSE]"),collapse="")))
      hcdata <- subcop(sub1) - subcop(subd) #recursive first differences
      out[i] <- min(hcdata) #find the smallest of the current i-th differences
    }
    out
  }
  annih <- function(hcdata) {
    sapply(
      1:dim,
      function(i) { #e.g. i=3,dim=4 yields hcdata[,,1,]
        sub <- rep.int(",",dim) #dim equals length(dim(hcdata))
        sub[i] <- 1
        max(abs(eval(parse(text = paste(c("hcdata[", sub, "]"),collapse=""))))) #find maximal deviation from 0 #max|C(x1,...xn,0,xm,...) - 0|
      }
    )
  }
  neutel <- function(hcdata) {
    sapply(
      1:dim,
      function(i) {
        sub <- as.character(dim(hcdata))
        sub[i] <- ""
        max(abs(eval(parse(text = paste(c("hcdata[", paste(sub,collapse=","), "]"),collapse="")))-axes[,i])) #max|C(1,...1,x,1,...1) - x|
      }
    )    
  }
  fun <- function(pars) {
    coparray <- cdfun(pars) 
    if(length(coparray)!=dagrid^dim || !all(is.finite(coparray))) stop(paste("Copula with parameter(s): ",paste(pars,collapse=" ")," led to errors.\n")) #test for propper length (nrows) to catch errors in cdf
    coparray <- array(coparray,rep.int(dagrid,dim))
    c(
      monot(coparray),
      -annih(coparray), #minus unary operator for evaluation compatibility with one-sided monot criterion
      -neutel(coparray)
    )
  }
  parsgrid <- as.matrix(expand.grid(parameters, KEEP.OUT.ATTRS = FALSE)) #all combinations
  result <- apply(parsgrid,1,fun); dim(result) <- c(dim,3,ncol(result)) #dim1~copuladimension,dim2~(monot,annih,neutel),dim3~parametersgrid
  ind <- which(result < 0 - tolerance, arr.ind=T, useNames=F) 
  output <- data.frame(
    dim=ind[,1], #to which variable the issue is related
    property=c("monot","annih","neutel")[ind[,2]], #which copula property is violated
    value=result[ind]*ifelse(ind[,2]==1,1,-1), #value of difference or deviation from theoretical value, that exceeds tolerance
    parsgrid[ind[,3],,drop=F] #actual parameters of copula related to the issue 
  )
  output <- list(
    is.copula = !as.logical(nrow(output)),
    issues = output
  )
  class(output) <- "isCopula"
  output
}
 
#method for printing isCopula results (does not work with output as data.frame; to fix)
print.isCopula <- function(x,...) {
  n <- nrow(x$issues)
  cat("\n Does the object appears to be a copula(?): ", x$is.copula, "\n")
  if(n>0) {
    cat("\n Showing", min(n,5), "of", n, "issues: \n\n")
    print.data.frame(x$issues[1:min(n,5),])
  }
  #invisible(x)
}