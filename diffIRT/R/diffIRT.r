 ########################################################################################
 #  diffIRT.R
 #
 #
 #  Created by Dylan Molenaar at the University of Amsterdam.
 #  Copyright 2013 Dylan Molenaar.
 #
 #  This file is part of the diffIRT package for R.
 #
 #  The diffIRT package is free software: you can redistribute it and/or modify
 #  it under the terms of the GNU General Public License as published by
 #  the Free Software Foundation, either version 2 of the License, or
 #  (at your option) any later version.
 #
 #  This program is distributed in the hope that it will be useful,
 #  but WITHOUT ANY WARRANTY; without even the implied warranty of
 #  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 #  GNU General Public License for more details.
 #
 #  You should have received a copy of the GNU General Public License
 #  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 #
 ########################################################################################
 
diffIRT=function(rt, x, model="D", constrain=NULL,start=NULL, se=F, control=list()){
  x=as.matrix(x)
  rt=as.matrix(rt)
  if(sum((dim(rt)==dim(x))*1)!=2)
    stop("The matrix of response times should be of same size as the matrix of responses\n")
  N=nrow(rt)
  nit=ncol(rt)
  if(toupper(model)!="D" & toupper(model)!="Q" )
      stop("The 'model' argument should be either 'D' for the D-diffusion IRT model or 'Q' for the Q-diffusion IRT model.")
   if(!is.null(start)){ 
     if(length(start)!=(3*nit+2))
      stop("If you provide starting values, this vector should be of length 3 x [number of items] + 2.\n")
     if(sum(1*na.omit(start[1:nit])>0)<length(na.omit(start[1:nit])))
      stop("Starting values for 'a[i]' should be strictly positive")
     if(sum(1*na.omit(start[(nit+1):(2*nit)])>0)<length(na.omit(start[(nit+1):(2*nit)])) & model=="Q")
      stop("Starting values for 'v[i]' should be strictly positive")
     if(sum(1*na.omit(start[(2*nit+1):(3*nit)])>0)<length(na.omit(start[(2*nit+1):(3*nit)])))
      stop("Starting values for 'Ter[i]' should be strictly positive")
     if(sum(1*na.omit(start[(3*nit+1):(3*nit+2)])>=0)<length(na.omit(start[(3*nit+1):(3*nit+2)])))
      stop("Starting values for 'omega[gamma]' and 'omega[theta]' should be positive")
  }
    if(!is.logical(se)) stop("the 'se' argument should be of type logical")
    
    if(is.character(constrain)){
      if(tolower(constrain)=="ai.equal") constrain=c(rep(1,nit),2:(nit+1),(nit+2):(2*nit+1),(2*nit+2),(2*nit+3))
      if(tolower(constrain[1])=="vi.equal") constrain=c(1:nit,rep(nit+1,nit),(nit+2):(2*nit+1),(2*nit+2),(2*nit+3))
      if(tolower(constrain[1])=="ter.equal") constrain=c(1:nit,(nit+1):(2*nit),rep(2*nit+1,nit),(2*nit+2),(2*nit+3))  
      if(tolower(constrain[1])=="scale.equal") constrain=c(1:nit,(nit+1):(2*nit),(2*nit+1):(3*nit),(3*nit+1),(3*nit+1))     
    }
  if(!is.null(constrain)){
    if(length(constrain)!=3*nit+2)
      stop("Your vector of constraints should be of length 3 x [number of items] + 2.\n")
    if(length(which(constrain==0))!=0){
      if(is.null(start))
        stop("When you provide fixed parameters in 'constrain', you should provide their values at the corresponding elements in 'start'")
      if(length(intersect(which(constrain==0), which(!is.na(start))))!=length(which(constrain==0)))
        stop("When you provide fixed parameters in 'constrain', you should provide their values at the corresponding elements in 'start'")
    }
    if(min(constrain[constrain!=0])!= 1)
       stop("Parameter numbering in 'constrain' should start with '1'")
    if(length(intersect(unique(constrain[constrain!=0]),1:max(constrain)))!= max(constrain))
        stop("Parameter numbering in 'constrain' is not valid. See 'help' for examples on parameter numbering")
    }
   if(toupper(model)=="Q") model=2  else model=1
   
   if(is.null(constrain)) constrain=1:(3*nit+2)
   s0=get.start(rt,x,constrain,start,model)
  if(!is.null(start)){
    start[c(1:nit,(2*nit+1):(3*nit+2))]=log(start[c(1:nit,(2*nit+1):(3*nit+2))])
    if(model==2) start[(nit+1):(2*nit)]=log(start[(nit+1):(2*nit)])
    s0[!is.na(start)]=start[!is.na(start)]
  }
  s0[c(3*nit+1,3*nit+2)]= start[c(3*nit+1,3*nit+2)]=s0[c(3*nit+1,3*nit+2)]*2
  if(!is.na(start[3*nit+1])) if(exp(start[3*nit+1])==0) start[3*nit+1]=-746
  if(!is.na(start[3*nit+2])) if(exp(start[3*nit+2])==0) start[3*nit+2]=-746
  
  npars=length(unique(constrain[constrain!=0]))
  s0=s0[constrain!=0]
  pars=s0[!duplicated(constrain[constrain!=0])]
  
  con=list(nq=c(7),method="BFGS",eps=.01,delta=1e-7,diagnostic=F,trace=0,fnscale=1,parscale=rep(1,npars),maxit=1999,reltol=sqrt(.Machine$double.eps))
  con[names(control)]=control
  if(length(con$nq)==0)
    stop("Invalid number of quadrature points provided")
  
  if(length(complete.cases(rt))!=N | length(complete.cases(x))!=N ) stop("One or more of the subjects has missings on all items.")
  rt[is.na(rt)]=-999
  rt[is.na(x)]=-999
  x[is.na(x)]=-999
  cat("Fitting", if(model==1) "the D-diffusion IRT model." else "the Q-diffusion IRT model.", "Please wait . . .","\n")
  
  hist=matrix(,length(con$nq)+1,npars+4)
  hist[1,1:(npars)]=pars
  sel=(1:(3*nit+2))[!duplicated(constrain)]
  if((0 %in% constrain)) sel=sel[sel!=which(constrain==0)[1]]
  if(model==1){ nms=c(paste(paste("a*[",1:nit,sep=""),"]",sep=""),
    paste(paste("v[",1:nit,sep=""),"]",sep=""),
    paste(paste("ter*[",1:nit,sep=""),"]",sep=""),"omega2A*","omega2V*")
    nms=c(nms[sel],"LL","converg","func","gradient")
  }
  if(model==2){ nms=c(paste(paste("a*[",1:nit,sep=""),"]",sep=""),
    paste(paste("v*[",1:nit,sep=""),"]",sep=""),
    paste(paste("ter*[",1:nit,sep=""),"]",sep=""),"omega2A*","omega2V*")
    nms=c(nms[sel],"LL","converg","func","gradient")
  }
  colnames(hist)=nms
  rownames(hist)=c("start",paste("nq",con$nq,sep=""))
    
  for(i in 1:length(con$nq)){
    nq=con$nq[i]
    gauss.quad(nq,kind="hermite")->Quad
    Quad$weights/sqrt(pi)->W
    Quad$nodes*sqrt(2)->A
    hess=F
    if(i==length(con$nq) & se==T) hess=T
    res=optim(pars,fn=likelihood,gr=deriv,rt=rt,x=x,N=N,nit=nit,model=model,A=A,W=W,nq=nq,eps=con$eps,delta=con$delta,start=start,constrain=constrain,method=con$method,
        control=list(trace=con$trace,fnscale=con$fnscale,parscale=con$parscale,maxit=con$maxit,reltol=con$reltol),hessian=hess)
    if(res$convergence!=0) stop("Convergence problems. Optim error code: ",res$convergence,"\n")
    else pars=res$par
    hist[i+1,1:npars]=pars
    hist[i+1,npars+1]=res$value
    hist[i+1,npars+2]=res$convergence
    hist[i+1,npars+3]=res$counts[1]
    hist[i+1,npars+4]=res$counts[2]
  }

  subjLL=calcSubjLL(pars,rt=rt,x=x,N=N,nit=nit,model=model,A=A,W=W,nq=nq,eps=con$eps,delta=con$delta,start=start,constrain=constrain)
  if(con$diagnostic){ 
    derivs=deriv(pars,rt=rt,x=x,N=N,nit=nit,model=model,A=A,W=W,nq=nq,eps=con$eps,delta=con$delta,start=start,constrain=constrain)
    cat("\nvector of first-order derivatives at parameter estimates:\n")
    print(derivs)
  } else derivs=NULL
  nr_fail=sum((subjLL==1e-300)*1)
  if(nr_fail!=0) cat("Warning: for",nr_fail,"subjects, the likelihood is numerically intractable. \nThis may have resulted in an instable solution. 
    Likelihood ratio tests cannot be trusted anymore. \nYou might consider using different starting values or you might check your data.\n")
  par=rep(999,3*nit+2)
  par[constrain!=0]=pars[constrain]
  par[constrain==0]=start[constrain==0]
  par.log=par
  par[c(1:nit,(2*nit+1):(3*nit+2))]=exp(par[c(1:nit,(2*nit+1):(3*nit+2))])
  if(model==2) par[(nit+1):(2*nit)]=exp(par.log[(nit+1):(2*nit)])
  par[c(3*nit+1,3*nit+2)]=sqrt(par[c(3*nit+1,3*nit+2)])
  if(se){
    var.log=diag(solve(.5*res$hessian))
    var=var.log
    hessian=res$hessian
    var[c(1:nit,(2*nit+1):(3*nit+2))]= par[c(1:nit,(2*nit+1):(3*nit+2))]^2 * var.log[c(1:nit,(2*nit+1):(3*nit+2))]
    if(model==2) var[(nit+1):(2*nit)]= par[(nit+1):(2*nit)]^2 * var.log[(nit+1):(2*nit)]
    var[c(3*nit+1,3*nit+2)]=.5/par[c(3*nit+1,3*nit+2)]*var[c(3*nit+1,3*nit+2)]
    std.err=rep(-999,3*nit+2)
    std.err[constrain!=0]=sqrt(var[constrain])
  } else { 
        std.err=NULL
        hessian=NULL
        }
  
  AIC=2*npars+res$value
  BIC=res$value+npars*log(N)
  sBIC=res$value+npars*log((N+2)/24)
  DIC=res$value+npars*log((N+2)/(2*pi))
  rt[rt==-999]=NA
  x[x==-999]=NA
  diffRES=list(model=model,N=N,nit=nit,rt=rt,score=x,start.val=s0[!duplicated(constrain[constrain!=0])],par=par,par.log=par.log,se=se,
  std.err=std.err,hessian=hessian,totLL=res$value, npars=npars,AIC=AIC,BIC=BIC,sBIC=sBIC,DIC=DIC,subjLL=subjLL,hist=hist,
  nq=nq,conv=res$convergence,nr_fail=nr_fail,gradient=derivs,constrain=constrain,control=con,W=W, A=A)
  class(diffRES)<-"diffIRT"
  return(diffRES)
}

