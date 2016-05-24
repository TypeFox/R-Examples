"svmpath" <-
  function(x, y, K = kernel.function(x, x,param.kernel=param.kernel), kernel.function=poly.kernel,param.kernel=1, trace = FALSE, plot.it = FALSE,eps = 
           1e-10, Nmoves = 3 * n, digits=6,lambda.min=1e-4,ridge=0,...)
{
### Function to compute the entire SVM path of solutions for a two-class
### classification problem
### Copyright Trevor Hastie, May 2003
###
### Initializations
### y must be -1 and 1
###
  this.call<-match.call()
  linear.plot=!is.null(attr(K,"linear"))
  if(plot.it&&(ncol(x) > 2))stop("Plotting only for 2-dim X")
  n <- length(y)
  yvals <- table(y)
  if(length(yvals) != 2)
    stop("SvmPath works with binary problems only")
  nplus <- yvals[2]
  nminus <- yvals[1]	### Initialize the sets of observations
  Right <- Elbow <- NULL
  Left <- seq(n)
  Kscript <- K * outer(y, y)	
### We start with a maximum of 2*n moves, but these can be increased
###
### Initializations of counters
  alpha <- matrix(1, n, Nmoves)
  alpha0 <- double(Nmoves)
  SumEps <- double(Nmoves)
  Elbow.list<-as.list(seq(Nmoves))
  Size.Elbow<-integer(Nmoves)
  Error<-integer(Nmoves)
  Step<-integer(2*Nmoves)
  Obs.step<-integer(2*Nmoves)
  Movefrom<-character(2*Nmoves)
  Moveto<-character(2*Nmoves)
  lambda <- double(Nmoves)
  Kstar<-matrix(0,1,1)
###Initialization of path
### Two cases nplus=nminus or else not
  if(nplus == nminus){
    init<-Balanced.Initialization(K, y,nplus,nminus)
    Elbow <- init$Elbow
    Left <- setdiff(Left, Elbow)
  }
  else{
    init<-Unbalanced.Initialization(K, y, nplus,nminus)
    Elbow<-init$Elbow
    Right<-init$Right
    Left<-init$Left
  }
  Elbow.list[[1]]<-Elbow
  lambda0 <- init$lambda
  alpha0[1] <- init$alpha0
  alpha[,1] <-init$alpha
  alpha00<-init$alpha00 #safekeeping
  Kstar<-UpdateKstar(Kstar,Kscript[Elbow,Elbow],NULL,y[Elbow])
  lambda[1] <- lambda0
  fl <- (K %*% (alpha[, 1] * y) +init$alpha0)/lambda0
  stats<-StatPath(y,fl,Elbow)
  SumEps[1]<-stats$margin
  Error[1]<-stats$error
  Size.Elbow[1]<-stats$selbow
  nobs<-seq(along=Elbow)
  Step[nobs]<-1
  Obs.step[nobs]<-Elbow
  Movefrom[nobs]<-" ";Moveto[nobs]<-"E"
  move.counter<-length(nobs)
  PrintPath(trace,1,Elbow,"E"," ",lambda0,digits,stats)
  k <- 1
  if(plot.it)
    SnapPath(k,x, y, fl, alpha[, k], alpha0[k], lambda[k], Elbow,kernel.function,param.kernel,linear.plot,...)
  while((k < Nmoves)&& (lambda[k]>lambda.min)) {
### Now we implement the updates in Section 4.0
    if(length(Elbow)==0){
### The elbow has become empty; need to resort to an initial condition
      if(sum(y[Left])!=0)stop("Unbalanced data in interior empty elbow situation")
      init<-Balanced.Initialization(K[Left,Left],y[Left],length(Left)/2)
      lambda0 <- init$lambda
      alpha0[k+1] <- init$alpha0
      Elbow <- Left[init$Elbow]
      Left <- setdiff(Left, Elbow)
      lambda[k+1] <- lambda0
      alpha[,k+1]<-alpha[,k]
      Kstar<-UpdateKstar(Kstar,Kscript[Elbow,Elbow],NULL,y[Elbow])
      fl<- (lambda[k]/lambda[k + 1]) * (fl + (alpha0[k+1]-alpha0[k])/lambda[k])
      movefrom<-" ";moveto<-"E";obs<-Elbow
    }
    else{
      bstar<-SolveKstar(Kstar,ridge=ridge,...)
      b0 <- bstar[1]
      b <- bstar[-1]
### Now find the first event
### Check for immobile margin
      
      gl <- K[, Elbow, drop=FALSE] %*% (y[Elbow] * b) + b0
      dl<-fl-gl
      immobile<-sum(abs(dl))/n < eps
### now check for exits from Elbow
      temp <-  - alpha[Elbow, k] + lambda[k] * b
      lambda.left<-(1+temp)/b
      lambda.left[abs(b)<eps]<- -1 #anything negative
      lambda.right<-temp/b
      lambda.right[abs(b)<eps]<- -1
      lambda01 <- c(lambda.right,lambda.left)
      lambda.exit <- max(lambda01[lambda01 < lambda[k] - eps])
### Check to see if we leave the margin when it is immobile
      if(immobile&(lambda.exit < eps))break
### Now check for entries

      if(!immobile){
      lambdai <- (lambda[k] * (dl))/(y - gl)
      lambdai[abs(y-gl)<eps]<- -Inf
      lambda.entry <- max(lambdai[lambdai < lambda[k] - eps])
    }
    else lambda.entry<- -1 #any negative will do
      lambda.max <- max(lambda.entry, lambda.exit)	
### update lambda, alphas and fit
      lambda[k + 1] <- lambda.max
      alpha[, k + 1] <- alpha[, k]
      alpha[Elbow, k + 1] <- alpha[Elbow, k] - (lambda[k] - 
                                                lambda.max) * b
      alpha0[k + 1] <- alpha0[k] - (lambda[k] - lambda[k + 1]) * b0
      fl <- (lambda[k]/lambda[k + 1]) * (dl) + gl	
### update active sets
      if(lambda.entry > lambda.exit) {
        
###point joins the elbow
        i <- match(lambda.entry, lambdai, 0)[1]
        obs<-i
###assumes for now there is only 1
        moveto<-"E"
        if(match(i, Left, FALSE)) {
          Left <- setdiff(Left, i)
          movefrom<-"L"
        }
        else {
          Right <- setdiff(Right, i)
          movefrom<-"R"
        }
        Kstar<-UpdateKstar(Kstar,Kscript[i, i],drop(Kscript[i, Elbow]),y[i])
        Elbow <- c(Elbow, i)
      }
      else {
###point(s) leaves the elbow; can be more than one
        movefrom<-"E"
        moveto<-NULL
        idrop<-Leaveright<-NULL
        i<-Elbow[abs(lambda.right-lambda.exit)< eps]
        if(length(i)>0){
          Leaveright<-rep(TRUE,length(i))
          idrop<-i
        }
        i<-Elbow[abs(lambda.left-lambda.exit)< eps]
        if(length(i)>0){
          Leaveright<-c(Leaveright,rep(FALSE,length(i)))
          idrop<-c(idrop,i)
        }
        obs<-idrop
        for(j in seq(along=idrop)){
          if(Leaveright[j]) {
            moveto<-c(moveto,"R")
            Right <- c(Right, idrop[j])
          }
          else {
            moveto<-c(moveto,"L")
            Left <- c(Left, idrop[j])
          }
          mi<-match(idrop[j],Elbow)
          Kstar<-DowndateKstar(Kstar,mi)
          Elbow <- Elbow[-mi]
        }
      }
    }
    k <- k + 1
    stats<-StatPath(y,fl,Elbow)
    SumEps[k]<-stats$margin
    Error[k]<-stats$error
    Size.Elbow[k]<-stats$selbow
    nobs<-seq(along=obs)
    Moveto[move.counter+nobs]<-moveto
    Movefrom[move.counter+nobs]<-movefrom
    Step[move.counter+nobs]<-k
    Obs.step[move.counter+nobs]<-obs
    move.counter<-move.counter+length(nobs)
    Elbow.list[[k]]<-Elbow
    PrintPath(trace,k,obs,moveto,movefrom,lambda[k],digits,stats)
    if(plot.it)
    SnapPath(k,x, y,fl, alpha[, k], alpha0[k], lambda[k], Elbow, kernel.function,param.kernel, linear.plot,...)
  }
  obj<-list(alpha=alpha[,seq(k),drop=FALSE],alpha0=alpha0[seq(k)],lambda=lambda[seq(k)],alpha00=alpha00,Error=Error[seq(k)],SumEps=SumEps[seq(k)],
            Size.Elbow=Size.Elbow[seq(k)],Elbow=Elbow.list[seq(k)],Moveto=Moveto[seq(move.counter)],Movefrom=Movefrom[seq(move.counter)],Obs.step=Obs.step[seq(move.counter)],Step=Step[seq(move.counter)],kernel=kernel.function,param.kernel=param.kernel,x=x,y=y,linear=linear.plot,call=this.call)
class(obj)<-"svmpath"
  obj
}

