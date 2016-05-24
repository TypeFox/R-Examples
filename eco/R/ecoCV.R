ecoX <- function(formula, Z, supplement = NULL, data = parent.frame(), 
                 nu0 = 4, S0 = 10, beta0 = 0, A0 = 100,		
                 grid = FALSE, parameter = FALSE,
                 n.draws = 5000, burnin = 0, thin = 5, verbose = TRUE){ 

  ## checking inputs
  if (burnin >= n.draws)
    stop("Error: n.draws should be larger than burnin")
  
  call <- match.call()

  ff <- as.formula(paste(call$Y, "~ -1 +", call$X))
  if (is.matrix(eval.parent(call$data)))
    data <- as.data.frame(data)
  X <- model.matrix(ff, data)
  Y <- model.response(model.frame(ff, data=data))
  
  ##survey data
  if (length(supplement) == 0) {
    survey.samp <- 0
    survey.data <- 0
    survey.yes<-0
  }
  else {
    survey.samp <- length(supplement[,1])
    survey.data <- as.matrix(supplement)
    survey.yes<-1
  }
  
  ind<-c(1:length(X))
  X1type<-0
  X0type<-0
  samp.X1<-0
  samp.X0<-0
  X1.W1<-0
  X0.W2<-0
  
  ##Xtype x=1
  X1.ind<-ind[along=(X==1)]
  if (length(X[X!=1])<length(X)){
    X1type<-1
    samp.X1<-length(X1.ind)
    X1.W1<-Y[X1.ind]
  }
  
  ##Xtype x=0
  X0.ind<-ind[along=(X==0)]
  if (length(X[X!=0])<length(X)){
    X0type<-1
    samp.X0<-length(X0.ind)
    X0.W2<-Y[X0.ind]
  }
  
  XX.ind<-setdiff(ind, union(X0.ind, X1.ind))
  X.use<-X[XX.ind]
  Y.use<-Y[XX.ind]

  order.old<-order(c(XX.ind, X0.ind, X1.ind))
  
  ## fitting the model
  n.samp <- length(Y.use)	 
  d <- cbind(X.use, Y.use)

  n.a <- floor((n.draws-burnin)/thin)
  n.par <- n.a
  n.w <- n.a * (n.samp+samp.X1+samp.X0) 
  unit.a <- 1
  unit.par <- 1
  unit.w <- (n.samp+samp.X1+samp.X0) 	
  Zmat<-Z%x%diag(1, 2)
  print(Zmat)
  Zp<-dim(Zmat)[2]
  cat("Zp", Zp)
  if (is.null(beta0)) beta0<-rep(0, Zp)
  if (is.null(A0)) A0<-diag(0.01, Zp)
  print(beta0)
  print(A0)
  n.a.b<-n.a*Zp
  n.a.V<-n.a*3
  res <- .C("cBaseecoZ", as.double(d), as.double(Zmat), as.integer(Zp),  
            as.integer(n.samp), as.integer(n.draws), as.integer(burnin), as.integer(thin),
            as.integer(verbose),
            as.integer(nu0), as.double(S0),
            as.double(beta0), as.double(A0),
            as.integer(survey.yes), as.integer(survey.samp), as.double(survey.data),
            as.integer(X1type), as.integer(samp.X1), as.double(X1.W1),
            as.integer(X0type), as.integer(samp.X0), as.double(X0.W2),
            as.integer(predict), as.integer(parameter), 
            pdSBeta=double(n.a.b),
            pdSSigma=double(n.a.V),
            pdSW1=double(n.w), pdSW2=double(n.w), 
            pdSWt1=double(n.w), pdSWt2=double(n.w), PACKAGE="eco")
  
  if (parameter) {
    beta.post <- matrix(res$pdSBeta, n.a, Zp, byrow=TRUE) 
    Sigma.post <- matrix(res$pdSSigma, n.a, 3, byrow=TRUE)
    colnames(Sigma.post) <- c("Sigma11", "Sigma12", "Sigma22")
  }
  W1.post <- matrix(res$pdSW1, n.a, unit.w, byrow=TRUE)[,order.old]
  W2.post <- matrix(res$pdSW2, n.a, unit.w, byrow=TRUE)[,order.old]
  
  res.out <- list(model="Normal prior", burnin=burnin, thin = thin, X=X, Y=Y,
                  nu0=nu0, A0=A0, beta0=beta0, S0=S0, call=call, beta.post=beta.post,
                  Sigma.post=Sigma.post, W1.post=W1.post, W2.post=W2.post)

  class(res.out) <- c("ecoCV", "eco")
  return(res.out)
}


