design.mat <-
function(mpars, moutput, par.reg, time.reg)
  {
    
  # PRELIMINARIES #!+
  Theta.mat <- t(mpars$par) # Parameter matrix Theta #!+
                            # [row, col] = [run index, parameter index]
  par.p     <- dim(Theta.mat)[1]      # Number of runs #!+
  par.m     <- dim(Theta.mat)[2]      # Number of parameters #!+
  par.n     <- length(moutput$t)  # Number of times #!+


  # CONSTRUCT MATRICES #!+
  # Data matrix #!+
  # Parameter index varies fastest, time index the slowest 
  Y.mat     <- matrix(t(moutput$out), par.p*par.n, 1 ) 

  # Design matrix #!+
  n1        <- as.matrix(seq(1, 1, length.out=par.n)) 
  p1        <- as.matrix(seq(1, 1, length.out=par.p)) 
  D.mat1    <- kronecker(n1, Theta.mat) 
  D.mat2    <- kronecker(moutput$t, p1) 
  D.mat     <- cbind(D.mat1, D.mat2) 

  # Regression matrix #!+
  X.mat     <- as.matrix(seq(1, 1, length.out=par.n*par.p)) 
  all.reg   <- c(par.reg, time.reg) # Logicals indicating which regressors to include 
  if (sum(all.reg) > 0) X.mat <- cbind(X.mat, D.mat[,all.reg])

  
  # OUTPUT #!+
  emul.mat  <- list(Theta.mat = Theta.mat, Y.mat = Y.mat, D.mat = D.mat,
                    X.mat = X.mat)
  emul.mat
  }
