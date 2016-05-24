dcc.results <- function(u, garch.para, dcc.para, h, model){
      nobs <- dim(u)[1]
      ndim <- dim(u)[2]
      In <- diag(ndim)
      param <- c(garch.para, In[lower.tri(In)])
      estimates <- p.mat(param, model, ndim)
      a <- estimates$a
      A <- estimates$A
      B <- estimates$B
      if(model=="diagonal"){
         parameters <- c(a, diag(A), diag(B), dcc.para)
      } else {
         parameters <- c(a, as.vector(A), as.vector(B), dcc.para)
      }
      d2lvdw <- d2lv(u, B, h, model=model)                              # nobs x npar.h^2
      dlvdw <- dlv(u, a, A, B, model=model)                             # npar.h x nobs

      O <- dlc(dcc.para, B, u, h, model=model)
      dlcdf <- O$dlc                                                    # nobs x 2
      d2lcdf <- O$d2lc                                                  # nobs x 2^2
      d2lcdfdw <- O$dfdwd2lc                                            # nobs x 2*npar.h

      G11 <- matrix(colMeans(d2lvdw), ncol=sqrt(dim(d2lvdw)[2]))        # npar.h x npar.h
      G22 <- matrix(colMeans(d2lcdf), 2, 2)                             # 2 x 2
      G21 <- matrix(colMeans(d2lcdfdw), nrow=2)                         # 2 x npar.h


      O11 <- dlvdw%*%t(dlvdw)/nobs                                      # npar.h x npar.h
      O22 <- t(dlcdf)%*%dlcdf/nobs                                      # 2 x 2
      O12 <- dlvdw%*%dlcdf/nobs                                         # npar.h x 2

      G. <- cbind(rbind(G11, G21), rbind(matrix(0, nrow=dim(G11)[1], ncol=2),G22))
      O. <- cbind(rbind(O11, t(O12)), rbind(O12, O22))
      invG. <- solve(G., tol=1e-50)
      # V <- invG.%*%O.%*%invG./nobs                                       covariance matrix
      V <- invG.%*%O.%*%t(invG.)/nobs                                      # covariance matrix
      se <- sqrt(diag(V))

      results <- rbind(parameters, se)
      
      
      name.a <- numeric(0)
      name.A <- numeric(0)
      name.B <- numeric(0)
   	for(i in 1:ndim){	         # index for column
   		name.a <- c(name.a, paste("a", paste(i), sep=""))
      		for(j in 1:ndim){	   # index for row
      			name.A <- c(name.A, paste("A",paste(paste(j), i, sep=""), sep=""))
      			name.B <- c(name.B, paste("B",paste(paste(j), i, sep=""), sep=""))
      		}
   	}
      if(model=="diagonal"){
         ind <- as.vector(diag(ndim))
      	name.A <- name.A[ind==1]
         name.B <- name.B[ind==1]
      } 
      colnames(results) <- c(name.a, name.A, name.B, "dcc alpha", "dcc beta")
      rownames(results) <- c("estimates", "std.err")
   results
}
