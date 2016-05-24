#*****************************************************************************************************************
# esitimating an (E)CCC-GARCH(1,1) model
      eccc.estimation <- function(a, A, B, R, dvar, model, method="BFGS"){
         dvar <- as.matrix(dvar)
         nobs <- dim(dvar)[1]
         ndim <- dim(dvar)[2]
         if(!is.matrix(A)||!is.matrix(B)){
            stop("A or B or both must be matrices")
         }
         if(model=="diagonal"){
          init <- c(sqrt(a), diag(sqrt(A)), diag(sqrt(B)), R[lower.tri(R)])
         } else {
          init <- c(sqrt(a), as.vector(sqrt(A)), as.vector(sqrt(B)), R[lower.tri(R)])
         }
         results <- optim(par=init, fn=loglik.eccc, method=method, dvar=dvar, model=model, control=list(maxit=10^5, reltol=1e-15))
      
         if(results$convergence != 0){
            cat("***********************************************************\n")
            cat("* Optimization is FAILED.                                 *\n")
           stop("* Fine tuning is required.                                *\n")
            cat("***********************************************************\n")
         }
         
         estimates <- p.mat(results$par, model=model, ndim=ndim)
         esta <- estimates$a
         estA <- estimates$A
         estB <- estimates$B
         estR <- estimates$R
         
         h <- vector.garch(dvar, esta, estA, estB)    # estimated volatility
         std.resid <- dvar/sqrt(h)                          # std. residuals

         grad <- analytical.grad(a=esta, A=estA, B=estB, R=estR, u=dvar, model=model)
         grad <- grad%*%t(grad)/nobs
         H    <- analytical.Hessian(a=esta, A=estA, B=estB, R = estR, u=dvar, model=model)/nobs
         invH <- solve(H)
         Avar <- (invH%*%grad%*%invH)/nobs
         rob.se <- sqrt(diag(Avar))                # robust standard errors
         out.se <- sqrt(diag(solve(grad)/nobs))                # standard errors based on the outer product of the gradient
         H.se   <- sqrt(diag(invH/nobs))                # standard errors based on the inverted Hessian

         std.error <- rbind(H.se, out.se, rob.se); rownames(std.error) <- c("Inv. Hessian", "Outer Prod.", "Robust")

         name.a <- numeric(0)
         name.A <- numeric(0)
         name.B <- numeric(0)
         name.R <- numeric(0)
         for(i in 1:ndim){          # index for column
            name.a <- c(name.a, paste("a", paste(i), sep=""))
               for(j in 1:ndim){    # index for row
                  name.A <- c(name.A, paste("A",paste(paste(j), i, sep=""), sep=""))
                  name.B <- c(name.B, paste("B",paste(paste(j), i, sep=""), sep=""))
                  name.R <- c(name.R, paste("R",paste(paste(j), i, sep=""), sep=""))
               }
         }
         names(esta) <- name.a
         if(model=="diagonal"){
            ind <- as.vector(diag(ndim))
            vecA <- diag(estA);        name.A <- name.A[ind==1]
            vecB <- diag(estB);        name.B <- name.B[ind==1]
         } else {
            vecA <- as.vector(estA)
            vecB <- as.vector(estB)
         }
         names(vecA) <- name.A
         names(vecB) <- name.B
         vecR <- estR[lower.tri(estR)];         name.R <- matrix(name.R, ndim, ndim);  name.R <- name.R[lower.tri(name.R)]
         names(vecR) <- name.R
         colnames(std.error) <- c(name.a, name.A, name.B, name.R)
         para.estimates <- c(esta, vecA, vecB, vecR)
         
      list(out=rbind(para.estimates, std.error), h=h, std.resid=std.resid, opt=results, para.mat=estimates)
      }

