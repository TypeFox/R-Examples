assign("graph.rbf",
       function (formula, data, eta.opt, rho.opt, n.neigh, func, np, x0=c(0.5,0.5), eta.dmax, rho.dmax, P.T=NULL, iter, ...)                    
       {
         if (!is.logical(P.T) & !is.null(P.T))
           stop(paste("P.T must be logical or null"))
         if (!is.logical(eta.opt))
           stop(paste("eta.opt must be logical"))
         if (!is.logical(rho.opt))
           stop(paste("rho.opt must be logical"))
         
         if (eta.opt==TRUE & rho.opt==FALSE) {
           Opt <- optimize(rbf.cv, c(1e-05, eta.dmax), formula=formula, data=data, 
                           rho=0, n.neigh=n.neigh, func=func, ...)
           Datos <- as.data.frame(matrix(NA, nrow = length(seq(0.01,
                                                               eta.dmax, length.out = np)), ncol = 2))
           eta <- seq(0.01, eta.dmax, length.out = np)
           colnames(Datos) <- c("P", "RMSPE")
           pb <- txtProgressBar(min = 0, max = np, char = "=", style = 3)
           for (i in 1:np) {
             Datos[i, 1] <- eta[i]
             Datos[i, 2] <- rbf.cv(formula, data, eta[i], rho=0, n.neigh, func)
             setTxtProgressBar(pb, i)
           } 
           close(pb)
           Table0 <- rbind(Datos, c(Opt$minimum, Opt$objective))
           orden <- order(Table0$P)
           Table <- Table0[orden, ]
           plot(Table, lty = 3, ylab = "RMSPE", col = 3, xlab = "ETA",
                type = "l")
           Optim <- Table[which.min(Table[, 2]), ]
           ifelse(P.T == TRUE, list(print(Table), cat("Optimal eta RBF: ",
                                                      func, "\n", "ETA  = ", Optim$P, "\n", "RMSPE   = ", Optim$RMSPE,
                                                      "\n")), list(cat("Optimal eta RBF: ", func, "\n", "ETA  = ",
                                                                       Optim$P, "\n", "RMSPE   = ", Optim$RMSPE, "\n")))   
         }
         
         if (rho.opt==TRUE & eta.opt==FALSE) {
           Opt <- optimize(rbf.cv, c(1e-05, rho.dmax), formula=formula, data=data, 
                           eta=1e-05, n.neigh=n.neigh, func=func, ...)
           Datos <- as.data.frame(matrix(NA, nrow = length(seq(0.01,
                                                               rho.dmax, length.out = np)), ncol = 2))
           rho <- seq(0.01, rho.dmax, length.out = np)
           colnames(Datos) <- c("P", "RMSPE")
           pb <- txtProgressBar(min = 0, max = np, char = "=", style = 3)
           for (i in 1:np) {
             Datos[i, 1] <- rho[i]
             Datos[i, 2] <- rbf.cv(formula, data, eta=0, rho[i], n.neigh, func) 
             setTxtProgressBar(pb, i)
           }
           close(pb)
           Table0 <- rbind(Datos, c(Opt$minimum, Opt$objective))
           orden <- order(Table0$P)
           Table <- Table0[orden, ]
           plot(Table, lty = 3, ylab = "RMSPE", col = 3, xlab = "RHO",
                type = "l")
           Optim <- Table[which.min(Table[, 2]), ]
           ifelse(P.T == TRUE, list(print(Table), cat("Optimal eta RBF: ",
                                                      func, "\n", "RHO  = ", Optim$P, "\n", "RMSPE   = ", Optim$RMSPE,
                                                      "\n")), list(cat("Optimal eta RBF: ", func, "\n", "RHO  = ",
                                                                       Optim$P, "\n", "RMSPE   = ", Optim$RMSPE, "\n")))
         }
         
         if (eta.opt==FALSE & rho.opt==FALSE){
           eta <- seq(0.01, eta.dmax, length.out = np)
           rho <- seq(0.01, rho.dmax, length.out = np)    
           grid.opt <- expand.grid(eta = eta, rho = rho)
           grid.rmspe <- as.data.frame(matrix(NA, nrow= nrow(grid.opt), ncol=3))
           colnames(grid.rmspe) <- c("eta","rho","rmspe")
           pb <- txtProgressBar(min = 0, max = nrow(grid.opt), char = "=", style = 3)
           for (i in 1:nrow(grid.opt)){
             grid.rmspe[i,3] <- rbf.cv(formula, data, eta=grid.opt[i,1], rho=grid.opt[i,2], n.neigh, func)
             grid.rmspe[,1:2] <- grid.opt
             grid.rmspe
             setTxtProgressBar(pb, i)
           }
           close(pb)
           coordinates(grid.rmspe) = c("eta", "rho")
           opt.table <- data.frame(coordinates(grid.rmspe),grid.rmspe@data)[which.min(data.frame(coordinates(grid.rmspe),grid.rmspe@data)[,3]),]
           gridded(grid.rmspe) <- TRUE
           p <- spplot(grid.rmspe, "rmspe", col.regions=heat.colors(100), cuts=60, cex.main=0.2, scales = list(draw =TRUE), xlab=expression(eta), ylab = expression(rho), key.space=list(space="right", cex=0.6))
           list(opt.table=opt.table, spplot=p)
         }
         
         else if (eta.opt==TRUE & rho.opt==TRUE) {
           Opt <- bobyqa(x0, rbf.cv1, lower=c(1e-05,0), upper=c(eta.dmax,rho.dmax), formula=formula,
                         data=data, n.neigh=n.neigh, func=func, control = list(maxfun=iter), ...)
           cat("Optimal eta RBF: ", func, "\n", "ETA  = ", Opt$par[1], "RHO  = ", Opt$par[2], "\n", "RMSPE   = ", Opt$fval, "\n")    
           list(Opt=Opt)
         }
       }
)