random.effects <- function(object, HPDinterval=FALSE, prob=0.95)
{
				if(!inherits(object, "Bayesthresh"))
								stop("Use an object of class Bayesthresh")
				tc <- (dim(object$compVar)[1])-1
				ret <- dim(coef(object))[1]
				aleat <- NULL
				if (tc == 1){
								aleat <- data.frame(object$EfRandom, row.names=object$NomesZ)
								class(aleat) <- c("random.effects", "data.frame","conf")
								if(HPDinterval==TRUE){
												aleat <- data.frame(object$EfRandom, row.names=object$NomesZ)
												HPD <- HPDinterval(mcmc(object$outp.mcmc[[1]][,c((ret+1):(ncol(object$outp.mcmc[[1]])))]), prob=prob)
												rownames(HPD) <- c(object$Nomes)
												hpd.data <- data.frame(aleat,HPD)
												class(hpd.data) <- c("random.effects","data.frame","hpd")
												return(hpd.data)
								}
								else
												return(aleat)
				}
				else {
								res.mcmc <- list()
								gera.out <- list()
								resmat <- object$outp.mcmc[[1]][,-(1:ret)]
								il <-  1
								ifc <- object$fl
								ic <- object$fl[1]
								for(i in 2:tc){
												ic[i] <- ifc[i]+ic[i-1]
												il[i] <- ic[i]-ifc[i]+1
								}
								for(i in 1:tc){
												aleat[[i]] <- data.frame(object$EfRandom[il[i]:ic[i],], row.names=object$NomesZ[il[i]:ic[i]])
												class(aleat[[i]]) <- "data.frame"
												if(HPDinterval==TRUE)
												{
																res.mcmc[[i]] <- HPDinterval(mcmc(resmat[,il[i]:ic[i]]),prob=prob)
																gera.out[[i]] <- data.frame(aleat[[i]],res.mcmc[[i]])
																class(gera.out[[i]]) <- "data.frame"
												}
								}
								if(HPDinterval==TRUE)
								{
												names(gera.out) <- rownames(object$compVar[1:(nrow(object$compVar)-1),])
												class(gera.out) <- c("random.effects","hpd")
												return(gera.out)
								}
								else{
								names(aleat) <- rownames(object$compVar[1:(nrow(object$compVar)-1),])
								class(aleat) <- c("random.effects","conf")
								return(aleat)
								}
				}
}


