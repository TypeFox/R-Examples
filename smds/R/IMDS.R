"IMDS" <-
function(IDM, p=2,eps= 1e-5 ,maxit =1000,model=c("sphere","box"),opt.method=c("MM","BFGS"), ini = "auto",report=100,grad.num=FALSE,rel=0,dil=1){
						nvec <- dim(IDM)
						n <- nvec[2]
						if(model[1]=="sphere"){
							if(ini[1]=="auto"){
								cmat <- (IDM[2,,] + IDM[1,,])/2
								iniX <- cmdscale(as.dist(cmat),k=p)
								inir <- runif(n,min = 0.1, max = 1)
								inixr <- c(iniX,sqrt(inir))
							}else{
								iniX <- ini[[1]]
								inir <- ini[[2]]
								inixr <- c(ini[[1]],sqrt(ini[[2]]))
							}
							if(opt.method[1]=="BFGS"){
								if(grad.num==FALSE){
									result.opt <- optim(par=inixr, fn = obj.sph, gr=Grad.sph, ldiss=IDM[1,,],udiss=IDM[2,,], P=p,method ="BFGS",control = list( trace = TRUE,abstol=eps,REPORT=report,maxit=maxit))
								}else if(grad.num==TRUE){
									result.opt <- optim(par=inixr, fn = obj.sph, ldiss=IDM[1,,],udiss=IDM[2,,], P=p,method ="BFGS",control = list( trace = TRUE,REPORT=report,abstol=eps,maxit=maxit))
								}else if(grad.num=="hyb"){
									tmp <- optim(par=inixr, fn = obj.sph, gr=Grad.sph, ldiss=IDM[1,,],udiss=IDM[2,,], P=p,method ="BFGS",control = list( trace = TRUE,REPORT=report,abstol=eps,maxit=maxit))
									result.opt <- optim(par=tmp$par, fn = obj.sph, ldiss=IDM[1,,],udiss=IDM[2,,], P=p,method ="BFGS",control = list( trace = TRUE,REPORT=report,abstol=eps,maxit=maxit))
								}
								X <- matrix(result.opt$par[1:(n*p)],nrow=n)
								r <- result.opt$par[(n*p)+1:n]^2
								str <- result.opt$value
								result <- list(X=X,r = r, str = str, str.vec=NA,IDM=IDM, EIDM = idistSph(X,r))
								class(result) <- c("imds","sph")
								return(result)
							}else if(opt.method[1]=="MM"){
								str.vec <- rep(-1,maxit+1)
								result.mmSph <- .C("mmSph",
																arg1 = as.double(iniX),
																arg2 = as.double(inir),
																arg3 = as.double(str.vec),
																arg4 = as.double( array( 0,dim=c(n,n) ) ),#ldist
																arg5 = as.double( array( 0,dim=c(n,n) ) ),#udist
																arg6 = as.double(IDM[1,,]),
																arg7 = as.double(IDM[2,,]),
																arg8 =as.double(eps),
																arg9 = as.integer(rel),
																arg10 = as.integer(dil),
																arg11 = as.integer(n),
																arg12 = as.integer(p),
																arg13 = as.integer(maxit),
																arg14 = as.integer(report)
																) 
								EIDM <- array(0,dim=c(2,n,n))
								EIDM[1,,] <- result.mmSph$arg4
								EIDM[2,,] <- result.mmSph$arg5
								X <- matrix(result.mmSph$arg1,nrow=n,ncol=2)
								tmp.vec <- result.mmSph$arg3
								str.vec <- tmp.vec[which(tmp.vec>-0.5)]
								str <- str.vec[length(str.vec)]
								result <- list(X=X,r= result.mmSph$arg2, str = str, str.vec=str.vec,IDM=IDM,EIDM = EIDM)
								class(result) <- c("imds","sph")
								return(result)
							}
						}else if(model[1]=="box"){
							if(ini[1]=="auto"){
								cmat <- (IDM[2,,] + IDM[1,,])/2
								iniX <- cmdscale(as.dist(cmat),k=p)
								iniR <- matrix(runif(n*p,min = 0.1, max = 1),nrow=n,ncol=p)
								inixr <- c(iniX,sqrt(iniR))
							}else{
								iniX <- ini[[1]]
								iniR <- ini[[2]]
								inixr <- c(ini[[1]],sqrt(ini[[2]]))
							}
							if(opt.method[1]=="BFGS"){
								if(grad.num==FALSE){
									result.opt <- optim(par=inixr, fn = obj.box, gr=Grad.box, ldiss=IDM[1,,],udiss=IDM[2,,], P=p,method ="BFGS",control = list( trace = TRUE,abstol=eps,REPORT=report,maxit=maxit))
								}else if(grad.num==TRUE){
									result.opt <- optim(par=inixr, fn = obj.box, ldiss=IDM[1,,],udiss=IDM[2,,], P=p,method ="BFGS",control = list( trace = TRUE,abstol=eps,REPORT=report,maxit=maxit))
								}else if(grad.num=="hyb"){
									tmp <- optim(par=inixr, fn = obj.box, gr=Grad.box, ldiss=IDM[1,,],udiss=IDM[2,,], P=p,method ="BFGS",control = list( trace = TRUE,abstol=eps,REPORT=report,maxit=maxit))
									result.opt <- optim(par=tmp$par, fn = obj.box, ldiss=IDM[1,,],udiss=IDM[2,,], P=p,method ="BFGS",control = list( trace = TRUE,abstol=eps,REPORT=report,maxit=maxit))
								}
								X <- matrix(result.opt$par[1:(n*p)],nrow=n)
								R <- matrix(result.opt$par[(n*p)+1:(n*p)]^2,nrow=n)
								str <- result.opt$value
								result <- list(X=X,R = R, str = str, str.vec=NA,IDM=IDM, EIDM = idistBox(X,R))
								class(result) <- c("imds","box")
								return(result)
							}else if(opt.method[1]=="MM"){
								str.vec <- rep(-1,maxit+1)
								result.mmBox <- .C("mmBox",
																arg1 = as.double(iniX),
																arg2 = as.double(iniR),
																arg3 = as.double(str.vec),
																arg4 = as.double( array( 0,dim=c(n,n) ) ),#ldist
																arg5 = as.double( array( 0,dim=c(n,n) ) ),#udist
																arg6 = as.double(IDM[1,,]),
																arg7 = as.double(IDM[2,,]),
																arg8 =as.double(eps),
																arg9 = as.integer(rel),
																arg10 = as.integer(dil),
																arg11 = as.integer(n),
																arg12 = as.integer(p),
																arg13 = as.integer(maxit),
																arg14 = as.integer(report)
																) 
								EIDM <- array(0,dim=c(2,n,n))
								EIDM[1,,] <- result.mmBox$arg4
								EIDM[2,,] <- result.mmBox$arg5
								X <- matrix(result.mmBox$arg1,nrow=n,ncol=2)
								R <- matrix(result.mmBox$arg2,nrow=n,ncol=2)
								tmp.vec <- result.mmBox$arg3
								str.vec <- tmp.vec[which(tmp.vec>-0.5)]
								str <- str.vec[length(str.vec)]
								result <- list(X=X,R= R, str = str, str.vec=str.vec,IDM=IDM,EIDM = EIDM)
								class(result) <- c("imds","box")
								return(result)
							}
						}
					}

