"SOE" <-
function(CM, N,p=2, c=0.1,maxit =1000,report=100, iniX = "rand",rnd=10000){
						if(iniX[1]=="rand"){
							iniX <- mvrnorm(N,rep(0,p),diag(p))
						}
						if(nrow(CM)>100000){
							sid <- sample(1:nrow(CM))
							PCM <- CM[sid[1:rnd],]
							inivecx <- as.vector(iniX)
							result.opt <- optim(par=inivecx, fn = Objt.SOE, gr= Grad.SOEobjt,cm=PCM,n=N,P=p,C=c,method ="BFGS",control = list( trace = TRUE,REPORT=report,maxit=maxit))
							X <- matrix(result.opt$par,nrow=N)
							str = result.opt$value
							return(list(X=X,str=str))
						}else{
							inivecx <- as.vector(iniX)
							result.opt <- optim(par=inivecx, fn = Objt.SOE, gr= Grad.SOEobjt,cm=CM,n=N,P=p,C=c,method ="BFGS",control = list( trace = TRUE,REPORT=report,maxit=maxit))
							X <- matrix(result.opt$par,nrow=N)
							str = result.opt$value
							return(list(X=X,str=str))
						}
					}

