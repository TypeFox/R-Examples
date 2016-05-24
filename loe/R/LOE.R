"LOE" <-
function(ADM, p=2, c=0.1,eps= 1e-5,maxit =1000,method=c("BFGS","SD","MM"), iniX = "auto",report=100,DEL=1,H=0.5){
						N <- nrow(ADM)
						diag(ADM) <- -1
						if(iniX[1]=="auto"){
							iniX <- spec.emb(A=ADM,p=p,norm=FALSE)
						}
						iniD <- as.matrix(dist(iniX))
						if(method[1]=="BFGS"){
							inivecx <- as.vector(iniX)
							result.opt <- optim(par=inivecx, fn = Objt.LOE, gr=Grad.LOEobjt, adm=ADM,P=p,C=c,method ="BFGS",control = list( trace = TRUE,REPORT=report,maxit=maxit))
							X <- matrix(result.opt$par,nrow=N)
							str = result.opt$value
							return(list(X=X,str=str))
						}else if(method[1]=="SD"){
							YY <- iniX
							inivecx <- as.vector(iniX)
							result.sdm <- SD.LOE(ini=iniX, adm=ADM,C=c,M=maxit,EPS=eps,del=DEL, h=H,report=report)
							tmpX <- matrix(result.sdm$arg1,N,p)
							tmpstr <- result.sdm$arg2
							return( list(X=tmpX,str=tmpstr[which(tmpstr>-0.5)]) )
						}else if(method[1]=="MM"){
							result.MM <- MM.LOE(ini=iniX, adm=ADM,C=c,M=maxit,EPS=eps,report=report)
							tmpX <- matrix(result.MM$arg1,N,p)
							tmpstr <- result.MM$arg2
							return( list(X=tmpX,str=tmpstr[which(tmpstr>-0.5)]) )
						}
					}

