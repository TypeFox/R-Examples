


marqLevAlg <- function(b,m=FALSE,fn,gr=NULL,hess=NULL,maxiter=500,epsa=0.001,epsb=0.001,epsd=0.01,digits=8,print.info=FALSE,blinding=TRUE,multipleTry=25){
	cl <- match.call()
	if (missing(m) & missing(b)) stop("The 'marqLevAlg' alogorithm needs a vector of parameters 'b' or his length 'm'")
	if(missing(m)) m <- length(b)	
	if(missing(b)) b <- rep(0.1,m)
	if(length(b) != m){
		if(length(b) < m){
			b.temp <-NULL
			b.temp <- c(b.temp,b,rep(b[1],(m-length(b))))
			b <- b.temp 
		}else{
			m <- length(b)
		}
	}  

	if(missing(fn)) stop("The argument 'funcpa' is missing.")

	funcpa <- function(b){-fn(b)}
	if(!is.null(gr)) grad <- function(b){-gr(b)}
	if(!is.null(hess)) hessian <- function(b){hess(b)}

	flush.console()
	ptm <- proc.time()
	cat("\n")
	cat("Be patient. The program is computing ...\n")
	
###initialisation
	binit <- b
	th <- 1e-5
	eps <- 1e-7
	nfmax <- m*(m+1)/2    
	ca <- epsa+1
	cb <- epsb+1
	rl1 <- -1.e+10    
	ni <- 0
	istop <- 0
	da <- 1E-2
	dm <- as.double(5)
	nql <- 1
	m1 <- m*(m+1)/2
	ep <- 1E-20
	delta <- rep(0,m)
	b1 <- rep(0,m)
	idpos <- 0
	v <- rep(0,m*(m+3)/2)
	ind.func1 <- 0
	fu <- rep(0,(m*(m+3)/2))
	gonflencountmax <- 10

## old parameters iteration -1
	old.b <- b
	old.rl <- 0
	old.ca <- 1
	old.cb <- 1
	old.dd <- 1	
## 	
	repeat{	
	
		if (sum(!is.finite(b))>0){

			cat("Probably too much accuracy requested...\n")
			cat("Last step values :\n")
			cat("      b :",round(old.b,digits),"\n")
			cat("      likelihood :",round(-old.rl,digits),"\n")
			cat("      Convergence criteria: parameters stability=", round(old.ca,digits), "\n")
			cat("                          : likelihood stability=", round(old.cb,digits), "\n") 
			cat("                          : best relative distance to maximum obtained (RDM)=", round(old.dd,digits), "\n")
			stop("")
			 
		}
		res.out.error <- list("old.b"=round(old.b,digits),"old.rl"=round(old.rl,digits),"old.ca"=round(old.ca,digits),"old.cb"=round(old.cb,digits),"old.dd"=round(old.dd,digits))
	
		if(missing(gr)){
			deriv <- deriva(b,funcpa)
			v <- deriv$v
			rl <- deriv$rl
			
			if((multipleTry > 1) & (ni ==0)){
				kk <- 0
				while(((kk < multipleTry) & (!is.finite(rl)))){
					kk <- kk + 1
					b <- b/2
					deriv <- deriva(b,funcpa)
					v <- deriv$v
					rl <- deriv$rl
				}
			} 
		}else{
			v <- NULL
			rl=funcpa(b)
			
			if(missing(hess)){
				deriv <- deriva_grad(b,grad)
				v <- c(v,deriv$hessian,grad(b))
			}else{
				tmp.hessian <- hessian(b) 
				if(is.matrix(tmp.hessian)){
					tmp.hessian <- tmp.hessian[upper.tri(tmp.hessian,diag=T)]
				}
				v <- c(v,tmp.hessian,grad(b))	
			}
			
		}
		if((sum(is.finite(b))==m) && !is.finite(rl)){
			cat("Problem of computation. Verify your function specification...\n")
			cat("Infinite likelihood with finite parameters : b=",round(old.b,digits),"\n")
			cat("      - Check the computation and the continuity,\n")
			cat("      - Check that you minimize the function.\n")
			stop("")
		
		}

		if(((sum(!is.finite(b)) > 0) || (sum(!is.finite(rl)) > 0)) && (ni==0)){
			cat("Problem of computation. Verify your function specification...\n")
			cat("First check b length in parameters specification.\n")
			stop("")
		}
		rl1 <- rl      
		dd <- 0 
			
		for(i in 1:m){
			for(j in i:m){
				ij <- (j-1)*j/2+i
				fu[ij]=v[ij]
			}
		}
		fu.int <- fu[1:(m*(m+1)/2)]
	
		dsinv <- .Fortran("dsinv",fu.out=as.double(fu.int),as.integer(m),as.double(ep),ier=as.integer(0),det=as.double(0),PACKAGE="marqLevAlg")
				
		ier <- dsinv$ier
		fu[1:(m*(m+1)/2)] <- dsinv$fu.out
		if (ier == -1){
			#dd <- epsd+1 #without dd approximation by Pierre Joly
			v_tmp <- v[(m*(m+1)/2+1):(m*(m+3)/2)]
			dd <- sum(v_tmp*v_tmp)
		}else{
			dd <- ghg(m,v,fu)$ghg/m
		}
		
        if(print.info){
		cat("------------------ iteration ",ni,"------------------\n")
		cat("Log_likelihood ",round(-rl,digits),"\n")
		cat("Convergence criteria: parameters stability=", round(ca,digits), "\n")
		cat("                    : likelihood stability=", round(cb,digits), "\n") 
		if (ier == -1){
			cat("                    : Matrix inversion for RDM failed \n")	
		}else{
			cat("                    : Matrix inversion for RDM successful \n")
		}
		cat("                    : relative distance to maximum(RDM)=", round(dd,digits), "\n")

		nom.par <- paste("parameter",c(1:m),sep="")
		id <- 1:m
		indice <- rep(id*(id+1)/2)
		Var <- fu[indice]
		SE <- sqrt(abs(Var))
		res.info <- data.frame("coef"=round(b,digits),"SE.coef"=round(SE,digits),"Var.coef"=round(Var,digits))
		rownames(res.info) <- nom.par
		cat("\n")
	}
		old.b <- b
		old.rl <- rl
		old.ca <- ca
		old.cb <- cb
		if(dd<=old.dd){old.dd <- dd}
		if((ca < epsa) & (cb < epsb) & (dd < epsd)){break}

		
		tr <- 0
		for(i in 1:m){
			ii <- i*(i+1)/2
			tr <- tr+abs(v[ii])
		}
		tr <- tr/m
		
		ncount <- 0
		ga <- 0.01

		fu <- v
		
		for(i in 1:m){
			ii <- i*(i+1)/2
			if (v[ii] != 0){
				fu[ii] <- v[ii]+da*((1.e0-ga)*abs(v[ii])+ga*tr)
			}else{
				fu[ii] <- da*ga*tr
			}
		}
			
		dchole <- .Fortran("dchole",fu=as.double(fu),as.integer(m),as.integer(nql),idpos=as.integer(0),PACKAGE="marqLevAlg")
		fu <- dchole$fu
		idpos <- dchole$idpos

		while(idpos != 0){
		
			ncount <- ncount + 1
			if((ncount <= 3) | (ga >= 1)){
				da <- da * dm
			}else{
				ga <- ga * dm
				if(ga > 1) ga <- 1
			}
			
			fu <- v
			
			for(i in 1:m){
				ii <- i*(i+1)/2
				if (v[ii] != 0){
					fu[ii] <- v[ii]+da*((1.e0-ga)*abs(v[ii])+ga*tr)
				}else{
					fu[ii] <- da*ga*tr
				}
			}
			
		dchole <- .Fortran("dchole",fu=as.double(fu),as.integer(m),as.integer(nql),idpos=as.integer(0),PACKAGE="marqLevAlg")
		
			idpos <- dchole$idpos
			fu <- dchole$fu
			if (ncount >= gonflencountmax){
				break
			} 

		}
		delta <- fu[(nfmax+1):(nfmax+m)]
		b1 <- b + delta
		rl <- funcpa(b1)
		
		if(blinding){
			if(is.na(rl)){
				cat("rl :",rl,"\n")
				rl <- -500000
			}
		}else{
			if(is.na(rl)){
				cat(" Probably wrong definition of function FN \n")
				cat("      ->  invalid number (infinite or NA)\n")
				cat("          value of function is :",round(-rl,digits),"\n")
				
				istop <- 4
				break
			}
		}
		if (rl1 < rl){
		
			if(da < eps){
				da <- eps
			}else{
				da <- da/(dm+2)
			}
			goto800 <- func1(b,rl1,rl,delta,ni,maxiter)
			ca <- goto800$ca
			cb <- goto800$cb
			b <- goto800$b
			ni <- goto800$ni
			ind.func1 <- 1
			if(ni >= maxiter){
				istop <- 2
				break
			}
		}else{
			maxt <- max(abs(delta)) 
	
			if(maxt == 0){
				vw <- th
			}else{
				vw <- th/maxt
			}
		
			step <- log(1.5)

			sears <- searpas(vw,step,b,delta,funcpa,res.out.error)

			fi <- sears$fi
			vw <- sears$vw
			rl <- -fi
			if(rl == -1.e9){
				istop <- 4
				break
			}
			delta <- vw*delta
			da <- (dm-3)*da
			goto800 <- func1(b,rl1,rl,delta,ni,maxiter)
 			ca <- goto800$ca
			cb <- goto800$cb
			b <- goto800$b
			ni <- goto800$ni
			   
			if(ni >= maxiter){
				istop <- 2
				break
			}
		}

	}

	if((istop %in% 2:4)==F) istop <- 1
	cost <- proc.time() - ptm
	result <- list(cl=cl,ni=ni,ier=ier,istop=istop,v=fu[1:(m*(m+1)/2)],fn.value=-rl,b=b,ca=ca,cb=cb,rdm=dd,time=round(cost[3],3))
	class(result) <- "marqLevAlg"
	

	cat("The program took",round(cost[3],3), "seconds \n")
	
	result
}























