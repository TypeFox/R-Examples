`mixed.mtc` <-
function (data.rec, data.don, match.vars, y.rec, z.don, method="ML", rho.yz=NULL, micro=FALSE, constr.alg="Hungarian")
{
 #   if(micro && (constr.alg=="Hungarian" || constr.alg=="hungarian")) require(clue)
#    if(micro && (constr.alg=="lpSolve" || constr.alg=="lpsolve")) require(lpSolve)
    
	nA <- nrow(data.rec)
	nB <- nrow(data.don)
	A.lab <- rownames(data.rec)
	if(is.null(A.lab)) A.lab <- 1:nA
	A.lab <- paste("A", A.lab, sep="=")
	B.lab <- rownames(data.don)
	if(is.null(B.lab)) B.lab <- 1:nB
	B.lab <- paste("B", B.lab, sep="=")

#	if(is.list(data.rec))  data.rec <- data.matrix(data.rec)
#	if(is.list(data.don))  data.don <- data.matrix(data.don)
    xrec <- data.rec[,match.vars, drop=FALSE]
    xrec <- fact2dummy(data=xrec, all=FALSE)
    data.rec <- cbind(xrec, data.rec[, y.rec])
    colnames(data.rec) <- c(colnames(xrec), y.rec)

    xdon <- data.don[,match.vars, drop=FALSE]
    xdon <- fact2dummy(data=xdon, all=FALSE)
    data.don <- cbind(xdon, data.don[, z.don])
    colnames(data.don) <- c(colnames(xdon), z.don)

    match.vars <- colnames(xrec)
	p.x <- length(match.vars)
	x.A <- matrix(c(data.rec[,match.vars]), ncol=p.x)
	pos.x.A <- match(match.vars, colnames(data.rec))
	x.B <- matrix(c(data.don[,match.vars]), ncol=p.x)
	pos.x.B <- match(match.vars, colnames(data.don))
	
# preparing dependent variables Y and Z (assumed to be
# continuous)
	y.lab <- y.rec
	y.A <- data.rec[ ,y.lab]
	
	z.lab <- z.don
	z.B <- data.don[ ,z.lab]
	
	p <- p.x+2
	vc <- matrix(NA, nrow=p, ncol=p) 
	v.names <- c(match.vars, y.lab, z.lab)
	dimnames(vc) <- list(v.names, v.names)
	pos.x <- 1:p.x
	pos.y <- p.x+1
	pos.z <- p.x+2
	if(method=="ML"){

		# regression in file B: Z vs. X
		lm.B <- lm(z.B ~ x.B)
		res.B <- residuals(lm.B)
		se.B <- sqrt(mean(res.B^2))  # ML estimate of V(Z|X)  
		  
		# regression in file A: Y vs. X
		lm.A <- lm(y.A ~ x.A)
		res.A <- residuals(lm.A)
		se.A <- sqrt(mean(res.A^2))  # ML estimate of V(Y|X)  

		# ML estimates for X variables
		mu.x <- colMeans(rbind(x.A, x.B))
		S.x <- var(rbind(x.A, x.B))*((nA+nB-1)/(nA+nB))

		# ML estimates for Y
		mu.y <- sum(coefficients(lm.A)*c(1,mu.x))
		S.yx <- coefficients(lm.A)[-1] %*% S.x
		S.y <- se.A^2 + S.yx %*% solve(S.x) %*%  t(S.yx)

		# ML estimates for Z
		mu.z <- sum(coefficients(lm.B)*c(1,mu.x))
		S.zx <- coefficients(lm.B)[-1] %*% S.x
		S.z <- se.B^2 + S.zx %*% solve(S.x) %*%  t(S.zx)
      
		# ML estimates for Y,Z, given the input value rho.yz for partial.cor[(Y,Z)|X]      
		if(is.null(rho.yz)) rho.yz <- 0 # CI assumption
        S.yzGx <- rho.yz * (se.A * se.B) # partial Cov[(Y,Z)|X]
		S.yz <- S.yzGx + S.yx %*% solve(S.x) %*% t(S.zx)
		S.zGyx <- se.B^2 - S.yzGx %*% solve(se.A^2) %*% t(S.yzGx)
		S.yGzx <- se.A^2 - S.yzGx %*% solve(se.B^2) %*% t(S.yzGx)
	
		vc[pos.x, pos.x] <- S.x
		vc[pos.x, pos.y] <- vc[pos.y, pos.x] <- S.yx
		vc[pos.x, pos.z] <- vc[pos.z, pos.x] <- S.zx
		vc[pos.y, pos.y] <- S.y
		vc[pos.z, pos.z] <- S.z
		vc[pos.y, pos.z] <- vc[pos.z, pos.y] <- S.yz
		
		# prepare for output
		mu <- c(mu.x, mu.y, mu.z)
		names(mu) <- v.names
		res.var <- c(S.yGzx=c(S.yGzx), S.zGyx=c(S.zGyx))
		fine <- list(start.prho.yz=rho.yz, mu=mu, vc=vc, cor=cov2cor(vc), res.var=res.var)
	}
	if(method=="MS"){
		# estimates for X variables
		mu.x <- colMeans(rbind(x.B,x.A) )
		S.x <- var(rbind(x.B, x.A))
		
		# estimates for Y
		mu.y <- mean(y.A)
		S.y <- var(y.A)
		S.xy <- var(x.A, y.A)

		# estimates for Z	
		mu.z <- mean(z.B)
		S.z <- var(z.B)
		S.xz <- var(x.B,z.B)
		
		#fills in the Var-Cov matrix
		vc[pos.x, pos.x] <- S.x
		vc[pos.x, pos.y] <- vc[pos.y, pos.x] <- S.xy
		vc[pos.x, pos.z] <- vc[pos.z, pos.x] <- S.xz
		vc[pos.y, pos.y] <- S.y
		vc[pos.z, pos.z] <- S.z

        # estimation of S.yz
        # step.1 checks if the input value for Cor(Y,Z), rho.yz, is admissible  
		if(p.x==1){
       		c.xy <- c(cor(x.A, y.A))
    		c.xz <- c(cor(x.B, z.B))
			low.c <- c.xy*c.xz - sqrt( (1-c.xy^2)*(1-c.xz^2) )
			up.c <-  c.xy*c.xz + sqrt( (1-c.xy^2)*(1-c.xz^2) )
            rho.yz.CI <- c.xy*c.xz
		}      
		else{
            eps <- 0.0001
            cc <- cov2cor(vc)
            rr <- seq(-1, 1, eps)
            k <- length(rr)
            vdet <- rep(0,k)
            for(i in 1:k){
                cc[pos.z, pos.y] <- cc[pos.y, pos.z] <- rr[i]
                vdet[i] <- det(cc)
            }
            cc.yz <- rr[vdet>=0]
            low.c <- min(cc.yz)
            up.c <- max(cc.yz)
            rho.yz.CI <- (1/2)*(low.c + up.c)
		}
        if(is.null(rho.yz)) rho.yz <- rho.yz.CI
        # step.2 checks whether the input value rho.yz for Cor(Y,Z) is addmisible. Otherwise takes the closest admissible value.
		cat("input value for rho.yz is", rho.yz, fill=TRUE)
		cat("low(rho.yz)=", low.c, fill=TRUE)
		cat("up(rho.yz)=", up.c, fill=TRUE)
 		sum.rho.yz <- c(start=rho.yz)
		if( (rho.yz<=up.c) && (rho.yz>=low.c) ) {
			cat("The input value for rho.yz is admissible", fill=TRUE)
		}      
		else{
			cat("Warning: value for rho.yz is not admissible: a new value is chosen for it", fill=TRUE)
			if(rho.yz > up.c) rho.yz <- up.c - 0.01
			if(rho.yz < low.c) rho.yz <- low.c + 0.01
			cat("The new value for rho.yz is ", rho.yz, fill=TRUE)
		}
		S.yz <- rho.yz * sqrt(S.y * S.z)
		sum.rho.yz <- c(sum.rho.yz, low.lim=low.c, up.lim=up.c, used=rho.yz)
		vc[pos.y, pos.z] <- vc[pos.z, pos.y] <- S.yz
		
		fi.3 <- rbind(c(S.xz, S.yz)) %*% solve(vc[c(pos.x,pos.y),c(pos.x,pos.y)]) %*% cbind(c(S.xz, S.yz))
		S.zGyx <- S.z - fi.3
		if(S.zGyx<0) S.zGyx <- 0
		fi.6 <- rbind(c(S.xy, S.yz)) %*% solve(vc[c(pos.x,pos.z),c(pos.x,pos.z)]) %*% cbind(c(S.xy, S.yz))
		S.yGzx <- S.y - fi.6
		if(S.yGzx<0) S.yGzx <- 0
	
	# prepare for output
		mu <- c(mu.x, mu.y, mu.z)
		names(mu) <- v.names
		res.var <- c(S.yGzx=c(S.yGzx), S.zGyx=c(S.zGyx))
		fi <- c(fi.6.y=fi.6, fi.3.z=fi.3)
		fine <- list(rho.yz=sum.rho.yz, mu=mu, vc=vc, cor=cov2cor(vc), phi=fi, res.var=res.var) 
	}
	
	if(micro){
		if(nA>nB) stop("The number of donors is less than the number of recipients")
		if(method=="ML"){
			# Prediction of the values of Z in file A
			z.pred <- mu.z + rbind(vc[pos.z,c(pos.x, pos.y)]) %*%  solve(vc[c(pos.x,pos.y),c(pos.x,pos.y)]) %*%  
																			t(cbind(sweep(x.A, 2, mu.x),y.A-mu.y))
			# Prediction of the values of Y in file B
			y.pred <- mu.y + rbind(vc[pos.y,c(pos.x,pos.z)]) %*% solve(vc[c(pos.x,pos.z),c(pos.x,pos.z)]) %*% 
																		t(cbind(sweep(x.B, 2, mu.x), z.B-mu.z))
			rSS <- vc[c(pos.y,pos.z), c(pos.y,pos.z)]
		}
		else if(method=="MS"){	
			# Prediction of the values of Z in file A
			B.zGxy <- t(rbind(c(S.xz, S.yz)) %*% solve(vc[c(pos.x,pos.y),c(pos.x,pos.y)]))
			z.pred <- mu.z + cbind( t(t(x.A)- mu.x), y.A-mu.y   ) %*% B.zGxy
		
			# Prediction of the values of Y in file B
			B.yGxz <- t(rbind(c(S.xy, S.yz)) %*% solve(vc[c(pos.x,pos.z),c(pos.x,pos.z)]))
			y.pred <- mu.y + cbind( t(t(x.B)-mu.x), z.B-mu.z   ) %*% B.yGxz
			
			S1 <- S2 <- vc
			S1[pos.z, pos.z] <- fi.3
			S2[pos.y, pos.y] <- fi.6
			SS <- S1+S2
			rSS	<- SS[c(pos.y,pos.z), c(pos.y,pos.z)]
		}
		irSS <- solve(rSS)
		z.ep <- c(z.pred) + rnorm(nA, 0, sqrt(S.zGyx)) 
		y.ep <- c(y.pred) + rnorm(nB, 0, sqrt(S.yGzx))
		new.B <- cbind(y.ep, z.B)
		madist <- matrix(0, nA, nB)
		for(i in 1:nA){
			new.A <- c(y.A[i], z.ep[i])
			madist[i,] <- sqrt(mahalanobis(new.B, new.A, irSS, inverted=TRUE))
		}
		dimnames(madist) <- list(A.lab, B.lab)
		# constrained nearest neighbour matching matching is performed
		if(constr.alg=="lpSolve" || constr.alg=="lpsolve"){
		    
		    if(nA==nB) appo <- lp.assign(cost.mat=madist)
		    else if(nA<nB){
		        r.sig <- rep("==", nA)
		        r.rhs <- rep(1, nA)
		        c.sig <- rep("<=", nB)
		        c.rhs <- rep(1, nB)
		        appo <- lp.transport(cost.mat=madist, row.signs=r.sig, row.rhs=r.rhs, col.signs=c.sig, col.rhs=c.rhs)
		    }   
		    sol <- appo$solution
		    ss <- c(t(sol))
		    cc <- c(t(col(sol)))
		    dist.rd <- madist[cbind(1:nA, cc[as.logical(ss)] )]
		    don.lab <- B.lab[c(cc[as.logical(ss)])]
		}
		
		# the function solve_LSAP in package clue is used
		else if(constr.alg=="hungarian" || constr.alg=="Hungarian"){
		    if(nA > nB) stop("It is required  that the no. of donors \n 
		                      is equal or greater than the no. of recipients")
		    sol <- solve_LSAP(x=madist, maximum=FALSE)
		    don.lab <- B.lab[as.integer(sol)]
		    dist.rd <- madist[cbind(A.lab, don.lab)]
        }
		
		rec.lab <- substring(A.lab, 3)
		don.lab <- substring(don.lab, 3)
		if(is.null(rownames(data.don))) {
			rec.lab <- as.numeric(rec.lab)
			don.lab <- as.numeric(don.lab)
		}
		mtc.ids <- cbind(rec.id=rec.lab, don.id=don.lab)
		fill.A <- cbind(data.rec, data.don[don.lab, z.lab])
		colnames(fill.A) <- c(colnames(data.rec), z.lab)
		
		#output
		fine$filled.rec <- fill.A
		fine$mtc.ids <- mtc.ids
		fine$dist.rd <- dist.rd
		
	}
	fine$call <- match.call()
  fine
}

