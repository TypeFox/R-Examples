RAMP <-
function(X, y, family = "gaussian", penalty = "LASSO", gamma, inter = TRUE, eps = 1e-15, tune='EBIC', lam.list, lambda.min.ratio, max.iter = 100, max.num, n.lambda = 100, ebic.gamma = 1, refit=TRUE, trace = FALSE) {
	if (penalty == "SCAD" & missing(gamma)) {
		gamma = 2.3
	}
	if (penalty == "MCP" & missing(gamma)) {
		gamma = 1.3
	}
	if (penalty == "LASSO" & !missing(gamma)) {
		stop("gamma is not needed")
	}

	n = dim(X)[1]
	p = dim(X)[2]
	oldX = X
	if (missing(max.num)) 
		max.num = p + 1
	else max.num = min(max.num, p + 1)

	if (missing(lambda.min.ratio)) 
		lambda.min.ratio = ifelse(n < p, 0.01, 1e-04)

	################prepare variables
	lambda = NULL
	beta.mat = NULL ##storing the beta coefficents along the path
	beta = rep(0, p) ##the current beta estimates including the interation terms
    
	index = NULL


	################standardize design matrix
	X = scale(X)
	xm0 = attr(X,"scaled:center")
	xsd0 = attr(X,"scaled:scale")
    

	#############lambda list
	if (family == "binomial") {
		max.lam = max(abs((1/n) * (t(X) %*% (y - mean(y)))))
		a0list = rep(log(mean(y)/(1 - mean(y))), n.lambda)
	}
	if (family == "poisson") {
		max.lam = max(abs((1/n) * (t(X) %*% (y - mean(y)))))
		a0list = rep(log(mean(y)), n.lambda)

	}

	if (family == "gaussian") {
		max.lam = max(abs((1/n) * (t(X) %*% (y - mean(y)))))
		a0list = rep(mean(y), n.lambda)

	}

	if (missing(lam.list)) {
		min.lam = max.lam * lambda.min.ratio
		lam.list = exp(seq(from = log(max.lam), to = log(min.lam), length.out = n.lambda))
	} else {
		lam.list = lam.list[lam.list <= max.lam]
		n.lambda = length(lam.list)
	}




    a0 = a0list[1]
	loglik.list = rep(0, n.lambda)
	cri.list = rep(0,n.lambda)
	AIC.list = rep(0,n.lambda)
	BIC.list = rep(0,n.lambda)
	EBIC.list = rep(0,n.lambda)
	df.list = rep(0,n.lambda)
	df.m.list = rep(0,n.lambda)
	df.i.list = rep(0,n.lambda)
	ind.list.inter = as.list(NULL)
	ind.list.main = as.list(NULL)
	ind.list.inter.xlab = as.list(NULL)
	beta.list.inter = as.list(NULL)
	beta.list.inter[[1]] = NULL
	ind.list.inter[[1]] = NULL
	ind.list.main[[1]] = NULL
	ind.list.inter.xlab[[1]] = NULL
	
	##############main part
	

	##expand X matrix dynamically each step by incorporating the new candidates of interactions,
	##which are the interactions of the current active variables minus the current interaction 
#effects, the interactions will play the same role as the main effect, except we keep track 
#of the list through the lambda sequence. Keep track of other indices as well.colnames(X)

	colnames(X) = paste("X", 1:p, sep = "")

	nonPen = rep(0, p)
	for (k in 1:n.lambda) {
		ind.list.inter.xlab[[k]] = 'None'
		#if(k==10)browser()
		#print(k)
		old.beta = beta ##back up the old beta
		nonPen = rep(0, p)
		
		break.ind = FALSE
        

		lam = lam.list[k]


		index = which(abs(beta[1:p]) > eps) ###selected main effects from last step
		if (length(index) > 0) { ###find the candidate interation effects following strong heredity
			aa = outer(index, index, f <- function(x, y) {
				paste("X", x, "X", y, sep = "")
			})
			aa[lower.tri(aa,diag=FALSE)] = NA
			bb = as.vector(aa)
			bb = bb[!is.na(bb)]
			curInter = colnames(X)[-(1:p)]
			candInter = setdiff(bb, curInter)	
			ncurInter = length(curInter)
			ncandInter = length(candInter)
						
		
			#cat("k=",k,'candidate interaction', candInter, '\n', sep=' ')
			if (ncurInter > 0) { ##active interaction terms, setting the penalty for the parents to 0
				for (indInter in 1:ncurInter) {
					pair = as.numeric(strsplit(curInter[indInter], "X")[[1]][2:3])
			nonPen[pair[1]] = 1
			nonPen[pair[2]] = 1
				}
				}
			#nonPen[index] = 1
				Xinter = NULL
			if (ncandInter > 0 && inter) { ##candidate interaction terms
			

				for (indInter in 1:ncandInter) {
					pair = as.numeric(strsplit(candInter[indInter], "X")[[1]][2:3])
					Xinter = cbind(Xinter, X[, pair[1]] * X[, pair[2]]) 
					
 
				}
				Xnew = cbind(X, Xinter)
				# allxm = c(xm,xinterm)
				# allxvar = c(xvar,xintervar)

				colnames(Xnew) = c(colnames(X), candInter)
				X = Xnew
				
				beta = c(beta, rep(0, ncandInter)) #expand the beta coefficent vector to include the candiate interaction terms.
			}
		}
		X[,1:p] = oldX
	    X = scale(X)
	    xm = attr(X,"scaled:center")
	    xvar = attr(X,"scaled:scale")
	   # xm[1:p] = m.m
	   # xvar[1:p] = m.var
		for (ite in 1:max.iter) {
			if (penalty == "LASSO") {
				if (ite == 1) {
					cd.temp1 = cd.lasso(X = X, y = y, a0 = a0, beta = beta, epsilon = eps, max.iter = 1, lambda = lam, 
						family = family, bInd = break.ind, nonPen = nonPen)
					a0 = cd.temp1$a0
					beta = cd.temp1$beta
					ind1 = which(abs(beta)>eps)
				}


				# #CD
				if (ite > 1) {
					cd.temp2 = cd.lasso(X = X[, ind1], y = y, a0 = a0, beta = beta[ind1], epsilon = eps, max.iter = max.iter, lambda = lam, family = family, bInd = break.ind, nonPen = nonPen[ind1[ind1 <= p]])
					a0list[k] = cd.temp2$a0
					beta2 = beta
					beta2[ind1] = cd.temp2$beta
				

					# ########redetect active set
					# check.beta=new.beta
					cd.temp3 = cd.lasso(X = X, y = y, a0 = a0, beta = beta2, epsilon = eps, max.iter = 1, lambda = 					lam,family = family, bInd = break.ind, nonPen = nonPen)
					a0list[k] = cd.temp3$a0
					beta = cd.temp3$beta
					ind3 = which(abs(beta)>eps)
					if (setequal(ind1,ind3)) {
						break
					}
					ind1 = ind3

				} ##END iter>1 

			} #END LASSO 


			#end check-cd-recheck
			} ##END iter
			
		#cat('k=',k,'iter num',ite,'\n',sep=' ')
			
		#}

		#end case if sparse enough
		#}

		ind.list.main[[k]] = which(abs(beta[1:p])> eps)
		ind.list.inter[[k]] = which(abs(beta[-(1:p)])> eps ) #record the interaction index pair location list
	    beta.ols = beta
        size.main = length(ind.list.main[[k]])
        size.inter = length(ind.list.inter[[k]])
        
        index = which(abs(beta)>eps)
        beta[-index] = 0
        if(size.inter>0){ ###if interaction effects are detected,
        #enforce the strong heredity in this step
        tmpindname = colnames(X)[ind.list.inter[[k]]+p]
    	cur.main.ind = intertomain(tmpindname)
    	if(trace==T)
    	cat('k=', k ,'Enforced main effects', setdiff(cur.main.ind,ind.list.main[[k]]), '\n')
      
    	ind.list.main[[k]] = union(ind.list.main[[k]],cur.main.ind)
    	index = union(index,cur.main.ind)
    	size.main = length(ind.list.main[[k]])
        }
        index = sort(index)
        
        df.list[k] = size.main+size.inter
        
        loglik.list[k] = 10^15 
        refit.beta = beta
        refit.a0 = a0
	    if(df.list[k]>0 && refit==TRUE &&length(index)<n/2 ){ ##update the estimate using the result from an OLS refit.
	        ols.fit = ols.refit(y, oldX, index[index<=p], colnames(X)[ind.list.inter[[k]]+p], family=family)
	        
	        refit.beta[index] = ols.fit$beta[-1]
	        refit.a0 = ols.fit$beta[1]
	        loglik.list[k] = ols.fit$loglik
	        
	        if(family=='binomial'&&min(abs(ols.fit$link))>10) break
	
	        }
	    
	  		if (length(beta) > p) {
			tmp = which(abs(beta[-(1:p)]) > eps)
			beta = beta[c(1:p,p+tmp)] ##update beta and X with only the main effect + important interaction effect
			refit.beta = refit.beta[c(1:p,p+tmp)]
	
			X = X[,c(1:p,p+tmp)]
			#ind.list.inter.xlab[[k]] = NULL
			if (length(tmp) > 0) { ##if interaction effects are selected. 
			    ind.list.inter.xlab[[k]] = colnames(X)[-(1:p)]
				beta.list.inter[[k]] = refit.beta[-(1:p)]
				}
		else {ind.list.inter.xlab[[k]] = 'None'}
		}
		
		a0list[k] = refit.a0
		beta.mat = cbind(beta.mat, refit.beta[1:p])
		#break criteria
		
	   

			AIC.list[k] = loglik.list[k] + 2*length(index)
			BIC.list[k] = loglik.list[k] + log(n)*length(index)
			EBIC.list[k]= loglik.list[k] + log(n)*length(index)+ 2*ebic.gamma*log(choose(p+size.main^2,size.main+size.inter))
			if(tune=='AIC') cri.list[k] = AIC.list[k]														
			if(tune=='BIC') cri.list[k] = BIC.list[k]						
			if(tune=='EBIC') cri.list[k] = EBIC.list[k]	
			
			df.m.list[k] = length(ind.list.main[[k]])
		    df.i.list[k] = length(beta)-p

			if (sum(beta != 0) >= min(n-1,40)) 
			break
			if (break.ind == TRUE) {
			print("Warning: Algorithm failed to converge for all values of lambda")
			break
		}

		
		if(trace==T){
        cat("k=",k,'current main effects', ind.list.main[[k]], '\n', sep=' ')
		cat("k=",k,'current interaction', ind.list.inter.xlab[[k]], '\n', sep=' ')
		}
		#end the outer loop: decrease in lambda	
				}
     cri.loc = which.min(cri.list[1:(k-1)])
     AIC.loc = which.min(AIC.list[1:(k-1)])
     BIC.loc = which.min(BIC.list[1:(k-1)])
     EBIC.loc = which.min(EBIC.list[1:(k-1)])

	lambda = lam.list[1:ncol(beta.mat)]
	
		val = list(a0 = a0list, beta.m.mat = beta.mat, beta.i.mat = beta.list.inter, beta.m = beta.mat[ind.list.main[[cri.loc]], cri.loc], beta.i= beta.list.inter[[cri.loc]], df=df.list[1:k],df.m=df.m.list[1:k],df.i=df.i.list[1:k], lambda = lambda[1:k], mainInd.list = ind.list.main, mainInd = ind.list.main[[cri.loc]], cri.list=cri.list, cri.loc = cri.loc,
interInd.list = ind.list.inter.xlab,interInd=ind.list.inter.xlab[[cri.loc]], family = family, X = oldX, y = y)
	class(val) = "RAMP"
	return(val)
}
