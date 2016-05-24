

.DrHSVM.predict <- function(object, newx, newy=NULL, newlam=NULL, eps = 1e-10, verbose=TRUE) {
# #merge the official version of this function with predict.lars 
# --> predict coeffs for new lambda point.
# if newlam exists --> Extract coefficients from a fitted DrHSVM model 
# Otherwise Predict fitted values for new samples. 
# 					In case of having labels for new samples, misclassification errors are also calculated  


	# interpolate  coeffs for a new point and returns the new coefs
	if (!is.null(newlam)){ 
		if (newlam < 0 ) stop("Error in .DrHSVM.predict: new.lam1 should be non negative. ")
		lambdas = object$ lambda1
		
		# check extreme cases:
		# 1. if the new lambda newlam is larger than the bigest lambda in the path ---> set it equal to bigest lambda
		# because the biggest lambda has all coef beta=0! 
		newlam[newlam > max(lambdas)] = max(lambdas)
		# 2. if the new lambda newlam is smaller than the smallest lambda in the path
		# -->  set it equal to smallest lambda
		if (verbose) print(paste("By prediction of Elastic Net coefs change the value of the new lambda from ", newlam, "to ", min(lambdas),
																", which is the smallest lambda value in the path." ))
		newlam[newlam < min(lambdas)] = min(lambdas)		
		
		if (newlam %in% object$lambda1){
			# new lam is one of the knicks
			newbeta<- object$beta[newlam == object$lambda1, ]
			newbeta0<- object$beta0[newlam == object$lambda1 ]
		
		} else{ # estimate coefs for a new point 
			
			# put the new point at the right place in the path
			# find neighbours for the new point
			tmp<-c(newlam, lambdas)
			rank(tmp)[1] # rank of the new lambda
			lam.l<-tmp[rank(tmp) == (rank(tmp)[1] - 1) ]
			lam.r<-tmp[rank(tmp) == (rank(tmp)[1] + 1) ]
			
			betas.l<-c(object$beta0[which(lambdas  == lam.l) ], 
								object$beta[which(lambdas  == lam.l),  ])
			betas.r<-c(object$beta0[which(lambdas  == lam.r) ], 
								object$beta[which(lambdas  == lam.r),  ])						
			names(betas.l)<-names(betas.r)<- c("b", colnames(newx))
			
			tmp.mat<-as.matrix(data.frame(lam.l, lam.r, betas.l, betas.r, newlam ))
			
			# use linear interpolation for each feature
			.lin.inter<-function(points ){
				# http://en.wikipedia.org/wiki/Linear_interpolation
				# decode
				x.l<- points[1]; x.r<- points[2]
				y.l<- points[3]; y.r<- points[4]
				x.new<- points[5]
				if ( x.new < x.l | x.new > x.r | x.l > x.r  ) stop ("linear interpolation: the new point should lie between 2 known points")
				
				return(y.new<- y.l + (x.new-x.l)*(y.r-y.l)/(x.r-x.l))
			} # end of  use linear interpolation for each feature
			
			betas.new<- apply(	tmp.mat, 1, .lin.inter  ) 
			
			#data.frame(lam.l, lam.r, betas.l, betas.r, newlam, betas.new=betas.new ) [1:20,]
			newbeta0<- betas.new[1]
			newbeta<- betas.new[-1]
		} # end of interpolate coefs for a new point
		
		return(list(newbeta0=newbeta0,newbeta=newbeta))
	
	} else {# if (is.null(newlam)){	
		### ELSE fit the path and return misclassification errors
		### Get rid of many zero coefficients
		# refine the path by adding the new point 'newlam', put the lambda at the right place into the path
		
		coef1 <- object$beta 
		c1 <- drop(rep(1, nrow(coef1)) %*% abs(coef1))   
		nonzeros <- c1 > eps
		coef1 <- cbind(object$beta0, coef1[, nonzeros])
		newx <- cbind(1, newx[, nonzeros, drop=FALSE])
		fit <- newx %*% t(coef1)
		predict <- sign(fit)
		if (is.null(newy))
			return(list(coef = coef1, fit = fit, err=NULL))
		
		if (length(newy) ==1 ){
			err<- as.numeric(drop(predict) != newy) 
		}else err <- apply(apply(predict, 2, FUN="!=", newy), 2, sum)/length(newy)
		
		return(list(fit = fit, err = err))
	}
}
