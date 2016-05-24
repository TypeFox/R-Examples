skktest <-
function(X,Y = NULL){
# X is a n1 x p matrix with columns corresponding to subjects.
# Y is a n2 x p matrix with columns corresponding to subjects.
Result = rep(0,2);
if (is.null(Y)){
	##One sample test on X
	q <- (colSums(X^2) != 0);
	p = sum(q);
	if (p > 1){
		X2 <- X[,q];
		Xbar <- colMeans(X2);
		diag.S <- apply(X2, 2, var);
		diag.S <- diag.S + (diag.S == 0);
		diag.S.inv <- (1/diag.S)*(diag.S >= .Machine$double.eps);
		R <- cor(X);
		Result[1] <- (n1*sum(Xbar^2*diag.S.inv) - (n1-1)*p/(n1-3))/sqrt(2*sum(R^2) - p^2/(n1-1));
		Result[2] <- 1 - pnorm(Result[1]);
	} else if (p <= 1){
		Result[1] <- 0;
		Result[2] <- 1;
	}
} else if (!is.null(Y)){
	if (ncol(X) != ncol(Y)){
		stop("Dimensions do not match");
	} 

	## DIMENSION
	p <- ncol(X);

	## SAMPLE SIZE
	n1 <- nrow(X);
	n2 <- nrow(Y);
	if (n1 <= 1|| n2 <= 1){
		stop("Minimum sample size required for both groups is 2.")
	} 

	## GO CASE BY CASE DEPENDING ON THE VALUES.
	## CASE 1 - SOME OF THE VALUES IN THE TEST ARE NAs. THEN RETURN (0,1).
	## CASE 2 - ALL VALUES OF BOTH X AND Y ARE ZEROS. THEN RETURN (0,1).
	## CASE 3 - ALL VALUES OF Y ARE ZERO - ONE SAMPLE TEST FOR X.
	## CASE 4 - ALL VALUES OF X ARE ZERO - ONE SAMPLE TEST FOR Y.
	## CASE 5 - STANDARD TEST WHEN NONE OF THE ABOVE CASES HOLD.


	if (sum(is.na(X)) + sum(is.na(Y)) > 0){
		Result[1] = 0;
		Result[2] = 1;
	} else if (sum(X^2) == 0 && sum(Y^2) == 0){
	Result[1] <- 0;
	Result[2] <- 1;
	} else if (sum(X^2) != 0 && sum(Y^2) == 0){
		##One sample test on X
		q <- (colSums(X^2) != 0);
		p = sum(q);
		if (p > 1){
			X2 <- X[,q];
			Xbar <- colMeans(X2);
			diag.S <- apply(X2, 2, var);
			diag.S <- diag.S + (diag.S == 0);
			diag.S.inv <- (1/diag.S)*(diag.S >= .Machine$double.eps);
			R <- cor(X);
			Result[1] <- (n1*sum(Xbar^2*diag.S.inv) - (n1-1)*p/(n1-3))/sqrt(2*sum(R^2) - p^2/(n1-1));
			Result[2] <- 1 - pnorm(Result[1]);
		} else if (p <= 1){
			Result[1] <- 0;
			Result[2] <- 1;
		}
	} else if (sum(X^2) == 0 && sum(Y^2) != 0){
		##One sample test on Y
		q <- (colSums(Y^2) != 0);
		p = sum(q);
		if (p > 1){
			Y2 <- Y[,q];
			Ybar <- colMeans(Y2);
			diag.S <- apply(Y2, 2, var);
			diag.S <- diag.S + (diag.S == 0);
			diag.S.inv <- (1/diag.S)*(diag.S >= .Machine$double.eps);
			R <- cor(X);
			Result[1] <- (n2*sum(Ybar^2*diag.S.inv) - (n2-1)*p/(n2-3))/sqrt(2*sum(R^2) - p^2/(n2-1));
			Result[2] <- 1 - pnorm(Result[1]);
		} else if (p <= 1){
			Result[1] <- 0;
			Result[2] <- 1;
		}
	} else {
		##Two sample test based on 2013 paper
		q <- (colSums(X^2) + colSums(Y^2) != 0);
		p = sum(q);
		X2 <- X[,q];
		Y2 <- Y[,q];
		if (p > 1){
			Xbar <- colMeans(X2);
			Ybar <- colMeans(Y2);
			S1 <- var(X2); S2 <- var(Y2);
			Ds1 <- S1*diag(rep(1,p)); Ds2 <- S2*diag(rep(1,p));
			D <- Ds1/n1 + Ds2/n2;
			Dinv <- sqrt(diag(1/(diag(D) + (diag(D) == 0)))*(D >= .Machine$double.eps));
			R <- Dinv%*%((S1/n1 + S2/n2)%*%Dinv);
			T1 <- sum((Xbar - Ybar)^2*(diag(Dinv))^2) - p;
			F11 <- Dinv%*%S1%*%Dinv;
			F1 <- (sum(F11^2) - (sum(diag(F11)))^2/(n1-1))/p;
			F21 <- Dinv%*%S2%*%Dinv;
			F2 <- (sum(F21^2) - (sum(diag(F21)))^2/(n2-1))/p;
			G <- sum(diag(F11%*%F21))/p;
			cpn <- (1 + sum(R^2)/(p^(3/2)));
			Result[1] <- T1/sqrt(p*(2*F1/(n1^2) + 2*F2/(n2^2) + 4*G/(n1*n2))*cpn);
			Result[2] <- 1 - pnorm(Result[1]);
		} else if (p <= 1){
			Result[1] <- 0;
			Result[2] <- 1;
		}
	}
}
return(Result);
}
