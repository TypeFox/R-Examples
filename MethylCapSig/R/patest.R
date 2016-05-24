patest <-
function(X,Y = NULL){
# X is a p x n1 matrix with columns corresponding to subjects.
# Y is a p x n2 matrix with columns corresponding to subjects.
Result = rep(0,2);
if (is.null(Y)){
	##One-sample test on X
	T1 <- 0; S <- 0;
	T2 <- 0;
	q <- 1:n1;
	for (i in q){
		for (j in q[-i]){
			diag.S <- apply(X[-c(i,j),], 2, var);
			diag.S <- diag.S + (diag.S == 0);
			diag.S.inv <- (1/diag.S)*(diag.S >= .Machine$double.eps);
			T1 <- T1 + sum(X[i,]*X[j,]*diag.S.inv);
			XMean1 <-  colMeans(X[-c(i,j),]);
			T2 <- T2 + sum(X[i,]*(X[j,] - XMean1)*diag.S.inv)*sum(X[j,]*(X[i,] - XMean1)*diag.S.inv);
		}
	}
	if (T2 != 0){
		Result[1] <- ((n1-5)/(n1*(n1-1)*(n1-3))*T1)/sqrt(2/(n1*(n1-1))*T2);
		Result[2] <- pnorm(Result[1],lower.tail = FALSE);
	}else if (T2 == 0){
		Result[1] <- 0;
		Result[2] <- 1;
	}
}else if (!is.null(Y)){
	if (ncol(X) != ncol(Y)){
		stop("Dimensions do not match");
	} 

	## DIMENSION
	p <- ncol(X);

	## SAMPLE SIZE
	n1 <- nrow(X);
	n2 <- nrow(Y);
	if (n1 <= 3 || n2 <= 3){
		stop("Minimum sample size required for both groups is 4.")
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
		##One-sample test on X
		T1 <- 0; S <- 0;
		T2 <- 0;
		q <- 1:n1;
		for (i in q){
			for (j in q[-i]){
				diag.S <- apply(X[-c(i,j),], 2, var);
				diag.S <- diag.S + (diag.S == 0);
				diag.S.inv <- (1/diag.S)*(diag.S >= .Machine$double.eps);
				T1 <- T1 + sum(X[i,]*X[j,]*diag.S.inv);
				XMean1 <-  colMeans(X[-c(i,j),]);
				T2 <- T2 + sum(X[i,]*(X[j,] - XMean1)*diag.S.inv)*sum(X[j,]*(X[i,] - XMean1)*diag.S.inv);
			}
		}
		if (T2 != 0){
			Result[1] <- ((n1-5)/(n1*(n1-1)*(n1-3))*T1)/sqrt(2/(n1*(n1-1))*T2);
			Result[2] <- pnorm(Result[1],lower.tail = FALSE);
		}else if (T2 == 0){
			Result[1] <- 0;
			Result[2] <- 1;
		}
	} else if (sum(X^2) == 0 && sum(Y^2) != 0){
		##One-sample test on Y
		T1 <- 0;
		T2 <- 0;
		q <- 1:n2;
		for (i in q){
			for (j in q[-i]){
			diag.S <- apply(Y[-c(i,j),], 2, var);
			diag.S <- diag.S + (diag.S == 0);
			diag.S.inv <- (1/diag.S)*(diag.S >= .Machine$double.eps);
			T1 <- T1 + sum(Y[,i]*Y[,j]*diag.S.inv);
			YMean1 <-  colMeans(Y[-c(i,j),]);
			T2 <- T2 + sum(Y[,i]*(Y[,j] - YMean1)*diag.S.inv)*sum(Y[,j]*(Y[,i] - YMean1)*diag.S.inv);
			}
		}
		if (T2 != 0){
			Result[1] <- (T1/(n1*(n1-1)))/sqrt(2*T2/(n1*(n1-1))^2);
			Result[2] <- pnorm(Result[1],lower.tail = FALSE);
		}else if (T2 == 0){
			Result[1] <- 0;
			Result[2] <- 1;
		}
	} else {
		##Two-sample test on X and Y
		T1 <- 0; T2 <- 0; T3 <- 0;  ## Numerator terms
		D1 <- 0; D2 <- 0; D3 <- 0;  ## Denominator terms
		S1 <- 0; S2 <- 0; S3 <- 0;  ## Covariance matrices
		q1 <- 1:n1; q2 <- 1:n2;
		for (i in q1){
			for (j  in q1[-i]){
				diag.S1 <- ((n1-3)*apply(X[-c(i,j),], 2, var) + (n2 - 1)*apply(Y, 2, var))/(n1 + n2 - 4);
				diag.S1 <- diag.S1 + (diag.S1 == 0);
				diag.S1.inv <- (1/diag.S1)*(diag.S1 >= .Machine$double.eps);
				T1 <- T1 + sum(X[i,]*X[j,]*diag.S1.inv);
				XMean1 <- colMeans(X[-c(i,j),]);
				D1 <- D1 + sum(X[i,]*(X[j,] - XMean1)*diag.S1.inv)*sum(X[j,]*(X[i,] - XMean1)*diag.S1.inv);
			}
		}
		for (i in q2){
			for (j in q2[-i]){
				diag.S2 <- ((n1 - 1)*apply(X, 2, var) + (n2 - 3)*apply(Y[-c(i,j),], 2, var))/(n1 + n2 - 4);
				diag.S2 <- diag.S2 + (diag.S2 == 0);
				diag.S2.inv <- (1/diag.S2)*(diag.S2 >= .Machine$double.eps);
				T2 <- T2 + sum(Y[i,]*Y[j,]*diag.S2.inv);
				YMean1 <- colMeans(Y[-c(i,j),]);
				D2 <- D2 + sum(Y[i,]*(Y[j,] - YMean1)*diag.S2.inv)*sum(Y[j,]*(Y[i,] - YMean1)*diag.S2.inv);
			}
		}		       
		for (i in q1){
			for (j in q2){
				diag.S3 <- ((n1 - 2)*apply(X[-i,], 2, var) + (n2 - 2)*apply(Y[-j,], 2, var))/(n1 + n2 - 4);
				diag.S3 <- diag.S3 + (diag.S3 == 0);
				diag.S3.inv <- (1/diag.S3)*(diag.S3 >= .Machine$double.eps);
				T3 <- T3 + sum(X[i,]*Y[j,]*diag.S3.inv);
				XMean1 <- colMeans(X[-i,]);
				YMean1 <- colMeans(Y[-j,]);
				D3 <- D3 + sum(X[i,]*(Y[j,] - YMean1)*diag.S3.inv)*sum(Y[j,]*(X[i,] - XMean1)*diag.S3.inv);
			}
		}
		Numer = (T1/(n1*(n1-1)) + T2/(n2*(n2-1)) - (2*T3)/(n1*n2));
		Denomin = (2*D1)/(n1*(n1-1))^2 + (2*D2)/(n2*(n2-1))^2 + (4*D3)/(n1*n2)^2;
		if (Denomin > 0){
			Result[1] <- Numer/sqrt(Denomin);
			Result[2] <- pnorm(Result[1], lower.tail = FALSE);               
		}
		if (Denomin <= 0){
			Result[1] <- NA;
			Result[2] <- 1;
		}
	}
}
return(Result);
}
