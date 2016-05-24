cqtest <-
function(X,Y = NULL){
# X is a n1 x p matrix with columns corresponding to subjects.
# Y is a n2 x p matrix with columns corresponding to subjects.
Result = rep(0,2);
if (is.null(Y)){
	## DO A ONE-SAMPLE TEST ON X.
	p <- ncol(X);
	n1 <- nrow(X);
	if (sum(is.na(X)) > 0 || sum(X^2) == 0){
		Result[1] <- 0;
		Result[2] <- 1;
	} else if (sum(X^2) > 0){
		##One-sample test on X
		T1 <- 0;
		T2 <- 0;
		q <- 1:n1;
		for (i in q){
			for (j in q[-i]){
				T1 <- T1 + sum(X[i,]*X[j,]);
				if (n1 > 3){
					AMean <- colMeans(X[-c(i,j),]);
				} else {AMean <- X[-c(i,j),];}
				T2 <- T2 + sum(X[j,]*(X[i,] - AMean))*sum(X[i,]*(X[j,] - AMean));
			}
		}
		if (T2 != 0){
			Result[1] <- (T1/(n1*(n1-1)))/sqrt(2*T2/(n1*(n1-1))^2);
			Result[2] <- pnorm(Result[1],lower.tail = FALSE);
		} else if (T2 == 0){
			Result[1] <- 0;
			Result[2] <- 1;
		}
	}
} else if (!is.null(Y)){		
	if (ncol(X) != ncol(Y)){
		stop("Dimensions do not match");
	} 

	## Dimension
	p <- ncol(X); 
	## Sample sizes
	n1 <- nrow(X);
	n2 <- nrow(Y);
	if (n1 <= 2 || n2 <= 2){
		stop("Minimum sample size required for both groups is 3.")
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
		T1 <- 0;
		T2 <- 0;
		q <- 1:n1;
		for (i in q){
			for (j in q[-i]){
				T1 <- T1 + sum(X[i,]*X[j,]);
				if (n1 > 3){
					AMean <- colMeans(X[-c(i,j),]);
				} else {AMean <- X[-c(i,j),];}
				T2 <- T2 + sum(X[j,]*(X[i,] - AMean))*sum(X[i,]*(X[j,] - AMean));
			}
		}
		if (T2 != 0){
			Result[1] <- (T1/(n1*(n1-1)))/sqrt(2*T2/(n1*(n1-1))^2);
			Result[2] <- pnorm(Result[1],lower.tail = FALSE);
		} else if (T2 == 0){
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
				T1 <- T1 + sum(Y[i,]*Y[j,]);
				if (n2 > 3){
				AMean <- colMeans(Y[-c(i,j),]);
				} else {AMean <- Y[-c(i,j),];}
				T2 <- T2 + sum(Y[j,]*(Y[i,] - AMean))*sum(Y[i,]*(Y[j,] - AMean));
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
		T1 <- 0; D1 <- 0;
		T2 <- 0; D2 <- 0;
		T3 <- 0; D3 <- 0;
		q1 <- 1:n1; q2 <- 1:n2;
		for (i in q1){
			for (j in q1[-i]){
				T1 <- T1 + sum(X[i,]*X[j,]);
				if (n1 > 3){
					AMean <- colMeans(X[-c(i,j),]);
				} else {AMean <- X[-c(i,j),];}
				D1 <- D1 + sum(X[j,]*(X[i,] - AMean))*sum(X[i,]*(X[j,] - AMean));
			}
		}
		for (i in q2){
			for (j in q2[-i]){
				T2 <- T2 + sum(Y[i,]*Y[j,]);
				if (n2 > 3){
				AMean <- colMeans(Y[-c(i,j),]);
				} else {AMean <- Y[-c(i,j),];}
				D2 <- D2 + sum(Y[j,]*(Y[i,] - AMean))*sum(Y[i,]*(Y[j,] - AMean));
			}
		}
		for (i in q1){
			for (j in q2){
				T3 <- T3 + sum(X[i,]*Y[j,]);
				AMean1 <- colMeans(X[-i,]);
				AMean2 <- colMeans(Y[-j,]);
				D3 <- D3 + sum(Y[j,]*(X[i,] - AMean1))*sum(X[i,]*(Y[j,] - AMean2));
			}
		}
		Numer <- (T1/(n1*(n1-1)) + T2/(n2*(n2-1)) - 2*T3/(n1*n2));
		Denomin <- sqrt(2*D1/(n1*(n1-1))^2 + 2*D2/(n2*(n2-1))^2 + 4*D3/(n1*n2)^2);
		if (Denomin != 0){
			Result[1] <- Numer/Denomin;
			Result[2] <- pnorm(Result[1],lower.tail=FALSE);
		}else if (Denomin == 0){
			Result[1] <- 0;
			Result[2] <- 1;
		}
	}
}
return(Result);
}

