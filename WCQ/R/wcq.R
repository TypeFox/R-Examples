onesamplemeantest <- function(X1) {
	#library(msm);
	n1 <- dim(X1)[2];
	p <- dim(X1)[1];
	
	#write("Chen-Qin Two Sample Mean Test:", "");
	#write(paste("X1_n =", n1), "");
	#write(paste("X2_n =", n2), "");
	#write(paste("p =", p), "");
	
	exp1 <- 0;
	T_n <- 0;
	F <- 0;
	pval <- 0;

	for (i in 1:n1) {	
		for (j in 1:n1) {
			if (i != j) {
				term1 <- t(X1[,i])%*%X1[,j];	
				exp1 <- exp1 + term1[1];
			}
		}
	}
	
	#write(paste("expression1 =", exp1), "");

	
	#write(paste("expression3 =", exp3), "");
	
	F_n <- (exp1/(n1*(n1-1)));
	#write(paste("T_n =", T_n), "");

	trS1 <- 0;
	sigma_n1 <- 0;
	
	for (j in 1:n1)	{
		for (k in 1:n1) {
			if (j != k) {
				tempMean <- (rowSums(X1) - X1[,j] - X1[,k])/(n1 - 2);
				term1 <- (X1[,j] - tempMean);
				term1 <- term1 %*% t(X1[,j]);
				term1 <- term1 %*% (X1[,k] - tempMean);
				term1 <- term1 %*% t(X1[,k]);
				trS1 <- trS1 + term1;
			}
		}
	}
	trS1 <- sum(diag(trS1));
	trS1 <- trS1/(n1 * (n1 - 1));

	
	sigma_n1 <- ((2 / (n1 * (n1 - 1))) * trS1);
	sigma_n1 <- sqrt(sigma_n1);
	
	
	Q_n <- F_n/sigma_n1;
	#pval <- ptnorm(Q_n, lower=0, upper=2);
	pval <- 1 - pnorm(Q_n);
	#pval <- pt(Q_n, df=10000);
	return(pval);
}

twosamplemeantest <- function(X1,X2) {
	#library(msm);
	n1 <- dim(X1)[2];
	n2 <- dim(X2)[2];
	p <- dim(X1)[1];
	
	#write("Chen-Qin Two Sample Mean Test:", "");
	#write(paste("X1_n =", n1), "");
	#write(paste("X2_n =", n2), "");
	#write(paste("p =", p), "");
	
	exp1 <- 0;
	exp2 <- 0;
	exp3 <- 0;
	T_n <- 0;
	F <- 0;
	pval <- 0;

	for (i in 1:n1) {	
		for (j in 1:n1) {
			if (i != j) {
				term1 <- t(X1[,i])%*%X1[,j];	
				exp1 <- exp1 + term1[1];
			}
		}
	}
	
	#write(paste("expression1 =", exp1), "");

	for (i in 1:n2) {
		for (j in 1:n2) {
			if (i != j) {
				term2 <- t(X2[,i])%*%X2[,j];	
				exp2 <- exp2 + term2[1];
			}
		}
	}
	
	#write(paste("expression2 =", exp2), "");
	
	for (i in 1:n1) {
		for (j in 1:n2) {
			term3 <- t(X1[,i])%*%X2[,j];
			exp3 <- exp3 + term3[1];
		}
	}
	
	#write(paste("expression3 =", exp3), "");
	
	T_n <- (exp1/(n1*(n1-1))) + (exp2/(n2*(n2-1))) - (2*(exp3/(n1*n2)));
	#write(paste("T_n =", T_n), "");

	trS1 <- 0;
	trS2 <- 0;
	trS1S2 <- 0;
	sigma_n1 <- 0;
	
	for (j in 1:n1)	{
		for (k in 1:n1) {
			if (j != k) {
				tempMean <- (rowSums(X1) - X1[,j] - X1[,k])/(n1 - 2);
				term1 <- (X1[,j] - tempMean);
				term1 <- term1 %*% t(X1[,j]);
				term1 <- term1 %*% (X1[,k] - tempMean);
				term1 <- term1 %*% t(X1[,k]);
				trS1 <- trS1 + term1;
			}
		}
	}
	trS1 <- sum(diag(trS1));
	trS1 <- trS1/(n1 * (n1 - 1));
	#write(paste("trS1 =", trS1), "");

	for (j in 1:n2)	{
		for (k in 1:n2) {
			if (j != k) {
				tempMean <- (rowSums(X2) - X2[,j] - X2[,k])/(n2 - 2);
				term2 <- (X2[,j] - tempMean);
				term2 <- term2 %*% t(X2[,j]);
				term2 <- term2 %*% (X2[,k] - tempMean);
				term2 <- term2 %*% t(X2[,k]);
				trS2 <- trS2 + term2;
			}
		}
	}
	trS2 <- sum(diag(trS2));
	trS2 <- trS2/(n2 * (n2 - 1));
	#write(paste("trS2 =", trS2), "");

	for (j in 1:n1)	{
		for (k in 1:n2) {
			tempMean1 <- (rowSums(X1) - X1[,j])/(n1 - 1);
			tempMean2 <- (rowSums(X2) - X2[,k])/(n2 - 1);
			term3 <- (X1[,j] - tempMean1);
			term3 <- term3 %*% t(X1[,j]);
			term3 <- term3 %*% (X2[,k] - tempMean2);
			term3 <- term3 %*% t(X2[,k]);
			trS1S2 <- trS1S2 + term3;
		}
	}
	trS1S2 <- sum(diag(trS1S2));
	trS1S2 <- trS1S2/(n1 * n2);
	#write(paste("trS1S2 =", trS1S2), "");
	
	sigma_n1 <- ((2 / (n1 * (n1 - 1))) * trS1) + 
				((2 / (n2 * (n2 - 1))) * trS2) + 
				((4 / (n1 * n2)) * trS1S2);
	sigma_n1 <- sqrt(sigma_n1);
	#write(paste("sigma_n1 =", sigma_n1), "");
	
	#psi_alpha <- qnorm(0.95);
	#write(paste("psi_alpha (at alpha = 0.05) =", psi_alpha), "");
	
	Q_n <- T_n/sigma_n1;
	#write(paste("Q_n =", Q_n), "");
	#pval <- ptnorm(Q_n, lower=0, upper=2);
	pval <- 1 - pnorm(Q_n);	
	#pval <- pt(Q_n, df=10000);
	return(pval);
}

wcq <- function(marker_data, trait_data, alleles) {
#	library(msm);
#	source("twosamplemeantest.R");
#	source("onesamplemeantest.R");
#	source("twosamplemeantest_qn.R");
#	source("onesamplemeantest_qn.R");
	num_samples <- dim(marker_data)[2];
	num_markers <- dim(marker_data)[1];
	num_traits <- dim(trait_data)[1];
	pval_matrix <- array(0, dim=c(num_markers, num_traits));

	#allele_sample_prob <- array(0, dim=c(length(alleles),num_samples));
	#for (i in 1:length(alleles)) {
	#	for (j in 1:num_samples) {
	#		allele_sample_prob[i,j] <- sum(marker_data[,j]==alleles[i]);
	#	}
	#}	
	#allele_sample_prob <- allele_sample_prob / num_markers;
		
	for (i in 1:num_markers) {		
		for (j in 1:num_traits) {
			current_trait <- trait_data[j,];
			mean_current_trait <- mean(current_trait);
			stdev_current_trait <- sqrt(var(current_trait));
			X_indices_1 <- which(marker_data[i,]==alleles[1], arr.ind=T);
			#current_trait[X_indices_1] <- current_trait[X_indices_1] * (allele_sample_prob[1,X_indices_1] * length(X_indices_1) / num_samples);
			X_indices_2 <- which(marker_data[i,]==alleles[3], arr.ind=T);
			#current_trait[X_indices_2] <- 0.75 * current_trait[X_indices_2];
			X_indices_tot <- c(X_indices_1, X_indices_2);
			#X_indices_tot <- c(X_indices_1);
			X <- as.vector(current_trait)[X_indices_tot];
			X <- X / pnorm(X, mean=mean_current_trait, sd=stdev_current_trait);
			#X <- X / pnorm(X);
			
			Y_indices <- which(marker_data[i,]==alleles[2], arr.ind=T);
			#Y_indices <- c(Y_indices, X_indices_2);
			#current_trait[Y_indices] <- current_trait[Y_indices] * (allele_sample_prob[2,Y_indices] * length(Y_indices) / num_samples);
			Y <- as.vector(current_trait)[Y_indices];
			Y <- Y / pnorm(Y, mean=mean_current_trait, sd=stdev_current_trait);
			#Y <- Y / pnorm(Y);
			
			pval <- 0;			
			X_arr <- array(X, dim=c(1,length(X_indices_tot)));
  			Y_arr <- array(Y, dim=c(1,length(Y_indices)));
			if ((length(X) > 1) && (length(Y) > 1)) {
				pval <- twosamplemeantest(X_arr,Y_arr);
				if (is.nan(pval)) { #Chen-Qin test "failed"
					if (length(X) > length(Y)) {
						X_arr_normalized <- X_arr - array(mean(X_arr), c(1, length(X_indices_tot)));
						pval <- onesamplemeantest(X_arr_normalized);
					} else if (length(Y) > length(X)){
						Y_arr_normalized <- Y_arr - array(mean(Y_arr), c(1, length(Y_indices)));
						pval <- onesamplemeantest(Y_arr_normalized);
					} else {
						XY_diff <- X_arr - Y_arr;
						pval <- onesamplemeantest(XY_diff);
					}
				}
			} else { #only one group has at least sample size = 2
				if (length(X) > 1) {
					X_arr_normalized <- X_arr - array(mean(X_arr), c(1, length(X_indices_tot)));
					pval <- onesamplemeantest(X_arr_normalized);
				} else if (length(Y) > 1){
					Y_arr_normalized <- Y_arr - array(mean(Y_arr), c(1, length(Y_indices)));
					pval <- onesamplemeantest(Y_arr_normalized);
				}
			}
			if (is.nan(pval)) { #if all possble tests fail
				pval <- 0;
			}
			# transform and express the result in terms of log(p)
			#pval <- 1 - pval;
			#if (pval != 0) {
			#	pval = -log10(pval);
			#}
			pval_matrix[i,j] <- pval;
		}
	}
	dimnames(pval_matrix)[[1]] <- dimnames(marker_data)[[1]];
	dimnames(pval_matrix)[[2]] <- dimnames(trait_data)[[1]];
	return(pval_matrix);
}

wcq.fdrc <- function(marker_data, trait_data, alleles, fdrc.method="none") {
	pval_matrix <- wcq(marker_data, trait_data, alleles);
	if (fdrc.method == "none"){
		return(pval_matrix);
	} else {
		adjusted_pval_matrix <- array(0, dim=dim(pval_matrix));
		for (i in 1:dim(pval_matrix)[2]) {
			adjusted_pval_matrix[,i] <- p.adjust(pval_matrix[,i], fdrc.method);		
		} 
		dimnames(adjusted_pval_matrix) <- dimnames(pval_matrix);
		return(adjusted_pval_matrix);
	}
}

detect.qtl <- function(marker_data, trait_data, alleles, fdrc.method="none", threshold=0.05) {
	pval_matrix <- wcq.fdrc(marker_data, trait_data, alleles, fdrc.method);
	qtl.list <- list();
	names(qtl.list) <- dimnames(trait_data)[[1]];
	for (i in 1:dim(trait_data)[1]) {
		qtl.list[[i]] <- pval_matrix[(pval_matrix[,i] < threshold),i];
	}
	return(qtl.list);
}