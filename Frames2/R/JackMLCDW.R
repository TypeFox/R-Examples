#' @name JackMLCDW
#' @aliases JackMLCDW
#' @title Confidence intervals for MLCDW estimator based on jackknife method
#' 
#' @description Calculates confidence intervals for MLCDW estimator using jackknife procedure
#' 
#' @usage JackMLCDW (ysA, ysB, pik_A, pik_B, domains_A, domains_B, xsA, xsB, x, 
#'  ind_sam, N_A, N_B, N_ab = NULL, met = "linear", conf_level, sdA = "srs", 
#'  sdB = "srs", strA = NULL, strB = NULL, clusA = NULL, clusB = NULL, 
#'  fcpA = FALSE, fcpB = FALSE)
#' @param ysA A data frame containing information about one or more factors, each one of dimension \eqn{n_A}, collected from \eqn{s_A}.
#' @param ysB A data frame containing information about one or more factors, each one of dimension \eqn{n_B}, collected from \eqn{s_B}.
#' @param pik_A A numeric vector of length \eqn{n_A} containing first order inclusion probabilities for units included in \eqn{s_A}.
#' @param pik_B A numeric vector of length \eqn{n_B} containing first order inclusion probabilities for units included in \eqn{s_B}.
#' @param domains_A A character vector of size \eqn{n_A} indicating the domain each unit from \eqn{s_A} belongs to. Possible values are "a" and "ab".
#' @param domains_B A character vector of size \eqn{n_B} indicating the domain each unit from \eqn{s_B} belongs to. Possible values are "b" and "ba".
#' @param xsA A numeric vector of length \eqn{n_A} or a numeric matrix or data frame of dimensions \eqn{n_A} x \eqn{m}, with \eqn{m} the number of auxiliary variables, containing auxiliary information in frame A for units included in \eqn{s_A}.
#' @param xsB A numeric vector of length \eqn{n_B} or a numeric matrix or data frame of dimensions \eqn{n_B} x \eqn{m}, with \eqn{m} the number of auxiliary variables, containing auxiliary information in frame B for units included in \eqn{s_B}.
#' @param x A numeric vector or length \eqn{N} or a numeric matrix or data frame of dimensions \eqn{N} x \eqn{m}, with \eqn{m} the number of auxiliary variables, containing auxiliary information for every unit in the population.
#' @param ind_sam A numeric vector of length \eqn{n = n_A + n_B} containing the identificators of units of the population (from 1 to \eqn{N}) that belongs to \eqn{s_A} or \eqn{s_B}
#' @param N_A A numeric value indicating the size of frame A
#' @param N_B A numeric value indicating the size of frame B
#' @param N_ab (Optional) A numeric value indicating the size of the overlap domain
#' @param met (Optional) A character vector indicating the distance that must be used in calibration process. Possible values are "linear", "raking" and "logit". Default is "linear".
#' @param conf_level A numeric value indicating the confidence level for the confidence intervals.
#' @param sdA (Optional) A character vector indicating the sampling design considered in frame A. Possible values are "srs" (simple random sampling without replacement), "pps" (probabilities proportional to size sampling), "str" (stratified sampling), "clu" (cluster sampling) and "strclu" (stratified cluster sampling). Default is "srs".
#' @param sdB (Optional) A character vector indicating the sampling design considered in frame B. Possible values are "srs" (simple random sampling without replacement), "pps" (probabilities proportional to size sampling), "str" (stratified sampling), "clu" (cluster sampling) and "strclu" (stratified cluster sampling). Default is "srs".
#' @param strA (Optional) A numeric vector indicating the stratum each unit in frame A belongs to, if a stratified sampling or a stratified cluster sampling has been considered in frame A.
#' @param strB (Optional) A numeric vector indicating the stratum each unit in frame B belongs to, if a stratified sampling or a stratified cluster sampling has been considered in frame B.
#' @param clusA (Optional) A numeric vector indicating the cluster each unit in frame A belongs to, if a cluster sampling or a stratified cluster sampling has been considered in frame A.
#' @param clusB (Optional) A numeric vector indicating the cluster each unit in frame B belongs to, if a cluster sampling or a stratified cluster sampling has been considered in frame B.
#' @param fcpA (Optional) A logic value indicating if a finite population correction factor should be considered in frame A. Default is FALSE.
#' @param fcpB (Optional) A logic value indicating if a finite population correction factor should be considered in frame B. Default is FALSE. 
#' @details Let suppose a non stratified sampling design in frame A and a stratified sampling design in frame B where frame has been divided into L strata and a sample of size \eqn{n_{Bl}} from the \eqn{N_{Bl}} composing the l-th stratum is selected
#' In this context, jackknife variance estimator of an estimator \eqn{\hat{Y}_c} is given by
#'  \deqn{v_J(\hat{Y}_c) = \frac{n_{A}-1}{n_{A}}\sum_{i\in s_A} (\hat{Y}_{c}^{A}(i) -\overline{Y}_{c}^{A})^2 + \sum_{l=1}^{L}\frac{n_{Bl}-1}{n_{Bl}}  \sum_{i\in s_{Bl}} (\hat{Y}_{c}^{B}(lj) -\overline{Y}_{c}^{Bl})^2}
#' with \eqn{\hat{Y}_c^A(i)} the value of estimator \eqn{\hat{Y}_c} after dropping \eqn{i-th} unit from \code{ysA} and \eqn{\overline{Y}_{c}^{A}} the mean of values \eqn{\hat{Y}_c^A(i)}.
#' Similarly, \eqn{\hat{Y}_c^B(lj)} is the value taken by \eqn{\hat{Y}_c} after dropping j-th unit of l-th from sample \code{ysB} and \eqn{\overline{Y}_{c}^{Bl}} is the mean of values \eqn{\hat{Y}_c^B(lj)}.
#' If needed, a finite population correction factor can be included in frames by replacing \eqn{\hat{Y}_{c}^{A}(i)} or \eqn{\hat{Y}_{c}^{B}(lj)} with \eqn{\hat{Y}_{c}^{A*}(i)= \hat{Y}_{c}+\sqrt{1-\overline{\pi}_A} (\hat{Y}_{c}^{A}(i) -\hat{Y}_{c})} or
#' \eqn{\hat{Y}_{c}^{B*}(lj)= \hat{Y}_{c}+\sqrt{1-\overline{\pi}_B} (\hat{Y}_{c}^{B}(lj) -\hat{Y}_{c})}, where \eqn{\overline{\pi}_A = \sum_{i \in s_A}\pi_{iA}/nA} and \eqn{\overline{\pi}_B = \sum_{j \in s_B}\pi_{jB}/nB}
#' A confidence interval for any parameter of interest, \eqn{Y} can be calculated, then, using the pivotal method.
#' @return A numeric matrix containing estimations of population total and population mean and their corresponding confidence intervals obtained through jackknife method.
#' @references Molina, D., Rueda, M., Arcos, A. and Ranalli, M. G. (2015)
#'  \emph{Multinomial logistic estimation in dual frame surveys}
#'  Statistics and Operations Research Transactions (SORT). To be printed.
#' @references Wolter, K. M. (2007)
#'  \emph{Introduction to Variance Estimation.}
#'  2nd Edition. Springer, Inc., New York.
#' @seealso \code{\link{MLCDW}}
#' @examples
#' data(DatMA)
#' data(DatMB)
#' data(DatPopM)
#' 
#' IndSample <- c(DatMA$Id_Pop, DatMB$Id_Pop)
#' N_FrameA <- nrow(DatPopM[DatPopM$Domain == "a" | DatPopM$Domain == "ab",])
#' N_FrameB <- nrow(DatPopM[DatPopM$Domain == "b" | DatPopM$Domain == "ab",])
#'
#' \donttest{
#' #Let obtain a 95% jackknife confidence interval for variable Feeding,
#' #supposing a pps sampling in frame A and a simple random sampling
#' #without replacement in frame B with no finite population correction
#' #factor in any frame.
#' JackMLCDW(DatMA$Prog, DatMB$Prog, DatMA$ProbA, DatMB$ProbB, DatMA$Domain, 
#' DatMB$Domain, DatMA$Read, DatMB$Read, DatPopM$Read, IndSample, N_FrameA, 
#' N_FrameB, conf_level = 0.95, sdA = "pps", sdB = "srs")
#' }
#' @export
JackMLCDW = function (ysA, ysB, pik_A, pik_B, domains_A, domains_B, xsA, xsB, x, ind_sam, N_A, N_B, N_ab = NULL, met = "linear", conf_level, sdA = "srs", sdB = "srs", strA = NULL, strB = NULL, clusA = NULL, clusB = NULL, fcpA = FALSE, fcpB = FALSE){

	ysA <- as.data.frame(ysA)
	ysB <- as.data.frame(ysB)
	xsA <- as.matrix(xsA)
	xsB <- as.matrix(xsB)
	x <- as.matrix(x)

	cl <- match.call()
	results <- list()
	N <- nrow(x)

	c <- ncol(ysA)

	estimation <- MLCDW(ysA, ysB, pik_A, pik_B, domains_A, domains_B, xsA, xsB, x, ind_sam, N_A, N_B, N_ab, met)[[2]] 

	for (l in 1:c) {

		ys <- factor(c(as.character(ysA[,l]),as.character(ysB[,l])))

		lev <- levels(ys)
		m <- length(lev)

		nA <- nrow(ysA)
		nB <- nrow(ysB)
		
		mat <- matrix(NA, 6, m)

		if (sdA == "str"){

			strataA <- unique(sort(strA))
			nhA <- table(strA)
			nstrataA <- length(nhA)
			YcstrataA <- matrix(0, nstrataA, m)
			nhA <- c(0, nhA)
			cnhA <- cumsum(nhA)

			for (i in 1:nstrataA){

				k <- 1
				YcA <- matrix(0, nhA[i+1], m)
				for (j in (cnhA[i]+1):cnhA[i+1]){

					YcA[k,] <- MLCDW(ysA[-j,], ysB, pik_A[-j], pik_B, domains_A[-j], domains_B, xsA[-j,], xsB, x, ind_sam[-j], N_A, N_B, N_ab, met)[[2]][[l]][1,]
					k <- k + 1
				}

				YcAMean <- matrix(colMeans(YcA), nhA[i+1], m, byrow = TRUE)

				if (fcpA)
					fA <- 1 - mean(pik_A[strA == strataA[i]])
				else
					fA <- 1

				YcstrataA[i,] <- (nhA[i+1] - 1) / nhA[i+1] * fA * colSums((YcA - YcAMean)^2)
			}
			vjA <- colSums(YcstrataA)	
		}
		else {

			if (sdA == "clu"){

				clustersA <- unique(clusA)
				probclustersA <- unique(data.frame(pik_A, clusA))[,1]
				nclustersA <- length(clustersA)
				if (nclustersA < 3)		
					stop("Number of clusters from frame A is less than 3. Variance cannot be computed.")
			
				YcA <- matrix(0, nclustersA, m)

				for (i in 1:nclustersA)

					YcA[i,] <- MLCDW(ysA[clusA %in% clustersA[-clustersA[i]],], ysB, pik_A[clusA %in% clustersA[-clustersA[i]]], pik_B, domains_A[clusA %in% clustersA[-clustersA[i]]], domains_B, xsA[clusA %in% clustersA[-clustersA[i]],], xsB, x, ind_sam[c(clusA %in% clustersA[-clustersA[i]], rep(TRUE, nB))], N_A, N_B, N_ab, met)[[2]][[l]][1,]
			
				YcAMean <- matrix(colMeans(YcA), nclustersA, m, byrow = TRUE)

				if (fcpA)
					fA <- 1 - mean(probclustersA)
				else
					fA <- 1

				vjA <- ((nclustersA - 1) / nclustersA) * fA * colSums ((YcA - YcAMean)^2)
			}
			else{
				if (sdA == "strclu"){

					strataA <- unique(sort(strA))
					nstrataA <- length(strataA)
					YcstrataA <- matrix(0, nstrataA, m)

					for (i in 1:nstrataA){

						clustersA <- unique(clusA[strA == strataA[i]])
						probclustersA <- unique(data.frame(pik_A[strA == strataA[i]], clusA[strA == strataA[i]]))[,1]
						nclustersA <- length(clustersA)
						if (nclustersA < 3)		
							stop("Number of clusters in any stratum from frame A is less than 3. Variance cannot be computed.")
						k <- 1
						YcA <- matrix(0, nclustersA, m)
						for (j in 1:nclustersA){

							YcA[k,] <- MLCDW(ysA[clusA %in% clustersA[-clustersA[j]],], ysB, pik_A[clusA %in% clustersA[-clustersA[j]]], pik_B, domains_A[clusA %in% clustersA[-clustersA[j]]], domains_B, xsA[clusA %in% clustersA[-clustersA[j]],], xsB, x, ind_sam[c(clusA %in% clustersA[-clustersA[j]], rep(TRUE, nrow(ysB)))], N_A, N_B, N_ab, met)[[2]][[l]][1,]
							k <- k + 1
						}

						YcAMean <- matrix(colMeans(YcA), nclustersA, m, byrow = TRUE)

						if (fcpA)
							fA <- 1 - mean(probclustersA)
						else
							fA <- 1
					
						YcstrataA[i,] <- (nclustersA - 1) / nclustersA * fA * colSums((YcA - YcAMean)^2)
					}
					vjA <- colSums(YcstrataA)
				}else{

					YcA <- matrix(0, nA, m)
		
					for (i in 1:nA)

						YcA[i,] <- MLCDW(ysA[-i,], ysB, pik_A[-i], pik_B, domains_A[-i], domains_B, xsA[-i,], xsB, x, ind_sam[-i], N_A, N_B, N_ab, met)[[2]][[l]][1,]

					YcAMean <- matrix(colMeans(YcA), nA, m, byrow = TRUE)
	
					if (fcpA)
						fA <- 1 - mean(pik_A)
					else
						fA <- 1

					vjA <- ((nA - 1) / nA) * fA * colSums ((YcA - YcAMean)^2)
				}
			}	
		}

		if (sdB == "str"){

			strataB <- unique(sort(strB))
			nhB <- table(strB)
			nstrataB <- length(nhB)
			YcstrataB <- matrix(0, nstrataB, m)
			nhB <- c(0,nhB)
			cnhB <- cumsum(nhB)

			for (i in 1:nstrataB){

				k <- 1
				YcB <- matrix(0, nhB[i+1], m, byrow = TRUE)
				for (j in (cnhB[i]+1):cnhB[i+1]){

					YcB[k,] <- MLCDW(ysA, ysB[-j,], pik_A, pik_B[-j], domains_A, domains_B[-j], xsA, xsB[-j,], x, ind_sam[-(nA + j)], N_A, N_B, N_ab, met)[[2]][[l]][1,]
					k <- k + 1
				}

				YcBMean <- matrix(colMeans(YcB), nhB[i+1], m, byrow = TRUE)

				if (fcpB)
					fB <- 1 - mean(pik_B[strB == strataB[i]])
				else
					fB <- 1

				YcstrataB[i,] <- (nhB[i+1] - 1) / nhB[i+1] * fB * colSums((YcB - YcBMean)^2)
			}
			vjB <- colSums(YcstrataB)	
		}
		else{
			if (sdB == "clu"){

				clustersB <- unique(clusB)
				probclustersB <- unique(data.frame(pik_B, clusB))[,1]
				nclustersB <- length(clustersB)
				if (nclustersB < 3)		
					stop("Number of clusters from frame B is less than 3. Variance cannot be computed.")
			
				YcB <- matrix(0, nclustersB, m)

				for (i in 1:nclustersB)

					YcB[i,] <- MLCDW(ysA, ysB[clusB %in% clustersB[-clustersB[i]],], pik_A, pik_B[clusB %in% clustersB[-clustersB[i]]], domains_A, domains_B[clusB %in% clustersB[-clustersB[i]]], xsA, xsB[clusB %in% clustersB[-clustersB[i]],], x, ind_sam[c(rep(TRUE, nA), clusB %in% clustersB[-clustersB[i]])], N_A, N_B, N_ab, met)[[2]][[l]][1,]
			
				YcBMean <- matrix(colMeans(YcB), nclustersB, m, byrow = TRUE)

				if (fcpB)
					fB <- 1 - mean(probclustersB)
				else
					fB <- 1

				vjB <- ((nclustersB - 1) / nclustersB) * fB * colSums ((YcB - YcBMean)^2)
			}
			else{
				if (sdB == "strclu"){

					strataB <- unique(sort(strB))
					nstrataB <- length(strataB)
					YcstrataB <- matrix(0, nstrataB, m)

					for (i in 1:nstrataB){

						clustersB <- unique(clusB[strB == strataB[i]])
						probclustersB <- unique(data.frame(pik_B[strB == strataB[i]], clusA[strB == strataB[i]]))[,1]
						nclustersB <- length(clustersB)
						if (nclustersB < 3)		
							stop("Number of clusters in any stratum from frame B is less than 3. Variance cannot be computed.")
						k <- 1
						YcB <- matrix(0, nclustersB, m)
						for (j in 1:nclustersB){

							YcB[k,] <- MLCDW(ysA, ysB[clusB %in% clustersB[-clustersB[j]],], pik_A, pik_B[clusB %in% clustersB[-clustersB[j]]], domains_A, domains_B[clusB %in% clustersB[-clustersB[j]]], xsA, xsB[clusB %in% clustersB[-clustersB[j]],], x, ind_sam[c(rep(TRUE, nrow(ysA)), clusB %in% clustersB[-clustersB[j]])], N_A, N_B, N_ab, met)[[2]][[l]][1,]
							k <- k + 1
						}

						YcBMean <- matrix(colMeans(YcB), nclustersB, m, byrow = TRUE)

						if (fcpB)
							fB <- 1 - mean(probclustersB)
						else
							fB <- 1
					
						YcstrataB[i,] <- (nclustersB - 1) / nclustersB * fB * colSums((YcB - YcBMean)^2)
					}
					vjB <- colSums(YcstrataB)
				}
				else{

					YcB <- matrix(0, nB, m)

					for (i in 1:nB)
						YcB[i,] <- MLCDW(ysA, ysB[-i,], pik_A, pik_B[-i], domains_A, domains_B[-i], xsA, xsB[-i,], x, ind_sam[-(nA + i)], N_A, N_B, N_ab, met)[[2]][[l]][1,]

					YcBMean <- matrix(colMeans(YcB), nB, m, byrow = TRUE)

					if (fcpB)
						fB <- 1 - mean(pik_B)
					else
						fB <- 1

					vjB <- ((nB - 1) / nB) * fB * colSums((YcB - YcBMean)^2)

				}
			}
		}

		VJack_Phat_MLCDW <- vjA + vjB

		mat[1,] <- estimation[[l]][1,]
		mat[2,] <- estimation[[l]][1,] + qnorm(1 - (1 - conf_level) / 2) * sqrt(VJack_Phat_MLCDW)
		mat[3,] <- estimation[[l]][1,] - qnorm(1 - (1 - conf_level) / 2) * sqrt(VJack_Phat_MLCDW)
		mat[4,] <- estimation[[l]][2,]
		mat[5,] <- estimation[[l]][2,] + qnorm(1 - (1 - conf_level) / 2) * sqrt(1/N^2 * VJack_Phat_MLCDW)
		mat[6,] <- estimation[[l]][2,] - qnorm(1 - (1 - conf_level) / 2) * sqrt(1/N^2 * VJack_Phat_MLCDW)
	
		results[[l]] <- mat
	}

	return(results)
}