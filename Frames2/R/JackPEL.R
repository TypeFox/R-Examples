#' @name JackPEL
#' @aliases JackPEL
#' @title Confidence intervals for the pseudo empirical likelihood estimator based on jackknife method
#' 
#' @description Calculates confidence intervals for pseudo empirical likelihood estimator using jackknife procedure
#' 
#' @usage JackPEL(ysA, ysB, piA, piB, domainsA, domainsB, N_A = NULL, N_B = NULL, 
#' N_ab = NULL, xsAFrameA = NULL, xsBFrameA = NULL, xsAFrameB = NULL, xsBFrameB = NULL, 
#' XA = NULL, XB = NULL, conf_level, sdA = "srs", sdB = "srs", strA = NULL, strB = NULL, 
#' clusA = NULL,clusB = NULL, fcpA = FALSE, fcpB = FALSE)
#' @param ysA A numeric vector of length \eqn{nA} or a numeric matrix or data frame of dimensions \eqn{nA} x \eqn{c} containing information about variable of interest from \eqn{s_A}.
#' @param ysB A numeric vector of length \eqn{nB} or a numeric matrix or data frame of dimensions \eqn{nB} x \eqn{c} containing information about variable of interest from \eqn{s_B}.
#' @param piA A numeric vector of length \eqn{nA} or a square numeric matrix of dimension \eqn{nA} containing first order or first and second order inclusion probabilities for units included in \eqn{s_A}.
#' @param piB A numeric vector of length \eqn{nB} or a square numeric matrix of dimension \eqn{nB} containing first order or first and second order inclusion probabilities for units included in \eqn{s_B}.
#' @param domainsA A character vector of size \eqn{nA} indicating the domain each unit from \eqn{s_A} belongs to. Possible values are "a" and "ab".
#' @param domainsB A character vector of size \eqn{nB} indicating the domain each unit from \eqn{s_B} belongs to. Possible values are "b" and "ba".
#' @param N_A (Optional) A numeric value indicating the size of frame A
#' @param N_B (Optional) A numeric value indicating the size of frame B
#' @param N_ab (Optional) A numeric value indicating the size of the overlap domain
#' @param xsAFrameA (Optional) A numeric vector of length \eqn{nA} or a numeric matrix or data frame of dimensions \eqn{nA} x \eqn{m_A}, with \eqn{m_A} the number of auxiliary variables in frame A, containing auxiliary information in frame A for units included in \eqn{s_A}.
#' @param xsBFrameA (Optional) A numeric vector of length \eqn{nB} or a numeric matrix or data frame of dimensions \eqn{nB} x \eqn{m_A}, with \eqn{m_A} the number of auxiliary variables in frame A, containing auxiliary information in frame A for units included in \eqn{s_B}. For units in domain \eqn{b}, these values are 0.
#' @param xsAFrameB (Optional) A numeric vector of length \eqn{nA} or a numeric matrix or data frame of dimensions \eqn{nA} x \eqn{m_B}, with \eqn{m_B} the number of auxiliary variables in frame B, containing auxiliary information in frame B for units included in \eqn{s_A}. For units in domain \eqn{a}, these values are 0.
#' @param xsBFrameB (Optional) A numeric vector of length \eqn{nB} or a numeric matrix or data frame of dimensions \eqn{nB} x \eqn{m_B}, with \eqn{m_B} the number of auxiliary variables in frame B, containing auxiliary information in frame B for units included in \eqn{s_B}.
#' @param XA (Optional) A numeric value or vector of length \eqn{m_A}, with \eqn{m_A} the number of auxiliary variables in frame A, indicating the population totals for the auxiliary variables considered in frame A.
#' @param XB (Optional) A numeric value or vector of length \eqn{m_B}, with \eqn{m_B} the number of auxiliary variables in frame B, indicating the population totals for the auxiliary variables considered in frame B.
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
#' In this context, jackknife variance estimator of a estimator \eqn{\hat{Y}_c} is given by
#'  \deqn{v_J(\hat{Y}_c) = \frac{n_{A}-1}{n_{A}}\sum_{i\in s_A} (\hat{Y}_{c}^{A}(i) -\overline{Y}_{c}^{A})^2 + \sum_{l=1}^{L}\frac{n_{Bl}-1}{n_{Bl}}  \sum_{i\in s_{Bl}} (\hat{Y}_{c}^{B}(lj) -\overline{Y}_{c}^{Bl})^2}
#' with \eqn{\hat{Y}_c^A(i)} the value of estimator \eqn{\hat{Y}_c} after dropping \eqn{i-th} unit from \code{ysA} and \eqn{\overline{Y}_{c}^{A}} the mean of values \eqn{\hat{Y}_c^A(i)}.
#' Similarly, \eqn{\hat{Y}_c^B(lj)} is the value taken by \eqn{\hat{Y}_c} after dropping j-th unit of l-th from sample \code{ysB} and \eqn{\overline{Y}_{c}^{Bl}} is the mean of values \eqn{\hat{Y}_c^B(lj)}.
#' If needed, a finite population correction factor can be included in frames by replacing \eqn{\hat{Y}_{c}^{A}(i)} or \eqn{\hat{Y}_{c}^{B}(lj)} with \eqn{\hat{Y}_{c}^{A*}(i)= \hat{Y}_{c}+\sqrt{1-\overline{\pi}_A} (\hat{Y}_{c}^{A}(i) -\hat{Y}_{c})} or
#' \eqn{\hat{Y}_{c}^{B*}(lj)= \hat{Y}_{c}+\sqrt{1-\overline{\pi}_B} (\hat{Y}_{c}^{B}(lj) -\hat{Y}_{c})}, where \eqn{\overline{\pi}_A = \sum_{i \in s_A}\pi_{iA}/nA} and \eqn{\overline{\pi}_B = \sum_{j \in s_B}\pi_{jB}/nB}
#' A confidence interval for any parameter of interest, \eqn{Y} can be calculated, then, using the pivotal method.
#' @return A numeric matrix containing estimations of population total and population mean and their corresponding confidence intervals obtained through jackknife method.
#' @references Wolter, K. M. (2007)
#'  \emph{Introduction to Variance Estimation.}
#'  2nd Edition. Springer, Inc., New York.
#' @seealso \code{\link{PEL}}
#' @export
JackPEL = function (ysA, ysB, piA, piB, domainsA, domainsB, N_A = NULL, N_B = NULL, N_ab = NULL, xsAFrameA = NULL, xsBFrameA = NULL, xsAFrameB = NULL, xsBFrameB = NULL, XA = NULL, XB = NULL, conf_level, sdA = "srs", sdB = "srs", strA = NULL, strB = NULL, clusA = NULL, clusB = NULL, fcpA = FALSE, fcpB = FALSE){

	message("Execution of this function may take several minutes. Please, wait.","\n",appendLF=FALSE)
  	flush.console()	

	cnames <- names(ysA)

	ysA <- as.matrix(ysA)
	ysB <- as.matrix(ysB)
	piA <- as.matrix(piA)
	piB <- as.matrix(piB)
	if (!is.null(xsAFrameA)){
		xsAFrameA <- as.matrix(xsAFrameA)
		xsBFrameA <- as.matrix(xsBFrameA)
	}
	if (!is.null(xsBFrameB)){
		xsAFrameB <- as.matrix(xsAFrameB)
		xsBFrameB <- as.matrix(xsBFrameB)
	}

	c <- ncol(ysA)
	results <- matrix(NA, nrow = 6, ncol = c)
	rownames(results) <- c("Total", "Jack Upper End", "Jack Lower End", "Mean", "Jack Upper End", "Jack Lower End")
	colnames(results) <- cnames

	estimation <- PEL(ysA, ysB, piA, piB, domainsA, domainsB, N_A, N_B, N_ab, xsAFrameA, xsBFrameA, xsAFrameB, xsBFrameB, XA, XB)
	size_estimation <- estimation[[2]][1,1] / estimation[[2]][2,1]

	if (!is.null(dim(drop(piA))))
		pikA <- diag(piA)
	else
		pikA <- piA

	if (!is.null(dim(drop(piB))))
		pikB <- diag(piB)
	else
		pikB <- piB

	if (sdA == "str"){

		strataA <- unique(strA)
		nhA <- table(strA)
		nstrataA <- length(nhA)
		YcstrataA <- matrix(0, nstrataA, c)
		nhA <- c(0, nhA)
		cnhA <- cumsum(nhA)

		for (i in 1:nstrataA){

			k <- 1
			YcA <- matrix(0, nhA[i+1], c)
			for (j in (cnhA[i]+1):cnhA[i+1]){

				if (!is.null(dim(drop(piA))))
					YcA[k,] <- PEL(ysA[-j,], ysB, piA[-j,-j], piB, domainsA[-j], domainsB, N_A, N_B, N_ab, xsAFrameA[-j,], xsBFrameA, xsAFrameB[-j,], xsBFrameB, XA, XB)[[2]][1,]
				else
					YcA[k,] <- PEL(ysA[-j,], ysB, piA[-j], piB, domainsA[-j], domainsB, N_A, N_B, N_ab, xsAFrameA[-j,], xsBFrameA, xsAFrameB[-j,], xsBFrameB, XA, XB)[[2]][1,]
				k <- k + 1
			}

			YcAMean <- matrix(colMeans(YcA), nhA[i+1], c, byrow = TRUE)

			if (fcpA)
				fA <- 1 - mean(pikA[strA == strataA[i]])
			else
				fA <- 1

			YcstrataA[i,] <- (nhA[i+1] - 1) / nhA[i+1] * fA * colSums((YcA - YcAMean)^2)
		}
		vjA <- colSums(YcstrataA)	
	}
	else {

		if (sdA == "clu"){

			clustersA <- unique(clusA)
			probclustersA <- unique(data.frame(pikA, clusA))[,1]
			nclustersA <- length(clustersA)
			if (nclustersA < 3)		
				stop("Number of clusters from frame A is less than 3. Variance cannot be computed.")
			
			YcA <- matrix(0, nclustersA, c)

			for (i in 1:nclustersA){
				if (!is.null(dim(drop(piA))))
					YcA[i,] <- PEL(ysA[clusA %in% clustersA[-clustersA[i]],], ysB, piA[clusA %in% clustersA[-clustersA[i]],clusA %in% clustersA[-clustersA[i]]], piB, domainsA[clusA %in% clustersA[-clustersA[i]]], domainsB, N_A, N_B, N_ab, xsAFrameA[clusA %in% clustersA[-clustersA[i]],], xsBFrameA, xsAFrameB[clusA %in% clustersA[-clustersA[i]],], xsBFrameB, XA, XB)[[2]][1,]
				else
					YcA[i,] <- PEL(ysA[clusA %in% clustersA[-clustersA[i]],], ysB, piA[clusA %in% clustersA[-clustersA[i]]], piB, domainsA[clusA %in% clustersA[-clustersA[i]]], domainsB, N_A, N_B, N_ab, xsAFrameA[clusA %in% clustersA[-clustersA[i]],], xsBFrameA, xsAFrameB[clusA %in% clustersA[-clustersA[i]],], xsBFrameB, XA, XB)[[2]][1,]
			}

			YcAMean <- matrix(colMeans(YcA), nclustersA, c, byrow = TRUE)

			if (fcpA)
				fA <- 1 - mean(probclustersA)
			else
				fA <- 1

			vjA <- ((nclustersA - 1) / nclustersA) * fA * colSums ((YcA - YcAMean)^2)
		}
		else{
			if (sdA == "strclu"){

				strataA <- unique(strA)
				nstrataA <- length(strataA)
				nhA <- table(strA)
				YcstrataA <- matrix(0, nstrataA, c)
				nhA <- c(0,nhA)
				cnhA <- cumsum(nhA)

				for (i in 1:nstrataA){

					clustersA <- unique(clusA[strA == strataA[i]])
					probclustersA <- unique(data.frame(pikA[strA == strataA[i]], clusA[strA == strataA[i]]))[,1]
					nclustersA <- length(clustersA)
					if (nclustersA < 3)		
						stop("Number of clusters in any stratum from frame A is less than 3. Variance cannot be computed.")
					k <- 1
					YcA <- matrix(0, nclustersA, c)
					for (j in 1:nclustersA){

						if (!is.null(dim(drop(piA))))
							YcA[k,] <- PEL(ysA[clusA %in% clustersA[-clustersA[j]],], ysB, piA[clusA %in% clustersA[-clustersA[j]],clusA %in% clustersA[-clustersA[j]]], piB, domainsA[clusA %in% clustersA[-clustersA[j]]], domainsB, N_A, N_B, N_ab, xsAFrameA[clusA %in% clustersA[-clustersA[j]],], xsBFrameA, xsAFrameB[clusA %in% clustersA[-clustersA[j]],], xsBFrameB, XA, XB)[[2]][1,]
						else
							YcA[k,] <- PEL(ysA[clusA %in% clustersA[-clustersA[j]],], ysB, piA[clusA %in% clustersA[-clustersA[j]]], piB, domainsA[clusA %in% clustersA[-clustersA[j]]], domainsB, N_A, N_B, N_ab, xsAFrameA[clusA %in% clustersA[-clustersA[j]],], xsBFrameA, xsAFrameB[clusA %in% clustersA[-clustersA[j]],], xsBFrameB, XA, XB)[[2]][1,]
						k <- k + 1
					}

					YcAMean <- matrix(colMeans(YcA), nclustersA, c, byrow = TRUE)

					if (fcpA)
						fA <- 1 - mean(probclustersA)
					else
						fA <- 1
					
					YcstrataA[i,] <- (nclustersA - 1) / nclustersA * fA * colSums((YcA - YcAMean)^2)
				}
				vjA <- colSums(YcstrataA)
			}else{
				nA <- nrow(ysA)
				YcA <- matrix(0, nA, c)
	
				for (i in 1:nA){
					if (!is.null(dim(drop(piA))))
						YcA[i,] <- PEL(ysA[-i,], ysB, piA[-i,-i], piB, domainsA[-i], domainsB, N_A, N_B, N_ab, xsAFrameA[-i,], xsBFrameA, xsAFrameB[-i,], xsBFrameB, XA, XB)[[2]][1,]
					else
						YcA[i,] <- PEL(ysA[-i,], ysB, piA[-i], piB, domainsA[-i], domainsB, N_A, N_B, N_ab, xsAFrameA[-i,], xsBFrameA, xsAFrameB[-i,], xsBFrameB, XA, XB)[[2]][1,]		
				}

				YcAMean <- matrix(colMeans(YcA), nA, c, byrow = TRUE)
	
				if (fcpA)
					fA <- 1 - mean(pikA)
				else
					fA <- 1

				vjA <- ((nA - 1) / nA) * fA * colSums ((YcA - YcAMean)^2)
			}	
		}
	}
	if (sdB == "str"){

		strataB <- unique(strB)
		nhB <- table(strB)
		nstrataB <- length(nhB)
		YcstrataB <- matrix(0, nstrataB, c)
		nhB <- c(0,nhB)
		cnhB <- cumsum(nhB)

		for (i in 1:nstrataB){

			k <- 1
			YcB <- matrix(0, nhB[i+1], c)
			for (j in (cnhB[i]+1):cnhB[i+1]){

				if (!is.null(dim(drop(piB))))
					YcB[k,] <- PEL(ysA, ysB[-j,], piA, piB[-j,-j], domainsA, domainsB[-j], N_A, N_B, N_ab, xsAFrameA, xsBFrameA[-j,], xsAFrameB, xsBFrameB[-j,], XA, XB)[[2]][1,]		
				else
					YcB[k,] <- PEL(ysA, ysB[-j,], piA, piB[-j], domainsA, domainsB[-j], N_A, N_B, N_ab, xsAFrameA, xsBFrameA[-j,], xsAFrameB, xsBFrameB[-j,], XA, XB)[[2]][1,]
				k <- k + 1
			}

			YcBMean <- matrix(colMeans(YcB), nhB[i+1], c, byrow = TRUE)

			if (fcpB)
				fB <- 1 - mean(pikB[strB == strataB[i]])
			else
				fB <- 1
			YcstrataB[i,] <- (nhB[i+1] - 1) / nhB[i+1] * fB * colSums((YcB - YcBMean)^2)
		}
		vjB <- colSums(YcstrataB)	
	}
	else{

		if (sdB == "clu"){

			clustersB <- unique(clusB)
			probclustersB <- unique(data.frame(pikB, clusB))[,1]
			nclustersB <- length(clustersB)
			if (nclustersB < 3)		
				stop("Number of clusters from frame B is less than 3. Variance cannot be computed.")		

			YcB <- matrix(0, nclustersB, c)

			for (i in 1:nclustersB){
				if (!is.null(dim(drop(piB))))
					YcB[i,] <- PEL(ysA, ysB[clusB %in% clustersB[-clustersB[i]],], piA, piB[clusB %in% clustersB[-clustersB[i]],clusB %in% clustersB[-clustersB[i]]], domainsA, domainsB[clusB %in% clustersB[-clustersB[i]]], N_A, N_B, N_ab, xsAFrameA, xsBFrameA[clusB %in% clustersB[-clustersB[i]],], xsAFrameB, xsBFrameB[clusB %in% clustersB[-clustersB[i]],], XA, XB)[[2]][1,]
				else
					YcB[i,] <- PEL(ysA, ysB[clusB %in% clustersB[-clustersB[i]],], piA, piB[clusB %in% clustersB[-clustersB[i]]], domainsA, domainsB[clusB %in% clustersB[-clustersB[i]]], N_A, N_B, N_ab, xsAFrameA, xsBFrameA[clusB %in% clustersB[-clustersB[i]],], xsAFrameB, xsBFrameB[clusB %in% clustersB[-clustersB[i]],], XA, XB)[[2]][1,]
			}

			YcBMean <- matrix(colMeans(YcB), nclustersB, c, byrow = TRUE)

			if (fcpB)
				fB <- 1 - mean(probclustersB)
			else
				fB <- 1

			vjB <- ((nclustersB - 1) / nclustersB) * fB * colSums ((YcB - YcBMean)^2)
		}
		else{

			if(sdB == "strclu"){

				strataB <- unique(strB)
				nstrataB <- length(strataB)
				nhB <- table(strB)
				YcstrataB <- matrix(0, nstrataB, c)
				nhB <- c(0, nhB)
				cnhB <- cumsum(nhB)

				for (i in 1:nstrataB){

					clustersB <- unique(clusB[strB == strataB[i]])
					probclustersB <- unique(data.frame(pikB[strB == strataB[i]], clusB[strB == strataB[i]]))[,1]
					nclustersB <- length(clustersB)
					if (nclustersB < 3)		
						stop("Number of clusters in any stratum from frame B is less than 3. Variance cannot be computed.")
					k <- 1
					YcB <- matrix(0, nclustersB, c)
					for (j in 1:nclustersB){

						if (!is.null(dim(drop(piB))))
							YcB[k,] <- PEL(ysA, ysB[clusB %in% clustersB[-clustersB[j]],], piA, piB[clusB %in% clustersB[-clustersB[j]],clusB %in% clustersB[-clustersB[j]]], domainsA, domainsB[clusB %in% clustersB[-clustersB[j]]], N_A, N_B, N_ab, xsAFrameA, xsBFrameA[clusB %in% clustersB[-clustersB[j]],], xsAFrameB, xsBFrameB[clusB %in% clustersB[-clustersB[j]],], XA, XB)[[2]][1,]
						else
							YcB[k,] <- PEL(ysA, ysB[clusB %in% clustersB[-clustersB[j]],], piA, piB[clusB %in% clustersB[-clustersB[j]]], domainsA, domainsB[clusB %in% clustersB[-clustersB[j]]], N_A, N_B, N_ab, xsAFrameA, xsBFrameA[clusB %in% clustersB[-clustersB[j]],], xsAFrameB, xsBFrameB[clusB %in% clustersB[-clustersB[j]],], XA, XB)[[2]][1,]
						k <- k + 1
					}

					YcBMean <- matrix(colMeans(YcB), nclustersB, c, byrow = TRUE)

					if (fcpB)
						fB <- 1 - mean(probclustersB)
					else
						fB <- 1
					
					YcstrataB[i,] <- (nclustersB - 1) / nclustersB * fB * colSums((YcB - YcBMean)^2)
				}
				vjB <- colSums(YcstrataB)
			}else{
				nB <- nrow(ysB)
				YcB <- matrix(0, nB, c)

				for (i in 1:nB){

					if (!is.null(dim(drop(piB))))
						YcB[i,] <- PEL(ysA, ysB[-i,], piA, piB[-i,-i], domainsA, domainsB[-i], N_A, N_B, N_ab, xsAFrameA, xsBFrameA[-i,], xsAFrameB, xsBFrameB[-i,], XA, XB)[[2]][1,]
					else
						YcB[i,] <- PEL(ysA, ysB[-i,], piA, piB[-i], domainsA, domainsB[-i], N_A, N_B, N_ab, xsAFrameA, xsBFrameA[-i,], xsAFrameB, xsBFrameB[-i,], XA, XB)[[2]][1,]
				}

				YcBMean <- matrix(colMeans(YcB), nB, c, byrow = TRUE)

				if (fcpB)
					if (!is.null(dim(drop(piB))))
						fB <- 1 - mean(pikB)
					else
						fB <- 1 - mean(piB)
				else
					fB <- 1

				vjB <- ((nB - 1) / nB) * fB * colSums((YcB - YcBMean)^2)
			}
		}
	}

	VJack_Yhat_PEL <- vjA + vjB

	results[1,] <- estimation[[2]][1,]
	results[2,] <- estimation[[2]][1,] + qnorm(1 - (1 - conf_level) / 2) * sqrt(VJack_Yhat_PEL)
	results[3,] <- estimation[[2]][1,] - qnorm(1 - (1 - conf_level) / 2) * sqrt(VJack_Yhat_PEL)
	results[4,] <- estimation[[2]][2,]
	results[5,] <- estimation[[2]][2,] + qnorm(1 - (1 - conf_level) / 2) * sqrt(1/size_estimation^2 * VJack_Yhat_PEL)
	results[6,] <- estimation[[2]][2,] - qnorm(1 - (1 - conf_level) / 2) * sqrt(1/size_estimation^2 * VJack_Yhat_PEL)

	return(results)
}