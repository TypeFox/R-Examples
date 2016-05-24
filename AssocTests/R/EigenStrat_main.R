##' Find the eigenvectors of the similarity matrix among the subjects
##' used for correcting for population stratification in the
##' population-based genetic association studies.
##'
##' Suppose that a total of \emph{n} cases and controls are randomly
##' enrolled in the source population and a panel of \emph{m}
##' single-nucleotide polymorphisms are genotyped. The genotype at a
##' marker locus is coded as 0, 1, or 2, with the value corresponding
##' to the copy number of risk alleles. All the genotypes are given in
##' the form of a \emph{m*n} matrix, in which the element in the
##' \emph{i}th row and the \emph{j}th column represents the genotype
##' of the \emph{j}th subject at the \emph{i}th marker. This function
##' calculates the top ten eigenvectors or the eigenvectors with
##' significant eigenvalues of the similarity matrix among the
##' subjects to infer the potential population structure. See also
##' \link{tw}.
##' @title EIGENSTRAT for correcting for population stratification
##' @param genoFile a txt file containing the genotypes (0, 1, 2, or
##' 9). The element of the file in Row \emph{i} and Column \emph{j}
##' represents the genotype at the \emph{i}th marker of the \emph{j}th
##' subject. 0, 1, and 2 denote the number of risk alleles, and 9
##' (default) is for the missing genotype.
##' @param outFile.Robj the name of an R object for saving the list of
##' the results which is the same as the return value of this
##' function. The default is "\code{out.list}".
##' @param outFile.txt a txt file for saving the eigenvectors
##' corresponding to the top 10 significant eigenvalues.
##' @param rm.marker.index a numeric vector for the indices of the
##' removed markers. The default is \code{NULL}.
##' @param rm.subject.index a numeric vector for the indices of the
##' removed subjects. The default is \code{NULL}.
##' @param miss.val the number representing the missing data in the
##' input data. The default is \code{9}. The element 9 for the missing data
##' in the \code{genoFile} should be changed according to the value of
##' \code{miss.val}.
##' @param num.splits the number of groups into which the markers are
##' split. The default is \code{10}.
##' @param topK the number of eigenvectors to return. If \code{NULL}, it is
##' calculated by the Tracy-Wisdom test. The default is \code{NULL}.
##' @param signt.eigen.level a numeric value which is the significance
##' level of the Tracy-Wisdom test. It should be \code{0.05}, \code{0.01}, \code{0.005}, or
##' \code{0.001}. The default is \code{0.01}.
##' @param signal.outlier logical. If \code{TRUE}, delete the outliers of the
##' subjects; otherwise, do not search for the outliers. The default
##' is \code{FALSE}.
##' @param iter.outlier a numeric value that is the iteration time for
##' finding the outliers of the subjects. The default is \code{5}.
##' @param sigma.thresh a numeric value that is the lower limit for
##' eliminating the outliers. The default is \code{6}.
##' @return \code{eigenstrat} returns a list, which contains the following components:
##' \tabular{llll}{
##' \code{num.markers} \tab \tab \tab the number of the markers excluding the removed markers.\cr
##' \code{num.subjects} \tab \tab \tab the number of the subjects excluding the outliers.\cr
##' \code{rm.marker.index} \tab \tab \tab the indices of the removed markers.\cr
##' \code{rm.subject.index} \tab \tab \tab the indices of the removed subjects.\cr
##' \code{TW.level} \tab \tab \tab the significance level of the Tracy-Wisdom test.\cr
##' \code{signal.outlier} \tab \tab \tab dealing with the outliers in the subjects or not.\cr
##' \code{iter.outlier} \tab \tab \tab the iteration time for finding the outliers.\cr
##' \code{sigma.thresh} \tab \tab \tab the lower limit for eliminating the outliers.\cr
##' \code{num.outliers} \tab \tab \tab the number of the outliers.\cr
##' \code{outliers.index} \tab \tab \tab the indices of the outliers.\cr
##' \code{num.used.subjects} \tab \tab \tab the number of the used subjects.\cr
##' \code{used.subjects.index} \tab \tab \tab the indices of the used subjects.\cr
##' \code{similarity.matrix} \tab \tab \tab the similarity matrix among the subjects.\cr
##' \code{eigenvalues} \tab \tab \tab the eigenvalues of the similarity matrix.\cr
##' \code{eigenvectors} \tab \tab \tab the eigenvectors corresponding to the eigenvalues.\cr
##' \code{topK} \tab \tab \tab the number of the significant eigenvalues.\cr
##' \code{TW.stat} \tab \tab \tab the observed values of the Tracy-Wisdom statistics.\cr
##' \code{topK.eigenvalues} \tab \tab \tab the significant eigenvalues.\cr
##' \code{topK.eigenvectors} \tab \tab \tab the eigenvectors corresponding to the significant eigenvalues.\cr
##' \code{runtime} \tab \tab \tab the running time of this function.
##' }
##' @author Lin Wang, Wei Zhang, and Qizhai Li.
##' @references AL Price, NJ Patterson, RM Plenge, ME Weinblatt, NA
##' Shadick, and D Reich. Principal Components Analysis Corrects for
##' Stratification in Genome-Wide Association Studies. \emph{Nature
##' Genetics}. 2006; 38(8): 904-909.
##' @references N Patterson, AL Price, and D Reich. Population
##' Structure and Eigenanalysis. \emph{PloS Genetics}. 2006; 2(12):
##' 2074-2093.
##' @references CA Tracy and H Widom. Level-Spacing Distributions and
##' the Airy Kernel. \emph{Communications in Mathematical
##' Physics}. 1994; 159(1): 151-174.
##' @examples
##' eigenstratG.eg <- matrix(rbinom(3000, 2, 0.5), ncol = 30)
##' write.table(eigenstratG.eg, file = "eigenstratG.eg.txt", quote = FALSE,
##'             sep = "", row.names = FALSE, col.names = FALSE)
##' eigenstrat(genoFile = "eigenstratG.eg.txt", outFile.Robj = "eigenstrat.result.list",
##'              outFile.txt = "eigenstrat.result.txt", rm.marker.index = NULL,
##'              rm.subject.index = NULL, miss.val = 9, num.splits = 10,
##'              topK = NULL, signt.eigen.level = 0.01, signal.outlier = FALSE,
##'              iter.outlier = 5, sigma.thresh = 6)
##' @export
eigenstrat <- function(genoFile = "eigenstratG.eg.txt", outFile.Robj = "out.list", outFile.txt = "out.txt", rm.marker.index = NULL, rm.subject.index = NULL,
                              miss.val = 9, num.splits = 10, topK = NULL, signt.eigen.level = 0.01, signal.outlier = FALSE, iter.outlier = 5, sigma.thresh = 6)
{
    timeX <- proc.time()

    SignEigenPoint <- 0
    if (signt.eigen.level == 0.05)
    {
        SignEigenPoint <- 0.9793
    }else if (signt.eigen.level == 0.01)
    {
        SignEigenPoint <- 2.0234
    }else if (signt.eigen.level == 0.005)
    {
        SignEigenPoint <- 2.4224
    }else if (signt.eigen.level == 0.001)
    {
        SignEigenPoint <- 3.2724
    }else
    {
        stop("signt.eigen.level must belong to the set { 0.05, 0.01, 0.005, 0.001 }\n")
    }

    # read file
    #xStr <- scan(file=genoFile, what='character')
    xStr <- readLines(con=genoFile)
    num.subjects <- nchar(xStr[1])
    num.original <- num.subjects
    if (is.null(rm.marker.index))
    {
        num.markers <- length(xStr)
    }else
    {
        xStr <- xStr[-rm.marker.index]
        num.markers <- length(xStr)
    }

    # calculate the lines of using scan each time
    if (num.splits==1)
    {
        each.lines <- num.markers
      	used.lines <- c(0,num.markers)
    }else
    {
        each.lines <- rep(0, num.splits)
        y <- floor(num.markers/num.splits)
        each.lines[1:(num.splits-1)] <- y
        each.lines[num.splits] <- num.markers - y*(num.splits-1)
        used.lines <- cumsum(c(0, each.lines))
    }

    # if need to eliminate subjects,
    if (!is.null(rm.subject.index))
    {
        num.subjects <- num.subjects-length(rm.subject.index)
    }
    #print(paste("num.subjects=",num.subjects,sep=""))

    # if need to eliminate outliers, generate outlier vector
    if (signal.outlier)
    {
        outlier <- rep(0, num.subjects)
    }

    indicator <- TRUE
    iter <- 0
    Q <- 1:num.subjects # deposit the index of individuals

    while (indicator)
    {
        if (signal.outlier) # eliminate outliers
      	{
	          addr <- which(outlier==1)
	          num.outlier <- length(addr)
	          n <- num.subjects-num.outlier
      	}else
	      {
	           n <- num.subjects
	      }

	      xM <- matrix(data=0, nrow=n, ncol=n)

	      for (i in 1:num.splits)
        {
	          w <- (used.lines[i]+1):used.lines[i+1]
      	    x <- matrix(unlist(lapply(xStr[w], Str2Num)), nrow=each.lines[i], ncol=num.original, byrow=TRUE)
            is.na(x[x==miss.val]) <- TRUE

	          # remove rm.subject.index
	          if (!is.null(rm.subject.index))
	          {
	              x <- x[,-rm.subject.index]
	          }
	          #print(x)
            # remove outlier
    	      if (signal.outlier)
	          {
	              if (num.outlier>0)
		            {
		                x <- x[,-addr] # eliminate outliers
		            }
	          }
	          x <- x/2
      	    z <- apply(x, 1, ModifyNormalization) # modify normalization by rows
	          rm(w, x)

            xM <- xM + (z %*% t(z))
	          rm(z)
        }

        xM <- xM/num.markers

	      #print(dim(xM))

	      # eigen analysis
        eig.list <- eigen(xM)
	      eVal <- eig.list$values
	      eVec <- eig.list$vectors

      	# eigen test
      	TW.stat <- NULL
	      if (is.null(topK))
	      {
	          a <- n - 1
	          wet <- tw(eVal[1:a], a, SignEigenPoint)
	          topK <- wet$SigntEigenL
	          TW.stat <- wet$statistic
	          if (topK==0)
            {
                topK <- 2
            }
	      }
	      else
	      {
	          a <- n - 1
	          wet <- tw(eVal[1:a], a, SignEigenPoint)
	          TW.stat <- wet$statistic
        }

        # stop the iteration
	      if (iter == iter.outlier)
	      {
	          indicator <- FALSE
	          break
	      }

      	# find outliers
      	#print("It is time to remove outliers")
      	if (signal.outlier)
	      {
	          norm.eVec <- abs(scale(eVec[,(1:topK)], center = TRUE, scale = TRUE)) # normalization as ususal
	          out.dex <- unique(which(norm.eVec >= sigma.thresh) %% n)
	          #out.dex <- unique(which(abs(eVec[,1:5]) >= 0.2) %% n)
            out.dex[out.dex==0] <- n

            if (length(out.dex) > 0)
	          {
	              iter <- iter + 1
	              if (num.outlier > 0)
	              {
	                  q1 <- Q[-addr]
		                outlier[q1[out.dex]] <- 1
	              }else
	              {
	                  outlier[out.dex] <- 1
	              }
	          }else
	          {
	              indicator <- FALSE
	              break
	          }
        }else # don't eliminate outliers
	      {
	          indicator <- FALSE
	          break
	      }

	      #print("------------------------------------------------------------------------")
	      #print(list(outlier=which(outlier==1),indicator=indicator,iter=iter,iter.outlier=iter.outlier))

    }

    #print("Eigen similarity finished")
    # output
    if (signal.outlier)
    {
        if (num.outlier < 1)
        {
	          outQ <- Q
	          outAddr <- NULL
        }else
        {
            outQ <- Q[-addr]
	          outAddr <- addr
        }
    }else
    {
        outQ <- Q
	      outAddr <- NULL
      	num.outlier <- NULL
    }

    # compute the executed time
    timeY <- proc.time() - timeX

    #print(paste("run time=", timeY[1], sep=""))

    # return a list
    res.list <- list(num.markers=num.markers,            num.subjects=num.subjects,
                     rm.marker.index=rm.marker.index,    rm.subject.index=rm.subject.index,
		                 TW.level=signt.eigen.level,         signal.outlier=signal.outlier,
		                 iter.outlier=iter.outlier,          sigma.thresh=sigma.thresh,
		                 num.outliers=num.outlier,           outliers.index=outAddr,
		                 num.used.subjects=n,                used.subjects.index=outQ,
		                 similarity.matrix=xM,     	         eigenvalues=eVal,
		                 eigenvectors=eVec,                  topK=topK,
                     TW.stat=TW.stat,
		                 topK.eigenvalues=eVal[1:topK],      topK.eigenvectors=eVec[,(1:topK)],
		                 runtime=timeY[3]
		                 )

    if (!is.null(outFile.Robj))
    {
        save(res.list, file=outFile.Robj)
    }
    if (!is.null(outFile.txt))
    {
         	write.table(eVec[,1:10], file=outFile.txt, row.names=F, col.names=F, sep="\t")
    }

    res.list
}

