##' Identify the clustered and continuous patterns of the genetic
##' variation using the PCoC, which calculates the principal
##' coordinates and the clustering of the subjects for correcting for
##' PS.
##'
##' The hidden population structure is a possible confounding effect
##' in the large-scale genome-wide association studies. Cases and
##' controls might have systematic differences because of the
##' unrecognized population structure. The PCoC procedure uses the
##' techniques from the multidimensional scaling and the clustering to
##' correct for the population stratification. The PCoC could be seen
##' as an extension of the EIGENSTRAT.
##' @title PCoC for correcting for population stratification
##' @param genoFile a txt file containing the genotypes (0, 1, 2, or
##' 9). The element of the file in Row \emph{i} and Column \emph{j}
##' represents the genotype at the \emph{i}th marker of the \emph{j}th
##' subject. 0, 1, and 2 denote the number of risk alleles, and 9
##' (default) is for the missing genotype.
##' @param outFile.txt a txt file for saving the result of this
##' function.
##' @param n.MonterCarlo the number of times for the Monter Carlo
##' procedure. The default is \code{1000}.
##' @param num.splits the number of groups into which the markers are
##' split. The default is \code{10}.
##' @param miss.val the number representing the missing data in the
##' input data. The default is \code{9}. The element 9 for the missing data
##' in the \code{genoFile} should be changed according to the value of
##' \code{miss.val}.
##' @return A list of \code{principal.coordinates} and
##' \code{cluster}. \code{principal.coordinates} is the principal
##' coordinates and \code{cluster} is the clustering of the
##' subjects. If the number of the clusters is only one,
##' \code{cluster} is omitted.
##' @author Lin Wang, Wei Zhang, and Qizhai Li.
##' @references Q Li and K Yu. Improved Correction for Population
##' Stratification in Genome-Wide Association Studies by Identifying
##' Hidden Population Structures. \emph{Genetic Epidemiology}. 2008;
##' 32(3): 215-226.
##' @references KV Mardia, JT Kent, and JM Bibby. Multivariate
##' Analysis. \emph{New York: Academic Press}. 1976.
##' @examples
##' pcocG.eg <- matrix(rbinom(4000, 2, 0.5), ncol = 40)
##' write.table(pcocG.eg, file = "pcocG.eg.txt", quote = FALSE,
##'        sep = "", row.names = FALSE, col.names = FALSE)
##' pcoc(genoFile = "pcocG.eg.txt", outFile.txt = "pcoc.result.txt",
##'        n.MonterCarlo = 50, num.splits = 10, miss.val = 9)
##' @export
pcoc <- function(genoFile="pcocG.eg.txt", outFile.txt="pcoc.result.txt", n.MonterCarlo = 1000, num.splits=10, miss.val=9)
{
    ## read genotype file
    xStr <- readLines(con=genoFile)
    num.subject <- nchar(xStr[1])
    num.marker <- length(xStr)

    ## calculate the lines of used to calculate similarity matrix each time
    if (num.splits==1)
    {
        num.lines <- num.marker
	      cum.lines <- c(0,num.marker)
    }else
    {
        num.lines <- rep(0, num.splits)
        y <- floor(num.marker/num.splits)
        num.lines[1:(num.splits-1)] <- y
        num.lines[num.splits] <- num.marker - y*(num.splits-1)
        cum.lines <- cumsum(c(0, num.lines))
    }

    ## calculate similarity matrix
    S <- matrix(data=0, nrow=num.subject, ncol=num.subject)

    for (i in 1:num.splits)
    {
        w <- (cum.lines[i]+1):cum.lines[i+1]
        x <- matrix(unlist(lapply(xStr[w], Str2Num)), nrow=num.lines[i], ncol=num.subject, byrow=TRUE)
	      is.na(x[x==miss.val])<-T

	      S <- S + SimilarityMatrix(x, num.lines[i])

	      rm(w,x)
    }

    S <- S/num.marker

    X <- MDS(S, num.subject, TopK=NULL, SignEigenPoint=2.0234)
    k <- FindCNumRandom(X, num.subject, kG=10, n.MonterCarlo)

    if (k>1)
    {
        F <- cluster::clara(X, k, samples=20, medoids.x=FALSE)$clustering
      	V <- list(principal.coordinates=X, cluster=factor(F))
      	V1 <- data.frame(X, factor(F))
    }
    else
    {
        V <- list(principal.coordinates=X)
        V1 <- data.frame(X)
    }

    write(t(V1), file=outFile.txt, ncolumns=ncol(V1), sep="\t")

    V
}

