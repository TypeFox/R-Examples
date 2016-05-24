#' An R package for running the GUESS program.  GUESS is a computationally optimised C++ implementation of a fully Bayesian variable selection approach that can analyse, in a genome-wide context, single and multiple responses in an integrated way. The program uses packages from the GNU Scientific Library (GSL) and offers the possibility to re-route computationally intensive linear algebra operations towards the Graphical Processing Unit (GPU) 
#' through the use of proprietary CULA-dense library.
#' The multi-SNP model of GUESS typically seeks for the best combinations of SNPs to predict the (possibly multivariate) outcome of interest.
#' In its current implementation, using its GPU capacities, GUESS is able to handle hundreds of thousands of predictors, which enables genome-wide 
#' sized datasets to be analysed. However, the use of GPU-based numerical libraries implies extensive data transfer between the memory/CPU and the GPU, which, in turn, can be computationally expensive. As a consequence, for smaller datasets (as the example provided in the package) for which the matrix operations are not rate-limiting, the CPU version of GUESS may be more efficient. 
#'  Hence, to ensure both an optimal use of the algorithm, and to enable running GUESS on non-CULA compatible systems, the call to GPU-based calculations within GUESS can easily be switched off with the argument \cite{CUDA} of the function \cite{\link{R2GUESS}}.
#'  A documentation of the C++ code is presenting at \url{http://www.bgx.org.uk/software/GUESS_Doc_short.pdf}. 
#' @name R2GUESS-package
#' @aliases R2GUESS-package
#' @docType package
#' @title Sparse Bayesian variable selection method for multiple correlated outcomes in a regression context.
#' @author Benoit Liquet, \email{benoit.liquet@@isped.u-bordeaux2.fr} , Marc Chadeau \email{m.chadeau@@imperial.ac.uk}, Leonardo Botollo \email{l.bottolo@@imperial.ac.uk}, Gianluca Campanelle \email{g.campanella11@@imperial.ac.uk}
#' @references
#' Bottolo L and Richardson S (2010). Evolutionary Stochastic Search for Bayesian model exploration. Bayesian Analysis 5(3), 583-618.
#' 
#' Petretto E, Bottolo L, Langley SR, Heinig M, McDermott-Roe C, Sarwar R, Pravenec M, Hubner N, Aitman TJ, Cook SA and Richardson S (2010).
#' New insights into the genetic control of gene expression using a Bayesian multi-tissue approach.
#' PLoS Comput. Biol., 6(4), e1000737.
#' @keywords package
#' @seealso \code{\link{R2GUESS}}, \code{\link{as.ESS.object}}, \code{\link{plotMPPI}}, \code{\link{plot.ESS}}
NULL
#' Data set relative to the regulation genetic of the experssion level of the gene HOPX from the Rats. 
#' 
#' This data set contains the gene expresion of the gene Hopx from 4 tissues of rats (Adrenal gland, Fat, Heart, Kidney).  
#' @name data.Y.Hopx
#' @docType data
#' @format A dataframe with 29 rows and 4 columns (ADR, FAT, Heart, Kidney)
NULL
#' Data set relative to the regulation genetic of the Rats. 
#' 
#' This data set contains the genetic variation (SNPs) of 29 rats.  
#' @name data.X
#' @docType data
#' @format A dataframe with 29 rows and 770 columns
NULL
#' MAP file relative to the genotype performed on the rats . 
#' 
#' This data set contains the information on the SNPs: SNP name (\code{SNPName}), Chromosome (\code{Chr}), start position (\code{Posn}), end position (\code{End}) 
#' @name MAP.file
#' @docType data
#' @format A data frame with 770 rows and 4 columns
NULL

