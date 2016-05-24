#' @import TeachingSampling
#' @export
#' 
#' @title
#' Sample Sizes in Two-Stage sampling Designs for Estimating Signle Proportions
#' @description 
#' This function computes a grid of possible sample sizes for estimating single proportions under two-stage sampling designs.
#' @return 
#' This function returns a grid of possible sample sizes. 
#' The first column represent the design effect,
#' the second column is the number of clusters to be selected, 
#' the third column is the number of units to be selected inside the clusters, 
#' and finally, the last column indicates the full sample size induced by this particular strategy.
#' @details
#' In two-stage (2S) sampling, the design effect is defined by
#' \deqn{DEFF = 1 + (m-1)\rho} 
#' Where \eqn{\rho} is defined as the intraclass correlation coefficient,  
#' m is the average sample size of units selected inside each cluster. 
#' The relationship of the full sample size of the two stage design (2S) with the 
#' simple random sample (SI) design is given by
#' \deqn{ n_{2S} =  n_{SI}*DEFF} 
#' @author Hugo Andres Gutierrez Rojas <hugogutierrez at usantotomas.edu.co>
#' @param N The population size.
#' @param p The value of the estimated proportion.
#' @param conf The statistical confidence. By default conf = 0.95. By default \code{conf = 0.95}.
#' @param me The maximun margin of error that can be allowed for the estimation.
#' @param M Number of clusters in the population.
#' @param by number: increment of the sequence in the grid.
#' @param rho The Intraclass Correlation Coefficient.
#' 
#' @references 
#' Gutierrez, H. A. (2009), \emph{Estrategias de muestreo: Diseno de encuestas y estimacion de parametros}. Editorial Universidad Santo Tomas
#' @seealso \code{\link{ICC}}
#' 
#' @examples
#' 
#' ss2s4p(N=100000, p=0.5, me=0.03, M=50, by=5, rho=0.01)
#' ss2s4p(N=100000, p=0.5, me=0.05, M=50, by=5, rho=0.1)
#' ss2s4p(N=100000, p=0.5, me=0.03, M=500, by=100, rho=0.2) 
#' 
#' ############################
#' # Example 2 with Lucy data #
#' ############################
#' 
#' data(Lucy)
#' attach(Lucy)
#' N <- nrow(Lucy)
#' p <- prop.table(table(SPAM))[1]
#' y <- as.double(SPAM)
#' cl <- Zone
#' 
#' rho <- ICC(y,cl)$ICC
#' M <- length(levels(Zone))
#' ss2s4p(N, 0.99, conf=0.9, me=0.03, M=5, by=1, rho=rho)

ss2s4p <- function(N, p, conf=0.95, me=0.03, M, by, rho){
  
  mseq <- seq(1, M, by=by)
  NIseq <- rep(NA, times=length(mseq))
  Deffseq <- rep(NA, times=length(mseq))
  n2seq <- rep(NA, times=length(mseq))
  
  for(i in 1: length(mseq)){
    Deffseq[i] = 1 + (mseq[i] - 1) * rho
    n2seq[i] = ss4p(N, p, DEFF=Deffseq[i], conf=conf, me=me)$n.me
    NIseq[i] <- ceiling(n2seq[i]/mseq[i])
  }
  
  result <- data.frame("Deff"=Deffseq, "NI"=NIseq, "m"=mseq, "n2s"=n2seq)
  result.adj <- result[(result$NI <= N/M),,]
  result.adj
}


