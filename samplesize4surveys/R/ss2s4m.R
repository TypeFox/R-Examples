#' @import TeachingSampling
#' @export
#' 
#' @title
#' Sample Sizes in Two-Stage sampling Designs for Estimating Signle Means
#' @description 
#' This function computes a grid of possible sample sizes for estimating single means under two-stage sampling designs.
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
#' @param mu The value of the estimated mean of a variable of interest.
#' @param sigma The value of the estimated standard deviation of a variable of interest.
#' @param conf The statistical confidence. By default conf = 0.95. By default \code{conf = 0.95}.
#' @param rme The maximun relative margin of error that can be allowed for the estimation.
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
#' ss2s4m(N=100000, mu=10, sigma=2, conf=0.95, rme=0.03, M=50, by=5, rho=0.01)
#' ss2s4m(N=100000, mu=10, sigma=2, conf=0.95, rme=0.03, M=50, by=5, rho=0.1)
#' ss2s4m(N=100000, mu=10, sigma=2, conf=0.95, rme=0.03, M=50, by=5, rho=0.2)
#' ss2s4m(N=100000, mu=10, sigma=2, conf=0.95, rme=0.05, M=50, by=5, rho=0.3)
#' 
#' ##########################################
#' # Almost same mean in each cluster       #
#' #                                        #
#' # - Heterogeneity within clusters        #
#' # - Homogeinity between clusters         #
#' #                                        #
#' #  Decision rule:                        #
#' #    * Select a lot of units per cluster #
#' #    * Select a few of clusters          #
#' ##########################################
#' 
#' # Population size
#' N <- 1000000
#' # Number of clusters in the population
#' NI <- 1000
#' # Number of elements per cluster
#' N/NI
#' 
#' # The variable of interest
#' y <- c(1:N)
#' # The clustering factor
#' cl <- rep(1:NI, length.out=N)
#' 
#' rho = ICC(y,cl)$ICC
#' rho
#' 
#' ss2s4m(N, mu=mean(y), sigma=sd(y), conf=0.95, rme=0.03, M=N/NI, by=10, rho)
#' 
#' 
#' ##########################################
#' # Very different means per cluster       #
#' #                                        #
#' # - Heterogeneity between clusters       #
#' # - Homogeinity within clusters          #
#' #                                        #
#' #  Decision rule:                        #
#' #    * Select a few of units per cluster #
#' #    * Select a lot of clusters          #
#' ##########################################
#' 
#' # Population size
#' N <- 1000000
#' # Number of clusters in the population
#' NI <- 1000
#' # Number of elements per cluster
#' N/NI
#' 
#' # The variable of interest
#' y <- c(1:N)
#' # The clustering factor
#' cl <- kronecker(c(1:NI),rep(1,N/NI))
#' 
#' rho = ICC(y,cl)$ICC
#' rho
#' 
#' ss2s4m(N, mu=mean(y), sigma=sd(y), conf=0.95, rme=0.03, M=N/NI, by=10, rho)

ss2s4m <- function(N, mu, sigma, conf=0.95, rme=0.03, M, by, rho){

  mseq <- seq(1,M, by=by)
  NIseq <- rep(NA, times=length(mseq))
  Deffseq <- rep(NA, times=length(mseq))
  n2seq <- rep(NA, times=length(mseq))
  
  for(i in 1: length(mseq)){
    Deffseq[i] = 1 + (mseq[i] - 1) * rho
    n2seq[i] = ss4m(N, mu=mu, sigma=sigma, DEFF=Deffseq[i], conf=conf, rme=rme)$n.rme
    NIseq[i] <- ceiling(n2seq[i]/mseq[i])
  }

  result <- data.frame("Deff"=Deffseq, "NI"=NIseq, "m"=mseq, "n2s"=n2seq)
  result.adj <- result[(result$NI <= N/M),]
  result.adj
}


