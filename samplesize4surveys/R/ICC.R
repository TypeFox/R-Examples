#' @import TeachingSampling
#' @export
#' 
#' @title
#' Intraclass Correlation Coefficient
#' @description 
#' This function computes the intraclass correlation coefficient.
#' @return 
#' The total sum of squares (TSS), the between sum of squqres (BSS), the within sum of squares (WSS) and the intraclass correlation coefficient.
#' @details
#' The intraclass correlation coefficient is defined as: 
#' \deqn{\rho = 1- \frac{m}{m-1} \frac{WSS}{TSS}} 
#' Where \eqn{m} is the average sample sie of units selected inside each sampled cluster. 
#' @author Hugo Andres Gutierrez Rojas <hugogutierrez at usantotomas.edu.co>
#' @param y The variable of interest. 
#' @param cl The variable indicating the membership of each element to a specific cluster. 
#' 
#' @references 
#' Gutierrez, H. A. (2009), \emph{Estrategias de muestreo: Diseno de encuestas y estimacion de parametros}. Editorial Universidad Santo Tomas
#' @seealso \code{\link{ss4p}}
#' 
#' @examples 
#' 
#' ##########################################
#' # Almost same mean in each cluster       #
#' #                                        #
#' # - Heterogeneity within clusters        #
#' # - Homogeinity between clusters         #
#' ##########################################
#' 
#' # Population size
#' N <- 100000
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
#' table(cl)
#' tapply(y, cl, FUN=mean)
#' boxplot(y~cl)
#' rho = ICC(y,cl)$ICC
#' rho
#' 
#' 
#' ##########################################
#' # Very different means per cluster       #
#' #                                        #
#' # - Heterogeneity between clusters       #
#' # - Homogeinity within clusters          #
#' ##########################################
#' 
#' # Population size
#' N <- 100000
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
#' table(cl)
#' tapply(y, cl, FUN=mean)
#' boxplot(y~cl)
#' rho = ICC(y,cl)$ICC
#' rho
#' 
#' ############################
#' # Example 1 with Lucy data #
#' ############################
#' 
#' data(Lucy)
#' attach(Lucy)
#' N <- nrow(Lucy)
#' y <- Income
#' cl <- Zone
#' ICC(y,cl)
#' 
#' ############################
#' # Example 2 with Lucy data #
#' ############################
#' 
#' data(Lucy)
#' attach(Lucy)
#' N <- nrow(Lucy)
#' y <- as.double(SPAM)
#' cl <- Zone
#' ICC(y,cl)


ICC <- function(y, cl){
  cl <- as.factor(cl)
  N <- length(y)
  Ni <- table(cl)
  t1 <- (tapply(y,cl,mean)-mean(y))^2
  TSS <- (N-1)*var(y)
  BSS <- sum(Ni*t1)
  WSS <- TSS - BSS
  ICC <- 1-(sum(Ni)/sum(Ni-1))*(WSS/TSS)
  res <- list(TSS=TSS, BSS=BSS, WSS=WSS, ICC=ICC)
  res
}

