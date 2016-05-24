#' Simulate an ordinal relational matrix
#' 
#' Simulates an ordinal relational matrix having a particular marginal
#' distribution
#' 
#' 
#' @usage simY_ord(EZ, rho, Y)
#' @param EZ square matrix giving the expected value of the latent Z matrix
#' @param rho scalar giving the within-dyad correlation
#' @param Y ordinal relational data matrix
#' @return a square matrix
#' @author Peter Hoff
#' @export simY_ord
simY_ord<-
function (EZ, rho,Y) 
{ 
    
    uY<-sort(unique(c(Y))) 
    FY<-table(c(Y)) ; FY<-FY/sum(FY) ; FY<-cumsum(FY) 
    ZS <- simZ(EZ, rho)
    diag(ZS) <- NA
    qZ <- quantile(ZS, FY[-length(FY)], na.rm = TRUE)
    YS <- ZS * 0 + max(uY)
    for (w in 1:(length(uY) - 1)) {
        YS[(YS == max(uY)) & ZS <= qZ[w]] <- uY[w]
    }
    YS
}



