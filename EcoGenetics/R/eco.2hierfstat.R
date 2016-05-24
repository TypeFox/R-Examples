#' Converting an ecogen genetic data frame into a hierfstat data frame
#' 
#' @description This function converts the genetic 
#' data of an ecogen object in a hierfstat data frame. 
#' @param eco Object of class "ecogen".
#' @param pop The name of the S slot column with the groups 
#' for the hierfstat data frame.
#' 
#' @examples
#' 
#' \dontrun{
#' 
#' data(eco.test)
#' hiereco <- eco.2hierfstat(eco, "pop")
#' require("hierfstat")
#' basic.stats(hiereco)
#' 
#' }
#' 
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @export


setGeneric("eco.2hierfstat", 
           function(eco, pop = NULL) {
             
             u <- eco@G
             
             grupo <- eco@S
             
             if(is.null(pop))
             {
               factord <- as.data.frame(rep(1, nrow(u)))
               cnom <- "pop"
               rnom <- rownames(eco@G)
               Gord <- u
             } else {
               
               pop <- match(pop, colnames(eco@S), nomatch = 0)
               pop <- pop[pop != 0]
               if(length(pop) == 0) {
                 stop("incorrect factor name")
               }
               orden <- order(eco@S[, pop])
               Gord <- u[orden,]
               factord <- eco@S[orden, pop]
               factord <- as.numeric(factord)
               cnom <- colnames(eco@S[pop])
               rnom <- rownames(eco@G)[orden]
             }
             
             datahier <- data.frame(factord, Gord)
             colnames(datahier)[1] <- cnom
             rownames(datahier) <- rnom
             datahier
             
             #class control
             clases <- character()
             j <- 1
             for(i in 2:ncol(datahier)) {
               clases[j] <- class(datahier[, i])
               j <- j + 1
             }
             if(any(clases != "numeric" | clases != "integer")) {
               datahier <- as.matrix(datahier)
               colhier <- ncol(datahier)
               rowhier <- nrow(datahier)
               datahier <- matrix(as.numeric(datahier), ncol = colhier, nrow= rowhier)
               datahier <- as.data.frame(datahier)
               datahier[, 1] <- as.factor(datahier[, 1])
             }
             
               rownames(datahier) <- rownames(eco@G)
               colnames(datahier)[1] <- "Pop"
               colnames(datahier)[-1] <- colnames(eco@G)
             
               datahier
                 
           })
