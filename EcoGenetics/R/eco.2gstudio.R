#' Converting a diploid ecogen genetic data frame into a gstudio object
#' 
#' @description This function converts the genetic 
#' data of an ecogen object in a gstudio data frame. 
#' @param eco Object of class "ecogen".
#' @param ... Further arguments passed to \code{\link[adegenet]{df2genind}}.
#' @param type The type of data passed to gstudio locus. 
#' Default is "separated" (data as microsatellites or individual haplotypes);
#' "aflp" for presence - absence data. 
#' 
#' @examples
#' 
#' \dontrun{
#' 
#' data(eco.test)
#' gsteco <- eco.2gstudio(eco, "separated")
#' gsteco
#' 
#' }
#' 
#' @author Leandro Roser \email{leandroroser@@ege.fcen.uba.ar}
#' @export

setGeneric("eco.2gstudio", 
           function(eco, type = "separated", ...) {
             
             
             if(type == "separated") {
               dat <- int.df2genind(eco@G)
               dat <- int.genind2df(dat, sep = ":")
               for(i in 1:ncol(dat)) {  
                 dat[, i] = gstudio::locus(dat[, i], type = "separated")
               }
             } else {
               dat<-eco@G
               for(i in 1:ncol(dat)) {
                 dat[, i] = gstudio::locus(dat[, i], type = type)
               }
             }
             
             dat
           })
