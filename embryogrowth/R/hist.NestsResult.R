#' hist.NestsResult shows the histogram of temperatures with set of nests
#' @title Show the histogram of temperatures with set of nests
#' @author Marc Girondot \email{marc.girondot@@u-psud.fr}
#' @return A list with an histogram object with information on histogram or 
#' NULL if no series was selected and the complete set of temperatures used.
#' @param x Results obtained after searchR
#' @param series Series to be used, logical (TRUE ou FALSE), numbers or names. If "all", all series are used.
#' @param ... Parameters used by hist function (example main="Title")
#' @description Show the histogram of temperatures with set of nests
#' hist(data)
#' @examples
#' \dontrun{
#' library(embryogrowth)
#' data(resultNest_4p)
#' h <- hist(resultNest_4p, series=c(1:5))
#' }
#' @method hist NestsResult
#' @export

hist.NestsResult <- function(x, series="all", ...) {

# Je prends seulement les donnees de temperature et j'envoie a hist.Nests
# Ce sera plus simple pour faire les mises a jour - 30/7/2012

# j'ai un objet de resultat
# je prends les donnees
nids <- x$data
class(nids) <- "Nests"

L <- modifyList(list(x=nids), list(...))
L <- modifyList(L, list(series=series))

a <- do.call(hist, L)

return(invisible(a))

}
