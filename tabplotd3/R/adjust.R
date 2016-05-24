#' adjustment of tableplot object, to be implemented in tabplot itself
#' @keywords internal 
adjust <- function(tp){
  atp <- list()
  varnms <- sapply(tp$columns, function(i) i$name)
  vars <- lapply( tp$columns
                  , function(i) { 
                       if (i$isnumeric){
                         return(list( mean = i$mean
                                 , compl = i$compl/100
                               ))
                       }
                      list( freq = i$freq, palet=i$palet, categories=i$categories)                
                    }
                )
  names(vars) <- varnms
  atp$vars <- vars
#   atp$dataset <- tp$dataset
  atp$sortCol = names(tp$columns)[tp$sortCol]
  atp$nBins <- tp$nBins
  atp
}