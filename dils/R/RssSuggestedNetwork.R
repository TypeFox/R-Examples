#' Suggest a network with imputed links
#'
#' A longer description of the function.  This can be perhaps
#' a paragraph, perhaps more than one.
#' 
#' @param g Object type, then description of \code{arg1}.
#' @param rss Object type, then description of \code{arg2}.
#' @param q.impute.above Object type, then description of \code{arg3}.
#' @return list
#' \item{g.imputed}{\code{igraph} contatining the original and the new links}
#' \item{g.new}{\code{igraph} containing only the new links}
#' \item{g.original}{original graph}
#' \item{q.impute.above}{quantile of RSS scores above which links should be imputed}
#' \item{frac.filled}{fraction of potential links that were actually filled with a new link}
#' @export
#' @seealso \code{\link{RelationStrengthSimilarity}}
#' @references
#' \url{http://www.haptonstahl.org/R}
#' @author Stephen R. Haptonstahl \email{srh@@haptonstahl.org}
#' @examples
#' g <- graph.atlas(128)
#' \dontrun{plot(g)}
#' 
#' suggested <- RssSuggestedNetwork(g, q.impute.above=.6)
#' \dontrun{plot(suggested$g.imputed)}
#' suggested$frac.filled
RssSuggestedNetwork <- function(g,
                                rss,
                                q.impute.above=.8) {
  # Guardians
  stopifnot(is(g, "igraph"),
            is(q.impute.above, "numeric"),
            1 == length(q.impute.above)
  )
  
  adj.original <- as.matrix(get.adjacency(g))
  
  if(missing(rss)) {
    rss <- RelationStrengthSimilarity(adj.original, 
                                      radius=2,
                                      directed=is.directed(g))
  } else {
    stopifnot(is(rss, "numeric"),
              is(rss, "matrix"),
              nrow(rss) == vcount(g),
              ncol(rss) == vcount(g))
  }
    
  # perform the function
  mode.g <- ifelse(is.directed(g), "directed", "undirected")
  
  rss.quantiles <- rss / max(rss)
  adj.threshhold <- ifelse(rss.quantiles > q.impute.above, 1, 0)
  adj.imputed <- (adj.original | adj.threshhold) + 0
  adj.new <- adj.imputed - adj.original
  
  g.imputed <- graph.adjacency(adj.imputed, mode=mode.g)
  g.new <- graph.adjacency(adj.new, mode=mode.g)
  
  frac.filled <- (graph.density(g.imputed) - graph.density(g)) / (1 - graph.density(g))
  
  # prepare and return the output
  out <- list(g.imputed=g.imputed,
              g.new=g.new,
              g.original=g,
              q.impute.above=q.impute.above,
              frac.filled=frac.filled)
  return(out)
}