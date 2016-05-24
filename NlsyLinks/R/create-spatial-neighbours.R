#' @name CreateSpatialNeighbours
#' @export
#' 
#' @title Distances between related family members, formated for spatial analysis.
#' @description This helper function formats the LinksPair datasets so it can be used in some types of spatial analyses. The \pkg{spdep} (Spatial Dependence) uses a sparse matrix (actually a \code{\link{data.frame}}) to represent neigbours.
#' @usage  CreateSpatialNeighbours(linksPairsDoubleEntered)
## CreateSpatialNeighbours79Gen2()
#' @param linksPairsDoubleEntered A \code{data.frame} containing the links, preferably created by a function like 
#' 
#' \code{\link{CreatePairLinksDoubleEntered}}.
#' 
#' @details There is one row per unique pair of subjects, \emph{respecting order}.  This has twice as many rows as \code{\link{Links79Pair}} and \code{\link{Links79PairExpanded}} (which have one row per unique pair of subjects, \emph{irrespective of order}).
#' 
#' \code{CreateSpatialNeighbours} accepts any paired relationships in a \code{data.frame}, as long as it contains the columns \code{SubjectTag_S1}, \code{SubjectTag_S2}, and \code{R}.  See \code{\link{Links79Pair}} for more details about these columns.
#' 
#' @return An S3 \code{spatial.neighbours} object to work with functions in the \pkg{spdep} package. 
#' 
#' \code{SubjectTag_S1} is renamed `\code{from}'. 
#' 
#' \code{SubjectTag_S2} is renamed `\code{to}'. 
#' 
#' \code{R} is renamed `\code{weight}'.
#' 
#' The attribute \code{region.id} specifies each unique SubjectTag.
#' 
#' The attribue \code{n} specifies the number of unique subjects.
#' 
#' @references 
#' 
#' Bard, D.E., Beasley, W.H., Meredith, K., & Rodgers, J.L. (2012). \href{http://link.springer.com/article/10.1007/s10519-012-9566-6}{\emph{Biometric Analysis of Complex NLSY Pedigrees: Introducing a Conditional Autoregressive Biometric (CARB) Mixed Model}}. Behavior Genetics Association 42nd Annual Meeting. [\href{https://r-forge.r-project.org/forum/forum.php?thread_id=4761&forum_id=4266&group_id=1330}{Slides}]
#' 
#' Bivand, R., Pebesma, E., & Gomez-Rubio, V. (2013). \href{http://link.springer.com/book/10.1007/978-1-4614-7618-4}{\emph{Applied Spatial Data Analysis with R.}} New York: Springer. (Especially Chapter 9.)
#' 
#' Banerjee, S., Carlin, B.P., & Gelfand, A.E. (2004). \href{http://books.google.com/books/about/Hierarchical_Modeling_and_Analysis_for_S.html?id=YqpZKTp-Wh0C}{\emph{Hierarchical Modeling and Analysis for Spatial Data}}. Boca Raton: CRC Press.
#' 
#' Lawson, A.B (2013). \href{http://books.google.com/books?id=g7RJEZb1umwC}{\emph{Bayesian Disease Mapping: Hierarchical Modeling in Spatial Epidemiology, Second Edition}}. Boca Raton: CRC Press.
#' 
#' The \pkg{spdep} package documentation: \href{http://cran.r-project.org/package=spdep}{spdep: Spatial dependence: weighting schemes, statistics and models}.
#' 
#' @author Will Beasley and  David Bard
#' @note Notice the British variant of 'neighbo\emph{u}rs' is used, to be consistent with the \code{spatial.neighbour} class in the  \href{http://cran.r-project.org/package=spdep}{\code{spdep}} package.
#' @examples
#' dsLinksAll <- Links79Pair
#' dsLinksGen1Housemates <- dsLinksAll[dsLinksAll$RelationshipPath=="Gen1Housemates", ]
#' dsLinksGen2Siblings <- dsLinksAll[dsLinksAll$RelationshipPath=="Gen2Siblings", ]
#' 
#' spGen1 <- CreateSpatialNeighbours(dsLinksGen1Housemates)
#' spGen2 <- CreateSpatialNeighbours(dsLinksGen2Siblings)
#' 
#' head(spGen2)
#' #Returns: 
#' #   from  to weight
#' #3   201 202   0.50
#' #6   301 302   0.50
#' #7   301 303   0.50
#' #9   302 303   0.50
#' #24  401 403   0.25
#' #28  801 802   0.50
#' 
#' table(spGen2$weight)
#' #Returns:
#' #0.25 0.375   0.5  0.75     1 
#' #3442   610  6997    12    27 
#' @keywords spatial analysis
#' 
CreateSpatialNeighbours <- function( linksPairsDoubleEntered )  {
  ValidatePairLinks(linksPairsDoubleEntered)
  
  ds <- base::subset(linksPairsDoubleEntered, select=c("SubjectTag_S1", "SubjectTag_S2", "R"))
  base::colnames(ds)[base::colnames(ds) == "SubjectTag_S1"] <- "from"
  base::colnames(ds)[base::colnames(ds) == "SubjectTag_S2"] <- "to"
  base::colnames(ds)[base::colnames(ds) == "R"] <- "weight"
  # summary(ds)
  
  base::class(ds) <- c("spatial.neighbour", base::class(ds))
  base::attr(ds, "region.id") <- base::unique(ds$from)
  base::attr(ds, "n") <- base::length(base::unique(ds$from)) 
  return( ds )
}
