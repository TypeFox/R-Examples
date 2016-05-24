#' Convert an object of class prevR into a data.frame.
#' 
#' This function merges the slots \code{clusters} et \code{rings} of 
#' a object of class \code{\link[=prevR-class]{prevR}}.
#' 
#' @param x object of class \code{\link[=prevR-class]{prevR}}.
#' @param ... not used, for compatibility with the generic method \code{\link[base]{as.data.frame}}.
#' @param N integer or list of integers setting elements of \code{rings} to extract.
#' @param R integer or list of integers setting elements of \code{rings} to extract.
#' @param clusters.only return only the slot \code{clusters} of \code{x}?
#' 
#' @return 
#' If \code{clusters.only = TRUE}, the function will return only the slot \code{clusters} of \code{x}.
#' 
#' Otherwise, slots \code{clusters} and \code{rings} of \code{x} will be merged in a unique data frame. 
#' The columns of \code{rings} will be renamed adding a suffix like \emph{.N300.RInf}.
#' 
#' \code{N} and \code{R} define the elements of \code{rings} to extract. If not specified (\code{NULL}), 
#' all the elements of \code{rings} will be included.
#' 
#' @seealso \code{\link[base]{as.data.frame}}\{\pkg{base}\}, \code{\link[prevR]{prevR-class}}.
#'
#' @examples 
#' str(fdhs)
#' str(as.data.frame(fdhs))
#' \dontrun{
#'  r.fdhs <- rings(fdhs, N=c(100,200,300))
#'  str(r.fdhs)
#'  str(as.data.frame(r.fdhs, clusters.only=TRUE))
#'  str(as.data.frame(r.fdhs))
#'  str(as.data.frame(r.fdhs, N=300))
#' }
#' 
#' @export
#' @aliases as.data.frame
#' @keywords manip

as.data.frame.prevR = function(x,..., N=NULL, R=NULL, clusters.only=FALSE){
  ###############################################################################################
  # Cette fonction renvoie un data frame qui est soit 
  #   le slot clusters (argument clusters.only = T)
  #   un merge des slots clusters et rings 
  #     (Dans ce cas ne sont transformes que les rings correspondant aux valeurs de N et R fournies en arguments)
  #     Si aucune valeur de N et de R n'est fournies , l'ensembes des slots rings sont assemblees (merge) avec clusters
  ###############################################################################################
  out = NULL
  if(!clusters.only && is.prevR(x,"rings")){
    rings = slot(x,"rings")
    if(is.null(N)) N = sapply(rings,function(x) x$N)
    if(is.null(R)) R = sapply(rings,function(x) x$R)
    .isInputOk.prevR(N = N, R = R)
    couples  = unique(data.frame(N = N, R = R, stringsAsFactors=F))
    ringsNames = paste("N",couples[,1],".R",couples[,2],sep="")
    ind = match(ringsNames, names(rings),nomatch = 0) 
    if(sum(ind!=0)) {
      ringsNames = names(rings)[ind] 
      ring       = slot(x,"rings")[[ringsNames[1]]]
      clusters   = slot(x,"clusters")
      out = merge(clusters,ring$estimates,by="id")
      if(length(ringsNames) > 1){
        for(one.ring in ringsNames[-1]){
          ring  = slot(x,"rings")[[one.ring]]
          out   = merge(out,ring$estimates,by="id")
        }
      }
      indId      = match("id",names(ring$estimate))
      namesCol   = c(names(clusters),outer(names(ring$estimate)[-indId],ringsNames,paste,sep="."))
      names(out) = namesCol
    }
  }
  if(clusters.only || !is.prevR(x,"rings")){
    out = slot(x,"clusters")
  }
  out
}

 
