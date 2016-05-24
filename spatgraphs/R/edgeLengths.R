#' Edge lengths
#'
#'
#' @param g sg-object
#' @param x point pattern
#' @param ... ignored

#' @export

edgeLengths<-function(g, x, ...)
{
  if(missing(x)) stop("Need 'x' for distances.")
  if(!is(g, "sg")) stop("Give sg object, from spatgraph-function.")
  res<-list()
  ivec <- jvec <- dvec <- NULL

  x <- sg_parse_coordinates(x)

  for(i in 1:g$N)
  {
    iedges<-g$edges[[i]]
    for(j in iedges)
    {
      ivec<-c(ivec, i)
      jvec<-c(jvec, j)
      d <- sqrt( sum((x[i,]-x[j,])^2) )
      dvec<-c(dvec, d)
      g$edges[[j]]<-setdiff(g$edges[[j]], i)
    }
  }

  res$i<-ivec
  res$j<-jvec
  res$d<-dvec
  res$n<-length(res$i)
  res
}
