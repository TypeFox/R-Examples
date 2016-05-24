############################### 'equality' is a structure which each element is [a unique named value, compatible with unlist()]
# equality is a list(varname1=value1,varname2=value2,...)
# uses rcdd if (include.Hrepr); otherwise depends only on convhulln (ie geometry) through redundant.addVeq
subHullWrapper <- function(vertices, equality, include.Hrepr=FALSE,precision="rational") {
  tmp <- vertices
  for (ii in seq(length(equality))) {
    value <- unlist(equality[ii]) ## named scalar
    ## special case where subHullWrapper should probably not have been called anyway:
    #     if (ncol(tmp)==1L) {
    #       if (colnames(vertices)!=names(value)) {stop("From subHullWrapper(): constraint (from 'equality') involves a variable not in 'vertices'.")}
    #       vrange <- range(vertices)
    #       if (vrange[1]>value || value> vrange[2]) {
    #         tmp <- NULL
    #       } else tmp <- matrix(value,ncol=1,dimnames=list(NULL,names(value)))
    #     } else
    tmp <- redundant.addVeq(vertices=tmp, value=value)
  }
  if ( ! is.null(tmp) && include.Hrepr ) {
    if (precision == "rational") {
      Hrepr <- q2d(scdd(d2q( cbind(0, 1, tmp) ), representation="V", roworder="maxcutoff")$output) ## FR->FR bottleneck
    } else {
      Hrepr <- try(scdd(cbind(0, 1, tmp) , representation="V", roworder="maxcutoff")$output,silent=TRUE)
      if (inherits(Hrepr,"try-error")) { ##  || class(Hrepr) != "matrix") { ?
        Hrepr <- q2d(scdd(d2q( cbind(0, 1, tmp) ), representation="V", roworder="maxcutoff")$output) ## FR->FR bottleneck
      }
    }
    colnames(Hrepr) <- c("eq", "b", colnames(tmp) )
    if (any(abs(Hrepr)>1e16)) { ## Problems sometimes happen despite the core rational computation
      ##  a more meainingful test would be to look whether (
      # check <- Hrepr[, colnames(tmp)] %*% t(tmp) -constraints[, 2]
      ## is widely < 0 for any vertex; Diagnostic plot:
      #  plot(tmp)
      #  lines(tmp[abs(check[1,])<1e-6,],col="red")  ## and so on for each line = each constraint
      ## patch deduced from a single case and clearly heuristic:
      Hrepr <- Hrepr[ ((apply(abs(Hrepr),1,max)<1e16)) <1e16,,drop=FALSE] ## remove suspect constraints
    }
  } else Hrepr <- NULL
  return(list(vertices=tmp, Hrepr=Hrepr))
}
####### Hrepr used in profileBySubHull -> generateInitpts -> rhull -> scdd.addHline
