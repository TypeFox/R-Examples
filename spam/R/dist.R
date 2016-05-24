# This is file ../spam/R/dist.R
# This file is part of the spam package, 
#      http://www.math.uzh.ch/furrer/software/spam/
# by Reinhard Furrer [aut, cre], Florian Gerber [ctb]
     



### in base:
# dist(x, method = "euclidean", diag = FALSE, upper = FALSE, p=2)
#
### in fields
# rdist( x1, x2)
# rdist.earth(x1, x2, miles = TRUE, R = NULL)
# fields.rdist.near( x1, x2, delta, max.points= NULL)
#
### in sp 
# spDistsN1(pts, pt, longlat=FALSE)
#
### in amap     ### nbproc 	integer, Number of subprocess for parallelization
# Dist(x, method = "euclidean", nbproc = 1, diag = FALSE, upper = FALSE)
#
### in argosfilter
# distance(lat1, lat2, lon1, lon2)  # gc between two pts in km
# distanceTrack(lat,lon)            # gc between pts  in km
#
### in proxy
# dist(x, y = NULL, method = NULL, ..., diag = FALSE, upper = FALSE,
#     pairwise = FALSE, by_rows = TRUE, convert_similarities = TRUE,
#     auto_convert_data_frames = TRUE)
#
### in RFOC
#  GreatDist(LON1, LAT1, LON2, LAT2, EARTHRAD= 6371)
#

nearest.dist <- function( x, y=NULL, method = "euclidean",
                         delta = 1,
                         upper = if(is.null(y)) FALSE else NULL,
                         p = 2, miles=TRUE, R=NULL
#                         eps =  NULL, diag = NULL
                         )
{
  # see help for exact parameter meaning

  # We always include all small distances. Hence, this function 
  #   works different than any other spam functions. An addititonal
  #   call to an as.spam would eliminate the small values. 
#  if (!is.null(diag)) warning("Argument 'diag' is deprecated")
#  if (!is.null(eps))  warning("Argument 'eps' is deprecated")
  
  if (!is.na(pmatch(method, "euclidian")))     method <- "euclidean"
  METHODS <- c("euclidean", "maximum", "minkowski", "greatcircle")
  method <- pmatch(method, METHODS)  # result is integer

  if (is.na(method))     stop("invalid distance method")

  if (method == 4) {
    if (is.null(R))
      p <- ifelse( miles,3963.34,6378.388)
    else {
      if (R <= 0)           stop("'R' should be postiive")
      p <- R
    }
    if (abs(delta)>180.1)  stop("'delta' should be smaller than 180 degrees.")
  }
  

  if (is.null(upper)) 
    part <- 0L
  else
    part <- ifelse(upper, 1L ,-1L)
  if (is.data.frame(x))  x <- as.matrix(x)
  if (is.list(x)) stop("'x' should be an array or matrix")
           # as.matrix( list() ) does not work
  if (!is.matrix(x)) x <- as.matrix(x)
  nd <- dim(x)[2]
  n1 <- dim(x)[1]

  if (!is.null(y)) {
    # we specify x and y:
    if (is.data.frame(y))  y <- as.matrix(y)
    if (is.list(x)) stop("'x' should be an array or matrix")
    if (!is.matrix(y)) y <- as.matrix(y)
    if (nd!=dim(y)[2]) stop("'x' and 'y' should have the same number of columns.")
    n2 <- dim(y)[1]
    mi <- min(n1,n2)
    ma <- max(n1,n2)
    nnz <- min(max(.Spam$nearestdistnnz[1],
                   ma*.Spam$nearestdistnnz[2]),
               (as.double(mi)*(mi+1)+(ma-mi)^2)/ ifelse( is.null(upper), 1, 2),
               2^31-2)
    # there is an as.double just in case that mi (and n1 below) is > 2^16
  } else {
    # x = y, i.e. proper distance matrix
    if (n1==1)         stop("More than a single point in 'x' is required.")
    if (method == 4) {
      p <- -p  # we save one variable...
    }
    y <- x
    n2 <- n1
    nnz  <- min(max(.Spam$nearestdistnnz[1],
                    n1*.Spam$nearestdistnnz[2]),
                (as.double(n1)*(n1+1))/ ifelse( is.null(upper), 1, 2),
                2^31-2)
  }
  repeat {
    d <- .Fortran("closestdist", nd, as.double(x), n1,  as.double(y), n2, 
                  part,
                  as.double(p[1]), method, 
                  as.double(abs( delta[1])),
                  colindices=vector("integer",nnz),
                  rowpointers=vector("integer",n1+1),
                  entries=vector("double",nnz),
                  nnz=as.integer(nnz),
                  iflag=as.integer(0),NAOK=.Spam$NAOK,
                  PACKAGE="spam")
    
    if (d$iflag==0) break else {
      if (nnz==2^31-2)
        stop("distance matrix is too dense (more than 2^31 entries).")
      nnz <-  min(2^31-2,nnz*.Spam$nearestdistincreasefactor*n1/(d$iflag-1))
      madens <- d$iflag
      on.exit(
              warning(paste("You ask for a 'dense' spase distance matrix, I require one more iteration.",
                            "\nTo avoid the iteration, increase 'nearestdistnnz' option to something like\n",
                            "'spam.options(nearestdistnnz=c(",d$nnz,",400))'\n(constructed ",madens,
                            " lines out of ",n1,").\n",sep=""), call. = TRUE)
              )
    }
  }


  dmat <- new("spam")
  slot(dmat,"entries",check=FALSE) <-     d$entries[1:d$nnz]
  slot(dmat,"colindices",check=FALSE) <-  d$colindices[1:d$nnz]
  slot(dmat,"rowpointers",check=FALSE) <- d$rowpointers
  slot(dmat,"dimension",check=FALSE) <-   as.integer(c(n1,n2))
  return( dmat)
}

# in fields:
# rdist <- function (x1, x2) 

spam_rdist <- function(x1, x2, delta = 1) 
     nearest.dist(x1, y=x2,   delta = delta,  upper = NULL)

# in fields:
# rdist.earth <- function (x1, x2, miles = TRUE, R = NULL) 
spam_rdist.earth <- function(x1, x2, delta=1, miles = TRUE, R = NULL)
    nearest.dist( x1, y=x2, method = "greatcircle",
                         delta = delta, miles=miles, R=R,  upper = NULL)

