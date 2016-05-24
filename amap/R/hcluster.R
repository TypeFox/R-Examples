## Hierarchical clustering
##
## Created       : 18/11/02
## Last Modified : Time-stamp: <2014-12-17 19:06:24 antoine>
##
## This function is a "mix" of function dist and function hclust.
##
## Author : Antoine Lucas
##



hclusterpar <- hcluster <- function (x, method = "euclidean", diag = FALSE, upper = FALSE, link = "complete", members = NULL, nbproc = 2, doubleprecision = TRUE)
{

  if(class(x) == "exprSet")
    x <- Biobase::exprs(x)

  ## take from dist
  if (!is.na(pmatch(method, "euclidian"))) 
    method <- "euclidean"
  METHODS <- c("euclidean", "maximum", "manhattan", "canberra", 
               "binary","pearson","correlation","spearman","kendall",
               "abspearson","abscorrelation")
  method <- pmatch(method, METHODS)
  if (is.na(method)) 
    stop("invalid distance method")
  if (method == -1) 
    stop("ambiguous distance method")
  N <- nrow(x <- as.matrix(x))
  

  
                                        #take from hclust
  METHODSLINKS <- c("ward", "single", "complete", "average", "mcquitty", 
                    "median", "centroid","centroid2")
  
  link <- pmatch(link, METHODSLINKS)
  if (is.na(link)) 
    stop("invalid clustering method")
  if (link == -1) 
    stop("ambiguous clustering method")
    if (N < 2) 
        stop("Must have n >= 2 objects to cluster")
  if (is.null(members)) 
    members <- rep(1, N)
  if (length(members) != N) 
    stop("Invalid length of members")
  n <- N

  precision <- 1
  if(doubleprecision)
    precision <- 2
  
  hcl <- .C("hcluster",
            x = as.double(x),
            nr = as.integer(n),
            nc = as.integer(ncol(x)),
            diag = as.integer(FALSE),
            method = as.integer(method), 
            iopt = as.integer(link),
            ia = integer(n),
            ib = integer(n),
            order = integer(n),
            crit = double(n),
            members = as.double(members),
            nbprocess  = as.integer(nbproc),
            precision  = as.integer(precision),
            res  = as.integer (1),
            NAOK=TRUE,
            PACKAGE= "amap")

  if(hcl$res == 2)
    stop("Cannot allocate memory")
  if(hcl$res == 3)
    stop("Missing values in distance Matrix")
  if(hcl$res == 1)
    stop("Error")


  tree <- list(merge = cbind(hcl$ia[1:(N - 1)],
                 hcl$ib[1:(N -  1)]),
               height = hcl$crit[1:(N - 1)],
               order = hcl$order, 
               labels = dimnames(x)[[1]],
               method = METHODSLINKS[link],
               call = match.call(),
               dist.method = METHODS[method]
               )


  class(tree) <- "hclust"
  tree
}
