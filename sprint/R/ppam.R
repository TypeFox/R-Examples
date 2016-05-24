##########################################################################
#                                                                        #
#  SPRINT: Simple Parallel R INTerface                                   #
#  Copyright © 2008,2009 The University of Edinburgh                     #
#                                                                        #
#  This program is free software: you can redistribute it and/or modify  #
#  it under the terms of the GNU General Public License as published by  #
#  the Free Software Foundation, either version 3 of the License, or     #
#  any later version.                                                    #
#                                                                        #
#  This program is distributed in the hope that it will be useful,       #
#  but WITHOUT ANY WARRANTY; without even the implied warranty of        #
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the          #
#  GNU General Public License for more details.                          #
#                                                                        #
#  You should have received a copy of the GNU General Public License     #
#  along with this program. If not, see <http://www.gnu.or/licenses/>.   #
#                                                                        #
##########################################################################

## Ensure consistent "diss.." class --- make "namespace-private-global" !
dissiCl <- c("dissimilarity", "dist")

ppam <- function(x, k, medoids = NULL, is_dist = inherits(x, "dist"),
                 cluster.only = FALSE, do.swap = TRUE, trace.lev = 0)
{	
	if (is_dist){
		xLab <- attr(x, "Labels")
	}else{
		xLab <- rownames(x)}
	
  # === CHECK DISTANCE DATA INPUT === #
  if(is.ff(x)) {
    ## check type of input ff object
    if(vmode(x) != "double")
      stop(..sprintMsg$error["non.square"])
    
    if(data.class(x) != "ff_matrix") {
      length_x <- attr(attr(x, "virtual"), "Length")
      n_rows <- sqrt(length_x)
      
      # check if the vector can form a square matrix
      if (length_x %% n_rows)
        stop(..sprintMsg$error["non.square"])
    } else {
      n_rows = dim.ff(x)[1]
      if (dim.ff(x)[1] != dim.ff(x)[2]) {
        stop(..sprintMsg$error["non.square"])
      }
    }
    filename = attr(attr(x, "physical"), "filename")
    if (!is.character(filename)) {
      stop(..sprintMsg$error["no.filename"])
    }
    
    MAP_FILE= TRUE
	  x = 0;

  } else {
    
    if(!is.numeric(x)) stop("x is not a numeric dataframe or matrix.")

    if(is_dist) 
      x = as.matrix(x)

    n_rows <- nrow(x)
    filename = ""
    MAP_FILE <- FALSE
	
  }
  
  if((k <- as.integer(k)) < 1 || k >= n_rows)
    stop(..sprintMsg$err["no.valid.k"])
  if(is.null(medoids))# default: using "build & swap to determine medoids"
    medID <- integer(k)# all 0 -> will be used as 'code' in C
  else {
    ## 'fixme': consider  sort(medoids) {and rely on it in ../src/pam.c }
    if(length(medID <- as.integer(medoids)) != k ||
       any(medID < 1) || any(medID > n_rows) || any(duplicated(medID)))
      stop("'medoids' must be NULL or vector of ",
           k, " distinct indices in {1,2, .., n}, n=", n_rows)
    ## use observation numbers  'medID' as starting medoids for 'swap' only
  }
  nisol <- integer(if(cluster.only) 1 else k)

  if(do.swap) nisol[1] <- 1L
  stopifnot(length(cluster.only) == 1,
            length(trace.lev) == 1)

  res <- .C("ppam",
              as.logical(MAP_FILE),
              x,
	      as.integer(n_rows),
	      as.integer(k),
	      as.character(filename),
	      integer(n_rows),		# nsend[]
	      logical(n_rows),		# nrepr[]
	      integer(if(cluster.only) 1 else n_rows), # nelem[]
	      double(n_rows),		# radus[]
	      double(n_rows),		# damer[]
	      avsil = double(n_rows),	# 'ttd'
	      double(n_rows),		# separ[]
	      ttsil = as.double(0),
	      obj = as.double(c(cluster.only, trace.lev)),# in & out!
	      med = medID,# in & out(if !cluster.only)
	      clu = integer(n_rows),
	      clusinf = if(cluster.only) 0. else matrix(0., k, 5),
	      silinf  = if(cluster.only) 0. else matrix(0., n_rows, 4),
	      isol = nisol,
	      PACKAGE = "sprint")

  #Compatibility with pam
  res$silinf[,1] = res$silinf[,1] + 1
  res$silinf[,2] = res$silinf[,2] + 1
  res$silinf[,4] = res$silinf[,4] + 1
	
  if(length(xLab) > 0)
      names(res$clu) <- xLab
  if(cluster.only)
      return(res$clu)

   ## Else, usually
  medID <- res$med+1
  if(any(medID < 0))
    stop("error from .C(\"ppam\", *): invalid medID's")
  sildim <- res$silinf[, 4]
	
# Return medoid dimnames if they exist
	if(length(xLab) > 0){
		res$med <- names(res$clu[medID])
	}
	else{
		res$med = medID
	}
	
  ## add dimnames to Fortran output
  names(res$obj) <- c("build", "swap")
  res$isol <- factor(res$isol, levels = 0:2, labels = c("no", "L", "L*"))
  names(res$isol) <- 1:k
  dimnames(res$clusinf) <- list(NULL, c("size", "max_diss", "av_diss",
                                        "diameter", "separation"))
  
  ## construct S object
  r <-
    list(medoids = res$med, id.med = medID, clustering = res$clu,
         objective = res$obj, isolation = res$isol,
         clusinfo = res$clusinf,
         silinfo = if(k != 1) {
           dimnames(res$silinf) <- list(sildim, c("cluster", "neighbor", "sil_width", ""))
           list(widths = res$silinf[, -4],
                clus.avg.widths = res$avsil[1:k],
                avg.width = res$ttsil)
         },
         call = match.call())
  class(r) <- c("pam", "partition")
  r
  
}

## non-exported:
.print.pam <- function(x, ...) {
    cat("Medoids:\n");		print(cbind(ID = x$id.med, x$medoids), ...)
    cat("Clustering vector:\n");	print(x$clustering, ...)
    cat("Objective function:\n");	print(x$objective, ...)
}


print.pam <- function(x, ...)
{
    .print.pam(x, ...)
    cat("\nAvailable components:\n")
    print(names(x), ...)
    invisible(x)
}

summary.pam <- function(object, ...)
{
    class(object) <- "summary.pam"
    object
}

print.summary.pam <- function(x, ...)
{
    .print.pam(x, ...)
    cat("\nNumerical information per cluster:\n"); print(x$clusinfo, ...)
    cat("\nIsolated clusters:\n L-clusters: ")
    print(names(x$isolation[x$isolation == "L"]), quote = FALSE, ...)
    cat(" L*-clusters: ")
    print(names(x$isolation[x$isolation == "L*"]), quote = FALSE, ...)
    if(length(x$silinfo) != 0) {
	cat("\nSilhouette plot information:\n")
	print(x$silinfo[[1]], ...)
	cat("Average silhouette width per cluster:\n")
	print(x$silinfo[[2]], ...)
	cat("Average silhouette width of total data set:\n")
	print(x$silinfo[[3]], ...)
    }
    if(!is.null(x$diss)) { ## Dissimilarities:
	cat("\n");			print(summary(x$diss, ...))
    }
    cat("\nAvailable components:\n");	print(names(x), ...)
    invisible(x)
}
