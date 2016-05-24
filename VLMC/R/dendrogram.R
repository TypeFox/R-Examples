####- vlmctree() & as.dendrogram() --- R-level recursive tree construction

vlmctree <- function(x)
{
  ## Purpose: Compute the Tree representation of a "vlmc" object (as R list).
  ## -------------------------------------------------------------------------
  ## Arguments: x: a "vlmc" object {usually a fit from vlmc(..)}.
  ## -------------------------------------------------------------------------
  ## Author: Martin Maechler, Date:  1 Apr 2000, 18:02
  if(!is.vlmc(x)) stop("first argument must be a \"vlmc\" object; see ?vlmc")
  vvec <- (x $ vlmc.vec)#.Alias
  k <- (x $ alpha.len)#.Alias
  if(vvec[1] != k) stop("invalid vlmc structure {alpha.len}")

  ## return vtree :
  .vvec2tree(vvec[-1], k = vvec[1], chk.lev = 0)
}

.vvec2tree <- function(vv, k, chk.lev)
{
  ## Purpose: RECURSIVELY construct tree from a vvec of a "vlmc" object
  ##	      *not* using alphabet, (just k = |alphabet|).
  ## Do as load_tree(.)	 {in ../src/io.c }
  ## -------------------------------------------------------------------------
  ## Author: Martin Maechler, Date:  1 Apr 2000, 18:11

  if((lev <- vv[1]) >= 0) { ## non-end
      if(lev != chk.lev)
	  stop(paste("malformed vlmc tree at level",chk.lev))

      ii <- 1:k
      node <- list(level = lev, count = vv[1 + ii], child = vector("list", k))
      node $ total <- sum(node $ count)

      vv <- vv[ - c(ii, k+1)]# the first 1..(k+1) ones
      for(i in ii) {
	  r <- .vvec2tree(vv, k=k, chk.lev = chk.lev+1)
          vv <-
              if(!is.null(r)) {
                  node$child[[i]] <- r[[1]]
                  r[[2]]
              } else vv[-1] ## child[i] remains NULL
      }
      node$level[2] <- ## parent level :
          if(all(sapply(node$child, is.null))) { ## this is a leaf
              node$child <- NULL
              0 # parent level
          } else 1L + max(sapply(node$child, function(n)
                                  if(is.null(n)) 0 else n$level[2]))
      node$level <- as.integer(node$level)
      if(lev > 0)
	  list(node, vv)
      else { ## toplevel
	  class(node) <- "vtree"
	  node
      }
  }
  ## else return NULL
}

str.vtree <- function(object, ...)
{
    ## Purpose: str method for "vtree" lists  [[recursive]]
    if(!is.null(lev <- object$level)) {
	nch <- length(object$child)
	cat(if(lev[1])
	    paste(rep(" ", lev[1]+1), collapse="..") else "'vtree':\n",
            paste("{", lev[2],"}", sep=""),
	    format(object$count),"|", object$total, "; ",
	    if(nch) paste(nch,"children") else "_leaf_",
	    "\n")
	for(ch in object$child)
	    str.vtree(ch, ...)
    }
}


###- as.dendrogram() method - in order to plot the context - tree

## Generic and hclust method are in
##  ~/R/D/r-devel/R/src/library/stats/R/dendrogram.R

### FIXME:
### =====

## (*) Add "midpoint" such that I can plot with center = FALSE
##     (I *can* plot, but it's not okay;   center = TRUE is okay)

## CONSIDER:  Make this a new class *inheriting* from dendrogram
##            The new class would have its own print and plot method...

as.dendrogram.vlmc <- function(object, ...)
{
    if(!is.vlmc(object))
        stop("first argument must be a \"vlmc\" object; see ?vlmc")
    vvec <- (object $ vlmc.vec)#.Alias
    k <- object $ alpha.len
    if(vvec[1] != k) stop("invalid vlmc structure {alpha.len}")
    p <- unname(object $ size["ord.MC"]) # maximal MC order
    abc <- alphabet(object)

    vv2dendro <- function(vv, cl.lev)
    {
        ## construct the nested list, level 'cl.lev' from 'vvec' -- recursively!
        if((lev <- vv[1]) >= 0) { ## non-end
            if(lev != cl.lev)
                stop(paste("malformed vlmc tree at level",cl.lev))

            ii <- 1:k
            node <- vector("list", k)   # k children
            names(node) <- 0:(k-1)
            count <- vv[1 + ii] 	# and their counts
            vv <- vv[ - c(ii, k+1)]     # drop the first 1..(k+1) ones
            for(i in ii) { ## extract child[i] (and *its* children)
                r <- vv2dendro(vv, cl.lev = cl.lev+1)
                ## downdating 'vv',  updating node[[i]]:
                vv <-
                    if(!is.null(r)) {
                        node[[i]] <- r[[1]]
                        attr(node[[i]], "edgetext") <- abc[i]

                        r[[2]]
                    }
                    else ## empty child[i]; drop NULL node[[i]] below
                        vv[-1]
            }
            ##- cat("lev=",lev,";", "count=",count,"  vv : \n"); str(vv)
            if(all(kids0 <- sapply(node, is.null))) { ## this is a leaf
                node <- sum(count)
                attr(node,"members") <- 1L
                attr(node,"leaf") <- TRUE
            }
            else { ## has at least one child;
                node[kids0] <- NULL ## drop the NULL ones (but remember which!)
                attr(node,"0-kids") <- (0:(k-1))[kids0]
                ## attr(node,"height") <- ## parent level :
                ##    1L + max(sapply(node, function(n) attr(n,"height")))
                attr(node,"members") <- ## parent level :
                    sum(sapply(node, function(n) attr(n,"members")))
            }
            attr(node,"height") <- p - lev
            ## keep the full count[] :
            attr(node,"count") <- count

            list(node, vv)
        }
        ## else vv[1] = -1 :  return NULL
    }

    r <- vv2dendro(vvec[-1], 0)[[1]]
    class(r) <- "dendrogram"
    r
}

