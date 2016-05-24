######################################################################################################################################
######################################################################################################################################

#Supplemental R code for plotting ICCs, TCCs, IIF, and TIFs from flexMIRT output files

#Authors: Li Cai & Carrie R. Houts
#Date: 10/12/2012

######################################################################################################################################
######################################################################################################################################

#The R programs provided are free software; you can redistribute and/or modify them under the terms of the GNU General Public License 
#(http://www.gnu.org/copyleft/gpl.html) as published by the Free Software Foundation; either version 2 of the License, or 
#(at your option) any later version.  The programs are distributed in the hope that they will be useful, but WITHOUT ANY WARRANTY;
#without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
#See the GNU General Public License for more details.

######################################################################################################################################
######################################################################################################################################

##' Read a flexMIRT PRM file
##'
##' Load the item parameters from a flexMIRT PRM file.
##'
##' @param fname file name
##' @return a list of groups each with item parameters and the latent distribution
##' @export
read.flexmirt <- function(fname) {
  groups <- list()

  ncol <- max(count.fields(fname, sep = "\t"))
  prm <- read.delim(fname,header=F, col.names=1:ncol)

  # find out how many rows
  nlines <- nrow(prm)

  # determine the number of groups
  ng <- 0
  for (i in 1:nlines) {
    if (prm[i,1] == 0) { # group
      ng <- ng+1
    }
  }

  for (g in 1:ng) {
    # select all items and this group
    index <- (prm[,3] == g)
    thisGroup <- prm[index,]

    g.name <- "?"
    g.dist <- list()
    g.label <- list()
    g.spec <- list()
    g.param <- list()

    for (i in 1:nrow(thisGroup)) {
      spec <- NULL
      param <- c()
      
      dims <- thisGroup[i,4]

      if (thisGroup[i,1] == 1) { # item
        if (thisGroup[i,5] == 1) { # 3PL
          # grab item parameters
          logitg <- thisGroup[i,7]
          c <- thisGroup[i,8]
          a <- thisGroup[i,9:(9+dims-1)]

          spec <- rpf.drm(factors = dims)
          param <- c(unlist(a), c, logitg, logit(1))
        }
        if (thisGroup[i,5] == 2) { # graded
          # grab item parameters
          nc <- thisGroup[i,6]
          c <- rep(0,nc-1)
          for (k in 1:(nc-1)) {
            c[k] <- thisGroup[i,6+k]
          }
          a <- thisGroup[i,(6+nc):(6+nc+dims-1)]

          spec <- rpf.grm(outcomes=nc, factors=dims)
          param <- c(unlist(a), c)
        }
        if (thisGroup[i,5] == 3) { # nominal
          offset <- 6
          # grab item parameters
          nc <- thisGroup[i,offset]
          offset <- offset+1
          Tmat <- thisGroup[i,offset]
          offset <- offset+1
          if (Tmat == 0) T <- 'trend'
          if (Tmat == 1) T <- 'id'
          if (Tmat == 2) { # user-defined
            v <- matrix(0,nc*(nc-1))
            for (k in 1:(nc*(nc-1))) {
              v[k] <- thisGroup[i,offset]
              offset <- offset+1
            }
            T <- matrix(c(v), nc, nc-1, byrow=TRUE)
          }
          # scoring fn
          alf <- rep(0,nc-1)
          for (k in 1:(nc-1)) {
            alf[k] <- thisGroup[i,offset]
            offset <- offset+1
          }
          # slope
          a <- thisGroup[i,offset]
          offset <- offset+1

          # intercept
          Lmat <- thisGroup[i,offset]
          offset <- offset+1
          if (Lmat == 0) L <- 'trend'
          if (Lmat == 1) L <- 'id'
          if (Lmat == 2) { # user-defined
            v <- matrix(0,nc*(nc-1))
            for (k in 1:(nc*(nc-1))) {
              v[k] <- thisGroup[i,offset]
              offset <- offset+1
            }
            L <- matrix(c(v), nc, nc-1, byrow=TRUE)
          }
          gam <- rep(0,nc-1)
          for (k in 1:(nc-1)) {
            gam[k] <- thisGroup[i,offset]
            offset <- offset+1
          }

          spec <- rpf.nrm(outcomes = nc, T.a = T, T.c = L)
          param <- c(a, alf, gam)
        }

        if (!is.null(spec)) {
          g.label[[i]] <- as.character(thisGroup[i,2])
          g.spec[[i]] <- spec
          g.param[[i]] <- param
        }
      } else { # group
        cov <- matrix(NA, nrow=dims, ncol=dims)
        col <- 6+dims
        for (rx in 1:dims) {
          for (cx in 1:rx) {
            cov[rx,cx] <- thisGroup[i,col]
            if (rx != cx) cov[cx,rx] <- thisGroup[i,col]
            col <- col+1
          }
        }
        g.dist <- list(mean=unlist(thisGroup[i,6:(6+dims-1)]), cov=cov)
        g.name <- thisGroup[i,2]
      }
    } # for every item

    names(g.spec) <- g.label
    pmat <- matrix(NA, ncol=length(g.param), nrow=max(vapply(g.param, length, 0)))
    for (i in 1:length(g.param)) {
      v <- g.param[[i]]
      pmat[1:length(v), i] <- v
    }
    colnames(pmat) <- g.label
    groups[[g]] <- list(spec=g.spec, param=pmat, mean=g.dist$mean, cov=g.dist$cov)
    names(groups)[[g]] <- as.character(g.name)
  }  # for every group

  groups
}

serialize.T <- function(spec, Tmat, T) {
  if (all(abs(T - Tnom.trend(spec@outcomes)) < 1e-4)) {
    c(0)
  } else if (Tmat == 'a' && all(abs(T - Tnom.ida(spec@outcomes)) < 1e-4)) {
    c(1)
  } else if (Tmat == 'c' && all(abs(T - Tnom.idc(spec@outcomes)) < 1e-4)) {
    c(1)
  } else {
    c(2, rep(0, dim(T)[1]), t(T))
  }
}

##' Write a flexMIRT PRM file
##'
##' Formats item parameters in the way that flexMIRT expects to read
##' them. Use \code{\link{read.flexmirt}} to see what shape the groups
##' parameter of this function should take.
##'
##' NOTE: Support for the graded response model may not be complete.
##'
##' @param groups a list of groups each with items and latent parameters
##' @param file the destination file name
##' @param fileEncoding how to encode the text file (optional)
write.flexmirt <- function(groups, file=NULL, fileEncoding="") {
  if (missing(file)) {
    stop("You must specify the destination file=")
  }
  if (file == "") 
    file <- stdout()
  else if (is.character(file)) {
    file <- if (nzchar(fileEncoding)) 
      file(file, "w", encoding = fileEncoding)
    else file(file, "w")
    on.exit(close(file))
  }
  else if (!isOpen(file, "w")) {
    open(file, "w")
    on.exit(close(file))
  }
  if (!inherits(file, "connection")) 
    stop("'file' must be a character string or connection")
  
  for (gx in 1:length(groups)) {
    grp <- groups[[gx]]
    for (ix in 1:length(grp$spec)) {
      spec <- grp$spec[[ix]]
      nfact <- spec@factors
      name <- names(grp$spec)[ix]
      iparam <- grp$param[,ix]
      if (is.null(name)) name <- paste("i",ix,sep="")
      if (class(spec) == "rpf.mdim.drm") {
        if (iparam[4] < 1) warnings("Nonzero upper bound ignored")
        if (iparam[3] == 0) {
          warnings("Guessing set to 0; using Graded(2) instead")
          cat(paste(1, name, gx, nfact, 2, 2, iparam[2], iparam[1], sep="\t"),
              file=file, fill=TRUE)
        } else {
          logitg <- iparam[3]
          cat(paste(1, name, gx, nfact, 1, 2, logitg, iparam[2], iparam[1], sep="\t"),
              file=file, fill=TRUE)
        }
      } else if (class(spec) == "rpf.mdim.grm") {
        if (spec@outcomes > 2) stop("Not implemented")
        cat(c(1, name, gx, nfact, 2, spec@outcomes, iparam[(nfact+1):length(iparam)],
                  iparam[1:nfact]), sep="\t", file=file, fill=TRUE)
      } else if (class(spec) == "rpf.mdim.nrm") {
        T.a <- serialize.T(spec, 'a', getT(spec, 0))
        T.c <- serialize.T(spec, 'c', getT(spec, 1))
        cat(paste(c(1, name, gx, nfact, 3, spec@outcomes, T.a,
                    iparam[(spec@factors+1):(spec@outcomes + spec@factors - 1)],
                    iparam[1:spec@factors],
                    T.c, iparam[(spec@factors + spec@outcomes):length(iparam)]), sep="\t"),
            file=file, fill=TRUE)
      } else {
        stop(paste("Not implemented for", class(spec)))
      }
    }
    tri <- rep(NA, nfact*(nfact+1)/2)
    tx <- 1
    for (rx in 1:nfact) {
      for (cx in 1:rx) {
        tri[tx] <- grp$cov[rx,cx]
        tx <- tx+1
      }
    }
    cat(c(0, names(groups)[gx], gx, nfact, 0, grp$mean, tri),
        sep="\t", file=file, fill=TRUE)
  }
}
