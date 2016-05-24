#
#  multinomRob
#
#  Walter R. Mebane, Jr.
#  Cornell University
#  http://macht.arts.cornell.edu/wrm1/
#  wrm1@macht.arts.cornell.edu
#
#  Jasjeet Singh Sekhon 
#  UC Berkeley
#  http://sekhon.polisci.berkeley.edu
#  sekhon@berkeley.edu
#
#  $Id: datamanip.R,v 1.7 2005/09/27 08:04:06 wrm1 Exp $
#

#Mapping from xvec to the beta.vector and back again
#forward==TRUE from xvec TO beta.vector
#forward==FALSE from beta.vector TO xvec
mnl.xvec.mapping <- function (forward=TRUE,base.xvec,work.xvec,beta.vector,
                             ncats,nvars.total) 
{
  n.ones <- sum(base.xvec==1);  # indicates unique parameters
  n.mults <- sum(base.xvec>1);  # indicates parameters constrained equal
  if (forward) {
    p <- 0;
    if (n.ones>0) {
      for (j in 1:ncats) {
        for (i in 1:nvars.total) {
          if (base.xvec[i,j] == 1) {
            p <- p + 1;
            beta.vector[p] <- work.xvec[i,j];
          } #end of if
        } #end of i
      } #end of j
    }
    if (n.mults>0) {
      nidxvals <- length(idxvals <- sort(unique(base.xvec[base.xvec>1])));
      for (k in 1:nidxvals) {
        for (j in 1:ncats) {
          if (any(jktest <- base.xvec[,j]==idxvals[k])) {
            beta.vector[p+k] <- work.xvec[jktest,j][1];
            break;
          }
        }
      }
    }
    return(beta.vector);
  } else {
    p <- 0;
    if (n.ones>0) {
      for (j in 1:ncats) {
        for (i in 1:nvars.total) {
          if (base.xvec[i,j] == 1) {
            p <- p + 1;
            work.xvec[i,j] <- beta.vector[p];
          } #end of if
        } #end of i
      } #end of j
    }
    if (n.mults>0) {
      nidxvals <- length(idxvals <- sort(unique(base.xvec[base.xvec>1])));
      for (k in 1:nidxvals) {
        for (j in 1:ncats) {
          if (any(jktest <- base.xvec[,j]==idxvals[k])) {
            work.xvec[jktest,j] <- beta.vector[p+k];
          }
        }
      }
    }
    return(work.xvec);
  } # end of else
} #end of mnl.xvec.mapping


############################################################################
## Create the jacstack (from tanh)                                     #
############################################################################    

# jacstack.function:  arrange regressors for computing the Jacobian matrix
##   produces jacstack:  array of regressors,
##      dim(jacstack) = c(n observations, n UNIQUE parameters, n categories)
jacstack.function <- function(X,tvars.unique,xvec) {
  xdim  <- dim(X)
  obs   <- xdim[1]
#  nvars <- xdim[2]
  ncats <- xdim[3]

  jacstack <- array(0,dim=c(obs,tvars.unique,ncats));
  n.ones <- sum(xvec==1);  # indicates unique parameters
  n.mults <- sum(xvec>1);  # indicates parameters constrained equal
  itu <- 0;
  if (n.ones>0) {
    for (j in 1:ncats) {
      nxj <- sum(xvec[,j]==1);
      if (nxj>0) {
        jacstack[1:obs, itu + 1:nxj, j] <- X[rep(TRUE,obs),xvec[,j]==1,j];
      }
      itu <- itu + nxj;
    }
  }
  if (n.mults>0) {
    nxvecrows <- dim(xvec)[1];
    nidxvals <- length(idxvals <- sort(unique(xvec[xvec>1])));
    for (k in 1:nidxvals) {
      for (j in 1:ncats) {
        if (any(jktest <- xvec[,j]==idxvals[k])) {
          kidx <- (1:nxvecrows)[jktest];  # xvec[,j] rows matching constraint k
          nxj <- sum(jktest);
          for (kk in 1:nxj) {
            jacstack[1:obs, itu + k, j] <-
              jacstack[1:obs, itu + k, j] + X[rep(TRUE,obs),kidx[kk],j];
          }
        }
      }
    }
  }
  return(jacstack)
} #jacstack.function

# jacstack.singles:  check jacstack array for regressors with a distinct value
#                    at only one observation
jacstack.singles <- function(jacstack) {
  nunique <- dim(jacstack)[2];

  jsingle <- rep(FALSE, nunique);
  for (i in 1:nunique) {
    if (length(unique(c(jacstack[,i,]))) == 2) {
      jsingle[i] <- ifelse(any(table(c(jacstack[,i,])) == 1), TRUE, FALSE)
    }
  }
  return(jsingle)
} #jacstack.singles.function
