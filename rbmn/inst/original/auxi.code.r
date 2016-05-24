
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
aux0 <- function(az)
#TITLE returns a convenient function
#DESCRIPTION \code{sqrt(1-sum(az^2))}.\cr Not intended for the
# final user
#DETAILS
#KEYWORDS 
#INPUTS
#{az} << numerical vector.>>
#[INPUTS]
#VALUE
# sqrt(1-sum(az^2))
#EXAMPLE
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 12_01_13
#REVISED 12_01_18
#--------------------------------------------
{
  sqrt(1-sum(az^2));
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
aux1 <- function(rho,inverse=FALSE)
#TITLE returns the generating matrix of a sequential chain
#DESCRIPTION From a vector of length \code{n-1} of correlations
# computes the generating matrix \code{n x n} of a sequential
# chain with the root in first position (centred and normalized).\cr
# This function is not intented for the end user: much better
# to use \code{chain2gema} which allows also non zero expectations
# and non unity variances.
#DETAILS
#KEYWORDS 
#INPUTS
#{rho} << Vector of the \code{n-1} correlations.>>
#[INPUTS]
#{inverse} <<Indicates if the order of rows and columns must be reverted.>>
#VALUE
# a square matrix, say \code{L} such
# that the vector \code{L%*%W} has identical
# distribution with the chain; \code{W} being
# a white noise.
#EXAMPLE
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 12_01_13
#REVISED 12_02_06
#--------------------------------------------
{
  # no check
  # completing the vector
  rho <- c(0,rho);
  # numbers of nodes
  nn <- length(rho);
  # defining an auxiliary function
  ax <- function(h,k) {
    res <- 0;
    if (h==k+1) {
      res <- 1;
    } else {
      res <- prod(rho[h:k]);
    }
    res;
  }
  # initialization
  res <- matrix(0,nn,nn);
  # filling
  res[1,1] <- 1;
  for (ii in r.bd(2,nn)) {
    for (jj in r.bd(1,ii)) {
      res[ii,jj] <- aux0(rho[jj])*ax(jj+1,ii);
    }
  }
  # inverting
  if (inverse) {
    res <- res[rev(r.bc(nn)),rev(r.bc(nn))];
  }
  # returning
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
aux2 <- function(chain,ordering=FALSE,parents=FALSE)
#TITLE arc orientations, ordering of a /chain/, the parents of 
# each node...
#DESCRIPTION returns a value for the
# \code{n-1} arcs of a chain: 1 if from
# \code{i} to \code{i+1} -1 if not when \code{!ordering}.\cr
# One of the possible topological order when \code{ordering}.\cr
# The parents when \code{parents}, that is \code{0} for root nodes,
# \code{-1} for colliders and the parent number for other nodes.\cr
# Not intended for the final user
#DETAILS
# \code{ordering} and \code{parents} cannot be both \code{TRUE}.
#KEYWORDS 
#INPUTS
#{chain} << The chain to consider.>>
#[INPUTS]
#{ordering} <<\code{TRUE} to get
# a possible topological order.>>
#{parents} <<\code{TRUE} to get the parents.>>
#VALUE
# A named vector of plus or minus ones when
#  \code{!ordering} if not a permutation of the
# nodes following a topological order.
#EXAMPLE
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 12_01_18
#REVISED 12_02_15
#--------------------------------------------
{
  # checking
  # <to be done>
  # getting the necessary information
  nn <- length(chain$names);
  aaa <- sort(r.numero(chain$roots,chain$names));
  if (length(aaa) > 1) {
    ccc <- sort(r.numero(chain$colliders,chain$names));
  } else {
    ccc <- numeric(0);
  }
  # a small check
  if (!(all(ccc-aaa[-1]<0)) |
      !(all(ccc-aaa[-length(aaa)]>0))
     ) {
    r.erreur(chain,message="roots/colliders misplaced");
  }
  ccc <- c(ccc,nn+1)
  if (!ordering) {
    # computing
    res <- rep(-1,nn-1);
    for (aa in r.bf(aaa)) {
      res[r.bd(aaa[aa],ccc[aa]-1)] <- 1;
    }
    # naming
    jun <- c("<-","<>","->")[2+res];
    nres <- paste(chain$names[-nn],jun,chain$names[-1],sep="");
    names(res) <- nres;
  } else {
    if (!parents) {
      # looking for a topological order
      possibles <- r.bc(nn+1);
      res <- numeric(0);
      for (nume in r.bf(aaa)) {
        anc <- aaa[nume];
        coi <- ccc[nume];
        ouanc <- which(anc==possibles);
        oucoi <- which(coi==possibles);
        gauche <- possibles[r.bd(1,ouanc)];
        if (nume > 1) { gauche <- c(gauche[1]-1,gauche);}
        droite <- possibles[r.bd(ouanc+1,oucoi-1)];
        possibles <- possibles[r.bd(oucoi+1,length(possibles))];
        res <- c(res,rev(gauche),droite);
      }
    }
  }
  if (parents) {
    # transforming arc orientations into parents
    rres <- rep(0,length(chain$names));
    names(rres) <- chain$names;
    for (no in r.bf(rres)) {
      if (no==1) {
        if (res[1] == -1) { rres[no] <- 2;}
      }
      if (no==nn) {
        if (res[nn] == 1) { rres[no] <- nn-1;}
      }
      if ((no>1)&(no<nn)) {
        nbp <- 0;
        if (res[no-1] == 1) {
          nbp <- 1;
          rres[no] <- no-1;
        }
        if (res[no] == -1) {
          nbp <- nbp+1;
          rres[no] <- no+1;
        }
        if (nbp == 2) {
          rres[no] <- -1;
        }
      }
    }
    res <- rres;
  }
  # returning
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
aux3 <- function(nam,ordering=NULL)
#TITLE gets the ordering of a vector of names
#DESCRIPTION 
# For most functions \code{print8...} to get
# the vector of integer indices (possibly with
# missing and repetitions.\cr
# Not intended for the final user
#DETAILS
#KEYWORDS 
#INPUTS
#{nam} << The names to consider.>>
#[INPUTS]
#{ordering} << the proposed order either
#              as indices of \code{names}
#              or components of \code{names}.
#              When null, the names as they are.
#VALUE
# A vector.
#EXAMPLE
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 12_01_19
#REVISED 12_01_19
#--------------------------------------------
{
  # checking
  # <to be done>
  # getting the necessary information
  nn <- length(nam);
  if (is.null(ordering)) {
    ordering <- r.bc(nn);
  } else {
    if (is.numeric(ordering)) {
      oo <- ordering;
      ordering <- numeric(0);
      for (rr in oo) {
        if (rr %in% r.bc(nn)) {
          ordering <- c(ordering,match(rr,r.bc(nn)));
        }
      }
    } else {
      if (is.character(ordering)) {
        oo <- ordering;
        ordering <- numeric(0);
        for (rr in oo) {
          if (rr %in% nam) {
            ordering <- c(ordering,match(rr,nam));
          }
        }
      } else {
        # not an accepted object, then default value
        ordering <- r.bc(nn);
      }
    }
  }
  # returning
  ordering;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
aux4 <- function(uu)
#TITLE gets conditional expectations and sd
#DESCRIPTION 
# from an pseudo /chain/ transforms \code{$mu}
# and \code{$sigma} from marginal to conditional
# ones.\cr
# Not intended for the final user
#DETAILS
#KEYWORDS 
#INPUTS
#{uu} << pseudo chain to consider.>>
#[INPUTS]
#VALUE
# Resulting chain
#EXAMPLE
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 12_01_31
#REVISED 12_02_02
#--------------------------------------------
{
    # finding the reference positions
    aaa <- sort(r.numero(uu$roots,uu$names));
    muma <- uu$mu;
    sima <- uu$sigma;
    oooo <- aux2(uu,TRUE);
    orro <- aux2(uu);
    avan <- c(0,orro)== 1;
    apre <- c(orro,0)==-1;
    for (nno in oooo) {
      # only non roots must be modified
      if (!(nno %in% aaa)) {
        # finding the parent(s) and the associated correlations
        if (avan[nno]) { papa <- nno-1; rho <- uu$corre[nno-1];}
        else { papa <- numeric(0); rho <- numeric(0);}
        if (apre[nno]) { papa <- c(papa,nno+1); rho <- c(rho,uu$corre[nno]);}
        if (length(papa)==0) { stop("transforming a chain");}
        # getting the new parameters
        cre <- uu$sigma[nno] * rho / sima[papa];
        uu$sigma[nno] <- uu$sigma[nno] * sqrt(1-sum(rho^2));
        uu$mu[nno] <- uu$mu[nno] - sum(cre*muma[papa]);
      }
    }
  # returning
  uu;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
aux5 <- function(rho)
#TITLE returns the product block from a vector
#DESCRIPTION 
# See code, quite straightforward\cr
# Not intended for the final user
#DETAILS
#KEYWORDS 
#INPUTS
#{rho} << numeric vector>>
#[INPUTS]
#VALUE
# Resulting matrix of correlations (if \code{rho} is)
# of size \code{length(rho)+1}.
#EXAMPLE
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 12_02_01
#REVISED 12_02_01
#--------------------------------------------
{
    #
    nn <- length(rho);
    res <- matrix(1,nn+1,nn+1);
    for (ii in r.bc(nn)) {
      for (jj in r.bd(ii+1,nn+1)) {
        res[ii,jj] <- prod(rho[r.bd(ii,jj-1)]);
        res[jj,ii] <- res[ii,jj];
      }
    }
    #
  # returning
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
aux6 <- function(rho)
#TITLE returns a convenient codiagonal matrix from a vector
#DESCRIPTION 
# See code, quite straightforward\cr
# Not intended for the final user
#DETAILS
#KEYWORDS 
#INPUTS
#{rho} << numeric vector with all square
#         values less than one.>>
#[INPUTS]
#VALUE
# Resulting matrix of size \code{length(rho)+1}.
#EXAMPLE
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 12_02_02
#REVISED 12_02_02
#--------------------------------------------
{
  #
  nn <- length(rho)+1;
  rb2 <- (1-rho^2);
  res <- matrix(0,nn,nn);
  if (nn > 1) {
    # first row
    res[1,1:2] <- c(1,-rho[1])/rb2[1];
    # last row
    res[nn,(nn-1):nn] <- c(-rho[nn-1],1)/rb2[nn-1];
    # other rows
    for (ii in r.bd(2,nn-1)) {
      res[ii,(ii-1):(ii+1)] <- c(-rho[ii-1]/rb2[ii-1],
                                 (1-rho[ii-1]^2*rho[ii]^2)/rb2[ii-1]/rb2[ii],
                                 -rho[ii]/rb2[ii]);
    }
  }
  # returning
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
aux7 <- function(rho,end)
#TITLE returns a 3x3 matrix precisionfrom colliding
# nodes of chain
#DESCRIPTION 
# See code, quite straightforward\cr
# Not intended for the final user
#DETAILS
#KEYWORDS 
#INPUTS
#{rho} << numeric(2) vector providing the 
#         two concerned correlations.>>
#{end} << logical(2) vector indicating
#         if the first and second node
#         of the block are one end of the chain.>>
#[INPUTS]
#VALUE
# Resulting matrix of size \code{3x3}.
#EXAMPLE
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 12_02_02
#REVISED 12_02_02
#--------------------------------------------
{
  #
  if (length(rho)!=2) { r.erreur(rho,message="must be of length 2");}
  if (sum(rho^2) >= 1) { r.erreur(rho,message="degenerate case");}
  #
  res <- matrix(1,3,3);
  res[2,1] <- res[1,2] <- -rho[1];
  res[3,2] <- res[2,3] <- -rho[2];
  res[3,1] <- res[1,3] <- rho[1]*rho[2];
  if (end[1]) {
    res[1,1] <- (1-rho[2]^2);
  } else {
    res[1,1] <- (1-rho[1]^4-rho[2]^2) / aux0(rho[1])^2;
  }
  if (end[2]) {
    res[3,3] <- (1-rho[1]^2);
  } else {
    res[3,3] <- (1-rho[2]^4-rho[1]^2) / aux0(rho[2])^2;
  }
  res <- res / aux0(rho)^2;
  # returning
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
aux8 <- function(rho)
#TITLE returns the 'conditioned' product block from a vector
#DESCRIPTION 
# See code, quite straightforward\cr
# Not intended for the final user
#DETAILS
# in fact, is equivalent to take the result of
# the function \code{aux5} and applied the
# conditional variance for the last.
#KEYWORDS 
#INPUTS
#{rho} << numeric vector>>
#[INPUTS]
#VALUE
# Resulting matrix of partial correlations (if \code{rho} is)
# of size \code{length(rho)}.
#EXAMPLE
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 12_02_07
#REVISED 12_02_07
#--------------------------------------------
{
    #
    nn <- length(rho);
    res <- matrix(1,nn,nn);
    for (ii in r.bc(nn)) {
      kon <- 1-prod(rho[r.bd(ii,nn)]^2);
      res[ii,ii] <- kon;
      for (jj in r.bd(1,ii-1)) {
        res[ii,jj] <- kon*prod(rho[r.bd(jj,ii-1)]);
        res[jj,ii] <- res[ii,jj];
      }
    }
    #
  # returning
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
