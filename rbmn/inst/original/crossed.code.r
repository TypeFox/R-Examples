
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
crossed4nbn1nbn <- function(nbn1,nbn2,
                        we1=rep(1,length(nbn1)),
                        we2=rep(1,length(nbn2)),
                        nona=as.vector(outer(names(nbn1),
                                             names(nbn2),paste,sep="_")))
#TITLE creates a crossed-nbn from two /nbn/s
#DESCRIPTION
# A crossed /nbn/ is a /nbn/ obtained when replacing
# each node of the first /nbn/ by the second /nbn/ and
# vice-versa.\cr
# Let \code{nn1/nn2} and \code{na1/na2} be the node and arc
# numbers of the two \code{nbn}s, the node number of the
# crossed \code{nbn} is \code{nn1*nn2} and its arc number
# is \code{nn1*na2+nn2*na1}.\cr
# The regression coefficients attributed to the crossed \code{nbn}
# are the products of the weights (\code{we1/we2}) and the regression
# coefficients of the initial \code{nbn}.
#DETAILS
# The \code{mu} coefficient is the sum of the two corresponding \code{mu}s
# of the generating \code{nbn}.\cr
# The \code{sigma} coefficient is the product of the two corresponding \code{sigma}s
# of the generating \code{nbn}.\cr
# The regression coefficient are directed inherited from the \code{nbn}
# which is duplicated with this arc.
#KEYWORDS 
#INPUTS
#{nbn1}<< The first generating /nbn/.>>
#{nbn2}<< The second generating /nbn/.>>
#[INPUTS]
#{we1}<< The weight to apply to the nodes of the first generating /nbn/.>>
#{we2}<< The weight to apply to the nodes of the second generating /nbn/.>>
#{nona}<< The node names to give to the crossed /nbn/, the nodes of the
# \code{nbn1} varying first.>>
#VALUE
# The resulting crossed \code{nbn} object.
#EXAMPLE
# print8nbn(crossed4nbn1nbn(rbmn0nbn.01,rbmn0nbn.04));
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 13_04_22
#REVISED 13_04_24
#--------------------------------------------
{
  # checking
  # < to be done >
  # constants
  nn1 <- length(nbn1); nn2 <- length(nbn2);
  if (length(nona) != nn1*nn2) {
    r.erreur(list(nn1,nn2,length(nona)),
           message="Bad correspondence with the proposed names"
          );
  }
  nuij <- function(ii,jj,nlig) {ii + nlig*(jj-1);}
  # initializing
  res <- vector("list",nn1*nn2);
  names(res) <- nona;
  # making each crossed node in turn
  for (ii in r.bc(nn1)) {
    ip1 <- match(nbn1[[ii]]$parents,names(nbn1));
    for (jj in r.bc(nn2)) {
      kk <- nuij(ii,jj,nn1);
      #cat(names(nbn1)[ii],names(nbn2)[jj],":",nona[kk],"\n");
      res[[kk]]$mu    <- nbn1[[ii]]$mu    + nbn2[[jj]]$mu;
      res[[kk]]$sigma <- nbn1[[ii]]$sigma * nbn2[[jj]]$sigma;
      ip2 <- match(nbn2[[jj]]$parents,names(nbn2));
      res[[kk]]$parents <- character(0);
      res[[kk]]$regcoef <- numeric(0);
      for (ij in r.bf(ip1)) {
        nupa <- nona[nuij(ip1[ij],jj,nn1)];
        res[[kk]]$parents <- c(res[[kk]]$parents,nupa);
        res[[kk]]$regcoef <- c(res[[kk]]$regcoef,we2[jj]*nbn1[[ii]]$regcoef[ij]);
      }
      for (ij in r.bf(ip2)) {
        nupa <- nona[nuij(ii,ip2[ij],nn1)];
        #cat(names(nbn1)[ii],names(nbn2)[ij],":",nupa,"\n");
        res[[kk]]$parents <- c(res[[kk]]$parents,nupa);
        res[[kk]]$regcoef <- c(res[[kk]]$regcoef,we1[ii]*nbn2[[jj]]$regcoef[ij]);
      }
      #pause(paste("<<",nona[kk],">>"));
    }
  }
  # returning
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
arcs4nbn1nbn <- function(nbn1,nbn2,type="a1",
                        nona=as.vector(outer(names(nbn1),
                                             names(nbn2),paste,sep="_")))
#TITLE returns the list of 'parallel' arcs of a crossed-nbn
#DESCRIPTION
# Returns a list of matrices with two columns (as needed by \code{estimate8constrainednbn})
# indicating corresponding arcs for each arcs/nodes of \code{nbn1} (or \code{nbn2}) of the 
# crossed /nbn/ obtained when crossing /nbn1/ and /nbn2/ with node names given by \code{nona}.
#DETAILS
#KEYWORDS 
#INPUTS
#{nbn1}<< The first generating /nbn/.>>
#{nbn2}<< The second generating /nbn/.>>
#[INPUTS]
#{type}<< Must be \code{"a1"} to indicate that the parallelism must be done for each arc of \code{nbn1}.
#         \code{"a2"} for each arc of \code{nbn2}.
#         Or \code{"n1"} for each node of \code{nbn1}.
#         Or \code{"n2"} for each node of \code{nbn2}.
#>>
#{nona}<< The node names to give to the crossed /nbn/, the nodes of the
# \code{nbn1} varying first.>>
#VALUE
# The resulting named (after node names) list of matrices.
#EXAMPLE
#  print(arcs4nbn1nbn(rbmn0nbn.01,rbmn0nbn.04));
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 13_04_26
#REVISED 13_04_26
#--------------------------------------------
{
  # checking
  # < to be done >
  # constants
  nn1 <- length(nbn1); nn2 <- length(nbn2);
  if (length(nona) != nn1*nn2) {
    r.erreur(list(nn1,nn2,length(nona)),
           message="Bad correspondence with the proposed names"
          );
  }
  nuij <- function(ii,jj,nlig) {ii + nlig*(jj-1);}
  fait <- FALSE;
  # 
  if (type=="n1") {
    # initializing
    res <- vector("list",nn1);
    names(res) <- names(nbn1);
    # making each crossed node in turn
    for (ii in r.bc(nn1)) {
      res[[ii]] <- matrix(NA,0,2);
      for (jj in r.bc(nn2)) {
        head <- nona[nuij(ii,jj,nn1)];
        ip2 <- match(nbn2[[jj]]$parents,names(nbn2));
        for (ij in r.bf(ip2)) {
          tail <- nona[nuij(ii,ip2[ij],nn1)];
          res[[ii]] <- rbind(res[[ii]],c(tail,head));
        }
      }
    }
    fait <- TRUE;
  }
  #
  if (type=="n2") {
    # initializing
    res <- vector("list",nn2);
    names(res) <- names(nbn2);
    # making each crossed node in turn
    for (ii in r.bc(nn2)) {
      res[[ii]] <- matrix(NA,0,2);
      for (jj in r.bc(nn1)) {
        head <- nona[nuij(jj,ii,nn1)];
        ip1 <- match(nbn1[[jj]]$parents,names(nbn1));
        for (ij in r.bf(ip1)) {
          tail <- nona[nuij(ip1[ij],ii,nn1)];
          res[[ii]] <- rbind(res[[ii]],c(tail,head));
        }
      }
    }
    fait <- TRUE;
  }
  #
  if (type=="a1") {
    # getting the arcs of nbn1
    arcs1 <- adja2arcs(adja4nbn(nbn1));
    # initializing
    res <- vector("list",nrow(arcs1));
    names(res) <- paste(arcs1[,1],arcs1[,2],sep="->");
    # making each crossed node in turn
    for (ii in r.bf(res)) {
      res[[ii]] <- matrix(NA,0,2);
      ii1 <- match(arcs1[ii,1],names(nbn1));
      ii2 <- match(arcs1[ii,2],names(nbn1));
      for (jj in r.bc(nn2)) {
        headn <- nona[nuij(ii1,jj,nn1)];
        tailn <- nona[nuij(ii2,jj,nn1)];
        res[[ii]] <- rbind(res[[ii]],c(headn,tailn));
      }
    }
    fait <- TRUE;
  }
  # 
  if (type=="a2") {
    # getting the arcs of nbn2
    arcs2 <- adja2arcs(adja4nbn(nbn2));
    # initializing
    res <- vector("list",nrow(arcs2));
    names(res) <- paste(arcs2[,1],arcs2[,2],sep="->");
    # making each crossed node in turn
    for (jj in r.bf(res)) {
      res[[jj]] <- matrix(NA,0,2);
      jj1 <- match(arcs2[jj,1],names(nbn2));
      jj2 <- match(arcs2[jj,2],names(nbn2));
      for (ii in r.bc(nn1)) {
        headn <- nona[nuij(ii,jj1,nn1)];
        tailn <- nona[nuij(ii,jj2,nn1)];
        res[[jj]] <- rbind(res[[jj]],c(headn,tailn));
      }
    }
    fait <- TRUE;
  }
  if (!fait) {
    stop(paste("'type' =",type,"not recognized by 'arcs4nbn1nbn'"));
  }
  # returning
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
adja2crossed <- function(adj1,adj2,
                        nona=as.vector(outer(dimnames(adj1)[[1]],
                                             dimnames(adj2)[[1]],
                                             paste,sep="_")))
#TITLE creates a crossed-adjacency matrix from two ones
#DESCRIPTION
# Like crossed4nbn1nbn but at the level of adjacency matrices. Must 
# be much efficient when regression coefficients are not needed.
#DETAILS
# Just two Kronecker products of matrices.
#KEYWORDS 
#INPUTS
#{adj1}<< The first adjacency matrix.>>
#{adj2}<< The second adjacency matrix.>>
#[INPUTS]
#{nona}<< The node names to give to the crossed /nbn/, the nodes of the
# \code{nbn1} varying first.>>
#VALUE
# The resulting crossed adjacency matrix.
#EXAMPLE
#  print(adja2crossed(rbmn0adja.01,rbmn0adja.01));
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 13_04_30
#REVISED 13_04_30
#--------------------------------------------
{
  # checking
  # < to be done >
  # constants
  nn1 <- nrow(adj1); nn2 <- nrow(adj2);
  if (length(nona) != nn1*nn2) {
    r.erreur(list(nn1,nn2,length(nona)),
           message="Bad correspondence with the proposed names"
          );
  }
  # computation
  res <- diag(nrow=nn2) %x% adj1 + adj2 %x% diag(nrow=nn1);
  dimnames(res) <- list(nona,nona);
  # returning
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
