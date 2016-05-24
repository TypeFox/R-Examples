
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
order4nbn <- function(nbn,ord=NULL)
#TITLE topological order of a /nbn/
#DESCRIPTION returns one of the orders of the nodes
# such as the parents of any node are less ranked than
# it when \code{is.null(ord)}. If not check that the proposed order
# is either a right permutation (\code{is.numeric(ord)}) or
# a vector of node names providing a topological order
# (\code{is.character(ord)}).
#DETAILS
# When \code{!is.null(ord)} the order must be an order, if
# not an error is issued.
#KEYWORDS 
#INPUTS
#{nbn} <<\code{nbn} object for which the order must be computed.>>
#[INPUTS]
#{ord} <<\code{NULL} or an order to test as a permutation or
#        a vector of names.>>
#VALUE
# a permutation vector of the nodes of the /nbn/
#        or a named list with the nodes not having
#        their parents before them.
#EXAMPLE
#  names(rbmn0nbn.04)[order4nbn(rbmn0nbn.04)];
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE 
#AUTHOR J.-B. Denis
#CREATED 12_01_14
#REVISED 13_07_05
#--------------------------------------------
{
  # basic constants
  nn <- length(nbn);
  nam <- names(nbn);
  # checking
  for (ii in r.bc(nn)) {
    if (length(union(nam,nbn[[ii]]$parents))!= nn) {
      r.erreur(list(nam,nbn[[ii]]$parents),
             message="Unknown Parents");
    }
  }
  # processing
  if (is.null(ord)) {
    # getting a topological order
    res <- numeric(0);
    while (length(res) < nn) {
      for (ii in r.bc(nn)) {
        if (!(ii %in% res)) {
          pare <- nbn[[ii]]$parents;
          # if its parents are already included, it can be added
          if (length(union(nam[res],pare))==length(res)) {
            res <- c(res,ii);
          }
        }
      }
    }
  } else {
    # checking the proposed order
    # translating 'ord' as numeric if character
    if (is.character(ord)) { nord <- r.numero(ord,nam);
    } else { nord <- ord;}
    # checking it is a permutation
    if (!all(sort(nord)==1:nn)) {
        r.erreur(ord,message="'ord' is not a permutation.")
    }
    # checking the order
    res <- TRUE; 
    parents <- character(0);
    for (ii in r.bc(nn)) {
      pare <- nbn[[nord[ii]]]$parents;
      if (length(union(pare,parents))>length(parents)) {
        res <- FALSE;
      }
    parents <- union(parents,nam[ii]);
    }
  }
  # returning
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
order4gema <- function(gema,ord=NULL)
#TITLE topological order of a /gema/
#DESCRIPTION returns one of the orders of the nodes
# such as the parents of any node are less ranked than
# it when \code{is.null(ord)}. If not check that the proposed order
# is either a right permutation (\code{is.numeric(ord)}) or
# a vector of node names providing a topological order
# (\code{is.character(ord)}).
#DETAILS
# When \code{!is.null(ord)} the order must be an order, if
# not an error is issued.
#KEYWORDS 
#INPUTS
#{gema} <<\code{gema} object for which the order must be computed.>>
#[INPUTS]
#{ord} <<\code{NULL} or an order to test as a permutation or
#        a vector of names.>>
#VALUE
# a permutation vector of the nodes of the /gema/
#        or a named list with the nodes not having
#        their parents before them. That is a topological order.
#EXAMPLE
#  names(rbmn0gema.04$mu)[order4gema(rbmn0gema.04)];
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE 'ord' option is bad
#AUTHOR J.-B. Denis
#CREATED 12_01_14
#REVISED 12_01_17
#--------------------------------------------
{
  # checking
  # < to be done >
  eps <- 10^-10;
  li <- gema$li;
  # number of nodes
  nn <- nrow(li);
  nam <- dimnames(li)[[1]];
  if (is.null(ord)) {
    # getting the order
    res <- numeric(0);
    while (length(res) < nn) {
      for (ii in r.bc(nn)) {
        if (!(ii %in% res)) {
          nbpar <- sum(abs(li[ii,])>eps) - 1;
          if (nbpar==0) {
            res <- c(res,ii);
            li[,ii] <- 0;
          }
        }
      }
    }
  } else {
    return("Sorry, the 'ord' option is wrong and must be corrected!");
    # checking the proposed permutation
    if (is.numeric(ord)) {
      if (length(union(ord,r.bc(nn)))!=nn) {
        r.erreur(ord,message="'ord' is not a permutation.")
      }
      ord <- nam[ord];
    } else {
      if (length(union(ord,nam))!=nn) {
        r.erreur(list(ord,nam),
               message="'ord' is not a permutation of 'nam'");
      }
    }
    # checking the order
    res <- vector("list",0); kk <- 0; nm <- character(0);
    parents <- character(0);
    for (ii in r.bc(nn)) {
      parents <- c(parents,nam[ii]);
      wpar <- which(abs(gema$li[ii,])>eps);
      pare <- nam[wpar];
      if (length(union(pare,parents))>length(parents)) {
        kk <- kk + 1;
        res[[kk]] <- pare;
        nm <- c(nm,nam[ii]);
        names(res) <- nm;
      }
    }
  }
  # returning
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
normalize8nbn <- function(nbn,mu=0,sigma=1)
#TITLE normalizes a /nbn/
#DESCRIPTION returns a \code{nbn} with a
# given expectation and variance through
# an transformation leaving the correlation
# unchanged.
#DETAILS
#KEYWORDS
#INPUTS
#{nbn}<< The \code{nbn} object to transform.>>
#[INPUTS]
#{mu} << Imposed expectations. When \code{NULL}
#        nothing is changed. When of length one,
#        this value is given to all the node
#        expectations. If not the complete vector
#        of expectations.
#{sigma} << The same as \code{mu} but for the
#        standard deviations.>>
#VALUE
# The transformed \code{nbn}.
#EXAMPLE
# print8nbn(normalize8nbn(rbmn0nbn.01));
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 11_12_19
#REVISED 11_12_19
#--------------------------------------------
{
  # number of nodes
  nn <- length(nbn);
  nna <- names(nbn);
  # developping the mu and sigma
  if (length(mu)==1) {
    mu <- rep(mu,nn);
  }
  if (length(sigma)==1) {
    sigma <- rep(sigma,nn);
  }
  # going to the gema form
  gema <- nbn2gema(nbn);
  # getting the mn form
  mn <- gema2mn(gema);  
  # computing the actual mu and sigma
  mu0 <- mn$mu;
  sigma0 <- sqrt(diag(mn$gamma));
  # transforming accordingly the gema
  if (!is.null(mu)) {
    gema$mu <- mu;
  }
  if (!is.null(sigma)) {
    if (!all(sigma>0)) {
      r.erreur(sigma,message="asked standard deviation(s) are not positive");
    }
    gema$li <- diag(sigma/sigma0,nrow=nn) %*% gema$li;
  }
  names(gema$mu) <- nna;
  dimnames(gema$li) <- list(nna,NULL);
  # back to the nbn
  res <- gema2nbn(gema);
  # returning
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
nbn4nbn <- function(nbn)
#TITLE From a /nbn/ computes the associated nbn1
#DESCRIPTION returns a /nbn/ object with the 
# same structure as \code{nbn} but all \code{$mu}
# are put to zero, all \code{$sigma} to one as well
# as \code{$regcof}.
#DETAILS
# These coefficient values allows the easy study of the   
# /nbn/ structure.
#KEYWORDS
#INPUTS
#{nbn}<< The \code{nbn} object to transform.>>
#[INPUTS]
#VALUE
# The resulting \code{nbn}.
#EXAMPLE
# print8nbn(nbn4nbn(rbmn0nbn.04));
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE To be integrated into 'normalize8nbn' since it is a kind of normalization
#AUTHOR J.-B. Denis
#CREATED 12_02_17
#REVISED 12_02_17
#--------------------------------------------
{
  # checking
  # <to be done>
  # modification
  for (no in r.bf(nbn)) {
    nbn[[no]]$mu <- 0;
    nbn[[no]]$sigma <- 1;
    nbn[[no]]$regcoef[] <- 1;
  }
  # returning
  nbn;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
adja4nbn <- function(nbn) 
#TITLE adjacency matrix of a /nbn/
#DESCRIPTION returns a dimnamed matrix
# indicating with 1 an arc from row to column nodes
# (0 everywhere else); i.e. the adjacency matrix.
#DETAILS
#KEYWORDS 
#INPUTS
#{nbn}<< The initial \code{nbn} object.>>
#[INPUTS]
#VALUE
# A dimnamed matrix
#EXAMPLE
# adja4nbn(rbmn0nbn.04);
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 12_07_27
#REVISED 12_07_27
#--------------------------------------------
{
  # checking
  # To be done
  # getting the parentship matrix
  nbno <- length(nbn);
  nbna <- names(nbn);
  res <- matrix(0,nbno,nbno);
  dimnames(res) <- list(from=nbna,to=nbna);
  for (nn in r.bf(nbn)) {
    if (length(nbn[[nn]]$parents) > 0) {
      res[match(nbn[[nn]]$parents,nbna),nn] <- 1;
    }
  }
  # returning
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
adja2nbn <- function(adja) 
#TITLE standardized /nbn/ from an adjacency matrix
#DESCRIPTION returns a \code{nbn} object
# with O/1 regression coefficients having \code{adja} as
# adjacency matrix.
#DETAILS
#KEYWORDS 
#INPUTS
#{adja}<< The initial adjacency matrix.>>
#[INPUTS]
#VALUE
# The corresponding standardized \code{nbn} object.
#EXAMPLE
# print8nbn(adja2nbn(adja4nbn(rbmn0nbn.03)));
#REFERENCE
#SEE ALSO adja4nbn
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 13_05_02
#REVISED 13_05_02
#--------------------------------------------
{
  # checking
  # To be done
  # getting the parentship matrix
  nbno <- nrow(adja);
  nbna <- dimnames(adja)[[1]];
  res <- vector("list",nbno);
  names(res) <- nbna;
  for (nn in r.bf(res)) {
    res[[nn]]$mu <- 0;
    res[[nn]]$sigma <- 1;
    res[[nn]]$parents <- nbna[which(adja[,nn]==1)];
    res[[nn]]$regcoef <- rep(1,length(res[[nn]]$parents));
  }
  # returning
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
rm8nd4adja <- function(adja,nodes)
#TITLE  removes somes nodes from an adjacency matrix
#DESCRIPTION 
# Eliminates from the adjacency matrix (\code{adja})all \code{nodes}
# not breaking the existing links.\cr
# Important: the node order in \code{adja} must be topological.
#DETAILS
# When a node is removed, all its parents become parent of its children.
#KEYWORDS utilities
#PKEYWORDS 
#INPUTS
#{adja} <<The relation matrix to be consider (same format as those
#        provided by the function \code{adja4nbn}. Must be in topological
#        order, roots first.>>
#{nodes} <<Numeric or character vector providing the node numbers
#          to use for the generation of the subset.>>
#[INPUTS]
#VALUE
# The reduced adjacency matrix.
#EXAMPLE
# rm8nd4adja(rbmn0adja.04,"1.1");
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE remove the topological order constraint
#AUTHOR J.-B. Denis
#CREATED 12_07_27
#REVISED 12_07_27
#--------------------------------------------
{
# checking
nbno <- nrow(adja);
# the topological order
masque <- outer(r.bc(nbno),r.bc(nbno),">=");
if (sum(masque*adja^2) > 0) {
  r.erreur(adja,message="Not a topological order");
}
# < to be finished >
# getting constants
nano <- dimnames(adja)[[1]];
if (is.character(nodes)) {
  nodes <- match(nodes,nano);
}
nodes <- sort(unique(nodes),decreasing=TRUE);
nbs <- length(nodes);
# degenerate case
if (nbno*nbs == 0) {
  return(matrix(0,0,0));
}
# removing starting from the leaves
for (nod in nodes) {
  nbc <- nrow(adja);
  des <- r.bc(nbno)[adja[nod,]!=0];
  asc <- r.bc(nbno)[adja[,nod]!=0];
  if (length(des)*length(asc)>0) {
    adja[asc,des] <- 1;
  }
  adja <- adja[-nod,]; adja <- adja[,-nod];
}
# returning
adja;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
rm8nd4nbn <- function(nbn,nodes)
#TITLE removes some nodes from a /nbn/
#DESCRIPTION returns a /nbn/ object deduced from
# an original /nbn/ by integrating on a given
# subset of nodes.
#DETAILS
# The transformation is made through the associated
# joint distributions for the probabilities and with
# the help of the function \code{rm8nd4adja} 
# for the relationships.
#KEYWORDS
#INPUTS
#{nbn}<< The \code{nbn} object to reduce.>>
#{nodes} <<\code{character} or \code{numeric} vector
#           giving the subset of nodes to remove.>>
#[INPUTS]
#VALUE
# The resulting \code{nbn}.
#EXAMPLE
# rm8nd4nbn(rbmn0nbn.04,"1.1"); 
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE to be done properly!
#AUTHOR J.-B. Denis
#CREATED 12_05_16
#REVISED 12_07_27
#--------------------------------------------
{
  # checking
  return("the presently proposed algorithm is false, see 'margina.r' to check it");
  # getting the ordering subset as a numeric
  if (is.numeric(nodes)) {
    ss <- sort(nodes);
    if ((max(ss) > length(nbn)) | (min(ss) < 1)) {
      r.erreur(list(length(nbn),nodes),"Not a subset (1)");
    }
  } else {
    ss <- sort(match(nodes,names(nbn)));
    if (length(ss)!=length(nodes)) {
      r.erreur(list(names(nbn),nodes),"Not a subset (2)");
    }
  }
  nbno <- length(nbn);
  # getting the structure
  pam <- adja4nbn(nbn);
  pam <- rm8nd4adja(pam,ss);
  # getting the parameters
  rr <- sort(setdiff(r.bc(nbno),ss));
  para <- gema2mn(nbn2gema(nbn));
  para <- list(mu=para$mu[rr],
               gamma=para$gamma[rr,rr,drop=FALSE]
              );
  # reconsistuting the nbn from pam and para
  res <- vector("list",length(rr));
  names(res) <- names(nbn)[rr];
  for (nn in r.bf(rr)) {
    # the parents
    rrpa <- which(pam[,nn]==1);
    res[[nn]]$parents <- names(res)[rrpa];
    # the parameters values
    if (length(rrpa)==0) {
      # no parents
      mumu <- para$mu[nn];
      gaga <- para$gamma[nn,nn];
      coco <- numeric(0);
    } else {
      # parents
      gvv <- solve(para$gamma[rrpa,rrpa,drop=FALSE]);
      guv <- para$gamma[nn,rrpa,drop=FALSE];
      mup <- matrix(para$mu[rrpa],length(rrpa));
      mumu <- para$mu[nn] - guv %*% gvv %*% mup;
      gaga <- guv %*% gvv %*% t(guv);
      coco <- guv %*% gvv;
    }
    res[[nn]]$mu <- as.vector(mumu);
    res[[nn]]$sigma <- as.vector(sqrt(gaga));
    res[[nn]]$regcoef <- as.vector(coco);
  }
  # returning
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
generate8nbn <- function(rnn=c(3,7),ppar=0.5,rreg=c(-1,1),
                         rmu=c(0,0),rsig=c(0,1),
                         nona=r.form3names(max(rnn)))
#TITLE returns a randomly built /nbn/ object.
#DESCRIPTION 
# To obtain systematic results, you have to call \code{set.seed}
# before hands.
#DETAILS
# Node numbers are uniformly drawn. Parent numbers are
# independently drawn from all ancestors with the probability
# associated to the considered node. Regression coefficient are uniformly drawn.
# Conditional expectations and standard deviations are uniformly drawn.\cr
# All range arguments can be given one value instead of two, to precise the
# unique value to use.
#KEYWORDS 
#INPUTS
#[INPUTS]
#{rnn} <<Range of the number of nodes.>>
#{ppar} <<Probabilities (not a range) of the parent occurrence for
# each ancestor of every node. Can be a vector,
# cycled as necessary.>>
#{rreg} <<Range of regression coefficients.>>
#{rmu}  <<Range of the conditional expectations.>>
#{rsig} <<Range of the conditional standard deviations.>>
#{nona} <<Proposed names for the maximum number of nodes, only the
#         necessary first ones will be used.>>
#VALUE
# a /nbn/ object, with nodes in topological order.
#EXAMPLE
# set.seed(1234)
# print8nbn(generate8nbn());
# print8nbn(generate8nbn());
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 13_04_21
#REVISED 13_04_24
#--------------------------------------------
{
  # checking
  # < to be done >
  # dealing with unique values
  if (length(rnn)==1) { rnn <- c(rnn,rnn);}
  if (length(rreg)==1) { rreg <- c(rreg,rreg);}
  if (length(rmu)==1) { rmu <- c(rmu,rmu);}
  if (length(rsig)==1) { rsig <- c(rsig,rsig);}
  # getting the node number.
  nn <- floor(runif(1,rnn[1],rnn[2]+1));
  # getting the node names
  nona <- nona[r.bc(nn)];
  # getting the probabilities of parents
  pp <- ppar;
  while (length(pp) < nn) { pp <- c(pp,ppar);}
  pp <- pp[1:nn];
  # building the /nbn/
  res <- vector("list",0);
  for (ii in r.bc(nn)) {
    mu <- runif(1,rmu[1],rmu[2]);
    sigma <- runif(1,rsig[1],rsig[2]);
    regcoef <- numeric(0);
    parents <- character(0);
    for (jj in r.bc(ii-1)) {
      if (rbinom(1,1,pp[ii]) > 0) {
        parents <- c(parents,nona[jj]);
        regcoef <- c(regcoef,runif(1,rreg[1],rreg[2]));
      }
    }
    res[[ii]] <- list(mu=mu,sigma=sigma,
                      parents=parents,
                      regcoef=regcoef);
  }
  names(res) <- nona;
  # returning
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
arc7nb4nbn <- function(nbn,each=FALSE) 
#TITLE returns the number(s) of arcs of a /nbn/
#DESCRIPTION returns the arc numbers of the node
# of /nbn/ object.
#DETAILS
# Parents associated with a zero regression coefficient
# are not excluded in the counting.
#KEYWORDS 
#INPUTS
#{nbn}<< The \code{nbn} object to consider.>>
#[INPUTS]
#{each} << When \code{TRUE}, returns a named vector of the
# number of parents of each node. If not the total number of arcs.>>
#VALUE
# Either a number or a named vector of numbers (names being the node names).
#EXAMPLE
# arc7nb4nbn(rbmn0nbn.05);
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE To complete with a broad description of properties
#       returning a list and associating a printing function.
#AUTHOR J.-B. Denis
#CREATED 13_04_22
#REVISED 13_04_22
#--------------------------------------------
{
  # checking
  # To be done
  # getting the arc number for each node
  nbno <- length(nbn);
  res <- rep(0,nbno);
  for (nn in r.bf(nbn)) {
    res[nn] <- length(nbn[[nn]]$parents);
  }
  # finalizing
  if (each) {
    names(res) <- names(nbn);
  } else {
    res <- sum(res);
  }
  # returning
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
