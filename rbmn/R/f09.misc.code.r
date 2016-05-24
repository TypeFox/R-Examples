
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
nb8bn <- function(n,label=FALSE)
#TITLE number of Bayesian networks
#DESCRIPTION
# returns the number of different Bayesian networks
# having \code{n} labelled or not nodes. Non labelled nodes means
# that nodes are exchangeable: \code{A -> B} is identical to
# \code{A <- B}.
#DETAILS
# When not labelled nodes, the results were proposed by Sloane in 'the on line
# encyclopedy of integer sequences' (http://oeis.org/A003087).
# For labelled nodes, just the application of the recursive formula of Robinson. 
#KEYWORDS 
#INPUTS
#{n} << number of nodes. Must be less or equal to 18.>>
#[INPUTS]
#{label} <<Indicates if the nodes must be considered as labelled or not.>>
#VALUE
# Number of Bayesian networks
#EXAMPLE
# nb8bn(5)
# nb8bn(5,TRUE);
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 12_06_13
#REVISED 13_04_29
#--------------------------------------------
{
  if (label) {
    n <- max(0,round(n));
    if (n > 18) {
      r.erreur(n,"'n' is too large",w=TRUE);
      return(NA);
    }
    res <- vector("numeric",19);
    res[1+ 0] <- 1;
    res[1+ 1] <- 1;
    res[1+ 2] <- 2;
    res[1+ 3] <- 6;
    res[1+ 4] <- 31;
    res[1+ 5] <- 302;
    res[1+ 6] <- 5984;
    res[1+ 7] <- 243668;
    res[1+ 8] <- 20286025;
    res[1+ 9] <- 3424938010;
    res[1+10] <- 1165948612902;
    res[1+11] <- 797561675349580;
    res[1+12] <- 1094026876269892596;
    res[1+13] <- 3005847365735456265830;
    res[1+14] <- 16530851611091131512031070;
    res[1+15] <- 181908117707763484218885361402;
    res[1+16] <- 4004495398476191849391903634065582;
    res[1+17] <- 176332675845335018307024273011267894000;
    res[1+18] <- 15530301094140830400618221389766986731287870;
    #
    res <- res[1+n];
  } else {
    n <- round(n);
    if ((n<0) | (n>20)) { r.erreur(n,"'n' must not be negative or too large");}
    #
    if (n == 0) {
      # initialization
      res <- 1;
    } else {
      # recurrence formula
      res <- 0;
      for (k in 1:n) {
        res <- res + (-1)^(k-1)*choose(n,k)*2^(k*(n-k))*Recall(n-k);
      }
    }
  }
  # returning
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
adja4three <- function(nona=LETTERS[1:3])
#TITLE Adjacency matrices of DAGs having three nodes
#DESCRIPTION
# Returns the list of the 25 adjacency matrices associated
# to DAGs comprising three nodes. The first character
# of the name components gives the number of arcs in the DAG.
#DETAILS
# Poor filling... 
#KEYWORDS 
#INPUTS
#[INPUTS]
#{nona} << The three node names.>>
#VALUE
# a named list having 25 components, each being a 3x3 matrix.
#EXAMPLE
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
  # initializing
  se <- "-";
  res <- vector("list",25);
  ama <- matrix(0,3,3,dimnames=list(nona,nona));
  for (ii in 1:25) { res[[ii]] <- ama; }
  kk <- 0;
  # filling the zero arc DAG
  kk <- kk + 1; names(res)[1] <- paste(0,1,sep=se);
  # filling the one arc DAGS
  aa <- matrix(c(1,2,2,1,1,3,3,1,2,3,3,2),
               ncol=2,byrow=TRUE);
  for (ii in r.bc(nrow(aa))) {
    kk <- kk+1;
    res[[kk]][aa[ii,1],aa[ii,2]] <- 1;
    names(res)[[kk]] <- paste(1,ii,sep=se);
  }
  # filling the two arcs DAGS
  aa <- matrix(c(1,2,1,3,
                 2,1,1,3,
                 1,2,3,1,
                 2,1,3,1,
                 2,1,2,3,
                 2,1,3,2,
                 1,2,2,3,
                 1,2,3,2,
                 3,1,3,2,
                 3,1,2,3,
                 1,3,3,2,
                 1,3,2,3
                 ),
               ncol=4,byrow=TRUE);
  for (ii in r.bc(nrow(aa))) {
    kk <- kk+1;
    res[[kk]][aa[ii,1],aa[ii,2]] <- 1;
    res[[kk]][aa[ii,3],aa[ii,4]] <- 1;
    names(res)[[kk]] <- paste(2,ii,sep=se);
  }
  # filling the three arcs DAGS
  aa <- matrix(c(1,2,1,3,2,3,
                 1,2,1,3,3,2,
                 2,1,2,3,1,3,
                 2,1,2,3,3,1,
                 3,1,3,2,1,2,
                 3,1,3,2,2,1
                 ),
               ncol=6,byrow=TRUE);
  for (ii in r.bc(nrow(aa))) {
    kk <- kk+1;
    res[[kk]][aa[ii,1],aa[ii,2]] <- 1;
    res[[kk]][aa[ii,3],aa[ii,4]] <- 1;
    res[[kk]][aa[ii,5],aa[ii,6]] <- 1;
    names(res)[[kk]] <- paste(3,ii,sep=se);
  }
  # returning
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
adja2arcs <- function(adj)
#TITLE Arc matrix from an adjacency matrix
#DESCRIPTION
# returns the arc matrix from an adjacency matrix.
#DETAILS
#KEYWORDS 
#INPUTS
#{adj} << The adjacency matrix.>>
#[INPUTS]
#VALUE
# a matrix with two columns ("from","to")
#EXAMPLE
# adja2arcs(rbmn0adja.02)
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
  # initializing
  nn <- nrow(adj); nona <- dimnames(adj)[[1]];
  res <- matrix(NA,0,2,dimnames=list(NULL,c("to","from")));
  # filling
  for (ii in r.bc(nn)) {
    ou <- which(adj[ii,]==1);
    res <- rbind(res,cbind(rep(nona[ii],length(ou)),nona[ou]));
  }
  # returning
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
