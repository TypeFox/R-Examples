
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
generate8chain <- function(rnn=c(3,7),proo=0.5,rcor=c(-1,1),
                         rmu=c(0,0),rsig=c(0,1),
                         nona=r.form3names(max(rnn)))
#TITLE generation of a /chain/ /nbn/
#DESCRIPTION [randomly] generates a /chain/ /nbn/.
#DETAILS
# Proposed ranges can be a unique value, implying no randomness
# in the value.\cr
# Roots are placed according to \code{proo} probabilities, then collider
# are placed in between with uniform probability on the possibles nodes.
#KEYWORDS 
#INPUTS
#[INPUTS]
#{rnn} <<Range of the number of nodes.>>
#{proo} <<Probabilit[y|ies] that the successive and acceptable nodes be colliders. Can be a vector.>>
#{rcor} <<Range of the correlations between neighbour nodes.>>
#{rmu}  <<Range of the expectations.>>
#{rsig} <<Range of the standard deviations.>>
#{nona} <<Proposed names for the maximum number of nodes, only the
#         necessary first ones will be used.>>
#VALUE
# A /chain/ coding list is returned.
#EXAMPLE
# set.seed(1234);
# print8chain(generate8chain());
# print8chain(generate8chain());
# print8chain(generate8chain(rnn=10,rcor=0.5));
# print8chain(generate8chain(rnn=10,rcor=0.5));
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 13_07_03
#REVISED 13_07_03
#--------------------------------------------
{
  # checking
  # < to be done >
  # dealing with unique values
  if (length(rnn)==1) { rnn <- c(rnn,rnn);}
  if (length(rcor)==1) { rcor <- c(rcor,rcor);}
  if (length(rmu)==1) { rmu <- c(rmu,rmu);}
  if (length(rsig)==1) { rsig <- c(rsig,rsig);}
  # getting the node number.
  nn <- floor(runif(1,rnn[1],rnn[2]+1));
  # getting the node names
  nona <- nona[r.bc(nn)];
  # getting the probabilities of downstream
  pp <- proo;
  while (length(pp) < nn) { pp <- c(pp,proo);}
  pp <- pp[1:nn];
  # building the /chain/
  # miscellaneous part
  res <- vector("list",0);
  res$names <- nona;
  res$mu <- runif(nn,rmu[1],rmu[2]);
  res$sigma <- runif(nn,rsig[1],rsig[2]);
  res$corre <- runif(nn-1,rcor[1],rcor[2]);
  # placing roots first
  root <- rbinom(nn,1,pp);
  for (ii in r.bc(nn-1)) {
    if (root[ii]==1) { root[ii+1] <- 0;}
  }
  # at least one root
  if (sum(root)==0) {
    root[sample.int(nn,1)] <- 1;
  }
  res$roots <- nona[root==1];
  root <- which(root==1);
  # placing colliders second
  res$colliders <- character(0);
  for (ii in r.bc(length(root)-1)) {
    r1 <- root[ii]; r2 <- root[ii+1];
    if (r2-r1<2) { stop("Erreur in generate8chain");}
    qui <- sample.int(r2-r1-1,1)+r1;
    res$colliders <- c(res$colliders,nona[qui]);
  }
  # returning
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
print8chain <- function(chain,digits=3)
#TITLE prints a /chain/ object
#DESCRIPTION prints a /chain/ object.
#DETAILS
# See \code{nbn2chain} code for some details about the
# definition of a /chain/.
#KEYWORDS 
#INPUTS
#{chain} << The \code{chain} object to print.>>  
#[INPUTS]
#{digits} << when not null, the number of digits for rounding the
# numerical values.>>
#VALUE
# nothing but something is printed
#EXAMPLE
# print8chain(rbmn0chain.01);
# print8chain(rbmn0chain.02);
# print8chain(rbmn0chain.03);
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 12_01_18
#REVISED 12_01_30
#--------------------------------------------
{
  # checking
  che <- check8chain(chain);
  if (is.character(che)) {
    print(che);
    stop("The provided 'chain' is not valid!");
  }
  # node number
  nn <- length(chain$names);
  orien <- aux2(chain)/2+1.5;
  cat("#-------------------------\n");
  for (node in r.bc(nn)) {
    nna <- chain$names[node];
    nona <- paste("(",nna,")",sep="");
    if (nna %in% chain$roots) {
      nona <- paste("<--(",nna,")-->",sep="");
    }
    if (nna %in% chain$colliders) {
      nona <- paste("-->)",nna,"(<--",sep="");
    }
    cat(r.form3justifie(nona,nbc=10,format=2));
    cat("  ",round(chain$mu[node],digits),"              (",
             round(chain$sigma[node],digits),")",sep="");
    cat("\n");
    if (node < nn) {
      cat(r.form3justifie(c("^","v")[orien[node]],nbc=10,format=2));
      cat("     <",round(chain$corre[node],digits),">",sep="");
      cat("\n");
    }
  }
  cat("#-------------------------\n");
  # returning
  invisible();
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
check8chain <- function(chain)
#TITLE checks a /chain/ object
#DESCRIPTION checks the consistency of \code{chain} as a /chain/ object
# issues a fatal error with some clues if inconsistent.
#DETAILS
# Looking a the code of this function provides a way to know which
# are the requirements of a /chain/ object.
#KEYWORDS 
#INPUTS
#{chain} << The \code{chain} object to check.>>  
#[INPUTS]
#VALUE
# \code{TRUE} or a \code{character} containing some clue
# about the discovered inconsistency.
#EXAMPLE
# check8chain(rbmn0chain.01);
# res <- check8chain(rbmn0adja.01);
# if (is.na(as.logical(res))) { print(res);}
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 13_07_26
#REVISED 13_07_26
#--------------------------------------------
{
  # checking the component names
  compo <- names(rbmn0chain.01);
  if (!setequal(compo,names(chain))) {
    return("The names of the chain are not 'names(rbmn0chain.01)'");
  }
  # checking the lengths
  nbno <- length(chain$names);
  res <- character(0);
  if (nbno != length(chain$mu)) { res <- c(res,"'$mu' has got a bad length");}
  if (nbno != length(chain$mu)) { res <- c(res,"'$sigma' has got a bad length");}
  if ((nbno-1) != length(chain$corre)) { res <- c(res,"'$corre' has got a bad length");}
  if ((nbno-1)/2 < length(chain$colliders)) { res <- c(res,"'$colliders' has got a bad length");}
  if ((nbno+1)/2 < length(chain$roots)) { res <- c(res,"'$roots' has got a bad length");}
  if (length(res) > 0) { return(res);}
  # checking the correlations
  if (any(abs(chain$corre)>1)) { return("Some correlation doesn't make sense");}
  # returning
  TRUE;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
order4chain <- function(chain,ord=NULL)
#TITLE returns a topological order of a /chain/ or
# checks a proposed order.
#DESCRIPTION From a \code{chain} object
# returns one of the possible topological orders,
# through a permutation when \code{is.null(ord)}.
# If not \code{ord} must be a proposed order to be
# checked given as a permutation if \code{is.numeric(ord)}
# or a vector of ordered names if \code{is.character(ord)}.
#DETAILS
# For the moment the \code{ord} option is
# bad and an error message is returned when used.
#KEYWORDS 
#INPUTS
#{chain}<< the \code{chain} object to be considered.>>
#[INPUTS]
#{ord} << Indicates what must be done. \code{NULL} to get a topological
# order associated to the chain otherwise a permutation to be checked as
# one of the possible topological orders of the chain.>>
#VALUE
# a permutation vector of the nodes of the /nbn/
#        or a named character with the nodes not having
#        their parents before them; when it is of
#        length zero this means that the check 
#        was successful.
#EXAMPLE
# order4chain(rbmn0chain.02);
# order4chain(rbmn0chain.02,order4chain(rbmn0chain.02));
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE Correct the 'ord' option!
#AUTHOR J.-B. Denis
#CREATED 12_01_31
#REVISED 12_01_31
#--------------------------------------------
{
  # checking
  # < to be done >
  # getting the permutation
  ooo <- aux2(chain,TRUE);
  if (is.null(ord)) {
    res <- ooo;
  } else {
    return("Sorry, the 'ord' option is wrong and must be corrected!");
    nam <- chain$names;
    nn <- length(nam);
    if (is.numeric(ord)) {
      # a numeric permutation is expected
      if (length(union(ord,r.bc(nn)))!=nn) {
        r.erreur(ord,message="'ord' is not a permutation.")
      }
    } else {
      # a list of names is expected
      if (length(union(ord,nam))!=nn) {
        r.erreur(list(ord,nam),
               message="'ord' is not a permutation of 'nam'");
      }
    ord <- r.numero(ord,nam);
    }
    # checking and storing inconsistencies
    res <- numeric(0); nm <- character(0);
    names(res) <- nm;
    etat <- state4chain(chain);
    for (uu in ord) {
      if (etat[uu] == "r") {
        if (uu > 1) {
	    if (etat[uu-1] == "t") { etat[uu-1] <- "r";}
	    if (etat[uu-1] == "c") { etat[uu-1] <- "t";}
        }
        if (uu < nn) {
	    if (etat[uu+1] == "t") { etat[uu+1] <- "r";}
	    if (etat[uu+1] == "c") { etat[uu+1] <- "t";}
        }
      } else {
	  res <- c(res,uu);
          nm <- c(nm,nam[uu]);
          names(res) <- nm;
      }
      etat[uu] <- "f";
    }
  }
  # returning
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
state4chain <- function(chain)
#TITLE returns the states of each node of a chain
#DESCRIPTION From a \code{chain} object
# returns a named character precising the role of each node:
# "r" for root, "c" for collider, "t" for transmitter and
# "l" for leaf.
#DETAILS
#KEYWORDS 
#INPUTS
#{chain}<< the \code{chain} object to be considered.>>
#[INPUTS]
#VALUE
# a character of the states named with node names.
#EXAMPLE
# state4chain(rbmn0chain.01);
# state4chain(rbmn0chain.03);
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 12_01_31
#REVISED 12_01_31
#--------------------------------------------
{
  # checking
  # < to be done >
  # identifying
  nn <- length(chain$names);
  etat <- rep("t",nn);
  etat[c(1,nn)] <- "l";
  if (length(chain$colliders) > 0) {
      etat[r.numero(chain$colliders,chain$names)] <- "c";
  }
  etat[r.numero(chain$roots,chain$names)] <- "r";
  # naming
  names(etat) <- chain$names;
  # returning
  etat;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
marginal4chain <- function(chain)
#TITLE returns marginal expectations and standard deviations of a chain
#DESCRIPTION From a \code{chain} object
# returns a list with two components: \code{$mu} and \code{$sigma}
# vectors of marginal expectations and standard deviations.\cr
#DETAILS
#KEYWORDS 
#INPUTS
#{chain}<< the \code{chain} object to be considered.>>
#[INPUTS]
#VALUE
# a list with the two components \code{$mu} and \code{$sigma}.
#EXAMPLE
# marginal4chain(rbmn0chain.02);
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 12_01_31
#REVISED 12_01_31
#--------------------------------------------
{
  # checking
  # < to be done >
  # getting a topological order
  # and the status of each node
    nn <- length(chain$names);
    ooo <- order4chain(chain);
    etat <- state4chain(chain);
  # following this order computing
  # progressively expectations and 
  # standard deviations
    mu <- sigma <- rep(NA,nn);
    orien <- aux2(chain);
    for (uu in ooo) {
	if (etat[uu]=="r") {
            # root, no modification
	    mu[uu] <- chain$mu[uu];
	    sigma[uu] <- chain$sigma[uu];
	} else {
          if (etat[uu]=="c") {
              # collider: the two neighbours are parents
              sigma[uu] <- chain$sigma[uu]/aux0(chain$corre[(uu-1):uu]);              
              mu[uu] <- chain$mu[uu] + 
                        chain$corre[uu-1]*sigma[uu]/sigma[uu-1]*mu[uu-1] +
                        chain$corre[uu]  *sigma[uu]/sigma[uu+1]*mu[uu+1];
	  } else {
              # other cases: only one parent
              if (uu>1) { if (orien[uu-1]==1) {
                  # parent upward
                  sigma[uu] <- chain$sigma[uu]/aux0(chain$corre[uu-1]);              
                  mu[uu] <- chain$mu[uu] + 
                            chain$corre[uu-1]*sigma[uu]/sigma[uu-1]*mu[uu-1];

	      }}
              if (uu<nn) { if (orien[uu]==-1) {
                  # parent downward
                  sigma[uu] <- chain$sigma[uu]/aux0(chain$corre[uu]);              
                  mu[uu] <- chain$mu[uu] + 
                            chain$corre[uu]*sigma[uu]/sigma[uu+1]*mu[uu+1];
	      }}
	  }
	}
    }
  if (any(is.na(mu))) {
      r.erreur(chain,message="either the chain or the algo");
  }
  # returning
  list(mu=mu,sigma=sigma);
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
chain2nbn <- function(chain)
#TITLE transforms a /chain/ to a /nbn/
#DESCRIPTION From a \code{chain} object
# returns the \code{nbn} translation.
#DETAILS
#KEYWORDS 
#INPUTS
#{chain}<< the \code{chain} object to be transformed.>>
#[INPUTS]
#VALUE
# The corresponding \code{nbn} object.
#EXAMPLE
# print8nbn(chain2nbn(rbmn0chain.02),ordering=names(rbmn0nbn.02));
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 12_01_19
#REVISED 13_07_04
#--------------------------------------------
{
  # checking
  # < to be done >
  # constants
  nna <- chain$names;
  nn <- length(nna);
  ooo <- aux2(chain,TRUE); # topological order
  rho <- aux2(chain); # arc orientations
  avan <- c(0,rho)== 1; # parent before
  apre <- c(rho,0)==-1; # parent after
  aaa <- sort(r.numero(chain$roots,chain$names));
  # transforming
  res <- vector("list",nn);
  sigma <- marginal4chain(chain)$sigma;
  mu <- chain$mu;
  for (nnu in ooo) {
    papa <- character(0); # parent specification
    coef <- numeric(0);   # regression coefficients
    # only non roots must be filled
    if (!(nnu %in% aaa)) {
      corri <- 1;
      if (avan[nnu] & apre[nnu]) {
        corri <- aux0(chain$corre[(nnu-1):nnu]);
      } else {
        if (avan[nnu]) { corri <- aux0(chain$corre[nnu-1]);}
        if (apre[nnu]) { corri <- aux0(chain$corre[nnu  ]);}
      }
      if (avan[nnu]) {
        # the node before is a parent
        papa <- nna[nnu-1];
        coef <- chain$sigma[nnu] / sigma[nnu-1] /
                corri * chain$corre[nnu-1];
      }
      if (apre[nnu]) {
        # the node after is a parent
        papa <- c(papa,nna[nnu+1]);
        coef <- c(coef,chain$sigma[nnu] / sigma[nnu+1] /
                       corri * chain$corre[nnu]);
      }
    }
    res[[nnu]]$parents <- papa;
    res[[nnu]]$regcoef <- coef;
    res[[nnu]]$mu <- chain$mu[nnu];             
    res[[nnu]]$sigma <- chain$sigma[nnu];
  }
  names(res) <- nna;
  # reordering the nodes
  ooo <- order4chain(chain);
  res <- res[ooo];
  # testing possible NaNs
  for (ii in r.bf(res)) {
    if (any(is.na(res[[ii]]$regcoef))) {
      stop("Dectected Numerical Inaccuracy in 'chain2nbn'");
    }
  }
  # returning
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
chain2gema <- function(chain)
#TITLE transforms a /chain/ to a /gema/
#DESCRIPTION From a \code{chain} object
# returns the \code{gema} using a direct formulae.\cr
# Much precised than to use the /nbn/ way.
#DETAILS
#KEYWORDS 
#INPUTS
#{chain}<< the \code{chain} object to be transformed.>>
#[INPUTS]
#VALUE
# The corresponding \code{gema} object.
#EXAMPLE
# identical(chain2gema(rbmn0chain.02)$mu,rbmn0gema.02$mu);
# print(chain2gema(rbmn0chain.02)$li-rbmn0gema.02$li);
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 12_02_06
#REVISED 12_02_06
#--------------------------------------------
{
  # checking
  # < to be done >
  # useful constants
  nna <- chain$names;
  nn <- length(nna);
  pare <- aux2(chain,parents=TRUE);
  # initializing
  res <- list(mu=marginal4chain(chain)$mu,
              li=diag(nrow=nn));
  names(res$mu) <- nna;
  dimnames(res$li) <- list(nna,NULL);
  # computing linear combination
  ## for correlations colliding structure
  for (no in r.bc(nn)) {
    if (pare[no] == -1) {
      # a colliding node
      res$li[no,no-1] <- chain$corre[no];
      res$li[no,no] <- aux0(chain$corre[c(no,no+1)]);
      res$li[no,no+1] <- chain$corre[no+1];
    } else {
      if (pare[no] > 0) {
        # a standard node
        if (pare[no] < no) { rho <- chain$corre[no-1];
                    } else { rho <- chain$corre[no];}
        res$li[no,pare[no]] <- rho;
        res$li[no,no] <- aux0(rho);
      }
    }
  }
  ## for correlations sequential segments
  if (length(chain$colliders)>0) {
    ccc <- r.numero(chain$colliders,nna);
  } else {
    ccc <- numeric(0);
  }
  aaa <- r.numero(chain$roots,nna);
  ava <- aaa[-length(aaa)];
  apr <- aaa[-1];
  for (cc in r.bf(ccc)) {
    pivo <- res$li[ccc[cc],ccc[cc]];
    ba <- ava[cc]:ccc[cc];
    bp <- ccc[cc]:apr[cc];
    res$li[ba,ba] <- aux1(chain$corre[ba[-length(ba)]],FALSE);
    res$li[bp,bp] <- aux1(chain$corre[bp[-length(bp)]],TRUE);
    res$li[ccc[cc],ccc[cc]] <- pivo;
  }
  ## without forgetting endind sequential segments
  if (aaa[1] > 1) {
    aa <- 1:aaa[1];
    res$li[aa,aa] <- aux1(chain$corre[1:(aaa[1]-1)],TRUE);
  }
  if (aaa[length(aaa)] < nn) {
    aa <- aaa[length(aaa)]:nn;
    res$li[aa,aa] <- aux1(chain$corre[aaa[length(aaa)]:(nn-1)],FALSE);
  }
  #
  ## adding standard deviations
  masig <- marginal4chain(chain)$sigma;
  res$li <- res$li * masig;
  # reordering the nodes
  ooo <- order4chain(chain);
  res$mu <- res$mu[ooo];
  res$li <- res$li[ooo,ooo];
  # returning
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
chain2correlation <- function(chain)
#TITLE computes the correlation matrix of a chain
#DESCRIPTION returns the correlation
# matrix of a /chain/ object.
#DETAILS
#KEYWORDS 
#INPUTS
#{chain}<< The chain object to consider.>>
#[INPUTS]
#VALUE
# The correlation matrix. It is not
# sorted to respect a topological order
# contrary to \code{chain2mn} function.
#EXAMPLE
# chain2correlation(rbmn0chain.03);
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 11_02_01
#REVISED 11_02_02
#--------------------------------------------
{
  # checking
  # <to be done>
  nn <- length(chain$mu);
  res <- matrix(0,nn,nn);
  dimnames(res) <- list(chain$names,chain$names);
  # finding the milestones
  etat <- state4chain(chain);
  sto <- sort(c(1,nn,which("c"==etat)));
  # introducing the blocks
  for (bb in r.bc(length(sto)-1)) {
    deb <- sto[bb]; fin <- sto[bb+1];
    res[r.bd(deb,fin),r.bd(deb,fin)] <- aux5(chain$corre[r.bd(deb,fin-1)]);
  }
  # returning
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
chain2mn <- function(chain,order=TRUE)
#TITLE computes the distribution of a chain
#DESCRIPTION returns the /mn/ object
# associated to a /chain/ object. Much better
# to use this function that the general function
# \code{nbn2mn} since exact formulae are applied.
#DETAILS
#KEYWORDS 
#INPUTS
#{chain}<< The chain object to consider.>>
#[INPUTS]
#{order} << Must a topological order be imposed?>>
#VALUE
# The resulting /mn/ object. Following the
# convention of \code{mn} objects, a topological
# order is given to it. This is necessary to retrieve
# the associate /nbn/.
#EXAMPLE
# print8mn(chain2mn(rbmn0chain.01));
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 11_12_05
#REVISED 13_07_30
#--------------------------------------------
{
  # checking
  # <to be done>
  # getting the marginal parameters
  margi <- marginal4chain(chain);
  # getting the multivariable parameters
  corre <- chain2correlation(chain);
  # computing the variance matrix
  gamma <- outer(margi$sigma,margi$sigma,"*") * corre;
  # naming the components
  names(margi$mu) <- chain$names;
  dimnames(gamma) <- list(chain$names,chain$names);
  # getting a topological order
  if (order) {
    ooo <- order4chain(chain);
  } else {
    ooo <- r.bf(margi$mu);
  }
  # returning
  list(mu=margi$mu[ooo],
       gamma=gamma[ooo,ooo]);
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
chain2pre <- function(chain,corre=FALSE)
#TITLE computes the precision of a chain
#DESCRIPTION returns the precision matrix
# of a chain, that is the inverse of its
# variance (correlation) matrix. Much better
# to use this function that 
# \code{solve(chain2mn(chain)$gamma)} since
# exact formulae are applied.
#DETAILS
#KEYWORDS 
#INPUTS
#{chain}<< The chain object to consider.>>
#[INPUTS]
#{corre} <<To get the inverse of the correlation
#          matrix instead of.>>
#VALUE
# A dimnamed matrix
#EXAMPLE
# chain2pre(rbmn0chain.02);
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 11_02_02
#REVISED 11_02_06
#--------------------------------------------
{
  # checking
  # <to be done>
  # building the inverse of the correlation matrix
  # constants
  nn <- length(chain$names);
  ## first a long sequential chain
  res <- aux6(chain$corre);
  ## then the 3x3 blocks for colliders
  coli <- sort(which("c"==state4chain(chain)));
  for (cc in coli) {
    ends <- c(cc==2,cc==nn-1);
    ou <- (cc-1):(cc+1);
    res[ou,ou] <- aux7(chain$corre[ou[-3]],ends);
  }
  # adding the variance components
  if (!corre) {
    sisi <- marginal4chain(chain)$sigma;
    res <- res / outer(sisi,sisi,"*");
  }
  # dimnaming
  dimnames(res) <- list(chain$names,chain$names);
  # returning
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
inout4chain <- function(chain)
#TITLE reduces a chain to its inputs and outputs
#DESCRIPTION From a \code{chain} returns the 
# reduced \code{chain} comprising only inputs 
# (that is root nodes) and outputs (that is
# colliders and ends which are not roots)
#DETAILS
#KEYWORDS 
#INPUTS
#{chain}<< The chain object to consider.>>
#[INPUTS]
#VALUE
# The resulting chain
#EXAMPLE
# print8chain(inout4chain(rbmn0chain.02));
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 11_02_07
#REVISED 11_02_07
#--------------------------------------------
{
  # checking
  # <to be done>
  # determining the node to keep
  eta <- state4chain(chain);
  kee <- chain$names[which(eta!="t")];
  nk <- length(kee);
  # building the resulting chain
  res <- vector("list",0);
  res$names <- kee;
  res$roots <- chain$roots;
  res$colliders <- chain$colliders;
  if (nk < 2) {
    res$mu <- chain$mu;
    res$sigma <- chain$sigma;
    res$corre <- chain$corre;
  } else {
    mn <- chain2mn(chain);
    mn$mu <- mn$mu[kee];
    mn$gamma <- mn$gamma[kee,kee,drop=FALSE];
    res$mu <- mn$mu;
    res$sigma <- sqrt(diag(mn$gamma));
    co <- cor4var(mn$gamma);
    res$corre <- diag(co[-1,-nk,drop=FALSE]);
    res <- aux4(res);
  }
  # returning
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
reverse8chain <- function(chain)
#TITLE reverses the nodes of a chain
#DESCRIPTION returns the chain obtained
# after reversing its node order
#DETAILS
#KEYWORDS 
#INPUTS
#{chain}<< The chain object to consider.>>
#[INPUTS]
#VALUE
# The resulting chain
#EXAMPLE
# print8chain(rbmn0chain.02);
# print8chain(reverse8chain(rbmn0chain.02));
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 11_02_07
#REVISED 11_02_07
#--------------------------------------------
{
  # checking
  # <to be done>
  # building it
  res <- vector("list",0);
  res$names <- rev(chain$names);
  res$roots <- chain$roots;
  res$colliders <- chain$colliders;
  res$mu <- rev(chain$mu);
  res$sigma <- rev(chain$sigma);
  res$corre <- rev(chain$corre);
  # returning
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
is8nbn8chain <- function(nbn,order=FALSE)
#TITLE Checks if a given /nbn/ is a /chain/
#DESCRIPTION returns \code{TRUE} [the order] or
# \code{FALSE} [NULL] according that \code{nbn}
# is a chain of not [according to \code{order}].
#DETAILS
#KEYWORDS 
#INPUTS
#{nbn}<< The nbn object to consider.>>
#[INPUTS]
#{order} << When \code{FALSE} the answer to the
# question is returned with \code{TRUE} or \code{FALSE}.\cr
# When \code{TRUE} the chain order
# of the nodes is returned if it is a /chain/
# else \code{NULL}.>>
#VALUE
# A \code{logical(1)} when \code{order} si \code{TRUE} if not
#  the resulting chain order versus NULL.
#EXAMPLE
#is8nbn8chain(rbmn0nbn.01);
#is8nbn8chain(rbmn0nbn.04);
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 11_02_14
#REVISED 11_02_14
#--------------------------------------------
{
  # checking
  # <to be done>
  # constants
  nn <- length(names(nbn));
  # degenerate case
  if (nn == 1) {
    if (order) {
      return(1);
    } else {
      return(TRUE);
    }
  }
  # symmetrical relationship matrix
  rela <- adja4nbn(nbn);
  rela <- rela + t(rela);
  # checking obvious cases
  res <- TRUE;
  nbr <- apply(rela,1,sum);
  if (sum(nbr == 1)!= 2) { res <- FALSE;}
  if (sum(nbr == 2)!= nn-2) { res <- FALSE;}
  if (!res) {
    if (order) {
      return(NULL);
    } else {
      return(FALSE);
    }
  }
  # looking for a chain between nodes
  anc <- which(nbr==1)[1];
  orde <- anc;
  mauvais <- FALSE;
  while ((length(orde)<nn) & !mauvais) {
    nou <- which(rela[anc,]==1);
    if (length(nou) != 1) {
      mauvais <- TRUE;
    } else {
      rela[anc,nou] <- rela[nou,anc] <- 0;
      orde <- c(orde,nou);
      anc <- nou;
    }
  }
  #
  if (mauvais) {
    if (order) {
      res <- NULL;
    } else {
      res <- FALSE;
    }
  } else {
    if (!order) {
      res <- TRUE;
    } else {
      res <- orde;
    }
  }
  # returning
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
nbn2chain <- function(nbn)
#TITLE transforms a /nbn/ into a /chain/
#DESCRIPTION returns the chain obtained
# from \code{nbn} which is supposed to a chain.
# If it is not a chain, an error is issued.
#DETAILS
# It is advised to use \code{is8nbn8chain} before
# calling this function.
#KEYWORDS 
#INPUTS
#{nbn}<< The /nbn/ object to consider.>>
#[INPUTS]
#VALUE
# The resulting chain
#EXAMPLE
# print8chain(nbn2chain(rbmn0nbn.02));
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 11_02_14
#REVISED 11_02_15
#--------------------------------------------
{
  # checking
  # <to be done for /nbn/>
  if (!is8nbn8chain(nbn)) {
    r.erreur(nbn,message="This /nbn/ is not a /chain/");
  }
  # getting one of the two possible chain orders
  ooo <- is8nbn8chain(nbn,TRUE);
  # building the chain except for correlations
  res <- vector("list",0);
  res$names <- names(nbn)[ooo];
  res$sigma <- sapply(nbn,function(a){a$sigma;})[ooo];
  res$mu    <- sapply(nbn,function(a){a$mu;})[ooo];
  nbpa <- sapply(nbn,function(a){length(a$parents);})
  nbpa <- nbpa[ooo];
  res$roots <- res$names[which(nbpa==0)];
  res$colliders <- res$names[which(nbpa==2)];
  res$corre <- rep(0.5,length(res$names)-1);
  # getting a topological order
  too <- order4chain(res);
  res$corre <- rep(0.5,length(res$names)-1);
  sig <- rep(NA,length(res$names));
  # adding the correlations
  if (length(res$names)>1) {
    for (nnu in r.bf(res$names)) {
      nnn <- too[nnu];
      nna <- res$names[nnn];
      if (!(nna %in% res$roots)) {
        papa <- nbn[[nna]]$parents;
        # computing the necessary standard deviations
        sb <- res$sigma[nnn];
        coco <- sasa <- numeric(0);
        for (pn in r.bf(papa)) {
          pp <- papa[pn];
          reco <- nbn[[nna]]$regcoef[pn];
          if (nnn > 1) {if (res$names[nnn-1]==pp) {
            sa <- sig[nnn-1];
            sasa <- c(sasa,sa);
            coco <- c(coco,reco);
          }}
          if (nnn < length(res$names)) {if (res$names[nnn+1]==pp) {
            sa <- sig[nnn+1];
            sasa <- c(sasa,sa);
            coco <- c(coco,reco);
          }}
        }
        # getting the marginal variance
        sb <- sqrt(sum(coco^2*sasa^2)+sb^2);
        sig[nnn] <- sb;
        # computing the correlation(s)
        for (pn in r.bf(papa)) {
          pp <- papa[pn];
          reco <- nbn[[nna]]$regcoef[pn];
          if (nnn > 1) {if (res$names[nnn-1]==pp) {
            sa <- sig[nnn-1];
            res$corre[nnn-1] <- reco * sa / sb;
          }}
          if (nnn < length(res$names)) {if (res$names[nnn+1]==pp) {
            sa <- sig[nnn+1];
            res$corre[nnn] <-  reco * sa / sb;
          }}
        }
      } else {
        # conditional and marginal variances are equal
        sig[nnn] <- res$sigma[nnn];
      }
    }
  }
  # returning
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
chain4chain <- function(chain,nodes,
                        condi=numeric(0),
                        value=rep(0,length(condi)))
#TITLE extracts a chain from a chain
#DESCRIPTION returns the chain obtained
# from \code{chain} retaining only nodes indicated
# by \code{nodes} and conditioned with nodes
# indicated in \code{condi}.
#DETAILS
# Integration is done for nodes not belonging to
# the extracted chain nor being in the conditioning
# subset. Then the distribution of the retained nodes
# is left identical to this in the initial chain.
#KEYWORDS 
#INPUTS
#{chain}<< The chain object to consider.>>
#[INPUTS]
#{nodes} << \code{numeric} (or \code{character}) vector
#         giving the numbers (or names) of the nodes
#         to be retained in the extracted chain.>>
#{condi} << \code{numeric} (or \code{character}) vector
#         giving the numbers (or names) of the
#         conditioning nodes for the extracted chain.>>
#{value} << Numerical values associated to \code{condi}.>>
#VALUE
# The resulting chain
#EXAMPLE
# chain4chain(rbmn0chain.02,c("a","d"),c("b"),12);
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE To be rewritten
#AUTHOR J.-B. Denis
#CREATED 11_02_08
#REVISED 11_02_08
#--------------------------------------------
{
  return("This function is under construction!");
  # checking
  # <to be done>
  # constants
  nn <- length(chain$names);
  # each case in turn
  fait <- FALSE;
  if (is.null(head) & is.null(tail)) {
    # no condition
    res <- chain;
    fait <- TRUE;
  }
  #
  if (!is.null(head) & !is.null(tail)) {
    # both conditioning
    if (nn < 3) {
      res <- NULL;
    } else {
      res <- "'a faire";
    fait <- TRUE;
    }
  }
  #
  if (!fait) {
    if (nn < 2) {
      res <- NULL;
    } else {
      if (is.null(tail)) { chain <- reverse8chain(chain);}
      # condioning has to be done for the last node
      # looking for the last root
      lro <- max(which(state4chain(chain)=="r"));
      ncor <- aux8(chain$corre[r.bd(lro,nn)]);
      ncor <- diag(ncor[-1,-(nn-lro),drop=FALSE]);
      nmu <- chain$mu[-nn];
      nsi <- chain$sigma[-nn];
      res <- list(names=chain$names[-nn],
                  roots=chain$roots,
                  colliders=chain$colliders,
                  mu=nmu,
                  sigma=nsi,
                  corre=ncor
                 );
      fait <- TRUE;
      if (is.null(tail)) { res <- reverse8chain(res);}
    }
  }      
  # returning
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
