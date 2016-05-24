
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
estimate8nbn <- function(nbn,data) 
#TITLE estimating the /nbn/ parameters
#DESCRIPTION 
# From a /nbn/ to describe the DAG, and a data.frame
# containing the necessary observations, returns the /nbn/ with
# all its parameters newly estimated.
#DETAILS
# No constraints are put on the parameters.
#KEYWORDS 
#INPUTS
#{nbn}<< The initial /nbn/.>>
#{data}<<The data frame comprising all /nbn/ nodes.>>
#[INPUTS]
#VALUE
# The resulting /nbn/ with the estimated parameters.
#EXAMPLE
# data(boco);
# print8nbn(rbmn0nbn.05);
# print8nbn(estimate8nbn(rbmn0nbn.05,boco));
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
  # To be done
  # constant
  nona <- names(nbn);
  # looping onto the node set
  for (nn in r.bf(nbn)) {
    YY <- nona[nn];
    if (length(nbn[[nn]]$parents)>0) {
      XX <- paste(nbn[[nn]]$parents,collapse="+");
    } else {
      XX <- "1";
    }
    FF <- paste(YY,"~",XX);
    flm <- lm(as.formula(FF),data=data);
    ano <- anova(flm);
    # 13_05_02 modif performed to obtain /bnlearn/ results
    #nbn[[nn]]$sigma <- sqrt(ano["Residuals",3]);
    nbn[[nn]]$sigma <- sqrt(ano["Residuals",2]/sum(ano[,1]));
    # end of 13_05_02 modification
    nbn[[nn]]$mu <- coefficients(flm)[1];
    nbn[[nn]]$parents <- names(coefficients(flm)[-1]);
    nbn[[nn]]$regcoef <- coefficients(flm)[-1];
    names(nbn[[nn]]$regcoef) <- NULL;
  }
  # returning
  nbn;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
estimate8arcs <- function(nbn,sarc,data)
#TITLE estimates a common regression coefficient to a set of arcs
#DESCRIPTION 
# Supposing known all parameters but one in a /nbn/ object, the function
# provides the esmate for the missing one which is supposed common to
# arcs described in the argument \code{sarc}. Constant term is included
# in the optimization.
#DETAILS
# This is mainly an ancillary function of \code{estimate8constrainednbn}\cr
# Two coefficient are updated: the regression coefficient and
# the estimate of the standard deviation.
#KEYWORDS 
#INPUTS
#{nbn} <<\code{nbn} object.>>
#{sarc} <<Matrix with two columns indicating the tails (1rst column) and the
#         heads (2d column) of the arcs having a common parameter. It is
#         checked that these arcs are indeed included in \code{nbn}.
#         Nodes must be indicated by their names (not their number).>>
#{data} <<Data frame to be used for the estimation. It must
#        comprise all necessary nodes (not only those involved
#        in \code{sarc} but also the remaining parents of \code{sarc[,2]}.
#        Usually, all used variables are centred but this is not
#        required.>>
#[INPUTS]
#VALUE
# the resulting /nbn/ object with the common parameter.
#EXAMPLE
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE 
#AUTHOR J.-B. Denis
#CREATED 13_04_24
#REVISED 13_04_25
#--------------------------------------------
{
  # checking
  # < to be done >
  # initializing the data frame
  don <- as.data.frame(matrix(NA,0,2));
  www <- numeric(0);
  nda <- nrow(data);
  # looping onto the set of arcs to build the estimating data frame
  # of the targetted parameter
  for (aa in r.bc(nrow(sarc))) {
    # getting the part to add
    tail <- sarc[aa,1]; head <- sarc[aa,2];
    noeud <- nbn[[head]];
    if (!(tail %in% noeud$parents)) {
      stop(paste(tail,"is not a parent of",head));
    }
    reco <- noeud$regcoef; names(reco) <- noeud$parents;
    reco <- reco[setdiff(noeud$parents,tail)];
    XX <- data[[tail]];
    YY <- data[[head]];
    for (pp in names(reco)) {
      YY <- YY - reco[pp]*data[[pp]];
    }
    newi <- paste(head,r.bc(nda),sep=".");
    dos <- matrix(c(XX,YY),ncol=2,
                  dimnames=list(newi,NULL));
    # adding it to the data.frame
    don <- rbind(don,dos);
    www <- c(www,rep(noeud$sigma^-2,nda));
  }
  dimnames(don)[[2]] <- c("X","Y");
  # estimating the parameter
  elm <- lm(Y ~ -1+X,data=don,weights=www);
  # reporting the estimation in the /nbn/
  for (aa in r.bc(nrow(sarc))) {
    # incorporating the new regression coefficient
    tail <- sarc[aa,1]; head <- sarc[aa,2];
    ou <- which(nbn[[head]]$parents==tail);
    nbn[[head]]$regcoef[ou] <- coefficients(elm);
    # computing the new residual sum of squares
    #cat("<",tail,"->",head,">\n"); browser();
    rss <- sum((
           data[,head] -
           as.matrix(data[,nbn[[head]]$parents,drop=FALSE]) %*%
           matrix(nbn[[head]]$regcoef,ncol=1))^2);
    # incorporating the estimation of the new sigma
    # 13_05_02 modif performed to obtain /bnlearn/ results
    #nbn[[head]]$sigma <- sqrt(rss/(nda-length(nbn[[head]]$parents)));
    nbn[[head]]$sigma <- sqrt(rss/(nda-1));
    # end of 13_05_02 modification
  }
  # returning
  nbn;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
estimate8constrainednbn <- function(nbn,sarc,data,
                                    imp=0,nite=10,eps=10^-5)
#TITLE estimates the parameters of a nbn with equality constraints
#DESCRIPTION 
# Estimations of the parameters of a /nbn/ is done when there
# are some equality constraints onto the regression coefficients.\cr
# Constant terms (\code{mu}) and conditional standard deviations (\code{sigma})
# are supposed independent (that is not constrained with equalities).\cr
# Equality constraints are given by \code{sarc}, a list of matrices
# with two columns, indicating each the series of arcs having the same
# regression coefficient.
#DETAILS
# Not linked regression coefficients doesn't require to be included in
# \code{sarc}, the function do it by itself.\cr
# The score to use to measure the differences between two successive
# estimations is not well established (see the code).
#KEYWORDS 
#INPUTS
#{nbn} <<\code{nbn} object.>>
#{sarc} <<List of Matrices with two columns indicating the tails (1rst column) and the
#         heads (2d column) of the arcs having a common parameter. It is
#         checked that these arcs are indeed included in \code{nbn}.
#         Nodes must be indicated by their names (not their number).>>
#{data} <<Data frame to be used for the estimation. It must
#        comprise all necessary nodes (not only those involved
#        in \code{sarc} but also the remaining parents of \code{sarc[,2]}.
#        Usually, all used variables are centred but this is not
#        required.>>
#[INPUTS]
#{imp}<<When \code{0} nothing displayed. When \code{1} the number of iterations is displayed.
#       When \code{2} the successive values of the criterion are also displayed. >>
#{nite}<<Maximum number of iterations.>>
#{eps}<<relative difference in successive scores needed to stop the iterations.>>
#VALUE
# The resulting /nbn/ object with the estimated parameters.
#EXAMPLE
# data(boco);
# print8nbn(rbmn0nbn.05);
# print8nbn(estimate8nbn(rbmn0nbn.05,boco));
# print8nbn(estimate8constrainednbn(rbmn0nbn.05,rbmn0crarc.05,boco));
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE 
#AUTHOR J.-B. Denis
#CREATED 13_04_25
#REVISED 13_05_21
#--------------------------------------------
{
  # checking
  # < to be done >
  sep <- " / ";
  # detecting and checking the proposed arcs for equality
  nona <- names(nbn);
  aaa <- adja4nbn(nbn);
  for (a1 in r.bf(sarc)) {
    ma <- sarc[[a1]];
    for (a2 in r.bc(nrow(ma))) {
      tail <- ma[a2,1]; head <- ma[a2,2];
      if (aaa[tail,head] == 0) {
        r.erreur(paste(tail,"->",head),message="This arc is not included in the /nbn/");
      } else {
        aaa[tail,head] <- aaa[tail,head] + 1;
      }
    }
  }
  # adding the remaining arcs alone in the list
  for (ii in r.bf(nona)) { for (jj in r.bf(nona)) {
    if (aaa[ii,jj] == 1) {
      sarc[[length(sarc)+1]] <- matrix(c(nona[ii],nona[jj]),1);
    }
  }}
  # centring every node
  mimi <- sapply(data,mean);
  for (nn in nona) {
    data[[nn]] <- data[[nn]] - mimi[nn];
  }
  # starting the estimation without constraint
  nbna <- estimate8nbn(nbn,data);
  nite <- round(max(nite,1));
  if (imp>1) { cat("estimate8constrainednbn:",sep,sep="");}
  # iterating for the estimation
  nbite <- 0; dsco <- 2*eps;
  while ((nbite < nite) & (dsco > eps)) {
    nbite <- nbite+1;
    # looping onto the parameters
    for (pp in r.bf(sarc)) {
      nbna <- estimate8arcs(nbna,sarc[[pp]],data);
    }
    # calculating the score
    if (nbite > 1) {
      dsco <- diff8nbn(nbna,nbnn);
      if (imp>1) { cat(dsco,sep,sep="");}
    }
    # updating the /nbn/
    nbnn <- nbna;
  }
  if (imp>1) { cat("\n");}
  if (imp>0) { cat("nb.ite =",nbite,"(",dsco,")\n");}
  # reintroducing the mean in every node
  for (nn in nona) {
    mu <- mimi[nn];
    mu <- mu - sum(mimi[nbnn[[nn]]$parents]*nbnn[[nn]]$regcoef);
    nbnn[[nn]]$mu <- mu;
  }
  # returning
  nbnn;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
