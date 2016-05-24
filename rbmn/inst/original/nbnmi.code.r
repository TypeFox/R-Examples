
#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
print8nbn <- function(nbn,what="pr",digits=3,ordering=NULL)
#TITLE print function for a /nbn/ object.
#DESCRIPTION prints a /nbn/ object.
#DETAILS
#KEYWORDS 
#INPUTS
#{nbn} <<\code{nbn} object to be printed.>>
#[INPUTS]
#{what} <<a \code{character(1)}; when comprising "p" the name of each node
#         with its parents are given, when comprising "r the formula
#         regression of each node is given with the node, when comprising
#         "m" the model is given.>>
#{digits} << when not null, the number of digits for rounding.>>
#{ordering} << Nodes are given following the indices of "ordering" if \code{numeric}
#         or the names if it is \code{character}. \code{NULL} means the
#         identity permutation. Repetitions or missing nodes are accepted.>>
#
#VALUE
# Nothing but but \code{nbn} is printed.
#EXAMPLE
# print8nbn(rbmn0nbn.01);
# print8nbn(rbmn0nbn.03,"pm",order=1:2)
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE add the option 'model'
#AUTHOR J.-B. Denis
#CREATED 12_01_14
#REVISED 13_07_04
#--------------------------------------------
{
  # checking
  che <- check8nbn(nbn);
  if (is.character(che)) {
    print(che);
    stop("The provided 'nbn' is not valid!");
  }
  # number of nodes
  nn <- length(nbn);
  # getting the ordering for the nodes
  ordering <- aux3(names(nbn),ordering);
  # printing each asked option
  if (r.expr3present("m",what)) {
    cat("\nModel:",string7dag4nbn(nbn),"\n\n")
  }
  pare <- r.expr3present("p",what);
  regr <- r.expr3present("r",what);
  if (pare | regr) {
    for (i in r.bf(ordering)) {
      ii <- ordering[i];
      if (i==1) {
        cat(r.form3justifie("Nodes",nbc=10,format=3,carac="="));
        cat("===");
        if (pare) {
          cat(r.form3justifie("[parents]",nbc=12,format=1,carac=" "));
        }
        if (regr) {
          cat(r.form3justifie("= Exp. (sd.dev)",nbc=22,format=1,carac=" "));
        }
        cat("\n");
        cat(rep("-",10+3+12+22),collapse="",sep="");
        cat("\n");
      }
      nbpa <- length(nbn[[ii]]$parents);
      cat(r.form3justifie(names(nbn)[ii],nbc=10,format=3,carac="-"));
      cat("---");
      if (pare) {
        if (nbpa == 0) {
          papa = "[-]";
        } else {
          papa <- paste(nbn[[ii]]$parents,collapse=",");
          papa <- paste("[",papa,"]",sep="");
        }
        cat(papa,sep="");
      }
      if (regr) {
        rere <- "  =";
        if ((nbn[[ii]]$mu != 0)| (nbpa==0)) {
          rere <- paste(rere,round(nbn[[ii]]$mu,digits));
          if (nbpa>0) {
            rere <- paste(rere,"+");
          }
        }
        for (pp in r.bc(nbpa)) {
          rere <- paste(rere,
                       paste(round(nbn[[ii]]$regcoef[pp],digits),"*",
                             nbn[[ii]]$parents[pp],sep="")
                      );
          if (pp < nbpa) {
            rere <- paste(rere,"+");
          }
        }
        rere <- paste(rere,
                      paste(" (",round(nbn[[ii]]$sigma,digits),")",sep="")
                     );
        cat(rere);
      }
      cat("\n");
    }
  }
  # returning
  invisible();
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
check8nbn <- function(nbn)
#TITLE checks a /nbn/ object
#DESCRIPTION checks the consistency of \code{nbn} as a /nbn/ object
# issues a fatal error with some clues if inconsistent.
#DETAILS
# Looking a the code of this function provides a way to know which
# are the requirements of a /chain/ object.
#KEYWORDS 
#INPUTS
#{nbn} << The \code{nbn} object to check.>>  
#[INPUTS]
#VALUE
# \code{TRUE} or a \code{character} containing some clue
# about the discovered inconsistency.
#EXAMPLE
# check8nbn(rbmn0nbn.01);
# res <- check8nbn(rbmn0adja.01);
# if (is.na(as.logical(res))) { print(res);}
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 13_07_26
#REVISED 13_07_30
#--------------------------------------------
{
  # initial checking
  if (!is.list(nbn)) {
    return("A /nbn/ must be a list");
  }
  # getting the node numbers
  nbno <- length(nbn);
  # checking the existence of node names
  nano <- names(nbn);
  if (length(unique(nano)) < nbno) {
    return("A /nbn/ must be a named list with different names");
  }
  # checking each node in turn
  res <- character(0);
  for (nn in r.bc(nbno)) {
    na <- names(nbn)[nn];
    if (is.null(nbn[[na]]$mu)) {
      res <- c(res,paste("node",na,"hasn't got a '$mu'"));
    } else {
      if (length(nbn[[nn]]$mu) != 1) {
        res <- c(res,paste("node",na,"has got a '$mu' with a length != to 1"));
      } else {
        if (!is.numeric(nbn[[nn]]$mu)) {
          res <- c(res,paste("node",na,"has got a '$mu' which is not numeric"));
        }
      }
    }
    if (is.null(nbn[[na]]$sigma)) {
      res <- c(res,paste("node",na,"hasn't got a '$sigma'"));
    } else {
      if (length(nbn[[nn]]$sigma) != 1) {
        res <- c(res,paste("node",na,"has got a '$sigma' with a length != to 1"));
      } else {
        if (!is.numeric(nbn[[nn]]$sigma)) {
          res <- c(res,paste("node",na,"has got a '$sigma' which is not numeric"));
        } else {
          if (nbn[[nn]]$sigma < 0) {
            res <- c(res,paste("node",na,"has got a '$sigma' which is negative"));
          }
        }
      }
    }
    if (length(nbn$parents) > 0) {
      papa <- nbn$parents;
      cpa <- TRUE;
      if (length(unique(papa))!=length(papa)) {
        res <- c(res,paste("node",na,"has got duplicated parents"));
        cpa <- FALSE;
      }
      if (length(union(papa,nano))!= length(nano)) {
        res <- c(res,paste("node",na,"has got parents which are not nodes"));
        cpa <- FALSE;
      }
      if (any(na==papa)) {
        res <- c(res,paste("node",na,"is itself its parent"));
        cpa <- FALSE;
      }
      if (cpa) {
        if (length(papa) != length(nbn$regcoef)) {
          res <- c(res,paste("node",na,"has got a '$parents' and '$regcoef' with different lengths"));
        }
        if (!is.numeric(nbn$regcoef)) {
          res <- c(res,paste("node",na,"has got a '$regcoef' which is not numeric"));
        }
      }
    }      
  }
  if (length(res) > 0) { return(res);}
  # returning
  TRUE;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
print8gema <- function(gema,what="ml",ordering=NULL,
                       digits=3,printed=TRUE)
#TITLE standard print function for a /gema/ object.
#DESCRIPTION prints a /gema/ object completely
# or a part of it according to \code{what} specification.
#DETAILS
#KEYWORDS 
#INPUTS
#{gema} <<\code{gema} object to be printed.>>
#[INPUTS]
#{what} <<a \code{character(1)}; when comprising "m" the 
#         expectations are printed, "l" the linear combinations
#         are printed.>>
#{ordering} << Nodes are given following the indices of "ordering" if \code{numeric}
#         or the names if it is \code{character}. \code{NULL} means the
#         identity permutation. Repetitions or missing nodes are accepted.>>
#{digits} << when not null, the number of digits for rounding.>>
#{printed} << \code{TRUE} to issue a printing, if not the prepared matrix
#           is returned.>>
#VALUE
# The \code{gema} is printed or a matrix having \code{nn x ?} is returned
#  binding which elements are precised in the argument \code{what}.
#EXAMPLE
# print8gema(rbmn0gema.01);
# print8gema(rbmn0gema.02,"m");
# print8gema(rbmn0gema.03,"l",digit=1);
# print8gema(rbmn0gema.04,printed=FALSE);
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 12_01_19
#REVISED 12_02_16
#--------------------------------------------
{
  # checking
  che <- check8gema(gema);
  if (is.character(che)) {
    print(che);
    stop("The provided 'gema' is not valid!");
  }
  # number of nodes
  nn <- length(gema$mu);
  nam <- names(gema$mu);
  # getting the ordering for the nodes
  ordering <- aux3(nam,ordering);
  nnr <- length(ordering);
  # initializing
  cnam <- character(0);
  res <- matrix(NA,nn,0);
  # printing each asked option
  if (r.expr3present("m",what)) {
    cnam <- c(cnam,"mu");
    res <- cbind(res,gema$mu);
  }
  if (r.expr3present("l",what)) {
    cnam <- c(cnam,paste("E",r.bc(nnr),sep=""));
    if (length(union(ordering,r.bc(nn))==nn)) {
      uop <- ordering;
    } else {
      uop <- r.bc(nn);
    }
    res <- cbind(res,gema$li[,uop]);
  }
  # permuting
  res <- res[ordering,,drop=FALSE];
  dimnames(res)[[2]] <- cnam;
  # rounding
  if (!is.null(digits)) {
    res <- round(res,digits);
  }
  # returning
  if (printed) {
    print(res);
    invisible();
  } else {
    return(res);
  }
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
check8gema <- function(gema)
#TITLE checks a /gema/ object
#DESCRIPTION checks the consistency of \code{gema} as a /gema/ object
# issues a fatal error with some clues if inconsistent.
#DETAILS
# Looking a the code of this function provides a way to know which
# are the requirements of a /chain/ object.
#KEYWORDS 
#INPUTS
#{gema} << The \code{gema} object to check.>>  
#[INPUTS]
#VALUE
# \code{TRUE} or a \code{character} containing some clue
# about the discovered inconsistency.
#EXAMPLE
# check8gema(rbmn0gema.01);
# res <- check8gema(rbmn0adja.01);
# if (is.na(as.logical(res))) { print(res);}
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 13_07_26
#REVISED 13_07_30
#--------------------------------------------
{
  # initial checking
  if (!is.list(gema)) {
    return("A /gema/ must be a list");
  }
  if (!setequal(names(gema),names(rbmn0gema.01))) {
    return(paste("Names of every /gema/ list must be:",paste(names(rbmn0gema.01),collapse=" ")));
  }
  # getting the node numbers
  nbno <- length(gema$mu);
  # checking the existence of node names
  nano <- names(gema$mu);
  if (length(unique(nano)) < nbno) {
    return("Proposed names for the nodes doesn't fit");
  }
  res <- character(0);
  if (!is.numeric(gema$mu)) {
    res <- c(res,"$mu must be numeric");
  }
  if (!is.numeric(gema$li)) {
    res <- c(res,"$li must be numeric");
  }
  if (!is.matrix(gema$li)) {
    res <- c(res,"$li must be a matrix");
  } else {
    if (!setequal(nano,dimnames(gema$li)[[1]])) {
      res <- c(res,"$li must have node names in rows"); 
    }
    if (nrow(gema$li)!=ncol(gema$li)) {
      res <- c(res,"$li must be a squared matrix"); 
    }
  }
  if (length(res) > 0) { return(res);}
  # returning
  TRUE;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
simulate8gema <- function(gema,nbs) 
#TITLE simulates from a /gema/ object
#DESCRIPTION returns a matrix of simulated
# values with the variable in columns and the 
# simulations in rows.
#DETAILS
# Just the application of the standard formula
# to a white noise. Variables names are taken
# from those of \code{gema$mu}, when these
# does not exist, standard ones are provided.
#KEYWORDS 
#INPUTS
#{gema}<< The \code{gema} object.>>
#{nbs}<< number of simulations to return.>>
#[INPUTS]
#VALUE
# A matrix of size : \code{nbs x length(gema$mu)}
#EXAMPLE
# simulate8gema(rbmn0gema.01,10);
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 12_04_27
#REVISED 12_12_03
#--------------------------------------------
{
  # number of variables and their names
  nbv <- length(gema$mu);
  # 
  if (is.null(names(gema$mu))) {
    va <- paste("V",as.character(r.bc(nbv)),sep="");
  } else {
    va <- names(gema$mu);
  }
  # number of simulations
  nbs <- round(max(0,nbs));
  # simulating
  if (nbv*nbs > 1) {
    res <- matrix(rnorm(nbs*nbv),nbv,nbs);
    res <- gema$mu + gema$li %*% res;
    res <- t(res);
  } else {
    res <- matrix(NA,nbs,nbv);
  }
  # adding the variable names
  dimnames(res) <- list(NULL,va);
  #
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
simulate8nbn <- function(nbn,nbs) 
#TITLE simulates from a /nbn/ object
#DESCRIPTION returns a matrix of simulated
# values with the variable in columns and the 
# simulations in rows.
#DETAILS
# Just the sequential simulations of the nodes
#KEYWORDS 
#INPUTS
#{nbn}<< The \code{nbn} object.>>
#{nbs}<< number of simulations to return.>>
#[INPUTS]
#VALUE
# A matrix of size : \code{nbs x length(nbn)}
#EXAMPLE
# simulate8nbn(rbmn0nbn.01,10);
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 12_04_27
#REVISED 12_12_03
#--------------------------------------------
{
  # number of variables and their names
  nbv <- length(nbn);
  va <- names(nbn);
  # number of simulations
  nbs <- round(max(0,nbs));
  # initializing
  res <- matrix(NA,nbs,nbv);
  dimnames(res) <- list(NULL,va);
  # simulating
  if (nbv*nbs > 1) {
    for (vv in r.bc(nbv)) {
      res[,vv] <- rnorm(nbs,nbn[[vv]]$mu,nbn[[vv]]$sigma);
      for (pp in r.bf(nbn[[vv]]$parents)) {
        pn <- nbn[[vv]]$parents[pp];
        res[,vv] <- res[,vv] + nbn[[vv]]$regcoef[pp] * res[,pn];
      }
    }
  }
  # returning
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
string7dag4nbn <- function(nbn,sep=";") 
#TITLE provides so-called string model of a /nbn/
#DESCRIPTION returns a \code{character(1)}
# describing the dag of the nbn under
# the string form.
#DETAILS
#KEYWORDS 
#INPUTS
#{nbn}<< The nbn.>>
#[INPUTS]
#{sep}<<Separation sign between parents after the
#       conditioning sign (\code{|}).>>
#VALUE
# A \code{character(1)}.
#EXAMPLE
# string7dag4nbn(rbmn0nbn.01);
# string7dag4nbn(rbmn0nbn.04,sep=",");
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE
#AUTHOR J.-B. Denis
#CREATED 12_11_29
#REVISED 12_11_29
#--------------------------------------------
{
  # initializing
  res <- character(0);
  # adding the nodes in turn
  for (nn in r.bf(nbn)) {
    res <- paste(res,"[",names(nbn)[nn],sep="");
    pp <- nbn[[nn]]$parents;
    if (length(pp) >= 1) {
      res <- paste(res,"|",pp[1],sep="");
      for (ip in r.bd(2,length(pp))) {
        res <- paste(res,sep,pp[ip],sep="");
      }
    }
    res <- paste(res,"]",sep="");
  }
  # returning
  res;
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>

#<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<<
diff8nbn <- function(nbn1,nbn2,type=1,scalar=TRUE)
#TITLE returns a score of the difference between two /nbn/s
#DESCRIPTION 
# Returns a positive scalar value measuring, in some way, the difference
# existing within two /nbn/s sharing the same structure.\cr
#DETAILS
# For \code{type==1} it is the canonical euclidian difference between
# all parameters, including the \code{sigma}.
# The score to use to measure the differences between two successive
# estimations is not well established (see the code).
#KEYWORDS 
#INPUTS
#{nbn1} <<First \code{nbn} object.>>
#{nbn2} <<Second \code{nbn} object.>>
#[INPUTS]
#{type}<<When 1, the score includes the difference between the sigmas.
#        When -1, sigmas are not taken into account.>>
#{scalar}<<When \code{TRUE} the squared norm is returned, if not the vector of difference.>>
#VALUE
# Either a scalar or a named vector (according to \code{scalar}).
#EXAMPLE
# diff8nbn(rbmn0nbn.01,rbmn0nbn.01);
# diff8nbn(rbmn0nbn.01,rbmn0nbn.01,scalar=FALSE);
#REFERENCE
#SEE ALSO
#CALLING
#COMMENT
#FUTURE 
#AUTHOR J.-B. Denis
#CREATED 13_04_25
#REVISED 13_07_30
#--------------------------------------------
{
  # checking
  # < to be done >
  nona <- names(nbn1);
  # constituting the vector of differences.
  didi <- numeric(0);
  for (nn in r.bf(nona)) {
    suff <- nona[nn];
    n0 <- length(didi);
    didi <- c(didi,nbn2[[nn]]$mu-nbn1[[nn]]$mu);
    nana <- paste(suff,"mu",sep=".");
    if (type > 0) {
      didi <- c(didi,nbn2[[nn]]$sigma-nbn1[[nn]]$sigma);
      nana <- c(nana,paste(suff,"sigma",sep="."));
    }
    if (length(nbn2[[nn]]$parents)>0) {
      didi <- c(didi,nbn2[[nn]]$regcoef-nbn1[[nn]]$regcoef);
      nana <- c(nana,paste(suff,nbn2[[nn]]$parents,sep="."));
    }
    names(didi)[r.bd(n0+1,length(didi))] <- nana;
  }
  # returning
  if (scalar) {
    return(sum(didi^2));
  } else {
    return(didi);
  }
  invisible();
}
#>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>>
