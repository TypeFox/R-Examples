## Library for evaluating various properties of liability networks

###########################################################################
#The following function computes the relative liabilities matrix
#from the total liabilities matrix and the external liabilities myexliab.
###########################################################################
computerelativeliab <-function(liabmat, myexliab){
    pimatrix<-liabmat*0;
    n<-dim(liabmat)[1];
    totalliab <-rowSums(liabmat) + myexliab;
    for (i in 1:n){
        if(totalliab[i]>0){pimatrix[i,] <- liabmat[i, ]/totalliab[i]; }
    }
    return(pimatrix);
}

###################################################################################################################
#The following function computes the clearing vector (using linear programming)
#for a network with n banks, relative liabilities matrix Mypi
#external assets mye and total liabilities mypbar; it returns the clearing vector.
#####################################################################################################################
clearingvector <-function(n, Mypi, mye, mypbar){
    pstar<-mat.or.vec(n, 1);
    default<-mat.or.vec(n, 1);
    MyI <-diag(n);
    MyA <-MyI - t(Mypi);
    bl<-mat.or.vec(n, 1);
    bu<-mypbar;
    f.obj<-mat.or.vec(n, 1)+1;
    f.con<-rbind(MyA, MyI, MyI);
    f.dir<-c(rep ("<=", n), rep (">=", n), rep ("<=", n));
    f.rhs<-c(mye, mat.or.vec(n, 1), mypbar);
    lp(direction="max", f.obj, f.con, f.dir, f.rhs)
    return(lp(direction="max", f.obj, f.con, f.dir, f.rhs)$solution);
}

######
## Computes Defaults and Clearing Vector with external liabilities.
#####
defaults_withexliab <-function(L, ea, el){
    Mypi <- computerelativeliab(L, el);
    Lbars <- rowSums(L)
    totalliab <-Lbars+el;
    n <- length(ea)
    clearingvec<-clearingvector(n, Mypi, ea,totalliab);
    list(defaultind=as.numeric((clearingvec+1e-7)<totalliab),clearingvec=clearingvec)
}

#' Default Cascade
#'
#' Computes bank defaults via the default cascade algorithm.
#'
#' @param L liability matrix
#' @param ea vector of external assets
#' @param el vector of external liabilites (default 0)
#' @param recoveryrate recovery rate in [0,1] (defaults to 0)
#'
#' @return vector indicating which banks default (1=default, 0= no default)
#'
#' @examples
#' ea <- c(1/2,5/8,3/4)
#' el <- c(3/2,1/2,1/2)
#' x <- 0.5
#' L <- matrix(c(0,x,1-x,1-x,0,x,x,1-x,0),nrow=3)
#' default_cascade(L,ea,el)
#' @export
default_cascade<-function(L, ea, el=0,recoveryrate=0){
    n <- length(ea)
    Lbars <- rowSums(L);
    Abars <- colSums(L);
    networth<-Abars+ea-Lbars-el;
    dind <-mat.or.vec(n, 1);
    if(any(networth<0)){
        dind[networth<0]<-1;
        vec1 <-rep(1 - recoveryrate,n);
        ndef <- sum(dind)
        for(nu in 1:(n-1)){
            loss <- (vec1*dind)%*%L
            dind <- as.numeric(networth<loss)
            ndefnew <- sum(dind)
            if (ndefnew==ndef)
                return(list(defaultind=dind));
            ndef <- ndefnew
        }
    }
    list(defaultind=dind);
}


#' @title Clearing Vector with Bankruptcy Costs
#'
#' @description
#' Computes bank defaults for the clearing vector approach without and
#' with bankruptcy costs (Eisenberg and Noe, 2001), (Rogers and Veraart, 2013).
#'
#' @details
#' Without bankruptcy costs the approach of Eisenberg and Noe (2001)
#' is used using a linear programme.  With bankruptcy costs, the
#' implementation is based on the Greatest Clearing Vector Algorithm (GA),
#' see Definition 3.6, Rogers & Veraart (2013).
#'
#' @param L Liabilities matrix
#' @param ea Vector of external assets
#' @param el Vector of external liabilites  (default 0)
#' @param alpha 1-proportional default costs on external assets in [0,
#' 1] (default to 1).
#' @param beta 1-proportional default costs on interbank assets in [0,
#' 1] (defaults to 1).
#' @return A list consisting of a vector indicating which banks
#' default (1=default, 0= no default) and the greatest clearing
#' vector.
#'
#' @references
#' Eisenberg, L. and Noe, T.H. (2001). Systemic risk in financial
#' systems. Management Science 47, 236--249.
#'
#' Rogers, L. C. G. and Veraart, L. A. M. (2013) Failure and Rescue in
#' an Interbank Network, Management Science 59 (4), 882--898.
#'
#' @examples
#' ea <- c(1/2,5/8,3/4)
#' el <- c(3/2,1/2,1/2)
#' x <- 0.5
#' L <- matrix(c(0,x,1-x,1-x,0,x,x,1-x,0),nrow=3)
#' default_clearing(L,ea,el)
#' default_clearing(L,ea,el, alpha=0.5, beta=0.7)
#' @export
default_clearing<-function(L, ea, el=0, alpha=1, beta=1){
    if (alpha==1&&beta==1)
        return(defaults_withexliab(L,ea,el))
    n <- length(ea)
    Lbars <- rowSums(L);
    Abars <- colSums(L);
    Mypi<-computerelativeliab(L, el);
    totalliab<-Lbars+el;
    networth<-Abars+ea-Lbars-el;
    dind <-mat.or.vec(n, 1);
    Defaults<-which(networth<0);
    Nondefaults<-which(networth>=0);
    dind[networth<0]<-1;

    defaultcount1<-sum(dind);
    defaultcount2<-0;
    Lambda<-totalliab;
    iterationcount<-0;
    while(defaultcount1>defaultcount2){
        Lambda[Nondefaults]<-totalliab[Nondefaults];
	tmpm<-diag(defaultcount1)-beta*t(Mypi[Defaults, Defaults]);
	if(sum(dind)==n){
            tmpv<-alpha*ea[Defaults];
	} else {
            tmpv<-alpha*ea[Defaults]+beta*t(t(totalliab[Nondefaults])%*%Mypi[Nondefaults, Defaults]);
	}
	Lambda[Defaults]<-solve(tmpm, tmpv);
        iterationcount<-iterationcount+1;
	v<-t(Mypi)%*%Lambda + ea - totalliab;
        dind[v<0]<-1;
	Defaults<-which(v<0);
        Nondefaults<-which(v>=0);
	defaultcount2<-defaultcount1;
        defaultcount1<-sum(dind);
    }
    list(defaultind=as.numeric((Lambda+1e-7)<totalliab),clearingvec=Lambda)
}


#' @title Default of Banks
#'
#' @description
#' Computes bank defaults based on a liabilities matrix and external assets and liabilities.
#'
#'
#' @param L liability matrix
#' @param ea vector of external assets
#' @param el vector of external liabilites.
#' @param method the method to be used. See Details.
#' @param ... Additional information for the various methods. See Details.
#'
#' @return A list with at least one element "defaultind", which is a
#' vector indicating which banks default (1=default, 0= no
#' default). Depending on the method, other results such as the
#' clearing vector may also be reported.
#' @seealso \code{\link{default_cascade}},  \code{\link{default_clearing}},
#'
#' @examples
#' ea <- c(1/2,5/8,3/4)
#' el <- c(3/2,1/2,1/2)
#' x <- 0.5
#' L <- matrix(c(0,x,1-x,1-x,0,x,x,1-x,0),nrow=3)
#' default(L,ea,el)
#' default(L,ea,el,"cascade")
#'
#' @export
default <- function(L, ea, el=0,method=c("clearing","cascade"),...){
    method <- match.arg(method)
    if (method=="clearing"){
        return(default_clearing(L,ea,el,...))
    }
    if (method=="cascade"){
        return(default_cascade(L,ea,el,...));
    }
}


