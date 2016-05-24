#******************************************************************************* 
#
# Bayesian Regression and Adaptive Sampling with Gaussian Process Trees
# Copyright (C) 2005, University of California
# 
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
# 
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# Lesser General Public License for more details.
# 
# You should have received a copy of the GNU Lesser General Public
# License along with this library; if not, write to the Free Software
# Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301  USA
#
# Questions? Contact Robert B. Gramacy (rbgramacy@ams.ucsc.edu)
#
#*******************************************************************************


## tgp:
##
## the master tgp R function which checks for valid inputs and
## calls the C-side via .C on those inputs -- and then calls the
## post-processing code accordingly

"tgp" <-
function(X, Z, XX=NULL, BTE=c(2000,7000,2), R=1, m0r1=FALSE, linburn=FALSE,
         params=NULL, itemps=NULL, pred.n=TRUE, krige=TRUE, zcov=FALSE,
         Ds2x=FALSE, improv=TRUE, sens.p=NULL, trace=FALSE, verb=1, rmfiles=TRUE)
{
  ## (quitely) double-check that tgp is clean before-hand
  tgp.cleanup(message="NOTICE", verb=verb, rmfiles=TRUE);

  ## what to do if fatally interrupted?
  on.exit(tgp.cleanup(verb=verb, rmfiles=rmfiles))

  ## check for two unsupported combinations of modeling options
  if(params$corr == "mrexpsep" && linburn)
    stop("Sorry, the linear burn-in is not available for corr=\"mrexpsep\"")
  if(params$corr == "mrexpsep" && !is.null(sens.p))
    stop("Sorry, sensitivity analysis is not available for corr=\"mrexpsep\"")

  ## get names
  Xnames <- names(X)
  response <- names(Z)
  
  ## check X and Z
  XZ <- check.matrix(X, Z) 
  X <- XZ$X; Z <- XZ$Z
  n <- nrow(X); d <- ncol(X)
  if(is.null(n)) stop("nrow(X) is NULL")

  ## check XX
  XX <- check.matrix(XX)$X
  if(is.null(XX)) { nn <- 0; XX <- matrix(0); nnprime <- 0 }
  else { 
    nn <- nrow(XX); nnprime <- nn 
    if(ncol(XX) != d) stop("mismatched column dimension of X and XX");
  }

  ## check that trace is true or false)
  if(length(trace) != 1 || !is.logical(trace))
    stop("trace argument should be TRUE or FALSE")
  else if(trace) {
    if(3*(10+d)*(BTE[2]-BTE[1])*R*(nn+1)/BTE[3] > 1e+7)
      warning(paste("for memory/storage reasons, ",
                    "trace not recommended when\n",
                    "\t 3*(10+d)*(BTE[2]-BTE[1])*R*(nn+1)/BTE[3]=",
                    3*(10+d)*(BTE[2]-BTE[1])*R*(nn+1)/BTE[3], " > 1e+7.\n",
                    "\t Try reducing nrow(XX)", sep=""), immediate.=TRUE)
  }

  ## check that pred.n, krige, and Ds2x is true or false
  if(length(pred.n) != 1 || !is.logical(pred.n))
    stop("pred.n should be TRUE or FALSE")
  if(length(krige) != 1 || !is.logical(krige))
    stop("krige should be TRUE or FALSE")
  if(length(zcov) != 1 || !is.logical(zcov))
    stop("zcov should be TRUE or FALSE")
  if(length(Ds2x) != 1 || !is.logical(Ds2x))
    stop("Ds2x should be TRUE or FALSE")

  ## check the form of the improv-power argument
  if(length(improv) == 2) { numirank <- improv[2]; improv <- improv[1] }
  else { numirank <- NULL }
  if(length(improv) != 1 || !(is.logical(improv) ||
             is.numeric(improv)) || (is.numeric(improv) && improv <= 0))
    stop(paste("improv [", improv,
               "] should be TRUE, FALSE, or a positive integer (power)", sep=""))
  g <- as.numeric(improv)

  ## check numirank, which is improv[2] in input
  if(is.null(numirank) && improv) numirank <- nn ## max(min(10, nn), 0.1*nn)
  else if(!is.null(numirank) && numirank > nn) stop("improv[2] must be <= nrow(XX)")
  else if(is.null(numirank)) numirank <- 0
  
  ## check for inconsistent XX and Ds2x/improv
  if(nn == 0 && (Ds2x || improv))
    warning("need to specify XX locations for Ds2x and improv")

  ## check the sanity of input arguments
  if(nn > 0 && sum(dim(XX)) > 0 && ncol(XX) != d) stop("XX has bad dimensions")
  if(length(Z) != n) stop("Z does not have length == nrow(Z)")
  if(BTE[1] < 0 || BTE[2] <= 0 || BTE[1] > BTE[2]) stop("bad B and T: must have 0<=B<=T")
  if(BTE[3] <= 0 || ((BTE[2]-BTE[1] != 0) && (BTE[2]-BTE[1] < BTE[3])))
    stop("bad E arg: if T-B>0, then must have T-B>=E")
  if((BTE[2] - BTE[1]) %% BTE[3] != 0) stop("E must divide T-B")
  if(R < 0) stop("R must be positive")
  
  ## deal with params
  if(is.null(params)) params <- tgp.default.params(d)

  ## check if X is of full rank
  if(params$meanfn == "linear" &&
     class(try(solve(t(X[,1:params$tree[5]]) %*% X[,1:params$tree[5]]),
               silent=TRUE)) == "try-error") {
    stop("X[,1:", params$tree[5], "]-matrix is not of full rank", sep="")
  }
  
  ## convert params into a double-vector for passing to C
  dparams <- tgp.check.params(params, d);
  if(is.null(dparams)) stop("Bad Parameter List")

  ## check starting importance-tempering inv-temp
  itemps <- check.itemps(itemps, params)
  
  ## might scale Z to mean of 0 range of 1
  if(m0r1) { Zm0r1 <- mean0.range1(Z); Z <- Zm0r1$X }
  else Zm0r1 <- NULL

  ## if performining a sensitivity analysis, set up XX ##
  if(!is.null(sens.p)) {
    if(nn > 0) warning("XX generated online in sensitivity analyses")
    nnprime <- 0
    sens.par <- check.sens(sens.p, d)
    nn <- sens.par$nn; nn.lhs <- sens.par$nn.lhs; XX <- sens.par$XX
    ngrid <- sens.par$ngrid; span <- sens.par$span
    MEgrid <- as.double(sens.par$MEgrid)
    if(verb >= 2) 
      cat(paste("Predict at", nn, "LHS XX locs for sensitivity analysis\n"))
  } else{ nn.lhs <- ngrid <- 0; MEgrid <- span <- double(0) }

  ## construct the set of candidate split locations
  Xsplit <- X
  if(is.null(sens.p) && nn > 0) Xsplit <- rbind(Xsplit, XX)

  ## for sens
  S = R*(BTE[2]-BTE[1])/BTE[3]
  
  # RNG seed
  state <- sample(seq(0,999), 3)
  
  ## run the C code
  ll <- .C("tgp",
           
           ## begin inputs
           state = as.integer(state),
           X = as.double(t(X)),
           n = as.integer(n),
           d = as.integer(d),
           Z = as.double(Z),
           XX = as.double(t(XX)),
           nn = as.integer(nn),
           Xsplit = as.double(t(Xsplit)),
           nsplit = as.integer(nrow(Xsplit)),
           trace = as.integer(trace),
           BTE = as.integer(BTE),
           R = as.integer(R),
           linburn = as.integer(linburn),
           zcov = as.integer(zcov),
           g = as.integer(c(g, numirank)),
           dparams = as.double(dparams),
           itemps = as.double(itemps),
           verb = as.integer(verb),
           tree = as.double(-1),
           hier = as.double(-1),
           MAP = as.integer(0),
           sens.ngrid = as.integer(ngrid),
           sens.span = as.double(span),
           sens.Xgrid = as.double(MEgrid),

           ## output dimensions for checking NULL
           pred.n = as.integer(pred.n),
           nnprime = as.integer(nnprime),
           krige = as.integer(krige),
           bDs2x = as.integer(Ds2x),
           bimprov = as.integer(as.logical(improv) * nnprime),
           
           ## begin outputs
           Zp.mean = double(pred.n * n),
           ZZ.mean = double(nnprime),
           Zp.km = double(krige * pred.n * n),
           ZZ.km = double(krige * nnprime),
           Zp.vark = double(krige * pred.n * n),
           ZZ.vark = double(krige * nnprime),
           Zp.q = double(pred.n * n),
           ZZ.q = double(nnprime),
           Zp.s2 = double(pred.n * (zcov*n^2 + (!zcov)*n)),
           ZZ.s2 = double(zcov*nnprime^2 + (!zcov)*nnprime),
           ZpZZ.s2 = double(pred.n * n * nnprime * zcov),
           Zp.ks2 = double(krige * pred.n * n),
           ZZ.ks2 = double(krige * nnprime),
           Zp.q1 = double(pred.n * n),
           Zp.med = double(pred.n * n),
           Zp.q2 = double(pred.n * n),
           ZZ.q1 = double(nnprime),
           ZZ.med = double(nnprime),
           ZZ.q2 = double(nnprime),
           Ds2x = double(Ds2x * nnprime),
           improv = double(as.logical(improv) * nnprime),
           irank = integer(as.logical(improv) * nnprime),
           ess = double(1 + itemps[1]*2),
           gpcs = double(4),
           sens.ZZ.mean = double(ngrid*d),
           sens.ZZ.q1 = double(ngrid*d),
           sens.ZZ.q2 = double(ngrid*d),
           sens.S = double(d*S*!is.null(sens.p)),
           sens.T = double(d*S*!is.null(sens.p)),
           
           ## end outputs
           PACKAGE = "tgp")

  ## all post-processing is moved into a new function so it
  ## can be shared by predict.tgp()
  ll <- tgp.postprocess(ll, Xnames, response, pred.n, zcov, Ds2x, improv,
                        sens.p, Zm0r1, params, rmfiles)
  
  return(ll)
}
