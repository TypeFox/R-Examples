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


## predict.tgp:
##
## generic the master tgp R function which takes most of its inputs
## from a tgp-object.  Most of the changeable outputs have to do with
## sampling from the posterior predictive distribution (hence a predict
## method).  It checks for valid inputs and then calls the C-side via .C
## on those inputs -- and then calls the post-processing code accordingly

"predict.tgp" <-
function(object, XX=NULL, BTE=c(0,1,1), R=1, MAP=TRUE, pred.n=TRUE, krige=TRUE,
         zcov=FALSE, Ds2x=FALSE, improv=FALSE, sens.p=NULL, trace=FALSE, verb=0, ...)
{
  ## (quitely) double-check that tgp is clean before-hand
  tgp.cleanup(message="NOTICE", verb=verb, rmfiles=TRUE);
  
  ## what to do if fatally interrupted?
  on.exit(tgp.cleanup(verb=verb, rmfiles=TRUE))

  if(object$params$corr == "mrexpsep" && !is.null(sens.p))
    stop("Sorry, sensitivity analysis is not yet available for corr=\"mrexpsep\"")

  ## get names
  Xnames <- names(object$X)
  response <- names(object$Z)

  ## check XX
  XX <- check.matrix(XX)$X
  if(is.null(XX)) { nn <- 0; XX<- matrix(0);  nnprime <- 0 }
  else { 
    nn <- nrow(XX); nnprime <- nn 
    if(ncol(XX) != object$d) stop("mismatched column dimension of object$X and XX");
  }

  ## check that pred.n, krige, MAP, and Ds2x is true or false
  if(length(pred.n) != 1 || !is.logical(pred.n))
    stop("pred.n should be TRUE or FALSE")
  if(length(krige) != 1 || !is.logical(krige))
    stop("krige should be TRUE or FALSE")
  if(length(zcov) != 1 || !is.logical(zcov))
    stop("zcov should be TRUE or FALSE")
  if(length(MAP) != 1 || !is.logical(MAP))
    stop("MAP should be TRUE or FALSE")
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
  if(is.null(numirank) && improv) numirank <- max(min(10, nn), 0.1*nn)
  else if(!is.null(numirank) && numirank > nn) stop("improv[2] must be <= nrow(XX)")
  else if(is.null(numirank)) numirank <- 0
  
  ## check for inconsistent XX and Ds2x/improv
  if(nn == 0 && (Ds2x || improv))
    warning("need to specify XX locations for Ds2x and improv")

  ## check the sanity of input arguments
  if(nn > 0 && sum(dim(XX)) > 0 && ncol(XX) != object$d) stop("XX has bad dimensions")
  if(BTE[1] < 0 || BTE[2] <= 0 || BTE[1] >= BTE[2]) stop("bad B and T: must have 0<=B<T")
  if(BTE[3] <= 0 || BTE[2]-BTE[1] < BTE[3]) stop("bad E arg: must have T-B>=E")
  
  ## might scale Z to mean of 0 range of 1
  if(object$m0r1) { Zm0r1 <- mean0.range1(object$Z); Z <- Zm0r1$X }
  else { Z <- object$Z; Zm0r1 <- NULL }

  ## get infor about the tree
  m <- which.max(object$posts$lpost)
  t2c <- tree2c(object$trees[[object$posts$height[m]]])
  
  # RNG seed
  state <- sample(seq(0,999), 3)

  ## get itemps from object, but set c0n0 <- c(0,0) 
  ## so no stochastic approx happens
  object$itemps$c0n0 <- c(0,0)
  itemps <- check.itemps(object$itemps, object$params)

  ## if performing a sensitivity analysis, set up XX 
  if(!is.null(sens.p)){
    nnprime <- 0
    if(nn > 0) warning("XX generated online in sensitivity analyses")
    sens.par <- check.sens(sens.p, object$d)
    nn <- sens.par$nn; nn.lhs <- sens.par$nn.lhs; XX <- sens.par$XX
    ngrid <- sens.par$ngrid; span <- sens.par$span
    MEgrid <- as.double(sens.par$MEgrid)
    if(verb >= 1) 
      cat(paste("Predict at", nn, "LHS XX locs for sensitivity analysis\n"))
  }
  else{ nn.lhs <- ngrid <- 0; MEgrid <- span <- double(0) }

  ## calculate the number of sampling rounds
  S = R*(BTE[2]-BTE[1])/BTE[3]

  ## run the C code
  ll <- .C("tgp",

           ## begin inputs
           state = as.integer(state),
           X = as.double(t(object$X)),
           n = as.integer(object$n),
           d = as.integer(object$d),
           Z = as.double(Z),
           XX = as.double(t(XX)),
           nn = as.integer(nn),
           Xsplit = as.double(t(object$Xsplit)),
           nsplit = as.integer(nrow(object$Xsplit)),
           trace = as.integer(trace),
           BTE = as.integer(BTE),
           R = as.integer(R),
           linburn = as.integer(FALSE),
           zcov = as.integer(zcov),
           g = as.integer(c(g, numirank)),
           dparams = as.double(object$dparams),
           itemps = as.double(itemps),
           verb = as.integer(verb),
           tree = as.double(c(ncol(t2c),t(t2c))),
           hier = as.double(object$posts[m,3:ncol(object$posts)]),
           MAP = as.integer(MAP),
           sens.ngrid = as.integer(ngrid),
           sens.span = as.double(span),
           sens.Xgrid = MEgrid,

           ## output dimensions for checking NULL
           pred.n = as.integer(pred.n),
           nnprime = as.integer(nnprime),
           krige = as.integer(krige),
           bDs2x = as.integer(Ds2x),
           improv = as.integer(as.logical(improv) * nnprime),

           ## begin outputs
           Zp.mean = double(pred.n * object$n),
           ZZ.mean = double(nnprime),
           Zp.km = double(krige * pred.n * object$n),
           ZZ.km = double(krige * nnprime),
           Zp.vark = double(krige * pred.n * object$n),
           ZZ.vark = double(krige * nnprime),
           Zp.q = double(pred.n * object$n),
           ZZ.q = double(nnprime),
           Zp.s2 = double(pred.n * (zcov*object$n^2) + (!zcov)*object$n),
           ZZ.s2 = double(zcov*nnprime^2 + (!zcov)*nnprime^2),
           ZpZZ.s2 = double(pred.n * object$n * nnprime * zcov),
           Zp.ks2 = double(krige * pred.n * object$n),
           ZZ.ks2 = double(krige * nnprime),
           Zp.q1 = double(pred.n * object$n),
           Zp.med = double(pred.n * object$n),
           Zp.q2 = double(pred.n * object$n),
           ZZ.q1 = double(nnprime),
           ZZ.med = double(nnprime),
           ZZ.q2 = double(nnprime),
           Ds2x = double(Ds2x * nnprime),
           improv = double(as.logical(improv) * nnprime),
           irank = integer(as.logical(improv) * nnprime),
           ess = double(1 + itemps[1]*2),
           gpcs = double(4),
           sens.ZZ.mean = double(ngrid*object$d),
           sens.ZZ.q1 = double(ngrid*object$d),
           sens.ZZ.q2 = double(ngrid*object$d),
           sens.S = double(object$d*S*!is.null(sens.p)),
           sens.T = double(object$d*S*!is.null(sens.p)),

           ## end outputs
           PACKAGE = "tgp")

  ## post-process before returning
  ll <- tgp.postprocess(ll, Xnames, response, pred.n, zcov, Ds2x, improv, sens.p,
                        Zm0r1, object$params, TRUE)
  return(ll)
}


## tree2c
##
## converts the list-and-data.frame style tree contained in
## the tgp-class object into a C-style double-vector so that
## the C-side can start from the MAP tree contained in the object

"tree2c" <-
function(t) {

  ## change var into a numeric vector
  var <- as.character(t$var)
  var[var == "<leaf>"] <- -1
  var <- as.numeric(var)
  
  ## to return
  tr <- data.frame(rows=t$rows, var=var)
  tr <- cbind(tr, t[,8:ncol(t)])

  ## order the rows by the row column
  o <- order(tr[,1])
  tr <- tr[o,]

  return(as.matrix(tr))
}

