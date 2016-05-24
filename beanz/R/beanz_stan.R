##----------------------------------------------------------------
##                  STAN MODEL
##----------------------------------------------------------------
#' Call STAN models
#'
#' Call STAN to draw posterior samples for Bayesian HTE models.
#'
#' @param mdls name of the Bayesian HTE model. The options are:
#'
#' \describe{
#'   \item{nse}{No subgroup effect model}
#'   \item{fs}{Full stratification model}
#'   \item{sr}{Simple regression model}
#'   \item{bs}{Basic shrinkage model}
#'   \item{srs}{Simple regression with shrinkage model}
#'   \item{ds}{Dixon-Simon model}
#'   \item{eds}{Extended Dixon-Simon model}
#' }
#'
#' @param dat.sub dataset with subgroup treatment effect summary data
#'
#' @param var.estvar column names in dat.sub that corresponds to treatment effect
#'     estimation and the estimated variance
#'
#' @param var.cov array of column names in dat.sub that corresponds to binary or
#'     ordinal baseline covaraites
#'
#' @param var.nom array of column names in dat.sub that corresponds to nominal
#'     baseline covariates
#'
#' @param lst.par.pri list of prior parameters for each model:
#'
#' \describe{
#'   \item{nse, fs}{\code{vtau}, \code{vrange} for (\eqn{\Delta_1, \Delta_2})}
#'   \item{sr}{\code{vtau}, \code{vgamma}, \code{vrange}}
#'   \item{bs, ds, eds}{\code{vtau}, \code{vw}, \code{vrange}}
#'   \item{srs}{\code{vtau}, \code{vw}, \code{vgamma}, \code{vrange}}
#' }
#'
#' @param ... options to call STAN sampling. These options include
#'     \code{chains}, \code{iter}, \code{warmup}, \code{thin}, \code{algorithm}.
#'     See \code{rstan::sampling} for details.
#'
#' @return A class \code{beanz.stan} list containing
#'  \describe{
#'    \item{mdl}{name of the Bayesian HTE model}
#'    \item{stan.rst}{raw \code{rstan} \code{sampling} results}
#'    \item{smps}{matrix of the posterior samples}
#'    \item{get.mus}{method to return the posterior sample of the subgroup treatment effects}
#'    \item{DIC}{DIC value}
#' }
#'
#' @examples
#' \dontrun{
#' var.cov    <- c("sodium", "lvef", "any.vasodilator.use");
#' var.resp   <- "y";
#' var.trt    <- "trt";
#' var.censor <- "censor";
#' resptype   <- "survival";
#' var.estvar <- c("Estimate", "Variance");
#'
#' subgrp.effect <- r.get.subgrp.raw(solvd.sub,
#'                                   var.resp   = var.resp,
#'                                   var.trt    = var.trt,
#'                                   var.cov    = var.cov,
#'                                   var.censor = var.censor,
#'                                   resptype   = resptype);
#'
#' rst.nse    <- call.stan("nse", dat.sub=subgrp.effect,
#'                          var.estvar = var.estvar, var.cov = var.cov,
#'                          lst.par.pri = list(vtau=1000, vrange=c(0,0)),
#'                          chains=1, iter=4000,
#'                          warmup=2000, thin=2, seed=1000);
#' rst.sr     <- call.stan("sr", dat.sub=subgrp.effect,
#'                         var.estvar=var.estvar, var.cov = var.cov,
#'                         lst.par.pri=list(vtau=1000, vgamma=1000, vrange=c(0,0)),
#'                         chains=1, iter=4000,
#'                         warmup=2000, thin=2, seed=1000)};
#' @export
#'
#'
#'
call.stan <- function(mdls = c("nse", "fs", "sr", "bs", "srs", "ds", "eds"),
                      dat.sub,
                      var.estvar,
                      var.cov,
                      lst.par.pri,
                      var.nom=NULL,
                      ...) {

    mdls <- match.arg(mdls);



    ##data manipulation
    eval(parse(text=paste("vname <- stan.model.", mdls,
               "(dat.sub, var.estvar, var.cov, var.nom)", sep="")));

    ##add prior data
    names(lst.par.pri) <- toupper(names(lst.par.pri));
    vname              <- c(vname, lst.par.pri);

    ##call stan
    stan.rst <- sampling(stanmodels[[mdls]], data=vname, ...);

    smps <- extract(stan.rst, permuted=FALSE);
    dic  <- get.dic(dat.sub[,var.estvar[1]], smps);

    ##return
    rst <- list(mdl      = get.mdl.name(mdls),
                stan.rst = stan.rst,
                smps     = smps,
                dic      = dic,
                get.mus=function() {

        n.mu <- stan.rst@par_dims$mu;
        mj   <- paste("mu[", 1:n.mu, "]", sep="");
        mus  <- smps[,1,mj];
        colnames(mus) <- paste("Subgroup", 1:n.mu);
        mus
    })

    class(rst) <- "beanz.stan";
    rst
}

##no subgroup effect model
stan.model.nse <- function(dat.sub, var.estvar, var.cov, var.nom) {
    ##data manipulation
    SIZE   <- nrow(dat.sub);
    Y      <- dat.sub[,var.estvar[1]];
    VY     <- dat.sub[,var.estvar[2]];
    vname  <- make.list(environment());
    vname
}

##full stratification
stan.model.fs <- function(dat.sub, var.estvar, var.cov, var.nom) {
    ##data manipulation
    SIZE  <- nrow(dat.sub);
    Y     <- dat.sub[,var.estvar[1]];
    VY    <- dat.sub[,var.estvar[2]];
    vname <- make.list(environment());
    vname
}

##simple regression
stan.model.sr <- function(dat.sub, var.estvar, var.cov, var.nom) {
    ##data manipulation
    Y      <- dat.sub[,var.estvar[1]];
    SIZE   <- nrow(dat.sub);
    VY     <- dat.sub[,var.estvar[2]];

    ##design matrix
    dx     <- dat.sub[,var.cov, drop=FALSE];
    dx     <- df.convert(dx, var.nom);
    fml    <- paste("~", paste(var.cov, collapse="+"), sep="");
    des.x  <- model.matrix(formula(fml), dx);
    X      <- des.x[,-1]; #remove intercept
    NX     <- ncol(X);

    vname  <- make.list(environment());
    vname
}

##basic shrinkage
stan.model.bs <- stan.model.fs;
##simple regression + shrinkage
stan.model.srs <- stan.model.sr;
##Dixon and Simon
stan.model.ds <- stan.model.sr;

##extended dixon and simon
stan.model.eds <- function(dat.sub, var.estvar, var.cov, var.nom) {
    ##data manipulation
    SIZE   <- nrow(dat.sub);
    Y      <- dat.sub[,var.estvar[1]];
    VY     <- dat.sub[,var.estvar[2]];

    ##design matrix
    dx     <- dat.sub[,var.cov, drop=FALSE];
    dx     <- df.convert(dx, var.nom);
    NTAU   <- ncol(dx);

    fml    <- paste("~", paste(var.cov, collapse="+"), "^", NTAU, sep="");
    des.x  <- model.matrix(formula(fml), dx);
    X      <- des.x[,-1]; #remove intercept
    NX     <- ncol(X);

    ##use the number of : for the order of interaction
    d.name <- colnames(X);
    TAUINX <- 1 + nchar(d.name) - nchar(gsub(":", "", d.name));

    ##return
    vname <- make.list(environment());
    vname
}

##convert dataframe into numeric values
##except the var.nom (nominal covariates)
df.convert <- function(dat, var.nom=NULL) {
    rst <- as.data.frame(lapply(dat, as.numeric));
    if (!is.null(var.nom)) {
        for (i in 1:length(var.nom)) {
            rst[,var.nom[i]] <- as.character(rst[,var.nom[i]]);
        }
    }
    rst
}

make.list <- function(par.env) {
    objlst  <- objects(par.env);
    rst     <- list();
    for (o in 1:length(objlst)) {
        curo <- objlst[o];
        if (toupper(curo) == curo) {
            rst[[curo]] <- get(curo, envir=par.env);
        }
    }
    rst
}

##-----------------------------------------------------------------------------------
##                             DIC
##-----------------------------------------------------------------------------------
get.log.ll <- function(y, mu, vs) {
    rst <- dnorm(y, mu, sqrt(vs), log=TRUE);
    sum(rst);
}

get.dic <- function(y, stan.smps) {
    ny    <- length(y);
    l.mus <- paste("mu[", 1:ny, "]", sep="");
    l.vs  <- paste("vs[", 1:ny, "]", sep="");

    mus   <- stan.smps[, 1, l.mus];
    vs    <- stan.smps[, 1, l.vs];

    mean.mu <- apply(mus, 2, mean);
    mean.vs <- 1/(apply(1/vs, 2, mean));
    d.theta.bar <- -2*get.log.ll(y, mean.mu, mean.vs);

    all.d   <- apply(cbind(mus, vs), 1,
                     function(x) {
                         rst <- -2*get.log.ll(y, x[1:ny], x[(ny+1):(2*ny)]);
                         rst
                     });
    bar.d   <- mean(all.d);

    ##dic
    DIC <- 2*bar.d - d.theta.bar;
    DIC
}

get.mdl.name <- function(mdl) {

    ALL.MODELS  <- c("No subgroup effect",
                     "Full stratification",
                     "Simple regression",
                     "Basic shrinkage",
                     "Regression and shrinkage",
                     "Dixon and Simon",
                     "Extended Dixon and Simon"
                     );
    STAN.NAME   <- c("nse", "fs", "sr", "bs", "srs", "ds", "eds");

    rst <- ALL.MODELS[which(mdl == STAN.NAME)];
}
