
##-----------------------------------------------------------------------------------
##                     comparison
##-----------------------------------------------------------------------------------

#' Comparison of posterior treatment effects
#'
#' Present the difference in the posterior treatment effects
#' between subgroups
#'
#'
#' @name comp.stan
#'
#' @return \code{r.summary.comp} generates a data frame with summary statistics
#'     of the difference of treatment effects between the selected subgroups.
#'     \code{r.plot.comp} generates the density plot of the difference in the
#'     posterior treatment effects between subgroups. \code{r.forest.comp}
#'     generates the forest plot of the difference in the posterior treatment
#'     effects between subgroups.
#'
#'
#' @examples
#'
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
#' rst.sr     <- call.stan("sr", dat.sub=subgrp.effect,
#'                         var.estvar=var.estvar, var.cov = var.cov,
#'                         lst.par.pri=list(vtau=1000, vgamma=1000, vrange=c(0,0)),
#'                         chains=1, iter=500,
#'                         warmup=100, thin=2, seed=1000);
#'
#' sel.grps <- c(1,4,5);
#' tbl.sub <- r.summary.comp(rst.sr, sel.grps=sel.grps);
#' r.plot.stan(rst.sr, sel.grps = sel.grps);
#' r.forest.stan(rst.sr, sel.grps = sel.grps);
#'
#'
#' @seealso \code{\link{call.stan}}
#'
#'
NULL


#' @rdname comp.stan
#'
#' @inheritParams r.summary.stan
#'
#' @export
#'
r.summary.comp <- function(stan.rst, sel.grps=NULL, cut=0, digits=3) {

    mus <- stan.rst$get.mus();

    sel.grps <- get.sel.subgrp(mus, sel.grps);
    if (length(sel.grps) < 2)
        return(NULL);

    mus  <- mus[,sel.grps, drop=FALSE];
    fsum <- function(x) {
        crst <- c(mean(x),
                  sd(x),
                  quantile(x, c(0.025, 0.25, 0.5, 0.75, 0.975)),
                  mean(x<cut));
        crst <- round(crst, digits);
    }

    rst <- NULL;
    for (i in 1:(ncol(mus)-1)) {
        for (j in (i+1):ncol(mus)) {
            cur.j <- sample(mus[,j], nrow(mus), TRUE);
            cur.i <- sample(mus[,i], nrow(mus), TRUE);
            rst   <- rbind(rst, fsum(cur.j-cur.i));
            rownames(rst)[nrow(rst)] <- paste("Subgroup ", sel.grps[j], "-", sel.grps[i], sep="");
        }
    }


    rst <- cbind(rownames(rst), rst);
    colnames(rst) <- c("Comparison", "Mean", "SD", "2.5%", "25%",
                       "Median", "75%", "97.5%", paste("Prob <", cut, sep=""));
    rst
}



#' @rdname comp.stan
#'
#' @inheritParams r.plot.stan
#'
#' @export
#'
r.plot.comp <- function(stan.rst, sel.grps=NULL, ...) {

    mus      <- stan.rst$get.mus();
    sel.grps <- get.sel.subgrp(mus, sel.grps);

    if (length(sel.grps) > 5 | length(sel.grps) < 2)
        return(NULL);

    mus     <- mus[,sel.grps, drop=FALSE];
    lst.den <- NULL;
    cmp.mus <- NULL;
    for (i in 1:(ncol(mus)-1)) {
        for (j in (i+1):ncol(mus)) {
            cur.j <- sample(mus[,j], nrow(mus), TRUE);
            cur.i <- sample(mus[,i], nrow(mus), TRUE);
            cmp.mus <- cbind(cmp.mus, cur.j-cur.i);
            lst.den[[length(lst.den)+1]]    <- density(cur.j - cur.i);
            names(lst.den)[length(lst.den)] <- paste("Subgroup ", sel.grps[j], "-", sel.grps[i], sep="");
        }
    }

    plot.densities(list(lst.den=lst.den, mus=cmp.mus),
                   main="Comparison of Subgroup Effects",
                   xlab="Difference of Treatment Effects",
                   ...);
}

#' @rdname comp.stan
#'
#' @inheritParams r.forest.stan
#'
#' @export
#'
r.forest.comp <- function(stan.rst, sel.grps=NULL, ..., quants=c(0.025,0.975)) {


    mus      <- stan.rst$get.mus();
    sel.grps <- get.sel.subgrp(mus, sel.grps);

    if (length(sel.grps) > 5 | length(sel.grps) < 2)
        return(NULL);

    mus      <- mus[,sel.grps];
    mu.qmean <- NULL;
    for (i in 1:(ncol(mus)-1)) {
        for (j in (i+1):ncol(mus)) {
            cur.j <- sample(mus[,j], nrow(mus), TRUE);
            cur.i <- sample(mus[,i], nrow(mus), TRUE);
            cur.m <- cur.j-cur.i;
            cq    <- quantile(cur.m, quants);

            mu.qmean <- rbind(mu.qmean, c(cq[1], mean(cq), cq[2]));
            rownames(mu.qmean)[nrow(mu.qmean)] <- paste("Subgroup ", sel.grps[j], "-", sel.grps[i], sep="");
        }
    }

    plot.forest(mu.qmean, main="Difference of Subgroup Effects Forest Plot", ...);
}


##-----------------------------------------------------------------------------------
##                     posterior summary
##-----------------------------------------------------------------------------------
#' Posterior subgroup treatment effects
#'
#' Present the posterior subgroup treatment effects
#'
#' @name summary.stan
#'
#' @return \code{r.summary.stan} generates a dataframe with summary statistics
#'     of the posterior treatment effect for the selected subgroups.
#'     \code{r.plot.stan} generates the density plot of the posterior treatment
#'     effects for the selected subgroups. \code{r.forest.stan}
#'     generates the forest plot of the posterior treatment
#'     effects.
#'
#'@examples
#' \dontrun{
#' sel.grps <- c(1,4,5);
#' tbl.sub <- r.summary.stan(rst.sr, ref.stan.rst=rst.nse, ref.sel.grps=1);
#' r.plot.stan(rst.sr, sel.grps = sel.grps, ref.stan.rst=rst.nse, ref.sel.grps=1);
#' r.forest.stan(rst.sr, sel.grps = sel.grps, ref.stan.rst=rst.nse, ref.sel.grps=1);}
#'
#' @seealso \code{\link{call.stan}}
#'
#'
NULL


#' @rdname summary.stan
#'
#' @param cut cut point to compute the probabiliby that the posterior subgroup
#'     treatment effects is below
#'
#' @param digits number of digits in the summary result table
#'
#' @inheritParams r.plot.stan
#'
#' @export
#'
r.summary.stan <- function(stan.rst, sel.grps=NULL, ref.stan.rst=NULL, ref.sel.grps=1, cut=0, digits=3) {

    s.mus <- get.all.mus(stan.rst, sel.grps, ref.stan.rst, ref.sel.grps);

    rst <- apply(s.mus, 2,
                 function(x) {
                     crst <- c(mean(x),
                               sd(x),
                               quantile(x, c(0.025, 0.25, 0.5, 0.75, 0.975)),
                               mean(x<cut)
                               );
                     crst <- round(crst, digits);
                 });

    rst <- t(rst);
    rst <- cbind(colnames(s.mus), rst);

    colnames(rst) <- c("Subgroup", "Mean", "SD", "2.5%", "25%", "Median", "75%",
                       "97.5%", paste("Prob < ", cut, sep=""));
    rst
}

#' @rdname summary.stan
#'
#' @param stan.rst a class \code{beanz.stan} object generated by
#'     \code{\link{call.stan}}
#'
#' @param sel.grps an array of subgroup numbers to be included in the summary results
#'
#' @param ref.stan.rst a class \code{beanz.stan} object from \code{\link{call.stan}} that
#'     is used as the reference
#'
#' @param ref.sel.grps subgroups from the reference model to be included in the
#'     summary table
#'
#' @param ... options for \code{plot} function
#'
#' @export
#'
r.plot.stan <- function(stan.rst, sel.grps=NULL,
                        ref.stan.rst=NULL, ref.sel.grps=1,
                        ... ) {

    s.mus          <- get.all.mus(stan.rst, sel.grps, ref.stan.rst, ref.sel.grps);
    lst.den        <- apply(s.mus, 2, density);
    names(lst.den) <- colnames(s.mus);

    plot.densities(list(mus=s.mus,lst.den=lst.den),
                   main="Posterior Distribution of Subgroup Effects",
                   xlab="Treatment Effect",
                   ...);
}

#' @rdname summary.stan
#' @inheritParams r.plot.stan
#' @param quants lower and upper quantiles of the credible intervals in the
#'     forest plot
#'
#' @export
#'
r.forest.stan <- function(stan.rst, sel.grps=NULL,
                        ref.stan.rst=NULL, ref.sel.grps=1,
                        ..., quants=c(0.025,0.975)) {

    s.mus    <- get.all.mus(stan.rst, sel.grps, ref.stan.rst, ref.sel.grps);
    mu.q     <- apply(s.mus, 2, quantile, quants);
    mu.mean  <- apply(s.mus, 2, mean);
    mu.qmean <- cbind(mu.q[1,], mu.mean, mu.q[2,]);
    rownames(mu.qmean) <- colnames(s.mus);

    plot.forest(mu.qmean, main="Subgroup Effects Forest Plot", ...);
}



##-----------------------------------------------------------------------------------
##                     traceplot
##-----------------------------------------------------------------------------------
#' Trace plot of \code{rstan} samples
#'
#' Trace plot of \code{rstan} samples for checking the convergence of the MCMC chains
#'
#' @param n.eachrow number of trace plot each row
#'
#' @param width plot size: width
#'
#' @param height plot size: height
#'
#' @inheritParams r.plot.stan
#'
#' @seealso \code{\link{call.stan}}
#'
#'
#' @export
#'
r.plot.trace <- function(stan.rst, n.eachrow=2, width=4, height=0.8) {

    stopifnot(class(stan.rst) == "beanz.stan");

    stan.samples <- stan.rst$smps;

    pars  <- dimnames(stan.samples)$parameters;
    n.par <- length(pars);
    n.row <- ceiling(n.par / n.eachrow);
    n.smp <- dim(stan.samples)[1];

    par(mfrow=c(n.row, n.eachrow), mar=c(0,0,2,0), pin=c(width, height));
    for (i in 1:n.par) {
        smps <- stan.samples[,1,i];
        plot(1:n.smp, smps, main=pars[i], type="l", xlab=NULL, ylab=NULL, axes=FALSE);
        box();
    }
}

##-----------------------------------------------------------------------------------
##                     SUMMARY REPORT
##-----------------------------------------------------------------------------------

#' Summary table of treatment effects
#'
#' Compare the DIC from different models and report the summary of treatment effects
#' based on the model with the smallest DIC value
#'
#' @param lst.stan.rst list of class \code{beanz.stan} results from
#'     \code{\link{call.stan}} for different models
#'
#' @inheritParams call.stan
#' @inheritParams r.summary.stan
#'
#' @return A dataframe with summary statistics of the model selected by DIC
#'
#' @export
#'
r.rpt.tbl <- function(lst.stan.rst, dat.sub, var.cov, cut=0, digits=3) {

    if (is.null(dat.sub) | is.null(lst.stan.rst))
        return(NULL);

    dic <- NULL;
    for (i in 1:length(lst.stan.rst)) {
        dic <- c(dic, lst.stan.rst[[i]]$dic);
    }

    min.mdl <- which.min(dic);
    mus     <- lst.stan.rst[[min.mdl]]$get.mus();
    rst     <- apply(mus, 2,
                     function(x) {
                         crst <- c(mean(x),
                                   sd(x),
                                   mean(x < cut)
                                   );
                         crst <- round(crst, digits);
                     });

    rst           <- t(rst);
    colnames(rst) <- c("Mean", "SD", paste("Prob < ", cut, sep=""));
    rst           <- cbind(dat.sub[, c("Subgroup", var.cov)] ,rst);

    rst
}

##-----------------------------------------------------------------------------------
##                             Frequentist GailSimon Test
##-----------------------------------------------------------------------------------
#' Gail-Simon Test
#'
#' Gail-Simon qualitative interaction test.
#'
#' @param effects subgroup treatment effects
#'
#' @param sderr standard deviation of the estimated treatment effects
#'
#' @param d clinically meaningful difference
#'
#'
#' @examples
#' \dontrun{
#' var.cov    <- c("sodium", "lvef", "any.vasodilator.use");
#' var.resp   <- "y";
#' var.trt    <- "trt";
#' var.censor <- "censor";
#' resptype   <- "survival";
#' subgrp.effect <- r.get.subgrp.raw(solvd.sub,
#'                                   var.resp   = var.resp,
#'                                   var.trt    = var.trt,
#'                                   var.cov    = var.cov,
#'                                   var.censor = var.censor,
#'                                   resptype   = resptype);
#'
#' gs.pval <- r.gailsimon(subgrp.effect$Estimate,
#'                        subgrp.effect$Variance); }
#'
#'
#' @export
#'
r.gailsimon <- function(effects, sderr, d=0) {
    d  <- abs(d);
    I  <- length(effects);
    Qm <- sum((effects > d) * ((effects-d)/sderr)^2);
    Qp <- sum((effects < -d) * ((effects+d)/sderr)^2);
    test.stats <- min(Qm, Qp);
    pval       <- sum(dbinom(1:(I-1), I-1, .5) * (1 - pchisq(test.stats, df=1:(I-1))));
    pval
}



##------------------------------------------------------------------
##------------tool kit----------------------------------------------
##------------------------------------------------------------------

plot.densities <- function(lst.mus.den,
                           main="", xlab="",
                           ylims=NULL,
                           xlims=NULL,
                           cols=c(rep(c("black", "red", "green", "brown", "gray", "yellow",
                                        "cyan", "blue", "ivory", "wheat"), 4),
                                  colors()),
                           ltys=c(rep(1:4, each=10),
                                  rep(1:6, each=100)),
                           quants=c(0.025,0.975)) {

    mus     <- lst.mus.den$mus;
    lst.den <- lst.mus.den$lst.den;

    ##get min max of x and y in dens
    f.tmp <- function(cur.den, cur.mu) {
        max.y <- max(cur.den$y);
        min.x <- quantile(cur.mu, quants[1]);
        max.x <- quantile(cur.mu, quants[2]);
        rst   <- c(max.y, min.x, max.x);
    }

    ##boundary
    lims <- NULL;
    for (i in 1:length(lst.den)) {
        cur.lim <- f.tmp(lst.den[[i]], mus[,i]);
        lims    <- cbind(lims, cur.lim);
    }

    if (is.null(ylims)) {
        max.y <- max(lims[1,]);
        ylims <- c(0,1.1*max.y);
    }

    if (is.null(xlims)) {
        min.x <- min(lims[2,]);
        max.x <- max(lims[3,]);
        xlims <- c(min.x,max.x)
    }

    plot(NULL, ylim=ylims, xlim=xlims,
         main=main,
         xlab=xlab,
         ylab="Density");

    ns <- length(lst.den);
    for (i in 1:ns) {
        lines(lst.den[[i]]$x,
              lst.den[[i]]$y,
              col=cols[i],  lty=ltys[i], lwd=2);
    }

    legend("topleft",
           legend=names(lst.den),
           lty=ltys[1:ns],
           col=cols[1:ns], bty="n", cex=1.2);
}




##plot forest plot
plot.forest <- function(quants,
                        cut=0,
                        y.bottom=0.1, y.top=0.95,
                        x.labs=0.2,
                        main="", xlab="") {

    ##where to put labels
    ng    <- nrow(quants);
    ys    <- seq(y.top, y.bottom, length.out=ng);
    lab.g <- rownames(quants);

    ##convert
    xrs <- range(c(quants, cut));
    xrs <- range(xrs, 0.1*(xrs[2]-xrs[1]), -0.1*(xrs[2]-xrs[1]));

    xs <- c(cut, seq(xrs[1], xrs[2], length.out=4));
    xs <- round(xs, digits=1);

    f.con <- function(x) {
        (x - xrs[1])/(xrs[2] - xrs[1])*(1-x.labs) + x.labs;
    }

    ##plot
    plot(NULL, ylim=c(0,1), xlim=c(0,1), main=main, xlab=xlab, ylab="", axes=FALSE);
    lines(f.con(c(cut,cut)), c(0,1), lty=2, lwd=2, col="gray");
    axis(1, at=f.con(xs), labels=xs);

    for (i in 1:nrow(quants)) {
        cur.q <- quants[i,];
        text(x.labs, ys[i], lab.g[i], pos=2);
        lines(f.con(cur.q[c(1,3)]), c(ys[i], ys[i]), lwd=2, lty=1, col="gray");
        points(f.con(cur.q[2]), ys[i], pch=23);
    }
}

get.sel.subgrp <- function(mus, sel.grps) {
    if (is.null(sel.grps)) {
        sel.grps <- 1:ncol(mus);
    } else {
        sel.grps <- as.numeric(sel.grps);
    }

    sel.grps;
}


##comibine all mus for presentation
get.all.mus <- function(stan.rst, sel.grps=NULL, ref.stan.rst=NULL, ref.sel.grps=1) {

    stopifnot(class(stan.rst) == "beanz.stan");
    stopifnot(is.null(ref.stan.rst) | class(ref.stan.rst) == "beanz.stan");

    mus      <- stan.rst$get.mus();
    sel.grps <- get.sel.subgrp(mus, sel.grps);
    s.mus    <- mus[, sel.grps, drop=FALSE];

    if (!is.null(ref.stan.rst)) {
        ref.mus <- ref.stan.rst$get.mus();
        if (is.null(ref.sel.grps))
            ref.sel.grps <- sel.grps;
        r.mus <- ref.mus[, ref.sel.grps, drop=FALSE];
        colnames(r.mus) <- paste(ref.stan.rst$mdl, "(", ref.sel.grps, ")", sep="");
        s.mus <- cbind(s.mus, r.mus);
    }

    s.mus
}
