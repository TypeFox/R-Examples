## ---- eval=T, echo=FALSE-------------------------------------------------
require(beanz);

## ---- eval=T, echo=TRUE--------------------------------------------------

var.cov    <- c("sodium", "lvef", "any.vasodilator.use");
var.resp   <- "y";
var.trt    <- "trt";
var.censor <- "censor";
resptype   <- "survival";

subgrp.effect <- r.get.subgrp.raw(solvd.sub,
                                  var.resp   = var.resp,
                                  var.trt    = var.trt,
                                  var.cov    = var.cov,
                                  var.censor = var.censor,
                                  resptype   = resptype);
print(subgrp.effect);


## ---- eval=T, echo=TRUE--------------------------------------------------

var.estvar <- c("Estimate", "Variance");

rst.nse <- call.stan("nse", dat.sub=subgrp.effect,
                     var.estvar = var.estvar, var.cov = var.cov,
                     lst.par.pri = list(vtau=1000, vrange=c(0,0)),
                     chains=1, iter=4000, warmup=2000, thin=2, seed=1000);

rst.sr  <- call.stan("sr", dat.sub=subgrp.effect,
                     var.estvar = var.estvar, var.cov = var.cov,
                     lst.par.pri = list(vtau=1000, vgamma=1000, vrange=c(0,0)),
                     chains=1, iter=4000, warmup=2000, thin=2, seed=1000);

rst.bs  <- call.stan("bs", dat.sub=subgrp.effect,
                     var.estvar = var.estvar, var.cov = var.cov,
                     lst.par.pri = list(vtau=1000, vw=100, vrange=c(-0.1,0.1)),
                     chains=1, iter=4000, warmup=2000, thin=2, seed=1000);


## ---- eval=T, echo=TRUE, fig.width=6, fig.height=5-----------------------
sel.grps <- c(1,4,5);
tbl.sub <- r.summary.stan(rst.sr, ref.stan.rst=rst.nse, ref.sel.grps=1);
print(tbl.sub);
r.plot.stan(rst.sr, sel.grps = sel.grps, ref.stan.rst=rst.nse, ref.sel.grps=1);
r.forest.stan(rst.sr, sel.grps = sel.grps, ref.stan.rst=rst.nse, ref.sel.grps=1);

## ---- eval=T, echo=TRUE, fig.width=6, fig.height=5-----------------------
tbl.sub <- r.summary.stan(rst.bs, ref.stan.rst=rst.nse, ref.sel.grps=1);
print(tbl.sub);
r.plot.stan(rst.bs, sel.grps = sel.grps, ref.stan.rst=rst.nse, ref.sel.grps=1);
r.forest.stan(rst.bs, sel.grps = sel.grps, ref.stan.rst=rst.nse, ref.sel.grps=1);

## ---- eval=T, echo=TRUE, fig.width=6, fig.height=5-----------------------
tbl.sub <- r.summary.comp(rst.sr, sel.grps=sel.grps);
print(tbl.sub);
r.plot.stan(rst.sr, sel.grps = sel.grps);
r.forest.stan(rst.sr, sel.grps = sel.grps);

## ---- eval=T, echo=TRUE, fig.width=6, fig.height=5-----------------------
tbl.sub <- r.summary.comp(rst.bs, sel.grps=sel.grps);
print(tbl.sub);
r.plot.comp(rst.bs, sel.grps = sel.grps);
r.forest.comp(rst.bs, sel.grps = sel.grps);

## ---- echo=TRUE----------------------------------------------------------
lst.rst     <- list(nse=rst.nse, sr=rst.sr, bs=rst.bs);
tbl.summary <- r.rpt.tbl(lst.rst, dat.sub = subgrp.effect, var.cov = var.cov);
print(tbl.summary);

## ---- eval=F-------------------------------------------------------------
#  run.beanz();

## ---- echo=T-------------------------------------------------------------
gs.pval <- r.gailsimon(subgrp.effect$Estimate, subgrp.effect$Variance);
print(gs.pval);

