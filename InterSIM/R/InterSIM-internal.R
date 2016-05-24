.onLoad <-
function (libname, pkgname) 
{
    utils::data("cov.expr", "cov.M", "cov.protein", "CpG.gene.map.for.DEG", "mean.expr", "mean.expr.with.mapped.protein", "mean.M", "mean.protein", "methyl.gene.level.mean", "protein.gene.map.for.DEP", "rho.expr.protein", "rho.methyl.expr", package=pkgname, envir=parent.env(environment()))
    utils::globalVariables(c("par", "rbinom"))
}
