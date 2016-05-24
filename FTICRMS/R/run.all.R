`run.all` <-
function(par.file = "parameters.RData", root.dir = "."){
    run.baselines(root.dir=root.dir, par.file=par.file, use.par.file=TRUE)
    run.peaks(root.dir=root.dir, par.file=par.file, use.par.file=TRUE)
    run.lrg.peaks(root.dir=root.dir, par.file=par.file, use.par.file=TRUE)
    run.strong.peaks(root.dir=root.dir, par.file=par.file, use.par.file=TRUE)
    run.cluster.matrix(root.dir=root.dir, par.file=par.file, use.par.file=TRUE)
    run.analysis(root.dir=root.dir, par.file=par.file, use.par.file=TRUE)
}
