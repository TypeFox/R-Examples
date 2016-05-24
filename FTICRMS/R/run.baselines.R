`run.baselines` <-
function(root.dir = ".", raw.dir, base.dir, overwrite = FALSE, use.par.file = FALSE,
        par.file = "parameters.RData", sm.par = 1e-11, sm.ord = 2, max.iter = 20,
        tol = 5e-8, sm.div = NA, sm.norm.by = c("baseline", "overestimate", "constant"),
        neg.div = NA, neg.norm.by = c("baseline", "overestimate", "constant"),
        rel.conv.crit = TRUE, zero.rm = TRUE, halve.search = FALSE){
    fail <- 0
    if(missing(base.dir)){base.dir <- paste(root.dir, "/Baselines", sep="")}
    if(missing(raw.dir)){raw.dir <- paste(root.dir, "/Raw_Data", sep="")}
    if(use.par.file){
        load(paste(root.dir, "/", par.file, sep=""))
        tmp <- match.call()
        tmp[[1]] <- as.name("list")
        tmp <- eval(tmp)
        if(length(tmp) > 0){
            for(i in 1:length(tmp)){
                assign(names(tmp)[i],tmp[[i]])
            }
        }
    }
    neg.norm.by <- match.arg(neg.norm.by)
    sm.norm.by <- match.arg(sm.norm.by)
    if(!file.exists(base.dir)){
        dir.create(base.dir)
    }
    for(i in list.files(raw.dir)){
        if(!file.exists(paste(base.dir, "/", sub("\\.txt$", ".RData", i), sep="")) || 
                overwrite){
            if(regexpr(",", readLines(paste(raw.dir, "/", i, sep=""), n=1)) != -1){ 
                spect <- read.csv(paste(raw.dir, "/", i, sep=""), header=FALSE)
            } else {
                spect <- read.table(paste(raw.dir, "/", i, sep=""), header=FALSE)
            }
            spect.base <- baseline(spect[,2], sm.par=sm.par, sm.ord=sm.ord, max.iter=max.iter,
                tol=tol, neg.div=neg.div, sm.div=sm.div, sm.norm.by=sm.norm.by, neg.norm.by=neg.norm.by,
                rel.conv.crit=rel.conv.crit, zero.rm=zero.rm, halve.search=halve.search)[[1]]
            names(spect) <- c("Freq", "Amp")
            save(spect,spect.base, file=sub("\\.txt$",".RData", paste(base.dir, "/", i, sep="")))
            rm(spect,spect.base)
        } else {
            fail <- fail + 1
        }
    }
    if(fail) {
        warning(paste(fail, "baseline file(s) already existed and overwrite = FALSE; those file(s) not updated"))
    }
}

