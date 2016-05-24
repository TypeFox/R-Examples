epcreg.baselearner.control <- function(baselearners=c("nnet","rf","svm","gbm","knn")
  , baselearner.configs=make.configs(baselearners, type="regression"), npart=1, nfold=5) {
    return (list(configs=baselearner.configs, npart=npart, nfold=nfold))
}

epcreg.integrator.control <- function(errfun=rmse.error, nfold=5, method=c("default")) {
  return (list(errfun=rmse.error, nfold=nfold, method=method))
}

epcreg.set.filemethod <- function(formula, data, instance.list, type="regression") FALSE # TODO: to be implemented

epcreg <- function(formula, data
  , baselearner.control=epcreg.baselearner.control()
  , integrator.control=epcreg.integrator.control()
  , ncores=1, filemethod=FALSE, print.level=1
  , preschedule = FALSE
  , schedule.method = c("random", "as.is", "task.length"), task.length) {
  if (integrator.control$method!="default") stop("invalid PCR integration method")
  ncores.max <- try(detectCores(),silent=T)
  mycall <- match.call()
  if (!inherits(ncores.max,"try-error")) ncores <- min(ncores,ncores.max)
  if (print.level>=1 && ncores>1) cat("running in parallel mode, using", ncores, "cores\n")

  # training a batch of baselearners
  partitions.bl <- generate.partitions(baselearner.control$npart, nrow(data), baselearner.control$nfold)
  my.instance.list <- make.instances(baselearner.control$configs, partitions.bl)

  # TODO: determining of base learner estimation objects must be saved to disk or not
  if (missing(filemethod)) filemethod <- epcreg.set.filemethod(formula, data, my.instance.list)
  
  if (print.level>=1) cat("CV training of base learners...\n")
  est.baselearner.cv.batch <- Regression.CV.Batch.Fit(my.instance.list, formula, data, ncores=ncores, filemethod=filemethod, print.level=print.level, preschedule = preschedule, schedule.method = schedule.method, task.length = task.length)
  if (print.level>=1) cat("finished CV training of base learners\n")
  Xcv <- est.baselearner.cv.batch@pred
  y <- data[,all.vars(formula)[1]] # TODO: more robust way of extracting y (here and inside base learner functions)
  
  partition.int <- generate.partition(nrow(data), integrator.control$nfold)
  my.integrator.config <- Regression.Integrator.PCR.SelMin.Config(errfun=integrator.control$errfun, partition=partition.int)
  est.integrator <- Regression.Integrator.Fit(my.integrator.config, X=Xcv, y=y, print.level=print.level)
  pred <- est.integrator@pred
  
  ret <- list(call=mycall, formula=formula, instance.list=my.instance.list, integrator.config=my.integrator.config, method=integrator.control$method
    , est=list(baselearner.cv.batch=est.baselearner.cv.batch, integrator=est.integrator)
    , y=y, pred=pred, filemethod=filemethod)
  class(ret) <- "epcreg"
  if (filemethod) class(ret) <- c(class(ret), "epcreg.file")
  return (ret)
}

predict.epcreg <- function(object, newdata=NULL, ncores=1, preschedule = TRUE, ...) {
  if (is.null(newdata)) return (object$pred)
  if (object$method=="default") {
    newpred.baselearner.cv.batch <- predict(object$est$baselearner.cv.batch, newdata, ncores=ncores, preschedule = preschedule, ...)
    newpred <- predict(object$est$integrator, Xnew=newpred.baselearner.cv.batch, config.list=object$est$baselearner.batch@config.list, ...)
  } else {
    stop("invalid PCR integration method")
  }
  return (as.numeric(newpred))
}

plot.epcreg <- function(x, ...) {
  errfun <- x$integrator.config@errfun 
  error <- errfun(x$pred, x$y)
  oldpar <- par(mfrow=c(1,2))
  plot(x$est$baselearner.cv.batch, errfun=errfun, ylim.adj = error)
  abline(h=error, lty=2)
  pcr.errors <- x$est$integrator@est$select@est$error
  #pcr.pred <- x$est$integrator@est$pcr@pred
  #pcr.errors <- apply(pcr.pred, 2, errfun, x$y)
  plot(pcr.errors, type="l", xlab="Number of Principal Components", ylab="CV Error", main="Integrator Performance")
  par(oldpar)
}

epcreg.save <- function(obj, file) {
  if (!identical(class(obj),c("epcreg","epcreg.file"))) stop("invalid object class (must be epcreg & epcreg.file)")
  if (missing(file)) stop("must provide file argument")
  tmpfiles <- obj$est$baselearner.cv.batch@tmpfiles
  tmpfile.new <- tempfile()
  save(obj, file=tmpfile.new, compress=F)
  all.files <- c(tmpfile.new, tmpfiles)
  all.files.basename <- basename(all.files)
  
  tmpdir <- paste0("./.", basename(tempfile("dir")),"/")
  dir.create(tmpdir)
  all.files.new <- paste0(tmpdir, all.files.basename)
  file.copy(all.files, all.files.new)
  meta <- list(filename.mainobj=all.files.basename[1], filenames.batchobj=all.files.basename[1+1:length(tmpfiles)])
  save(meta, file=paste0(tmpdir, "meta"), compress=FALSE)

  tar(file, files=tmpdir, compression="gzip")
  unlink(tmpdir, recursive=TRUE)
}

epcreg.load <- function(file) {
  filepaths <- untar(file, list=T)
  basenames <- basename(filepaths)
  dirnames <- dirname(filepaths)
  if (length(unique(dirnames))>1) stop("unexpected multiple directories in tar filepaths")
  metafile.index <- which(basenames=="meta")
  
  extdir <- dirnames[1] # this is where untar will extract the files to
  untar(file)
  meta <- NULL # to overcome codetools error: "no visible binding for meta"
  load(filepaths[metafile.index]) # this will load "meta"
  mainfile.index <- which(basenames==meta$filename.mainobj)
  load(filepaths[mainfile.index]) # this will load "obj"
  if (!identical(class(obj),c("epcreg","epcreg.file"))) stop("invalid object class (must be epcreg & epcreg.file)")
  
  basenames.ordered <- basename(obj$est$baselearner.cv.batch@tmpfiles)
  #if (!identical(sort(basenames.ordered),sort(basenames[-c(metafile.index,mainfile.index)]))) stop("basenames mismatch")
  filepaths.ordered <- paste(extdir, basenames.ordered, sep="/")
  
  # copy batch files to new tempfiles in tempdir
  tmpfiles.new <- tempfile(rep("file", length(filepaths.ordered)))
  # replaced file.rename with file.copy and unlink to handle cross-device cases (where . and R tmp directories are on different devices)
  file.copy(from=filepaths.ordered, to=tmpfiles.new)
  unlink(filepaths.ordered, recursive=TRUE)
  
  unlink(extdir, recursive=TRUE)
  
  obj$est$baselearner.cv.batch@tmpfiles <- tmpfiles.new
  n.instance <- length(obj$est$baselearner.cv.batch@instance.list@instances)
  for (i in 1:n.instance) {
    partid <- obj$est$baselearner.cv.batch@instance.list@instances[[1]]@partid
    nfold <- length(unique(obj$est$baselearner.cv.batch@instance.list@partitions[,partid]))
    for (j in 1:nfold) {
      obj$est$baselearner.cv.batch@fitobj.list[[i]]@fitobj.list[[j]]@est <- tmpfiles.new[obj$est$baselearner.cv.batch@tmpfiles.index.list$start[i]+j-1]
    }
  }

  return (obj)
}

print.epcreg <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
}

summary.epcreg <- function(object, ...) {
  #summary.baselearner <- NULL
  #summary.integrator <- summary(object$est$integrator)
  #summary.integrator <- NULL
  #ret <- list(baselearner=summary.baselearner, integrator=summary.integrator)
  
  n.instance <- length(object$instance.list@instances)
  maxpc <- object$est$integrator@est$pcr@sweep.list[[1]]@config@n
  index.min <- object$est$integrator@est$select@est$index.min
  error.min <- object$est$integrator@est$select@est$error.min
  tvec <- object$est$baselearner.cv.batch@tvec
  
  ret <- list(n.instance=n.instance, maxpc=maxpc, index.min=index.min, error.min=error.min, tvec = tvec)
  
  class(ret) <- "summary.epcreg"
  return (ret)
}

print.summary.epcreg <- function(x, ...) {
  cat("number of base learner instances:", x$n.instance, "\n")
  cat("maximum number of PC's considered:", x$maxpc, "\n")
  cat("optimal number of PC's:", x$index.min, "\n")
  cat("minimum error:", x$error.min, "\n")
}


#  ret <- list(call=mycall, formula=formula, instance.list=my.instance.list, integrator.config=my.integrator.config, method=integrator.control$method
#    , est=list(baselearner.cv.batch=est.baselearner.cv.batch, integrator=est.integrator)
#    , y=y, pred=pred, filemethod=filemethod)



