ecv.regression.baselearner.control <- function(baselearners=c("nnet","rf","svm","gbm","knn")
  , baselearner.configs=make.configs(baselearners, type="regression"), npart=1, nfold=5) {
    return (list(configs=baselearner.configs, npart=npart, nfold=nfold))
}

ecv.regression.integrator.control <- function(errfun=rmse.error, method=c("default")) {
  return (list(errfun=rmse.error, method=method))
}

ecv.set.filemethod <- function(formula, data, instance.list, type="regression") FALSE # TODO: to be implemented

ecv.regression <- function(formula, data
  , baselearner.control=ecv.regression.baselearner.control()
  , integrator.control=ecv.regression.integrator.control()
  , ncores=1, filemethod=FALSE, print.level=1) {
  if (integrator.control$method!="default") stop("invalid CV integration method")
  ncores.max <- try(detectCores(),silent=T)
  mycall <- match.call()
  if (!inherits(ncores.max,"try-error")) ncores <- min(ncores,ncores.max)
  if (print.level>=1 && ncores>1) cat("running in parallel mode, using", ncores, "cores\n")

  # training a batch of baselearners
  partitions.bl <- generate.partitions(baselearner.control$npart, nrow(data), baselearner.control$nfold)
  my.instance.list <- make.instances(baselearner.control$configs, partitions.bl)

  # TODO: determining of base learner estimation objects must be saved to disk or not
  if (missing(filemethod)) filemethod <- ecv.set.filemethod(formula, data, my.instance.list)
  
  if (print.level>=1) cat("CV training of base learners...\n")
  est.baselearner.cv.batch <- Regression.CV.Batch.Fit(my.instance.list, formula, data, ncores=ncores, filemethod=filemethod, print.level=print.level)
  if (print.level>=1) cat("finished CV training of base learners\n")
  Xcv <- est.baselearner.cv.batch@pred
  y <- data[,all.vars(formula)[1]] # TODO: more robust way of extracting y (here and inside base learner functions)
  
  # trainig on full dataset is needed for a subset of methods
  if (print.level>=1) cat("full training of base learners...\n")
  #Regression.Batch.Fit <- function(config.list, formula, data, ncores=1, filemethod=FALSE, print.level=1)
  est.baselearner.batch <- Regression.Batch.Fit(baselearner.control$configs, formula, data, ncores=ncores, filemethod=filemethod, print.level=print.level)
  if (print.level>=1) cat("finished full training of base learners\n")
  Xfull <- est.baselearner.batch@pred

  my.integrator.config <- Regression.Select.MinAvgErr.Config(errfun=integrator.control$errfun, instance.list=my.instance.list)
  est.integrator <- Regression.Select.Fit(my.integrator.config, X=Xcv, y=y, print.level=print.level)
  #pred <- predict(est.integrator, Xnew=Xfull, config.list=baselearner.control$configs)
  pred <- est.integrator@pred
  
  ret <- list(call=mycall, formula=formula, instance.list=my.instance.list, integrator.config=my.integrator.config, method=integrator.control$method
    , est=list(baselearner.cv.batch=est.baselearner.cv.batch, baselearner.batch=est.baselearner.batch, integrator=est.integrator)
    , y=y, pred=pred, filemethod=filemethod)
  class(ret) <- "ecv.regression"
  if (filemethod) class(ret) <- c(class(ret), "ecv.file")
  return (ret)
}

print.ecv.regression <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
}

predict.ecv.regression <- function(object, newdata=NULL, ncores=1, ...) {
  if (is.null(newdata)) return (object@pred)
  if (object$method=="default") {
    newpred.baselearner.batch <- predict(object$est$baselearner.batch, newdata, ncores=ncores, ...)
    newpred <- predict(object$est$integrator, Xnew=newpred.baselearner.batch, config.list=object$est$baselearner.batch@config.list, ...)
  } else {
    stop("invalid CV integration method")
  }
  return (newpred)
}

summary.ecv.regression <- function(object, ...) {
  #summary.baselearner <- summary(object$est$baselearner.cv.batch, errfun=object$integrator.config@errfun) # not implemented yet in EnsembleBase
  summary.baselearner <- NULL
  summary.integrator <- summary(object$est$integrator)
  ret <- list(baselearner=summary.baselearner, integrator=summary.integrator)
  class(ret) <- "summary.ecv.regression"
  return (ret)
}

print.summary.ecv.regression <- function(x, ...) {
  #cat("### baselearner summary ###\n")
  #print(x$baselearner)
  cat("### integrator summary ###\n")
  print(x$integrator)
}

plot.ecv.regression <- function(x, ...) {
  errfun <- x$integrator.config@errfun 
  error <- errfun(x$pred, x$y)
  #plot(x$est$baselearner.batch, errfun=x$integrator.config@errfun)
  plot(x$est$baselearner.cv.batch, errfun=x$integrator.config@errfun)
  abline(h=error, lty=2)
  #plot(x$est$baselearner.cv.batch, errfun=x$integrator.config@errfun)
  #abline(h=est$est$integrator@est$error.opt, lty=2)
  #cat("full error:", error, "\n")
  #cat("cv error:", est$est$integrator@est$error.opt, "\n")
}

## determine if save and load methods are generic enough to be applicable to ALL integrators
## if yes, we should move these functions to EnsembleBase to make them available to other derivative packages

ecv.save <- function(obj, file) {
  if (!identical(class(obj),c("ecv.regression","ecv.file"))) stop("invalid object class (must be ecv.regression & ecv.file)") # must be edited to incorporate non-regression models
  if (missing(file)) stop("must provide file argument")
  tmpfiles <- obj$est$baselearner.cv.batch@tmpfiles
  tmpfiles.full <- obj$est$baselearner.batch@tmpfiles
  tmpfile.new <- tempfile()
  save(obj, file=tmpfile.new, compress=F)
  all.files <- c(tmpfile.new, tmpfiles, tmpfiles.full)
  all.files.basename <- basename(all.files)
  
  tmpdir <- paste0("./.", basename(tempfile("dir")),"/")
  dir.create(tmpdir)
  all.files.new <- paste0(tmpdir, all.files.basename)
  file.copy(all.files, all.files.new)
  meta <- list(filename.mainobj=all.files.basename[1], filenames.batchobj=all.files.basename[1+1:length(tmpfiles)], filenames.batchobj.full=1+length(tmpfiles)+1:length(tmpfiles.full))
  save(meta, file=paste0(tmpdir, "meta"), compress=FALSE)

  tar(file, files=tmpdir, compression="gzip")
  unlink(tmpdir, recursive=TRUE)
}

ecv.load <- function(file) {
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
  if (!identical(class(obj),c("ecv.regression","ecv.file"))) stop("invalid object class (must be ecv.regression & ecv.file)") # extend later to allow for non-regression models
  
  basenames.ordered <- basename(obj$est$baselearner.cv.batch@tmpfiles)
  basenames.ordered.full <- basename(obj$est$baselearner.batch@tmpfiles)
  #if (!identical(sort(basenames.ordered),sort(basenames[-c(metafile.index,mainfile.index)]))) stop("basenames mismatch")
  filepaths.ordered <- paste(extdir, basenames.ordered, sep="/")
  filepaths.ordered.full <- paste(extdir, basenames.ordered.full, sep="/")
  
  # copy batch files to new tempfiles in tempdir
  tmpfiles.new <- tempfile(rep("file", length(filepaths.ordered)))
  # replaced file.rename with file.copy and unlink to handle cross-device cases (where . and R tmp directories are on different devices)
  file.copy(from=filepaths.ordered, to=tmpfiles.new)
  unlink(filepaths.ordered, recursive=TRUE)
  
  tmpfiles.new.full <- tempfile(rep("file", length(filepaths.ordered.full)))
  # replaced file.rename with file.copy and unlink to handle cross-device cases (where . and R tmp directories are on different devices)
  file.copy(from=filepaths.ordered.full, to=tmpfiles.new.full)
  unlink(filepaths.ordered.full, recursive=TRUE)

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
  
  obj$est$baselearner.batch@tmpfiles <- tmpfiles.new.full
  n.config <- length(obj$est$baselearner.batch@config.list)
  for (i in 1:n.config) {
    obj$est$baselearner.batch@fitobj.list[[i]]@est <- tmpfiles.new.full[i]
  }
  
  return (obj)
}



