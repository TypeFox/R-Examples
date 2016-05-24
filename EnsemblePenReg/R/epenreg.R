epenreg.baselearner.control <- function(baselearners=c("nnet","rf","svm","gbm","knn")
  , baselearner.configs=make.configs(baselearners, type="regression"), npart=1, nfold=5) {
    return (list(configs=baselearner.configs, npart=npart, nfold=nfold))
}

epenreg.integrator.control <- function(errfun=rmse.error, alpha=1.0, n=100, nfold=5, method=c("default")) {
  return (list(errfun=rmse.error, alpha=alpha, n=n, nfold=nfold, method=method))
}

epenreg.set.filemethod <- function(formula, data, instance.list, type="regression") FALSE # TODO: to be implemented

epenreg <- function(formula, data
  , baselearner.control=epenreg.baselearner.control()
  , integrator.control=epenreg.integrator.control()
  , ncores=1, filemethod=FALSE, print.level=1) {
  if (integrator.control$method!="default") stop("invalid PenReg integration method")
  ncores.max <- try(detectCores(),silent=T)
  mycall <- match.call()
  if (!inherits(ncores.max,"try-error")) ncores <- min(ncores,ncores.max)
  if (print.level>=1 && ncores>1) cat("running in parallel mode, using", ncores, "cores\n")

  # training a batch of baselearners
  partitions.bl <- generate.partitions(baselearner.control$npart, nrow(data), baselearner.control$nfold)
  my.instance.list <- make.instances(baselearner.control$configs, partitions.bl)

  # TODO: determining of base learner estimation objects must be saved to disk or not
  if (missing(filemethod)) filemethod <- epenreg.set.filemethod(formula, data, my.instance.list)
  
  if (print.level>=1) cat("CV training of base learners...\n")
  est.baselearner.cv.batch <- Regression.CV.Batch.Fit(my.instance.list, formula, data, ncores=ncores, filemethod=filemethod, print.level=print.level)
  if (print.level>=1) cat("finished CV training of base learners\n")
  Xcv <- est.baselearner.cv.batch@pred
  y <- data[,all.vars(formula)[1]] # TODO: more robust way of extracting y (here and inside base learner functions)
  
  partition.int <- generate.partition(nrow(data), integrator.control$nfold)
  my.integrator.config <- Regression.Integrator.PenReg.SelMin.Config(errfun=integrator.control$errfun, partition=partition.int, n=integrator.control$n, alpha=integrator.control$alpha)
  est.integrator <- Regression.Integrator.Fit(my.integrator.config, X=Xcv, y=y, print.level=print.level)
  pred <- est.integrator@pred
  
  ret <- list(call=mycall, formula=formula, instance.list=my.instance.list, integrator.config=my.integrator.config, method=integrator.control$method
    , est=list(baselearner.cv.batch=est.baselearner.cv.batch, integrator=est.integrator)
    , y=y, pred=pred, filemethod=filemethod)
  class(ret) <- "epenreg"
  if (filemethod) class(ret) <- c(class(ret), "epenreg.file")
  return (ret)
}

predict.epenreg <- function(object, newdata=NULL, ncores=1, ...) {
  if (is.null(newdata)) return (object$pred)
  if (object$method=="default") {
    newpred.baselearner.cv.batch <- predict(object$est$baselearner.cv.batch, newdata, ncores=ncores, ...)
    newpred <- predict(object$est$integrator, Xnew=newpred.baselearner.cv.batch, ...)
  } else {
    stop("invalid PenReg integration method")
  }
  return (as.numeric(newpred))
}

plot.epenreg <- function(x, ...) {
  errfun <- x$integrator.config@errfun 
  error <- errfun(x$pred, x$y)
  oldpar <- par(mfrow=c(1,2))
  plot(x$est$baselearner.cv.batch, errfun=errfun)
  abline(h=error, lty=2)
  pcr.errors <- x$est$integrator@est$select@est$error
  alpha <- x$est$integrator@est$penreg@sweep.list[[1]]@config@alpha
  lambda <- x$est$integrator@est$penreg@sweep.list[[1]]@config@lambda
  plot(lambda, pcr.errors, type="l", xlab="Shrinkage Parameter (lambda)", ylab="CV Error", main=paste0("Integrator Performance (alpha=", alpha, ")"), log="x")
  par(oldpar)
}

epenreg.save <- function(obj, file) {
  if (!identical(class(obj),c("epenreg","epenreg.file"))) stop("invalid object class (must be epenreg & epenreg.file)")
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

epenreg.load <- function(file) {
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
  if (!identical(class(obj),c("epenreg","epenreg.file"))) stop("invalid object class (must be epenreg & epenreg.file)")
  
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


print.epenreg <- function(x, ...) {
  cat("Call:\n")
  print(x$call)
}

summary.epenreg <- function(object, ...) {
  n.instance <- length(object$instance.list@instances)
  nlambda <- object$est$integrator@est$penreg@sweep.list[[1]]@config@n
  index.min <- object$est$integrator@est$select@est$index.min
  lambda.opt <- object$est$integrator@est$penreg@sweep.list[[1]]@config@lambda[index.min]
  error.min <- object$est$integrator@est$select@est$error.min
  ret <- list(n.instance=n.instance, nlambda=nlambda, index.min=index.min, lambda.opt=lambda.opt, error.min=error.min)
  
  class(ret) <- "summary.epenreg"
  return (ret)
}

print.summary.epenreg <- function(x, ...) {
  cat("number of base learner instances:", x$n.instance, "\n")
  cat("number of lambda values considered:", x$nlambda, "\n")
  cat("optimal lambda:", x$lambda.opt, "\n")
  cat("\toptimal index:", x$index.min, "\n")
  cat("minimum error:", x$error.min, "\n")
}




