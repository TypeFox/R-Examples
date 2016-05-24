bootstrap.grow.prune <-
function(data, method=c("marginal", "gamma.frailty", "exp.frailty"),
                                 col.time, col.status, col.id, col.split.var, col.ctg=NULL,
                                 minsplit=20, min.nevents=3, max.depth=10, mtry=length(col.split.var),
                                 cont.split=c("distinct","percentiles"), delta=0.05, nCutPoints=50,
                                 B=30, LeBlanc=TRUE, min.boot.tree.size=1, details=FALSE){
  method<-match.arg(method,c("marginal", "gamma.frailty", "exp.frailty"))
  cont.split<-match.arg(cont.split,c("distinct", "percentiles"))

  call <- match.call(); out <- match.call(expand.dots = FALSE)
  out$boot.tree <- out$boot.prune <- out$... <- NULL
  time.start <- date()
  tree0 <- grow.MST(dat=data, test=NULL, method=method,
                    col.time=col.time, col.status=col.status, col.id=col.id, col.split.var=col.split.var, col.ctg=col.ctg,
                    minsplit=minsplit, min.nevents=min.nevents, max.depth=max.depth, mtry=mtry,
                    cont.split=cont.split, delta=delta, nCutPoints=nCutPoints, details=details)
  if (NROW(tree0)==1) {print(tree0); stop("Your initial tree is a NULL tree! There is no need for tree size selection.")}
  # print(tree0);
  prune0 <- prune.size(tree0); 
  boot.tree <- list(tree0); boot.prune <- list(prune0) 
  b <- 1
  while (b <= B) {
    # SAMPLING OBSERVATION
    samp <- sample(1:nrow(data), size=nrow(data), replace=TRUE) 
    dat <- data[samp, ]
    dat.oob <- data[-unique(samp),]
    n.oob <- nrow(dat.oob); # print(n.oob)
    if (LeBlanc) tree <- grow.MST(dat=dat, test=data, method=method,
                                  col.time=col.time, col.status=col.status, col.id=col.id, col.split.var=col.split.var, col.ctg=col.ctg,
                                  minsplit=minsplit, min.nevents=min.nevents, max.depth=max.depth, mtry=mtry,
                                  cont.split=cont.split, delta=delta, nCutPoints=nCutPoints, details=details)
    else  tree <- grow.MST(dat=dat, test=dat.oob, method=method,
                           col.time=col.time, col.status=col.status, col.id=col.id, col.split.var=col.split.var, col.ctg=col.ctg,
                           minsplit=minsplit, min.nevents=min.nevents, max.depth=max.depth, mtry=mtry,
                           cont.split=cont.split, delta=delta, nCutPoints=nCutPoints, details=details)
    if (details) print(tree)
    if (nrow(tree)> min.boot.tree.size) {
      if((b %% 5)==0){print(paste0("Current Bootstrap Sample: ", b))}
      boot.tree <- c(boot.tree, list(tree))
      prune <- prune.size(tree); # print(prune)
      boot.prune <- c(boot.prune, list(prune))
      b <- b+1
    }
  }
  time.end <- date()
  print(paste("The Start and End time for ", B, "bootstrap runs is:"))
  print(rbind(time.start, time.end))
  out$boot.tree <- boot.tree
  out$boot.prune <- boot.prune
  # THE INITIAL LARGE TREE
  out$initial.tree <- tree0
  out
}
