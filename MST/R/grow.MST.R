grow.MST <-
function(dat, test=NULL, method=c("marginal", "gamma.frailty", "exp.frailty"),
                     col.time, col.status, col.id, col.split.var, col.ctg=NULL,
                     minsplit=20, min.nevents=3, max.depth=10, mtry=length(col.split.var),
                     cont.split=c("distinct","percentiles"), delta=0.05, nCutPoints=50, details=FALSE){
  method<-match.arg(method,c("marginal", "gamma.frailty", "exp.frailty"))
  cont.split<-match.arg(cont.split,c("distinct", "percentiles"))

  out <- list.nd <- list.test <- temp.list <- temp.test <- temp.name <- NULL
  list.nd <- list(dat); 
  if (!is.null(test)) list.test <- list(test)
  name <- 1
  while (length(list.nd)!=0) {
    for (i in 1:length(list.nd)){
      if (!is.null(dim(list.nd[[i]])) && nrow(list.nd[[i]]) > 1){ 
        test0 <- NULL
        if (!is.null(test)) test0 <- list.test[[i]]
        split <- partition.MST(dat=list.nd[[i]], test=test0, name=name[i], method=method,
                               col.time=col.time, col.status=col.status, col.id=col.id, col.split.var=col.split.var, col.ctg=NULL, 
                               minsplit=minsplit, min.nevents=min.nevents, max.depth=max.depth, mtry=mtry,
                               cont.split=cont.split, delta=delta, nCutPoints=nCutPoints, details=details)
        # print(split$info)
        out <- rbind(out, split$info)
        if (!is.null(split$left) && is.null(test)) {
          temp.list <- temp.test <- c(temp.list, list(split$left, split$right))
          temp.name <- c(temp.name, split$name.l, split$name.r)
        } else if (!is.null(split$left) && !is.null(test) && !is.null(split$left.test)) {
          temp.list <- c(temp.list, list(split$left, split$right))
          temp.name <- c(temp.name, split$name.l, split$name.r)
          temp.test <- c(temp.test, list(split$left.test, split$right.test))
        }		
      }
    }
    list.nd <- temp.list; list.test <- temp.test; name <- temp.name
    temp.list <- temp.test <- temp.name <- NULL
  }
  out$node <- as.character(out$node)
  out <- out[order(out$node), ]
  out
}
