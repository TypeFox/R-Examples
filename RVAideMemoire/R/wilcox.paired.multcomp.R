wilcox.paired.multcomp <-
function(formula,data,p.method="fdr") {
  if (missing(formula)) {stop("formula missing")}
  if ((length(formula)!=3) || (length(formula[[3]])!=3) || (formula[[3]][[1]]!=as.name("|")) ||
    (length(formula[[3]][[2]])!=1) || (length(formula[[3]][[3]])!=1)) {stop("incorrect specification for formula")}
  formula[[3]][[1]] <- as.name("+")
  m <- match.call()
  m$formula <- formula
  if (is.matrix(eval(m$data,parent.frame()))) {m$data <- as.data.frame(m$data)}
  m[[1]] <- as.name("model.frame")
  m$p.method <- NULL
  mf <- eval(m,parent.frame())
  mf <- droplevels(mf[complete.cases(mf),])
  dname <- paste(names(mf)[1]," by ",names(mf)[2],", block = ",names(mf)[3],sep="")
  resp <- mf[,1]
  fact <- mf[,2]
  block <- mf[,3]
  tab <- data.frame(fact,block,resp)
  tab <- tab[order(tab$fact),]
  method <- "Wilcoxon signed rank test"
  fun.p <- function(i,j) {
    test <- suppressWarnings(wilcox.test(tab$resp[as.integer(tab$fact)==i],tab$resp[as.integer(tab$fact)==j],paired=TRUE))
    test$p.value
  }
  comp <- pairwise.table(fun.p,levels(fact),p.adjust.method=p.method)
  result <- list(data.name=dname,method=method,p.adjust.method=p.method,p.value=comp)
  class(result) <- "pairwise.htest"
  return(result)
}
