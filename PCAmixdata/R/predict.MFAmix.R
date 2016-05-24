predict.MFAmix<-function (object, data, groups,name.groups,...) 
{
  mfa <- object
  if (!inherits(mfa, "MFAmix")) 
    stop("use only with \"MFAmix\" objects")
  
  n <- nrow(data)
  nbr.groups <- length(unique(groups))
  Lst.groups <- splitgroups(base = data, groups = groups, name.groups = name.groups)
  long.groups <- sapply(Lst.groups, ncol)
  typ.groups <- unlist(sapply(Lst.groups, splitmix)[3, ])
  DATA.ord <- data.frame(matrix(NA, ncol = ncol(data), nrow = nrow(data)))
  init <- 0
  for (g in 1:nbr.groups) {
    DATA.ord[, c((1 + init):(init + ncol(Lst.groups[[g]])))] <- Lst.groups[[g]]
    init <- init + ncol(Lst.groups[[g]])
  }
  colnames(DATA.ord) <- unlist(sapply(Lst.groups, colnames))
  rownames(DATA.ord) <- rownames(data)
  groups.ord <- NULL
  for (g in 1:nbr.groups) {
    groups.ord <- c(groups.ord, rep(g, ncol(Lst.groups[[g]])))
  }
  groups <- groups.ord
  data <- DATA.ord
  
  X.quanti<-splitmix(data)$X.quanti
  X.quali<-splitmix(data)$X.quali
  
  
  rec <- recod(X.quanti, X.quali,rename.level=TRUE)
  Y <- rec$Y
  n <- rec$n
  beta <- mfa$global.pca$coef
  
  ncol_tot<-length(beta[[1]])-1
  #test if the subsample of individual to predict has all the levels of the categ var
  if (ncol_tot!=ncol(Y)){
    Ymodif<-matrix(0,nrow=n,ncol=ncol_tot)
    colname.part<-colnames(Y)
    colname.tot<-names(beta[[1]][-1,])
    colnames(Ymodif)<-colname.tot
    Ymodif[,colname.part]<-Y[,colname.part]
    Y<-Ymodif  
  }
  
  if (!is.null(X.quanti)) {
    label <- rownames(X.quanti)
    n1 <- nrow(X.quanti)
    p1 <- ncol(X.quanti)
    if (p1 != mfa$global.pca$rec$p1) 
      stop("The number of numerical variables in data must be the same than in the learning set")
  }
  if (!is.null(X.quali)) {
    label <- rownames(X.quali)
    n2 <- nrow(X.quali)
    p2 <- ncol(X.quali)
    if (p2 != mfa$global.pca$rec$p2) 
      stop("The number of categorical variables in data must be the same than in the learning set")
  }
  if (!is.null(X.quanti) && !is.null(X.quali)) {
    if (n1 != n2) 
      stop("The number of objects in X.quanti and X.quali must be the same")
    if (sum(rownames(X.quali) != rownames(X.quanti)) != 0) 
      stop("The names of the objects in X.quanti and X.quali must be the same")
  }
  scores <- matrix(, n, length(beta))
  for (g in 1:length(beta)) scores[, g] <- Y %*% beta[[g]][-1] + 
    beta[[g]][1]
  colnames(scores) <- paste("dim", 1:length(beta), sep = "")
  rownames(scores) <- label
  return(scores)
}


