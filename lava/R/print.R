###{{{ print.lvm

##' @export
`print.lvm` <-
function(x, ...) {
  res <- NULL
  myhooks <- gethook("print.hooks")
  for (f in myhooks) {
      res <- do.call(f, list(x=x,...))
  }
  if (is.null(res)) {
    k <- length(vars(x))
    L <- rep(FALSE,k); names(L) <- vars(x); L[latent(x)] <- TRUE
    cat("Latent Variable Model\n") ##;" \n\twith: ", k, " variables.\n", sep="");
    if (k==0)
      return()
    ff <- formula(x,char=TRUE,all=TRUE)
    R <- Rx <- c()
    exo <- exogenous(x)
    for (f in ff) {
      oneline <- as.character(f);
      y <- strsplit(f,"~")[[1]][1]
      y <- trim(y)
      {
        col1 <- as.character(oneline)          
        D <- attributes(distribution(x)[[y]])$family
        col2 <- x$attributes$type[y]
        if (is.null(col2) || is.na(col2)) col2 <- "gaussian"
        if (!is.null(D$family)) {
            col2 <- paste0(D$family)
        }
        if (!is.null(D$link)) col2 <- paste0(col2,"(",D$link,")")
        if (!is.null(D$par)) col2 <- paste0(col2,"(",paste(D$par,collapse=","),")")
        if (is.list(distribution(x)[[y]]) && is.vector(distribution(x)[[y]][[1]])) col2 <- "fixed"
        if (L[y]) col2 <- paste0(col2,", latent")
        if (y%in%exo) {
            Rx <- rbind(Rx,c(col1,col2))
        } else {
            R <- rbind(R,c(col1,col2))
        }
      }
    }
    if (length(R)>0) {
        rownames(R) <- paste(" ",R[,1]," "); colnames(R) <- rep("",ncol(R))
        print(R[,2,drop=FALSE],quote=FALSE,...)
    }
    if (length(Rx)>0) {
        cat("\nExogenous variables:")
        rownames(Rx) <- paste(" ",Rx[,1]," "); colnames(Rx) <- rep("",ncol(Rx))
        print(Rx[,2,drop=FALSE],quote=FALSE,...)
    }
    

##      oneline <- as.character(f);
##      cat(as.character(oneline),"\n")

  }
  cat("\n")
  invisible(x)
}

###}}} print.lvm

###{{{ print.lvmfit

##' @export
`print.lvmfit` <-
function(x,level=2,labels=FALSE,...) {
    print(CoefMat(x,labels=labels,level=level,...),quote=FALSE,right=TRUE)
    minSV <- attr(vcov(x),"minSV")
    if (!is.null(minSV) && minSV<1e-12) {
        warning("Small singular value: ", format(minSV))
    }
    pseudo <- attr(vcov(x),"pseudo")
    if (!is.null(pseudo) && pseudo) warning("Singular covariance matrix. Pseudo-inverse used.")
    invisible(x)
}

###}}} print.lvmfit

###{{{ print.lvmfit.randomslope

##' @export
print.lvmfit.randomslope <- function(x,labels=FALSE,level=2,...) {
  print(CoefMat(x,labels=labels,level=level,...),quote=FALSE,right=TRUE)
  invisible(x)
}

###}}}

###{{{ print.multigroupfit

##' @export
print.multigroupfit <- function(x,groups=NULL,...)  {
  if (is.null(groups)) {
    if (x$model$missing) {
      modelclass <- attributes(x$model0)$modelclass
      nmis <- attributes(x$model0)$nmis
      orggroup <- unique(modelclass)
      groupn <- unlist(lapply(orggroup,function(i) sum(modelclass==i)))
      cumsumgroup <- cumsum(c(0,groupn))
      groups <- unlist(lapply(orggroup,function(i)
                              which.min(nmis[which(modelclass==i)])+cumsumgroup[i])) ##  groups with max. number of variables
      for (i in seq_len(length(groups))) {
        if (nmis[groups[i]]>0) warning("No complete cases in group ",i,". Showing results of group with max number of variables. All coefficients can be extracted with 'coef'. All missing pattern groups belonging to this sub-model can be extracted by calling: coef(..., groups=c(",paste(which(modelclass==i),collapse=","),"))")
      }
      if (!is.null(x$model$mnameses))
        x$model$names <- x$model$mnames
    } else {
      groups <- seq_len(length(x$model$lvm))
    }
  }
  res <- coef(x,level=2,groups=groups,...)
  counter <- 0
  dots <- list(...)
  dots$groups <- groups
  level <- if (is.null(dots$level)) {
    dots$level <- 2
##    dots$level <- ifelse("lvmfit.randomslope"%in%class(x),2,9)
  }
  myargs <- c(list(x=x), dots)
  myargs$groups <- groups
  CC <- do.call("CoefMat.multigroupfit",myargs)
  for (cc in res) {
    counter <- counter+1
    cat(rep("_",52),"\n",sep="")
    cat("Group ", counter, sep="")
    myname <- x$model$names[counter]
    if (!is.null(myname) && !is.na(myname))
      cat(": ",myname,sep="")
    if (!x$model$missing) cat(" (n=",nrow(Model(x)$data[[groups[counter]]]), ")", sep="")
    cat("\n")
    print(CC[[counter]],quote=FALSE,right=TRUE)
  }
  cat("\n")
  invisible(x)
}

###}}} print.multigroupfit

###{{{ print.multigroup

##' @export
print.multigroup <- function(x,...) {
  cat("\n")
  cat("Number of groups:", x$ngroup, "\n")
  cat("Number of free parameters (not counting mean-parameters):", x$npar,"\n")
##  cat("Parameter-vector:", unlist(x$parlist), "\n\n")
  cat("Number of free mean parameters:", length(grep("m",x$mean)),"\n")
##  cat("Mean-vector:", x$mean, "\n\n")
  invisible(x)
}

###}}} print.multigroup

###{{{ printmany

printmany <- function(A,B,nspace=1,name1=NULL,name2=NULL,digits=3,rownames=NULL,emptystr=" ",bothrows=TRUE,print=TRUE,...) {
  cA <- colnames(A); cB <- colnames(B)
  A <- format(A, digits=digits)
  B <- format(B, digits=digits)
  nA <- nrow(A); nB <- nrow(B)
  if (nrow(A)<nrow(B)) {
    rA <- rownames(A)
    A <- rbind(A, matrix("", nrow=nB-nA, ncol=ncol(A)))
  }
  if (nrow(B)<nrow(A)) {
    rB <- rownames(B)
    B <- rbind(B, rep("", nrow=nA-nB, ncol=ncol(B)))
  }
  if (!is.null(rownames) & length(rownames)==nrow(A))
    rownames(A) <- rownames(B) <- rownames
  res <- cbind(A, matrix("", nrow=nrow(A), ncol=nspace)); colnames(res) <- c(colnames(A), rep(emptystr,nspace))
  if (!is.null(name1)) {
    oldname <- colnames(res)
    res <- cbind(rep("",nrow(res)), rownames(res), res); colnames(res) <- c(name1,"",oldname)
    rownames(res) <- rep("",nrow(res))
  }
  if (!is.null(name2)) {
    oldname <- colnames(res)
    res <- cbind(res,rep("",nrow(res))); colnames(res) <- c(oldname,name2)
  }
  if (!identical(rownames(A),rownames(B)) | bothrows)
    res <- cbind(res, rownames(B))
  res <- cbind(res, B)
  if (print) print(res, quote=FALSE,...)
  invisible(res)
}

###}}} printmany
