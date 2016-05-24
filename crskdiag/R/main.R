Crsk <- function(t,ic) {cbind(t,ic)}

diag_crr <- function(formula, data, test=c("lin","prop"),
                     Nit=20,n.sim=1000,n.plot=10,seed=NULL,minor_included=1) {
  
  # Pre-processing for formula
  Call <- match.call(expand.dots = FALSE)
  indx <- match(c("formula","data"),names(Call),nomatch=0)
  if(indx[1]==0) stop("A formula argument is required.")
  temp <- Call[c(1,indx)]
  temp[[1]] <- as.name("model.frame") 
  
  special <- c("strata","cluster","tt")
  if(missing(data)) temp$formula <- terms(formula,special)
  else temp$formula <- terms(formula,special,data=data)
  
  if(is.R()) m <- eval(temp,parent.frame())  #obs with missing is removed here
  else m <- eval(temp,sys.parent())

  if(nrow(m)==0) stop("No (non-missing) observations")
  
  Y <- model.extract(m,"response")
  t <- Y[,1]
  ic <- Y[,2]
  
  # Factor  - class variables
  fac.indx <- NULL
  for(i in 2:dim(m)[2]) if(is.factor(m[,i])) fac.indx <- c(fac.indx,i)
  
  fac.len <- NULL
  if(length(fac.indx)>0) {
    fac.label <- NULL
    for(i in 1:length(fac.indx)) {
      tmp.string <- gsub("factor\\(","",colnames(m)[fac.indx[i]])
      tmp.string <- gsub("\\)","",tmp.string)
      tmp.string <- paste(tmp.string," : ",sep="")
      fac.label <- c(fac.label,tmp.string)
    }
    
    fac.mat <- model.matrix(~factor(m[,fac.indx[1]]))[,-1]
    fac.len <- ncol(as.matrix(fac.mat))
    if(fac.len == 1 && is.vector(fac.mat)) {
      fac.mat <- matrix(fac.mat,ncol=1)
      i <- 1
      colnames(fac.mat) <- paste("factor(m[, fac.indx[",i,"]])",
                                 strsplit(fac.label[1]," ")[[1]][1],sep="")
    }
    if(length(fac.indx)>1) {
      for(i in 2:length(fac.indx)) {
        tmp <- model.matrix(~factor(m[,fac.indx[i]]))[,-1]
        fac.mat <- cbind(fac.mat,tmp)
        fac.len <- c(fac.len,ncol(as.matrix(tmp)))
        if(colnames(fac.mat)[dim(fac.mat)[2]]=="tmp") {
          tmp.string <- paste("factor(m[, fac.indx[",i,"]])",
                              strsplit(fac.label[i]," ")[[1]][1],sep="")
          colnames(fac.mat)[dim(fac.mat)[2]] <- tmp.string
        }
      }
    }
    
    exp.label <- NULL
    for(i in 1:length(fac.label)) {
      exp.label <- c(exp.label,rep(fac.label[i],fac.len[i]))
    }
    
    fac.names <- colnames(fac.mat)
    for(i in 1:length(fac.names)) {
      name.split <- strsplit(fac.names[i],")")
      name.common <- paste(name.split[[1]][1],")",sep="")
      name.common <- gsub("\\(","\\\\\\(",name.common)
      name.common <- gsub("\\)","\\\\\\)",name.common)
      name.common <- gsub("\\[","\\\\\\[",name.common)
      name.common <- gsub("\\]","\\\\\\]",name.common)
      
      colnames(fac.mat)[i] <- gsub(name.common,exp.label[i],fac.names[i])
    }
    
    m <- cbind(m,fac.mat)
    m <- m[,-fac.indx]
  }
  
  z <- as.matrix(m[,2:dim(m)[2]],nc=dim(m[,2:dim(m)[2]])[2])
  colnames(z) <- colnames(m)[-1]
  
  n.total <- length(t)
  n.missing <- 0
  var.name <- colnames(z) 
  miss.index <- NULL
  
#   # Missing values removal
#   for(i in 1:n.total){
#     if(is.na(t[i]) || is.na(ic[i]) || any(is.na(z[i,]))) {
#         n.missing <- n.missing+1
#         miss.index <- c(miss.index,i)
#     }     
#   }
#   if(is.null(miss.index)==FALSE){
#     t <- t[-miss.index]
#     ic <- ic[-miss.index]
#     z <- z[-miss.index,]
#   }
#   n <- n.total-n.missing
  
#   if(fail.code > 1) {
#     ic <- ifelse(ic==1, 999, ic)  
#     ic <- ifelse(ic==fail.code, 1, ic)
#     ic <- ifelse(ic==999, fail.code, ic)
#   }
  
  # Fitting the model
  if(test=="lin") {
    res <- diag_lin(t,ic,z,n.total,Nit,n.sim,n.plot,seed,minor_included)
    return(res)
  }
  if(test=="prop") {
    res <- diag_prop(t,ic,z,n.total,Nit,n.sim,n.plot,seed,minor_included)
    return(res)
  }
}