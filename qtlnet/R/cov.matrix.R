cov.formula <- function(addcov=NULL, intcov=NULL, nQ, cross.type = "f2")
{
  if(nQ > 0){
    Qnames <- paste("add",1:nQ,sep="")
    if(cross.type == "f2")
      Qnames <- as.vector(rbind(Qnames, paste("dom",1:nQ,sep="")))
  }
  if(!is.null(addcov) & is.null(intcov)){ 
    if(nQ > 0){
      myform <- as.formula(paste("y ~ ", paste(c(addcov,Qnames), collapse = "+")))
    }
    else{
      myform <- as.formula(paste("y ~ ", paste(addcov, collapse = "+")))
    }
  }
  if(!is.null(intcov)){
    le <- length(intcov)	
    intaddcov <- unique(c(intcov,addcov))
    if(nQ > 0){
      intQnames <- c()
      for(i in 1:le){
        intQnames <- c(intQnames,paste(intaddcov[i], Qnames, sep=":"))
      }
      myform <- as.formula(paste("y ~ ", paste(c(intaddcov,Qnames,intQnames), collapse = "+")))
    }
    else{
      myform <- as.formula(paste("y ~ ", paste(intaddcov, collapse = "+")))
    }
  }
  if(is.null(addcov) & is.null(intcov)){
    if(nQ > 0){
      myform <- as.formula(paste("y ~ ", paste(Qnames, collapse = "+")))
    }
    else{
      myform <- as.formula("y ~ 1")
    }
  }	
  myform
}
#######################################################################
create.cov.matrix <- function(cross, cov.names)
{
  if(is.numeric(cov.names))
    cov.names <- names(cross$pheno)[cov.names]
  
  if(length(cov.names)){
    myformula <- formula(paste("~", paste(cov.names, collapse = "+")))
    model.matrix(myformula, cross$pheno)[,-1, drop = FALSE]
  }
  else
    NULL
}
#######################################################################
pull.pheno.null <- function(cross, cols) {
  if(is.null(cols))
    NULL
  else {
    if(length(cols))
      cross$pheno[, cols, drop = FALSE]
    else
      NULL
  }
}
