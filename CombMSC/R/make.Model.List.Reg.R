`make.Model.List.Reg` <-
function(fram, max.Size=min(8, dim(fram)[2]), no.Intercepts = FALSE, GLM = FALSE){
  if(!is.data.frame(fram)) stop ("fram must be a data.frame")
  if(max.Size > dim(fram)[2]) stop ("max.Size cannot be larger than the number of variables in fram")
  num.Vars <- dim(fram)[2]
  var.Names <- colnames(fram)
  subset.List <- list()
  for(j in 1:max.Size){
    # if condition necessary since if j==num.Vars, there is only ONE subset of size j
    # and hence we can't use apply on it
    ifelse(j==num.Vars, subset.List[[j]] <- paste(subsets(num.Vars,j,var.Names), collapse="+"),
    subset.List[[j]] <- apply(subsets(num.Vars, j, var.Names),1,function(x)paste(x,collapse="+")))}

  all.Subsets <- unlist(subset.List)
  all.Subsets <- c("1", all.Subsets)
  if(no.Intercepts) all.Subsets <- c(all.Subsets, paste(all.Subsets, " + 0", sep = ''))
  all.Formulas <- sapply(all.Subsets, function(x)paste("y ~ ", x, sep=''))
  names(all.Formulas) <- NULL
  all.Formulas <- lapply(all.Formulas, function(x) x <- as.formula(x))
  for(i in 1:length(all.Formulas)){
    class(all.Formulas[[i]]) <- c("lmFormula", "formula")
    if(GLM) class(all.Formulas[[i]]) <- c("glmFormula", "lmFormula", "formula")
    }
  all.Formulas}

