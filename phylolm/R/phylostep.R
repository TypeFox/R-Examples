################################################
### Stepwise selection for phylolm using AIC
################################################

phylostep <- function(formula, starting.formula = NULL, data=list(), phy, 
                    model=c("BM","OUrandomRoot","OUfixedRoot","lambda","kappa","delta","EB","trend"),
                    direction = c("both", "backward", "forward"), trace = 2,
                    lower.bound=NULL, upper.bound=NULL, starting.value=NULL,
                    k=2, ...) 
{
  ## initialize  
  model = match.arg(model)	
  direction = match.arg(direction)  
  fit.full = phylolm(formula, data, phy, model, lower.bound, upper.bound, starting.value, ...)
  response = fit.full$formula[[2]] # name of the response
  covariates = names(fit.full$coefficients) # name of the covariates
  if (length(data) == 0) {
    data = as.data.frame(fit.full$X)
  } # check if there is no data structure
  
  if (is.null(covariates)) stop("Covariates has no name.")
  p = length(covariates)
  if (p==1) stop("Your full model only has intercept. Model selection is not needed.")
  plm.full = rep(1,p)

  ## create a formula from current model
  create.formula <- function(plm) {
    on = which(plm==1)
    str = paste(response," ~ 1",sep="")
    for (i in 1:length(on))
      if (i>1) str = paste(str," + ",covariates[on[i]],sep="")
    return(str)
  }
  
  ## fit phylolm
  fit <- function(plm) {
    return(phylolm(create.formula(plm), data, phy, model, 
                   lower.bound, upper.bound, starting.value, ...))
  }
  
  ## plm.current is a binary vector of length p
  ## where 0 at i-th position means excluding i-th covariate
  
  if (direction == "forward") {
    plm.current = c(1,rep(0,p-1)) # only intercept
    fit.current = fit(plm.current)
  } else {
    plm.current = plm.full # all covariates
    fit.current = fit.full
  }
  
  if (!is.null(starting.formula)) {
    fit.current = phylolm(starting.formula, data, phy, model, lower.bound, upper.bound, starting.value, ...)
    covariates.current = names(fit.current$coefficients)
    plm.current = rep(0,p)
    position = match(covariates.current,covariates)
    if (any(is.na(position))) stop("The starting model is not a submodel of the full model.")
    plm.current[position] = 1
  }

  if (trace>0) {
    cat("----------\n")
    cat(paste("Starting model: ",create.formula(plm.current),"\n",sep=""))
    cat(paste("Direction: ",direction,"\n",sep=""))
    cat(paste("AIC(k=",k,"): ",AIC(fit.current,k),"\n",sep=""))
#     cat("----------\n")
  }
    
  flag = 0 # flag of termination
  count = 1
  while (flag == 0) {
    flag = 1
    plm.best = plm.current
    fit.best = fit.current
    
    ### to preserve hierarchy priciple
    terms.add = add.scope(formula(create.formula(plm.current)), 
                                  formula(create.formula(plm.full)))
    terms.drop = drop.scope(formula(create.formula(plm.current)))
    
    for (i in 2:p) {
      
      plm.propose = plm.current
      do.update = FALSE
      if ((plm.current[i]==1)&&(direction %in% c("both", "backward"))
          &&(covariates[i] %in% terms.drop)) {
        plm.propose[i] = 0  # remove i-th covariate
        do.update = TRUE
      }
      if ((plm.current[i]==0)&&(direction %in% c("both", "forward"))
          &&(covariates[i] %in% terms.add)) {
        plm.propose[i] = 1 # add i-th covariate
        do.update = TRUE
      }
        
      ## check if proposed model is better
      if (do.update){
        fit.propose = fit(plm.propose)
        if (trace>1) {
          cat(paste0("\tProposed: ",create.formula(plm.propose),"\n"))
          cat(paste0("\tAIC(k=",k,"): ",AIC(fit.propose,k),"\n"))
        }
        if (AIC(fit.propose,k) < AIC(fit.best,k)) {
          plm.best = plm.propose
          fit.best = fit.propose
          flag = 0
        }
      }
    }
    
    ## Set current model as the best model
    plm.current = plm.best
    fit.current = fit.best
    
    if (trace>0) {
      cat("----------\nStep ",count,"\n",sep="")
      cat(paste("Current model: ",create.formula(plm.current),"\n",sep=""))
      cat(paste("AIC(k=",k,"): ",AIC(fit.current,k),"\n",sep=""))
    }
    count = count + 1
    }
  if (trace>0) cat("---END\n")
  return(fit.current)
}  
