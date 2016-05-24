SetStaticVariables <-
function(dataset, covariates) {
  # if-test to preserve class, avoids possible numerical covariates as factors in lm
  if(length(covariates)==1) {
    dataset[, covariates] <- na.locf(dataset[, covariates])
  } 
  if(length(covariates)>=2) {
    classes <- sapply(dataset[, covariates], class)
    u.class <- unique(classes)
    for(cla in u.class) {
      do.me <- covariates[which(classes==cla)]
      dataset[, do.me] <- na.locf(dataset[, do.me])
    }
  }
  return(dataset)
}
