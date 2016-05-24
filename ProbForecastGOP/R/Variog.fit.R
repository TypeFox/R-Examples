Variog.fit <-
function(emp.variog,variog.model="exponential",max.dist.fit=NULL,init.val=NULL,fix.nugget=FALSE){
# INPUT CHECK
### Here there should be the check on whether emp.variog is an object output of the function emp.variog
init.var <- emp.variog$mar.var
emp.variog.mod <- list(bin.midpoints=emp.variog$bin.midpoints,number.pairs=emp.variog$number.pairs,empir.variog=emp.variog$empir.variog)
# Here we do the input check on the rest of the input
# default
if(missing(max.dist.fit))
  max.dist.fit <- NULL
if(missing(fix.nugget))
  fix.nugget <- FALSE
## check on the variog.model
l.variog <- length(variog.model)
case.var <- 0
if(l.variog==0){
   variog.model <- "exponential"
   case.var <- 1
}
if(l.variog==1){
  if(is.character(variog.model)=="TRUE"){
    if(variog.model=="exponential"){
        case.var <- 1}
    if(variog.model=="spherical"){
        case.var <- 1}
    if(variog.model=="whittlematern" | variog.model=="matern"){
        case.var <- 1}
    if(variog.model=="gencauchy"){
        case.var <- 1}
    if(variog.model=="gauss"){
        case.var <- 1}
  }
}
if(case.var==0){
   stop("Incorrect variogram model specification")
}
# check on the max.dist.fit and the initial values
l.max.dist.fit <- length(max.dist.fit)
if(l.max.dist.fit > 1){
  stop("max.dist.fit should be numeric field: not a vector")
}
if(l.max.dist.fit==1 & is.numeric(max.dist.fit)==FALSE){
  stop("max.dist.fit should be a numeric field")
}
if(l.max.dist.fit==1 & is.numeric(max.dist.fit)==TRUE){
  if(max.dist.fit < 0){
    stop("max.dist.fit should be positive number")
  }
  if(max.dist.fit > max(emp.variog.mod$bin.midpoints)){
    stop("max.dist.fit should be less or equal than max.dist")
  }
}
 
## Input check on the initial values for the parameters
l.init.val <- length(init.val)
if(l.init.val > 0 & l.init.val <3 ){
  stop("init.val should be equal to NULL or to a vector of at least length 3")
}
if(l.init.val ==3 & (sum(is.numeric(init.val)==rep("TRUE",3)) < l.init.val)){
  stop("The initial values should be numeric entries")  
}
if(l.init.val > 3 & (sum(is.numeric(init.val)==rep("TRUE",3))== l.init.val)){
  if(l.init.val==4 & variog.model!="matern"){
    stop("Incorrect number of initial values")
    }
  if(l.init.val==5 & variog.model!="gencauchy"){
    stop("Incorrect number of initial values")
    }
  if(l.init.val > 5){
    stop("Incorrect number of initial values")
    }
  if(init.val[1] <0){
    stop("The nugget effect cannot be negative")
  }
  if(sum(init.val[2:l.init.val]>0)<(l.init.val-1)){
    stop("Initial values for all the parameters, but the nugget effect, should be positive numbers")
    }
}
## Input check for the fix.nugget field
l.nug <- length(fix.nugget)
if(l.nug==0){
  fix.nugget <- "FALSE"}
l.nug <- length(fix.nugget)
if(l.nug==1){
  if(fix.nugget!=TRUE & fix.nugget!=FALSE){
   stop("fix.nugget should be either equal to TRUE or FALSE")
   }
}
if(l.nug==2){
  if(fix.nugget[1]!=TRUE){
    stop("Invalid input for the fix.nugget field")
  }
  if(is.numeric(fix.nugget[2])==FALSE){
    stop("The second entry of the fix.nugget field should be a numeric field")
  }
  if(fix.nugget[2] <0){
    stop("The second entry of the fix.nugget field should be a non-negative number")
  }
  if(l.init.val >0){
    if(fix.nugget[2]!=init.val[1]){
       fix.nugget[2] <- init.val[1]
    }
  }   
}
if(l.nug >2){
  stop("fix.nugget is either a character field or a 2x1 field")
}
param.est <- round(model.fit(variog.model,emp.variog.mod,max.dist.fit,init.val,init.var,fix.nugget),3)
if(variog.model=="matern"){variog.model <- "whittlematern"}
ifelse((variog.model=="whittlematern" | variog.model=="gencauchy"),
   output <- 
list(model=variog.model,nugget=param.est[1],variance=param.est[2],range=param.est[3],additional.par=param.est[-seq(1:3)]),
   output <- 
list(model=variog.model,nugget=param.est[1],variance=param.est[2],range=param.est[3]))
   return(output)
}

