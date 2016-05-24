#' A harvested classification tree 
#' 
#' Basic function to apply the harvest algorithm to the training data set, computing whether we can harvest any nodes based on the classic classification tree algorithm.  
#' @param rpart.object classification result of training data from traditional classification tree(rpart function).
#' @param data original training data where 'y' stores classmembership
#' @param varname the name of each explaanatory variables
#' @param sig significance level (default 0.95)
#' @return the list of orginial result of classification, likelihood improvment and harvested classification result.
#' @export
#' @import stats

harfunc <- function(rpart.object, data, varname,sig=0.95) 
{ 
  rownn <- as.numeric(rownames(rpart.object$frame)[rpart.object$where]) 
  if(length(rownn)!=nrow(data)){
    delete.index <- (1:nrow(data))[which(!(1:nrow(data)) %in% names(rpart.object$where))]
    data <- data[-delete.index,] 
  }
  # varname <- all.vars(rpart.object$call$formula)[-1] 
  rname <- all.vars(rpart.object$call$formula)[1]    
  colnames(data)[colnames(data)==rname]<-"y"  
  predited.value <- attr(rpart.object,"ylevels")
  active_name <- predited.value[2]
  active_index <- which(data$y==active_name)
  nonact_index <- which(data$y!=active_name)
  data$y <- as.numeric(data$y)
  data[active_index,"y"] <- 1
  data[nonact_index,"y"] <- 0  
  newdata <- cbind(data,rownn) 
  termnodes <- rpart.object$frame[rpart.object$frame$var=="<leaf>",c(2,4,5)]
  term.nodes <- termnodes[,-3]
  
  term.nodes2 <- numeric(length(term.nodes[,1]))
  nodename <- as.numeric(dimnames(term.nodes)[[1]])
  
  for(i in 1:length(term.nodes2)){
    term.nodes2[i]=sum(newdata[,1]==1 & newdata[,length(newdata)]==nodename[i])
  }
  term.nodes[,2]=term.nodes2
  
  
  nodenumb <- rownames(term.nodes) 
  nci.try <- vector("list", length(nodenumb))
  for (i in 1:length(nodenumb)){
    pathr <- path.rpart.new( rpart.object, as.numeric(nodenumb)[i])
    nci.try[[i]] <- list(bounds=extrule(pathr, varname),total=term.nodes[i,1],active=term.nodes[i,2],label=nodenumb[i]) 
  }
  original_set <- nci.try
  original_set <- lapply(original_set,function(x){
    names(x)[3] <- eval(parse(text='active_name'))
    x})
  modif.b <-  vector("list", length(nci.try))
  
  
  cutoff=numeric(length(varname))
  for(i in 1:length(varname)){
    cutoff[i]=-0.5*qchisq(p=sig,i)
  }
  
  
  for(i in 1:length(nci.try)){
    noden <- nci.try[[i]] 
    ruleset <- rulesets(noden, newdata, varname, nodenumb) 
    
    modif <- change.log(ruleset, noden, nci.try, newdata)
    modif.b[[i]] <- drulenew(ruleset, modif, nci.try, sig=sig, newdata, noden,cutoff) 
  } 
  har.rule <- vector('list', 1)
  modif.a <- modif.b
  nn <- length(nci.try)
  
  for (i in 1:nn){
    print(i)
    har.r <- iterlog(modif.a, nci.try, newdata, varname) 
    
    if (har.r$delta != -Inf & har.r$flag ==0){
      logvec <- logivec(har.r$har.rule, har.r$label, nodenumb, newdata)
      harv <- newdata[logvec==T,] 
      har.rule[[i]] <- list(rule=har.r$har.rule, total=nrow(harv), active=sum(harv$y), logchange=har.r$delta)
      newdata <- newdata[logvec==F,] 
      if (nrow(newdata)==0 || sum(newdata$rownn %in% nodenumb)==0) break
      nci.try <- simu.new(nci.try, harv)  
      modif.a <- vector("list", length(nci.try))
      nn=length(nci.try)
      for (j in 1:nn){
        noden <- nci.try[[j]]
        ruleset <- rulesets(noden, newdata, varname,nodenumb)
        modif <- change.log(ruleset, noden, nci.try, newdata)
        modif.a[[j]] <- drulenew(ruleset, modif, nci.try, sig=sig, newdata, noden,cutoff) 
      }
    } 
    
    else if(har.r$flag == 1){
      har.rule[[i]] <- list(rule=har.r$har.rule, total=nrow(newdata), active=sum(newdata[,1]), logchange=har.r$delta)
      
      break
    }
    else break
  }
  har.rule <- lapply(har.rule,function(x){
    names(x)[3] <- eval(parse(text='active_name'))
    x})
  
  return(list(original_set=original_set, modif=modif.b,har.rule=har.rule))
}
