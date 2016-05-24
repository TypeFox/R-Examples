"vim.interactions" <-
function(x, pred.data, x.pair, w, learner=c("trees"),verbose=FALSE,nseq=20)  	# use weights for samples 
{
if(!learner%in%c("tree","mars","tensor.product","linear.regression"))
stop("learner", learner, "is not implemented yet\n")
if(any(x.pair[,1]==x.pair[,2]))
stop("The paired interaction should have different input variable index (names)\n")


  if(learner=="linear.regression") pred.data <- pred.data[,-1] #remove the first columen, which has 1s for intercept
  n.preds <- ncol(pred.data)
  if(learner=="linear.regression") n.preds <- -0.5+0.5*sqrt(1+8*n.preds)  #For linear regression, the first column is the intercept; p+p*(p-1)/2 = d for main + interaction terms
  pred.names <- colnames(pred.data)[1:n.preds]
  cross.tab <- matrix(0,ncol=n.preds,nrow=n.preds)
  dimnames(cross.tab) <- list(pred.names,pred.names)


   data <- pred.data
   for (s in 1:nrow(x.pair)){
    i <- x.pair[s,1]
    if (is.vector(data[,i])) {  # create a sequence through the range
       x.var <- seq(min(data[,i],na.rm=T),max(data[,i],na.rm=T),length = nseq)
       }
    else {                      # otherwise set up simple factor variable
       x.var <- factor(names(table(data[,i])),levels = levels(data[,i]))
       }
    x.length <- length(x.var)
    if(verbose)
    cat(i,"")

    j <- x.pair[s,2]     
    if(verbose)
    cat(j,"")
      if (is.vector(data[,j])) {
        y.var <- seq(min(data[,j],na.rm=T),max(data[,j],na.rm=T),length = nseq)
      }
      else {
        y.var <- factor(names(table(data[,j])),levels = levels(data[,j]))
      }
      y.length <- length(y.var)


      tmp <- expand.grid(list(x.var,y.var))
      pred.frame <- matrix(NA,nrow=nrow(tmp),ncol=n.preds)
      pred.frame[,i] <- tmp[[1]]
      pred.frame[,j] <- tmp[[2]]
      nreplace <- ifelse(learner=="tensor.product",floor(nrow(pred.frame)*3/5),nrow(pred.frame))

      for (k in 1:n.preds) {
        if (k != i & k != j) {
          if (is.vector(data[,k])) {  # either with the mean
            pred.frame[1:nreplace,k] <- mean(data[,k],na.rm=T)
          }
          else {   # or the most common factor level
            temp.table <- sort(table(data[,k]),decreasing = TRUE)
            pred.frame[1:nreplace,k] <- rep(names(temp.table)[1],x.length * y.length)
            pred.frame[1:nreplace,k] <- as.factor(pred.frame[,k])  #perhaps needs more work here
          }
        }
      }
      colnames(pred.frame) <- pred.names
      if(learner=="tree")
      prediction <- predict.gbm(x,as.data.frame(pred.frame),n.trees = x$n.trees, type="link")
      else if(learner=="tensor.product" || learner=="mars")
      prediction <- predict(x,newdata=pred.frame)
           else {x.data <- pred.frame   #for linear regression with interactions
                for(jj in 1:(n.preds-1))
                  for (kk in (jj+1):n.preds)
                    x.data <- cbind(x.data,x.data[,jj]*x.data[,kk])
                x.data <- cbind(1,x.data) # add the first column for intercept; at least for glmboost, the prediction is based on column location, not on column names.
                prediction <- predict(x,newdata=as.matrix(x.data))
}
      if (missing(w)) w <- rep(1, nrow(pred.frame))
        interaction.test.model <- lm(prediction[1:nreplace] ~ as.factor(pred.frame[1:nreplace,i]) + as.factor(pred.frame[1:nreplace,j]), 
          weights = w[1:nreplace])
  
        
      interaction.flag <- mean(resid(interaction.test.model)^2)
      cross.tab[i,j] <- interaction.flag
      if(verbose) cat(" ",interaction.flag,"\n")
    }   # end of j loop

  search.index <- ((n.preds^2) + 1) - rank(cross.tab, ties.method = "first")

  n.important <- max(2,round(0.1 * ((n.preds^2)/2),0))
  var1.names <- rep(" ",n.important)
  var1.index <- rep(0,n.important)
  var2.names <- rep(" ",n.important)
  var2.index <- rep(0,n.important)
  int.size <- rep(0,n.important)

  for (i in 1:n.important) {
  
    index.match <- match(i,search.index)

    j <- min(n.preds,trunc(index.match/n.preds) + 1)
    var1.index[i] <- j
    var1.names[i] <- pred.names[j]
    k <- index.match%%n.preds
    if (k > 0) {   #only do this if k > 0 - otherwise we have all zeros from here on 
      var2.index[i] <- k
      var2.names[i] <- pred.names[k]
      int.size[i] <- cross.tab[k,j]
    }
  }
  cross.tab <- round(100*cross.tab/sum(cross.tab,na.rm=TRUE),2)
  rank.list <- data.frame(var1.index,var1.names,var2.index,var2.names,int.size)
  return(list(rank.list = rank.list, interactions = cross.tab))
}

