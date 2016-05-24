#
# gbm.interactions version 2.9
#
# j. leathwick, j. elith - May 2007
#
# functions assesses the magnitude of 2nd order interaction effects 
# in gbm models fitted with interaction depths greater than 1
# this is achieved by:
#   1. forming predictions on the linear scale for each predictor pair;
#   2. fitting a linear model that relates these predictions to the predictor
#        pair, with the the predictors fitted as factors;
#   3. calculating the mean value of the residuals, the magnitude of which
#        increases with the strength of any interaction effect;
#   4. results are stored in an array;
#   5. finally, the n most important interactions are identified,
#        where n is 25% of the number of interaction pairs;




gbm.interactions <- function(gbm.object, 
   use.weights = FALSE,   	# use weights for samples 
   mask.object)       		# a gbm object describing sample intensity 
{

    if (! requireNamespace('gbm') ) { stop ('you need to install the gbm package to run this function') }
	requireNamespace('splines')
	
  gbm.call <- gbm.object$gbm.call
  n.trees <- gbm.call$best.trees
  depth <- gbm.call$interaction.depth
  gbm.x <- gbm.call$gbm.x
  n.preds <- length(gbm.x)
  pred.names <- gbm.object$gbm.call$predictor.names
  cross.tab <- matrix(0,ncol=n.preds,nrow=n.preds)
  dimnames(cross.tab) <- list(pred.names,pred.names)

  if (use.weights) mask.trees <- mask.object$gbm.call$best.trees

  message("gbm.interactions - version 2.9\nCross tabulating interactions for gbm model with ", n.preds, " predictors")

  data <- gbm.call$dataframe[,gbm.x]  

  for (i in 1:(n.preds - 1)) {  # step through the predictor set

    if (is.vector(data[,i])) {  # create a sequence through the range
       x.var <- seq(min(data[,i],na.rm=T),max(data[,i],na.rm=T),length = 20)
       }
    else {                      # otherwise set up simple factor variable
       x.var <- factor(names(table(data[,i])),levels = levels(data[,i]))
       }
    x.length <- length(x.var)

    message(i, " ", appendLF=FALSE)

    for (j in (i+1):n.preds) { #create vector or factor data for second variable
      
      if (is.vector(data[,j])) {
        y.var <- seq(min(data[,j],na.rm=T),max(data[,j],na.rm=T),length = 20)
      }
      else {
        y.var <- factor(names(table(data[,j])),levels = levels(data[,j]))
      }
      y.length <- length(y.var)

# and now make a temporary data frame

      pred.frame <- expand.grid(list(x.var,y.var))
      names(pred.frame) <- c(pred.names[i],pred.names[j])

      n <- 3 # and add the balance of the variables to it

      for (k in 1:n.preds) {
        if (k != i & k != j) {
          if (is.vector(data[,k])) {  # either with the mean
            pred.frame[,n] <- mean(data[,k],na.rm=T)
          }
          else {   # or the most common factor level
            temp.table <- sort(table(data[,k]),decreasing = TRUE)
            pred.frame[,n] <- rep(names(temp.table)[1],x.length * y.length)
            pred.frame[,n] <- as.factor(pred.frame[,n])
          }
          names(pred.frame)[n] <- pred.names[k]
          n <- n + 1
        }
      }
#
# form the prediction
#
      prediction <- gbm::predict.gbm(gbm.object,pred.frame,n.trees = n.trees, type="link")

      if (use.weights) {
        point.prob <- gbm::predict.gbm(mask.object[[1]],pred.frame, n.trees = mask.trees, type="response")
        interaction.test.model <- lm(prediction ~ as.factor(pred.frame[,1]) + as.factor(pred.frame[,2]), 
          weights = point.prob)
      }
  
      else {

        interaction.test.model <- lm(prediction ~ as.factor(pred.frame[,1]) + as.factor(pred.frame[,2]))
      }
        
      interaction.flag <- round(mean(resid(interaction.test.model)^2) * 1000,2)

      cross.tab[i,j] <- interaction.flag

    }   # end of j loop
  }  # end of i loop

# create an index of the values in descending order

  search.index <- ((n.preds^2) + 1) - rank(cross.tab, ties.method = "first")

  n.important <- max(2,round(0.1 * ((n.preds^2)/2),0))
  var1.names <- rep(" ",n.important)
  var1.index <- rep(0,n.important)
  var2.names <- rep(" ",n.important)
  var2.index <- rep(0,n.important)
  int.size <- rep(0,n.important)

  for (i in 1:n.important) {
  
    index.match <- match(i,search.index)

    j <- trunc(index.match/n.preds) + 1
    var1.index[i] <- j
    var1.names[i] <- pred.names[j]

    k <- index.match%%n.preds
    if (k > 0) {   #only do this if k > 0 - otherwise we have all zeros from here on 
      var2.index[i] <- k
      var2.names[i] <- pred.names[k]

      int.size[i] <- cross.tab[k,j]
    }

  }

  rank.list <- data.frame(var1.index,var1.names,var2.index,var2.names,int.size)

  message("")
  return(list(rank.list = rank.list, interactions = cross.tab, gbm.call = gbm.object$gbm.call))
}



