##############################################################
#Example 4
#Gianola et al. (2011).
#Warning, it will take a while, substitute the FALSE
#statement with the TRUE statement if you really want to run the example

#This example uses the formula interface for the fitting
if(FALSE)
{
  #Load the Jersey dataset
  data(Jersey)
  
  #Fit the model with the FULL DATA
  #Formula interface
  out=brnn(pheno$yield_devMilk~G,neurons=2,verbose=TRUE)
  
  #Obtain predictions and plot them against fitted values
  plot(pheno$yield_devMilk,predict(out))
  
  #Predictive power of the model using the SECOND set for 10 fold CROSS-VALIDATION
  data=pheno
  data$X=G
  data$partitions=partitions
  
  #Fit the model for the TESTING DATA
  out=brnn(yield_devMilk~X,
           data=subset(data,partitions!=2),neurons=2,verbose=TRUE)
           
  #Plot the results
  #Predicted vs observed values for the training set
  par(mfrow=c(2,1))
  plot(out$y,predict(out),xlab=expression(hat(y)),ylab="y")
  cor(out$y,predict(out))
  
  #Predicted vs observed values for the testing set
  yhat_R_testing=predict(out,newdata=subset(data,partitions==2))
  ytesting=pheno$yield_devMilk[partitions==2]
  plot(ytesting,yhat_R_testing,xlab=expression(hat(y)),ylab="y")
  cor(ytesting,yhat_R_testing)
}

#This example uses the default method for the call
if(FALSE)
{
  #Load the Jersey dataset
  data(Jersey)
  
  #Fit the model with the FULL DATA
  out=brnn(y=pheno$yield_devMilk,x=G,neurons=2,verbose=TRUE)

  #Obtain predictions and plot them against fitted values
  plot(pheno$yield_devMilk,predict(out))

  #Predictive power of the model using the SECOND set for 10 fold CROSS-VALIDATION
  index=partitions==2
  Xtraining=G[!index,]
  ytraining=pheno$yield_devMilk[!index]
  Xtesting=G[index,]
  ytesting=pheno$yield_devMilk[index]

  #Fit the model for the TESTING DATA
  out=brnn(y=ytraining,x=Xtraining,neurons=2,verbose=TRUE)

  #Plot the results
  #Predicted vs observed values for the training set
  par(mfrow=c(2,1))
  yhat_R_training=predict(out)
  plot(ytraining,yhat_R_training,xlab=expression(hat(y)),ylab="y")
  cor(ytraining,yhat_R_training)
  
  #Predicted vs observed values for the testing set
  yhat_R_testing=predict(out,Xtesting)
  plot(ytesting,yhat_R_testing,xlab=expression(hat(y)),ylab="y")
  cor(ytesting,yhat_R_testing)
}
