#Example 5
#Warning, it will take a while

#Load the Jersey dataset
data(Jersey)

#Predictive power of the model using the SECOND set for 10 fold CROSS-VALIDATION
data=pheno
data$G=G
data$D=D
data$partitions=partitions

#Fit the model for the TESTING DATA for Additive + Dominant
out=brnn_extended(yield_devMilk ~ G | D,
				  data=subset(data,partitions!=2),
				  neurons1=2,neurons2=2,epochs=100,verbose=TRUE)

#Plot the results
#Predicted vs observed values for the training set
par(mfrow=c(2,1))
yhat_R_training=predict(out)
plot(out$y,yhat_R_training,xlab=expression(hat(y)),ylab="y")
cor(out$y,yhat_R_training)
  
#Predicted vs observed values for the testing set
newdata=subset(data,partitions==2,select=c(D,G))
ytesting=pheno$yield_devMilk[partitions==2]
yhat_R_testing=predict(out,newdata=newdata)
plot(ytesting,yhat_R_testing,xlab=expression(hat(y)),ylab="y")
cor(ytesting,yhat_R_testing)
