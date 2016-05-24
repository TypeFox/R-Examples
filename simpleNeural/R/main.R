nnCostFunction = function(Theta1,Theta2,X,y,lambda) {
  m = nrow(X);
  # Feedforward the neural network
  predictions = sN.MLPpredict(nnModel=list(Theta1=Theta1,Theta2=Theta2),X=X,raw=TRUE);
  
  # put observed outcome into a matrix so that we can compare it to the prediction matrix
  # ie convert y=2 into sth like y=[0 1 0 0]
  # cf http://stackoverflow.com/questions/3963565/matlab-how-to-assign-values-on-the-diagonal
  Y = zeros(max(y)+1,length(y)); # max(y) is to determine the number of classes. +1 is because we start classes at 0
  startIdx = seq(from=0,to=(max(y)+1)*length(y)-1,by=max(y)+1);
  Y[startIdx+t(y)+1] = 1;
  Y = t(Y);
  
  # compute cost
  # without regularization
  # J=1/m*sum(sum(-Y*log(predictions)-(1-Y)*log(1-predictions)));
  # with regularization (NB: we remove first column of Thetas because no penalization on the +1 unit)
  J=1/m*sum(sum(-Y*log(predictions)-(1-Y)*log(1-predictions))) + lambda/(2*m)*(sum(Theta1[,-1]^2) + sum(Theta2[,-1]^2));
  
  # now compute gradients
  X2 = cbind(ones(m, 1), X); # add +1 unit
  Z2 = X2%*%t(Theta1);
  A2 = sigmoid(Z2);
  
  A2 = cbind(ones(nrow(A2), 1), A2); # add the +1 unit
  Z3 = A2%*%t(Theta2);
  A3 = sigmoid(Z3);
  
  Delta3 = A3-Y;
  Delta2 = Delta3%*%Theta2; # just the first step to compute Delta2
  Delta2 = Delta2[,-1] * sigmoidGradient(Z2); # remove delta corresponding to bias unit, plus first step
  
  Theta2_grad=1/m*(t(Delta3)%*%A2);
  Theta1_grad=1/m*(t(Delta2)%*%X2);
  # with regularization (again, don't forget to ignore the +1 units)
  Theta1_grad[,-1]=Theta1_grad[,-1]+lambda/m*Theta1[,-1];
  Theta2_grad[,-1]=Theta2_grad[,-1]+lambda/m*Theta2[,-1];
  
  return(list(J=J,grad=list(Theta1_grad=Theta1_grad,Theta2_grad=Theta2_grad)));
}


# A dummy exported function
# @param hellow - say hello
# @return A matrix of zeros with nrows rows and ncols columns
# @export
# helloWorld = function(hellow='Hello World') {
#   return(hellow);
# }

#' Trains a multilayer perceptron with 1 hidden layer
#' @description Trains a multilayer perceptron with 1 hidden layer and a sigmoid activation function,
#' using backpropagation and gradient descent.
#' Don't forget to normalize the data first - sN.normalizeDF(), provided in the package, can be used to do so.
#' @param X Matrix of predictors
#' @param y Vector of output (the ANN learns y=ANN(X)).
#' Classes should be assigned an integer number, starting at 0 for the first class.
#' @param hidden_layer_size Number of units in the hidden layer
#' @param it Number of iterations for the gradient descent.
#' The default value of 50 may be a little low in some cases. 100 to 1000 are generally sensible values.
#' @param lambda Penalization for model coefficients (regularization parameter)
#' @param alpha Speed multiplier (learning rate) for gradient descent
#' @return The coefficients of the MLP, in a list (Theta1 between input and hidden layers, Theta2 between hidden and output layers)
#' @examples
#' # NB: the provided examples are just here to help use the package's functions.
#' # In real use cases you should perform a proper validation (cross-validation,
#' # external validation data...)
#' data(UCI.BCD.Wisconsin);
#' X=as.matrix(sN.normalizeDF(as.data.frame(UCI.BCD.Wisconsin[,3:32])));
#' y=as.matrix(UCI.BCD.Wisconsin[,2]);
#' myMLP=sN.MLPtrain(X=X,y=y,hidden_layer_size=20,it=50,lambda=0.5,alpha=0.5);
#' myPrediction=sN.MLPpredict(nnModel=myMLP,X=X,raw=TRUE);
#' #library('verification');
#' #roc.area(y,myPrediction[,2]);
#' @references M.W Gardner, S.R Dorling, Artificial neural networks (the multilayer perceptron)-
#' a review of applications in the atmospheric sciences, Atmospheric Environment, Volume 32,
#' Issues 14-15, 1 August 1998, Pages 2627-2636, ISSN 1352-2310, doi: 10.1016/S1352-2310(97)00447-0
#' [\url{http://www.sciencedirect.com/science/article/pii/S1352231097004470}]
#' @references Jain, A.K.; Jianchang Mao; Mohiuddin, K.M., "Artificial neural networks: a tutorial,"
#' Computer , vol.29, no.3, pp.31,44, Mar 1996. doi: 10.1109/2.485891
#' [\url{http://ieeexplore.ieee.org/stamp/stamp.jsp?tp=&arnumber=485891&isnumber=10412}]
#' @references Rumelhart, David E., Geoffrey E. Hinton, and R. J. Williams.
#' "Learning Internal Representations by Error Propagation". David E. Rumelhart, James L. McClelland, and
#' the PDP research group (editors).
#' Parallel distributed processing: Explorations in the microstructure of cognition, Volume 1: Foundations. MIT Press, 1986.
#' @export
sN.MLPtrain = function (X, y, hidden_layer_size=5, it=50, lambda=0.5, alpha=0.5) {
  X=as.matrix(X);
  y=as.matrix(y,length(y),1);
  input_layer_size=ncol(X);
  num_labels=max(y)+1;
  
  Theta1=randInitializeWeights(L_in=input_layer_size, L_out=hidden_layer_size);
  Theta2=randInitializeWeights(L_in=hidden_layer_size, L_out=num_labels);
  for (step in 1:it) {
    thisCost=nnCostFunction(Theta1,Theta2,X,y,lambda);
    Theta1=Theta1-alpha*thisCost[['grad']][['Theta1_grad']];
    Theta2=Theta2-alpha*thisCost[['grad']][['Theta2_grad']];
    #cat(thisCost[['J']],"\n");
  }
  
  return(list(Theta1=Theta1,Theta2=Theta2));
}

#' Runs a multilayer perceptron
#' @param nnModel A list containing the coefficients for the MLP (as produced with sN.MLPtrain())
#' @param X Matrix of predictors
#' @param raw If true, returns score of each output option. If false, returns the output option with highest value.
#' @return The predicted values obtained by the MLP
#' @examples
#' data(UCI.transfusion);
#' X=as.matrix(sN.normalizeDF(as.data.frame(UCI.transfusion[,1:4])));
#' y=as.matrix(UCI.transfusion[,5]);
#' myMLP=sN.MLPtrain(X=X,y=y,hidden_layer_size=4,it=50,lambda=0.5,alpha=0.5);
#' myPrediction=sN.MLPpredict(nnModel=myMLP,X=X,raw=TRUE);
#' #library('verification');
#' #roc.area(y,myPrediction[,2]);
#' @export
sN.MLPpredict = function(nnModel,X,raw=FALSE) {  
  m = nrow(X);
  
  X = cbind(ones(m, 1), X); # Add ones to the X data matrix
  h1 = sigmoid(X%*%t(nnModel$Theta1)); # Compute hidden layer
  h1 = cbind(ones(nrow(h1),1),h1); # Add ones to the hidden layer
  h2 = sigmoid(h1%*%t(nnModel$Theta2));
  if(raw) return(h2);
  # Returns which column (=which output category) is most likely. NB: for 2-class classification probably this could/should be optimized to only 1 column
  # http://stackoverflow.com/questions/743622/finding-row-index-containing-maximum-value-using-r
  return(apply(h2,1,which.max)-1);
}
  