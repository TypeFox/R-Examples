
#' Trains an online sequential extreme learning machine with random weights
#' @param p dataset used to perform the training of the model
#' @param y classes vector for classiication or regressors for regression
#' @param Elm_Type select if the ELM must perform a "regression" or "classification"
#' @param nHiddenNeurons number of neurons in the hidden layer
#' @param ActivationFunction "rbf" for radial basis function with Gaussian kernels , "sig" for sigmoidal fucntion, "sin" for sine function, "hardlim" for hard limit function
#' @param N0 size of the first block to be processed
#' @param Block size of each chunk to be processed at each step
#' @importFrom stats model.frame model.matrix model.response runif
#' @return returns all the parameters used in the function, the weight matrix, the labels for the classification, the number of classes found, the bias, the beta activation function and the accuracy on the trainingset
#' @examples
#' x = runif(100, 0, 50)
#' y = sqrt(x)
#' train = data.frame(y,x)
#' train = data.frame(preProcess(train))
#' OSelm_train.formula(y~x, train, "regression", 100, "hardlim", 10, 10)
#' @references
#' [1] N.-Y. Liang, G.-B. Huang, P. Saratchandran, and N. Sundararajan, 'A Fast and Accurate On-line Sequential Learning Algorithm for Feedforward Networks' IEEE Transactions on Neural Networks, vol. 17, no. 6, pp. 1411-1423, 2006
#' @export
OSelm_training = function(p, y,Elm_Type, nHiddenNeurons, ActivationFunction, N0, Block){
  p = matrix(p, ncol = ncol(p))
  #tv.t = testing[,1]
  #tv.p = testing[,-1]
  t = as.matrix(as.numeric(y), ncol = ncol(y))
  #tv.t = as.numeric(tv.t)

  nTrainingData = dim(p)[1]
  #nTestingData = dim(tv.p)[1]
  nInputNeurons= dim(p)[2]


  ###classification###
  if(Elm_Type == "classification"){
    #label = unique(c(t,tv.t))
    label = unique(t)
    nClass = nOutputNeurons = length(label)

    temp_T = matrix(0, nrow = nTrainingData, ncol = nClass)
    for(i in 1:nTrainingData){
      for(j in 1:nClass){
        if(label[j] == t[i]){
          break;
        }
      }
      temp_T[i,j]=1;
    }
    t=temp_T*2-1;

#     temp_TV_T = matrix(0, nrow = nTestingData, ncol= nClass)
#     for(i in 1:nTestingData){
#       for(j in 1:nClass){
#         if(label[j] == tv.t[i]){
#           break;
#         }
#       }
#       temp_TV_T[i,j]=1;
#     }
#     tv.t=temp_TV_T*2-1;
  }

  #start time

  p0 = matrix(p[1:N0,], ncol = dim(p)[2])
  t0 = t[1:N0,]

  iw = matrix(runif(nHiddenNeurons*nInputNeurons), ncol=nInputNeurons)*2-1
  switch(tolower(ActivationFunction),
         rbf={
           bias = runif(nHiddenNeurons);
           H0 = RBFun(p0,iw,bias); #todo
         },
         sig={
           bias = runif(nHiddenNeurons)*2-1;
           H0 = sigActFun(p0,iw,bias); #todo
         },
         sin={
           bias = runif(nHiddenNeurons)*2-1;
           H0 = sinActFun(p0,iw,bias); #todo
         },
         hardlim={
           bias = runif(nHiddenNeurons)*2-1;
           H0 = hardLimActFun(p0,iw,bias); #todo
         })
  M = myPinv(t(H0) %*% H0)
  beta = myPinv(H0) %*% t0
  rm(p0,t0,H0)

  for(n in seq(N0, nTrainingData, by = Block)){
    if ((n+Block-1) > nTrainingData){
     pn = matrix(p[n:nTrainingData,], ncol = dim(p)[2])
     tn = t[n:nTrainingData]  #check, must be nTestData
     Block = dim(pn)[1]
    }
    else{
      pn = matrix(p[n:(n+Block-1),], ncol = dim(p)[2])
      tn = t[n:(n+Block-1),]
    }
    switch(tolower(ActivationFunction),
           rbf={
             H = RBFun(pn,iw,bias); #todo
           },
           sig={
             H = sigActFun(pn,iw,bias); #todo
           },
           sin={
             H = sinActFun(pn,iw,bias); #todo
           },
           hardlim={
             H = hardLimActFun(pn,iw,bias); #todo
           })
    M = M - M %*% t(H) %*% ((diag(Block) + H %*% M %*% t(H))%^%(-1)) %*% H %*% M
    beta = beta + M %*% t(H) %*% (tn - H %*% beta)
  }

  #stop time
  #clear Pn Tn H M;
  switch(tolower(ActivationFunction),
         rbf={
           HTrain = RBFun(p, iw, bias); #todo
         },
         sig={
           HTrain = sigActFun(p, iw, bias); #todo
         },
         sin={
           HTrain = sinActFun(p, iw, bias); #todo
         },
         hardlim={
           HTrain = hardLimActFun(p, iw, bias); #todo
         })
  Y=HTrain %*% beta;
  #clear HTrain;

  if(Elm_Type == "regression"){
    ### Calculate RMSE in the case of REGRESSION ###
    TrainingAccuracy=sqrt(mean(sum((t - Y)^2)))
    return(list(ActivationFunction=ActivationFunction, Elm_Type = Elm_Type,iw=iw, bias=bias, beta=beta, TrainingAccuracy=TrainingAccuracy))
  }
  else if(Elm_Type == "classification"){
    ### Calculate correct classification rate in the case of CLASSIFICATION ###
    MissClassificationRate_Training=0;
    label_index_expected = array()
    label_index_actual = array()
    for(i in 1:nTrainingData){
       label_index_expected[i] = which.max(t[i,])
       label_index_actual[i] = which.max(Y[i,])
       #if(label_index_expected != label_index_actual){
        # MissClassificationRate_Training=MissClassificationRate_Training+1;
       #}
    }
    t = table(observed = label_index_expected, predicted = label_index_actual)
    accuracy.reg <- sum(diag(t)) / sum(t)
    #TrainingAccuracy=1-MissClassificationRate_Training/nTrainingData
    return(list(ActivationFunction=ActivationFunction, Elm_Type = Elm_Type,iw=iw, label=label, nClass=nClass,bias=bias, beta=beta, TrainingAccuracy=accuracy.reg))
  }
}


#' Prediction function for the ELM model generated with the elm_training() function
#' @param model the output of the elm_training() function
#' @param test dataset used to perform the testing of the model, the first column must be the column to be fitted for the regression or the labels for the classification
#' @return returns the accuracy on the testset
#' @examples
#' x = runif(100, 0, 50)
#' y = sqrt(x)
#' train = data.frame(y,x)
#' train = data.frame(preProcess(train))
#' model = OSelm_train.formula(y~x, train, "regression", 100, "hardlim", 10, 10)
#' #' x = runif(100, 0, 50)
#' y = sqrt(x)
#' test = data.frame(y,x)
#' test = data.frame(preProcess(train))
#' accuracy = predict_elm(model, test)
#' @references
#' [1] N.-Y. Liang, G.-B. Huang, P. Saratchandran, and N. Sundararajan, "A Fast and Accurate On-line Sequential Learning Algorithm for Feedforward Networks" IEEE Transactions on Neural Networks, vol. 17, no. 6, pp. 1411-1423, 2006
#' @export
predict_elm = function(model, test){
  if(!is.null(model$formula)){
    mf <- model.frame(formula = model$formula, data = test)
    x <- model.matrix(attr(mf, "terms"), data = mf)
    y <- model.response(mf)
  }

  tv.p = matrix(x, ncol = ncol(x))
  #tv.t = testing[,1]
  #tv.p = testing[,-1]
  tv.t = as.matrix(as.numeric(y), ncol = ncol(y))

  #tv.t = test[,1]
  #tv.p = matrix(test[,-1], ncol = dim(test)[2]-1)
  #tv.t = as.numeric(tv.t)
  nTestingData = dim(tv.p)[1]
  ActivationFunction = model$ActivationFunction
  iw = model$iw
  bias = model$bias
  beta = model$beta
  Elm_Type = model$Elm_Type
  nClass = model$nClass
  label = model$label
  if(Elm_Type == "classification"){
    temp_TV_T = matrix(0, nrow = nTestingData, ncol= nClass)
    for(i in 1:nTestingData){
      for(j in 1:nClass){
        if(label[j] == tv.t[i]){
          break;
        }
      }
      temp_TV_T[i,j]=1;
    }
    tv.t=temp_TV_T*2-1;
  }


  ##### Performance Evaluation #####

  #start time#
  switch(tolower(ActivationFunction),
         rbf={
           HTest = RBFun(tv.p, iw, bias) #todo
         },
         sig={
           HTest = sigActFun(tv.p, iw, bias); #todo
         },
         sin={
           HTest = sinActFun(tv.p, iw, bias); #todo
         },
         hardlim={
           HTest = hardLimActFun(tv.p, iw, bias); #todo
         })

  TY=HTest %*% beta;
  #clear HTest;
  #stop time#
  #TestingTime=end_time_test-start_time_test

  if(Elm_Type == "regression"){
    ### Calculate RMSE in the case of REGRESSION ###
    TestingAccuracy=sqrt(mean(sum((tv.t - TY)^2)))
    return(list(predicted = TY, testAccuracy = TestingAccuracy))
  }
  else if(Elm_Type == "classification"){
    ### Calculate correct classification rate in the case of CLASSIFICATION ###
    MissClassificationRate_Testing=0;
    label_index_expected = array()
    label_index_actual = array()
    for(i in 1:nTestingData){
      label_index_expected[i] = which.max(tv.t[i,])
      label_index_actual[i] = which.max(TY[i,])
      #if(label_index_expected != label_index_actual){
      #  MissClassificationRate_Testing=MissClassificationRate_Testing+1;
      #}
    }
    t = table(observed = label_index_expected, predicted = label_index_actual)
    accuracy.reg <- sum(diag(t)) / sum(t)
    #TestingAccuracy=1-MissClassificationRate_Testing/nTestingData
    return(list(prediction=label_index_actual,accuracy.reg=accuracy.reg))
  }
}
