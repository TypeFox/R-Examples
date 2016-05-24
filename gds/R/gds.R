#' Descriptive statistics of grouped data: with the help of this package we calculate mean, median, mode, variance, standard deviation, coefficient of variance, quartiles, IQR, skewness, and kurtosis of grouped data.
#'@param ll A data vector to store lower limit of the classes
#'@param ul A data vector to store upper limit of the classes
#'@param freq A data vector to store the frequencies of the corresponding classes
#'
#'@examples
#' gds(c(10,20,30,40,50),c(20,30,40,50,60),c(7,13,23,20,8))
#'
#'
#' @references
#' 1. Gupta, S.P., and Gupta, M.P. (2005) Business statistics, Sultan Chand and Sons educational publishers, New Delhi.
#'
#' 2. Levine, D.M., Krehbiel,  T.C.,  Bereson, M.L. and Viswanathan, P.K. (2011) Business statistics: a first course, 5th edition, Pearson.
#'
#' 3. Langford, E. (2006) Quartiles in Elementary Statistics, Journal of Statistics Education Volume 14, Number 3.
#'
#' 4. Das, N. G. (2010) Statistical Methods- Combined Edition (Volumes I & II), Tata McGraw Hill Education Private Limited, New Delhi.
#'
#'@return gmean, gmedian, gmode, gvar, gstdev, gcv, gq1, gq2, gq3, gIQR, g1, g2
#'@export

gds<-function(ll,ul,freq){

  input<-data.frame(ll,ul,freq)
  colnames(input)<-c("lowerLimit","upperLimit","frequency")
  cat("The input data is: \n"); print(input); cat("\n")

  cat("For the given grouped data \n")
  ###############################################
  # The grouped mean calculation
  input1 <- cbind(input, (input[,1] + input[,2])/2) #identification of mid point
  gmean<-sum(input1[,3]*input1[,4])/sum(input1[,3]) # grouped mean calculation
  cat("The mean is: ",gmean)
  cat("\n")

  ###############################################
  # The grouped median calculation
  input2 <- cbind(input, cumsum(input[,3]))		# display the entered data
  # Identification of median class
  mClass<-input2[input2[,4]==input2[,4][input2[,4]>sum(input2[,3])/2][1],]
  # Lower limit of the median class
  L<-mClass[1,1]
  # Median class frequency
  f<-mClass[1,3]
  # Class interval
  i<-(mClass[1,2]-mClass[1,1])
  # To identification of the p.c.f
  rowLength<-matrix(dim(input2[input2[,4]<sum(input2[,3])/2,]))[1,1]
  # The p.c.f value
  pcf<-input2[which(input2[,4]==mClass[,4])-1,][,4]
  N<- sum(input2[ ,3]) # Total frequency
  gmedian<-L + (N/2 - pcf)*i/f # Median calculation
  cat("The median is: ",gmedian)
  cat("\n")

  ###############################################
  # The grouped mode calculation
  input3<-input
  modalClass<-input3[which(input3[,3]==max(input3[,3])),] # modal class identification
  preModalClass<-(input3[which(input3[,3]==max(input3[,3]))-1,]) # pre modal class
  postModalClass<-(input3[which(input3[,3]==max(input3[,3]))+1,]) # post modal class
  D1<-modalClass[,3]-preModalClass[,3]
  D2<-modalClass[,3]-postModalClass[,3]
  L<-modalClass[,1] # identification of lower limit of modal class
  i<-modalClass[,2]-modalClass[,1] # modal class interval
  gmode = L + D1*i/(D1+D2) # mode calculation
  cat("The mode is: ",gmode)
  cat("\n")

  ###############################################
  # The grouped standard deviation calculation
  input4<-input
  x<-(input4[,2]+ input4[,1])/2
  n<-sum(input4[,3]) # for sample std deviation n<-sum(input4[,3]) - 1
  xMean= sum(input4[,3]*x)/n
  gvar<-sum(input4[,3]*(x-xMean)^2)/n
  gstdev<-(gvar)^(1/2)
  cat("The variance is: ", gvar)
  cat("\n")
  cat("The standard deviation is: ", gstdev)
  cat("\n")

  # The coefficients of variation or relative standard deviation
  gcv<-(gstdev/gmean)*100
  cat("The coefficient of variance (around the mean) is: ", gcv);
  cat("\n")

  ###############################################
  input5<-input
  input5 <- cbind(input5, cumsum(input5[,3]))
  #k=1 for Q1, k=2 for Q2, and k=3 for Q3
  for(k in 1:3)
  {
    mClass<-input5[input5[,4]==input5[,4][input5[,4]>sum(input5[,3])*k/4][1],]
    L<-mClass[1,1] 					# to get lower limit of the kth quartile
    f<-mClass[1,3] 					# to get frequency the kth quartile
    i<-(mClass[1,2]-mClass[1,1]) 			# to get class interval of the kth quartile
    N<- sum(input5[ ,3]) 				# total frequency
    pcf<-input5[which(input5[,4]==mClass[,4])-1,][,4] 	# to get pcf of the kth quartile
    quantile<-L + (N*k/4 - pcf)*i/f 			# Quartiles calculation

    if(k ==1){
      cat("The 1st quartile is: ", quantile) 			# to print the result
      gq1<-quantile
      cat("\n")}

    if(k ==2){
      cat("The 2nd quartile is: ", quantile) 			# to print the result
      gq2<-quantile
      cat("\n")}

    if(k ==3){
      cat("The 3rd quartile is: ", quantile)			# to print the result
      gq3<-quantile }
  }
  gIQR<- (gq3-gq1)
  cat("\nThe inter-quartile range is: ",gIQR)
  cat("\n")

  ###############################################
  #Calculation of skewness and kutosis using group moments
  input6<-input
  x<-(input6[,2]+input6[,1])/2 	# identify the mid point
  N<-sum(input6[,3]) 		# total number of data
  f<-input6[,3] 			# frequency of each class
  xMean= sum(input6[,3]*x)/N 	# mean
  mu1<-sum(f*(x-xMean))/N # mu1
  mu2<- sum(f*(x-xMean)^2)/N  # mu2
  mu3<- sum(f*(x-xMean)^3)/N  # mu3
  mu4<- sum(f*(x-xMean)^4)/N  # mu4
  b1<-mu3^2/mu2^3 #

  # The skewness calculation
  g1<-b1^(1/2) #  The skewness

  # The kurtosis calculation
  b2<-mu4/mu2^2 #
  g2<- b2-3 # the kurtosis
  cat("The skewness is: ", g1); cat("\n")
  cat("The kurtosis is: ", g2); cat("\n")


  ###############################################
  # returning section
  glist <- list("mean" = gmean, "median" = gmedian, "mode" = gmode, "var" = gvar, "stdDev" = gstdev, "cv" = gcv, "quartile1"= gq1, "quartile2"= gq2, "quartile3"= gq3, "IQR" = gIQR, "skewness"= g1, "kurtosis"= g2)
  return(glist)

}
