#' Incorporating Prior Knowledge
#' 
#' Learns a function using prior information.
#' 
#' @param priorData A \code{"numeric"} array which contains the values of the variable we have information apriori about.
#' @param data A \code{"numeric"} array which contains the values to fit.
#' @param s A \code{"numeric"} coefficient which fixes the confidence of the prior knowledge 
#' we are going to introduce. By default it is \code{NULL}, only we must modify it if we want 
#' to incorporate prior information to the fits.
#' @param POTENTIAL_TYPE A \code{"character"} string specifying the posibles potential
#' types, must be one of \code{"MOP"} or \code{"MTE"}.
#' @param domain A \code{"numeric"} array which contains the limits to defined the data function.
#' By default it is the range of the data.
#' @param coeffversion A \code{"numeric"} value between \code{1--4} which contains the used version for computing the coefficients of the linear opinion pool to
#' combine the prior function and the data function. By default \code{coeffversion = "4"} is used, so the combination
#' depends on the goodness of the model versus another random model.
#' @param restrictDomain This argument lets us choose if the domain is used joining both domains,
#' the prior one and the data domain or trimming them. By default \code{TRUE} is used, so 
#' the domain will be trimmed.
#' @param maxParam A \code{"numeric"} value which indicates the maximum number of coefficients in the function. By default it is \code{NULL}; 
#' if not, the output is the function which gets the best BIC with at most this number of parameters.
#' @return A list with the elements
#' \item{coeffs}{An \code{"numeric"} array with the two coefficients of the linear opinion pool}
#' \item{posteriorFunction}{The final function after combining.}
#' \item{priorFunction}{The fit of the prior data.}
#' \item{dataFunction}{The fit of the original data.}
#' \item{rangeNewPriorData}{A \code{"numeric"} vector which contains the final domain where the functions are defined.}
#' @seealso \link{getCoefficients}
#' @export
#' @examples
#' 
#' ## Data
#' X <- rnorm(15)
#' 
#' ## Prior Data
#' priordata <- rnorm(5000)
#' 
#' ## Test data
#' test <- rnorm(1000)
#' testData <- test[test>=min(X)&test<=max(X)]
#' 
#' ## Learning
#' type <- "MOP" ## type <- "MTE"
#' confident <- 3 ## confident <- 1,2,...,length(X)
#' f <- learnMoTBFpriorInformation(priorData = priordata, data = X, s = confident,
#' POTENTIAL_TYPE = type)
#' attributes(f)
#' 
#' ## Log-likelihood
#' sum(log(as.function(f$dataFunction)(testData)))
#' sum(log(as.function(f$posteriorFunction)(testData))) ## best loglikelihood
#' 
#'
learnMoTBFpriorInformation <- function(priorData, data, s, POTENTIAL_TYPE, domain=range(data), coeffversion=4, restrictDomain=TRUE, maxParam=NULL)
{  
  ## Learning prior function
  fPI <- univMoTBF(priorData,POTENTIAL_TYPE, range(priorData, domain))
  
  ## Computing the new domain
  if (is.null(restrictDomain)||(restrictDomain==FALSE)) rangeNewPriorData <- range(priorData)
  else rangeNewPriorData <- newRangePriorData(fPI,priorData, length(data), domain, s, POTENTIAL_TYPE)
  
  ## Normalizing the prior function
  iP <- integrate(as.function(fPI), min(domain, rangeNewPriorData), max(domain, rangeNewPriorData))$value
  if(POTENTIAL_TYPE=="MOP"){
    f=coef(fPI)/iP
    fPI <- asMOPString(f)
    fPI <- motbf(fPI)
  } else {
    f <- coef(fPI)/iP
    fPI <- asMTEString(f)
    fPI <- motbf(fPI)
  }
  
  ## Learning data function
  fD <- univMoTBF(data,POTENTIAL_TYPE, range(rangeNewPriorData, domain), maxParam=maxParam)
  
  ## Learning posterior function
  coeffs <- getCoefficients(fPI, rangeNewPriorData, fD, data, domain, coeffversion)
  
  if(coeffs[1]==0){
      fX <- fD
  } else if (coeffs[2]==0){
      fX <- fPI
  } else {
    fPI1 <- coef(fPI); fD1=coef(fD)
    maxLength <- max(length(fPI1), length(fD1))
    if(length(fD1)<maxLength){
      fD1 <- c(fD1, rep(0, maxLength-length(fD1)))
    } else {
      fPI1 <- c(fPI1, rep(0,maxLength-length(fPI1)))
    }  
    fX <- fPI1*coeffs[1]+fD1*coeffs[2]
      
    if(POTENTIAL_TYPE=="MOP") fX <- asMOPString(fX) else fX <- asMTEString(fX)
    fX <- motbf(fX)  
  }
  
  return(list(coeffs=coeffs, posteriorFunction=fX, 
              priorFunction=fPI, dataFunction=fD, 
              domain=range(domain, rangeNewPriorData)))
}

#' Redefining the Domain
#' 
#' Computes the new domain of two datasets. 
#' 
#' @param fPI The fitted function to the prior data of class \code{"motbf"}.
#' @param priorData A \code{"numeric"} array with the values we want to include as prior information.
#' @param N A \code{"numeric"} value which is the size of the data.
#' @param domain A \code{"numeric"} array with the limits where defining the data function.
#' @param s A \code{"numeric"} value which is the confident of the expert in his information. It is between 0 and the data size.
#' @param POTENTIAL_TYPE A \code{"character"} string giving the potential of the model, i.e. \code{"MOP"} if the basis functions are polynomials,
#' or \code{"MTE"} if they are exponentials.
#' @return A \code{"numeric"} array which contains the new domain of the prior function.
#' @export
#' @examples
#' 
#' ## Data
#' X <- rnorm(15)
#' 
#' ## Prior Data
#' priordata <- rnorm(5000)
#' 
#' ## Learning
#' type = "MTE" ## type <- "MOP"
#' fPrior <- univMoTBF(priordata, POTENTIAL_TYPE = type)
#' 
#' ## New range
#' confident <- 5 ## confident <- 1,2,...,length(X)
#' domain <- range(X)
#' N <- length(X)
#' newRange <- newRangePriorData(fPrior, priorData = priordata, N = N,
#' domain = domain, s = confident, POTENTIAL_TYPE = type)
#' newRange
#' 
newRangePriorData <- function(fPI, priorData, N, domain, s, POTENTIAL_TYPE)
{
  if(s<N){
    coeff <- 1-s/N
    diffmin <- 0;  diffmax <- 0
    
    if((min(priorData)<min(domain))) diffmin <- min(domain) - min(priorData)
    
    if((max(priorData)>max(domain))) diffmax <- max(priorData) - max(domain)
    
    if((diffmax==0)&&(diffmin==0)) {
      ## Prior data range is inside of data range
      return(range(priorData))
    } else{
      
      ## Left tail
      if(diffmin!=0){
        coeff1 <- coeff*diffmin/(diffmax + diffmin)
        fx <- coef(fPI)
        CDF <- integralMoTBF(fPI)
        min1 <- as.function(CDF)(min(domain, priorData))
        int <- integrate(as.function(fPI), min(priorData, domain), min(domain))$value
        y1 <- int*coeff1
        if(-min1-y1>=0) sign = "+" else sign = ""
        CDF1 <- noquote(paste(CDF,sign,-min1-y1 , sep=""))
        CDF1 <- motbf(CDF1)
        q1 <- uniroot(as.function(CDF1), range(domain, priorData))$root
      }else {
        q1 <- min(priorData)
      }
      
      ## Right tail
      if(diffmax!=0){
        coeff2 <- coeff*diffmax/(diffmax + diffmin)
        CDF <- integralMoTBF(fPI)
        min2 <- as.function(CDF)(min(domain, priorData))
        int <- integrate(as.function(fPI), min(domain, priorData), max(domain))$value
        y2 <- (1-int)*coeff2
        if(-min2-(1-y2)>=0) sign = "+" else sign = ""
        CDF1 <- noquote(paste(CDF,sign,-min2-(1-y2), sep=""))
        CDF1 <- motbf(CDF1)
        q2 <- uniroot(as.function(CDF1), range(priorData, domain))$root
      } else {
        q2 <- max(priorData)
      }  

      ## New Prior domain
      newDomain <- c(q1,q2)
    }
    return(range(newDomain))
  } else{
    ## Parameter s is greater than the data size
    return(range(priorData))
  }
}

#' Get the Coefficients
#' 
#' Compute the coefficients for the linear opinion pool
#' 
#' @param fPI The fitted function to the prior data of class \code{"motbf"}.
#' @param rangeNewPriorData An array of length two with the new domain of the prior function.
#' @param fD The fitted function to the original data of class \code{"motbf"}.
#' @param data A \code{"numeric"} array which contains the values to fit.
#' @param domain A \code{"numeric"} array with the limits where defining the data function.
#' @param coeffversion A \code{"numeric"} value between \code{1--4} which contains the used version for computing the coefficients in the linear 
#' opinion pool to combine the prior function and the data function. By default \code{coeffversion = "4"} is used, so the combination
#' depends on the goodness of the model versus another random positive MoTBF model.
#' @details \code{coeffversion} can be:
#' \code{"1"} coef1 and coef2 are the sum of the probabilities of one of the function over the sum of all probabilities, respectively;
#' \code{"2"} coef1 and coef2 are the solution of a linear optimization problem which tries to maximize the sum 1 for each row of probabilities;
#' \code{"3"} coef1 and coef2 are the difference of the log-likelihood of the evaluated model and a random uniform model over the sum of both differences, respectively;
#' \code{"4"} coef1 and coef2 are the difference of the log-likelihood of the evaluated model and a ramdom positive MoTBF model over the sum of both differences, respectively.
#' @return A \code{"numeric"} value of length 2 giving the coefficients which are the weigth of the two function to combine.
#' @seealso \link{learnMoTBFpriorInformation}
#' @importFrom lpSolve lp
#' @export
#' @examples
#' 
#' ## Data
#' X <- rnorm(15)
#' 
#' ## Prior Data
#' priordata <- rnorm(5000)
#' 
#' ## Learning
#' confident <- 5
#' type <- "MOP"
#' f <- learnMoTBFpriorInformation(priorData = priordata, data = X, s = confident,
#' POTENTIAL_TYPE = type)
#' attributes(f)
#'  
#' ## Coefficients: linear opinion pool
#' getCoefficients(fPI = f$priorFunction, rangeNewPriorData = f$domain, fD = f$dataFunction, 
#' data = X, domain = range(X), coeffversion = 4)
#' 
#' getCoefficients(fPI = f$priorFunction, rangeNewPriorData = f$domain, fD = f$dataFunction, 
#' data = X, domain = range(X), coeffversion = 1)
#' 
#' getCoefficients(fPI = f$priorFunction, rangeNewPriorData = f$domain, fD = f$dataFunction, 
#' data = X, domain = range(X), coeffversion = 3)
#'
#' getCoefficients(fPI = f$priorFunction, rangeNewPriorData = f$domain, fD = f$dataFunction, 
#' data = X, domain = range(X), coeffversion = 2)
#' 
getCoefficients=function(fPI, rangeNewPriorData, fD, data, domain, coeffversion)
{
  ## Density values for both functions
  pfPI <- as.function(fPI)(data)
  pfD <- as.function(fD)(data)
  
  switch(coeffversion,   

## coeffversion <- 1     
{coef1 <- sum(pfPI)/(sum(pfPI,pfD)); coef2 <- sum(pfD)/(sum(pfPI,pfD))},

## coeffversion <- 2    
{XX <- matrix(c(pfPI, pfD), ncol=2); maxL=c()
 for(i in 1:nrow(XX)){
   maxLL <- lp("max", objective.in=XX[i,], const.mat=matrix(c(1,1), nrow=1), const.rhs=c(1), const.dir=c("="))$solution
   maxL <- matrix(rbind(maxL,maxLL), ncol=2)
 }
 coef1 <- mean(maxL[,1])
 coef2 <- 1-coef1},

## coeffversion <- 3    
{probUnif <- dunif(data, min(data, rangeNewPriorData), max(data, rangeNewPriorData));
 a <- sum(log(pfPI))-sum(log(probUnif)); if(a<0) a=0;
 b <- sum(log(pfD))-sum(log(probUnif)); if(b<0) b=0;
 coef1 <- a/(a+b); coef2 <- b/(a+b)},

## coeffversion <- 4
{NonNormalisedRandomMoTBF <- getNonNormalisedRandomMoTBF(length(coef(fPI))-1)   
 a <- sum(log(pfPI))-UpperBoundLogLikelihood(NonNormalisedRandomMoTBF,data, min(domain, rangeNewPriorData),max(domain, rangeNewPriorData)); if(a<0) a=0;
 b <- sum(log(pfD))- UpperBoundLogLikelihood(NonNormalisedRandomMoTBF,data, min(domain, rangeNewPriorData),max(domain, rangeNewPriorData)); if(b<0) b=0;
 coef1 <- a/(a+b); coef2 <- b/(a+b)}
)

return(c(coef1, coef2))
}

#' Ramdom MoTBF
#' 
#' Get an MoTBF function ensuring positives values without being
#' a normalized function.
#' 
#' @param degree A \code{"numeric"} value containing the degree of the random function.
#' @param POTENTIAL_TYPE A \code{"character"} string specifying the posibles potential
#' types, must be one of \code{"MOP"} or \code{"MTE"}. 
#' @return A \code{"numeric"} vector of length 2 giving the coefficients.
#' @export
#' @examples
#' 
#' getNonNormalisedRandomMoTBF(8, POTENTIAL_TYPE = "MOP")
#' getNonNormalisedRandomMoTBF(11, POTENTIAL_TYPE = "MTE")
#' 
getNonNormalisedRandomMoTBF <- function(degree, POTENTIAL_TYPE="MOP") 
{
  if(degree%%2!=0) degree <- degree-1
  parameters=c()
  if(POTENTIAL_TYPE=="MOP"){
    for(i in 1:((degree/2)+1)) parameters <- c(parameters,0.5,0)
    pol <- asMOPString(parameters)
    pol <- motbf(pol)
  } else {
    parameters <- rep(0.5, degree/2+1)
    pol <- asMTEString(parameters)
    pol <- motbf(pol)
  }
  return(pol)
}

#' Upper Bound Loglikelihood
#' 
#' Computes an upper bound of the expected loglikelihood of a dataset given a randomly 
#' generated MoTBF density.
#' 
#' @param f A function to evaluate of class \code{"character"}, \code{"motbf"} or others.
#' @param data A \code{"numeric"} array which contains the values to evaluate.
#' @param min A \code{"numeric"} value giving the lower limit of the function.
#' @param max A \code{"numeric"} value giving the highter limit of the function.
#' @return A \code{"numeric"} value which is the log-likelihood of the evaluated ramdom function.
#' @seealso \link{getNonNormalisedRandomMoTBF}
#' @export
#' @examples
#' 
#' data <- rnorm(20)
#' f <- getNonNormalisedRandomMoTBF(degree = 8, POTENTIAL_TYPE = "MOP")
#' UpperBoundLogLikelihood(f, data, min = -2.5, max = 3.2)
#' 
#' data <- rexp(20)
#' f <- getNonNormalisedRandomMoTBF(degree = 8, POTENTIAL_TYPE = "MTE")
#' UpperBoundLogLikelihood(f, data, min = 0, max = 5)
#' 
UpperBoundLogLikelihood <- function(f,data, min, max){
  n <- length(data)
  k <- integrate(as.function(f),min,max)$value
  prob <- as.function(f)(data)
  return(sum(log(prob))-n*log(k))
}

#' Multivariate Normal Sample
#' 
#' Generate a multivariate normal data vector taking into account the real data and the 
#' relationships with other variables in the dataset.
#' 
#' @param n A \code{"numeric"} value which is the size of the prior data to generate.
#' @param dataParents A data set of class \code{"data.frame"} giving the data of the set of coditional parent variables.
#' @param dataChild A \code{"numeric"} vector containing the original data of the child variable.
#' @return A \code{"numeric"} vector giving the prior data values.
#' @seealso \link{generateNormalPriorData}
#' @export
#' @examples
#' 
#' ## Data
#' data(ecoli)
#' data <- ecoli[,-c(1,9)] ## remove sequece.name and class
#' 
#' ## DAG
#' dag <- LearningHC(data)
#' plot(dag)
#' getChildParentsFromGraph(dag)
#' 
#' ## 1. Random sample
#' parents <- "mcg"
#' child <- "alm1"
#' n <- 1000
#' rnormMultiv(n, dataParents = data.frame(data[,parents]), dataChild = data[,child])
#' 
#' ## 2. Random sample
#' parents <- "alm1"
#' child <- "aac"
#' n <- 256
#' rnormMultiv(n, dataParents = data.frame(data[,parents]), dataChild = data[,child])
#' 
rnormMultiv <- function(n, dataParents, dataChild)
{
  ma <- as.data.frame(cbind(dataParents, dataChild))
  mu <- colMeans(ma); covar <- cov(ma)
  
  ## Covariate matrix
  cc <- covar[ncol(ma),c(1:(ncol(ma)-1))] 
  
  ## Computing the alphas for the linear regression
  ro <- sapply(1:(ncol(ma)-1), function(i) cc[i]/diag(covar)[i])
  alpha <- c(); for(i in 1:length(ro)) alpha <- paste(alpha,"-ro[", i, "]*mu[",i, "]",sep="")
  
  ## Computing the means for the linear regression
  muZ <- c(); for(i in 1:length(ro)) muZ <- paste(muZ,"+ro[",i,"]*dataParents[,",i,"]", sep="")
  muZ <- paste("mu[", ncol(ma), "]", alpha, muZ, sep="")
  
  Z <- rnorm(n,eval(parse(text=muZ)), sd(dataChild))
  return(Z) 
}

#' Prior Data
#' 
#' Generate a prior dataset taking in to account the relationships
#' between varibles inside a given network.
#' 
#' @param graph A network of the class \code{"bn"}, \code{"graphNEL"} or \code{"network"}.
#' @param data A datase of class \code{"data.frame"} containing the continuous variables of the dataset.
#' @param size A \code{"numeric"} value indicating the number of records to generate for each variable in the dataset.
#' @param means A \code{"numeric"} vector with the averiage of each variable. The names of the vector must be the name
#' of the variables of which the information is given a priori by the expert.
#' @param deviations A \code{"numeric"} vector with the desviations of each variable. The names of the vector must be the name
#' of the variables of which the information is given a priori by the expert. By default it is \code{NULL} and the desviations
#' of the given data are taken.
#' @seealso \link{rnormMultiv}
#' @return A normal prior data set of class \code{"data.frame"}.
#' @export
#' @examples
#' 
#' ## Data
#' data(ecoli)
#' data <- ecoli[,-c(1,9)] ## remove sequece.name and class
#' X <- TrainingandTestData(data, percentage_test = 0.95)
#' Xtraining <- X$Training
#' Xtest <- X$Test
#' 
#' ## DAG
#' dag <- LearningHC(data)
#' plot(dag)
#' 
#' ## Means and desviations
#' colnames(data)
#' 
#' m <- sapply(data, mean)
#' m <- m[-which(is.na(m))]
#' names(m)
#' 
#' d <- sapply(data, sd)
#' d <- d[-which(is.na(m))]
#' names(d)
#' 
#' ## Prior Dataset
#' n <- 5600
#' priorData <- generateNormalPriorData(dag, data = Xtraining, size = n, means = m)
#' summary(priorData)
#' ncol(priorData)
#' nrow(priorData)
#' class(priorData)
#' 
generateNormalPriorData <- function(graph, data, size, means, deviations=NULL)
{
  if(is.null(deviations)){
    options(warn=-1)
    for(i in 1:length(means)){
      if(is.na(means[i])) deviations[i] <- NA
      else deviations[i] <- sd(data[,i])
    }
    #deviations <- sapply(data, sd)
  }
  columnsNames <- colnames(data); XX=c()
  for(i in 1:length(columnsNames)){ 
    if(is.na(means[i])||is.na(deviations[i])) X <- NA
    else X <- rnorm(size, means[i], deviations[i])
    XX <- matrix(cbind(XX,X), nrow=size)
  }
  XX <- as.data.frame(XX)
  colnames(XX) <- columnsNames
  noNormal <- colnames(XX[sapply(XX[1,], is.na)])
  
  childrenAndParents <- getChildParentsFromGraph(graph,columnsNames)
  YY <- c()
  for(i in 1:length(columnsNames)){
    if(length(childrenAndParents[[i]])==1){
      XX[,i] <- XX[,i]
    } else {
      child <- childrenAndParents[[i]][1]
      dataChild <- cbind(XX[,child])
      if(child%in%noNormal){
        XX[,i] <- cbind(XX[,child])
        next
      }
      parents <- childrenAndParents[[i]][2:length(childrenAndParents[[i]])]
      if(all(parents%in%noNormal)){
        XX[,i] <- cbind(XX[,child])
        next
      }
      if(any(parents%in%noNormal)){
        parents <- parents[-which(parents%in%noNormal)]
      }
      dataParents <- as.data.frame(XX[,parents])
      XX[,i] <- rnormMultiv(size,dataParents, dataChild) 
    }
    YY <- XX   
  }
  pos <- which(sapply(YY[1,], is.na))
  if(length(pos)!=0) YY <- YY[,-pos]
  return(as.data.frame(YY))
}

