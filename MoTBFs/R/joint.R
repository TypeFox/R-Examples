#' Learning Joint Functions 
#' 
#' Two functions for learning joint MoTBFs. The first one, \code{parametersJointMoTBF()},
#' gets the parameters by solving a quadratic optimization problem, minimizing
#' the error of the empirical cumulative joint values versus the estimated ones.
#' Finally it derives these parameters in order to obtain parameters of the joint density function.
#' The second one, \code{jointMoTBF()}, fixes the equation of the joint function using 
#' the previously learned parameters and converting this \code{"character"} string in an 
#' object of class \code{"jointmotbf"}.
#' 
#' @name jointmotbf.learning
#' @rdname jointmotbf.learning
#' @param X A dataset of class \code{"data.frame"}.
#' @param ranges A \code{"numeric"} matrix by columns with the ranges of the varibles where fitting the function.
#' @param dimensions A \code{"numeric"} vector with the number of parameters of each varible.
#' @param object A list with the output of the function \code{parametersJointMoTBF()}.
#' @return
#' \code{parametersJointMoTBF()} returns a list with the elements: \bold{Parameters}, which contains the computed 
#' coefficients of the resulting function; 
#' \bold{Dimension} which is a \code{"numeric"} vector containing the number 
#' of coefficients used for each variable; 
#' \bold{Range} contains a \code{"numeric"} matrix with the domain of each variable by columns;
#' \bold{Iterations} contains a number of iterations needed to solve the problem;
#' \bold{Time} contains the time that the functions spent to solve the problem.
#' 
#' \code{jointMoTBF()} retunrs an object of class \code{"jointmotbf"}; It is a list whose unique visible element
#' is the mathematical expression, furthermore it contains the other elements of the output of the
#' \code{parametersJointMoTBF()} function.
#' @examples
#' ## 1. EXAMPLE 
#' require(mnormt) ## Load the package to generate a multinormal dataset
#' ## Dataset
#' Mean <- 0; nVar <- 2; ro <- 0
#' varcov <- matrix(c(rep(c(1, rep(ro, nVar)), nVar-1),1), nrow=nVar)
#' means <- rep(Mean, nVar)
#' X <- rmnorm(100,means,varcov) 
#' data <- standardizeDataset(data.frame(X))
#' 
#' ## Joint learnings
#' dim <- c(2,3)
#' param <- parametersJointMoTBF(X = data, dimensions = dim)
#'
#' param$Parameters
#' length(param$Parameters) ## prod(dim)
#' param$Dimension
#' param$Range
#' 
#' P <- jointMoTBF(param)
#' P
#' attributes(P)
#' class(P)
#' 
#' ###############################################################################
#' ## MORE EXAMPLES ##############################################################
#' ###############################################################################
#' ## ## Dataset
#' ## Mean <- 1; nVar <- 3; ro <- 0.5
#' ## varcov <- matrix(c(rep(c(1, rep(ro, nVar)), nVar-1),1), nrow=nVar)
#' ## means <- rep(Mean, nVar)
#' ## X <- rmnorm(200,means,varcov) 
#' ## data <- standardizeDataset(data.frame(X))
#' ## 
#' ## ## Joint learnings
#' ## dim <- c(3,2,4,2)
#' ## param <- parametersJointMoTBF(X = data, dimensions = dim)
#' ## 
#' ## param$Parameters
#' ## length(param$Parameters) ## prod(dim)
#' ## param$Dimension
#' ## param$Range
#' ## param$Time
#' ## 
#' ## P <- jointMoTBF(param)
#' ## P
#' ## attributes(P)
#' ## class(P)
#' ###############################################################################
#' ###############################################################################

#' @rdname jointmotbf.learning
#' @export
parametersJointMoTBF=function(X, ranges=NULL, dimensions=NULL)
{
  tm <- Sys.time()
  if(is.null(ranges)) ranges <- sapply(X, range)
  if(is.null(dimensions)){
    Fx <- lapply(X, univMoTBF, "MOP")
    ll <- list(); for(i in 1:length(Fx)) ll[[length(ll)+1]] <- length(coef(Fx[[i]]))
  } else {
    ll <- list(); for(i in 1:length(dimensions)) ll[[length(ll)+1]] <- dimensions[i]
  }
  
  dim <- c(); pos <- which(c(length(ll[[1]]), length(ll[[2]]))==min(length(ll[[1]]), length(ll[[2]])))
  if(pos[1]==1) {pos1 <- 1; pos2 <- 2}
  if(pos[1]!=1) {pos1 <- 2; pos2 <- 1}
  for(i in 1:length(ll[[pos1]])) for(j in 1:length(ll[[pos2]])) dim <- c(dim,list(c(ll[[pos1]][i],ll[[pos2]][j])))
  if(pos1==2) for(i in 1:length(dim)) dim[[i]] <- dim[[i]][length(dim[[i]]):1]
  
  size <- prod(sapply(1:length(ll), function(i) length(ll[[i]])))
  si <- lapply(1:ncol(X),function(i) rep(ll[[i]][1], size))
  dim <- c(); for(i in 1:length(si)) dim <- cbind(dim,unlist(si[[i]]))
  dim <- lapply(1:nrow(dim), function(i) dim[i,])
  
  ## Create the grid: depending on the n.record and the n.variables
  npointsgrid <- 100
  if(nrow(X)>npointsgrid){
    num <- round((npointsgrid^2)^(1/ncol(X)))
    ran <- sapply(X, range)
    d <- unique(as.vector(sapply(1:ncol(ran), function(i) which(X[,i]%in%ran[,i]))))
    if((num-length(d)<=0)) x1 <- sort(X[d,1])
    else{
      x <- X[-d,]
      p <- sample(1:nrow(x), (num-length(d)))
      p <- c(p,d)
      x1 <- sort(X[p,1])
    }
    p <- sapply(1:length(x1), function(i) which(x1[i]==X[,1]))
    if(class(p)=="list") p <- unlist(unique(p))
    x <- X[p,]
    rownames(x) <- 1:nrow(x)
  }else{
    if(ncol(X)==2){
      x1 <- sort(X[,1])
      p <- sapply(1:length(x1), function(i) which(x1[i]==X[,1]))
      if(class(p)=="list") p <- unlist(unique(p))
      x <- X[p,]
    } else{
      num <- round((npointsgrid^2)^(1/ncol(X)))
      ran <- sapply(X, range)
      d <- unique(as.vector(sapply(1:ncol(ran), function(i) which(X[,i]%in%ran[,i]))))
      if((num-length(d)<=0)) x1 <- sort(X[d,1])
      else{
        x <- X[-d,]
        p <- sample(1:nrow(x), (num-length(d)))
        p <- c(p,d)
        x1 <- sort(X[p,1])
      }
      p <- sapply(1:length(x1), function(i) which(x1[i]==X[,1]))
      if(class(p)=="list") p <- unlist(unique(p))
      x <- X[p,]
      rownames(x) <- 1:nrow(x)
    }
  }
  x <- expand.grid(x)
  
  ## Cumulative densities
  y <- jointCDF(X,x)
  
  P <- c()
  for(s in 1:length(dim)){
    Xt <- ranges; pos <- coefExpJointCDF(dim[[s]])
    n <- nrow(x); nterms <- length(pos)
    
    x1 <- c()
    for(j in 1:ncol(X)){
      xx=rep(1, n)
      for(i in 1:nterms){
        xi <- cbind((x[,j]^pos[[i]][j]))
        xx <- cbind(xx,xi)
      }
      x1[[length(x1)+1]] <- xx
    }
    xx <-Reduce("*", x1)
    
    ## Dmat
    XX <- t(xx)%*%xx 
    
    ## dvec
    Xy <- t(xx)%*%y 
    
    ## Constrains
    xNew <- list()
    for( j in 1:ncol(Xt)){
      xnew <- seq(min(Xt[,j]), max(Xt[,j]), length=round(60000^(1/ncol(Xt))))
      if(xnew[length(xnew)]!=max(Xt[,j])) xnew <- c(xnew,max(Xt[,j]))
      xNew[[length(xNew)+1]] <- xnew
    }
    xN <- expand.grid(xNew)

    ma22 <- c()
    for(i in 1:ncol(X)){
      d <- sapply(1:nterms, function(j) pos[[j]][i]*(xN[,i]^(pos[[j]][i]-1)))
      tf <- sapply(1:length(pos), function(j) pos[[j]][i]==0)
      if(any(tf)) d[,which(tf)] <- 0
      d <- cbind(rep(0, nrow(xN)),d)
      ma22[[length(ma22)+1]] <- t(d)
    }
    ma2 <- Reduce("*", ma22)

    ma3 <- c()
    for(j in 1:ncol(Xt)){
      ma <- 0
      for(i in 1:nterms) ma <- rbind(ma,(max(Xt[,j])^pos[[i]][j]) - (min(Xt[,j])^pos[[i]][j]))
      ma3 <- cbind(ma3,ma)
    }
    for(i in 1:(ncol(X)-1)) ma3[,i+1] <- ma3[,i] * ma3[,i+1]
    ma33 <- cbind(ma3[,i+1])
    
    ## Amat
    AA <- cbind(ma33, ma2)
    
    ## bvec
    B <- c(1,rep(1.0E-5,ncol(ma2)))
    
    tr <- tryCatch(solve.QP(XX, Xy, AA, B, meq=1), error = function(e) NULL)
    finaltm <- Sys.time() - tm
    if(is.null(tr)==T){
      return(NULL)
    }else{    
      soluc <- tr 
      parameters <- soluc$solution
    }
    
    ## Parameters PDF
    a <- sapply(pos, prod)
    param <- parameters[-1]*a
    parameters <- param[which(param!=0)]
    
    P <- list(Parameters = parameters, 
              Dimension = dim[[s]], Range = Xt, 
              Iterations = tr$iterations[1], 
              Time=finaltm)
  }
  
  return(P)
}

#' @rdname jointmotbf.learning
#' @export
jointMoTBF <- function(object){
  parameters <- object$Parameters
  dimensions <- object$Dimension
  v <- letters[24:(24+length(dimensions)-1)]
  v[is.na(v)] <- letters[1:length(v[is.na(v)])]
  if(length(parameters)==1) {
    s <- paste(parameters[1], "+0", sep="")
    for(i in 1:length(v)) s <- paste(s, "*", v[i], sep="")
    P <- list(Function = s, Domain = object$Range)
    P <- jointmotbf(P)
    return(P)
  }
  pos <- coefExpJointCDF(dimensions)
  t <- list()
  for(i in 1:length(pos)) if(prod(pos[[i]])!=0) 
    t[[length(t)+1]] <- pos[[i]] else next
  pos <- t
  
  if(any(dimensions==1)) posvardim1 <- which(dimensions==1)
  
  if(parameters[1]==0) s <- c() else s <- parameters[1]
  for(i in 2:length(parameters)){
    if(is.null(s)){
      s <- paste(s, parameters[i],sep="")
    } else{
      s <- paste(s, ifelse(parameters[i]>=0, "+", ""), parameters[i],sep="")
      if(i != length(parameters)) for(p in 1:length(v)) s <- paste(s, ifelse((pos[[i]][p]-1)==0, "", paste("*",v[p], ifelse((pos[[i]][p]-1)==1, "", paste("^", pos[[i]][p]-1, sep="")), sep="")),sep="")
      if(i == length(parameters)){
        for(p in 1:length(v)) s <- paste(s, ifelse((pos[[i]][p]-1)==0, ifelse(any(p==posvardim1),
                                                                              paste("*",v[p],"^",pos[[i]][p]-1, sep=""),""), 
                                                   paste("*",v[p], ifelse((pos[[i]][p]-1)==1, "",
                                                                          paste("^", pos[[i]][p]-1, sep="")), sep="")),sep="")
      }
    }
  }
  P <- noquote(s)
  P <- list(Function = P, Domain = object$Range,
            Iterations = object$Iterations, Time = object$Time)
  P <- jointmotbf(P)
  return(P)
}


#' Degree Function
#'
#' Compute the degree for each term of a cumulative joint function.
#'
#' @param dimensions A \code{"numeric"} vector including the number of parameters of each variable.
#' @return A list with n element. Each element contains a \code{numeric} vector with the degree for 
#' each variable and each term of the cumulative joint function.
#' @export
#' @examples
#'
#' ## Dimension of the joint PDF of 2 variables
#' dim <- c(4,5) 
#' ## Potentials of each term of the CDF
#' c <- coefExpJointCDF(dim)
#' length(c) + 1 ## prod(dim+1) ## plus 1 because of the constant coefficient
#'
#' ## Dimension of the joint density function of 2 variables
#' dim <- c(5,5,3)
#' ## Potentials of the cumulative function
#' coefExpJointCDF(dim)
#'
coefExpJointCDF <- function(dimensions)
{
  t <- lapply(1:length(dimensions), function(i) 0:dimensions[i])
  l <- sapply(1:length(t), function(i) length(t[[i]]))
  mm <- c()
  for(i in 1:length(t)){
    m <- c()
    for(j in 1:length(t[[i]])) m <- c(m,rep(t[[i]][j], prod(l[-c(1:i)])))
    if(i!=1) m <- rep(m, l[i-1])
    mm <- cbind(mm, m)
  }
  mm <- mm[-1,]
  colnames(mm) <- NULL
  if(is.matrix(mm)) pos <- lapply(1:nrow(mm), function(i) mm[i,])
  else pos <-list(as.numeric(mm))
  return(pos)
}


#' Cumulative Joint Distribution
#' 
#' Functions to compute the multivariate cumulative distribution. It uses a grid data
#' over the ranges of variables in the dataset to calculate the cumulative values.
#' 
#' @name jointCDF
#' @rdname jointCDF
#' @param data The dataset as data.frame class.
#' @param grid a data frame with the select rows to include in the grid
#' @param nrows Number of rows.
#' @details \code{jointCDF()} is the main function which uses the internal
#' function \code{posGrid()}.
#' @examples
#' ## 1. EXAMPLE
#' ## Dataset with 2 variables
#' X <- data.frame(rnorm(100), rnorm(100))
#' 
#' ## Grid dataset
#' gridX <- sapply(1:ncol(X), function(i) seq(min(X[,i]), max(X[,i]), length=10))
#' gridPoints <- expand.grid(gridX) ## expand gridX
#' ncol(gridPoints)
#' nrow(gridPoints)
#' 
#' ## Joint cumulative values
#' jointCDF(data = X, grid = gridPoints)
#' 
#' ## 2. EXAMPLE
#' ## Dataset with 3 variables
#' X <- data.frame(rnorm(100), rexp(100), rchisq(100, df=2))
#' 
#' ## Grid dataset
#' gridX <- sapply(1:ncol(X), function(i) seq(min(X[,i]), max(X[,i]), length=10))
#' gridPoints <- expand.grid(gridX) ## expand gridX
#' ncol(gridPoints)
#' nrow(gridPoints)
#' 
#' ## Joint cumulative values
#' jointCDF(data = X, grid = gridPoints)
#' 

#' @rdname jointCDF
#' @export
jointCDF <- function(data, grid) unlist(lapply(1:nrow(grid), posGrid, data, grid))

#' @rdname jointCDF
#' @export
posGrid <- function(nrows, data, grid)
{
  p <- lapply(1:ncol(grid), function(i) which(data[,i]<=grid[nrows,i]))
  pp <- c()
  for(i in 1:length(p)){
    pp <- c(pp,p[[i]])
    if(i!=1) pp <- pp[duplicated(pp, last=T)]
  }
  return(length(pp)/nrow(data))
}


#' Coerce a \code{"jointmotbf"} Object to a Function
#'
#' Takes a \code{"jointmotbf"} object and contructs an \R function to evaluate it at multidimensional points.
#' 
#' @param x An object of class \code{"joinmotbf"}.
#' @param \dots Further arguments to be passed to or from the method. Not necessary for this method.
#' @details This is an \code{S3} method for the generic function \link{as.function}.
#' @return It returns a function to evaluate an object of class \code{"jointmotbf"}.
#' @seealso \link{parametersJointMoTBF} and \link{jointMoTBF}
#' @export
#' @examples
#' ## 1.EXAMPLE
#' ## Dataset
#' data <- data.frame(X = rnorm(100), Y = rexp(100))
#' 
#' ## Joint function
#' dim <- c(3,2)
#' param <- parametersJointMoTBF(data, dimensions = dim)
#' P <- jointMoTBF(param)
#' density <- as.function(P)(data[,1], data[,2])
#' density
#' 
#' ## Log-likelihood
#' sum(log(density))
#' 
#' #############################################################################
#' ## MORE EXAMPLES ############################################################
#' #############################################################################
#' ## ## Dataset
#' ## data <- data.frame(X = rnorm(100), Y = rexp(100), Z = rnorm(100))
#' ##
#' ## ## Joint function
#' ## dim <- c(2,3,4)
#' ## param <- parametersJointMoTBF(data, dimensions = dim)
#' ## P <- jointMoTBF(param)
#' ## density <- as.function(P)(data[,1], data[,2], data[,3])
#' ## density
#' ##
#' ## ## Log-likelihood
#' ## sum(log(density))
#' #############################################################################
#' #############################################################################
#' 
as.function.jointmotbf <- function(x, ...)
{
  P <- x[[1]]
  v <- nVariables(P)
  f <- function(a,b,c,d,e,f,g,h,i,j,k,l,m,n,o,p,q,r,
                s,t,u,v,w,x,y,z) eval(parse(text=P))
  formals(f) <- formals(f)[1:length(v)]
  names(formals(f)) <- v
  return(f)  
}

#' Extract Coefficients  of a \code{"jointmotbf"} Object
#' 
#' Extracts the parameters of the learned multivariate mixtures of truncated basis
#' functions. 
#' 
#' @param object An MoTBF function.
#' @param \dots Other arguments, unnecessary for this function.
#' @return A \code{"numeric"} vector with the parameters of the function.
#' @seealso \link{parametersJointMoTBF} and \link{jointMoTBF}
#' @export
#' @examples
#' require(mnormt) ## Load the package to generate a multinormal dataset
#' ## Dataset
#' Mean <- 0; nVar <- 2; ro <- 0
#' varcov <- matrix(c(rep(c(1, rep(ro, nVar)), nVar-1),1), nrow=nVar)
#' means <- rep(Mean, nVar)
#' X <- rmnorm(100,means,varcov) 
#' data <- standardizeDataset(data.frame(X))
#' 
#' ## Joint function
#' dim <-c(2,4)
#' param <- parametersJointMoTBF(data, dimensions = dim)
#' P <- jointMoTBF(param)
#' P$Time
#' 
#' ## Coefficients
#' coef(P)
#'
#' #############################################################################
#' ## MORE EXAMPLES ############################################################
#' #############################################################################
#' ## ## Dataset
#' ## Mean <- 0; nVar <- 3; ro <- 0
#' ## varcov <- matrix(c(rep(c(1, rep(ro, nVar)), nVar-1),1), nrow=nVar)
#' ## means <- rep(Mean, nVar)
#' ## X <- rmnorm(100,means,varcov) 
#' ## data <- standardizeDataset(data.frame(X))
#' ## 
#' ## ## Joint function
#' ## dim <-c(2,4,3)
#' ## param <- parametersJointMoTBF(data, dimensions = dim)
#' ## P <- jointMoTBF(param)
#' ## P$Time
#' ##
#' ## ## Coefficients
#' ## coef(P)
#' ##############################################################################
#' ##############################################################################
coef.jointmotbf <- function(object, ...) coeffMOP(object)


#' Number of Variables in a Joint Function
#' 
#' Compute the number of variables which are in a 'jointmotbf' object.
#' 
#' @param P An \code{"motbf"} object or a \code{"jointmotbf"} object.
#' @return A \code{"character"} vector with the names of variables of the function.
#' @export
#' @examples
#' # 1. EXAMPLE
#' ## Dataset
#' require(mnormt) ## Load the package to generate a multinormal dataset
#' ## Dataset
#' Mean <- 0; nVar <- 2; ro <- 0
#' varcov <- matrix(c(rep(c(1, rep(ro, nVar)), nVar-1),1), nrow=nVar)
#' means <- rep(Mean, nVar)
#' X <- rmnorm(100,means,varcov) 
#' data <- standardizeDataset(data.frame(X))
#' 
#' ## Joint function
#' dim <-c(3,2)
#' param <- parametersJointMoTBF(data, dimensions = dim)
#' P <- jointMoTBF(param)
#' P
#' 
#' ## Variables
#' nVariables(P)
#' 
#' ##############################################################################
#' ## MORE EXAMPLES #############################################################
#' ##############################################################################
#' ## ## Dataset
#' ## Mean <- 0; nVar <- 3; ro <- 0
#' ## varcov <- matrix(c(rep(c(1, rep(ro, nVar)), nVar-1),1), nrow=nVar)
#' ## means <- rep(Mean, nVar)
#' ## X <- rmnorm(100,means,varcov) 
#' ## data <- standardizeDataset(data.frame(X))
#' ##
#' ## ## Joint function
#' ## dim <- c(2,1,3)
#' ## param <- parametersJointMoTBF(data, dimensions = dim)
#' ## P <- jointMoTBF(param)
#' ##
#' ## ## Variables
#' ## nVariables(P)
#' #############################################################################
#' #############################################################################
nVariables=function(P)
{
  Letters <- letters[c(24:26,1:23)]
  string <- as.character(P)
  t <- strsplit(string, split="-", fixed = T)[[1]]
  t2 <- c()
  for(i in 1:length(t)){
    t1 <- strsplit(t[i], split="+", fixed = T, perl = FALSE, useBytes = FALSE)[[1]]
    t2 <- c(t2,t1)
  }
  t <- unlist(sapply(1:length(t2), function(i) strsplit(t2[i], split="*", fixed = T, perl = FALSE, useBytes = FALSE)[[1]]))
  
  options(warn=-1)
  variables <- t[is.na(as.numeric(t))]
  t <- unlist(strsplit(variables, split="^", fixed = T, perl = FALSE, useBytes = FALSE))
  t <- unlist(strsplit(t, split="(", fixed = T, perl = FALSE, useBytes = FALSE))
  variables <- t[is.na(as.numeric(t))]
  variables <- unique(variables)
  variables <- Letters[Letters%in%variables]
  return(variables)
}

#' Dimension of Functions
#' 
#' Gets the dimension of \code{"motbf"} and \code{"jointmotbf"} functions.
#' 
#' @param P An object of class \code{"motbf"} and subclass 'mop' or \code{"jointmotbf"}.
#' @return Dimension of the function.
#' @seealso \link{univMoTBF} and \link{jointMoTBF}
#' @export
#' @examples
#' ## 1. EXAMPLE 
#' ## Data
#' X <- rnorm(2000)
#' 
#' ## Univariate function
#' subclass <- "MOP"
#' f <- univMoTBF(X, POTENTIAL_TYPE = subclass)
#' dimensionFunction(f)
#' 
#' ## 2. EXAMPLE 
#' ## Dataset with 2 variables
#' X <- data.frame(rnorm(100), rnorm(100))
#' 
#' ## Joint function
#' dim <- c(2,3)
#' param <- parametersJointMoTBF(X, dimensions = dim)
#' P <- jointMoTBF(param)
#' 
#' ## Dimension of the joint function
#' dimensionFunction(P)
#'
dimensionFunction <- function(P){
  string <- as.character(P)
  nVar <- nVariables(P)
  
  if(length(nVar)==1){
    str <- strsplit(string, split="^", fixed=T)[[1]]
    return(as.numeric(str[length(str)])+1)
  }
  dimensions <- c()
  for(i in 1:length(nVar)){
    elements <- strsplit(string, split=NULL)[[1]]
    pos <- grep(nVar[i], elements)
    nele <- elements[(pos[length(pos)]+1):(pos[length(pos)]+3)]
    if(is.na(as.numeric(nele[3]))) nele <- nele[-3]
    if(length(which(nele%in%"^"))!=0){
      #nele[1]=="^")
      pos <- as.numeric(nele[c(2,3)])
      if(all(!is.na(pos))){
        num <- c()
        for(h in 2:length(nele)) num <- paste(num,nele[h], sep="")
        dimensions=c(dimensions, as.numeric(num))
      }else dimensions=c(dimensions, as.numeric(nele[2]))
    } else dimensions=c(dimensions, 1)
  }
  return(dimensions+1)
}

#' Integral Joint MoTBF
#' 
#' Integrate a \code{"jointmotbf"} object over an non defined domain. It is able to
#' get the integral of a joint function over a set of variables or over all
#' the variables of the function.
#' 
#' @param P A \code{"jointmotbf"} object.
#' @param var A \code{"character"} vector containing the name of the variables where integrating the joint function.
#' Instead the names it can contain the position of the variables.
#' By default it's \code{NULL} then all the variables are taken.
#' @return A multiintegral of a joint function of class \code{"jointmotbf"}.
#' @export
#' @examples
#' ## 1. EXAMPLE
#' ## Dataset with 2 variables
#' X <- data.frame(rnorm(100), rnorm(100))
#' 
#' ## Joint function
#' dim <- c(2,3)
#' param <- parametersJointMoTBF(X, dimensions = dim)
#' P <- jointMoTBF(param)
#' 
#' ## Integral
#' integralJointMoTBF(P)
#' integralJointMoTBF(P, var="x")
#' integralJointMoTBF(P, var="y")
#' 
#' ##############################################################################
#' ## MORE EXAMPLES #############################################################
#' ##############################################################################
#' ## ## Dataset with 3 variables
#' ## X <- data.frame(rnorm(50), rnorm(50), rnorm(50))
#' ##
#' ## ## Joint function
#' ## dim <- c(2,1,3)
#' ## param <- parametersJointMoTBF(X, dimensions = dim)
#' ## P <- jointMoTBF(param)
#' ## 
#' ## ## Integral
#' ## integralJointMoTBF(P)
#' ## integralJointMoTBF(P, var="x")
#' ## integralJointMoTBF(P, var="y")
#' ## integralJointMoTBF(P, var="z")
#' ## integralJointMoTBF(P, var=c("x","z"))
#' ##############################################################################
#' ##############################################################################
integralJointMoTBF <- function(P, var=NULL)
{
  Letters <- c(letters[24:26], letters[1:23])
  if(is.numeric(var)) var <- Letters[var]
  if(is.null(var)) var <- nVariables(P)
  nVar <- nVariables(P)
  for(h in 1:length(var)){
    string <- as.character(P)
    parameters <- coef(P)
    f1 <- substr(string, 1, 1)
    t <- strsplit(string, split="-", fixed = T)[[1]]
    for(i in 1:length(t)) t[i] <- paste("-", t[i], sep="")##le volvemos a añadir el simbolo negativo
    if(f1!=substr(t[1], 1, 1)) t[1] <- substr(t[1], 2, nchar(t[1]))
    
    t2 <- c()
    for(i in 1:(length(t))){
      t1 <- strsplit(t[i], split="+", fixed = T, perl = FALSE, useBytes = FALSE)[[1]]
      t2 <- c(t2,t1)
    }
    
    pos <- grep("e", t2)
    if(length(pos)!=0){
      for(i in pos){
        t2[i] <- paste(t2[i],t2[i+1], sep="")
        t2[i+1] <- NA
      }
    }
    str <- t2[!is.na(t2)]
    
    s <- strsplit(str, split=as.character(parameters), fixed=T)
    for (i in 1:length(s)){
      if(is.na(s[[i]][2])) s[[i]][2]= paste("*", var[h], sep="")
      else{
        t <- strsplit(s[[i]][2], split=NULL)[[1]]
        if(!any(t%in%var[h])){
          pVar <- which(nVar%in%var[h])
          pVars <- which(nVar%in%t)
          if(any(pVar<pVars)){
            d <- which(nVar[pVars[which(pVar<pVars)[1]]]==t)
            t <- append(t, paste("*", var[h], sep=""), d-2)
            f <- c()
            for(j in 1:length(t))  f <- paste(f, t[j], sep="")
            s[[i]][2] <- f
          }else{
            s[[i]][2] <- paste(s[[i]][2], "*", var[h], sep="")
          }
        }else{
          p <- which(t%in%var[h])
          if(is.na(as.numeric(t[p+2]))){
            t <- append(t, paste("^" , 2, sep=""),after=p)
            parameters[i] <- parameters[i]/2
          }else{
            t[p+2]=as.numeric(t[p+2])+1
            parameters[i] <- parameters[i]/as.numeric(t[p+2])
          } 
          f <- c()
          for(j in 1:length(t))  f <- paste(f, t[j], sep="")
          s[[i]][2] <- f
        }
        
      }
    }
    if(any(parameters[2:length(parameters)]>=0)) parameters[parameters[2:length(parameters)]>=0]
    parameters[2:length(parameters)][parameters[2:length(parameters)]>=0] = 
      paste("+",parameters[2:length(parameters)][parameters[2:length(parameters)]>=0], sep="")
    
    s <- unlist(s)
    s[which(s%in%"")]=parameters
    str <- c()
    for(i in 1:length(s)) str <- paste(str, s[i], sep="")
    
    str <- noquote(str)
    if(length(nVariables(str))>1) P <- jointmotbf(str)
    else P <- motbf(str)
  }
  
  return(P)
}


#' Evaluate a Joint Function
#' 
#' Evaluates a \code{"jointmotbf"} object at a specific point.
#' 
#' @param P A \code{"jointmotbf"} object.
#' @param values A list with the name of the variables equal to the values to be evaluate.
#' @return If all the variables in the equation are evaluated then the return is a \code{"numeric"} value,
#' if not an \code{"motbf"} object or a \code{"jointmotbf"} object is returned.
#' @export
#' @examples
#' #' ## 1. EXAMPLE
#' ## Dataset with 2 variables
#' X <- data.frame(rnorm(100), rexp(100))
#' 
#' ## Joint function
#' dim <- c(3,3) # dim <- c(5,4)
#' param <- parametersJointMoTBF(X, dimensions = dim)
#' P <- jointMoTBF(param)
#' P
#' 
#' ## Evaluation
#' nVariables(P)
#' val <- list(x = -1.5, y = 3)
#' evalJointFunction(P, values = val)
#' val <- list(x = -1.5)
#' evalJointFunction(P, values = val)
#' val <- list(y = 3)
#' evalJointFunction(P, values = val)
#' 
#' ##############################################################################
#' ## MORE EXAMPLES #############################################################
#' ############################################################################## 
#' ## ## Dataset with 3 variables
#' ## X <- data.frame(rnorm(100), rexp(100), rnorm(100, 1))
#' ##
#' ## ## Joint function
#' ## dim <- c(2,1,3)
#' ## param <- parametersJointMoTBF(X, dimensions = dim)
#' ## P <- jointMoTBF(param)
#' ## P
#' ##
#' ## ## Evaluation
#' ## nVariables(P)
#' ## val <- list(x = 0.8, y = -2.1, z = 1.2)
#' ## evalJointFunction(P, values = val)
#' ## val <- list(x = 0.8, z = 1.2)
#' ## evalJointFunction(P, values = val)
#' ## val <- list(y = -2.1)
#' ## evalJointFunction(P, values = val)
#' ## val <- list(y = -2.1)
#' ## evalJointFunction(P, values = val)
#' ##############################################################################
#' ##############################################################################
evalJointFunction <- function(P, values)
{
  nVar <- nVariables(P)
  if(is.motbf(P)) return(as.numeric(as.function(P)(values[[1]])))
  string <- as.character(P)
  parameters <- coef(P)
  f1 <- substr(string, 1, 1)
  t <- strsplit(string, split="-", fixed = T)[[1]]
  for(i in 1:length(t)) t[i] <- paste("-", t[i], sep="")##le volvemos a añadir el simbolo negativo
  if(f1!=substr(t[1], 1, 1)) t[1] <- substr(t[1], 2, nchar(t[1]))
  
  t2 <- c()
  for(i in 1:(length(t))){
    t1 <- strsplit(t[i], split="+", fixed = T, perl = FALSE, useBytes = FALSE)[[1]]
    t2 <- c(t2,t1)
  }
  
  pos <- grep("e", t2)
  if(length(pos)!=0){
    for(i in pos){
      t2[i] <- paste(t2[i],t2[i+1], sep="")
      t2[i+1] <- NA
    }
  }
  str <- t2[!is.na(t2)]
  
  s <- strsplit(str, split=as.character(parameters), fixed=T)
  var <- names(values)
  for(i in 1:length(s)){
    for(j in 1:length(var)){
      if(length(grep(var[j],s[[i]]))==0){
        if(j==1){s[[i]][1]=parameters[i]; if(s[[i]][1]>0) s[[i]][1]=paste("+", s[[i]][1], sep="")}
      }else{
        if(j==1){s[[i]][1]=parameters[i]; if(s[[i]][1]>0) s[[i]][1]=paste("+", s[[i]][1], sep="")}
        if(!is.na(s[[i]][2])){
          s1=strsplit(s[[i]][2], split=NULL)[[1]]
          if(any(s1%in%var[j])){
            if(!is.na(s1[which(s1==var[j])+1])){
              if(s1[which(s1==var[j])+1]=="^"){
                exp=as.numeric(s1[which(s1==var[j])+2])
                s[[i]][1]=as.numeric(s[[i]][1])*(values[[j]]^exp)
                if(s[[i]][1]>0) s[[i]][1]=paste("+", s[[i]][1], sep="")
                s1[which(s1==var[j])+1]=""
                s1[which(s1==var[j])-1]=""
                s1[which(s1==var[j])+2]=""
                s1[s1%in%var[j]]=""
                s2=c(); for(h in 1:length(s1)) s2=paste(s2, s1[h], sep="")
                s[[i]][2]=s2
              }else{
                s[[i]][1]=as.numeric(s[[i]][1])*values[[j]]
                if(s[[i]][1]>0) s[[i]][1]=paste("+", s[[i]][1], sep="")
                s1[which(s1%in%var[j])-1]=""
                s1[which(s1%in%var[j])]=""
                s2=c(); for(h in 1:length(s1)) s2=paste(s2, s1[h], sep="")
                s[[i]][2]=s2
              }
            }else{
              s[[i]][1]=as.numeric(s[[i]][1])*values[[j]]
              if(s[[i]][1]>0) s[[i]][1]=paste("+", s[[i]][1], sep="")
              s1[which(s1%in%var[j])-1]=""
              s1[which(s1%in%var[j])]=""
              s2=c(); for(h in 1:length(s1)) s2=paste(s2, s1[h], sep="")
              s[[i]][2]=s2
            }
          } 
        }
      }
    }
  }
  
  if(is.na(s[[1]][2])) s[[1]][2]=""
  exp <- sapply(1:length(s), function(i) s[[i]][2])
  exp <- unique(exp)  
  s <- unlist(s)
  param <- c()
  for(i in 1:length(exp)) param <- c(param,sum(as.numeric(s[which(s==exp[i])-1])))
  sign <- ""
  if(length(param)>1) for(i in 2:length(param)) if(param[i]>0) sign <- c(sign,"+") else sign <- c(sign,"")
  str <- c()
  for(i in 1:length(exp)) str <- paste(str, sign[i], param[i], exp[i], sep="")
  
  if(all(names(values)==nVar)){
    cat("The method 'as.function()' can be used. \n")
    return(eval(parse(text=str)))
  } else{
    f <- noquote(str)
    l <- length(nVariables(f))
    if(l==1) f <- motbf(f)
    if(l>1) f <- jointmotbf(f)
    return(f)
  }  
}

#'Marginal Joint MoTBF
#'
#'Computes the marginal functions of a \code{"jointmotbf"}
#'object.
#'
#'@param P A joint function of class \code{"jointmotbf"}.
#'@param var The \code{"numeric"} position or the \code{"character"} name of the marginal variable.
#'@return The marginal of a \code{"jointmotbf"} function. The result is an object of class \code{"motbf"}.
#'@seealso \link{jointMoTBF} and \link{evalJointFunction}
#'@export
#'@examples
#' ## 1. EXAMPLE 
#' ## Dataset with 2 variables
#' X <- data.frame(rnorm(100), rnorm(100))
#' 
#' ## Joint function
#' dim <- c(4,3)
#' param <- parametersJointMoTBF(X, dimensions = dim)
#' P <- jointMoTBF(param)
#' P
#' 
#' ## Marginal
#' marginalJointMoTBF(P, var = "x")
#' marginalJointMoTBF(P, var = 2)
#' 
#' ##############################################################################
#' ## MORE EXAMPLES #############################################################
#' ##############################################################################
#' ## require(mnormt) ## Load the package to generate a multinormal dataset
#' ## ## Dataset with 5 variables
#' ## Mean <- 0; nVar <- 5; ro <- 0.3
#' ## varcov <- matrix(c(rep(c(1, rep(ro, nVar)), nVar-1),1), nrow=nVar)
#' ## means <- rep(Mean, nVar)
#' ## X <- rmnorm(100,means,varcov) 
#' ## data <- standardizeDataset(data.frame(X))
#' ##
#' ## ## Joint function
#' ## dim <- c(2,1,3,1,2)
#' ## param <- parametersJointMoTBF(data, dimensions = dim)
#' ## P <- jointMoTBF(param)
#' ## nVariables(P)
#' ##
#' ## ## Marginal
#' ## marginalJointMoTBF(P, var="x")
#' ## marginalJointMoTBF(P, var="y")
#' ## marginalJointMoTBF(P, var="z")
#' ## marginalJointMoTBF(P, var="a")
#' ## marginalJointMoTBF(P, var="b")
#' #############################################################################
#' #############################################################################
marginalJointMoTBF <- function(P, var)
{
  Letters <- c(letters[24:26], letters[1:23])
  if(is.numeric(var)) var <- Letters[var] else var <- var
  varN <- nVariables(P); noVar=varN[!varN%in%var]
  l <- length(noVar)
  for(j in 1:l){
    varN <- nVariables(P)
    posnoVar <- which(varN%in%noVar[j])
    posvar <- which(varN%in%var)
    
    min <- lapply(posnoVar, function(i) P$Domain[1,i])
    names(min) <- noVar[j]
    max <- lapply(posnoVar, function(i) P$Domain[2,i])
    names(max) <- noVar[j]
    
    iP <- integralJointMoTBF(P, noVar[j])
    minF <- evalJointFunction(iP, min)
    maxF <- evalJointFunction(iP, max)
    if(is.numeric(minF)&&is.numeric(maxF)){
      P <- list(Function=(maxF-minF), Subclass="mop", Domain= P$Domain[,which(!varN%in%noVar[j])])
      P <- motbf(P)
      return(P)
    }
    
    ## Min part
    nVar <- nVariables(minF)
    string <- as.character(minF)
    f1 <- substr(string, 1, 1)
    t <- strsplit(string, split="-", fixed = T)[[1]]
    for(i in 1:length(t)) t[i] <- paste("-", t[i], sep="")
    if(f1!=substr(t[1], 1, 1)) t[1] <- substr(t[1], 2, nchar(t[1]))
    if(t[1]=="-") t=t[-1]
    t2 <- c()
    for(i in 1:(length(t))){
      t1 <- strsplit(t[i], split="+", fixed = T, perl = FALSE, useBytes = FALSE)[[1]]
      t2 <- c(t2,t1)
    }
    pos <- grep("e", t2)
    if(length(pos)!=0){
      for(i in pos){
        t2[i] <- paste(t2[i],t2[i+1], sep="")
        t2[i+1] <- NA
      }
    }
    str <- t2[!is.na(t2)]
    coefMIN <- coef(minF)
    
    sMIN <- unlist(strsplit(str, split=coefMIN, fixed=T))
    varseq=(unique(sMIN))[-1]
    sMIN[which(sMIN=="")]=coefMIN
    paramMIN <-c()
    for(i in 1:length(varseq)){
      pos <- which(sMIN%in%varseq[i])
      paramMIN <- c(paramMIN, sum(as.numeric(sMIN[pos-1])))
      sMIN <- sMIN[-c(pos, pos-1)]
    }
    paramMIN <- c(sum(as.numeric(sMIN)), paramMIN)
    
    ## Max part
    nVar <- nVariables(maxF)
    string <- as.character(maxF)
    f1 <- substr(string, 1, 1)
    t <- strsplit(string, split="-", fixed = T)[[1]]
    for(i in 1:length(t)) t[i] <- paste("-", t[i], sep="")
    if(f1!=substr(t[1], 1, 1)) t[1] <- substr(t[1], 2, nchar(t[1]))
    if(t[1]=="-") t=t[-1]
    t2 <- c()
    for(i in 1:(length(t))){
      t1 <- strsplit(t[i], split="+", fixed = T, perl = FALSE, useBytes = FALSE)[[1]]
      t2 <- c(t2,t1)
    }
    pos <- grep("e", t2)
    if(length(pos)!=0){
      for(i in pos){
        t2[i] <- paste(t2[i],t2[i+1], sep="")
        t2[i+1] <- NA
      }
    }
    str <- t2[!is.na(t2)]
    coefMAX <- coef(maxF)
    sMAX <- unlist(strsplit(str, split=coefMAX, fixed=T))
    varseq=(unique(sMAX))[-1]
    sMAX[which(sMAX=="")]=coefMAX
    paramMAX <-c()
    for(i in 1:length(varseq)){
      pos <- which(sMAX%in%varseq[i])
      paramMAX <- c(paramMAX, sum(as.numeric(sMAX[pos-1])))
      sMAX <- sMAX[-c(pos, pos-1)]
    }
    paramMAX <- c(sum(as.numeric(sMAX)), paramMAX)
    
    ## New parameters
    newparam <- paramMAX - paramMIN
    
    sign=""
    for(v in 2:length(newparam)) if(newparam[v]>0) sign=c(sign,"+") else sign=c(sign,"")
    st <- newparam[1]
    for(v in 2:length(newparam)) st <- paste(st, sign[v],newparam[v], varseq[v-1], sep="")
    if(length(nVariables(st))>1){
      P <- list(Function=noquote(st), Domain= P$Domain[,which(!varN%in%noVar[j])])
      P <- jointmotbf(P)
    }else{
      P <- list(Function=noquote(st), Subclass="mop", Domain= P$Domain[,which(!varN%in%noVar[j])])
      P <- motbf(P)
      if(length(unique(coeffPol(P)))==1){
        st <- paste(sum(coef(P)), ifelse(unique(coeffPol(P))==0, "",paste("*", var, ifelse(unique(coeffPol(P))==1, "", paste("^", unique(coeffPol(P)), sep="")), sep="")), sep="")
        P <- list(Function=noquote(st), Subclass="mop", Domain= P$Domain)
        P <- motbf(P)
      }
    }
  }
  return(P)
}

#' Bidimensional Plots for \code{'jointmotbf'} Objects
#' 
#' Draws the perpective and the contour plots for joint functions.
#' 
#' @param x An object of class \code{'jointmotbf'}.
#' @param type A \code{"character"} string: by default it is \code{"contour"}, the alternative one is \code{"perspective"}.
#' @param ranges A \code{"numeric"} matrix containing the domain of the variables by columns.
#' @param orientation A \code{"numeric"} vector indicating the perpective of the plot in degrees.
#' By default it is \code{(5,-30)}.
#' @param ticktype A \code{"character"} string: by default it is \code{"simple"} which draws just an arrow parallel to the axis 
#' to indicate direction of increase; \code{"detailed"} draws normal ticks.
#' @param data A dataset of class \code{"data.frame"}. Only two columns are allowed, it means two variables.
#' It is only necessary to draw the points over the main plot. By defailt it is \code{NULL}.
#' @param filled A logical argument; it is only used if \code{type = "contour"}
#' is active; by default is \code{TRUE}  so filled contours are plotted.
#' @param \dots Further arguments to be passed as for \link{plot}.
#' @method plot jointmotbf
#' @return A plot of the joint function.
#' @seealso \link{jointMoTBF}
#' @export
#' @examples
#' ## 1 .EXAMPLE
#' ## Dataset
#' X <- data.frame(rnorm(500), rnorm(500))
#' 
#' ## Joint function
#' dim <- c(3,3) ## dim <- c(5,5)
#' param1 <- parametersJointMoTBF(X, dimensions = dim)
#' P <- jointMoTBF(param1)
#' P
#' 
#' ## Plots
#' plot(P)
#' plot(P, type = "perspective", orientation = c(90,0))
#' 
#' #############################################################################
#' ## MORE EXAMPLES ############################################################
#' #############################################################################
#' ## ## Dataset
#' ## X <- data.frame(rnorm(200,2), rexp(200, 1))
#' ##
#' ## ## Joint function
#' ## dim <- c(4,5)
#' ## param2 <- parametersJointMoTBF(X, dimensions = dim)
#' ## P <- jointMoTBF(param2)
#' ## P
#' ##
#' ## ## Plots
#' ## plot(P)
#' ## plot(P, filled = FALSE, data = X)
#' ## plot(P, type = "perspective", orientation = c(10,180))
#' #############################################################################
#' #############################################################################
plot.jointmotbf <- function(x, type="contour", ranges=NULL, orientation=c(5,-30), 
                            data=NULL, filled=TRUE, ticktype="simple", ...)
{
  if(length(nVariables(x))!=2) return(cat("It is not possible plotting a joint function with more than 2 variables."))
  if(is.null(ranges)) ranges <- x$Domain
  #{ranges <- x$Domain; ranges[1,] <- ranges[1,]+2E-1; ranges[2,] <- ranges[2,]-2E-1}
  X <- seq(min(ranges[,1]), max(ranges[,1]), length = 20)
  Y <- seq(min(ranges[,2]), max(ranges[,2]), length = 25)
  Z <- outer(X, Y, as.function(x))
  Z[Z<=0]=1.0E-10
  
  if(type == "perspective"){    
    nrz <- nrow(Z)
    ncz <- ncol(Z)
    ## Create a function interpolating colors in the range of specified colors
    jet.colors <- colorRampPalette( c("blue", "green") )
    
    ## Generate the desired number of colors from this palette
    nbcol <- 100; color <- jet.colors(nbcol)
    
    ## Compute the z-value at the facet centres
    zfacet <- Z[-1, -1] + Z[-1, -ncz] + Z[-nrz, -1] + Z[-nrz, -ncz]
    
    ## Recode facet z-values into color indices
    facetcol <- cut(zfacet, nbcol)
    persp(X, Y, Z, col = color[facetcol], phi = orientation[1], theta = orientation[2], ticktype=ticktype) 
  }
  if(type == "contour"){
    if(!filled){
      if(!is.null(data)) plot(data, xlab = "", ylab = "", xlim=ranges[,1], ylim=ranges[,2])
      if(is.null(data)) plot(NULL, xlim=ranges[,1], ylim=ranges[,2], xlab="", ylab="")
      contour(X,Y,Z, col=terrain.colors(15), add=T)
      abline(h =round(ranges[,1])[1]:round(ranges[,1])[2] , v =round(ranges[,2])[1]:round(ranges[,2])[2],
             col = "gray", lty = 2, lwd = 0.1)   
    }
    if(filled){
      zlim <- range(Z, finite=TRUE)
      zlim[1] <- 0
      nlevels <- 20
      levels <- pretty(zlim, nlevels)
      nlevels <- length(levels)
      color <- colorRampPalette(c("#FFFFD9", "#EDF8B1", "#C7E9B4", "#7FCDBB", "#41B6C4", "#1D91C0", "#225EA8", "#253494", "#081D58"))
      filled.contour(X,Y,Z,nlevels=nlevels,levels=levels,
                     las=1,
                     col=color(nlevels),                           
                     #col=rainbow(nlevels, alpha=0.8),                           
                     plot.title = {title(xlab = "X", ylab = "Y")},
      )
      if(!is.null(data)){
        mar.orig <- par("mar")
        w <- (3 + mar.orig[2]) * par("csi") * 2.54
        layout(matrix(c(2, 1), ncol = 2), widths = c(1, lcm(w)))
        points(data)
        par(mfrow=c(1,1))
      }
    }
  }
}
