#' Learning Conditional Functions
#' 
#' Collection of functions for learning conditional MoTBFs,
#' computing the internal BIC, selecting the parents that get 
#' the best BIC value, and other internal functions needed to get the 
#' final conditionals.
#' 
#' @name conditionalmotbf.learning
#' @rdname conditionalmotbf.learning
#' 
#' @param data A dataset of class \code{"data.frame"}.
#' @param nameParents A \code{"character"} vector with the names of the parent variables.
#' @param nameChild The name of the child variable as a \code{"character"} string.
#' @param domainChild An \code{"numeric"} array with the range of the child variable.
#' @param domainParents A matrix of class \code{"matrix"} with the range of the parent 
#' variables or an \code{"numeric"} array if there is a only one parent.
#' @param numIntervals A \code{"numeric"} value indicating the maximum number of intervals in which we 
#' want to split the domain of the parent variables.
#' @param POTENTIAL_TYPE A \code{"character"} string specifying the posibles potential
#' types, must be one of \code{"MOP"} or \code{"MTE"}. You can specify just the initial letter.
#' @param maxParam A \code{"numeric"} value which indicate the maximum number of coefficients in the function. By default it is \code{NULL}; 
#' if not, the output is the function which gets the best BIC with at most this number of parameters.
#' @param s A \code{"numeric"} coefficient which fixes the confidence of the prior knowledge 
#' we are going to introduce. By default it is \code{NULL}, only we must modify it if we want 
#' to incorporate prior information to the fits.
#' @param priorData Prior dataset with values of the variables we have information apriori about.
#' This dataset must be of \code{"data.frame"} class.
#' @param conditionalfunction The output of the function \code{learn.tree.Intervals}.
#' @param mm One of the inputs and the output of the recursive function \code{"conditional"}.
#' @return The main function \code{conditionalMethod} returns a list with the name of the parents, 
#' the different intervals and  the fitted functions.
#' @details The main function \code{conditionalMethod} fits truncated functions for a variable which depends
#' on others variables. The domain of the parent variables is splitted
#' in different intervals and univariate functions are fitted in these 
#' ranges. The rest of the described functions are internal ones of the main function.
#' @seealso \link{printConditional}
#' @examples
#' ## Dataset
#' X <- rnorm(1000)
#' Y <- rbeta(1000, shape1 = abs(X)/2, shape2 = abs(X)/2)
#' Z <- rnorm(1000, mean = Y)
#' data <- data.frame(X = X, Y = Y, Z = Z)
#' 
#' ## Conditional Method
#' parents <- c("X","Y")
#' child <- "Z"
#' intervals <- 2
#'
#' potential <- "MTE"
#' fMTE <- conditionalMethod(data, nameParents = parents, nameChild = child, 
#' numIntervals = intervals, POTENTIAL_TYPE = potential)
#' printConditional(fMTE)
#' 
#' ##############################################################################
#' ## potential <- "MOP"
#' ## fMOP <- conditionalMethod(data, nameParents = parents, nameChild = child, 
#' ## numIntervals = intervals, POTENTIAL_TYPE = potential, maxParam = 15)
#' ## printConditional(fMOP)
#' ##############################################################################
#' 
#' ##############################################################################
#' ## Internal functions: Not needed to run #####################################
#' ##############################################################################
#' ## domainP <- range(data[,parents])
#' ## domainC <- range(data[, child])
#' ## t <- conditional(data, nameParents = parents, nameChild = child, 
#' ## domainParents = domainP, domainChild = domainC, numIntervals = intervals, 
#' ## mm = NULL, POTENTIAL_TYPE = potential)
#' ## printConditional(t)
#' ## selection <- select(data, nameParents = parents, nameChild = child, 
#' ## domainParents = domainP, domainChild = domainC, numIntervals = intervals, 
#' ## POTENTIAL_TYPE = potential)
#' ## parent1 <- selection$parent; parent1
#' ## domainParent1 <- range(data[,parent1])
#' ## treeParent1 <- learn.tree.Intervals(data, nameParents = parent1, 
#' ## nameChild = child, domainParents = domainParent1, domainChild = domainC, 
#' ## numIntervals = intervals, POTENTIAL_TYPE = potential)
#' ## BICscoreMoTBF(treeParent1, data, nameParents = parent1, nameChild = child)
#' ###############################################################################
#' ###############################################################################

#' @rdname conditionalmotbf.learning
#' @export
conditionalMethod <- function(data, nameParents, nameChild, numIntervals, POTENTIAL_TYPE, maxParam=NULL, s=NULL, priorData=NULL)
{
  if((POTENTIAL_TYPE=="MOP")||(POTENTIAL_TYPE=="MTE")){
    data <- newData(data, nameChild, nameParents)
    
    ## Domains
    if(is.character(data[,nameChild])) domainChild <- discreteVariablesStates(nameChild, data)[[1]]$states
    else domainChild <- range(data[,nameChild])
    
    if(is.numeric(data[,nameParents])) domainParents <- range(data[,nameParents])
    else{
      domainParents <- lapply(1:length(nameParents), function(i) if(is.numeric(data[,nameParents[i]])) range(data[,nameParents[i]]) else discreteVariablesStates(nameParents[i], data)[[1]]$states)
      names(domainParents) <- nameParents
    }
    
    ## Recursive process
    mm <- c()
    mm <- conditional(data, nameParents, nameChild, domainChild, domainParents, numIntervals, mm, POTENTIAL_TYPE, maxParam, s, priorData)
    return(mm)
  }else{
    return(cat("Unknown method, please use MOP or MTE"))
  }
}

#'@rdname conditionalmotbf.learning
#'@export
conditional <- function(data, nameParents, nameChild, domainChild, domainParents, numIntervals, mm, POTENTIAL_TYPE, maxParam=NULL, s=NULL, priorData=NULL)
{  
  ## select the parent who get the best BIC score when its domain is splitted
  f <- select(data, nameParents, nameChild, domainChild, domainParents, numIntervals, POTENTIAL_TYPE, maxParam, s, priorData)
  nameParents <- nameParents[which(nameParents!=f$parent)]
  
  for(i in 1:length(f$t)){
    m <- list(parent=f$parent, interval=f$t[[i]]$interval, Px=f$t[[i]]$Px)
    mm[[length(mm)+1]] <- m
    if(is.numeric(data[,f$parent])){
      dataInterval <- splitdata(data, f$parent, f$t[[i]]$interval[1], f$t[[i]]$interval[2])
      if(nrow(dataInterval)>20||length(nameParents)==0){
        if(length(nameParents)==0){
          mm <- mm
        } else {
          mm[[length(mm)]]$Px <- NULL
          
          ## An recursive process
          mm <- conditional(dataInterval, nameParents, nameChild,domainChild, domainParents, numIntervals, mm, POTENTIAL_TYPE, maxParam, s, priorData)
        }
      }
    } else{
      dataInterval <- subset(data,(data[,f$parent]==f$t[[i]]$interval))
      if(nrow(dataInterval)>20||(length(nameParents)==0)){
        if(length(nameParents)==0){
          mm <- mm
        } else {
          mm[[length(mm)]]$Px <- NULL
          mm <- conditional(dataInterval, nameParents, nameChild,domainChild, domainParents, numIntervals, mm, POTENTIAL_TYPE, maxParam, s, priorData)
        }
      }
    }
  }
  return(mm)
}

#'@rdname conditionalmotbf.learning
#'@export
select <- function(data, nameParents, nameChild, domainChild, domainParents, numIntervals, POTENTIAL_TYPE, maxParam=NULL, s=NULL, priorData=NULL)
{
  bestbic <- -10^10; bestvalues <- 0; Bic <- 0
  for(i in 1:length(nameParents)){
    t <- learn.tree.Intervals(data, nameParents[i], nameChild, domainParents, domainChild, numIntervals, POTENTIAL_TYPE, maxParam, s, priorData)
    Bic <- BICscoreMoTBF(t, data, nameParents[i], nameChild)
    values <- list(parent=nameParents[i], t=t)
    
    if(bestbic<Bic){
      bestbic <- Bic
      bestvalues <- values
    }
  }
  return(bestvalues)      
}

#'@rdname conditionalmotbf.learning
#'@export
learn.tree.Intervals <- function(data, nameParents, nameChild, domainParents, domainChild, numIntervals, POTENTIAL_TYPE, maxParam=NULL, s=NULL, priorData=NULL)
{
  X <- data[, nameParents]
  Y <- data[, nameChild]
  if(is.numeric(domainParents)) Xrange <- domainParents
  else Xrange <- domainParents[nameParents][[1]]
  
  NY <- length(Y)
  if(is.numeric(X)){
    priorD <- priorData[,nameChild]
    bestb <- min(X)-1
    B <- quantileIntervals(X, numIntervals)
    Pp1 <- list(); points <- c()
    for(i in 1:length(B)){
      Yf <- Y[which(X>bestb)]
      Xf <- X[which(X>bestb)]
      if(length(Y)==0) next
      if(is.character(Y)){
        cutPoints <- discreteVariablesStates(nameChild, data)[[1]]$states
        pD <- probDiscreteVariable(cutPoints, Yf)
        pos <- which(domainChild%in%cutPoints)
        coeff <- rep(0,length(domainChild))
        pD$coeff <- replace(coeff, pos, pD$coeff)
        names(pD$coeff) <- domainChild
        pD$sizeDataLeaf <- replace(coeff, pos, pD$sizeDataLeaf)
        prob <- list(pD)
        bestBIC <-  getBICDiscreteBN(prob)
      } else{
        if(is.null(priorData)) P <- univMoTBF(Yf, POTENTIAL_TYPE, domainChild, maxParam=maxParam)
        else P <- learnMoTBFpriorInformation(priorD, Yf, s, POTENTIAL_TYPE, domainChild, maxParam=maxParam)$posteriorFunction 
        bestBIC <- BICMoTBF(P, Yf)
      }
      b <- B[i]
      Xl <- Xf[Xf<=b]; Xl <- sort(Xl)
      Yl <- Yf[which(Xf<=b)]
      if(length(Yl)<=5) break
      if(length(Yl)==0) next
      if(max(Yl)==min(Yl)) next
      
      ## discrete child variable
      if(is.character(Yl)){
        cutPoints <- unique(Yl)
        pD <- probDiscreteVariable(cutPoints, Yl)
        pos <- which(domainChild%in%cutPoints)
        coeff <- rep(0,length(domainChild))
        pD$coeff <- replace(coeff, pos, pD$coeff)
        names(pD$coeff) <- domainChild
        pD$sizeDataLeaf <- replace(coeff, pos, pD$sizeDataLeaf)
        Px1 <- pD
      } else{
        if(is.null(priorD)) Px1 <- univMoTBF(Yl, POTENTIAL_TYPE, domainChild, maxParam=maxParam) 
        else Px1 <- learnMoTBFpriorInformation(priorD, Yl, s, POTENTIAL_TYPE, domainChild, maxParam=maxParam)$posteriorFunction 
      }
      
      Xr <- Xf[Xf>b]; Xr <- sort(Xr)
      Y22 <- Yf[which(Xf>b)]; Y22 <- sort(Y22)
      if(length(Y22)==0) next
      if(max(Y22)==min(Y22)) next
      
      if(is.character(Y22)){
        cutPoints <- unique(Y22)
        pD <- probDiscreteVariable(cutPoints, Y22)
        pos <- which(domainChild%in%cutPoints)
        coeff <- rep(0,length(domainChild))
        pD$coeff <- replace(coeff, pos, pD$coeff)
        names(pD$coeff) <- domainChild
        pD$sizeDataLeaf <- replace(coeff, pos, pD$sizeDataLeaf)
        Px2 <- pD
      } else{
        if(is.null(priorD)) Px2 <- univMoTBF(Y22, POTENTIAL_TYPE, domainChild, maxParam=maxParam)
        else Px2 <- learnMoTBFpriorInformation(priorD, Y22, s, POTENTIAL_TYPE, domainChild, maxParam=maxParam)$posteriorFunction 
      }
      
      if(is.character(Y)){
        DiscreteBN <-  list(Px1, Px2)
        BICT <-  getBICDiscreteBN(DiscreteBN)
      } else {
        PX <- list(Px1,Px2)
        multiY <- list(Yl, Y22)
        BICT <- BICMultiFunctions(PX, multiY)
      }
      
      if(is.na(BICT)) next ###quitar
      
      if(BICT>bestBIC){
        Pp1[[length(Pp1)+1]] <- Px1
        Pp2 <- Px2
        bestb <- b
        points <- c(points,bestb)
      }
    }
    
    if(is.null(points)){
      Y <- Y[which(X==X)]
      if(is.character(Y)){
        cutPoints <- unique(Y)
        pD <- probDiscreteVariable(cutPoints, Y)
        pos <- which(domainChild%in%cutPoints)
        coeff <- rep(0,length(domainChild))
        pD$coeff <- replace(coeff, pos, pD$coeff)
        names(pD$coeff) <- domainChild
        pD$sizeDataLeaf <- replace(coeff, pos, pD$sizeDataLeaf)
        prob <- list(pD)
        v <- c()
        values <- list(Px=prob, interval=c(Xrange[1], Xrange[2]))
      } else{
        v <- c()
        if(is.null(priorData)) P <- univMoTBF(Y, POTENTIAL_TYPE, domainChild, maxParam=maxParam)
        else P <- learnMoTBFpriorInformation(priorD, Y, s, POTENTIAL_TYPE, domainChild, maxParam=maxParam)$posteriorFunction 
        values <- list(Px=P, interval=c(Xrange[1], Xrange[2]))
      }
      v[[length(v)+1]] <- values
      return(v)
    } else{
      if(length(points)==1){
        v <- c(); X <- data[,nameParents]; X <- X[X<=points]
        values <- list(Px=Pp1[[1]], interval=c(Xrange[1], points))
        v[[length(v)+1]] <- values

        X <- data[, nameParents]; X <- X[X>points]
        values <- list(Px=Pp2, interval=c(points, Xrange[2]))
        v[[length(v)+1]] <- values
      } else{
        v <- c(); X <- data[,nameParents]
        dom <- c(Xrange[1], points, Xrange[2])
        Pp <- Pp1; Pp[[length(Pp)+1]] <- Pp2
        for(i in 1:length(Pp)){
          values <- list(Px=Pp[[i]], interval=c(dom[i], dom[i+1]))
          v[[length(v)+1]] <- values
        }
      }
      
      return(v) 
    }
  } else {
    B <- discreteVariablesStates(nameParents, data)[[1]]$states
    v <- c()
    for(i in 1:length(B)){
      Y1 <- Y[which(X==B[i])]
      if(is.character(Y1)){
        cutPoints <- unique(Y1)
        pD <- probDiscreteVariable(cutPoints, Y1)
        pos <- which(domainChild%in%cutPoints)
        coeff <- rep(0,length(domainChild))
        pD$coeff <- replace(coeff, pos, pD$coeff)
        names(pD$coeff) <- domainChild
        pD$sizeDataLeaf <- replace(coeff, pos, pD$sizeDataLeaf)
        values <- list(Px=pD, interval=B[i])
        v[[length(v)+1]] <- values
        
      } else{
        if(is.null(priorData)){
          P <- univMoTBF(Y1, POTENTIAL_TYPE, domainChild, maxParam=maxParam) 
        }else{
          priorChild <- priorData[,nameChild]
          if(ncol(priorData)!=1){
            priorParent <- priorData[,nameParents]
            priorD <- priorChild[which(priorParent==B[i])]
            P <- learnMoTBFpriorInformation(priorD, Y1, s, POTENTIAL_TYPE, domainChild)$posteriorFunction 
          } else {
            P <- univMoTBF(Y1, POTENTIAL_TYPE, domainChild, maxParam=maxParam)
          }
        }
        
        values <- list(Px=P, interval=B[i])
        v[[length(v)+1]] <- values
      } 
    }
    return(v)
  }
}

#'@rdname conditionalmotbf.learning
#'@export
BICscoreMoTBF <- function(conditionalfunction, data, nameParents, nameChild)
{
  X <- data[, nameParents]
  Y <- data[, nameChild]
  
  ## Compute the density values
  valuesT <- c(); DiscreteBN <- c(); nlevel <- c()
  for(i in 1:length(conditionalfunction)){
    if(is.numeric(X)){
      
      ##Continuous parent
      if(is.numeric(Y)){
        
        ## Continuous child
        domain <- Y[which((X>=conditionalfunction[[i]]$interval[1])&(X<=conditionalfunction[[i]]$interval[2]))]
        if(is.motbf(conditionalfunction[[i]]$Px)) values <- as.function(conditionalfunction[[i]]$Px)(domain)
        else values <- NULL
        if((length(values)!=length(domain))&&(is.numeric(conditionalfunction[[i]]$Px))) values <- rep(values, length(domain))
        valuesT <- c(valuesT,values)
      } else{
        
        ## Discrete child
        if(length(conditionalfunction[[i]]$Px)==1) DiscreteBN[[length(DiscreteBN)+1]] <- conditionalfunction[[i]]$Px[[i]]
        else DiscreteBN[[length(DiscreteBN)+1]] <- conditionalfunction[[i]]$Px
      } 
    } else 
      
      ## Discrete Parent
      if(is.numeric(Y)){ 
        
        ## Continuous child
        domain <- Y[which(X==conditionalfunction[[i]]$interval)] 
        if(is.motbf(conditionalfunction[[i]]$Px)) values <- as.function(conditionalfunction[[i]]$Px)(domain)
        else values <- NULL
        if((length(values)!=length(domain))&&(is.numeric(conditionalfunction[[i]]$Px))) values <- rep(values, length(domain))
        valuesT <- c(valuesT,values)
      } else {
        
        ## Discrete child
        if(length(conditionalfunction[[i]]$Px)==1) DiscreteBN[[length(DiscreteBN)+1]] <- conditionalfunction[[i]]$Px[[i]]
        else DiscreteBN[[length(DiscreteBN)+1]] <- conditionalfunction[[i]]$Px
      }
  }
  
  ## Compute the BIC
  if(is.numeric(X)){
    if(is.numeric(Y)){
      sizePx=c()
      for(i in 1:length(conditionalfunction)){
        if(!is.motbf(conditionalfunction[[i]]$Px)) tam <- 0
        else tam <- length(coef(conditionalfunction[[i]]$Px))
        sizePx <- c(sizePx, tam)
      }
      sizeP <- sum(sizePx)
      BiC <- sum(log(valuesT))-1/2*(sizeP+(length(conditionalfunction)-1))*log(length(valuesT))
    } else {
      BiC <- getBICDiscreteBN(DiscreteBN)
    }
  } else {
    if(is.numeric(Y)){
      sizePx=c()
      for(i in 1:length(conditionalfunction)){
        if(!is.motbf(conditionalfunction[[i]]$Px)) tam <- 0
        else tam <- length(coef(conditionalfunction[[i]]$Px))
        sizePx <- c(sizePx, tam)
      }
      sizeP <- sum(sizePx)
      valuesT <- valuesT[valuesT>0]
      BiC <- sum(log(valuesT))-1/2*(sizeP+(length(conditionalfunction)-1))*log(length(valuesT))
    } else {
      BiC <- getBICDiscreteBN(DiscreteBN)
    }
  }
  return(BiC)
}


#' BIC for Multiple Functions
#' 
#' Compute the BIC score for more than one function
#' 
#' @param Px A list with 'n' elements. Each element contains a \code{"motbf"} function.
#' @param X A list with 'n' elements. Each one contains a \code{"numeric"} vector with the 
#' values of the data set which correspond to the diferent functions.
#' @return The \code{"numeric"} BIC value.
#' @seealso \link{univMoTBF}
#' @export
#' @examples
#' ## Data
#' X <- rnorm(500)
#' Y <- rnorm(500, mean=1)
#' data <- data.frame(X=X, Y=Y)
#' ## Data as a "list"
#' Xlist <- sapply(data, list)
#' 
#' ## Learning as a "list"
#' Plist <- lapply(data, univMoTBF, POTENTIAL_TYPE="MOP")
#' Plist
#' 
#' ## BIC value
#' BICMultiFunctions(Px=Plist, X=Xlist)
#' 
BICMultiFunctions <- function(Px, X){
  coeffs <- unlist(sapply(1:length(Px), function(i) coef(Px[[i]])))
  totalValues <- unlist(sapply(1:length(Px), function(i) as.function(Px[[i]])(X[[i]])))
  BiC <- sum(log(totalValues))-(1/2*((length(coeffs)+length(Px))*log(length(totalValues))))
  return(BiC)
}

#' Plots for Conditional Functions
#' 
#' Get the graphical result of an MoTBF conditional function of two variables, i.e. a parent and his child.
#' 
#' @param conditionalFunction the output of the conditionalMethod. A list with the the interval of the parent
#' and the final MoTBF density function.
#' @param data The dataset used.
#' @param nameChild Name of the child variable in the conditional function. By default is \code{NULL}.
#' @param points Logical value. If \code{TRUE} the points of the data are overplotted.
#' @param color By default \code{NULL}, a selection of colors of the color palette of \R is used.
#' @details If the number of parents is bigger than one, then the message 
#' "It is not possible plotting the conditional function." is shown.
#' @return A plot of the conditional function.
#' @seealso \link{conditionalMethod}
#' @export
#' @examples
#' ## Data
#' X <- rnorm(1000)
#' Y <- rnorm(1000, mean=X)
#' data <- data.frame(X=X,Y=Y)
#' cov(data)
#' 
#' ## Conditional Learning
#' parent <- "X"
#' child <- "Y"
#' intervals <- 5
#' potential <- "MTE"
#' P <- conditionalMethod(data, nameParents=parent, nameChild=child, 
#' numIntervals=intervals, POTENTIAL_TYPE=potential)
#' plotConditional(conditionalFunction=P, data=data)
#' plotConditional(conditionalFunction=P, data=data, points=TRUE)
#' 
plotConditional <- function(conditionalFunction, data, nameChild=NULL, points=FALSE, color=NULL)
{
  ## Define the color
  if(is.null(color)) color <- colorRampPalette(c("#FFFFD9", "#EDF8B1", "#C7E9B4", "#7FCDBB", "#41B6C4", "#1D91C0", "#225EA8", "#253494", "#081D58"))
  
  nameParent <- conditionalFunction[[1]]$parent
  if(length(nameParent)==1){ 
    if(is.null(nameChild)){
      if(ncol(data)==2) nameChild <- colnames(data)[which(colnames(data)!=nameParent)]
      else return(cat("The name of the child variable is needed because the number of columns in the dataset is bigger than two."))
    }else{
      if(ncol(data)==2)
        if(colnames(data)[which(colnames(data)!=nameParent)]!=nameChild)
          return(cat("The name of the child variable has not been found in the dataset."))
    }
    X <- data[, nameParent]
    Y <- data[, nameChild]
    
    xgrid <- seq(min(Y),max(Y),length.out=50) 
    ygrid <- seq(min(X),max(X),length.out=50)
    griddata <- as.data.frame(expand.grid(xgrid,ygrid))
    D <- c()
    for(i in 1:length(conditionalFunction)){
      j <- which((griddata[,2]>=conditionalFunction[[i]]$interval[1])&(griddata[,2]<=conditionalFunction[[i]]$interval[2]))
      gdi <- griddata[j,]
      D <- c(D,as.function(conditionalFunction[[i]]$Px)(gdi[,1]))
    }
    d <- matrix(D, nrow = length(xgrid), ncol = length(ygrid), byrow=F)
    zlim <- range(d, finite=TRUE)
    zlim[1] <- 0
    nlevels <- 20
    levels <- pretty(zlim, nlevels)
    nlevels <- length(levels)

    filled.contour(x = ygrid, y = xgrid, z = t(d), nlevels = nlevels,
                   levels = levels, col = color(nlevels),                           
                   plot.title = {title(xlab = "X", ylab = "Y")},     
    )
    if(points){
      mar.orig <- par("mar")
      w <- (3 + mar.orig[2]) * par("csi") * 2.54
      layout(matrix(c(2, 1), ncol = 2), widths = c(1, lcm(w)))
      points(data[,c(nameParent, nameChild)])
      par(mfrow=c(1,1))
    }
  }else{
    return(cat("It is not possible plotting the conditional function."))
  }
}

#' Prints Conditional Functions
#' 
#' Prints the result of an MoTBF conditional function of two variables, i.e. a parent and his child.
#' 
#' @param conditionalFunction the output of the function \code{conditionalMethod}. A list with the interval of the parent
#' and the final \code{"motbf"} density function fitted in each interval.
#' @return The results  of the conditional function are shown.
#' @seealso \link{conditionalMethod}
#' @export
#' @examples
#' ## Data
#' X <- rexp(500)
#' Y <- rnorm(500, mean=X)
#' data <- data.frame(X=X,Y=Y)
#' cov(data)
#' 
#' ## Conditional Learning
#' parent <- "X"
#' child <- "Y"
#' intervals <- 5
#' potential <- "MOP"
#' P <- conditionalMethod(data, nameParents=parent, nameChild=child, 
#' numIntervals=intervals, POTENTIAL_TYPE=potential)
#' printConditional(P)
#' 
printConditional <- function(conditionalFunction)
{
  for(i in 1:length(conditionalFunction)){
    if((conditionalFunction[[i]]$parent)==(conditionalFunction[[1]]$parent)){
      if(length(conditionalFunction[[i]])==3){
        if(is.character(conditionalFunction[[i]]$interval)){
          cat("Parent:", conditionalFunction[[i]]$parent, "   \t Range =", paste("\"",conditionalFunction[[i]]$interval,"\"", sep=""),"\n")
        } else {
          cat("Parent:", conditionalFunction[[i]]$parent, "   \t Range:", conditionalFunction[[i]]$interval[1], "<",conditionalFunction[[i]]$parent,"<", conditionalFunction[[i]]$interval[2], "\n")
        }
        print(              conditionalFunction[[i]]$Px)
      } else{
        if(is.character(conditionalFunction[[i]]$interval)){
          cat("Parent:", conditionalFunction[[i]]$parent, "   \t Range =", paste("\"",conditionalFunction[[i]]$interval,"\"", sep=""),"\n")
        } else {
          cat("Parent:", conditionalFunction[[i]]$parent, "   \t Range:", conditionalFunction[[i]]$interval[1], "<",conditionalFunction[[i]]$parent,"<", conditionalFunction[[i]]$interval[2], "\n")
        }
      }
    }else{
      if(length(conditionalFunction[[i]])==3){
        if(is.character(conditionalFunction[[i]]$interval)){
          cat("Parent:", conditionalFunction[[i]]$parent, "   \t Range =", paste("\"",conditionalFunction[[i]]$interval,"\"", sep=""),"\n")
        } else {
          cat("Parent:", conditionalFunction[[i]]$parent, "   \t Range:", conditionalFunction[[i]]$interval[1], "<",conditionalFunction[[i]]$parent,"<", conditionalFunction[[i]]$interval[2], "\n")
        }
        print(              conditionalFunction[[i]]$Px)
      } else{
        if(is.character(conditionalFunction[[i]]$interval)){
          cat("Parent:", conditionalFunction[[i]]$parent, "   \t Range =", paste("\"",conditionalFunction[[i]]$interval,"\"", sep=""),"\n")
        } else {
          cat("Parent:", conditionalFunction[[i]]$parent, "   \t Range:", conditionalFunction[[i]]$interval[1], "<",conditionalFunction[[i]]$parent,"<", conditionalFunction[[i]]$interval[2], "\n")
        }
      }
    }
  }
}
