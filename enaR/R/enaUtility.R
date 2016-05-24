#' enautility --- utility analysis of a flow network
#' INPUT = network object
#' OUTPUT = list of utility statistics
#'
#' M. Lau | July 2011
#' ---------------------------------------------------

enaUtility <- function(x, type=c('flow','storage'),
                       eigen.check=TRUE,balance.override=FALSE,tol=10){
                                        #Missing Data Check
  if (any(is.na(x%v%'storage'))){
    warning('This function requires quantified storage values.')
  }else{
    orient <- get.orient()
                                        #Check for network class
    if (class(x) != 'network'){warning('x is not a network class object')}

                                        #set default for type == 'flow'
    if (length(type) > 1){type <- 'flow'}

                                        #Check for balancing
    if (balance.override){}else{
      if (any(list.network.attributes(x) == 'balanced') == FALSE){x%n%'balanced' <- ssCheck(x)}
      if (x%n%'balanced' == FALSE){warning('Model is not balanced'); stop}
    }
                                        #unpack data from x
    Flow <- t(as.matrix(x, attrname = 'flow')) #flows
    input <- x%v%'input' #inputs
    stor <- x%v%'storage' #storage values
    I <- Flow*0;diag(I) <- 1 #create the identity matrix
    T. <- input + apply(Flow,1,sum)
    FD <- Flow;
    diag(FD) <- -T.
                                        #
    if (type == 'flow'){
                                        #flow utilities
      D <- ginv(diag(T.)) %*% (FD - t(FD))
      rownames(D) <- colnames(D) <- colnames(Flow)

      if (eigen.check & abs(eigen(D)$values[1]) > 1){
        print(paste('Largest eigen value of D > 1:',eigen(D)$values[1]),quote=FALSE)
        out <- NA
      }else{
        U <- ginv(I-D) #non-dimensional integral flow utility
        U <- round(U,tol)
        Y <- diag(T.) %*% U #dimensional integral flow utility
                                        #indices
        synergism.F <- bcratio(Y) #flow benefit cost ratio (calls other function) (Synergism)
        mutualism.F <- bcratio(sign(Y)) # flow ratio of positive to negative signs )


        ns <- cbind('lam1D'=abs(eigen(D)$values[1]),'synergism.F' = synergism.F,'mutualism.F'=mutualism.F)
                                        #re-orient
        if (orient == 'rc'){
          D <- t(D)
          U <- t(U)
          Y <- t(Y)
        }else{}
        out <- list('D'=D,'U'=U,'Y'=Y,'ns'=ns) #pack output

      }

    }else if (type == 'storage'){
                                        #storage utilities
      x <- stor
      DS <- ginv(diag(x)) %*% (FD - t(FD))
      rownames(DS) <- colnames(DS) <- colnames(Flow)

      if (eigen.check & abs(eigen(DS)$values[1]) > 1){
        print(paste('Largest eigen value of DS > 1:',eigen(DS)$values[1]),quote=FALSE)
        out <- NA
      }else{

        US <- ginv(I - DS)
        US <- round(US,tol)
        YS <- diag(T.) %*% US
                                        #indices
        synergism.S <- bcratio(YS) #storage benefit cost ratio (calls other function) (Synergism)
        mutualism.S <- bcratio(sign(YS)) #storage ratio of positive to negative signs (Y/abs(Y) == sign of Y)

        ns <- cbind('lam1DS'=abs(eigen(DS)$values[1]),'synergism.S' = synergism.S,'mutualism.S'=mutualism.S)
                                        #re-orient
        if (orient == 'rc'){
          DS <- t(DS)
          US <- t(US)
          YS <- t(YS)
        }else{}
        out <- list('DS'=DS,'US'=US,'YS'=YS,'ns'=ns) #package output
      }
    }
                                        #labeling
    if (length(out)>1){
      for (i in 1:(length(out)-1)){
        if (class(out[[i]])=='matrix'){
          rownames(out[[i]]) <- colnames(out[[i]]) <- colnames(Flow)
        }
      }
    }
                                        #output
    return(out)
  }
}
