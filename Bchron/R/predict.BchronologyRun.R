predict.BchronologyRun = function(object, newPositions, newPositionThicknesses=NULL, maxExtrap=500, ...) {
  # This function takes a BchronologyRun object and produces new age predictions based on the values given in newPositions. If thicknesses are given as well it will produce values for that too
  # The output is a matrix of values where the number of rows is the number of stored iterations in object, and the number of columns is the number of positions given (averaged over thicknesses)
  
  # Check that, if thicknesses are given, they are the same length as newPositions
  if(!is.null(newPositionThicknesses)) {
    if(length(newPositionThicknesses)!=length(newPositions)) stop("newPositionThicknesses and newPositions must be of same length")
  }
  
  # Need some of the C functions for prediction
  predictInterp = function(alpha,lambda,beta,predictPositions,diffPositionj,currPositionsj,currPositionsjp1,thetaj,thetajp1) {
    return(.C('predictInterp',as.double(alpha),as.double(lambda),as.double(beta),as.double(predictPositions),as.integer(length(predictPositions)), as.double(diffPositionj), as.double(currPositionsj),as.double(currPositionsjp1), as.double(thetaj), as.double(thetajp1),as.double(rep(0,length(predictPositions))))[11][[1]])
  }
  predictExtrapUp = function(alpha,lambda,beta,predictPositions,currPositions1,theta1,maxExtrap,extractDate) {
    return(.C('predictExtrapUp',as.double(alpha),as.double(lambda),as.double(beta),as.double(predictPositions),as.integer(length(predictPositions)), as.double(currPositions1), as.double(theta1), as.integer(maxExtrap), as.double(extractDate),as.double(rep(0,length(predictPositions))))[10][[1]])
  }
  predictExtrapDown = function(alpha,lambda,beta,predictPositions,currPositionsn,thetan,maxExtrap) {
    return(.C('predictExtrapDown',as.double(alpha),as.double(lambda),as.double(beta),as.double(predictPositions),as.integer(length(predictPositions)), as.double(currPositionsn), as.double(thetan), as.integer(maxExtrap), as.double(rep(0,length(predictPositions))))[9][[1]])
  }
  
  # Get some useful things to start of
  nSamples = length(object$mu)
  out = matrix(ncol=length(newPositions),nrow=nSamples)
  p = 1.2
  alpha = (2-p)/(p-1)
  oldPositions = object$positions/object$positionScaleVal
  
  # Now loop through all the values in newPositions
  pb = utils::txtProgressBar(min = 1, max = nSamples*length(newPositions), style = 3,width=60,title='Calculating predictions... ')
  for(j in 1:length(newPositions)) {
    for(i in 1:nSamples) {
      utils::setTxtProgressBar(pb, (j-1)*nSamples+i)
      # Sample a current values of newPositions
      if(is.null(newPositionThicknesses)) {
        currPosition = newPositions[j]/object$positionScaleVal
      } else {
        currPosition = stats::runif(1,newPositions[j]/object$positionScaleVal-0.5*newPositionThicknesses[j]/object$positionScaleVal,newPositions[j]/object$positionScaleVal+0.5*newPositionThicknesses[j]/object$positionScaleVal)
      }
      
      # Get sedimentation rate parameters
      lambda = (object$mu[i]^(2-p))/(object$psi[i]*(2-p))
      beta = 1/(object$psi[i]*(p-1)*(object$mu[i]^(p-1)))
      
      # Find out if the currPosition is inside one of the old ones
      inside = findInterval(currPosition,oldPositions)
      
      # Now run appropriate interpolation/extrapolation routine
      if(inside %in% 1:(length(oldPositions)-1)) {
        diffPosition = oldPositions[inside+1]-oldPositions[inside]
        out[i,j] = round(predictInterp(alpha,lambda,beta,currPosition,diffPosition,oldPositions[inside],oldPositions[inside+1],object$theta[i,inside]/object$ageScaleVal,object$theta[i,inside+1]/object$ageScaleVal),3)
      } else if(inside==0) {
        out[i,j] = round(predictExtrapUp(alpha,lambda,beta,currPosition,oldPositions[1],object$theta[i,1]/object$ageScaleVal,maxExtrap,object$extractDate/object$ageScaleVal),3)
      } else if(inside==length(oldPositions)) {
        out[i,j] = round(predictExtrapDown(alpha,lambda,beta,currPosition,oldPositions[length(oldPositions)],object$theta[i,length(oldPositions)]/object$ageScaleVal,maxExtrap),3)
      } else {
        stop("Error in newPositions")
      }
    }
  }
  colnames(out) = paste('Pos',newPositions,sep='')
  return(out*object$ageScaleVal)
}