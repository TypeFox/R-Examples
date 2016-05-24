kernTest <-
function(kernType, numIn=4, tieParamNames=list(), testindex=NULL) {
  numData = 20

  kern <- kernCreate(numIn, kernType)
##   if exists([kern.type 'KernSetIndex'])==2 
##     for i = 1:length(kern.comp)
##       if rand(1)>0.5
##         indices = randperm(numIn);
##         indices = indices(1:ceil(rand(1)*numIn));
##         kern = kernSetIndex(kern, i, indices);
##       end
##     end
##   end
  if (length(dim(numIn)) > 0) {
    if (is.list(numIn)) {
      x <- numIn[[1]]
      x2 <- numIn[[2]]
    } else {
      x <- numIn
      x2 <- numIn
    }
  } else {
    if ( 'positiveTime' %in% names(kern) && kern$positiveTime ) {
      # For convolutional kernels starting at t=0 it does not make sense to use
      # negative inputs...
      x <- abs(matrix(rnorm(numData * numIn), numData))
      x2 <- abs(matrix(rnorm(numData/2 * numIn), numData/2))
    }
    else {
      x <- matrix(rnorm(numData * numIn), numData)
      x2 <- matrix(rnorm(numData/2 * numIn), numData/2)
    }
  }

#   kern <- modelTieParam(kern, tieParamNames)
  # Set the parameters randomly.
  params <- kernExtractParam(kern)
  params <- rnorm(length(params))/sqrt(rnorm(length(params))^2)
  kern <- kernExpandParam(kern, params)

  # Test for positive definiteness
  K <- kernCompute(kern, x);
  e <- tryCatch(eigen(K, only.values=TRUE)$values, error=function(e) NaN)
  if (!all(is.finite(e))) {
    cat('NaNs in the kernel, unable to test.\n')
    kernDisplay(kern)
    return(kern)
  }
  else if (!all(is.numeric(e))) {
    cat('Kernel has imaginary eigenvalues(!?!).\n')
  }
  else if (min(e) > 0)
    cat('The kernel is positive definite.\n')
  else
    cat(sprintf('The kernel is not positive definite: max eig %g, min eig %g\n',
                max(e), min(e)))

  if (!is.null(testindex)) {
    covGrad <- array(0, dim(K))
    covGrad[testindex[1], testindex[2]] = 1
    covGrad <- covGrad + t(covGrad)
  }
  else
    covGrad <- array(1, dim(K))
  epsilon <- 1e-6
  params <- kernExtractParam(kern)
  origParams <- params
  Lplus <- 0*params
  Lminus <- 0*params
  for ( i in seq(1, length(params)) ) {
    params <- origParams
    params[i] <- origParams[i] + epsilon
    kern <- kernExpandParam(kern, params)
    Lplus[i] <- sum(kernCompute(kern, x) * covGrad)
    params[i] <- origParams[i] - epsilon
    kern <- kernExpandParam(kern, params)
    Lminus[i] <- sum(kernCompute(kern, x) * covGrad)
  }
  params <- origParams
  kern <- kernExpandParam(kern, params)
  names <- names(kernExtractParam(kern, only.values=FALSE))
  gLDiff <- .5*(Lplus - Lminus)/epsilon
  g <- kernGradient(kern, x, covGrad)

  paramMaxDiff <- max(abs(gLDiff-g))
  if (paramMaxDiff > 2*epsilon) {
    l <- 0
    for (i in seq(1, length(names))) {
      if (l < length(names[[i]]))
        l <- length(names[[i]]);
    }
  
    cat(rep(' ', l), '\tanalytic   diffs     delta\n')
    for (i in seq(1, length(names))) {
      spaceLen = l - length(names[[i]]);
      space = rep(' ', spaceLen)
      cat(space, names[[i]],
          sprintf(':\t%4.6g\t%4.6g\t%4.6g\n', g[i], gLDiff[i], gLDiff[i] - g[i]))
    }
  }
## try 
##   Lplus = zeros(size(x));
##   Lminus = zeros(size(x));
##   gx = zeros(size(x));
##   origX = x;
##   for i = 1:size(x, 1)
##     for j = 1:size(x, 2)
##       x = origX;
##       x(i, j) = origX(i, j) + epsilon;
##       K = kernCompute(kern, x);
##       Lplus(i, j) =  full(sum(sum(K)));
##       LplusDiag(i, j) = full(trace(K));
##       x(i, j) = origX(i, j) - epsilon;
##       K = kernCompute(kern, x);
##       Lminus(i, j) = full(sum(sum(K)));
##       LminusDiag(i, j) = full(trace(K));
##     end
##     x = origX;
##     gx(i, :) = 2*sum(kernGradX(kern, x(i, :), x), 1);
##     gxDiag(i, :) = kernDiagGradX(kern, x(i, :));
##   end

##   gXDiff = .5*(Lplus - Lminus)/epsilon;
##   xMaxDiff = max(max(abs(gx-gXDiff)));
  
##   if xMaxDiff > 2*epsilon
##     fprintf('gX\n')
##     disp(gx)
##     fprintf('gXDiff\n')
##     disp(gXDiff)
##   end
  
##   gXDiagDiff = .5*(LplusDiag - LminusDiag)/epsilon;
##   xDiagMaxDiff = max(max(abs(gxDiag-gXDiagDiff)));
  
##   if xDiagMaxDiff > 2*epsilon
##     fprintf('gxDiag\n')
##     disp(gxDiag)
##     fprintf('gXDiagDiff\n')
##     disp(gXDiagDiff)
##   end
## catch
##   fprintf('kernGradX has an error.\n')
##   warning(lasterr)
##   xMaxDiff = 0;
##   xDiagMaxDiff = 0;
## end

## K = kernCompute(kern, x);
## traceK =  full(trace(K));
## K2 = kernDiagCompute(kern, x);
## traceK2 = full(sum(K2));
## traceDiff = traceK - traceK2; 
## %if abs(traceDiff) > 2*epsilon,
## %  fprintf('kernDiagCompute is not in sync with kernCompute.\n')
## %  fprintf('diag(kernCompute)\tkernDiagCompute')
## %  disp([diag(K), K2])
## %end

## covGrad = ones(size(kernCompute(kern, x, x2)));
## epsilon = 1e-6;
## params = kernExtractParam(kern) * paramPack;
## origParams = params;
## Lplus = zeros(size(params));
## Lminus = zeros(size(params));
## for i = 1:length(params);
##   params = origParams;
##   params(i) = origParams(i) + epsilon;
##   kern = kernExpandParam(kern, params * paramExpand);
##   Lplus(i) = full(sum(sum(kernCompute(kern, x, x2))));
##   params(i) = origParams(i) - epsilon;
##   kern = kernExpandParam(kern, params * paramExpand);
##   Lminus(i) = full(sum(sum(kernCompute(kern, x, x2))));
## end
## params = origParams;
## kern = kernExpandParam(kern, params * paramExpand);
## [void, names] = kernExtractParam(kern);
## names(toRemove) = [];
## gL2Diff = .5*(Lplus - Lminus)/epsilon;
## g = kernGradient(kern, x, x2, covGrad) * paramExpand;

## param2MaxDiff = max(max(abs(gL2Diff-g)));
## if param2MaxDiff > 2*epsilon
##   l = 0;
##   for i = 1:length(names)
##     if l < length(names{i})
##       l = length(names{i});
##     end
##   end
  
##   fprintf([char(repmat(32, 1, l)) '\tanalytic   diffs     delta\n']);
##   for i = 1:length(names)
##     spaceLen = l - length(names{i});
##     space = char(repmat(32, 1, spaceLen));
##     fprintf([space names{i} ':\t%4.6g\t%4.6g\t%4.6g\n'], ...
##             g(i), gL2Diff(i), gL2Diff(i) - g(i));
##   end
##   pause(0);
## end


## fprintf('Trace max diff: %2.6g.\n', traceDiff);
## fprintf('Param max diff: %2.6g.\n', paramMaxDiff)
## fprintf('Param X2 max diff: %2.6g.\n', param2MaxDiff)
## fprintf('X max diff: %2.6g.\n', xMaxDiff)
## fprintf('XDiag max diff: %2.6g.\n', xDiagMaxDiff)
## fprintf('\n');

## if nargout > 0
##   kernRet = kern;
## else
##   kernDisplay(kern);
## end

  kernDisplay(kern)
  return(kern)
# We don't test kernCompute(kern, x, x2) here at all!
}
