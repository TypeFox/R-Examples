"preProcess" <-
function (data, sample = 1, sample_time = 1, sample_lambda = 1, 
    sel_time = vector(), sel_lambda = vector(), baselinetime = vector(), 
    baselinelambda = vector(), scalx = NULL, scalx2 = NULL, 
    sel_lambda_ab = vector(), sel_time_ab = vector(), rm_x2 = vector(), 
    rm_x = vector(), svdResid = list(), numV = 0, sel_special = list(),
    doubleDiff = FALSE, doubleDiffFile = "doubleDiff.txt") 
{
    ## Note that the options are _not_ compatible with each other
    ## call preProcess repeatedly for consistency
  
  dataN <- sample_sel(data = data, sample = sample, sample_time = 
                      sample_time, sample_lambda = sample_lambda, 
                      sel_time = sel_time, sel_lambda = sel_lambda, 
                      sel_lambda_ab = sel_lambda_ab,
                      sel_time_ab = sel_time_ab, sel_special = sel_special)
  if (length(baselinelambda) != 0) 
    dataN <- baseCorlambda(dataN, baselinelambda)
  if (length(baselinetime) != 0) 
    dataN <- baseCortime(dataN, baselinetime)
  if (length(rm_x2) != 0) { 
    dataN@psi.df <- dataN@psi.df[,-rm_x2] 
    dataN@x2 <-  dataN@x2[-rm_x2]
    dataN@nl <- dataN@nl - length(rm_x2)
  }
  if (length(rm_x) != 0) { 
    dataN@psi.df <- dataN@psi.df[-rm_x,] 
    dataN@x <-  dataN@x[-rm_x]
    dataN@nt <- dataN@nt - length(rm_x)
  }
  if (!is.null(scalx)) 
        dataN@x <- dataN@x * scalx
  if (!is.null(scalx2)) 
        dataN@x2 <- dataN@x2 * scalx2[1] + scalx2[2]
  if (numV > 0) {
    if(numV == 1)
      subtr <- as.matrix( as.matrix(svdResid$left[,1:numV]) * svdResid$value[1:numV]) %*% t(as.matrix(svdResid$right[1:numV,])) 
    else {
      leftscaled <- svdResid$left[,1:numV]
	    for(i in 1:numV) 
              leftscaled[,i] <- leftscaled[,i] * svdResid$value[i]
      subtr <- leftscaled %*% svdResid$right[1:numV,]
    }
    if(svdResid$weight) 
		subtr <- subtr / svdResid$weightM
    dataN@psi.df <- dataN@psi.df - subtr
  }
  if(doubleDiff) 
    dataN <- doubleDiff(dataN, file=doubleDiffFile)
  dataN@datCall <- append(data@datCall, match.call())
  dataN
}

