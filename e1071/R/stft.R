stft <- function(X, win=min(80,floor(length(X)/10)), 
                 inc=min(24, floor(length(X)/30)), coef=64, 
		 wtype="hanning.window")
  {
    numcoef <- 2*coef
    if (win > numcoef)
      {
	win <- numcoef
	cat ("stft: window size adjusted to", win, ".\n")
      }
    numwin <- trunc ((length(X) - win) / inc)

    ## compute the windows coefficients
    wincoef <- eval(parse(text=wtype))(win)

    ## create a matrix Z whose columns contain the windowed time-slices
    z <- matrix (0, numwin + 1, numcoef)
    y <- z
    st <- 1
    for (i in 0:numwin)
      {
	z[i+1, 1:win] <- X[st:(st+win-1)] * wincoef
	y[i+1,] <- fft(z[i+1,])
	st <- st + inc
      }

    Y<- list (values = Mod(y[,1:coef]), windowsize=win, increment=inc,
		  windowtype=wtype)
    class(Y) <- "stft"
    return(Y)
  }



