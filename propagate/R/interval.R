interval <- function(df, expr, seq = 10, plot = TRUE)
{
  if (!is.data.frame(df)) df <- as.data.frame(df)
  
  ## check for matching column/variable names
  VARS <- all.vars(expr)
  m <- match(VARS, colnames(df))
  if (any(is.na(m))) 
    stop("Variable names of input dataframe and expression do not match!")
  if (length(unique(m)) != length(m)) 
    stop("Some variable names are repetitive!")
  
  ## convert to list
  LIST <- as.list(df)
  
  ## create sequence list items
  if (is.numeric(seq)) LIST <- lapply(LIST, function(x) seq(x[1], x[2], length.out = seq))
    
  ## add 0 if x in [-a, b] 
  LIST <- lapply(LIST, function(x) if (x[1] < 0 & x[length(x)] > 0) c(x[x < 0], 0, x[x > 0]) else x)  
  
  ## create combinations grid
  GRID <- expand.grid(LIST) 
  
  ## evaluate combinations
  EVAL <- try(eval(expr, GRID), silent = TRUE) 
  if (inherits(EVAL, "try-error")) {
    print("Using 'vectorized' evaluation gave an error. Switching to 'row-wise' evaluation...")
    EVAL <- apply(GRID, 1, function(x) eval(expr, envir = as.list(x)))     
  }
  
  ## plot evalutions
  if (plot) {
    plot(1:length(EVAL), EVAL, type = "l", xlab = "Evaluation Number", ylab = "Evaluation Value")
    abline(h = c(min(EVAL, na.rm = TRUE), max(EVAL, na.rm = TRUE)), col = "blue", lwd = 1.5)
  }
    
  ## min/max evaluation
  OUT <- c(min(EVAL, na.rm = TRUE), max(EVAL, na.rm = TRUE))  
    
  class(OUT) <- "interval"
  OUT
}






