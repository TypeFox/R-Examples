createVennData <-
function(x, Cols = NULL, Splits = c(0, 0, 0), Labels = NULL, type = c('count','percent'), ToSkip = '000')
{ ### Generates data to feed into plotVenn3d() only.

  if (ncol(x) < 3) stop("Input data 'x' must have at least 3 columns.")

  # Determine which columns of data 'x' to use
  if (!is.null(Cols)) if (length(Cols) < 3) Cols <- NULL
  if (is.null(Cols)) Cols <- colnames(x)[1:3]
  
  # Split data 'x' into groups by 'Splits'
  out <- c(
    length(which(x[,Cols[1]] > Splits[1] & x[,Cols[2]] > Splits[2] & x[,Cols[3]] > Splits[3]  )), #1
    length(which(x[,Cols[1]] < Splits[1] & x[,Cols[2]] > Splits[2] & x[,Cols[3]] > Splits[3]  )),
    length(which(x[,Cols[1]] > Splits[1] & x[,Cols[2]] < Splits[2] & x[,Cols[3]] > Splits[3]  )), #3
    length(which(x[,Cols[1]] < Splits[1] & x[,Cols[2]] < Splits[2] & x[,Cols[3]] > Splits[3]  )),
    length(which(x[,Cols[1]] > Splits[1] & x[,Cols[2]] > Splits[2] & x[,Cols[3]] < Splits[3]  )), #5
    length(which(x[,Cols[1]] < Splits[1] & x[,Cols[2]] > Splits[2] & x[,Cols[3]] < Splits[3]  )),
    length(which(x[,Cols[1]] > Splits[1] & x[,Cols[2]] < Splits[2] & x[,Cols[3]] < Splits[3]  )), #7
    length(which(x[,Cols[1]] < Splits[1] & x[,Cols[2]] < Splits[2] & x[,Cols[3]] < Splits[3]  )))
  
  # Re-write as percents
  if (strtrim(tolower(type),1) == 'p') out <- round(out / sum(out) * 100, 2)
  
  # Name grouping vector
  names(out) <- c("111","011","101","001","110","010","100",'000')

  # There is always one group which will not be plotted - eliminate ToSkip group.
  # Generally, '000' should be eliminated.
  out <- out[!names(out) %in% ToSkip]

  # Specify labels for the 3 variables
  if (!is.null(Labels)) if (length(Labels) != 3) Labels <- NULL
  if (is.null(Labels)) Labels <- Cols

  return(list(x = out, labels = Labels))
}
