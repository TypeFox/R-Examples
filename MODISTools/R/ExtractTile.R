ExtractTile <- 
function(Data, Rows, Cols, Grid = FALSE)
{
  if(!is.object(Data)) stop("Data input should be an R object - numeric vector, matrix, or data frame.")
  
  if(!is.vector(Data) & !is.matrix(Data) & !is.data.frame(Data)){
    stop("Data should be a vector (one tile), or a matrix/data.frame (multiple tiles).")
  }
  
  if(is.vector(Data)) Data <- matrix(Data, nrow = 1, ncol = length(Data))
  if(is.data.frame(Data)) Data <- as.matrix(Data)
  
  if(!is.numeric(Data)) stop("Data is not numeric class: should be MODIS data only to extract a nested subset.")
  
  if(ncol(Data) <= 1) stop("Not enough pixels (columns) found to extract a subset.")
  
  if(!is.numeric(Rows) | !is.numeric(Cols)) stop("Rows and Cols should be both be numeric class - two integers.")
  if(length(Rows) != 2 | length(Cols) != 2) stop("Rows and Cols input must both be a vector of integers, with two elements.")
  if(abs(Rows[1] - round(Rows[1])) > .Machine$double.eps^0.5 | 
       abs(Rows[2] - round(Rows[2])) > .Machine$double.eps^0.5 |
       abs(Cols[2] - round(Cols[2])) > .Machine$double.eps^0.5 |
       abs(Cols[2] - round(Cols[2])) > .Machine$double.eps^0.5){
    stop("Size input must be integers.")
  }
  
  if((Rows[1] %% 2) != 1 | (Cols[1] %% 2) != 1) stop("The dimensions from any tile downloaded should be odd numbered")
  
  # Check Rows & Cols [1] == ncol Data, i.e. the length of data in a tile fits a matrix of dim Rows[1] & Cols[1]
  if(ncol(Data) != length(matrix(nrow=Rows[1], ncol=Cols[1]))) stop("Tile size of Data does not match Rows and Cols input.")
  
  if(((Rows[2] * 2) + 1) >= Rows[1] & ((Cols[2] * 2) + 1) >= Cols[1]) stop("Tile size requested is not smaller than Data.")
  
  if(!is.logical(Grid)) stop("Grid should be logical, to specify the format of the output.")
  #####
  
  # Get Data into a workable format and work out the subscripts of the nested subset.
  full.tile <- apply(Data, 1, function(x) list(matrix(x, nrow = Rows[1], ncol = Cols[1], byrow = TRUE)))
  centre <- c(ceiling(nrow(full.tile[[1]][[1]]) / 2), ceiling(ncol(full.tile[[1]][[1]]) / 2))
  row.range <- (centre[1] - Rows[2]):(centre[1] + Rows[2])
  col.range <- (centre[2] - Cols[2]):(centre[2] + Cols[2])
  
  # Put output in either array or matrix format.
  if(Grid){
    res <- array(dim = c( ((Rows[2] * 2) + 1), ((Cols[2] * 2) + 1), nrow(Data)))
    for(i in 1:nrow(Data)) res[ , ,i] <- full.tile[[i]][[1]][row.range,col.range]
  } else if(!Grid){
    res <- matrix(nrow = nrow(Data), ncol = length(matrix(nrow=((Rows[2] * 2) + 1), ncol = ((Cols[2] * 2) + 1))))
    for(i in 1:nrow(Data)) res[i, ] <- as.vector(full.tile[[i]][[1]][row.range,col.range])
  }
  
  return(res)
}