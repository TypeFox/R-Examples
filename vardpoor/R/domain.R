 
domain <- function(Y, D, dataset=NULL) { 

   if(!is.null(dataset)) {
        dataset <- data.table(dataset)
        if (min(Y %in% names(dataset))!=1) stop("'Y' does not exist in 'dataset'!")
        if (min(Y %in% names(dataset))==1) Y <- dataset[, Y, with=FALSE] 
    
    if (!is.null(D)) {
        if (min(D %in% names(dataset))!=1) stop("'D' does not exist in 'data'!")
        if (min(D %in% names(dataset))==1) D <- dataset[, D, with=FALSE] }
   }

  name.Y <- substitute(Y)
  name.D <- substitute(D)
  
  # Y
  Y <- data.table(Y, check.names = TRUE)
  if (!all(sapply(Y, is.numeric))) stop(name.Y, " must be numerical")
  if (any(is.na(Y))) stop(name.Y, " has unknown values")
  n <- nrow(Y)

  # D
  D <- data.table(D, check.names = FALSE)
  if (any(duplicated(names(D))))
    stop(name.D, " has duplicate column names: ",
         paste(names(D)[duplicated(names(D))], collapse = ", "))
  if (nrow(D) != n) stop(name.Y, " and ", name.D ," have different row count")
  D <- D[, lapply(.SD, as.character), .SDcols = names(D)]

  Dom_agg <- unique(D)
  setkeyv(Dom_agg, names(Dom_agg))
  
  i <- k <- NULL 	
  domen <- foreach(i = 1:ncol(Y), .combine = data.table) %:%
    foreach(k = 1:nrow(Dom_agg), .combine = data.table) %do%
      ifelse(rowSums(D == Dom_agg[k, ][rep(1, n), ]) == ncol(D), Y[[i]], 0)
  
  if (!is.data.table(domen)) domen <- data.table(domen)
  namesD <- function(Y, D) {
    h <- vector(mode = "character", length = nrow(Dom_agg))
    for (i in 1:nrow(Dom_agg)) {
      cc <- paste(names(D), Dom_agg[i, ], sep = ".")
      h[i] <- paste(cc, collapse = "__")
    }
    foreach(i = 1:ncol(Y), .combine = c) %do% paste(names(Y)[i], h, sep="__")
  }
  setnames(domen, namesD(Y, D))
  domen <- data.table(domen, check.names=TRUE)
  return(domen)
}
