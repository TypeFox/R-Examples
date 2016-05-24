MODISTimeSeries <-
function(Dir, Band, Simplify = FALSE)
{
  # DEFINE
  NUM_METADATA_COLS <- 10
  WHICH_ID <- 6
  
  if(!file.exists(Dir)) stop("Character string input for Dir argument does not resemble an existing file path.")
  
  file.set <- list.files(path = Dir, pattern = ".asc")
  
  file.ids <- sapply(file.path(Dir, file.set), function(x) 
                    any(grepl(pattern = Band, x = read.csv(file = x, header = FALSE, as.is = TRUE)[ ,WHICH_ID]))
              )
  file.set <- file.set[file.ids]
  
  if(length(file.set) < 1) stop("No downloaded files found in the requested directory for the requested data band.")
  
  data.collector <- vector(mode = "list", length = length(file.set))
  ts.row.names <- vector(mode = "list", length = length(file.set))
  ts.col.names <- vector(mode = "list", length = length(file.set))
  nrow.recorder <- ncol.recorder <- numeric(length = length(file.set))
  
  for(i in 1:length(file.set)){
    data.file <- read.csv(file.path(Dir, file.set[i]), header = FALSE, as.is = TRUE)
    names(data.file) <- c("nrow", "ncol", "xll", "yll", "pixelsize", "row.id", "product.code", "MODIS.acq.date",
                          "where", "MODIS.proc.date", 1:(ncol(data.file) - NUM_METADATA_COLS))
    data.file <- data.file[grepl(pattern = Band, x = data.file$row.id), ]
    
    data.collector[[i]] <- as.matrix(data.file[ ,(NUM_METADATA_COLS+1):ncol(data.file)])
    
    nrow.recorder[i] <- nrow(as.matrix(data.file[ ,(NUM_METADATA_COLS+1):ncol(data.file)]))
    ncol.recorder[i] <- ncol(as.matrix(data.file[ ,(NUM_METADATA_COLS+1):ncol(data.file)]))
    
    ts.col.names[[i]] <- paste(unique(data.file$where), "_pixel", 1:ncol.recorder[i], sep = "")
    ts.row.names[[i]] <- data.file$MODIS.acq.date
    colnames(data.collector[[i]]) <- ts.col.names[[i]]
    rownames(data.collector[[i]]) <- ts.row.names[[i]]
  }
  
  if(!Simplify) return(data.collector)
  
  if(!all(sapply(data.collector, nrow) == sapply(data.collector, nrow)[1])){
    cat('Simplify == TRUE, but not all tiles have the same number of rows so cannot be\n',
        'simplified into one matrix. Returning data as an array instead.\n', sep = '')
    return(data.collector)
  } else {
    res <- matrix(nrow = max(nrow.recorder), ncol = sum(ncol.recorder))
    rownames(res) <- ts.row.names[[which(nrow.recorder == max(nrow.recorder))[1]]]
    colnames(res) <- unlist(ts.col.names)
    for(j in 1:length(data.collector)){
      res[1:nrow.recorder[j],(sum(1, ncol.recorder[1:j]) - ncol.recorder[j]):sum(ncol.recorder[1:j])] <- 
        as.matrix(data.collector[[j]])
    }
    return(res)
  }
  
}