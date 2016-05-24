#Insert column in data frame at specific position
# Version:             0.1
# Date:         2015-03-01
#	Author:             T.S.

inbind <- function(frame,pos_after,data_vec){

  df_width <- NCOL(data_vec)
  if(is.vector(frame)){
    warning("Vectors will be coerced to data.frame.")
    frame <- data.frame(origvec=frame)
  }
  
  newframe <- data.frame(matrix(nrow=NROW(frame),ncol=NCOL(frame)+df_width))
  oldnames <- names(frame)
  
  ### Warnings & stops
  
  if(NROW(data_vec) != NROW(newframe)){
    stop("Wrong dimensions.")
  }
  
  if(is.matrix(frame)){
    warning("Matrix will be coerced to data.frame.")
  }
  
  if(is.null(oldnames)){
    warning("No column names provided. Create own.")
    oldnames <- paste0("X",1:NCOL(frame))
  }
  
  newnames <- character(NCOL(newframe))
  
  
  ## Main loop
  i <- 0
  while(i < NCOL(newframe)){
    i <- i + 1
    if(i <= pos_after){
      newframe[,i] <- frame[,i]
      newnames[i]  <- oldnames[i]
      } else if (i == pos_after+1){
        newframe[,i:(i+df_width-1)] <- data_vec
        if(!is.null(names(data_vec))){
          newnames[i:(i+df_width-1)] <- colnames(data_vec)
          } else {
            newnames[i:(i+df_width-1)] <- paste0("inserted",i:(i+df_width-1))
          }
        i <- (i+df_width-1)
        } else {
          newframe[,i] <- frame[,i-df_width]
          newnames[i] <- oldnames[i-df_width]
        }
  }
  
  names(newframe) <- newnames
  return(newframe)
}