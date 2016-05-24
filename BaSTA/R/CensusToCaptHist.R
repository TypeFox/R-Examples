CensusToCaptHist <-
    function(ID, d, dformat = "%Y", timeInt = "Y") {
  # Check data
  if(class(ID) != "character") {
    ID <- as.character(ID)
  } 
  if (is.numeric(d)) {
    if (length(which(round(d) != d)) > 0) {
      stop("Please provide integer values or Date class ", 
          "values for agrument 'd'.", call. = FALSE)
    } else {
      int <- d
    }
  } else if (is.character(d) | class(d) == "Date") {
    if (is.character(d)) {
      d <- as.Date(d, format = dformat)
      if (length(which(is.na(d)))) {
        stop("Wrong 'dformat' argument or wrong 'd' values.", 
            call. = FALSE)
      }
    }
    if (timeInt == "Y"){
      int <- as.numeric(format(d, format = "%Y"))
    } else if (timeInt == "M") {
      int <- as.numeric(format(d, format = "%m")) + 
          12 * (as.numeric(format(d, format = "%Y")) - 
            min(as.numeric(format(d, format = "%Y"))))
    } else if (timeInt == "D" | timeInt == "W") {
      jul <- julian(d, origin = min(d)) + 1
      if (timeInt == "W"){
        int <- ceiling(jul / 7)
      } else {
        int <- jul
      }
    }
  } else {
    stop("Wrong class for argument 'd'. Values\n", 
        "need to be of class 'integer', 'character' or 'Date'.", 
        call. = FALSE)
  }  
  
  # Add a dummy individual seen in all years,
  # to cope with years where there are no recaptures
  dint <- min(int):max(int)
  ID  <- c(ID, rep("XdummyX", length(dint)))
  int  <- c(int, dint)
  
  mat <- as.matrix(table(ID, int))
  mat[mat > 0] <- 1
  ID <- rownames(mat)
  mat <- as.data.frame(cbind(as.factor(ID), as.matrix(mat)))
  colnames(mat) = c("ID", colnames(mat)[-1])
  
  # Remove the dummy row
  mat <- mat[-which(rownames(mat) == "XdummyX"), ]
  return(mat)
}