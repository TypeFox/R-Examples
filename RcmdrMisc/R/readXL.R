readXL <- function(file, rownames=FALSE, header=TRUE, na="", sheet=1, 
                   stringsAsFactors=default.stringsAsFactors()){
  data <- readxl::read_excel(path=file, sheet=sheet, col_names=header, na=na)
  class(data) <- "data.frame"
  if (rownames){
    rownames(data) <- data[, 1]
    data[[1]] <- NULL
  }
  colnames(data) <- make.names(colnames(data), unique=TRUE)
  if (stringsAsFactors){
    char <- sapply(data, class) == "character"
    for (var in which(char)){
      data[[var]] <- factor(data[[var]])
    }
  }
  data
}