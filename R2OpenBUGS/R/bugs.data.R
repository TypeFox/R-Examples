"bugs.data" <- 
function(data, dir = getwd(), digits = 5, data.file = "data.txt"){
  if(is.numeric(unlist(data)))
            write.datafile(lapply(data, formatC, digits = digits, format = "E"), 
                file.path(dir, data.file))
  else {
            data.list <- lapply(as.list(data), get, pos = parent.frame(2))
            names(data.list) <- as.list(data)
            write.datafile(lapply(data.list, formatC, digits = digits, format = "E"), 
                file.path(dir, data.file))
  }
  return(data.file)
}
