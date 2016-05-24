cdup <-
function(){setwd(sub("(.*)/.+$","\\1",getwd()));as.matrix(list.files())}

