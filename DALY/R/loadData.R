## Helper function for 'readDALYdata'
## Load data into 'DALY' database

loadData <-
function(type, dst, str, data, readFile, x, y, O, D){
  if(type == "outcome")
    if(!is.na(as.character(readFile)))
      tclvalue(data) <- as.character(readFile)

  if(type == "data"){
    if(!is.na(as.character(readFile[8+O*(x-1)+D*(y-1)+1, 2])))
      tclvalue(dst) <- as.character(readFile[8+O*(x-1)+D*(y-1)+1, 2])
    if(!is.na(as.character(readFile[8+O*(x-1)+D*(y-1)+1, 3])))
      tclvalue(str) <- as.character(readFile[8+O*(x-1)+D*(y-1)+1, 3])
    if(30/(D-1) == 2)
      for (a in seq(15))
        for (b in seq(2))
         if(!is.na(as.double(as.character(readFile[9+O*(x-1)+D*(y-1)+a, b]))))
            data[[ceiling(a/3), 2*(a-1)%%3+b]] <-
            as.double(as.character(readFile[9+O*(x-1)+D*(y-1)+a, b]))
    if(30/(D-1) == 6)
      for (a in seq(5))
        for (b in seq(6))
         if(!is.na(as.double(as.character(readFile[9+O*(x-1)+D*(y-1)+a, b]))))
            data[[a, b]] <-
            as.double(as.character(readFile[9+O*(x-1)+D*(y-1)+a, b]))
    }
}