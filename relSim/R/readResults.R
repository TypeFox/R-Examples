readResults = function(N = 0, rel = "UN", gzip = TRUE,
                       strPath = "", strVer = "", fileName = NULL){

    if(is.null(fileName)){
        if(nchar(strVer)==0){
            fileName = sprintf("results-sim-%s-%d.csv", rel, N)
        }else{
            fileName = sprintf("results-sim-%s-%d-%s.csv", rel, N, strVer)
        }
    }

    if(nchar(strPath)>0){
        fileName = paste(strPath, fileName, sep = "")
    }

    if(gzip & !grepl("gz$", fileName))
        fileName = paste(fileName, ".gz", sep = "")

    res = matrix(scan(fileName, sep = ','), ncol = 3, byrow = T)

    return(data.frame(sib = res[,1], pc = res[,2], ibs = res[,3]))
}
