SweaveBook = function(idx = 0:6,
                      fileList = paste(paste("Chapter",idx, "/", sep = ""),
                                       paste("ch",idx,".rnw",sep=""),
                                       sep = "")){
    for(f in fileList){
        Sweave(f)
    }
}
