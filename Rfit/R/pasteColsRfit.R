#  Renamed pasteCols of library plotrix written by Jim Lemon et. al.
#  June 2011 under GPL 2

pasteColsRfit<- function (x, sep = "") {
    pastestring <- paste("list(", paste("x", "[", 1:dim(x)[1], 
        ",]", sep = "", collapse = ","), ")", sep = "")
    return(do.call(paste, c(eval(parse(text = pastestring)), 
        sep = sep)))
}

