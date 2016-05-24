allPerm = function(mcObj){
    if(class(mcObj) != "mc")
      stop("mcObject must be of class mc")
        
    r = mcObj$mc$allPerm()
    x = unlist(r)
    
    return(matrix(mcObj$elements[x], ncol = mcObj$length, byrow = TRUE))
}
