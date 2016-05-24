"vec2mat" <- function(fct, no)
{
    parName <- names(formals(fct)[no])
    if (is.na(parName)) {stop("Argument number does not exist")}
    parName1 <- paste(parName, "[", sep = "")
    parName2 <- paste(parName, "[,", sep = "")

#    bodyStr <- as.character(body(fct))
#    if (bodyStr[1] == "{") {bodyStr <- bodyStr[-1]}  #  else {bodyStr <- bodyStr[1]}
    bodyStr <- deparse(body(fct))
    if (bodyStr[1] == "{") {bodyStr <- paste(head(tail(bodyStr, -1), -1), collapse = "")}  #  else {bodyStr <- bodyStr[1]}
    
    lenbs <- length(bodyStr)
    bsList <- list()
#    options(warn = -1)
    for (i in 1:lenbs)
    {
        tempText <- gsub(parName1, parName2, bodyStr[i], fixed = TRUE)
        bsList[[i]] <- paste(tempText, ";", sep = "")
    }
#    options(warn = 0)
    bodyStr2 <- paste("{", paste(as.vector(unlist(bsList)), collapse = ""), "}")
    bodyStr3 <- parse(text = bodyStr2)

    newfct <- fct
    body(newfct) <- bodyStr3  # bodyStr2

    return(list(newfct, bodyStr2, parName))
}

"nParm" <- function(bodyStr)
{
    gregObj <- gregexpr("(\\[,){1}[[:digit:]]+\\]{1}", bodyStr)[[1]]
    posVec <- gregObj
    lenVec <- attr(gregObj, "match.length")

    lenpv <- length(posVec)
    numVec <- rep(0, lenpv)
    for (i in 1:lenpv)
    {
        numVec[i] <- as.numeric(substr(bodyStr, posVec[i] + 2, posVec[i] + lenVec[i] - 2))
    }
    return(length(unique(numVec)))
}

"fParm" <- function(fct, no, fixed)
{
    v2m <- vec2mat(fct, no)
    if (all(is.na(fixed))) {return(v2m[[1]])}    
    
    bodyStr <- v2m[[2]]

    fStr <- paste("{ f <- c(", paste(fixed[!is.na(fixed)], collapse = ", "), ");")
    bodyStr <- paste(fStr, substr(bodyStr, 2, nchar(bodyStr)))
    parName <- v2m[[3]]

    gregObj <- gregexpr("(\\[,){1}[[:digit:]]+\\]{1}", bodyStr)[[1]]
    posVec <- gregObj
    lenVec <- attr(gregObj, "match.length")

    lenpv <- length(posVec)
    numVec <- rep(0, lenpv)
    
    realPos <- rep(NA, length(fixed))
    realPos[is.na(fixed)] <- 1:sum(is.na(fixed)) 
    
    fixPos <- rep(NA, length(fixed))
    fixPos[!is.na(fixed)] <- 1:sum(!is.na(fixed))    
    
#    options(warn = -1)    
    for (i in 1:lenpv)
    {
        numVec[i] <- as.numeric(substr(bodyStr, posVec[i] + 2, posVec[i] + lenVec[i] - 2))
        
        if (is.na(realPos[numVec[i]]))
        {
            inStr0 <- paste("f[", as.character(fixPos[numVec[i]]), "]", sep = "")
           
            lenBl <- nchar(parName) + lenVec[numVec[i]] - nchar(inStr0)
            inStr1 <- paste(rep(" ", nchar(parName) + lenVec[numVec[i]] - nchar(inStr0)), collapse = "")
            inStr2 <- paste(inStr0, inStr1, sep = "")
  
            substr(bodyStr, posVec[i] - nchar(parName), posVec[i] + lenVec[i] - 1) <- inStr2



#           inStr0 <- as.character(fixed[numVec[i]])
#           
#           lenBl <- nchar(parName) + lenVec[numVec[i]] - nchar(inStr0)
#           if (lenBl > 0) 
#           {
#               inStr1 <- paste(rep(" ", lenBl), collapse = "")
#               inStr2 <- paste(inStr0, inStr1, sep = "")
#           } else {
#               inStr2 <- inStr0
#           }         
#           substr(bodyStr, posVec[i] - nchar(parName), posVec[i] + lenVec[i] - 1) <- inStr2

        } else {
        
            inStr3 <- as.character(realPos[numVec[i]])
            
            ## In case the number changes from 10 to 9, 100 to 99 and so on
            if (nchar(inStr3) < nchar(as.character(numVec[i])) )
            {
                numSpaces <- nchar(as.character(numVec[i])) - nchar(inStr3)
                tempStr <- paste(rep(" ", numSpaces), collapse = "")
            
                inStr3 <- paste(tempStr, inStr3, sep = "")
            }
            substr(bodyStr, posVec[i] + 2, posVec[i] + lenVec[i] - 3) <- inStr3           
        }
    }
#    options(warn = 0)  
    # to avoid warnings when the replacement is shorter or longer: 1 or 10000 to replace x[,1]
    

    bodyStr2 <- parse(text = bodyStr)
    newfct <- fct
    body(newfct) <- bodyStr2

    return(newfct)
}


fct2list <- function(fct, no)
{
    v2m <- vec2mat(fct, no) 
    list(v2m[[1]], NULL, letters[1:nParm(v2m[[2]])])
}
