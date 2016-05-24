"repChar" <- function(str, names, fixed, keep)  # used in 'mixdrc'
{
    if (is.null(fixed)) {fixed <- rep(NA, length(names))}

    "replaceChar" <- function(str, names, fixed, keep, sep=c(",", ";"))
    {
        lenK <- length(keep)
        lenN <- length(names)
        strVal <- str
    
        cutFrom <- rep(0, lenK)
        cutTo <- rep(0, lenK)
        keep2 <- rep("", lenK)
        for (i in 1:lenK)
        {

            cutFrom[i] <- regexpr(keep[i], strVal)  # matchVec[i]
            cutTo[i] <- cutFrom[i] + attr(regexpr(keep[i], strVal), "match.length") - 1
        
            keep2[i] <- paste(rep(sep[i], nchar(keep[i])), collapse="")
            substr(strVal, cutFrom[i], cutTo[i]) <- keep2[i]
        }
        #print(strVal)
    
        for (i in 1:lenN)
        {
            if (!is.na(fixed[i]))
            {
                strVal <- gsub(names[i], as.character(fixed[i]), strVal)
            }
        }
    
        for (i in 1:lenK)
        {
            cutFrom[i] <- regexpr(keep2[i], strVal)  # matchVec[i]
            cutTo[i] <- cutFrom[i] + attr(regexpr(keep2[i], strVal), "match.length") - 1    
            substr(strVal, cutFrom[i], cutTo[i]) <- keep[i]
        }
        return(strVal)
    }



    "buildFct" <- function(bodyStr, names, fixed)
    {
        argNames <- paste(names[is.na(fixed)], collapse=",")
        headerStr <- paste("function(DOSE," , argNames, "){(")
    
        fctStr <- paste(headerStr, bodyStr, "^lambda - 1)/lambda}")
        
        
        formStr <- paste("formula(respVar ~ opfct(doseVar,", argNames, "))")
        
#        print(fctStr)
        return(list(fctStr, formStr))        
#        return(eval(parse(text=fctStr)))
    }

    bodyS <- replaceChar(str, names, fixed, keep)
    return(buildFct(bodyS, names, fixed))
}
