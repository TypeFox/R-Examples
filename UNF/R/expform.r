.expform <- function(z, digits = 7) {
    y <- sprintf(paste('%+1.',digits,'e',sep=''), z)
    plus <- grep('e+', y, fixed=TRUE)
    neg <- grep('e-', y, fixed=TRUE)
    nas <- seq_along(y)[-c(plus,neg)]
    
    b <- numeric(length=length(z))
    b[plus] <- sapply(y[plus], function(i) strsplit(i, 'e+', fixed=TRUE)[[1]][1])
    b[neg] <- sapply(y[neg], function(i) strsplit(i, 'e-', fixed=TRUE)[[1]][1])
    b[c(plus,neg)] <- sapply(b[c(plus,neg)], function(i) {
        sp <- strsplit(as.character(i),'.',fixed=TRUE)[[1]]
        dec <- gsub('0+$','',sp[2])
        paste(sp[1], if(is.na(dec) | digits==1) '' else dec, sep='.')
    })
    
    e <- numeric(length=length(z))
    if(length(plus))
        e[plus] <- as.numeric(sapply(y[plus], function(i) strsplit(i, 'e+', fixed=TRUE)[[1]][2]))
    if(length(neg))
        e[neg] <- as.numeric(sapply(y[neg], function(i) strsplit(i, 'e-', fixed=TRUE)[[1]][2]))
    
    char <- character(length=length(y))
    char[plus] <- ifelse(e[plus]==0, paste(b[plus],'e+',sep=''), paste(b[plus],'e+',e[plus],sep=''))
    char[neg] <- ifelse(e[neg]==0, paste(b[neg],'e-',sep=''), paste(b[neg],'e-',e[neg],sep=''))
    char[nas] <- NA
    return(char)
}
