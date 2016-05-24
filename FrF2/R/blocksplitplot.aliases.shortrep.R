## utility functions in support of valid and readable 
## alias info in the summary for block and splitplot designs

recalc.alias.block <- function(dia,leg){
## dia is the aliased element of a blocked or split plot design without the legend element
## leg is the separate legend element
    hilf <- lapply(dia, function(obj) strsplit(obj, ":", fixed=TRUE))
    leg <- lapply(leg, function(obj) unlist(strsplit(obj, "=", fixed=TRUE)))
    leg.c <- sapply(leg, function(obj) obj[1])
    leg.n <- sapply(leg, function(obj) obj[2])
    hilf <- lapply(hilf, function(obj) lapply(obj, function(obj2) sort(leg.c[which(leg.n %in% obj2)])))
    hilf <- lapply(hilf, function(obj) obj[ord(data.frame(sapply(obj, length)))])
    unlist(sort(unlist(sapply(hilf, function(obj) {
        ## vector of elements to be equated
        hilf <- paste(sapply(obj, function(obj2) paste(obj2, collapse="")), collapse="=")
        hilf[ord(data.frame(nchar(hilf), hilf))]
        }
        ))))
}

struc.aliased <- function(struc, nk, order){
        ## nk distinguishes Letters (nk<=50) or F1, F2, ...
        struc <- sapply(struc, function(obj) {
            if (length(grep("^-", obj)) == 1) {
                obj <- gsub("-", "~", obj)
                obj <- gsub("=", "=-", obj)
                obj <- gsub("=-~", "=", obj)
                obj <- gsub("~", "", obj)
            }
            obj
        })
        names(struc) <- NULL
        if (nk <= 50) {
            wme <- grep("^[[:alpha:]]=[[:alpha:][:punct:]]*", 
                struc)
            wme2 <- grep("^[[:alpha:]]{2}=[[:alpha:][:punct:]]*", 
                struc)
            if (order == 3) 
                wme3 <- grep("^[[:alpha:]]{3}=[[:alpha:][:punct:]]*", 
                  struc)
        }
        else {
            wme <- grep("^F[[:digit:]]+=F[[:digit:][:punct:]]*", 
                struc)
            wme2 <- grep("^F[[:digit:]]+:F[[:digit:]]+=F[[:digit:][:punct:]]*", 
                struc)
            if (order == 3) 
                wme3 <- grep("^F[[:digit:]]+:F[[:digit:]]+:F[[:digit:]]+=F[[:digit:][:punct:]]*", 
                  struc)
        }
        if (order == 2) 
            aus <- list(main = sort(struc[wme]), fi2 = sort(struc[wme2]))
        else aus <- list(main = sort(struc[wme]), fi2 = sort(struc[wme2]), 
            fi3 = sort(struc[wme3]))
    aus
}