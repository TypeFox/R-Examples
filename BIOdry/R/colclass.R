colclass <- structure(function#Column classes
### Column names in multilevel data frames are grouped into classes according to three level types: \code{numeric}, temporal, or \code{factor}. 
(
    rd, ##<<\code{data.frame} object with factor-level columns.
    as.list = TRUE ##<<\code{logical}. List the classes. If
                   ##\code{TRUE} then column names in data are
                   ##grouped in the form of a list.
) {
    
    tim <- c('second','hour','day','week',
             'month','year','time','millennium')
    
    n. <- sapply(rd,is.numeric)
    t. <- tim%in%names(rd)#;names(t.) <- tim
    tmp <- tim[t.]
    num <- names(n.)[n.]      
    fac <- names(n.)[!n.]
    num <- num[!num%in%tmp]
    fac <- fac[!fac%in%tmp]
    orv <- list(num = num,
                tmp = tmp, fac = fac)
    if(!as.list){
        orv <- c(orv, recursive = TRUE)
        names(orv) <- NULL}
    
    
    return(orv)
### \code{character} vector of the column names, or \code{list} of grouped
### column names.
} , ex=function() {
    ##Multilevel data frame of tree-ring widths:
    data(Prings05,envir = environment())
    ## getting grouped column names    
    classes <- colclass(Prings05)
    str(classes)
})
