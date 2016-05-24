splitFrame <- structure(function#Splitting of multilevel data
### This function splits a Multilevel data frame into the specified level.
(
    rd, ##<<\code{data.frame} object with factor-level columns.
    lv = 1 ##<< \code{Numeric} or \code{character}. Position number in
           ##the factor-level columns of \code{rd}, or correspondant
           ##column name to split the data. If the spliting column is
           ##not a factor, the character name of the column should be used.
) {

    levs <- colclass(rd,TRUE)[['fac']]
    if(is.numeric(lv)) lv <- levs[lv]
    
    ff <- function(x,lv,as.num = TRUE){
        x[,lv] <- as.character(x[,lv])
        if(as.num)
        x[,lv] <- as.numeric(x[,lv])
        return(x)}
    if(is.null(lv))
        lv <- colclass(rd,TRUE)[['fac']][1]
    nf <- is.numeric(rd[,lv])
    if(nf)
    rd <- ff(rd,lv,as.num = FALSE)
    nmx <- names(rd)
    n. <- sapply(rd,is.numeric)
    nux <- names(rd)[n.]
    lt <- rev(names(rd)[!n.])
    lt <- lt[1:grep(lv,lt)]
    rds <- split(rd,rd[,lt],drop = TRUE)
    if(nf) rds <- Map(function(x)ff(x,lv),rds)
    return(rds)    
### \code{list} of \code{data.frame} objects.
} ,
ex=function() {
    ##Ring data frame:
    ##Multilevel data frame of tree-ring widths:
    data(Prings05,envir = environment())
    
    ## split the multilevel data into its second factor-level column:
    spl <- splitFrame(Prings05,2)
    str(spl)
    ## split the data into the factor-level: 'year':
    spl <- splitFrame(Prings05,'year')
    str(spl)
})
