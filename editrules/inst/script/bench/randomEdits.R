# functions to generate some random catgegorical edits.


# make a set expression
setexpr <- function(var,cats){
    paste(var,'%in% c(',paste(cats,collapse=","),')')
}

# declare a subset FALSE
makeedit <- function(setexpr){
    paste('if (',setexpr,') FALSE')
}


# generate m random edits
genedits <- function(m,group){
    require(editrules)
    # use Boskovitz group A or B edits
    if ( group=='A' ){
        dk <- c(
            x.2.1 = 2, 
            x.2.2 = 2, 
            x.3.1 = 3, 
            x.3.2 = 3, 
            x.4.1 = 4, 
            x.4.2 = 4, 
            x.5.1 = 5, 
            x.5.2 = 5, 
            x.6.1 = 6, 
            x.6.2 = 6, 
            x.7.1 = 7, 
            x.7.1 = 7,
            x.8.1 = 8
            )
    } else if (group=='B') {
        dk <- c(
            x.2.1 = 2, 
            x.2.2 = 2, 
            x.3.1 = 3, 
            x.4.2 = 4, 
            x.5.1 = 5, 
            x.6.1 = 6, 
            x.7.1 = 7, 
            x.8.1 = 8, 
            x.9.1 = 9,
            x.10.1 = 10
            )
    }
    vars <- names(dk)
    dm <- c()
    for ( v in vars ) dm <- c(dm,setexpr(v,1:dk[v]))
    
    # generate random edits
    n <- length(dk)
    edits <- character(m)
    for ( i in 1:m ){
        invar <- sample(2:n,1)
        ivars <- sample(vars,invar,replace=FALSE)
        if ( i > 1 && !ivars1[1] %in% ivars ) ivars <- c(ivars,ivars[1])
        ivars1 <- ivars
        s <- character()
        for ( iv in ivars ){
            incat <- sample(1:(dk[iv]-1),1)
            icats <- sample(1:dk[iv],incat)
            s <- c(s,setexpr(iv,icats))
        }
        edits[i] <- makeedit(paste(s,collapse=" & "))
    }
    editarray(c(dm,edits))
}









