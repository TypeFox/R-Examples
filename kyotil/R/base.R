cbind.uneven=function(li) {
    # bind a list of data frame or named vector that are not of the same length
    allnames=lapply(li, rownames)
    alllen=lapply(allnames, length)
    nams = allnames[[which.max(alllen)]]
    nams= c(nams, setdiff(unique(unlist(allnames)), nams)) # append additional names
    
    res=NULL
    for (a in li){
        p=ncol(a)
        toadd = matrix(NA, nrow=length(nams), ncol=p, dimnames=list(nams,NULL))
        toadd[rownames(a),]=as.matrix(a)
        if (is.null(res)) {
            res=as.data.frame(toadd, stringsAsFactors=FALSE) # if stringsAsFactors is here in cbind, it won't work, that is why we have to do a if on is.null(res)
        } else {
            res=as.data.frame(cbind(res,toadd, stringsAsFactors=FALSE))
        }
    }
    res
}



# returns binary representation of an integer
binary<-function(i) if (i) paste(binary(i %/% 2), i %% 2, sep="") else "" 

# returns binary representatin of an integer with leading 0, the length of string is n
binary2<-function(i, n) {
    a<-2^((n-1):0)
    b<-2*a
    sapply(i,function(x) paste(as.integer((x %% b)>=a),collapse=""))
} 

unix=function (){
    substr(Sys.getenv("R_HOME") , 1,1)=="/"
}


#mysystem can call any exe file
mysystem = function (cmd, ...) {
    system ( paste(Sys.getenv("COMSPEC")," /c ",cmd) , invisible =TRUE, intern=FALSE, ...)
}


# convert temp from f to c
f2c = function (f) {
    (f-32)*5/9
}

# convert a factor to integer using its value, e.g. 1 to 1, 2 to 2
ftoi = function (f) {
    as.integer (as.character (f) )
}


# DONT CHANGE THIS, USED By RUMINEX!
# if ret.mat is set to true, always return a matrix
# in the output, each row corresponds to one element of X, instead of each column
mysapply=function (X, FUN, ..., simplify = TRUE, USE.NAMES = TRUE, ret.mat=TRUE) 
{
    if (is.null(names(X)) & is.numeric(X)) names(X)=X%+%""
    FUN <- match.fun(FUN)
    answer <- lapply(X, FUN, ...)
    if (USE.NAMES && is.character(X) && is.null(names(answer))) 
        names(answer) <- X
    if (simplify && length(answer) && length(common.len <- unique(unlist(lapply(answer, 
        length)))) == 1) {
         if (common.len >= 1) 
            if (common.len == 1 & !ret.mat) 
                unlist(answer, recursive = FALSE)
            else 
                t(array(unlist(answer, recursive = FALSE), dim = c(common.len, 
                    length(X)), dimnames = if (!(is.null(n1 <- names(answer[[1]])) & 
                    is.null(n2 <- names(answer)))) 
                    list(n1, n2)))
        else t(answer)
    }
    else t(answer)
}
## test mysapply
#sapply(1:3, function (i) if (i==2) rep(NA,2) else 1:3 )
#mysapply(1:3, function (i) if (i==2) rep(NA,2) else 1:3 )


# This function process all columns of x together instead of processing them one at a time
# FUN can return an array or a list. It does not have to return a scalar. This
#    saves from having to redo grouping for every col that has to be returned
#    also this eliminates the necessity to process a column of x at a time
myaggregate = function (x, by, FUN, new.col.name="aggregate.value", ...) 
{
    if (!is.data.frame(x)) 
        x <- as.data.frame(x)
    if (!is.list(by)) 
        stop(paste(sQuote("by"), "must be a list"))
    if (is.null(names(by))) 
        names(by) <- paste("Group", seq(along = by), sep = ".")
    else {
        nam <- names(by)
        ind <- which(nchar(nam) == 0)
        names(by)[ind] <- paste("Group", ind, sep = ".")
    }
    
    #original
    #y <- lapply(x, tapply, by, FUN, ..., simplify = FALSE)
    #modified
    z=mytapply(x, by, FUN, ...)
    
    #original
    #if (any(sapply(unlist(y, recursive = FALSE), length) > 1)) 
    #    stop(paste(sQuote("FUN"), "must always return a scalar"))
    #z <- y[[1]]
    
    d <- dim(z)
    w <- NULL
    for (i in seq(along = d)) {
        j <- rep.int(rep.int(seq(1:d[i]), prod(d[seq(length = i - 
            1)]) * rep.int(1, d[i])), prod(d[seq(from = i + 1, 
            length = length(d) - i)]))
        w <- cbind(w, dimnames(z)[[i]][j])
    }
    w <- w[which(!unlist(lapply(z, is.null))), , drop = FALSE]
        
    #original
    #y <- data.frame(w, lapply(y, unlist, use.names = FALSE))
    #modified
    if(nrow(w)!=nrow(matrix(unlist(z), nrow=nrow(w), byrow=TRUE))) stop("SOMETHING WRONG IN myaggregate")
    y <- data.frame(w, matrix(unlist(z), nrow=nrow(w), byrow=TRUE))
    #original
    #names(y) <- c(names(by), names(x))
    #modified
    names(y) <- c(names(by), new.col.name)
    y
}

# This function can handle X being matrix instead of just a vector
mytapply = function (X, INDEX, FUN = NULL, ..., simplify = TRUE) 
{
    FUN <- if (!is.null(FUN)) 
        match.fun(FUN)
    if (!is.list(INDEX)) 
        INDEX <- list(INDEX)
    nI <- length(INDEX)
    namelist <- vector("list", nI)
    names(namelist) <- names(INDEX)
    extent <- integer(nI)
    
    #original
    #nx <- length(X)
    #modified
    nx = ifelse(!(is.data.frame(X) | is.matrix(X)), length(X), length(X[,1]) )
    
    one <- as.integer(1)   
    group <- rep.int(one, nx)
    ngroup <- one
    for (i in seq(INDEX)) {
        index <- as.factor(INDEX[[i]])
        if (length(index) != nx) 
            stop("arguments must have same length")
        namelist[[i]] <- levels(index)
        extent[i] <- nlevels(index)
        group <- group + ngroup * (as.integer(index) - one)
        ngroup <- ngroup * nlevels(index)
    }
    if (is.null(FUN)) 
        return(group)
    ans <- lapply(split(X, group), FUN, ...)
    index <- as.numeric(names(ans))
    if (simplify && all(unlist(lapply(ans, length)) == 1)) {
        ansmat <- array(dim = extent, dimnames = namelist)
        ans <- unlist(ans, recursive = FALSE)
    }
    else {
        ansmat <- array(vector("list", prod(extent)), dim = extent, 
            dimnames = namelist)
    }
    names(ans) <- NULL
    ansmat[index] <- ans
    ansmat
}

# category.var is 
myreshapewide=function(formula, dat, idvar=NULL){
    tmp = as.character(formula)
    category.var=tmp[3]
    outcome.var=tmp[2]
    
    if (is.null(idvar)) {
        idvar=setdiff(names(dat), c(category.var,outcome.var))
        # if any of idvar has NA then it is a problem if not treated, here we opt to remove such columns
        idvar=idvar[sapply(idvar, function (idvar.) all(!is.na(dat[,idvar.])) )]
    } else {
        # remove those columns with changing values within an id
        # need [,-(1:length(idvar)),drop=FALSE] because the first columns are idvar
        tmp=apply(aggregate (x=dat[,!names(dat) %in% c(idvar,category.var,outcome.var),drop=FALSE], by=dat[,names(dat) %in% idvar,drop=FALSE], function(y) length(rle(y)$values)==1)[,-(1:length(idvar)),drop=FALSE], 
            2, all)
        varying.var=names(tmp)[!tmp]
        dat=dat[,!names(dat) %in% setdiff(varying.var, c(category.var,outcome.var)),drop=FALSE]
    }
    reshape(dat, direction="wide", v.names=outcome.var, timevar=category.var, idvar=idvar )
}

# keep column names
# return data frame if input is data frame
myscale=function(x){
    flag=is.data.frame(x)
    nam=dimnames(x)[[2]]
    aux=scale(x)
    dimnames(aux)[[2]]=nam
    if(flag) {
        tmp=as.data.frame(aux)
        mostattributes(tmp)=attributes(x)
        attr(tmp,"scaled:center")=attr(aux,"scaled:center")
        attr(tmp,"scaled:scale")=attr(aux,"scaled:scale")
#        
#        # add old attr
#        tmp1=attributes(aux)
#        for(a in tmp1){
#            #print(a)
#            print(names(a))
#            print("a")
#            #attr(tmp,"scaled:scale")=attr(x,"scaled:scale")
#        }
        aux=tmp
        names(aux)=names(x)
    }
    aux
}

meanmed=function(x, na.rm=FALSE){
    c(mean=mean(x, na.rm=na.rm), sd=sd(x, na.rm=na.rm), median=median(x, na.rm=na.rm), iqr=IQR(x, na.rm=na.rm), var=var(x, na.rm=na.rm))
}
read.tsv=function (file, header = TRUE, sep = "\t", ...) {
    read.csv(file, header = header, sep = sep, ...)
}

read.sv=function (file, header = TRUE, ...) {
    sep=","
    if (getExt(file)=="tsv") sep="\t"
    read.csv(file, header = header, sep = sep, ...)
}

keepWarnings <- function(expr) { 
    localWarnings <- list() 
    value <- withCallingHandlers(expr, 
    warning = function(w) { 
    localWarnings[[length(localWarnings)+1]] <<- w 
    invokeRestart("muffleWarning") 
    }) 
    list(value=value, warnings=localWarnings) 
} 

# make table that shows both counts/frequency and proportions
# style 1: count only; 2: count + percentage; 3: percentage only
table.prop=function (x,y=NULL,digit=1,style=2,whole.table.add.to.1=FALSE) {
    if (is.null(y)) {
        tbl=table(x)
        whole.table.add.to.1=TRUE  # to trick the computation of prop
    } else {
        tbl=table(x,y)
    }
    if (whole.table.add.to.1) prop = tbl/sum(tbl) else prop = apply (tbl, 2, prop.table)
    if (style==2) {
        res = tbl %+% " (" %+% round(prop*100,digit) %+% ")"
    } else if (style==3) {
        res = prop*100 # no need to do formating here
    } else res=tbl
    res=matrix(res, nrow=nrow(tbl))    
    dimnames(res)=dimnames(tbl)
    names(dimnames(res))=NULL
    
    if (is.null(y)) res=res[,1]
    
    res
}


# from Thomas, found on R user group
methods4<-function(classes, super=FALSE, ANY=FALSE){ 
    if (super) classes<-unlist(sapply(classes, function(cl) getAllSuperClasses(getClass(cl)))) 
    if (ANY) classes<-c(classes,"ANY") 
    gens<-getGenerics()@.Data 
    sigs<-lapply(gens, function(g) linearizeMlist(getMethods(g))@classes) 
    names(sigs)<-gens@.Data 
    sigs<-lapply(sigs, function(gen){ gen[unlist(sapply(gen, function(sig) any(sig %in% classes)))]}) 
    sigs[sapply(sigs,length)>0] 
} 
 

# p1 and p2 are two points. return y that corresponds to x on the line between p1 and p2
interpolate=function(pt1, pt2, x){
    slope=(pt2-pt1)[2]/(pt2-pt1)[1]
    intercept = pt1[2]-slope*pt1[1]
    intercept + slope * x    
}
