## There are only two functions that set the class of an object:
## permutation(), which uses class(P) <- "permutation", and cycle(),
## which uses class(x) <- "cycle".


## a *cyclist* is a list of cycles: list(1:4,8:9) is a cyclist; but
## this is informal.  The cycles are notionally distinct.

## a *cycle* object is a list of cyclists.


"word" <- function(M){
    ## takes a matrix and returns a word object; silently coerces to
    ## integer first
    stopifnot(is.matrix(M))
    storage.mode(M) <- "integer"
    if(nrow(M)>0){
        jj <- apply(M,1, singleword_valid)
    }
        class(M) <- c("permutation", "word")  # this is the *only*
                                              # time an object is
                                              # coerced to class
                                              # permutation or word.
        return(M)
}

"cycle" <- function(x){
    ## Function cycle() takes a list whose elements are lists whose
    ## elements are vectors (which are disjoint cycles); and returns
    ## an object of class "cycle".  It nicifies its input (eg removes
    ## length-1 cycles) before returning it.

    ## A use-case might be
    ## cycle(list(list(c(1,2,4),c(3,6)),list(c(1,2),c(3,4,5,6,7))))

    jj <- unlist(lapply(x, cyclist_valid))
        
    if(all(sapply(jj,isTRUE))){
        x <- lapply(x,nicify_cyclist)    
        class(x) <- c("permutation", "cycle")  ## NB this is the
                                               ## *only* place that
        ## class "cycle" is
                                               ## assigned to an
                                               ## object
        return(x)
    } else {
          stop(jj)
      }
}

is.id <- function(x){ UseMethod("is.id",x) }

is.id_single_cycle <- function(x){ is.null(unlist(x)) }

is.id.cycle <- function(x){ unlist(lapply(x,is.id_single_cycle)) }
is.id.word <- function(x){
    if(length(x)==0){return(logical(0))}
    if(size(x)==0){return(rep(TRUE,length(x)))}
    jj <- as.matrix(x)
    apply(jj == col(jj),1,all)
}

is.id.list <- function(x){ length(unlist(x))==0 }  # use for cyclists.

is.word <- function(x){ inherits(x, "word") }

is.cycle <- function(x){ inherits(x,"cycle") }

is.permutation <- function(x){ inherits(x,"permutaton") }

as.matrix.word <- function(x,...){unclass(x)}

names.word <- function(x){rownames(x)}

"names<-.word" <- function(x,value){
    rownames(x) <- value
    return(x)
}

"[.word" <- function(x, ...){
    x <- unclass(x)
    word(x[...,,drop=FALSE])
}

"[<-.word" <- function(x, index, value){
    out <- t(as.matrix(x))
    value <- t(as.matrix(as.word(value,size(x))))
    out[,index] <- value
    return(word(t(out)))
}

"[.cycle" <- function(x,...){
    x <- unclass(x)
    cycle(x[...])
}

#"[<-.cycle" <- function(x, index, value){
#    x <- unclass(x)
#    x[index] <- unclass(as.cycle(value))
#    browser()
#    return(cycle(x))  # sic -- not as.cycle(x), because x is a list of cyclists here.
#}

"c.word" <- function(...){
    a <- list(...)
    if(!all(unlist(lapply(a,is.word)))){
        stop("all arguments must be the same class")
    } else {
          n <- max(unlist(lapply(a,size)))
          a <- lapply(a,"size<-",n)
          word(do.call("rbind", a))
      }
}

"c.cycle" <- function(...){
    if(!all(unlist(lapply(list(...),is.cycle)))){
        stop("all arguments must be the same class")
    } else {
          return(cycle(unlist(list(...),recursive=FALSE)))
      }
}

addcols <- function(M,n){

    ##takes a matrix and adds columns [corresponding to fixed
    ## elements] so the returned value has 'n' columns.  Used by
    ## as.word(), so cannot coerce output to class word.

    if(nrow(M)==0){return(matrix(integer(0),0,n))}
    nm <- ncol(M)
    if(n>=nm){
        return(cbind(M,matrix(seq(from=1+nm,len=n-nm),nrow(M),n-nm,byrow=TRUE)))
    } else {
        stop("cannot remove columns")
    }
}

as.word <- function(x,n=NULL){

    ## This function is the user-friendly way to create a word object
    ## (compare word(), which is not terribly friendly).  Function
    ## as.word() does its best to coerce its argument to a word.
    ## Argument 'n' cannot act to reduce the size of the word, only
    ## increase it.  If you want to reduce the size, use trim() or
    ## tidy().  This function does not call word() except directly
    ## (e.g. it does not call size<-.word(), as this would give a
    ## recursion).

    if(is.word(x)){
        if(missing(n)){
           return(x)
        } else {
           size(x) <- n
          return(x)
        }
    } else if(is.cycle(x)){
        return(cycle2word(x,n)) 
    } else if(!is.numeric(x)){
        stop("can only coerce numeric objects to word")
    } else if(is.matrix(x)) {
        if(missing(n)){n <- ncol(x)}
        return(word(addcols(x,n)))
    } else if(is.vector(x)){
        if(missing(n)){n <- length(x)}
        return(word(addcols(t(x),n)))
    } else {
        warning("cannot coerce to class word")
        return(NA)
    }
}

print.word <- function(x, ...){  # contortions needed because x might have zero columns
    given <- x
    x <- unclass(x)
    if(is.null(rownames(x)) & length(x)>0){
        rownames(x) <- paste("[",seq_len(nrow(x)),"]",sep="")
    }
    if(ncol(x)>0){
        colnames(x) <- paste("{",seq_len(ncol(x)),"}",sep="")
    } else {
        cat(" {}")
    }
    jj <- x        
    jj[jj==col(jj)] <- '.'
    print(noquote(jj))
    return(invisible(given))
}

as.cycle <- function(x){   # does its best to coerce to cycle form.
                           # Takes character strings and permutation
                           # matrices

    if(missing(x)){
        return(id)
    } else if(is.cycle(x)){
        return(x)
    } else if(is.character(x)){
        return(char2cycle(x))
    } else if(is.vector(x)){
        return(cycle(list(list(x))))
    } else if(is.list(x) & all(unlist(lapply(x,is.vector)))){ # a cyclist
        return(cycle(list(x)))
    } else if(is.matrix(x)){  # includes words
        if(nrow(x)==0){return(nullcycle)}
        out <- apply(word(x),1,vec2cyclist_single)
        names(out) <- rownames(x)
        return(cycle(out))
    } else {
       stop("not recognised")
    }
}

cyc_len <- function(n){as.cycle(seq_len(n))}

char2cyclist_single <- function (x){

    if(all(unlist(strsplit(x,"")) != ",")) {#no commas anywhere
        commas <- ""
    } else {
        commas <- ","
    }
                                        
        jj <- lapply(
            strsplit(
                gsub(
                    "\\(", "",
                    unlist(strsplit(gsub(" ","",x),")"))
                    ), commas
                ),as.numeric)

    return(jj)
}

char2cycle <- function(char){
    ## char2cycle(c("(1,4)(6,7)","(3,4,2)(8,19)", "(56)","(12345)(78)","(78)"))
    cycle(sapply(char,char2cyclist_single,simplify=FALSE))
}

cycle2word <- function(x,n=NULL){  # cycle2word(as.cycle(1:5))
    if(is.null(n)){
        if(all(is.id(x))){
          n <- 0
        } else {
          n <- max(unlist(x,recursive=TRUE))}
    }
    word(do.call("rbind",lapply(x,cyclist2word_single,n=n)))
}

cyclist2word_single <- function(cyc,n){     #converts a cyclist to a single
                                        #permutation (vecor):
                                        #cycle2word_single(list(c(1,4,3),c(7,8)))

    if(length(unlist(cyc))==0){ return(seq_len(n)) }  # checking for the identity
    maxn <-  max(unlist(cyc,recursive=TRUE))

    if(missing(n)){
        n <- maxn
    } else {
        if(n<maxn){
            stop("supplied value of 'n' is too small")
        }
    }
        
    out <- seq_len(n)
    for(i in seq_along(cyc)){
        out[cyc[[i]]] <- shift(out[cyc[[i]]],-1)
    }
    return(out)
}

print.cycle <- function(x,...){  # x is a cycle.  Use case: print(cycle(list(x,y,z)))

    if((length(unlist(x))>0)){
        if(max(unlist(x,recursive=TRUE)) > 9){
            comma <- TRUE
        } else {
            comma <- FALSE
        }
    } else {
        comma <- FALSE
    }
    
    out <- unlist(lapply(x,as.character_cyclist,comma=comma))
    if(is.null(out)){
        cat("cycle(0)\n")
        return(out)
    } else {
        return(invisible(print(noquote(out))))
    }
}

as.character_cyclist <- function(y,comma=TRUE){
    
    ## Use case:
    ## as.character_cyclist(list(1:4,10:11,20:33)) x is a cyclist;
    ## as.character_cyclist(list(c(1,5,4),c(2,9)))
    ## as.character_cyclist(list(c(1,5,4),c(2,9)),comma=TRUE)
    
    if(length(y)==0){return("()")}
    
    if(comma){s <- ","} else {s <- ""}
    paste(sapply(y,function(u){paste(paste("(",paste(u,collapse=s),sep=""),")",sep="")}),collapse="")
}

as.character.cycle <- function(x,...){
    stopifnot(is.cycle(x))
    unlist(lapply(x,function(x){as.character_cyclist(x)}))
}

standard_cyclist <- function(x,n=NULL){

    ## standard representation as defined by Stanley, p30.  NB
    ## standard_cyclist() retains length-one cycles (compare
    ## nicify_cyclist(), which does not).

    ## standard_cyclist(list(c(4, 6), c(7), c(2, 5, 1), c(8, 3)))

    xvec <- unlist(x,recursive=TRUE)
    if(is.null(n)){n <- max(xvec)}
    jj <- seq_len(n)
    nicify_cyclist(c(x,as.list(jj[!(jj %in% xvec)])),rm1=FALSE,smallest_first=FALSE)
}

standard <- function(cyc,n=NULL){

    ## Take an object of class cycle and returns a list of cyclists.
    ## NB does not return a cycle object because cycle() calls
    ## nicify().
   
    ## Use-cases:
    ## standard(c(as.cycle(1:3),as.cycle(2:3) + as.cycle(6:7)))
    
    cyc <- as.cycle(cyc)
    xvec <- unlist(cyc,recursive=TRUE)
    if(is.null(n)){n <- max(xvec)}
    
    lapply(cyc,standard_cyclist,n=n)
}

fbin_single <- function(vec){  # takes a vector: fbin_single(sample(9))
    cycle(list(split(vec,cumsum(vec == cummax(vec)))))
}

fbin <- function(W){  # use-case: fbin(rperm(30,9))
    W <- as.matrix(as.word(W))
    cycle(unlist(apply(W,1,fbin_single),recursive=FALSE))
}

fbin_inv <- function(cyc){  # use-case: fbin_inv(as.cycle(rperm(30,9)))
    cyc <- as.cycle(cyc)
    f <- function(x){c(x,recursive=TRUE)}
    word(do.call("rbind",lapply(standard(cyc),f)))
}

nicify_cyclist <- function(x,rm1=TRUE,smallest_first=TRUE){  # needs rm1 argument for
                                         # standard_cyclist()

    ## takes a *cyclist* and puts it in a nice form (does not alter
    ## the permutation).  Note that nicify_cyclist() removes
    ## length-one cycles (compare standard_cyclist(), which does not).

    ## NB: nicify_cyclist() is called automatically by cycle().
    ## Remember that nicify_cyclist() takes a cyclist!


    ## use-cases:
    ## nicify_cyclist(list(c(4, 6), c(7), c(2, 5, 1), c(8, 3)),smallest_first=FALSE,rm1=FALSE)
    ## nicify_cyclist(list(c(4, 6), c(7), c(2, 5, 1), c(8, 3)),smallest_first=FALSE,rm1=TRUE )
    ## nicify_cyclist(list(c(4, 6), c(7), c(2, 5, 1), c(8, 3)),smallest_first=TRUE ,rm1=FALSE)
    ## nicify_cyclist(list(c(4, 6), c(7), c(2, 5, 1), c(8, 3)),smallest_first=TRUE ,rm1=TRUE )
   
    
    if(isTRUE(rm1)){  # remove singletons
        x <- remove_length_one(x)  
    }
    if(smallest_first){
        f <- which.min
    } else {
        f <- which.max
    }

    out <- lapply(x,function(o){shift(o,1-f(o))})
    order_wanted <- order(sapply(out,function(o){o[1]}))
    out <- 
        do.call("list", sapply(order_wanted, function(i){out[[i]]},simplify=FALSE))  # sort it by first [that is, the largest] element
    return(out)
}

#nicify <- function(x){
#    cycle(lapply(x,nicify_cyclist))
#}

remove_length_one <- function(x){
    x[unlist(lapply(x,function(u){length(u)>1}))]
}

vec2cyclist_single_cpp <- function(p){
    stop("vec2cyclist_single_cpp() not written yet")
}

vec2cyclist_single <- function(p){

    ## converts a (vector!) permutation into cycle form and returns a
    ## *list*.  Just a list! A cyclist! The elements of this list are
    ## the [disjoint] cycles.  Note the redundancies inherent:
    ## firstly, because the cycles commute, their order is immaterial
    ## (and a list is ordered); and secondly, the cycles themselves
    ## are invariant under cyclic permutation.  Heigh hoo.

    ## test: 793586142 -> (17)(458)(29)(3)
    ## as.cycle(c(7,9,3,5,8,6,1,4,2))

    n <- length(p)  #NB min(p) = 1 (not 0, off-by-one)
    out <- list()
    not_done <- rep(TRUE,length(p)) 
    while(any(not_done)){
        f <- min(which(not_done))  # first in bracket
       neew <- u <- f
       not_done[neew] <- FALSE
       while(u != (neew <- p[neew])){
           not_done[neew] <- FALSE
           f <- c(f,neew)
       }
        if(length(f)>1){out <- c(out,list(f))}
    }
    out    # NB a list whose elements are vectors which represent the cycles
}

inverse <- function(x){ UseMethod("inverse",x) }

inverse_word_single <- function(W){
    W[W] <- seq_along(W)
    return(W)
}

inverse_cyclist_single <- function(cyc){  # takes a cyclist, returns a cyclist
    ## use case:  inverse_cyclist_single(list(c(4, 6), c(2, 5, 1), c(8, 3)))

    lapply(cyc,function(o){c(o[1],rev(o[-1]))})
}

inverse.word <- function(x){   # takes a word, returns a word.  inverse.word(rperm(8,5))
    x <- as.word(x)
    word(t(apply(x,1,inverse_word_single)))
}

inverse.cycle <- function(x){ cycle(lapply(x,inverse_cyclist_single)) }

rperm <- function(n,r,moved=NA){
    if(is.na(moved)){
        return(word(matrix(replicate(n,sample(seq_len(r))),n,r,byrow=TRUE)))
    } else {
        f <- function(moved){
            out <- seq_len(r)
            jj <- sample(r,moved)
            out[jj] <- sample(jj)
            return(out)
        }
        return(as.word(matrix(replicate(n,f(moved)),n,r,byrow=TRUE)))
    }
}

"shape" <- function(x,drop=TRUE,id1=TRUE){
    x <- as.cycle(x)
    out <- lapply(x,shape_cyclist,id1=id1)
    if(drop & (length(x)==1)){
        out <- unlist(out)
    }
    return(out)
}

shape_cyclist <- function(cyc,id1=TRUE){  # use case: shape_cyclist(list(1:4,8:9))
    out <- unlist(lapply(cyc,length))
    if(id1 & is.null(out)){
        return(1)
    } else {
        return(out)
    }
}

shapepart_cyclist <- function(cyc,n=NULL){
    nmax <-  max(unlist(cyc,recursive=TRUE))
    if(is.null(n)){ n <- nmax } 
    if(n<nmax){stop("value of n too small")}
    out <- rep(0,n)
    for(i in seq_along(cyc)){out[cyc[[i]]] <- i}
    out[out==0] <- seq(from=max(out)+1,len=sum(out==0))
    return(out)
}

shapepart <- function(x){
    x <- as.cycle(x)
    out <- do.call("cbind",lapply(x,shapepart_cyclist,n=size(x)))
    colnames(out) <- names(x)
    as.partition(out)
}

size <- function(x){ UseMethod("size",x) }

size.word <- function(x){ # size(word(
    ncol(as.matrix(x))
}

size.cycle <- function(x){
    if(all(is.id(x))){return(0)}
    max(unlist(x,recursive=TRUE))
}

"size<-" <- function(x, value){ UseMethod("size<-") }

"size<-.word" <- function(x, value){   # trim it down then build
                                        # it up if necessary;
                                        # compare addcols() which
                                        # works purely on
                                        # matrices.
    stopifnot(is.word(x))
    return(word(addcols(trim(x),value)))
}

"size<-.cycle" <- function(x, value){
    stop("you cannot alter the size of a cycle")
}

"length<-.permutation" <- function(x,value){
    stop("you cannot change the length of a permutation")
}

length.word <- function(x){ nrow(x) }

trim <- function(x){
#    stop("problems: trim(as.word(1:6)) should return the empty word, but doesn't")
    stopifnot(is.word(x))
    if(length(x)==0){return(nullword)}
    if(all(is.id(x))){return(x)}
    y <- as.matrix(as.word(x))
    n <- ncol(y)
    jj <- apply(y,2,max)
    fix <- (jj==apply(y,2,min)) & (jj==seq_len(n))
    if(fix[n]){
        lose <- sum(cumprod(as.numeric(rev(fix)))>0)
        return(word(y[,seq_len(n-lose),drop=FALSE]))
    } else {
        return(x)
    }
}
        
fixed <- function(x){ UseMethod("fixed",x) }    

fixed.word <- function(x){  # fixed(word(t(c(2,3,5,4,1))))
    x <- as.matrix(x)
    if(nrow(x)>0){
        return(apply(x== col(x),2,all))
    } else {
        return(logical(0))
    }
}

fixed.cycle <- function(x){  # fixed(as.cycle(1:3) + as.cycle(10:11))
    n <- size(x)
    jj <- unlist(x,recursive=TRUE)
    !(seq_len(n) %in% jj)
}

"tidy" <- function(x){
    x <- as.word(x)
    x <- as.matrix(x)[,!fixed(x),drop=FALSE]
    if(nrow(x)==0){return(nullword)}
    n <- seq_len(ncol(x))
    u <- unique(sort(x))
    ind <- 0*n
    ind[u] <- n
    x[] <- ind[x]
    word(x)
}

rep.permutation <- function(x, ...){
    u <- seq(length.out = length(x))
    return(x[rep(u, ...)])
}

sgn <- function(x){
    x <- as.cycle(x)
    unlist(lapply(shape(as.cycle(x)),function(o){ifelse(is.null(o),1,1-2*sum(o-1)%%2)}))
}

is.even <- function(x){sgn(x)==1}

are_conjugate_single <- function(a,b){  # difficulties arise with the identity.
    stopifnot((length(a)==1) & (length(b)==1))
    if(is.id(a) & is.id(b)){
       return(TRUE)
    } else if(xor(is.id(a), is.id(b))){
       return(FALSE)
    } else {
       return(identical(unname(sort(shape(a))),unname(sort(shape(b)))))
    }
}

are_conjugate <- function(x,y){
    jj <- helper(x,y)
    apply(jj,1,function(ind){are_conjugate_single(x[ind[1]], y[ind[2]])})
}

as.function.permutation <- function(x,...){
    a <- NULL # to suppress the warning about 'no visible binding for global variable 'a' 
    x <- as.matrix(as.word(x))
    as.function(alist(a=, x[,a]))
}
           
#commutator_single <- function(x,y,n){
#    out <- inverse(x)*inverse(y)*x*y
#    if(!missing(n)){size(out) <- n}
#    return(out)
#}

commutator <- function(x,y){
    n <- max(size(x),size(y))
    jj <- helper(x,y)
    e1 <- as.matrix(as.word(x,n))
    e2 <- as.matrix(as.word(y,n))

    ##    f <- function(ind){commutator_single(e1[ind[1]],e2[ind[2]],n=n)}

    f <- function(ind){
        j1 <- e1[ind[1],]
        j2 <- e2[ind[2],]
        word_prod_single(
            word_prod_single(
                word_prod_single(
                    inverse_word_single(j1),
                    inverse_word_single(j2)
                    ),j1),j2)
    }
    
    return(as.word(t(apply(jj,1,f))))
}

permorder <- function(x,singly=TRUE){
    jj <- shape(x,id1=TRUE,drop=FALSE)
    f <- function(n){mLCM(c(1,n))} # needed because mLCM(5) fails
    if(singly){
        return(unlist(lapply(jj,f)))
    } else {
        return(f(as.vector(unlist(jj)))) # as.vector() needed to return an unnamed integer
    }
}

is.derangement <- function(x){
    x <- as.word(x)
    n <- seq_len(size(x))
    if(size(x)==0){   # identity element is not a derangment
        out <- rep(FALSE,length(x))
        names(out) <- names(x)
        return(out)
    }
    apply(as.matrix(x),1,function(u){all(u!=n)})
}

permprod <- function(x){
    out <- id
    x <- as.word(x)
    for(i in seq_along(x)){
        out <- out*x[i]
    }
    return(out)
}    

"get1" <- function(x,drop=TRUE){
    out <- lapply(as.cycle(x),function(u){unlist(lapply(u,min))})
    if(drop & (length(x)==1)){ out <- out[[1]] }
    return(out)
}

"get_cyc" <- function(x,elt){
    f <- function(u){ u[unlist(lapply(u,function(v){elt %in% v}))] }
    cycle(lapply(as.cycle(x), f))
}

"orbit_single" <- function(c1,n1){  # c1 is a cyclist, n1 an integer vector
    unlist(c1[which(unlist(lapply(c1,function(x){n1 %in% x})))])
}

"orbit" <- function(cyc,n){
    cyc <- as.cycle(cyc)
    jj <- helper(cyc,n)
    apply(jj,1,function(ind){orbit_single(unlist(unclass(cyc[ind[1]]),recursive=FALSE),n[ind[2]])})
}
   
