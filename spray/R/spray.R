`is_valid_spray` <- function(L){
    stopifnot(is.list(L))
    stopifnot(length(L)==2)
    if(is.empty(L)){
        return(TRUE)
    } else {
        stopifnot(is.matrix(L[[1]]))
        stopifnot(is.vector(L[[2]]))
        stopifnot(nrow(L[[1]])==length(L[[2]]))
        return(TRUE)
    }
}

`is.spray` <- function(S){inherits(S,"spray")}

`spraymaker` <- function(L,addrepeats=FALSE){    # formal; this is the *only* way to create a spray object; L = list(M,x)
    stopifnot(is_valid_spray(L))
    if(is.empty(L)){
        out <- L
    } else {
        M <- L[[1]]
        x <- L[[2]]
        if(!addrepeats & anyDuplicated(M)){ 
            stop("repeated indices; yet argument 'addrepeats' is FALSE")
        }
        out <- spray_maker(M,x)   # see RcppExports.R; strips out zeros and puts index rows in (some inscrutable) order
    }
    class(out) <- "spray"  # this is the *only time class<-() is called
    return(out)
}

`is.empty` <- function(L){ is.null(L[[1]]) & is.null(L[[2]])}
`is.zero` <- is.empty

`spray` <- function(M,x,addrepeats=FALSE){
    if(is.vector(M)){M <- rbind(M)}
    if(inherits(M,"partition")){M <- t(M)} # saves typing t() in 'spray(parts(4))' which drives me nuts
    M <- as.matrix(M) # e.g., coerces expand.grid(...) to a matrix
    if(missing(x)){x <- rep(1,nrow(M))}
    if(length(x)==1){x <- rep(x,nrow(M))}
    return(spraymaker(list(M,x),addrepeats=addrepeats))
}

`index` <- function(S){S[[1]]}    # these two functions are the only
`value` <- function(S){S[[2]]}    # 'accessor' functions in the package

`as.spray` <- function(arg1, arg2, addrepeats=FALSE, offbyone=FALSE){  # tries quite hard to coerce things to a spray
    if(is.spray(arg1)){
        return(arg1)  
    } else if(is.list(arg1)){
        return(spraymaker(arg1),addrepeats=addrepeats)
    } else if(is.matrix(arg1) & !missing(arg2)){
        return(spraymaker(list(arg1,arg2),addrepeats=addrepeats))
    } else if(is.array(arg1)){
        ind <- which(arg1!=0,arr.ind=TRUE)
        if(offbyone){ ind <- ind-1 }
        return(spray(ind, arg1[arg1!=0]))
    } else {
        stop("as.spray() needs a list or an array")
    }  
}

`arity` <- function(S){
    ncol(index(S))
}

setGeneric("dim")
setGeneric("deriv")

`dim.spray` <- function(x){
    ind <- index(x)
    if(any(ind<0)){
        stop("dim() does not make sense for sprays with negative indices")
    } else {
        return(apply(ind,2,max))
    }
}

`as.array.spray` <- function(x, offbyone=FALSE, compact=FALSE, ...){
    if(compact){ind <- sweep(ind,2,apply(ind,2,min)) }
    dS <- dim(x) 
    ind <- index(x)
    if(offbyone) {
      ind <- ind+1
      dS <- dS + 1
    }

    if(any(ind<0)){
        stop("There are negative index elements")
    } else if(any(ind==0)){
        stop("There are zero index elements")
    }
    
    out <- array(0,dS)
    out[ind] <- value(x)
    dimnames(out) <- names(ind)
    return(out)
}


#`spray_missing_accessor` <- function(S,dots){
#  stop("not yet implemented")
#  browser()



  ## Not implemented.
  
  ## the problem is illustrated by removing the stop() statement
  ## above, and changing the first line of the body of `[.spray`() to
  ## something like spray_missing_accessor(S,dots=match.call()).
  
  ## Then loading the package and typing 'S <- spray(diag(1:6)) ;
  ## S[1,,1:3]' on the commandline.  Then you enter the browser;
  ## 'dots' is a 5-element list.  The first two are the function name
  ## and S.  Strip these two out. Then the first element seems to be
  ## just a numeric "1".  The second element, corresponding to the
  ## thing between the two consecutive commas, is some sort of ghostly
  ## not NULL, not a list, nothing that I can make sense of.  For
  ## example, 'dput(dots[[2]])' gives an empty line; 'dput(dots[2])'
  ## gives structure(list(), .Names = ""); but note that
  ## structure(list(), .Names = "") will not parse on the commandline.
  ## Therefore one cannot use constructions such as
  ## 'identical(dots[2], structure(list(), .Names = ""))'.
  ##
  ## This might be surmountable, but for the following: dots[[3]]
  ## prints "1:3" on the commandline.  So typing '1 + dots[[3]]' gives
  ## an error (specifically, "non-numeric argument to binary
  ## operator").

  ## This is why Jonty had to resort to peculiar tricks with
  ## eval.parent() in Oarray.

  ## ... quite apart from all this, there is another reason why I have
  ## not implemented things like S[1,,1:3].  Such constructions simply
  ## do not fit nicely into the STL map methods in spray_ops.cpp.  At
  ## least I can't see any nice way.  You'd need some sort of partial
  ## matching system

#}




`[.spray` <- function(S, ...,drop=FALSE){
    dots <- list(...)
    first <- dots[[1]]
    if(is.matrix(first)){
        M <- first
    } else if(is.spray(first)){
        M <- index(first)
    } else {
        M <- as.matrix(expand.grid(dots))
    }
    if(ncol(M) != arity(S)){
        stop("incorrect number of dimensions specified")
    }
    jj <- spray_accessor(index(S),value(S),M) # returns a NumericVector
    if(drop){
        return(jj)
    } else {
        return(spray(M,jj))
    }
}  

`[<-.spray` <- function(S, index, ..., value){

    if(missing(index)){ # S[] <- something
        return(spray(spray::index(S),value))
    }
    
    if(is.matrix(index)){
        M <- index
    } else if(is.spray(index)){
        M <- spray::index(index)
    } else {
        M <- as.matrix(expand.grid(c(list(index), list(...))))
    }
    if(ncol(M) != arity(S)){
        stop("incorrect number of dimensions specified")
    }

    if(length(value)==1){value <- rep(value,nrow(M))}
    stopifnot(length(value)==nrow(M))
    return(spraymaker(spray_setter(spray::index(S),spray::value(S),M,value)))
}

`as.function.spray` <- function(x,...){
    function(X){
        if(!is.matrix(X)){X  <- rbind(X)}
        stopifnot(ncol(X)==arity(x))
        jj <- matrix(1,nrow(X),nrow(index(x)))
        for(i in seq_len(arity(x))){ # iterate over index columns
            jj <- jj * outer(X[,i],index(x)[,i],"^")
        }
        return(rowSums(sweep(jj,2,value(x),"*")))
    }
}

`constant` <- function(S,drop=FALSE){ # returns 'the constant (term) of S'
  S[t(rep(0,arity(S))),drop=drop]
}

`constant<-` <- function(S,value){ "[<-"(S,t(rep(0,arity(S))),value=value)}

`product` <- function(power){spray(rbind(power))}
`homog` <- function(d,power=1){spray(partitions::compositions(power,d))}
`linear` <- function(x,power=1){spray(diag(power,length(x)),x)}
`lone` <- function(n,d=n){
  jj <- rep(0,d)
  jj[n] <- 1
  return(product(jj))
}

`xyz` <- function(d){product(rep(1,d))}

`one` <- function(d){
  if(is.spray(d)){d <- arity(d)}
  return(spray(t(rep(0,d))))
}

`ooom` <- function(S,n){
  out <- one(S)
  for(i in seq_len(n)){
    out <- 1  + S*out
  }
  return(out)
}

`print_spray_matrixform` <- function(S){
    if(is.empty(S)){
        cat('the NULL polynomial\n')
    } else {
        jj <-
            data.frame(index(S),symbol= " = ", val=round(value(S),getOption("digits")))
        colnames(jj) <- c(rep(" ",arity(S)+1),'val')
        print(jj,row.names=FALSE)
    }
    return(invisible(S))
}

`print_spray_polyform` <- function(S){

    vector_of_3vars <- c("x","y","z")
    multiply_symbol1 <- "*"
    multiply_symbol2 <- "*"

    if(is.empty(S)){
    cat('the NULL polynomial\n')
    return(S)
  }

  if(arity(S) <= 3){
     variables <- vector_of_3vars
    } else {
     variables <- paste("x",seq_len(arity(S)),sep="")
    }

  if(is.empty(constant(S))){
    string <- ""
  } else {
    string <- constant(S,drop=TRUE)
  }
  
  ind <- index(S)
  val <- value(S)

    for(i in seq_len(nrow(ind))){ 
        coeff <- val[i]
        coeff_abs <- abs(coeff)
        coeff_sign <- sign(coeff)
        if(all(ind[i,]==0)){next}  # constant term
        if(coeff_abs == 1){
            term <- ""
        } else {
            term <- paste(round(coeff_abs, getOption("digits")), "",multiply_symbol1, sep="")
        }
        
        if(coeff_sign == 1){
            term <- paste(" +", term, sep="")
        } else {
            term <- paste(" -", term, sep="")
        }

        jj <- ind[i,]
        vars <- variables[which(jj !=0)]
        powers <- jj[jj!=0]
        lv <- length(vars)
        if(lv==0){ # constant; already printed
            ignore <- 3  # do nothing
        } else if(lv==1){   # just one variable
            if(powers==1){
                term <- paste(term, vars, sep="") 
            } else {
                term <- paste(term, vars, "^", powers, sep="")
            }
        } else {   ## length(vars) > 1
            for(i in seq_len(lv-1)){
                if(powers[i] == 1){
                    term <- paste(term, vars[i], multiply_symbol2,sep="")
                } else {
                    term <- paste(term, vars[i],"^",powers[i],multiply_symbol2,sep="") 
                }
            }
            
            if(powers[lv] == 1){
                term <- paste(term, vars[lv],sep="") 
            } else {
                term <- paste(term, vars[lv],"^",powers[lv],sep="") 
            }
        }
        string <- paste(string, term, sep="")

      }  # row of index loop closes
    
    string <- paste(string,'\n',sep="")
    cat(string)
    return(invisible(S))
}
  
`print.spray` <- function(x, ...){
  if(isTRUE(getOption("polyform",default=FALSE))){
    out <- print_spray_polyform(x)
  } else {
    out <- print_spray_matrixform(x)
  }
  return(out)
}

`asum` <- function(S, dims,drop=TRUE, ...) {UseMethod("asum")}

#      M <- matrix(sample(0:2,5*29,rep=T),ncol=5)

`process_dimensions` <- function(S,dims){
  if(is.logical(dims)){ dims <- which(dims)  }
  stopifnot(all(dims <= arity(S)))
  stopifnot(all(dims >0))
  stopifnot(all(dims == round(dims)))
  return(dims)
}

`asum.spray` <- function(S, dims, drop=TRUE, ...){   # off-by-one dealt with in C++
  dims <- process_dimensions(S,dims)
  jj <- spraymaker(spray_asum_include(index(S),value(S),dims))
  if(drop){
    return(spray(index(jj)[,-dims,drop=FALSE],value(jj),addrepeats=TRUE))
  } else {
    return(jj)
  }
}

`asum_inverted` <- function(S, dims){   # off-by-one dealt with in C++
  dims <- process_dimensions(S,dims)
  return(spraymaker(spray_asum_exclude(index(S),value(S),dims)))
}

`subs` <- function(S,dims,x){
  dims <- process_dimensions(S,dims)
  return(spray(index(S)[,-dims,drop=FALSE], drop(index(S)[,dims,drop=FALSE]%*%x), addrepeats=TRUE))
}

`deriv.spray` <- function(expr, i, derivative=1, ...){
  S <- as.spray(expr)
  orders <- rep(0,arity(S))
  orders[i] <- derivative
  return(aderiv(S,orders))
}

`aderiv` <- function(S,orders){
  orders <- round(orders)
  stopifnot(length(orders) == arity(S))
  stopifnot(all(orders >= 0))
  return(spraymaker(spray_deriv(index(S),value(S),orders)))
}  
