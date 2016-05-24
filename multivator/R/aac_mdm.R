setClass("mdm", # "mdm" == "multivariate design matrix"
         representation = representation(
           xold  = "matrix",
           types = "factor"
           )
         )

"mdm" <- function(xold,type){new("mdm", xold=xold, types=types)}

setGeneric("xold", function(x){standardGeneric("xold")})
setGeneric("xold<-", function(x,value){standardGeneric("xold<-")})
setGeneric("levels",function(x){standardGeneric("levels")})
### types() and names() already generic (in mhp.R)


setMethod("xold","mdm",function(x){x@xold})
setMethod("types","mdm",function(x){x@types})
### No occurrences of "@" below this line.


setMethod("types<-","mdm",function(x,value){
  mdm(xold=xold(x),types=value)
} )


setMethod("xold<-","mdm",function(x,value){
  jj <- xold(x)
  jj[] <- value
  return(mdm(jj,types(x)))
} )


setMethod("levels","mdm",function(x){levels(types(x))})
setMethod("levels<-","mdm",function(x,value){
  jj <- types(x)
  levels(jj) <- value
  return(mdm(xold=xold(x),types=jj))
} )


setMethod("names","mdm",function(x){colnames(xold(x))})

setMethod("names<-","mdm",function(x,value){
  jj <- xold(x)
  colnames(jj) <- value
  return(mdm(jj,types(x)))
})


".mdm_valid" <- function(object){
  if(nrow(xold(object)) !=  length(types(object))){
    return("length of types not equal to nrow(xold)")
  } else {
    return(TRUE)
  }
}

setValidity("mdm", .mdm_valid)

"mdm" <- function(xold, types){ 
  new("mdm" , xold=xold,types=as.factor(types))
# This is the only place new("mdm", ...) is called
}

"as.mdm" <- function(x,  ...){
    n <- ncol(x)
    mdm(x[,-n], types=as.factor(x[,n,drop=TRUE]))
}

is.mdm <- function(x){inherits(x,"mdm")}

setAs("mdm","matrix", function(from){
  cbind(xold(from),type=as.numeric(types(from)))
} )

setGeneric("as.matrix")
setMethod("as.matrix","mdm", function(x){as(x,"matrix")})

setAs("mdm","list",function(from){
  jj <- by(xold(from),types(from), function(x){x})
  attributes(jj)$class <- NULL
  attributes(jj)$call <- NULL
  attributes(jj)$dimnames <- list(levels(from))
  return(jj)
 } )

setGeneric("as.list")
setMethod("as.list","mdm", function(x){as(x,"list")})

".mdm_print" <- function(x, ...){
  data.frame(xold(x),type=types(x))
}
    
"print.mdm" <- function(x, ...){
  jj <- .mdm_print(x, ...)
  print(jj)
  return(invisible(jj))
}

setMethod("show", "mdm", function(object){print.mdm(object)})

setAs("mdm","data.frame",function(from){
  data.frame(xold(from),type=types(from))
})

setGeneric("as.data.frame")
setMethod("as.data.frame",signature=c("mdm","missing","missing"),function(x,row.names=NULL,optional=TRUE, ...){as(x,"data.frame")})

setGeneric("rownames")
setMethod("rownames","mdm",function(x, do.NULL=TRUE,prefix="row"){rownames(xold(x))})

setGeneric("rownames<-")
setMethod("rownames<-","mdm",function(x, value){
  jj <- xold(x)
  rownames(jj) <- value
  return(mdm(jj,types(x)))
} )

setMethod("[",signature(x="mdm"),
          function(x,i,j,drop=FALSE){
            if(missing(i)){
              return(xold(x)[,j,drop=drop])
            }
            if(missing(j)){
              j <- TRUE
            }
            return(mdm(xold=xold(x)[i,j,drop=drop],types=types(x)[i,drop=FALSE]))
          } )
          
setReplaceMethod("[",signature(x="mdm"),
                 function(x,i,j,value){
                   jj <- xold(x)
                   if(missing(j)){
                     jj[i,] <- value
                   } else {
                     jj[i,j] <- value
                   }
                   return(mdm(xold=jj, types=types(x)))
                   } )


setGeneric("nrow")
setGeneric("ncol")
setGeneric("dim" )

setMethod("nrow",signature=c("mdm"),function(x){nrow(xold(x))})
setMethod("ncol",signature=c("mdm"),function(x){ncol(xold(x))})
setMethod("dim" ,signature=c("mdm"),function(x){ dim(xold(x))})


setAs("mdm","mhp",function(from){
  levs <- levels(from)
  nams <- names(from)
  M <- diag(nrow=length(levs))
  B <- array(rep(diag(nrow=length(nams)),length(levs)),c(length(nams),length(nams),length(levs)))
  return(mhp(M,B,levels=levs,names=nams))
} )

setGeneric("as.mhp", function(x){standardGeneric("as.mhp")})
setMethod("as.mhp","mdm",function(x){as(x,"mhp")})

setMethod("head",signature="mdm",function(x,n=6,...){
  mdm(head(xold(x),n=n,...) , types=factor(head(types(x),n=n,...)))
} )

setMethod("tail",signature="mdm",function(x,n=6,...){
  mdm(tail(xold(x),n=n,...) , types=factor(tail(types(x),n=n,...)))
} )

# following lines copied from the Brobdingnag c() and .cPair() functions:
setGeneric(".mdm_rbind_pair", function(x,y,deparse.level){standardGeneric(".mdm_rbind_pair")})
setMethod (".mdm_rbind_pair", c("mdm", "mdm"), function(x,y,deparse.level){.mdm_rbind(x,y,deparse.level)})
setMethod (".mdm_rbind_pair", c("mdm", "ANY"), function(x,y,deparse.level){.mdm_rbind_error(x,y)})
setMethod (".mdm_rbind_pair", c("ANY", "mdm"), function(x,y,deparse.level){.mdm_rbind_error(x,y)})
setMethod (".mdm_rbind_pair", c("ANY", "ANY"), function(x,y,deparse.level){.mdm_rbind_error(x,y)})

".mdm_rbind_error" <- function(x,y){
  stop("an mdm object may only be rbinded to another mdm object")
}

".mdm_rbind" <- function(x,y,deparse.level){
  stopifnot(compatible(x,y))
  jj <- as.factor(c(types(x),types(y)))
  levels(jj) <- levels(x)     # identical to levels(y) ... as per compatible() above
  mdm(rbind(xold(x),xold(y),deparse.level=deparse.level),jj)
}

# Thanks to Martin Morgan for supplying the following setGeneric():
setGeneric("rbind",
    function(..., deparse.level=1) standardGeneric("rbind"),
    signature = "...")

setMethod("rbind", signature="mdm", function(x, ..., deparse.level=1) {
  if(nargs()<4)
    .mdm_rbind_pair(x, ..., deparse.level=deparse.level)
  else
    .mdm_rbind_pair(x, Recall(..., deparse.level=deparse.level))
})
          
