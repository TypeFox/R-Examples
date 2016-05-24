
.matrix2dataFrame2 <- function (x, at, restoreAll = TRUE){

    d <- dimnames(x)
    k <- length(d[[2]])
    w <- vector("list", k)
    names(w) <- d[[2]]
    
    for (i in 1:k) {
      a   <- at[[i]]
      xi  <- x[,i]
      if (a$sm=='character'){
        storage.mode(xi) <- 'integer'
        attributes(xi) <- a$atr
        xi  <- as.character(xi)
      } else {
        storage.mode(xi) <- a$sm
        attributes(xi)   <- a$atr
      }
      w[[i]] <- xi
    }
    
    val <-structure(w, class = 'data.frame', row.names = d[[1]])
    return(val)

}


.subsAttr2<-function (x){

  at<-lapply(x, function(y){
    ## print("--------------")
    sm  <- storage.mode(y)
    if (sm=='character')
      a   <- attributes(factor(y))
    else
      a   <- attributes(y)
    v <- list(sm=sm,atr=a)
    v
  })
  return(at)
}



.asNumericMatrix2 <- function (x){
  a <- attributes(x)
  k <- length(a$names)
  
  idx <- which(lapply(x,class)=="character")
  if (length(idx)>0){
    for (j in idx){
      x[,j] <- as.factor(x[,j])
    }
  }
  
  y <- matrix(unlist(x), ncol = k, dimnames = list(a$row.names,  a$names))
  y
}




























# .matrix2dataFrame2 <- function (x, at, restoreAll = TRUE)
# {
#   at <<- at
#   d <- dimnames(x)
#   k <- length(d[[2]])
#   w <- vector("list", k)
#   nam <<- names(w) <- d[[2]]
#   sm <- storage.mode(x) ##;print(sm)

#   for (i in 1:k) {
#     a <- at[[nam[i]]];
#     cat("name:", nam[[i]],"\n");    print(a)
#     if (!length(a))
#       next
#     xi <- x[, i]
#     names(xi) <- NULL
    
#     ## SHD Add:
#     ## Handle characters
#     ## OLD
#     # if (!is.null(a$storage.mode) && a$storage.mode == "character"){
# #       xi <- a$flevels[xi]
# #       a$flevels <- NULL
# #     }

#     ## NEW
#     print(storage.mode(a))
#     if (!is.null(storage.mode(a)) && storage.mode(a) == "character"){
#       print("char")
#       print(a$flevels)
#       xi <- a$flevels[xi]
#       a$flevels <- NULL
#     }

    
#     ## Handle times
#     ## OLD
# #     if (!is.null(a$class)){
# #       if (identical(a$class, c("POSIXt", "POSIXct"))){
# #         if (a$storage.mode != sm)
# #           storage.mode(xi) <- a$storage.mode
# #         ##a$storage.mode <- NULL
# #         attributes(xi) <- a       
# #       }          
# #     }

#     ## NEW
#     if (!is.null(class(a))){
#       if (identical(class(a), c("POSIXt", "POSIXct"))){
#         if (storage.mode(a) != sm)
#           storage.mode(xi) <- storage.mode(a)
#         ##a$storage.mode <- NULL
#         attributes(xi) <- a       
#       }          
#     }

    
#     ## End !
    
#     if (restoreAll) {
#       ## OLD
#       ##if (a$storage.mode != sm)
#       ##  storage.mode(xi) <- a$storage.mode
#       ##a$storage.mode <- NULL

#       ## NEW
#       if (storage.mode(a) != sm)
#         storage.mode(xi) <- storage.mode(a)
#       storage.mode(a) <- NULL

#       attributes(xi) <- a
#     }
#     else {
#       a<<-a
#       ##print("Here")
#       ##if (length(l <- a$label)){ ## OLD
#       if (length(l <- label(a))){  ## NEW
#         ## print("---label"); print(l)
#         label(xi) <- l
#       }
#       ##if (length(u <- a$units)){  ## OLD
#       if (length(u <- units(a))){   ## NEW
#         ##print("---units"); print(u)
#         units(xi) <- u
#       }
#       ##if (length(lev <- a$levels)){ ## OLD
#       if (length(lev <- levels(a))){   ## NEW
#         ##print(xi);
#         ##print(a$levels)
#         ##xi <<- xi
#         ##lev <<- lev
#         xi <- factor(xi, 1:length(lev), lev)
#         ##print(xi)
#       }
#     }
#     w[[i]] <- xi
#   }
#   structure(w, class = "data.frame", row.names = d[[1]])
# }


# # .subsAttr2<-function (x)
# # {
# #     g <- function(y) {
# #       a <- attributes(y)
# #       print("in g")
# #       print(storage.mode(y))
# #       print(a)
# #       a$dim <- a$names <- a$dimnames <- a$flevels <- NULL
# #       ##storage.mode(a) <- storage.mode(y)
# #       a$storage.mode <- storage.mode(y)   ## IKKE ændret ....
# #       ## Add:
# #       if (length(class(y))==1 && class(y)=="character"){
# #         print("a char...")
# #         a$flevels <- levels(as.factor(y))   ## NEW
# #         ##a$flevels <- levels(as.factor(y)) ## OLD
# #       }
# #       ## End!
# #       a
# #     }
# #     if (is.list(x)){
# #       ##cat("calling 1\n")
# #       sapply(x, g)
# #     } else {
# #       ##cat("calling 2\n")
# #       g(x)
# #     }
# #   }


# .subsAttr2<-function (x)
# {
#     g <- function(y) {
#       a <- attributes(y)
#       print("in g")
#       print(storage.mode(y))
#       print(a)
#       a$dim <- a$names <- a$dimnames <- a$flevels <- NULL
#       ##storage.mode(a) <- storage.mode(y)
#       a$storage.mode <- storage.mode(y)   ## IKKE ændret ....
#       ## Add:
#       if (length(class(y))==1 && class(y)=="character"){
#         print("a char...")
#         a$flevels <- levels(as.factor(y))   ## NEW
#         ##a$flevels <- levels(as.factor(y)) ## OLD
#       }
#       ## End!
#       a
#     }
#     if (is.list(x)){
#       ##cat("calling 1\n")
#       sapply(x, g)
#     } else {
#       ##cat("calling 2\n")
#       g(x)
#     }
#   }


# .asNumericMatrix2 <- function (x)
# {
#     a <- attributes(x)
#     k <- length(a$names)

#     idx <- which(lapply(x,class)=="character")
#     if (length(idx)>0){
#       for (j in idx){
#         x[,j] <- as.factor(x[,j])
#       }
#     }
    

# #     val <- lapply(x, function(xx){
# #       if (class(xx)=="character")
# #         as.factor(xx)
# #       else
# #         xx
# #     })
# #     x<- val


#     y <- matrix(unlist(x), ncol = k, dimnames = list(a$row.names,  a$names))
#     #y <- matrix(unlist(val), ncol = k, dimnames = list(a$row.names,  a$names))
#     ## End

#     # y <- matrix(as.numeric(unlist(x)), ncol = k, dimnames = list(a$row.names,  a$names))
#     #if (storage.mode(y) == "character")
#     #    warning("x had at least one character vector")
#     y
# }



### SHD re-implementation
