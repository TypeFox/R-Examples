##################################################
## This tree reads a list
offspring <- function(path=character(0), lst, ...) {
  if(length(path))
    obj <- lst[[path]]
    else
      obj <- lst
  nms <- names(obj)
  hasOffspring <- sapply(nms, function(i) {
    newobj <- obj[[i]]
    is.recursive(newobj) && !is.null(names(newobj))
    })
  data.frame(comps=nms, hasOffspring=hasOffspring, ## fred=nms,
             stringsAsFactors=FALSE)
}
l <- list(a="1", b= list(a="21", b="22", c=list(a="231")))

w <- gwindow("Tree test")
t <- gtree(offspring=offspring, offspring.data=l, cont=w)


##################################################
## This tree looks at recursive objects
describe <- function(x) UseMethod("describe")
describe.default <- function(x) sprintf("An object with class %s", class(x)[1])
describe.integer <- function(x) sprintf("An integer with %s value%s", length(x), ifelse(length(x) > 1, "s", ""))
describe.numeric <- function(x) sprintf("A numeric with %s value%s", length(x), ifelse(length(x) > 1, "s", ""))
describe.factor <- function(x) sprint("A factor with %s level%s", length(levels(x)), ifelse(length(levels(x)) > 1, "s", ""))

offspring <- function(path, obj) {
  if(length(path) > 0)
    x <- obj[[path]]
  else
    x <- obj

  nms <- names(x)
  recursive <- sapply(x, function(i) {
    is.recursive(i) &&
    !is.null(attr(i, "names")) &&
    length(i) > 0
    })
  descr <- sapply(x, describe)
  
  data.frame(Variable=nms, offspring=recursive, Description=descr, stringsAsFactors=FALSE)
}

l <- lm(mpg ~ wt, mtcars)

w <- gwindow("test")
gtree(offspring=offspring, offspring.data=l, cont=w)


##################################################


offspring <- function(path, ...) {
  data.frame(key=letters[1:2], offspring=c(TRUE, TRUE), icon=c("ok", "quit"), descr=paste("Level ", length(path), "letter", letters[1:2], sep=" "),
             stringsAsFactors=FALSE)
}

w <- gwindow("Testing")
tr <- gtree(offspring=offspring, icon.col=3, cont=w)



## svalue
ind <- c(1,2,1,2)
svalue(tr, index=TRUE) <- ind
expect_equal(svalue(tr, index=TRUE), ind)

## svalue<- index=TRUE

## multiple

## 
