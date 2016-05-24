#------------------------------------------------------------------------------
# .format.perc : for formatting percentages
#------------------------------------------------------------------------------
### code borrowed from non-exported code from confint.default from package stats
.format.perc <- function (probs, digits)
    paste(format(100 * probs, trim = TRUE, scientific = FALSE, digits = digits),
    "%")

#------------------------------------------------------------------------------
### for distrXXX pre 2.5
#------------------------------------------------------------------------------


if(packageVersion("distr")<"2.5"){

.inArgs <- function(arg, fct)
          {as.character(arg) %in% names(formals(fct))}

.fillList <- function(list0, len = length(list0)){
            if(len == length(list0))
               return(list0)
            i <- 0
            ll0 <- length(list0)
            li0 <- vector("list",len)
            if(ll0)
            while(i < len){
               j <- 1 + ( i %% ll0)
               i <- i + 1
               li0[[i]] <- list0[[j]]
            }
           return(li0)
}

.ULC.cast <- function(x){
         if( is(x,"AbscontDistribution"))
             x <- as(as(x,"AbscontDistribution"), "UnivarLebDecDistribution")
         if(is(x,"DiscreteDistribution"))
             x <- as(as(x,"DiscreteDistribution"), "UnivarLebDecDistribution")
         if(!is(x,"UnivarLebDecDistribution"))
            x <- as(x,"UnivarLebDecDistribution")
         return(x)
}

.isEqual <- function(p0, p1, tol = min( getdistrOption("TruncQuantile")/2,
                                          .Machine$double.eps^.7
                                          ))
                abs(p0-p1)< tol

.isIn <- function(p0, pmat, tol = min( getdistrOption("TruncQuantile")/2,
                                          .Machine$double.eps^.7
                                          ))
                  {list1 <- lapply(1:nrow(pmat), function(x){
                            (p0+tol > pmat[x,1]) & (p0-tol < pmat[x,2]) })
                   apply(matrix(unlist(list1), ncol = nrow(pmat)), 1, any)}


.isEqual01<- function(x) .isEqual(x,0)|.isEqual(x,1)

.presubs <- function(inp, frompat, topat){
### replaces in an expression or a string all frompat patterns to topat patterns

logic <- FALSE
inCx <- sapply(inp,
   function(inpx){
      inC <- deparse(inpx)
      l <- length(frompat)
      for(i in 1:l)
         { if (is.language(topat[[i]])){
               totxt <- deparse(topat[[i]])
               totxt <- gsub("expression\\(", "\", ", gsub("\\)$",", \"",totxt))
               if (length(grep(frompat[i],inC))) logic <<- TRUE
               inC <- gsub(frompat[i],totxt,inC)
           }else inC <- gsub(frompat[i], topat[[i]], inC)
         }
      return(inC)
    })
if(length(grep("expression",inCx))>0)
   inCx <- gsub("expression\\(", "", gsub("\\)$","",inCx))
if (length(inCx) > 1) {
   inCx <- paste(inCx, c(rep(",", length(inCx)-1), ""),
                 sep = "", collapse = "\"\\n\",")
   if ( any(as.logical(c(lapply(inp,is.language)))) | logic )
      inCx <- paste("expression(paste(", gsub("\\\\n"," ", inCx), "))", sep ="")
   else
      inCx <- paste("paste(",inCx,")", sep ="")
}else inCx <- paste("expression(paste(",inCx,"))",sep="")
outC <- eval(parse(text = eval(inCx)))
return(outC)
}

.DistrCollapse <- function(support, prob,
                              eps = getdistrOption("DistrResolution")){
    supp <- support
    prob <- as.vector(prob)
    suppIncr <- diff(c(supp[1]-2*eps,supp)) < eps
    groups <- cumsum(!suppIncr)
    prob <- as.vector(tapply(prob, groups, sum))
    supp <- as.vector(tapply(supp, groups, quantile, probs = 0.5, type = 1))
           ### in order to get a "support member" take the leftmost median
    return(list(supp = supp, prob = prob))
#    newDistribution <- DiscreteDistribution(supp=supp,prob=prob)
#    return(newDistribution)
}

.makeLenAndOrder <- function(x,ord){
   n <- length(ord)
   x <- rep(x, length.out=n)
   x[ord]
}

}

if(packageVersion("distrMod")<"2.5"){
.isUnitMatrix <- function(m){
### checks whether m is unit matrix
              m.row <- nrow(m)
              isTRUE(all.equal(m, diag(m.row), check.attributes = FALSE))
              }

.deleteDim <- function(x){
     attribs <- attributes(x)
     attribs$dim <- NULL
     attribs$dimnames <- NULL
     attributes(x) <- attribs
     x
     }

}


