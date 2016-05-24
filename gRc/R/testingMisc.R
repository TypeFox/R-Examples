
.addStat <- function(tab,n=NULL,direction="drop"){
  tab$p   <- with(tab, 1-pchisq(X2, df))  
  switch(direction,
         "drop"={
           tab$aic <- with(tab, -X2+2*df)
           if (!is.null(n))
             tab$bic <- with(tab, -X2+log(n)*df)  
         },
         "add"={
           tab$aic <- with(tab, X2-2*df)
           if (!is.null(n))
             tab$bic <- with(tab, X2-log(n)*df)  
           

         })
  tab
}


.nestedList2df <- function(nlist){
  nc    <- unique(sapply(nlist, length))
  x     <- lapply(1:nc, function(i) sapply(nlist, "[[", i))
  x2    <- as.data.frame(x)
  names(x2) <- paste("V", 1:nc, sep='')
  x2
}


.findWald <- function(i2,j2, b,V){
  wald <- as.numeric((b[i2] - b[j2])^2 / (V[i2,i2]+V[j2,j2]-2*V[i2,j2]))
  wald 
}
