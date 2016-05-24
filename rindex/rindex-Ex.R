pkgname <- "rindex"
source(file.path(R.home("share"), "R", "examples-header.R"))
options(warn = 1)
library('rindex')

assign(".oldSearch", search(), pos = 'CheckExEnv')
cleanEx()
nameEx("index")
### * index

flush(stderr()); flush(stdout())

### Name: index
### Title: Indexing vectors
### Aliases: index c.index indexAddTree indexDelTree indexAutobatch rindex
###   c.rindex rindexAddTree rindexDelTree rindexAutobatch
### Keywords: misc database

### ** Examples

  #library(rindex)
  x <- sample(c(rep("c", 5), paste("d", 1:5, sep=""), c(letters[c(1,2,5,6)], NA)))

  cat("\n")
  cat("creating an index from atomic vector\n")
  i <- index(x)
  i
  cat("creating an index by combining indices\n")
  i2 <- c(i,i)
  i2
  cat("if the index (or index$tree) is removed, the C-tree is removed at the next garbage collection\n")
  cat("the index tree can also removed and created explicitely\n")
  i <- indexDelTree(i)
  i
  i <- indexAddTree(i, batch=3)
  print(i, tree=TRUE)
  indexNodes(i)
  indexBytes(i)

  cat("\n")
  cat("extracting the original vector\n")
  i[]
  cat("subsetting works as expected\n")
  i[1:3]
  cat("accessing the sorted data is much faster\n")
  sort(i)[1:3]
  cat("accessing the ordering is also much faster (order.index is not dispatched since order is not yet generic)\n")
  order.index(i)[1:3]

  identical(is.na(i),is.na(x))
  identical(length(i),length(x))

  cat("\n")
  cat("LOW LEVEL SEARCH returns position in SORTED VECTOR\n")
  cat("low level search for position of lowest instance of value in the index\n")
  indexFind(i, "c")
  cat("low level search for position of highest instance of value in the index\n")
  indexFind(i, "c", findlow=FALSE)
  cat("low level search for position of lowest instance beginning like search value\n")
  indexFindlike(i,"d")
  cat("low level search for position of highest instance beginning like search value\n")
  indexFindlike(i,"d",findlow=FALSE)

  cat("\n")
  cat("MID LEVEL SEARCH also returns position in SORTED VECTOR\n")
  cat("mid level search for a set of values\n")
  indexMatch(i,c("c","f"), findlow=TRUE)  # giving parameter findlow= suppresses the warning issued on non-unique indices
  sort(i)[indexMatch(i,c("c","f"), findlow=TRUE)]
  i[indexMatch(i,c("c","f"), findlow=TRUE, what="pos")]
  indexMatch(i,c("c","f"), findlow=TRUE, what="val")
  indexMatch(i,c("c","f"), findlow=TRUE, what="pos")
  indexMatch(i,c("c","f"), findlow=FALSE, what="pos")

  cat("mid level search for interval of values\n")
  indexFindInterval(i,"b","f")
  cat("by default the searched endpoints are included\n")
  sort(i)[indexFindInterval(i,"b","f")]
  cat("but they can be excluded\n")
  sort(i)[indexFindInterval(i,"b","f",high.include=FALSE)]
  cat("by default the searched endpoints need not to be present\n")
  sort(i)[indexFindInterval(i,"a1","e1")]
  cat("but this can be required\n")
  sort(i)[indexFindInterval(i,"a1","e1",low.exact=TRUE)]
  cat("each of the searched endpoints can be defined via indexFindlike\n")
  sort(i)[indexFindInterval(i,"c","d",FUN=indexFindlike)]


  cat("\n")
  cat("HIGH LEVEL SEARCH returns POSITION(s) IN ORIGINAL VECTOR but in SEQUENCE OF INDEX\n")
  indexEQ(i,"d3")
  indexNE(i,"d3")
  indexLT(i,"d3")
  indexLE(i,"d3")
  indexGT(i,"d3")
  indexGE(i,"d3")
  cat("searching for several values returns a list\n")
  indexEQ(i,c("b","c","z",NA))
  indexEQ(i,c("b","c","z",NA), what="val")


  cat("\n")
  cat("HIGH LEVEL OPERATORS returns TRUE/FALSE AT ORIGINAL VECTOR POSITIONS\n")
  i=="d3"
  i!="d3"
  i<"d3"
  i<="d3"
  i>"d3"
  i>="d3"

  cat("HIGH LEVEL match.index and  behave as expected\n")
  match(c("b","c","z",NA), x)
  match.index(c("b","c","z",NA), i)

## Don't show: 

  cat("\n")
  cat("creating an rindex from atomic vector\n")
  ri <- rindex(x)
  ri
  cat("creating an rindex by combining indices\n")
  ri2 <- c(ri,ri)
  ri2
  cat("if the rindex (or rindex$tree) is removed, the C-tree is removed at the next garbage collection\n")
  cat("the rindex tree can also removed and created explicitely\n")
  ri <- rindexDelTree(ri)
  ri
  ri <- rindexAddTree(ri, batch=3)
  print(ri, tree=TRUE)
  rindexNodes(ri)
  rindexBytes(ri)

  cat("\n")
  cat("extracting the original vector\n")
  ri[]
  cat("subsetting works as expected\n")
  ri[1:3]
  cat("accessing the sorted data is much faster\n")
  sort(ri)[1:3]
  cat("accessing the ordering is also much faster (order.rindex is not dispatched since order is not yet generic)\n")
  order.rindex(ri)[1:3]

  identical(is.na(ri),is.na(x))
  identical(length(ri),length(x))

  cat("\n")
  cat("LOW LEVEL SEARCH returns position in SORTED VECTOR\n")
  cat("low level search for position of lowest instance of value in the rindex\n")
  rindexFind(ri, "c")
  cat("low level search for position of highest instance of value in the rindex\n")
  rindexFind(ri, "c", findlow=FALSE)
  cat("low level search for position of lowest instance beginning like search value\n")
  rindexFindlike(ri,"d")
  cat("low level search for position of highest instance beginning like search value\n")
  rindexFindlike(ri,"d",findlow=FALSE)

  cat("\n")
  cat("MID LEVEL SEARCH also returns position in SORTED VECTOR\n")
  cat("mid level search for a set of values\n")
  rindexMatch(ri,c("c","f"), findlow=TRUE)  # giving parameter findlow= suppresses the warning issued on non-unique indices
  sort(ri)[rindexMatch(ri,c("c","f"), findlow=TRUE)]
  ri[rindexMatch(ri,c("c","f"), findlow=TRUE, what="pos")]
  rindexMatch(ri,c("c","f"), findlow=TRUE, what="val")
  rindexMatch(ri,c("c","f"), findlow=TRUE, what="pos")
  rindexMatch(ri,c("c","f"), findlow=FALSE, what="pos")

  cat("mid level search for interval of values\n")
  rindexFindInterval(ri,"b","f")
  cat("by default the searched endpoints are included\n")
  sort(ri)[rindexFindInterval(ri,"b","f")]
  cat("but they can be excluded\n")
  sort(ri)[rindexFindInterval(ri,"b","f",high.include=FALSE)]
  cat("by default the searched endpoints need not to be present\n")
  sort(ri)[rindexFindInterval(ri,"a1","e1")]
  cat("but this can be required\n")
  sort(ri)[rindexFindInterval(ri,"a1","e1",low.exact=TRUE)]
  cat("each of the searched endpoints can be defined via rindexFindlike\n")
  sort(ri)[rindexFindInterval(ri,"c","d",FUN=rindexFindlike)]

  cat("\n")
  cat("HIGH LEVEL SEARCH returns POSITION(s) IN ORIGINAL VECTOR but in SEQUENCE OF INDEX\n")
  rindexEQ(ri,"d3")
  rindexNE(ri,"d3")
  rindexLT(ri,"d3")
  rindexLE(ri,"d3")
  rindexGT(ri,"d3")
  rindexGE(ri,"d3")
  cat("searching for several values returns a list\n")
  rindexEQ(ri,c("b","c","z",NA))
  rindexEQ(ri,c("b","c","z",NA), what="val")

  cat("\n")
  cat("HIGH LEVEL OPERATORS returns TRUE/FALSE AT ORIGINAL VECTOR POSITIONS\n")
  ri=="d3"
  ri!="d3"
  ri<"d3"
  ri<="d3"
  ri>"d3"
  ri>="d3"

  cat("HIGH LEVEL match AND  behave as expected\n")
  match(c("b","c","z",NA), x)
  match.rindex(c("b","c","z",NA), ri)
## End Don't show

## Not run: 
##D 
##D    cat("function timefactor helps with timing\n")
##D 
##D    n <-  1000000
##D    x <- sample(1:n)
##D    names(x) <- paste("a", x, sep="")
##D    d <- data.frame(x=as.vector(x), row.names=names(x))
##D 
##D    nsub  <- 100
##D    i <- sample(1:n, nsub)
##D    ni <- names(x)[i]
##D 
##D    ind <- index(names(x), verbose=TRUE)
##D    ind
##D 
##D    # test vectors
##D    cat("character susetting is by magnitude slower than integer subsettting\n")
##D    timefactor( x[ni] , x[i] , 10, 10000)
##D    cat("character susetting is approx as slow as matching\n")
##D    timefactor( x[ni] , x[match(ni,names(x))] , 10, 10)
##D    cat("for small fractions of n indexing is much faster\n")
##D    timefactor( x[ni] , x[indexMatch(ind,ni)] , 10, 100)
##D 
##D    # test dataframes
##D    cat("character susetting is by magnitude slower than integer subsettting\n")
##D    timefactor( d[ni,] , d[i,] , 1, 100)
##D    cat("obvious implementation problem (in R-2.5.1 subsetting is much slower than via matching)\n")
##D    timefactor( d[ni,] , d[match(ni,rownames(d)),] , 1, 1)
##D    cat("for small fractions of n indexing is much faster\n")
##D    timefactor( d[ni,] , d[indexMatch(ind,ni),] , 1, 10)
##D 
## End(Not run)

## Don't show: 

  #library(rindex)
  x <- c(rep("c", 5), paste("d", 1:5, sep=""), c(letters[c(1,2,5,6)], NA))[
   c(6L, 14L, 12L, 8L, 1L, 11L, 5L, 4L, 10L, 3L, 2L, 15L, 7L, 9L, 13L) ]

  i <- index(x)
  ri <- rindex(x)
  stopifnot(identical(indexNodes(i), rindexNodes(ri)))
  i$tree <- NULL
  ri$tree <- NULL
  stopifnot(identical(unclass(i),unclass(ri)))

  i <- index(x)
  i2 <- c(i,i)
  ri <- rindex(x)
  ri2 <- c(ri,ri)
  stopifnot(identical(indexNodes(i2), rindexNodes(ri2)))
  i2$tree <- NULL
  ri2$tree <- NULL
  stopifnot(identical(unclass(i2),unclass(ri2)))

  stopifnot(identical(length(i),length(x)))
  stopifnot(identical(length(ri),length(x)))
  stopifnot(identical(i[],x))
  stopifnot(identical(ri[],x))
  stopifnot(identical(is.na(i),is.na(x)))
  stopifnot(identical(is.na(ri),is.na(x)))

  success <- TRUE

  success <- success && binregtest(
    sort
  , sort
  , PAR1=list(i)
  , PAR2=list(x)
  , na.last=list(missing, FALSE, TRUE)
  , decreasing=list(missing, FALSE, TRUE)
  , COMP=identical
  , NAME="COMPARE sort BETWEEN index AND original vector"
  )

  success <- success && binregtest(
    sort
  , sort
  , PAR1=list(ri)
  , PAR2=list(x)
  , na.last=list(missing,FALSE, TRUE)
  , decreasing=list(missing,FALSE,TRUE)
  , COMP=identical
  , NAME="COMPARE sort BETWEEN rindex AND original vector"
  )

  success <- success && binregtest(
    order.index
  , order
  , PAR1=list(i)
  , PAR2=list(x)
  , na.last=list(missing,FALSE, TRUE)
  , decreasing=list(missing,FALSE) # ,TRUE
  , COMP=identical
  , NAME="COMPARE order BETWEEN non-unique index AND original vector"
  )

  success <- success && binregtest(
   order.rindex
  , order
  , PAR1=list(ri)
  , PAR2=list(x)
  , na.last=list(missing,FALSE, TRUE)
  , decreasing=list(missing,FALSE) # ,TRUE
  , COMP=identical
  , NAME="COMPARE order BETWEEN non-unique rindex AND original vector"
  )

  success <- success && binregtest(
    order.index
  , order
  , PAR1=list(index(unique(x)))
  , PAR2=list(unique(x))
  , na.last=list(missing,FALSE, TRUE)
  , decreasing=list(missing,FALSE,TRUE)
  , COMP=identical
  , NAME="COMPARE order BETWEEN unique index AND original vector"
  )

  success <- success && binregtest(
    order.rindex
  , order
  , PAR1=list(rindex(unique(x)))
  , PAR2=list(unique(x))
  , na.last=list(missing,FALSE, TRUE)
  , decreasing=list(missing,FALSE,TRUE)
  , COMP=identical
  , NAME="COMPARE order BETWEEN unique rindex AND original vector"
  )

  stopifnot(identical(indexFind(i, "c"),as.integer(c(0,3))))
  stopifnot(identical(rindexFind(ri, "c"),as.integer(c(0,3))))
  stopifnot(identical(indexFind(i, "c", findlow=FALSE),as.integer(c(0,7))))
  stopifnot(identical(rindexFind(ri, "c", findlow=FALSE),as.integer(c(0,7))))

  stopifnot(identical(indexFindlike(i,"d"), as.integer(c(0,8))))
  stopifnot(identical(rindexFindlike(ri,"d"), as.integer(c(0,8))))
  stopifnot(identical(indexFindlike(i,"d",findlow=FALSE), as.integer(c(0,12))))
  stopifnot(identical(rindexFindlike(ri,"d",findlow=FALSE), as.integer(c(0,12))))

  stopifnot(identical(indexMatch(i,c("c","f"), findlow=TRUE), as.integer(c(3,14))))
  stopifnot(identical(sort(i)[indexMatch(i,c("c","f"), findlow=TRUE)], c("c","f")))
  stopifnot(identical(i[indexMatch(i,c("c","f"), findlow=TRUE, what="pos")], c("c","f")))
  stopifnot(identical(indexMatch(i,c("c","f"), findlow=TRUE, what="val"), c("c","f")))
  stopifnot(identical(indexMatch(i,c("c","f"), findlow=TRUE, what="pos"), as.integer(c(5,2))))
  stopifnot(identical(indexMatch(i,c("c","f"), findlow=FALSE, what="pos"), as.integer(c(11,2))))

  success <- success && binregtest(
    indexMatch
  , rindexMatch
  , PAR1=list(i)
  , PAR2=list(ri)
  , x = list(c("c","f"))
  , findlow=list(missing, FALSE, TRUE)
  , what = list(missing, "ind", "pos", "val")
  , COMP=identical
  , NAME="COMPARE indexMatch with rindexMatch"
  )


  stopifnot(identical(indexFindInterval(i,"","z"), 1:14))
  stopifnot(identical(indexFindInterval(i,"","z", low.include=FALSE, high.include=FALSE), 1:14))
  stopifnot(identical(indexFindInterval(i,"","z", low.include=TRUE, high.include=TRUE), 1:14))
  stopifnot(identical(indexFindInterval(i,"","z", low.exact=TRUE), integer()))
  stopifnot(identical(indexFindInterval(i,"","z", high.exact=TRUE), integer()))

  stopifnot(identical(indexFindInterval(i,"b","f"), 2:14))
  stopifnot(identical(indexFindInterval(i,"b","f", low.include=FALSE), 3:14))
  stopifnot(identical(indexFindInterval(i,"b","f", high.include=FALSE), 2:13))
  stopifnot(identical(indexFindInterval(i,"b","f", low.include=TRUE), 2:14))
  stopifnot(identical(indexFindInterval(i,"b","f", high.include=TRUE), 2:14))

  stopifnot(identical(indexFindInterval(i,"a1","e1"), 2:13))
  stopifnot(identical(indexFindInterval(i,"a1","e1", low.include=FALSE, high.include=FALSE), 2:13))
  stopifnot(identical(indexFindInterval(i,"a1","e1", low.include=TRUE, high.include=TRUE), 2:13))
  stopifnot(identical(indexFindInterval(i,"a1","e1", low.exact=TRUE), integer()))
  stopifnot(identical(indexFindInterval(i,"a1","e1", high.exact=TRUE), integer()))

  stopifnot(identical(indexFindInterval(i,"","", FUN=indexFindlike), 1:14))
  stopifnot(identical(indexFindInterval(i,"c","d", FUN=indexFindlike), 3:12))
  stopifnot(identical(indexFindInterval(i,"c","d", FUN=indexFindlike, low.include=FALSE), 8:12))
  stopifnot(identical(indexFindInterval(i,"c","d", FUN=indexFindlike, high.include=FALSE), 3:7))
  stopifnot(identical(indexFindInterval(i,"c","d", FUN=indexFindlike, low.include=FALSE, high.include=FALSE), integer()))
  stopifnot(identical(indexFindInterval(i,"c","d", FUN=indexFindlike, low.exact=TRUE, high.exact=TRUE), 3:12))
  stopifnot(identical(indexFindInterval(i,"c","d", highFUN=indexFindlike), 3:12))
  stopifnot(identical(indexFindInterval(i,"c","d", lowFUN=indexFindlike), 3:7))

  success <- success && binregtest(
    indexFindInterval
  , rindexFindInterval
  , PAR1=list(i)
  , PAR2=list(ri)
  , list("b")
  , list("f")
  , low.include=list(missing, FALSE, TRUE)
  , high.include=list(missing, FALSE, TRUE)
  , low.exact=list(missing, FALSE, TRUE)
  , high.exact=list(missing, FALSE, TRUE)
  , lowFUN=list(missing, function(obj, ...)if (inherits(obj, "index")) indexFind(obj, ...) else rindexFind(obj, ...), function(obj, ...)if (inherits(obj, "index")) indexFindlike(obj, ...) else rindexFindlike(obj, ...))
  , highFUN=list(missing, function(obj, ...)if (inherits(obj, "index")) indexFind(obj, ...) else rindexFind(obj, ...), function(obj, ...)if (inherits(obj, "index")) indexFindlike(obj, ...) else rindexFindlike(obj, ...))
  , COMP=identical
  , NAME="COMPARE indexFindInterval with rindexFindInterval"
  )


  stopifnot(identical(indexEQ(i,"c"),c(5L, 7L, 8L, 10L, 11L)))
  stopifnot(identical(indexNE(i,"c"),c(6L, 3L, 1L, 13L, 4L, 14L, 9L, 15L, 2L)))
  stopifnot(identical(indexLT(i,"c"),c(6L, 3L)))
  stopifnot(identical(indexLE(i,"c"),c(6L, 3L, 5L, 7L, 8L, 10L, 11L)))
  stopifnot(identical(indexGT(i,"c"),c(1L, 13L, 4L, 14L, 9L, 15L, 2L)))
  stopifnot(identical(indexGE(i,"c"),c(5L, 7L, 8L, 10L, 11L, 1L, 13L, 4L, 14L, 9L, 15L, 2L)))
  stopifnot(identical(indexEQ(i,c("b","c","z",NA)), list(3L, c(5L, 7L, 8L, 10L, 11L), integer(0), 12L)))
  stopifnot(identical(indexEQ(i,c("b","c","z",NA), what="val"), list("b", c("c", "c", "c", "c", "c"), character(0), NA_character_)))

  for (j in list(list(indexEQ,rindexEQ), list(indexNE,rindexNE), list(indexLT,rindexLT), list(indexLE,rindexLE), list(indexGT,rindexGT), list(indexGE,rindexGE))){
   success <- success && binregtest(
      j[[1]]
    , j[[2]]
    , PAR1=list(i)
    , PAR2=list(ri)
    , list("c")
    , low.exact=list(missing, FALSE, TRUE)
    , high.exact=list(missing, FALSE, TRUE)
    , lowFUN=list(missing, function(obj, ...)if (inherits(obj, "index")) indexFind(obj, ...) else rindexFind(obj, ...), function(obj, ...)if (inherits(obj, "index")) indexFindlike(obj, ...) else rindexFindlike(obj, ...))
    , highFUN=list(missing, function(obj, ...)if (inherits(obj, "index")) indexFind(obj, ...) else rindexFind(obj, ...), function(obj, ...)if (inherits(obj, "index")) indexFindlike(obj, ...) else rindexFindlike(obj, ...))
    , COMP=identical
    , NAME="COMPARE indexXX with rindexXX"
    )
  }

  stopifnot(identical(i=="c",x=="c"))
  stopifnot(identical(i!="c",x!="c"))
  stopifnot(identical(i<"c",x<"c"))
  stopifnot(identical(i<="c",x<="c"))
  stopifnot(identical(i>"c",x>"c"))
  stopifnot(identical(i>="c",x>="c"))

  stopifnot(identical(ri=="c",x=="c"))
  stopifnot(identical(ri!="c",x!="c"))
  stopifnot(identical(ri<"c",x<"c"))
  stopifnot(identical(ri<="c",x<="c"))
  stopifnot(identical(ri>"c",x>"c"))
  stopifnot(identical(ri>="c",x>="c"))

  stopifnot(identical(match.index(c("b","c","z",NA), i), match(c("b","c","z",NA), x)))
  stopifnot(identical(match.rindex(c("b","c","z",NA), ri), match(c("b","c","z",NA), x)))

  stopifnot(success)

  cat("8 warnings 'indexMatch used with non-unique index' are expected\n")

  
## End Don't show




cleanEx()
nameEx("indexDemoClose")
### * indexDemoClose

flush(stderr()); flush(stdout())

### Name: indexDemoClose
### Title: Demo functions for creating and removing external pointers
### Aliases: indexDemoOpen indexDemoClose
### Keywords: misc

### ** Examples

ptr <- indexDemoOpen()
rm(ptr)
gc()

ptr <- indexDemoOpen()
indexDemoClose(ptr)
rm(ptr)
gc()



### * <FOOTER>
###
cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
quit('no')
