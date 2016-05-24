# /*
# R-Code for matching and other functions based on hashing
# S3 atomic 64bit integers for R
# (c) 2012 Jens Oehlschägel
# Licence: GPL2
# Provided 'as is', use at your own risk
# Created: 2011-12-11
# Last changed:  2011-12-11
# */

#! \name{benchmark64}
#! \alias{benchmark64}
#! \alias{optimizer64}
#! \title{
#!  Function for measuring algorithmic performance \cr 
#!  of high-level and low-level integer64 functions
#! }
#! \description{
#!  \code{benchmark64} compares high-level integer64 functions against the integer functions from Base R \cr
#!  \code{optimizer64} compares for each high-level integer64 function the Base R integer function with several low-level integer64 functions with and without caching \cr
#! }
#! \usage{
#! benchmark64(nsmall = 2^16, nbig = 2^25, timefun = repeat.time
#! )
#! optimizer64(nsmall = 2^16, nbig = 2^25, timefun = repeat.time
#! , what = c("match", "\%in\%", "duplicated", "unique", "unipos", "table", "rank", "quantile")
#! , uniorder = c("original", "values", "any")
#! , taborder = c("values", "counts")
#! , plot = TRUE
#! )
#! }
#! \arguments{
#!   \item{nsmall}{ size of smaller vector }
#!   \item{nbig}{ size of larger bigger vector }
#!   \item{timefun}{ a function for timing such as \code{\link[bit]{repeat.time}} or \code{\link{system.time}} }
#!   \item{what}{
#!  a vector of names of high-level functions
#! }
#!   \item{uniorder}{
#!  one of the order parameters that are allowed in \code{\link{unique.integer64}} and \code{\link{unipos.integer64}}
#! }
#!   \item{taborder}{
#!  one of the order parameters that are allowed in \code{\link{table.integer64}} 
#! }
#!   \item{plot}{
#!  set to FALSE to suppress plotting 
#! }
#! }
#! \details{
#!  \code{benchmark64} compares the following scenarios for the following use cases: 
#!  \tabular{rl}{
#!   \bold{scenario name} \tab \bold{explanation} \cr
#!   32-bit  \tab applying Base R function to 32-bit integer data \cr
#!   64-bit \tab applying bit64 function to 64-bit integer data (with no cache) \cr
#!   hashcache \tab dito when cache contains \code{\link{hashmap}}, see \code{\link{hashcache}} \cr
#!   sortordercache \tab dito when cache contains sorting and ordering, see \code{\link{sortordercache}} \cr
#!   ordercache \tab dito when cache contains ordering only, see \code{\link{ordercache}} \cr
#!   allcache \tab dito when cache contains sorting, ordering and hashing \cr
#!  }
#!  \tabular{rl}{
#!   \bold{use case name} \tab \bold{explanation} \cr
#!   cache         \tab filling the cache according to scenario \cr
#!   match(s,b)    \tab match small in big vector \cr
#!   s \%in\% b      \tab small \%in\% big vector \cr
#!   match(b,s)    \tab match big in small vector \cr
#!   b \%in\% s      \tab big \%in\% small vector \cr
#!   match(b,b)    \tab match big in (different) big vector \cr
#!   b \%in\% b      \tab big \%in\% (different) big vector \cr
#!   duplicated(b) \tab duplicated of big vector \cr
#!   unique(b)     \tab unique of big vector \cr
#!   table(b)      \tab table of big vector \cr
#!   sort(b)       \tab sorting of big vector \cr
#!   order(b)      \tab ordering of big vector \cr
#!   rank(b)       \tab ranking of big vector \cr
#!   quantile(b)   \tab quantiles of big vector \cr
#!   summary(b)    \tab summary of of big vector \cr
#!   SESSION       \tab exemplary session involving multiple calls (including cache filling costs) \cr
#!  }
#!  Note that the timings for the cached variants do \emph{not} contain the time costs of building the cache, except for the timing of the exemplary user session, where the cache costs are included in order to evaluate amortization. 
#! }
#! \value{
#!  \code{benchmark64} returns a matrix with elapsed seconds, different high-level tasks in rows and different scenarios to solve the task in columns. The last row named 'SESSION' contains the elapsed seconds of the exemplary sesssion.
#!  \cr
#!  \code{optimizer64} returns a dimensioned list with one row for each high-level function timed and two columns named after the values of the \code{nsmall} and \code{nbig} sample sizes. Each list cell contains a matrix with timings, low-level-methods in rows and three measurements \code{c("prep","both","use")} in columns. If it can be measured separately, \code{prep} contains the timing of preparatory work such as sorting and hashing, and \code{use} contains the timing of using the prepared work. If the function timed does both, preparation and use, the timing is in \code{both}.  
#! }
#! \author{
#!  Jens Oehlschlägel <Jens.Oehlschlaegel@truecluster.com>
#! }
#! \seealso{
#!  \code{\link{integer64}}
#! }
#! \examples{
#! message("this small example using system.time does not give serious timings\n
#! this we do this only to run regression tests")
#! benchmark64(nsmall=2^7, nbig=2^13, timefun=function(expr)system.time(expr, gcFirst=FALSE))
#! optimizer64(nsmall=2^7, nbig=2^13, timefun=function(expr)system.time(expr, gcFirst=FALSE)
#! , plot=FALSE
#! )
#!\dontrun{
#! message("for real measurement of sufficiently large datasets run this on your machine")
#! benchmark64()
#! optimizer64()
#!}
#! message("let's look at the performance results on Core i7 Lenovo T410 with 8 GB RAM")
#! data(benchmark64.data)
#! print(benchmark64.data)
#! 
#! matplot(log2(benchmark64.data[-1,1]/benchmark64.data[-1,])
#! , pch=c("3", "6", "h", "s", "o", "a") 
#! , xlab="tasks [last=session]"
#! , ylab="log2(relative speed) [bigger is better]"
#! )
#! matplot(t(log2(benchmark64.data[-1,1]/benchmark64.data[-1,]))
#! , type="b", axes=FALSE 
#! , lwd=c(rep(1, 14), 3)
#! , xlab="context"
#! , ylab="log2(relative speed) [bigger is better]"
#! )
#! axis(1
#! , labels=c("32-bit", "64-bit", "hash", "sortorder", "order", "hash+sortorder")
#! , at=1:6
#! )
#! axis(2)
#! data(optimizer64.data)
#! print(optimizer64.data)
#! oldpar <- par(no.readonly = TRUE)
#! par(mfrow=c(2,1))
#! par(cex=0.7)
#! for (i in 1:nrow(optimizer64.data)){
#!  for (j in 1:2){
#!    tim <- optimizer64.data[[i,j]]
#!   barplot(t(tim))
#!   if (rownames(optimizer64.data)[i]=="match")
#!    title(paste("match", colnames(optimizer64.data)[j], "in", colnames(optimizer64.data)[3-j]))
#!   else if (rownames(optimizer64.data)[i]=="\%in\%")
#!    title(paste(colnames(optimizer64.data)[j], "\%in\%", colnames(optimizer64.data)[3-j]))
#!   else
#!    title(paste(rownames(optimizer64.data)[i], colnames(optimizer64.data)[j]))
#!  }
#! }
#! par(mfrow=c(1,1))
#!}
#! \keyword{ misc }

#! \name{benchmark64.data}
#! \alias{benchmark64.data}
#! \docType{data}
#! \title{
#!  Results of performance measurement on a Core i7 Lenovo T410 8 GB RAM under Windows 7 64bit
#! }
#! \description{
#!   These are the results of calling \code{\link{benchmark64}}
#! }
#! \usage{data(benchmark64.data)}
#! \format{
#!   The format is:
#!  num [1:16, 1:6] 2.55e-05 2.37 2.39 1.28 1.39 ...
#!  - attr(*, "dimnames")=List of 2
#!   ..$ : chr [1:16] "cache" "match(s,b)" "s \%in\% b" "match(b,s)" ...
#!   ..$ : chr [1:6] "32-bit" "64-bit" "hashcache" "sortordercache" ...
#! }
#! \examples{
#! data(benchmark64.data)
#! print(benchmark64.data)
#! matplot(log2(benchmark64.data[-1,1]/benchmark64.data[-1,])
#! , pch=c("3", "6", "h", "s", "o", "a")
#! , xlab="tasks [last=session]"
#! , ylab="log2(relative speed) [bigger is better]"
#! )
#! matplot(t(log2(benchmark64.data[-1,1]/benchmark64.data[-1,]))
#! , axes=FALSE
#! , type="b"
#! , lwd=c(rep(1, 14), 3)
#! , xlab="context"
#! , ylab="log2(relative speed) [bigger is better]"
#! )
#! axis(1
#! , labels=c("32-bit", "64-bit", "hash", "sortorder", "order", "hash+sortorder")
#! , at=1:6
#! )
#! axis(2)
#! }
#! \keyword{datasets}


#! \name{optimizer64.data}
#! \alias{optimizer64.data}
#! \docType{data}
#! \title{
#!  Results of performance measurement on a Core i7 Lenovo T410 8 GB RAM under Windows 7 64bit
#! }
#! \description{
#!   These are the results of calling \code{\link{optimizer64}}
#! }
#! \usage{data(optimizer64.data)}
#! \format{
#!   The format is:
#! List of 16
#!  $ : num [1:9, 1:3] 0 0 1.63 0.00114 2.44 ...
#!   ..- attr(*, "dimnames")=List of 2
#!   .. ..$ : chr [1:9] "match" "match.64" "hashpos" "hashrev" ...
#!   .. ..$ : chr [1:3] "prep" "both" "use"
#!  $ : num [1:10, 1:3] 0 0 0 1.62 0.00114 ...
#!   ..- attr(*, "dimnames")=List of 2
#!   .. ..$ : chr [1:10] "\%in\%" "match.64" "\%in\%.64" "hashfin" ...
#!   .. ..$ : chr [1:3] "prep" "both" "use"
#!  $ : num [1:10, 1:3] 0 0 0.00105 0.00313 0.00313 ...
#!   ..- attr(*, "dimnames")=List of 2
#!   .. ..$ : chr [1:10] "duplicated" "duplicated.64" "hashdup" "sortorderdup1" ...
#!   .. ..$ : chr [1:3] "prep" "both" "use"
#!  $ : num [1:15, 1:3] 0 0 0 0.00104 0.00104 ...
#!   ..- attr(*, "dimnames")=List of 2
#!   .. ..$ : chr [1:15] "unique" "unique.64" "hashmapuni" "hashuni" ...
#!   .. ..$ : chr [1:3] "prep" "both" "use"
#!  $ : num [1:14, 1:3] 0 0 0 0.000992 0.000992 ...
#!   ..- attr(*, "dimnames")=List of 2
#!   .. ..$ : chr [1:14] "unique" "unipos.64" "hashmapupo" "hashupo" ...
#!   .. ..$ : chr [1:3] "prep" "both" "use"
#!  $ : num [1:13, 1:3] 0 0 0 0 0.000419 ...
#!   ..- attr(*, "dimnames")=List of 2
#!   .. ..$ : chr [1:13] "tabulate" "table" "table.64" "hashmaptab" ...
#!   .. ..$ : chr [1:3] "prep" "both" "use"
#!  $ : num [1:7, 1:3] 0 0 0 0.00236 0.00714 ...
#!   ..- attr(*, "dimnames")=List of 2
#!   .. ..$ : chr [1:7] "rank" "rank.keep" "rank.64" "sortorderrnk" ...
#!   .. ..$ : chr [1:3] "prep" "both" "use"
#!  $ : num [1:6, 1:3] 0 0 0.00189 0.00714 0 ...
#!   ..- attr(*, "dimnames")=List of 2
#!   .. ..$ : chr [1:6] "quantile" "quantile.64" "sortqtl" "orderqtl" ...
#!   .. ..$ : chr [1:3] "prep" "both" "use"
#!  $ : num [1:9, 1:3] 0 0 0.00105 1.17 0 ...
#!   ..- attr(*, "dimnames")=List of 2
#!   .. ..$ : chr [1:9] "match" "match.64" "hashpos" "hashrev" ...
#!   .. ..$ : chr [1:3] "prep" "both" "use"
#!  $ : num [1:10, 1:3] 0 0 0 0.00104 1.18 ...
#!   ..- attr(*, "dimnames")=List of 2
#!   .. ..$ : chr [1:10] "\%in\%" "match.64" "\%in\%.64" "hashfin" ...
#!   .. ..$ : chr [1:3] "prep" "both" "use"
#!  $ : num [1:10, 1:3] 0 0 1.64 2.48 2.48 ...
#!   ..- attr(*, "dimnames")=List of 2
#!   .. ..$ : chr [1:10] "duplicated" "duplicated.64" "hashdup" "sortorderdup1" ...
#!   .. ..$ : chr [1:3] "prep" "both" "use"
#!  $ : num [1:15, 1:3] 0 0 0 1.64 1.64 ...
#!   ..- attr(*, "dimnames")=List of 2
#!   .. ..$ : chr [1:15] "unique" "unique.64" "hashmapuni" "hashuni" ...
#!   .. ..$ : chr [1:3] "prep" "both" "use"
#!  $ : num [1:14, 1:3] 0 0 0 1.62 1.62 ...
#!   ..- attr(*, "dimnames")=List of 2
#!   .. ..$ : chr [1:14] "unique" "unipos.64" "hashmapupo" "hashupo" ...
#!   .. ..$ : chr [1:3] "prep" "both" "use"
#!  $ : num [1:13, 1:3] 0 0 0 0 0.32 ...
#!   ..- attr(*, "dimnames")=List of 2
#!   .. ..$ : chr [1:13] "tabulate" "table" "table.64" "hashmaptab" ...
#!   .. ..$ : chr [1:3] "prep" "both" "use"
#!  $ : num [1:7, 1:3] 0 0 0 2.96 10.69 ...
#!   ..- attr(*, "dimnames")=List of 2
#!   .. ..$ : chr [1:7] "rank" "rank.keep" "rank.64" "sortorderrnk" ...
#!   .. ..$ : chr [1:3] "prep" "both" "use"
#!  $ : num [1:6, 1:3] 0 0 1.62 10.61 0 ...
#!   ..- attr(*, "dimnames")=List of 2
#!   .. ..$ : chr [1:6] "quantile" "quantile.64" "sortqtl" "orderqtl" ...
#!   .. ..$ : chr [1:3] "prep" "both" "use"
#!  - attr(*, "dim")= int [1:2] 8 2
#!  - attr(*, "dimnames")=List of 2
#!   ..$ : chr [1:8] "match" "\%in\%" "duplicated" "unique" ...
#!   ..$ : chr [1:2] "65536" "33554432"
#! }
#! \examples{
#! data(optimizer64.data)
#! print(optimizer64.data)
#! oldpar <- par(no.readonly = TRUE)
#! par(mfrow=c(2,1))
#! par(cex=0.7)
#! for (i in 1:nrow(optimizer64.data)){
#!  for (j in 1:2){
#!    tim <- optimizer64.data[[i,j]]
#!   barplot(t(tim))
#!   if (rownames(optimizer64.data)[i]=="match")
#!    title(paste("match", colnames(optimizer64.data)[j], "in", colnames(optimizer64.data)[3-j]))
#!   else if (rownames(optimizer64.data)[i]=="\%in\%")
#!    title(paste(colnames(optimizer64.data)[j], "\%in\%", colnames(optimizer64.data)[3-j]))
#!   else
#!    title(paste(rownames(optimizer64.data)[i], colnames(optimizer64.data)[j]))
#!  }
#! }
#! par(mfrow=c(1,1))
#! }
#! \keyword{datasets}


benchmark64 <- function(nsmall=2^16, nbig=2^25, timefun=repeat.time)
{
 
 message('\ncompare performance for a complete sessions of calls')
 s <- sample(nbig, nsmall, TRUE)
 b <- sample(nbig, nbig, TRUE)
 b2 <- sample(nbig, nbig, TRUE)
 
 tim1 <- double(6)
 names(tim1) <- c("32-bit","64-bit","hashcache","sortordercache","ordercache","allcache")

 s <- as.integer(s)
 b <- as.integer(b)
 b2 <- as.integer(b2)
 
 i <- 1
 for (i in 1:6){
  message("\n=== ", names(tim1)[i], " ===")
  
  if (i==2){
   s <- as.integer64(s)
   b <- as.integer64(b)
   b2 <- as.integer64(b2)
  }

  tim1[i] <- 0
 
  tim1[i] <- tim1[i] + timefun({
   switch(as.character(i)
   , "3" = {hashcache(s); hashcache(b); hashcache(b2)}
   , "4" = {sortordercache(s); sortordercache(b); sortordercache(b2)}
   , "5" = {ordercache(s); ordercache(b); ordercache(b2)}
   , "6" = {hashcache(s); hashcache(b); hashcache(b2);sortordercache(s); sortordercache(b); sortordercache(b2)}
   )
  })[3]
 
  message('check data range, mean etc.')
  tim1[i] <- tim1[i] + timefun({
   summary(b)
  })[3]
  message('get all percentiles for plotting distribution shape')
  tim1[i] <- tim1[i] + timefun({
   quantile(b, probs=seq(0, 1, 0.01))
  })[3]
  message('list the upper and lower permille of values')
  tim1[i] <- tim1[i] + timefun({
   quantile(b, probs=c(0.001, 0.999))
   sort(b, na.last=NA)
  })[3]
  message('OK, for some of these values I want to see the complete ROW, so I need their positions in the data.frame')
  tim1[i] <- tim1[i] + timefun({
   if(i==1)order(b) else order.integer64(b)
  })[3]
  message('check if any values are duplicated')
  tim1[i] <- tim1[i] + timefun({
   any(duplicated(b))
  })[3]
  message('since not unique, then check distribution of frequencies')
  tim1[i] <- tim1[i] + timefun({
   if(i==1)tabulate(table(b, exclude=NULL)) else tabulate(table.integer64(b, return='list')$counts)
  })[3]
  message("OK, let's plot the percentiles of unique values versus the percentiles allowing for duplicates")
  tim1[i] <- tim1[i] + timefun({
   quantile(b, probs=seq(0, 1, 0.01))
   quantile(unique(b), probs=seq(0, 1, 0.01))
  })[3]
  message('check whether we find a match for each fact in the dimension table')
  tim1[i] <- tim1[i] + timefun({
   all(if(i==1) b %in% s else "%in%.integer64"(b, s))
  })[3]
  message('check whether there are any dimension table entries not in the fact table')
  tim1[i] <- tim1[i] + timefun({
   all(if(i==1) s %in% b else "%in%.integer64"(s, b))
  })[3]
  message('check whether we find a match for each fact in a parallel fact table')
  tim1[i] <- tim1[i] + timefun({
   all(if(i==1) b %in% b2 else "%in%.integer64"(b, b2))
  })[3]
  message('find positions of facts in dimension table for joining')
  tim1[i] <- tim1[i] + timefun({
   if(i==1) match(b, s) else match.integer64(b, s)
  })[3]
  message('find positions of facts in parallel fact table for joining')
  tim1[i] <- tim1[i] + timefun({
   if(i==1) match(b, b2) else match.integer64(b, b2)
  })[3]
  message('out of curiosity: how well rank-correlated are fact and parallel fact table?')
  tim1[i] <- tim1[i] + timefun({
   if (i==1){
    cor(rank(b, na.last="keep"), rank(b2, na.last="keep"), use="na.or.complete")
   }else{
    cor(rank.integer64(b), rank.integer64(b2), use="na.or.complete")
   }
  })[3]
  
  remcache(s)
  remcache(b)
  remcache(b2)
  
  print(round(rbind(seconds=tim1, factor=tim1[1]/tim1), 3))
 
 }

        # 32-bit         64-bit      hashcache sortordercache     ordercache       allcache 
       # 196.510          8.963          8.242          5.183         12.325          6.043 
        # 32-bit         64-bit      hashcache sortordercache     ordercache       allcache 
         # 1.000         21.924         23.842         37.913         15.944         32.519 

   
 message("\nnow let's look more systematically at the components involved")
 s <- sample(nbig, nsmall, TRUE)
 b <- sample(nbig, nbig, TRUE)
 b2 <- sample(nbig, nbig, TRUE)
 
 tim2 <- matrix(0, 15, 6)
 dimnames(tim2) <- list(c("cache", "match(s,b)", "s %in% b", "match(b,s)", "b %in% s", "match(b,b)", "b %in% b", "duplicated(b)", "unique(b)", "table(b)", "sort(b)", "order(b)", "rank(b)", "quantile(b)", "summary(b)")
 , c("32-bit","64-bit","hashcache","sortordercache","ordercache","allcache"))

 s <- as.integer(s)
 b <- as.integer(b)
 b2 <- as.integer(b2)
 
 i <- 1
 for (i in 1:6){
  if (i==2){
   s <- as.integer64(s)
   b <- as.integer64(b)
   b2 <- as.integer64(b2)
  }
 
  if (i>2)message(colnames(tim2)[i], " cache")
  tim2["cache",i] <- timefun({
   switch(as.character(i)
   , "3" = {hashcache(s); hashcache(b); hashcache(b2)}
   , "4" = {sortordercache(s); sortordercache(b); sortordercache(b2)}
   , "5" = {ordercache(s); ordercache(b); ordercache(b2)}
   , "6" = {hashcache(s); hashcache(b); hashcache(b2);sortordercache(s); sortordercache(b); sortordercache(b2)}
   )
  })[3]
 
  message(colnames(tim2)[i], " match(s,b)")
  tim2["match(s,b)",i] <- timefun({
   if (i==1) match(s, b) else match.integer64(s, b)
  })[3]
 
  message(colnames(tim2)[i], " s %in% b")
  tim2["s %in% b",i] <- timefun({
   if (i==1) s %in% b else "%in%.integer64"(s,b)
  })[3]
 
  message(colnames(tim2)[i], " match(b,s)")
  tim2["match(b,s)",i] <- timefun({
   if (i==1) match(b, s) else match.integer64(b, s)
  })[3]
 
  message(colnames(tim2)[i], " b %in% s")
  tim2["b %in% s",i] <- timefun({
   if (i==1) b %in% s else "%in%.integer64"(b,s)
  })[3]
 
  message(colnames(tim2)[i], " match(b,b)")
  tim2["match(b,b)",i] <- timefun({
   if (i==1) match(b, b2) else match.integer64(b, b2)
  })[3]
 
  message(colnames(tim2)[i], " b %in% b")
  tim2["b %in% b",i] <- timefun({
   if (i==1) b %in% b2 else "%in%.integer64"(b,b2)
  })[3]
 
  message(colnames(tim2)[i], " duplicated(b)")
  tim2["duplicated(b)",i] <- timefun({
   duplicated(b)
  })[3]
 
  message(colnames(tim2)[i], " unique(b)")
  tim2["unique(b)",i] <- timefun({
   unique(b)
  })[3]
 
  message(colnames(tim2)[i], " table(b)")
  tim2["table(b)",i] <- timefun({
   if(i==1) table(b) else table.integer64(b, return='list')
  })[3]
 
  message(colnames(tim2)[i], " sort(b)")
  tim2["sort(b)",i] <- timefun({
   sort(b)
  })[3]
 
  message(colnames(tim2)[i], " order(b)")
  tim2["order(b)",i] <- timefun({
   if(i==1) order(b) else order.integer64(b)
  })[3]
 
  message(colnames(tim2)[i], " rank(b)")
  tim2["rank(b)",i] <- timefun({
   if(i==1) rank(b) else rank.integer64(b)
  })[3]
 
  message(colnames(tim2)[i], " quantile(b)")
  tim2["quantile(b)",i] <- timefun({
   quantile(b)
  })[3]
 
  message(colnames(tim2)[i], " summary(b)")
  tim2["summary(b)",i] <- timefun({
   summary(b)
  })[3]
  
  remcache(s)
  remcache(b)
  remcache(b2)
  
  tim3 <- rbind(tim2, SESSION=tim1)
  #tim2 <- tim2[,1]/tim2
  
  cat("seconds")
  print(round(tim3, 3))
  cat("factor")
  print(round(tim3[,1]/tim3, 3))
 
 }


 
               # 32-bit 64-bit hashcache sortordercache ordercache allcache
# cache           0.000  0.000     0.775          1.330      6.500    2.660
# match(s,b)      0.820  0.218     0.004          0.025      0.093    0.004
# s %in% b        0.810  0.234     0.003          0.022      0.093    0.003
# match(b,s)      0.450  0.228     0.232          0.224      0.224    0.226
# b %in% s        0.510  0.226     0.224          0.222      0.218    0.222
# match(b,b)      2.370  0.870     0.505          0.890      0.880    0.505
# b %in% b        2.350  0.850     0.480          0.865      0.870    0.483
# duplicated(b)   0.875  0.510     0.141          0.116      0.383    0.117
# unique(b)       0.930  0.555     0.447          0.156      0.427    0.450
# table(b)      110.340  0.725     0.680          0.234      0.575    0.202
# sort(b)         2.440  0.400     0.433          0.072      0.460    0.069
# order(b)       12.780  0.680     0.615          0.036      0.036    0.035
# rank(b)        13.480  0.860     0.915          0.240      0.545    0.246
# quantile(b)     0.373  0.400     0.410          0.000      0.000    0.000
# summary(b)      0.645  0.423     0.427          0.016      0.016    0.016
# TOTAL         149.173  7.179     6.291          4.448     11.320    5.239
              # 32-bit  64-bit hashcache sortordercache ordercache allcache
# cache              1   1.062     0.000          0.000      0.000    0.000
# match(s,b)         1   3.761   230.420         32.475      8.843  217.300
# s %in% b           1   3.462   234.090         36.450      8.735  237.386
# match(b,s)         1   1.974     1.940          2.009      2.009    1.991
# b %in% s           1   2.257     2.277          2.297      2.339    2.297
# match(b,b)         1   2.724     4.693          2.663      2.693    4.693
# b %in% b           1   2.765     4.896          2.717      2.701    4.862
# duplicated(b)      1   1.716     6.195          7.572      2.283    7.500
# unique(b)          1   1.676     2.082          5.972      2.180    2.067
# table(b)           1 152.193   162.265        471.538    191.896  546.238
# sort(b)            1   6.100     5.631         33.822      5.304   35.534
# order(b)           1  18.794    20.780        357.840    354.297  366.950
# rank(b)            1  15.674    14.732         56.167     24.734   54.797
# quantile(b)        1   0.933     0.911        804.907    806.027  810.133
# summary(b)         1   1.524     1.512         39.345     39.345   39.345
# TOTAL              1  20.778    23.712         33.534     13.177   28.476 

  tim3
}


optimizer64 <- function(nsmall=2^16, nbig=2^25, timefun=repeat.time
, what=c("match","%in%","duplicated","unique","unipos","table","rank","quantile")
, uniorder = c("original", "values", "any")
, taborder = c("values", "counts")
, plot = TRUE
)
{
 uniorder <- match.arg(uniorder)
 taborder <- match.arg(taborder)
 ret <- vector("list", 2*length(what))
 dim(ret) <- c(length(what), 2L)
 dimnames(ret) <- list(what, c(nsmall, nbig))
 
 if (plot){
  oldpar <- par(no.readonly = TRUE)
  on.exit(par(oldpar))
  par(mfrow=c(2,1))
 }
 
 if ("match" %in% what){
  message("match: timings of different methods")
  N1 <- c(nsmall, nbig)
  N2 <- c(nbig, nsmall)
  for (i in seq_along(N1)){
   n1 <- N1[i]
   n2 <- N2[i]
   x1 <- c(sample(n2, n1-1, TRUE), NA)
   x2 <- c(sample(n2, n2-1, TRUE), NA)
   tim <- matrix(0, 9, 3)
   dimnames(tim) <- list(c("match","match.64","hashpos","hashrev","sortorderpos","orderpos","hashcache","sortorder.cache","order.cache"), c("prep","both","use"))

   tim["match","both"] <- timefun({
    p <- match(x1, x2)
   })[3]
   x1 <- as.integer64(x1)
   x2 <- as.integer64(x2)

   tim["match.64","both"] <- timefun({
    p2 <- match.integer64(x1, x2)
   })[3]
   stopifnot(identical(p2, p))

   tim["hashpos","prep"] <- timefun({
    h2 <- hashmap(x2)
   })[3]
   tim["hashpos","use"] <- timefun({
    p2 <- hashpos(h2, x1)
   })[3]
   stopifnot(identical(p2, p))
   
   tim["hashrev","prep"] <- timefun({
    h1 <- hashmap(x1)
   })[3]
   tim["hashrev","use"] <- timefun({
    p1 <- hashrev(h1, x2)
   })[3]
   stopifnot(identical(p2, p))
   
   tim["sortorderpos","prep"] <- system.time({
    s2 <- clone(x2)
    o2 <- seq_along(x2)
    ramsortorder(s2, o2, na.last=FALSE)
   })[3]
   tim["sortorderpos","use"] <- timefun({
    p2 <- sortorderpos(s2, o2, x1)
   })[3]
   stopifnot(identical(p2, p))
   
   tim["orderpos","prep"] <- timefun({
    o2 <- seq_along(x2)
    ramorder(x2, o2, na.last=FALSE)
   })[3]
   tim["orderpos","use"] <- timefun({
    p2 <- orderpos(x2, o2, x1, method=2)
   })[3]
   stopifnot(identical(p2, p))
   
   hashcache(x2)
   tim["hashcache","use"] <- timefun({
    p2 <- match.integer64(x1, x2)
   })[3]
   stopifnot(identical(p2, p))
   remcache(x2)
   
   sortordercache(x2)
   tim["sortorder.cache","use"] <- timefun({
    p2 <- match.integer64(x1, x2)
   })[3]
   stopifnot(identical(p2, p))
   remcache(x2)
   
   ordercache(x2)
   tim["order.cache","use"] <- timefun({
    p2 <- match.integer64(x1, x2)
   })[3]
   stopifnot(identical(p2, p))
   remcache(x2)

   if (plot){
    barplot(t(tim))
    n <- format(c(n1, n2))
    title(paste("match", n[1], "in", n[2]))
   }
   
   ret[["match",as.character(n1)]] <- tim
  }
 }

 if ("%in%" %in% what){
  message("%in%: timings of different methods")
  N1 <- c(nsmall, nbig)
  N2 <- c(nbig, nsmall)
  for (i in seq_along(N1)){
   n1 <- N1[i]
   n2 <- N2[i]
   x1 <- c(sample(n2, n1-1, TRUE), NA)
   x2 <- c(sample(n2, n2-1, TRUE), NA)
   tim <- matrix(0, 10, 3)
   dimnames(tim) <- list(c("%in%","match.64","%in%.64","hashfin","hashrin","sortfin","orderfin","hash.cache","sortorder.cache","order.cache"), c("prep","both","use"))

   tim["%in%","both"] <- timefun({
    p <- x1 %in% x2
   })[3]
   x1 <- as.integer64(x1)
   x2 <- as.integer64(x2)

   tim["match.64","both"] <- timefun({
    p2 <- match.integer64(x1,x2, nomatch = 0L) > 0L
   })[3]
   stopifnot(identical(p2, p))

   tim["%in%.64","both"] <- timefun({
    p2 <- "%in%.integer64"(x1,x2) # this is using the custom version
   })[3]
   stopifnot(identical(p2, p))

   tim["hashfin","prep"] <- timefun({
    h2 <- hashmap(x2)
   })[3]
   tim["hashfin","use"] <- timefun({
    p2 <- hashfin(h2, x1)
   })[3]
   stopifnot(identical(p2, p))
   
   tim["hashrin","prep"] <- timefun({
    h1 <- hashmap(x1)
   })[3]
   tim["hashrin","use"] <- timefun({
    p1 <- hashrin(h1, x2)
   })[3]
   stopifnot(identical(p2, p))
   
   tim["sortfin","prep"] <- timefun({
    s2 <- clone(x2)
    ramsort(s2, na.last=FALSE)
   })[3]
   tim["sortfin","use"] <- timefun({
    p2 <- sortfin(s2, x1)
   })[3]
   stopifnot(identical(p2, p))
   
   tim["orderfin","prep"] <- timefun({
    o2 <- seq_along(x2)
    ramorder(x2, o2, na.last=FALSE)
   })[3]
   tim["orderfin","use"] <- timefun({
    p2 <- orderfin(x2, o2, x1)
   })[3]
   stopifnot(identical(p2, p))
   
   hashcache(x2)
   tim["hash.cache","use"] <- timefun({
    p2 <- "%in%.integer64"(x1, x2)
   })[3]
   stopifnot(identical(p2, p))
   remcache(x2)
   
   sortordercache(x2)
   tim["sortorder.cache","use"] <- timefun({
    p2 <- "%in%.integer64"(x1, x2)
   })[3]
   stopifnot(identical(p2, p))
   remcache(x2)
   
   ordercache(x2)
   tim["order.cache","use"] <- timefun({
    p2 <- "%in%.integer64"(x1, x2)
   })[3]
   stopifnot(identical(p2, p))
   remcache(x2)
   
   if (plot){
    barplot(t(tim))
    n <- format(c(n1, n2))
    title(paste(n[1], "%in%", n[2]))
   }
    
   ret[["%in%",as.character(n1)]] <- tim
  }
 }
 if ("duplicated" %in% what){
  message("duplicated: timings of different methods")
  N <- c(nsmall, nbig)
  for (i in seq_along(N)){
   n <- N[i]
   x <- c(sample(n, n-1, TRUE), NA)
   tim <- matrix(0, 10, 3)
   dimnames(tim) <- list(c("duplicated","duplicated.64","hashdup","sortorderdup1","sortorderdup2","orderdup1","orderdup2"
    ,"hash.cache","sortorder.cache","order.cache")
   , c("prep","both","use"))

   tim["duplicated","both"] <- timefun({
    p <- duplicated(x)
   })[3]
   x <- as.integer64(x)

   tim["duplicated.64","both"] <- timefun({
    p2 <- duplicated(x)
   })[3]
   stopifnot(identical(p2, p))

   tim["hashdup","prep"] <- timefun({
    h <- hashmap(x)
   })[3]
   tim["hashdup","use"] <- timefun({
    p2 <- hashdup(h)
   })[3]
   stopifnot(identical(p2, p))
   
   tim["sortorderdup1","prep"] <- timefun({
    s <- clone(x)
    o <- seq_along(x)
    ramsortorder(s, o, na.last=FALSE)
    nunique <- sortnut(s)[1]
   })[3]
   tim["sortorderdup1","use"] <- timefun({
    p2 <- sortorderdup(s, o, method=1)
   })[3]
   stopifnot(identical(p2, p))
    
   tim["sortorderdup2","prep"] <- tim["sortorderdup1","prep"]
   tim["sortorderdup2","use"] <- timefun({
    p2 <- sortorderdup(s, o, method=2)
   })[3]
   stopifnot(identical(p2, p))
    
   tim["orderdup1","prep"] <- timefun({
    o <- seq_along(x)
    ramorder(x, o, na.last=FALSE)
    nunique <- ordernut(x,o)[1]
   })[3]
   tim["orderdup1","use"] <- timefun({
    p2 <- orderdup(x, o, method=1)
   })[3]
   stopifnot(identical(p2, p))
   
   tim["orderdup2","prep"] <- tim["orderdup1","prep"]
   tim["orderdup2","use"] <- timefun({
    p2 <- orderdup(x, o, method=2)
   })[3]
   stopifnot(identical(p2, p))
   
   hashcache(x)
   tim["hash.cache","use"] <- timefun({
    p2 <- duplicated(x)
   })[3]
   stopifnot(identical(p2, p))
   remcache(x)

   sortordercache(x)
   tim["sortorder.cache","use"] <- timefun({
    p2 <- duplicated(x)
   })[3]
   stopifnot(identical(p2, p))
   remcache(x)

   ordercache(x)
   tim["order.cache","use"] <- timefun({
    p2 <- duplicated(x)
   })[3]
   stopifnot(identical(p2, p))
   remcache(x)
   
   if (plot){
    barplot(t(tim), cex.names=0.7)
    title(paste("duplicated(",n,")", sep=""))
   }
   
   ret[["duplicated",as.character(n)]] <- tim
  }
 }
 if ("unique" %in% what){
  message("unique: timings of different methods")
  N <- c(nsmall, nbig)
  for (i in seq_along(N)){
   n <- N[i]
   x <- c(sample(n, n-1, TRUE), NA)
   tim <- matrix(0, 15, 3)
   dimnames(tim) <- list(
   c("unique","unique.64","hashmapuni","hashuni","hashunikeep","sortuni","sortunikeep","orderuni","orderunikeep","hashdup","sortorderdup"
    ,"hash.cache","sort.cache","sortorder.cache","order.cache")
   , c("prep","both","use"))

   tim["unique","both"] <- timefun({
    p <- unique(x)
   })[3]
   x <- as.integer64(x)
   p <- as.integer64(p)
   if (uniorder=="values")
    ramsort(p, na.last=FALSE)

   tim["unique.64","both"] <- timefun({
    p2 <- unique(x, order=uniorder)
   })[3]
   if (uniorder!="any")
    stopifnot(identical.integer64(p2, p))

   tim["hashmapuni","both"] <- timefun({
    p2 <- hashmapuni(x)
   })[3]
   if (uniorder=="original")
    stopifnot(identical.integer64(p2, p))
   
   tim["hashuni","prep"] <- timefun({
    h <- hashmap(x)
    # for(r in 1:r)h <- hashmap(x, nunique=h$nunique)
   })[3]
   tim["hashuni","use"] <- timefun({
    p2 <- hashuni(h)
   })[3]
   if (uniorder=="values")
    stopifnot(identical.integer64(sort(p2, na.last=FALSE), p))
   
   tim["hashunikeep","prep"] <- tim["hashuni","prep"] 
   tim["hashunikeep","use"] <- timefun({
    p2 <- hashuni(h, keep.order=TRUE)
   })[3]
   if (uniorder=="original")
    stopifnot(identical.integer64(p2, p))

   tim["sortuni","prep"] <- timefun({
    s <- clone(x)
    ramsort(s, na.last=FALSE)
    nunique <- sortnut(s)[1]
   })[3]
   tim["sortuni","use"] <- timefun({
    p2 <- sortuni(s, nunique)
   })[3]
   if (uniorder=="values")
    stopifnot(identical.integer64(sort(p2, na.last=FALSE), p))
   
   tim["sortunikeep","prep"] <- timefun({
    s <- clone(x)
    o <- seq_along(x)
    ramsortorder(s, o, na.last=FALSE)
    nunique <- sortnut(s)[1]
   })[3]
   tim["sortunikeep","use"] <- timefun({
    p2 <- sortorderuni(x, s, o, nunique)
   })[3]
   if (uniorder=="original")
    stopifnot(identical.integer64(p2, p))
    
   tim["orderuni","prep"] <- timefun({
    o <- seq_along(x)
    ramorder(x, o, na.last=FALSE)
    nunique <- ordernut(x,o)[1]
   })[3]
   tim["orderuni","use"] <- timefun({
    p2 <- orderuni(x, o, nunique)
   })[3]
   if (uniorder=="values")
    stopifnot(identical.integer64(sort(p2, na.last=FALSE), p))
   
   tim["orderunikeep","prep"] <- tim["orderuni","prep"]
   tim["orderunikeep","use"] <- timefun({
    p2 <- orderuni(x, o, nunique, keep.order=TRUE)
    nunique <- ordernut(x,o)[1]
   })[3]
   if (uniorder=="original")
    stopifnot(identical.integer64(p2, p))

   tim["hashdup","prep"] <- tim["hashuni","prep"]
   tim["hashdup","use"] <- timefun({
    p2 <- x[!hashdup(h)]
   })[3]
   if (uniorder=="original")
    stopifnot(identical.integer64(p2, p))

   tim["sortorderdup","prep"] <- tim["sortunikeep","prep"]
   tim["sortorderdup","use"] <- timefun({
    p2 <- x[!sortorderdup(s, o)]
   })[3]
   if (uniorder=="original")
    stopifnot(identical.integer64(p2, p))

   
   hashcache(x)
   tim["hash.cache","use"] <- timefun({
    p2 <- unique(x, order=uniorder)
   })[3]
   if (uniorder!="any")
    stopifnot(identical.integer64(p2, p))
   remcache(x)

   sortcache(x)
   tim["sort.cache","use"] <- timefun({
    p2 <- unique(x, order=uniorder)
   })[3]
   if (uniorder!="any")
    stopifnot(identical.integer64(p2, p))
   remcache(x)

   sortordercache(x)
   tim["sortorder.cache","use"] <- timefun({
    p2 <- unique(x, order=uniorder)
   })[3]
   if (uniorder!="any")
    stopifnot(identical.integer64(p2, p))
   remcache(x)

   ordercache(x)
   tim["order.cache","use"] <- timefun({
    p2 <- unique(x, order=uniorder)
   })[3]
   if (uniorder!="any")
    stopifnot(identical.integer64(p2, p))
   remcache(x)
   
   if (plot){
    barplot(t(tim), cex.names=0.7)
    title(paste("unique(",n,", order=",uniorder,")", sep=""))
   }
   
   ret[["unique",as.character(n)]] <- tim
  }
 }
 if ("unipos" %in% what){
  message("unipos: timings of different methods")
  N <- c(nsmall, nbig)
  for (i in seq_along(N)){
   n <- N[i]
   x <- c(sample(n, n-1, TRUE), NA)
   tim <- matrix(0, 14, 3)
   dimnames(tim) <- list(
   c("unique","unipos.64","hashmapupo","hashupo","hashupokeep","sortorderupo","sortorderupokeep","orderupo","orderupokeep","hashdup","sortorderdup"
    ,"hash.cache","sortorder.cache","order.cache")
   , c("prep","both","use"))

   tim["unique","both"] <- timefun({
    unique(x)
   })[3]
   x <- as.integer64(x)

   tim["unipos.64","both"] <- timefun({
    p <- unipos(x, order=uniorder)
   })[3]

   tim["hashmapupo","both"] <- timefun({
    p2 <- hashmapupo(x)
   })[3]
   if (uniorder=="original")
    stopifnot(identical(p2, p))
   
   tim["hashupo","prep"] <- timefun({
    h <- hashmap(x)
    # if nunique is small we could re-build the hashmap at a smaller size
    # h <- hashmap(x, nunique=h$nunique)
   })[3]
   tim["hashupo","use"] <- timefun({
    p2 <- hashupo(h)
   })[3]
   if (uniorder=="values")
    stopifnot(identical(sort(p2, na.last=FALSE), sort(p, na.last=FALSE)))
   
   tim["hashupokeep","prep"] <- tim["hashupo","prep"] 
   tim["hashupokeep","use"] <- timefun({
    p2 <- hashupo(h, keep.order=TRUE)
   })[3]
   if (uniorder=="original")
    stopifnot(identical(p2, p))

   
   tim["sortorderupo","prep"] <- timefun({
    s <- clone(x)
    o <- seq_along(x)
    ramsortorder(s, o, na.last=FALSE)
    nunique <- sortnut(s)[1]
   })[3]
   tim["sortorderupo","use"] <- timefun({
    p2 <- sortorderupo(s, o, nunique)
   })[3]
   if (uniorder=="values")
    stopifnot(identical(p2, p))
    
   tim["sortorderupokeep","prep"] <- timefun({
    s <- clone(x)
    o <- seq_along(x)
    ramsortorder(s, o, na.last=FALSE)
    nunique <- sortnut(s)[1]
   })[3]
   tim["sortorderupokeep","use"] <- timefun({
    p2 <- sortorderupo(s, o, nunique, keep.order=TRUE)
   })[3]
   if (uniorder=="original")
    stopifnot(identical(p2, p))
    
   tim["orderupo","prep"] <- timefun({
    o <- seq_along(x)
    ramorder(x, o, na.last=FALSE)
    nunique <- ordernut(x,o)[1]
   })[3]
   tim["orderupo","use"] <- timefun({
    p2 <- orderupo(x, o, nunique)
   })[3]
   if (uniorder=="values")
    stopifnot(identical(p2, p))
   
   tim["orderupokeep","prep"] <- tim["orderupo","prep"]
   tim["orderupokeep","use"] <- timefun({
    p2 <- orderupo(x, o, nunique, keep.order=TRUE)
    nunique <- ordernut(x,o)[1]
   })[3]
   if (uniorder=="original")
    stopifnot(identical(p2, p))

   tim["hashdup","prep"] <- tim["hashupo","prep"]
   tim["hashdup","use"] <- timefun({
    p2 <- (1:n)[!hashdup(h)]
   })[3]
   if (uniorder=="original")
    stopifnot(identical(p2, p))

   tim["sortorderdup","prep"] <- tim["sortorderupokeep","prep"]
   tim["sortorderdup","use"] <- timefun({
    p2 <- (1:n)[!sortorderdup(s, o)]
   })[3]
   if (uniorder=="original")
    stopifnot(identical(p2, p))
   
   hashcache(x)
   tim["hash.cache","use"] <- timefun({
    p2 <- unipos(x, order=uniorder)
   })[3]
   if (uniorder!="any")
    stopifnot(identical(p2, p))
   remcache(x)

   sortordercache(x)
   tim["sortorder.cache","use"] <- timefun({
    p2 <- unipos(x, order=uniorder)
   })[3]
   if (uniorder!="any")
    stopifnot(identical(p2, p))
   remcache(x)

   ordercache(x)
   tim["order.cache","use"] <- timefun({
    p2 <- unipos(x, order=uniorder)
   })[3]
   if (uniorder!="any")
    stopifnot(identical(p2, p))
   remcache(x)
   
   if (plot){
    barplot(t(tim), cex.names=0.7)
    title(paste("unipos(",n,", order=",uniorder,")", sep=""))
   }
   
   ret[["unipos",as.character(n)]] <- tim
  }
 }
 if ("table" %in% what){
  message("table: timings of different methods")
  N <- c(nsmall, nbig)
  for (i in seq_along(N)){
   n <- N[i]
   x <- c(sample(1024, n-1, TRUE), NA)
   tim <- matrix(0, 13, 3)
   dimnames(tim) <- list(c("tabulate","table","table.64","hashmaptab","hashtab","hashtab2","sorttab","sortordertab","ordertab","ordertabkeep"
    ,"hash.cache","sort.cache","order.cache")
   , c("prep","both","use"))

   tim["tabulate","both"] <- timefun({
    tabulate(x)
   })[3]
   
   tim["table","both"] <- timefun({
    p <- table(x, exclude=NULL)
   })[3]
   p <- p[-length(p)]
   
   x <- as.integer64(x)

   tim["table.64","both"] <- timefun({
    p2 <- table.integer64(x, order=taborder)
   })[3]
   p2 <- p2[-1]
   stopifnot(identical(p2, p))

   tim["hashmaptab","both"] <- timefun({
    p <- hashmaptab(x)
   })[3]
   
   tim["hashtab","prep"] <- timefun({
    h <- hashmap(x)
   })[3]
   tim["hashtab","use"] <- timefun({
    p2 <- hashtab(h)
   })[3]
   stopifnot(identical(p2, p))
   
   tim["hashtab2","prep"] <- tim["hashtab","prep"] + timefun({
    h <- hashmap(x, nunique=h$nunique)
   })[3]
   tim["hashtab2","use"] <- timefun({
    p2 <- hashtab(h)
   })[3]
   
   sortp <- function(p){
    s <- p$values
    o <- seq_along(s)
    ramsortorder(s,o, na.last=FALSE)
    list(values=s, counts=p$counts[o])
   }
   p <- sortp(p)
   p2 <- sortp(p2)
   stopifnot(identical(p2, p))
   
   tim["sorttab","prep"] <- timefun({
    s <- clone(x)
    ramsort(s, na.last=FALSE)
    nunique <- sortnut(s)[1]
   })[3]
   tim["sorttab","use"] <- timefun({
    p2 <- list(values=sortuni(s, nunique), counts=sorttab(s, nunique))
   })[3]
   stopifnot(identical(p2, p))
    
   tim["sortordertab","prep"] <- timefun({
    s <- clone(x)
    o <- seq_along(x)
    ramsortorder(s, o, na.last=FALSE)
    nunique <- sortnut(s)[1]
  	})[3]
			tim["sortordertab","use"] <- timefun({
				p2 <- list(values=sortorderuni(x, s, o, nunique), counts=sortordertab(s, o))
			})[3]
			p2 <- sortp(p2)
			stopifnot(identical(p2, p))
				
			tim["ordertab","prep"] <- timefun({
				o <- seq_along(x)
				ramorder(x, o, na.last=FALSE)
				nunique <- ordernut(x, o)[1]
			})[3]
			tim["ordertab","use"] <- timefun({
				p2 <- list(values=orderuni(x, o, nunique), counts=ordertab(x, o, nunique))
			})[3]
			stopifnot(identical(p2, p))
				
			tim["ordertabkeep","prep"] <- tim["ordertab","prep"] 
			tim["ordertabkeep","use"] <- timefun({
				p2 <- list(values=orderuni(x, o, nunique, keep.order=TRUE), counts=ordertab(x, o, nunique, keep.order=TRUE))
			})[3]
			p2 <- sortp(p2)
			stopifnot(identical(p2, p))
			
			hashcache(x)
			tim["hash.cache","use"] <- timefun({
				p <- table.integer64(x, order=taborder)
			})[3]
			remcache(x)

			sortordercache(x)
			tim["sort.cache","use"] <- timefun({
				p2 <- table.integer64(x, order=taborder)
			})[3]
			stopifnot(identical(p2, p))
			remcache(x)

			ordercache(x)
			tim["order.cache","use"] <- timefun({
				p2 <- table.integer64(x, order=taborder)
			})[3]
			stopifnot(identical(p2, p))
			remcache(x)
			
			if (plot){
				barplot(t(tim), cex.names=0.7)
				title(paste("table.integer64(",n,", order=",taborder,")", sep=""))
			}
			
			ret[["table",as.character(n)]] <- tim
		}
	}
	if ("rank" %in% what){
		message("rank: timings of different methods")
		N <- c(nsmall, nbig)
		for (i in seq_along(N)){
			n <- N[i]
			x <- c(sample(n, n-1, TRUE), NA)
			tim <- matrix(0, 7, 3)
			dimnames(tim) <- list(c("rank","rank.keep","rank.64","sortorderrnk","orderrnk"
				,"sort.cache","order.cache")
			, c("prep","both","use"))

			tim["rank","both"] <- timefun({
				rank(x)
			})[3]
			
			tim["rank.keep","both"] <- timefun({
				p <- rank(x, na.last="keep")
			})[3]
			
			x <- as.integer64(x)

			tim["rank.64","both"] <- timefun({
				p2 <- rank.integer64(x)
			})[3]
			stopifnot(identical(p2, p))
				
			tim["sortorderrnk","prep"] <- timefun({
				s <- clone(x)
				o <- seq_along(x)
				na.count <- ramsortorder(s, o, na.last=FALSE)
			})[3]
			tim["sortorderrnk","use"] <- timefun({
				p2 <- sortorderrnk(s, o, na.count)
			})[3]
			stopifnot(identical(p2, p))
				
			tim["orderrnk","prep"] <- timefun({
				o <- seq_along(x)
				na.count <- ramorder(x, o, na.last=FALSE)
			})[3]
			tim["orderrnk","use"] <- timefun({
				p2 <- orderrnk(x, o, na.count)
			})[3]
			stopifnot(identical(p2, p))
				
			sortordercache(x)
			tim["sort.cache","use"] <- timefun({
				p2 <- rank.integer64(x)
			})[3]
			stopifnot(identical(p2, p))
			remcache(x)

			ordercache(x)
			tim["order.cache","use"] <- timefun({
				p2 <- rank.integer64(x)
			})[3]
			stopifnot(identical(p2, p))
			remcache(x)
			
			if (plot){
				barplot(t(tim), cex.names=0.7)
				title(paste("rank.integer64(",n,")", sep=""))
			}
			
			ret[["rank",as.character(n)]] <- tim
		}
	}
	if ("quantile" %in% what){
		message("quantile: timings of different methods")
		N <- c(nsmall, nbig)
		for (i in seq_along(N)){
			n <- N[i]
			x <- c(sample(n, n-1, TRUE), NA)
			tim <- matrix(0, 6, 3)
			dimnames(tim) <- list(c("quantile","quantile.64","sortqtl","orderqtl"
				,"sort.cache","order.cache")
			, c("prep","both","use"))

			tim["quantile","both"] <- timefun({
				p <- quantile(x, type=1, na.rm=TRUE)
			})[3]
			p2 <- p
			p <- as.integer64(p2)
			names(p) <- names(p2)
			
			x <- as.integer64(x)

			tim["quantile.64","both"] <- timefun({
				p2 <- quantile(x, na.rm=TRUE)
			})[3]
			stopifnot(identical(p2, p))
				
			tim["sortqtl","prep"] <- timefun({
				s <- clone(x)
				na.count <- ramsort(s, na.last=FALSE)
			})[3]
			tim["sortqtl","use"] <- timefun({
				p2 <- sortqtl(s, na.count, seq(0, 1, 0.25))
			})[3]
			stopifnot(identical(unname(p2), unname(p)))
				
			tim["orderqtl","prep"] <- timefun({
				o <- seq_along(x)
				na.count <- ramorder(x, o, na.last=FALSE)
			})[3]
			tim["orderqtl","use"] <- timefun({
				p2 <- orderqtl(x, o, na.count, seq(0, 1, 0.25))
			})[3]
			stopifnot(identical(unname(p2), unname(p)))
				
			sortordercache(x)
			tim["sort.cache","use"] <- timefun({
				p2 <- quantile(x, na.rm=TRUE)
			})[3]
			stopifnot(identical(p2, p))
			remcache(x)

			ordercache(x)
			tim["order.cache","use"] <- timefun({
				p2 <- quantile(x, na.rm=TRUE)
			})[3]
			stopifnot(identical(p2, p))
			remcache(x)
			
			if (plot){
				barplot(t(tim), cex.names=0.7)
				title(paste("quantile(",n,")", sep=""))
			}
			
			ret[["quantile",as.character(n)]] <- tim
		}
	}

	ret
	
}


#! \name{match.integer64}
#! \alias{match.integer64}
#! \alias{\%in\%.integer64}
#! \title{
#! 64-bit integer matching
#! }
#! \description{
#! \code{match} returns a vector of the positions of (first) matches of its first argument in its second. 
#! 
#! \code{\%in\%} is a more intuitive interface as a binary operator, which returns a logical vector indicating if there is a match or not for its left operand. 
#! 
#! }
#! \usage{
#! \method{match}{integer64}(x, table, nomatch = NA_integer_, nunique = NULL, method = NULL, ...)
#! \method{\%in\%}{integer64}(x, table, ...)
#! }
#! \arguments{
#!   \item{x}{
#! 	integer64 vector: the values to be matched, optionally carrying a cache created with \code{\link{hashcache}}
#! }
#!   \item{table}{
#! 	integer64 vector: the values to be matched against, optionally carrying a cache created with \code{\link{hashcache}} or \code{\link{sortordercache}}
#! }
#!   \item{nomatch}{
#!   the value to be returned in the case when no match is found. Note that it is coerced to integer.
#! }
#!   \item{nunique}{
#! 	NULL or the number of unique values of table (including NA). Providing \code{nunique} can speed-up matching when \code{table} has no cache. Note that a wrong nunique can cause undefined behaviour up to a crash.
#! }
#!   \item{method}{
#! 	NULL for automatic method selection or a suitable low-level method, see details
#! }
#!   \item{\dots}{
#! ignored
#! }
#! }
#! \details{
#!   These functions automatically choose from several low-level functions considering the size of \code{x} and \code{table} and the availability of caches. 
#! 
#! 
#!   Suitable methods for \code{\%in\%.integer64} are \code{\link{hashpos}} (hash table lookup), \code{\link{hashrev}} (reverse lookup), \code{\link{sortorderpos}} (fast ordering) and \code{\link{orderpos}} (memory saving ordering).
#!   Suitable methods for \code{match.integer64} are \code{\link{hashfin}} (hash table lookup), \code{\link{hashrin}} (reverse lookup), \code{\link{sortfin}} (fast sorting) and \code{\link{orderfin}} (memory saving ordering).
#! }
#! \value{
#!   A vector of the same length as \code{x}.
#! 
#!   \code{match}: An integer vector giving the position in \code{table} of
#!   the first match if there is a match, otherwise \code{nomatch}.
#! 
#!   If \code{x[i]} is found to equal \code{table[j]} then the value
#!   returned in the \code{i}-th position of the return value is \code{j},
#!   for the smallest possible \code{j}.  If no match is found, the value
#!   is \code{nomatch}.
#! 
#!   \code{\%in\%}: A logical vector, indicating if a match was located for
#!   each element of \code{x}: thus the values are \code{TRUE} or
#!   \code{FALSE} and never \code{NA}.
#! }
#! \author{
#! 	Jens Oehlschlägel <Jens.Oehlschlaegel@truecluster.com>
#! }
#! \seealso{
#! 	\code{\link{match}}
#! }
#! \examples{
#! x <- as.integer64(c(NA, 0:9), 32)
#! table <- as.integer64(c(1:9, NA))
#! match.integer64(x, table)
#! "\%in\%.integer64"(x, table)
#!
#! x <- as.integer64(sample(c(rep(NA, 9), 0:9), 32, TRUE))
#! table <- as.integer64(sample(c(rep(NA, 9), 1:9), 32, TRUE))
#! stopifnot(identical(match.integer64(x, table), match(as.integer(x), as.integer(table))))
#! stopifnot(identical("\%in\%.integer64"(x, table), as.integer(x) \%in\% as.integer(table)))
#! 
#! \dontrun{
#! 	message("check when reverse hash-lookup beats standard hash-lookup")
#! 	e <- 4:24
#! 	timx <- timy <- matrix(NA, length(e), length(e), dimnames=list(e,e))
#! 	for (iy in seq_along(e))
#! 	for (ix in 1:iy){
#! 		nx <- 2^e[ix]
#! 		ny <- 2^e[iy]
#! 		x <- as.integer64(sample(ny, nx, FALSE))
#! 		y <- as.integer64(sample(ny, ny, FALSE))
#! 		#hashfun(x, bits=as.integer(5))
#! 		timx[ix,iy] <- repeat.time({
#! 		hx <- hashmap(x)
#! 		py <- hashrev(hx, y)
#! 		})[3]
#! 		timy[ix,iy] <- repeat.time({
#! 		hy <- hashmap(y)
#! 		px <- hashpos(hy, x)
#! 		})[3]
#! 		#identical(px, py)
#! 		print(round(timx[1:iy,1:iy]/timy[1:iy,1:iy], 2), na.print="")
#! 	}
#!
#! 	message("explore best low-level method given size of x and table")
#! 	B1 <- 1:27
#! 	B2 <- 1:27
#! 	tim <- array(NA, dim=c(length(B1), length(B2), 5)
#!  , dimnames=list(B1, B2, c("hashpos","hashrev","sortpos1","sortpos2","sortpos3")))
#! 	for (i1 in B1)
#! 	for (i2 in B2)
#! 	{
#! 	  b1 <- B1[i1]
#! 	  b2 <- B1[i2]
#! 	  n1 <- 2^b1
#! 	  n2 <- 2^b2
#! 	  x1 <- as.integer64(c(sample(n2, n1-1, TRUE), NA))
#! 	  x2 <- as.integer64(c(sample(n2, n2-1, TRUE), NA))
#! 	  tim[i1,i2,1] <- repeat.time({h <- hashmap(x2);hashpos(h, x1);rm(h)})[3]
#! 	  tim[i1,i2,2] <- repeat.time({h <- hashmap(x1);hashrev(h, x2);rm(h)})[3]
#! 	  s <- clone(x2); o <- seq_along(s); ramsortorder(s, o)
#! 	  tim[i1,i2,3] <- repeat.time(sortorderpos(s, o, x1, method=1))[3]
#! 	  tim[i1,i2,4] <- repeat.time(sortorderpos(s, o, x1, method=2))[3]
#! 	  tim[i1,i2,5] <- repeat.time(sortorderpos(s, o, x1, method=3))[3]
#! 	  rm(s,o)
#! 	  print(apply(tim, 1:2, function(ti)if(any(is.na(ti)))NA else which.min(ti)))
#! 	}
#! }
#! }
#! \keyword{manip}
#! \keyword{logic}


match.integer64 <- function(x, table, nomatch = NA_integer_, nunique=NULL, method=NULL, ...){
  stopifnot(is.integer64(x) &&  is.integer64(table))  # xx TODO
  c <- cache(table)
  if (is.null(method)){
    if (is.null(c)){
			nx <- length(x)
			if (is.null(nunique))
				nunique <- length(table)
			btable <- as.integer(ceiling(log2(nunique*1.5)))
			bx <- as.integer(ceiling(log2(nx*1.5)))
			if (bx<=17 && btable>=16){
				method <- "hashrev"
			}else{
				method <- "hashpos"
			}
	}else{
		if (exists("hashmap", envir=c, inherits=FALSE)){
			method <- "hashpos"
		}else if (exists("sort", envir=c, inherits=FALSE) && exists("order", envir=c, inherits=FALSE) && (length(table)>length(x) || length(x)<4096)){
			method <- "sortorderpos"
		}else if (exists("order", envir=c, inherits=FALSE) && (length(table)>length(x) || length(x)<4096)){
			method <- "orderpos"
		}else{
			nx <- length(x)
			if (is.null(nunique)){
			  if (exists("nunique", envir=c, inherits=FALSE))
				nunique <- c$nunique
			  else
				nunique <- length(table)
			}
			btable <- as.integer(ceiling(log2(nunique*1.5)))
			bx <- as.integer(ceiling(log2(nx*1.5)))
			if (bx<=17 && btable>=16){
				method <- "hashrev"
			}else{
				method <- "hashpos"
			}
		}
	}
  }
  switch(method
  , hashpos={
			if (is.null(c) || !exists("hashmap", envir=c, inherits=FALSE)){
				if (exists("btable", inherits=FALSE))
					h <- hashmap(table, hashbits=btable)
				else{
					if (is.null(nunique))
						nunique <- c$nunique
					h <- hashmap(table, nunique=nunique)
				}
			}else
				h <- c
			p <- hashpos(h, x, nomatch=nomatch)
    }
  , hashrev={
		c <- cache(x)
		if (is.null(c) || !exists("hashmap", envir=c, inherits=FALSE)){
				if (exists("bx", inherits=FALSE))
					h <- hashmap(x, bits=bx)
				else{
					if (is.null(nunique))
						nunique <- c$nunique
					h <- hashmap(x, nunique=nunique)
				}
			}else
				h <- c
		p <- hashrev(h, table, nomatch=nomatch)
    }
  , sortorderpos={
		if (is.null(c) || !exists("sort", c) || !exists("order", c)){
			s <- clone(table)
			o <- seq_along(s)
			ramsortorder(s, o, na.last=FALSE)
		}else{
			s <- get("sort", c)
			o <- get("order", c)
		}
		p <- sortorderpos(s, o, x, nomatch=nomatch)
    }
  , orderpos={
		if (is.null(c) || !exists("order", c)){
			o <- seq_along(s)
			ramorder(table, o, na.last=FALSE)
		}else{
			o <- get("order", c)
		}
		p <- orderpos(table, o, x, nomatch=nomatch)
    }
  , stop("unknown method")
  )
  p
}


"%in%.integer64" <- function(x, table, ...){
  stopifnot(is.integer64(x) &&  is.integer64(table))  # xx TODO
	nunique <- NULL
	method <- NULL
  c <- cache(table)
  if (is.null(method)){
    if (is.null(c)){
			nx <- length(x)
			if (is.null(nunique))
				nunique <- length(table)
			btable <- as.integer(ceiling(log2(nunique*1.5)))
			bx <- as.integer(ceiling(log2(nx*1.5)))
			if (bx<=17 && btable>=16){
				method <- "hashrin"
			}else{
				method <- "hashfin"
			}
	}else{
		if (exists("hashmap", envir=c, inherits=FALSE)){
			method <- "hashfin"
		}else if (exists("sort", envir=c, inherits=FALSE) && (length(table)>length(x) || length(x)<4096)){
			method <- "sortfin"
		}else if (exists("order", envir=c, inherits=FALSE) && (length(table)>length(x) || length(x)<4096)){
			method <- "orderfin"
		}else{
			nx <- length(x)
			if (is.null(nunique)){
			  if (exists("nunique", envir=c, inherits=FALSE))
				nunique <- c$nunique
			  else
				nunique <- length(table)
			}
			btable <- as.integer(ceiling(log2(nunique*1.5)))
			bx <- as.integer(ceiling(log2(nx*1.5)))
			if (bx<=17 && btable>=16){
				method <- "hashrin"
			}else{
				method <- "hashfin"
			}
		}
	}
  }
  switch(method
  , hashfin={
		if (is.null(c) || !exists("hashmap", envir=c, inherits=FALSE)){
			if (exists("btable", inherits=FALSE))
				h <- hashmap(table, hashbits=btable)
			else{
				if (is.null(nunique))
					nunique <- c$nunique
				h <- hashmap(table, nunique=nunique)
			}
		}else
			h <- c
		p <- hashfin(h, x)
    }
  , hashrin={
		c <- cache(x)
		if (is.null(c) || !exists("hashmap", envir=c, inherits=FALSE)){
				if (exists("bx", inherits=FALSE))
					h <- hashmap(x, bits=bx)
				else{
					if (is.null(nunique))
						nunique <- c$nunique
					h <- hashmap(x, nunique=nunique)
				}
		}else
			h <- c
		p <- hashrin(h, table)
    }
  , sortfin={
		if (is.null(c) || !exists("sort", c)){
			s <- clone(table)
			ramsort(s, na.last=FALSE)
		}else{
			s <- get("sort", c)
		}
		p <- sortfin(s, x)
    }
  , orderfin={
		if (is.null(c) || !exists("order", c)){
			o <- seq_along(s)
			ramorder(table, o, na.last=FALSE)
		}else{
			o <- get("order", c)
		}
		p <- orderfin(table, o, x)
    }
  , stop("unknown method")
  )
  p
}

#! \name{duplicated.integer64}
#! \alias{duplicated.integer64}
#! \title{Determine Duplicate Elements of integer64}
#! \description{
#!   \code{duplicated()} determines which elements of a vector or data frame are duplicates
#!   of elements with smaller subscripts, and returns a logical vector
#!   indicating which elements (rows) are duplicates.
#! }
#! \usage{
#! \method{duplicated}{integer64}(x, incomparables = FALSE, nunique = NULL, method = NULL, \dots)
#! }
#! \arguments{
#!   \item{x}{a vector or a data frame or an array or \code{NULL}.}
#!   \item{incomparables}{ignored}
#!   \item{nunique}{
#! 	NULL or the number of unique values (including NA). Providing \code{nunique} can speed-up matching when \code{x} has no cache. Note that a wrong nunique can cause undefined behaviour up to a crash.
#! }
#!   \item{method}{
#! 	NULL for automatic method selection or a suitable low-level method, see details
#! }
#!   \item{\dots}{ignored}
#! }
#! \details{
#!   This function automatically chooses from several low-level functions considering the size of \code{x} and the availability of a cache. 
#! 
#!   Suitable methods are \code{\link{hashdup}} (hashing), \code{\link{sortorderdup}} (fast ordering) and \code{\link{orderdup}} (memory saving ordering).
#! }
#! \value{
#!     \code{duplicated()}: a logical vector of the same length as \code{x}.  
#! }
#! \author{
#! 	Jens Oehlschlägel <Jens.Oehlschlaegel@truecluster.com>
#! }
#! \seealso{ \code{\link{duplicated}}, \code{\link{unique.integer64}}  }
#! \examples{
#! x <- as.integer64(sample(c(rep(NA, 9), 1:9), 32, TRUE))
#! duplicated(x)
#! 
#! stopifnot(identical(duplicated(x),  duplicated(as.integer(x))))
#! }
#! \keyword{logic}
#! \keyword{manip}
#! 

duplicated.integer64 <- function(x
, incomparables = FALSE  # dummy parameter
, nunique = NULL
, method = NULL
, ...
){
  stopifnot(identical(incomparables, FALSE))
  c <- cache(x)
  if (is.null(nunique) && !is.null(c))
	nunique <- c$nunique
  if (is.null(method)){
    if (is.null(c)){
		if (length(x)>5e7)
			method <- "sortorderdup"
		else
			method <- "hashdup"
	}else{
		if (exists("sort", envir=c, inherits=FALSE) && exists("order", envir=c, inherits=FALSE))
			method <- "sortorderdup"
		else if (exists("hashmap", envir=c, inherits=FALSE))
			method <- "hashdup"
		else if (exists("order", envir=c, inherits=FALSE))
			method <- "orderdup"
		else if (length(x)>5e7)
			method <- "sortorderdup"
		else
			method <- "hashdup"
	}
  }
  switch(method
  , hashdup={
		if (is.null(c) || !exists("hashmap", envir=c, inherits=FALSE))
			h <- hashmap(x, nunique=nunique)
		else
			h <- c
		p <- hashdup(h)
    }
  , sortorderdup={
		if (is.null(c) || !exists("sort", c, inherits=FALSE) || !exists("order", c, inherits=FALSE)){
			s <- clone(x)
			o <- seq_along(s)
			ramsortorder(s, o, na.last=FALSE)
		}else{
			s <- get("sort", c, inherits=FALSE)
			o <- get("order", c, inherits=FALSE)
		}
		p <- sortorderdup(s, o)
    }
  , orderdup={
		if (is.null(c) || !exists("order", c, inherits=FALSE)){
			o <- seq_along(s)
			ramorder(x, o, na.last=FALSE)
		}else{
			o <- get("order", c, inherits=FALSE)
		}
		p <- orderdup(x, o)
    }
  , stop("unknown method", method)
  )
  p
}


#! \name{unique.integer64}
#! \alias{unique.integer64}
#! \title{Extract Unique Elements from integer64}
#! \description{
#!   \code{unique} returns a vector like \code{x} but with duplicate elements/rows removed.
#! }
#! \usage{
#! \method{unique}{integer64}(x, incomparables = FALSE, order = c("original","values","any")
#! , nunique = NULL, method = NULL, \dots)
#! }
#! \arguments{
#!   \item{x}{a vector or a data frame or an array or \code{NULL}.}
#!   \item{incomparables}{ignored}
#!   \item{order}{The order in which unique values will be returned, see details}
#!   \item{nunique}{
#! 	NULL or the number of unique values (including NA). Providing \code{nunique} can speed-up matching when \code{x} has no cache. Note that a wrong nunique can cause undefined behaviour up to a crash.
#! }
#!   \item{method}{
#! 	NULL for automatic method selection or a suitable low-level method, see details
#! }
#!   \item{\dots}{ignored}
#! }
#! \details{
#!   This function automatically chooses from several low-level functions considering the size of \code{x} and the availability of a cache. 
#!   Suitable methods are \code{\link{hashmapuni}} (simultaneously creating and using a hashmap)
#! , \code{\link{hashuni}} (first creating a hashmap then using it)
#! , \code{\link{sortuni}} (fast sorting for sorted order only)
#! , \code{\link{sortorderuni}} (fast ordering for original order only) 
#! and \code{\link{orderuni}} (memory saving ordering).
#! \cr
#! The default \code{order="original"} returns unique values in the order of the first appearance in \code{x} like in \code{\link{unique}}, this costs extra processing. 
#! \code{order="values"} returns unique values in sorted order like in \code{\link{table}}, this costs extra processing with the hash methods but comes for free. 
#! \code{order="any"} returns unique values in undefined order, possibly faster. For hash methods this will be a quasi random order, for sort methods this will be sorted order.
#! }
#! \value{
#!   For a vector, an object of the same type of \code{x}, but with only
#!   one copy of each duplicated element.  No attributes are copied (so
#!   the result has no names).
#! }
#! \author{
#! 	Jens Oehlschlägel <Jens.Oehlschlaegel@truecluster.com>
#! }
#! \seealso{
#!   \code{\link{unique}} for the generic, \code{\link{unipos}} which gives the indices of the unique
#!   elements and \code{\link{table.integer64}} which gives frequencies of the unique elements.
#! }
#! \examples{
#! x <- as.integer64(sample(c(rep(NA, 9), 1:9), 32, TRUE))
#! unique(x)
#! unique(x, order="values")
#! 
#! stopifnot(identical(unique(x),  x[!duplicated(x)]))
#! stopifnot(identical(unique(x),  as.integer64(unique(as.integer(x)))))
#! stopifnot(identical(unique(x, order="values")
#! ,  as.integer64(sort(unique(as.integer(x)), na.last=FALSE))))
#! }
#! \keyword{manip}
#! \keyword{logic}


unique.integer64 <- function(x
, incomparables = FALSE  # dummy parameter
, order = c("original","values","any")
, nunique = NULL
, method = NULL
, ...
){
  stopifnot(identical(incomparables, FALSE))
  order <- match.arg(order)
  c <- cache(x)
  keep.order <- order == "original"
  if (is.null(nunique) && !is.null(c))
	nunique <- c$nunique
  if (is.null(method)){
    if (is.null(c)){
		if (order=="values")
			method <- "sortuni"
		else
			method <- "hashmapuni"
	}else{
		switch(order
		, "original" = {
			if (exists("hashmap", envir=c, inherits=FALSE))
				method <- "hashuni"
			else if (exists("order", envir=c, inherits=FALSE)){
				if (exists("sort", envir=c, inherits=FALSE))
					method <- "sortorderuni"
				else
					method <- "orderuni"
			}else
				method <- "hashmapuni"
		}
		, "values" = {
			if (exists("sort", envir=c, inherits=FALSE))
				method <- "sortuni"
			else if (exists("order", envir=c, inherits=FALSE))
				method <- "orderuni"
			else if (exists("hashmap", envir=c, inherits=FALSE) && c$nunique<length(x)/2)
				method <- "hashuni"
			else
				method <- "sortuni"
		}
		, "any" = {
			if (exists("sort", envir=c, inherits=FALSE))
				method <- "sortuni"
			else if (exists("hashmap", envir=c, inherits=FALSE))
				method <- "hashuni"
			else if (exists("order", envir=c, inherits=FALSE))
				method <- "orderuni"
			else
				method <- "sortuni"
		}
		)
	}
  }
  switch(method
  , hashmapuni={
		p <- hashmapuni(x, nunique=nunique)
    }
  , hashuni={
		if (is.null(c) || !exists("hashmap", envir=c, inherits=FALSE))
			h <- hashmap(x, nunique=nunique)
		else
			h <- c
		p <- hashuni(h, keep.order=keep.order)
		if (order=="values")
			ramsort(p, na.last=FALSE)
    }
  , sortuni={
		if (is.null(c) || !exists("sort", c, inherits=FALSE)){
			s <- clone(x)
			ramsort(s, na.last=FALSE)
		}else{
			s <- get("sort", c, inherits=FALSE)
		}
		if (is.null(nunique))
			nunique <- sortnut(s)[1]
		p <- sortuni(s, nunique)
    }
  , sortorderuni={
		if (is.null(c) || !exists("sort", c, inherits=FALSE) || !exists("order", c, inherits=FALSE)){
			s <- clone(x)
			o <- seq_along(x)
			ramsortorder(s, o, na.last=FALSE)
		}else{
			s <- get("sort", c, inherits=FALSE)
			o <- get("order", c, inherits=FALSE)
		}
		if (is.null(nunique))
			nunique <- sortnut(s)[1]
		p <- sortorderuni(x, s, o, nunique)
    }
  , orderuni={
		if (is.null(c) || !exists("order", c, inherits=FALSE)){
			o <- seq_along(x)
			ramorder(x, o, na.last=FALSE)
		}else{
			o <- get("order", c, inherits=FALSE)
		}
		if (is.null(nunique))
			nunique <- ordernut(x, o)[1]
		p <- orderuni(x, o, nunique, keep.order=keep.order)
    }
  , stop("unknown method", method)
  )
  p
}


#! \name{unipos}
#! \alias{unipos}
#! \alias{unipos.integer64}
#! \title{Extract Positions of Unique Elements}
#! \description{
#!   \code{unipos} returns the positions of those elements returned by \code{\link{unique}}.
#! }
#! \usage{
#! unipos(x, incomparables = FALSE, order = c("original","values","any"), \dots)
#! \method{unipos}{integer64}(x, incomparables = FALSE, order = c("original","values","any")
#! , nunique = NULL, method = NULL, \dots)
#! }
#! \arguments{
#!   \item{x}{a vector or a data frame or an array or \code{NULL}.}
#!   \item{incomparables}{ignored}
#!   \item{order}{The order in which positions of unique values will be returned, see details}
#!   \item{nunique}{
#! 	NULL or the number of unique values (including NA). Providing \code{nunique} can speed-up when \code{x} has no cache. Note that a wrong nunique can cause undefined behaviour up to a crash.
#! }
#!   \item{method}{
#! 	NULL for automatic method selection or a suitable low-level method, see details
#! }
#!   \item{\dots}{ignored}
#! }
#! \details{
#!   This function automatically chooses from several low-level functions considering the size of \code{x} and the availability of a cache. 
#!   Suitable methods are \code{\link{hashmapupo}} (simultaneously creating and using a hashmap)
#! , \code{\link{hashupo}} (first creating a hashmap then using it)
#! , \code{\link{sortorderupo}} (fast ordering) 
#! and \code{\link{orderupo}} (memory saving ordering).
#! \cr
#! The default \code{order="original"} collects unique values in the order of the first appearance in \code{x} like in \code{\link{unique}}, this costs extra processing. 
#! \code{order="values"} collects unique values in sorted order like in \code{\link{table}}, this costs extra processing with the hash methods but comes for free. 
#! \code{order="any"} collects unique values in undefined order, possibly faster. For hash methods this will be a quasi random order, for sort methods this will be sorted order.
#! }
#! \value{
#!   an integer vector of positions
#! }
#! \author{
#! 	Jens Oehlschlägel <Jens.Oehlschlaegel@truecluster.com>
#! }
#! \seealso{
#!   \code{\link{unique.integer64}} for unique values and \code{\link{match.integer64}} for general matching.
#! }
#! \examples{
#! x <- as.integer64(sample(c(rep(NA, 9), 1:9), 32, TRUE))
#! unipos(x)
#! unipos(x, order="values")
#! 
#! stopifnot(identical(unipos(x),  (1:length(x))[!duplicated(x)]))
#! stopifnot(identical(unipos(x),  match.integer64(unique(x), x)))
#! stopifnot(identical(unipos(x, order="values"),  match.integer64(unique(x, order="values"), x)))
#! stopifnot(identical(unique(x),  x[unipos(x)]))
#! stopifnot(identical(unique(x, order="values"),  x[unipos(x, order="values")]))
#! }
#! \keyword{manip}
#! \keyword{logic}


unipos <- function(x, incomparables = FALSE, order = c("original","values","any"), ...)UseMethod("unipos")
unipos.integer64 <- function(x
, incomparables = FALSE  # dummy parameter
, order = c("original","values","any")
, nunique = NULL
, method = NULL
, ...
){
  stopifnot(identical(incomparables, FALSE))
  order <- match.arg(order)
  c <- cache(x)
  keep.order <- order == "original"
  if (is.null(nunique) && !is.null(c))
	nunique <- c$nunique
  if (is.null(method)){
    if (is.null(c)){
		if (order=="values")
			method <- "sortorderupo"
		else
			method <- "hashmapupo"
	}else{
		switch(order
		, "original" = {
			if (exists("hashmap", envir=c, inherits=FALSE))
				method <- "hashupo"
			else if (exists("order", envir=c, inherits=FALSE)){
				if (exists("sort", envir=c, inherits=FALSE))
					method <- "sortorderupo"
				else
					method <- "orderupo"
			}else
				method <- "hashmapupo"
		}
		, "values" = {
			if (exists("order", envir=c, inherits=FALSE)){
				if (exists("sort", envir=c, inherits=FALSE))
					method <- "sortorderupo"
				else
					method <- "orderupo"
			}else if (exists("hashmap", envir=c, inherits=FALSE) && c$nunique<length(x)/2)
				method <- "hashupo"
			else
				method <- "sortorderupo"
		}
		, "any" = {
			if (exists("sort", envir=c, inherits=FALSE) && exists("order", envir=c, inherits=FALSE))
				method <- "sortorderupo"
			else if (exists("hashmap", envir=c, inherits=FALSE))
				method <- "hashupo"
			else if (exists("order", envir=c, inherits=FALSE))
				method <- "orderupo"
			else
				method <- "sortorderupo"
		}
		)
	}
  }
  switch(method
  , hashmapupo={
		p <- hashmapupo(x, nunique=nunique)
    }
  , hashupo={
		if (is.null(c) || !exists("hashmap", envir=c, inherits=FALSE))
			h <- hashmap(x, nunique=nunique)
		else
			h <- c
		p <- hashupo(h, keep.order=keep.order)
		if (order=="values"){
			s <- x[p]
			ramsortorder(s, p, na.last=FALSE)
		}
    }
  , sortorderupo={
		if (is.null(c) || !exists("sort", c, inherits=FALSE) || !exists("order", c, inherits=FALSE)){
			s <- clone(x)
			o <- seq_along(x)
			ramsortorder(s, o, na.last=FALSE)
		}else{
			s <- get("sort", c, inherits=FALSE)
			o <- get("order", c, inherits=FALSE)
		}
		if (is.null(nunique))
			nunique <- sortnut(s)[1]
		p <- sortorderupo(s, o, nunique, keep.order=keep.order)
    }
  , orderupo={
		if (is.null(c) || !exists("order", c, inherits=FALSE)){
			o <- seq_along(x)
			ramorder(x, o, na.last=FALSE)
		}else{
			o <- get("order", c, inherits=FALSE)
		}
		if (is.null(nunique))
			nunique <- ordernut(x, o)[1]
		p <- orderupo(x, o, nunique, keep.order=keep.order)
    }
  , stop("unknown method", method)
  )
  p
}



#! \name{table.integer64}
#! \title{Cross Tabulation and Table Creation for integer64}
#! \alias{table.integer64}
#! 
#! \concept{counts}
#! \concept{frequencies}
#! \concept{occurrences}
#! \concept{contingency table}
#! 
#! \description{
#!   \code{table.integer64} uses the cross-classifying integer64 vectors to build a contingency
#!   table of the counts at each combination of vector values.
#! }
#! \usage{
#! table.integer64(\dots
#! , return = c("table","data.frame","list")
#! , order = c("values","counts")
#! , nunique = NULL
#! , method = NULL
#! , dnn = list.names(...), deparse.level = 1
#! ) 
#! }
#! \arguments{
#!   \item{\dots}{one or more objects which can be interpreted as factors
#!     (including character strings), or a list (or data frame) whose
#!     components can be so interpreted.  (For \code{as.table} and
#!     \code{as.data.frame}, arguments passed to specific methods.)}
#!   \item{nunique}{
#! 	NULL or the number of unique values of table (including NA). Providing \code{nunique} can speed-up matching when \code{table} has no cache. Note that a wrong nunique can cause undefined behaviour up to a crash.
#! }
#!   \item{order}{
#! 	By default results are created sorted by "values", or by "counts"
#! }
#!   \item{method}{
#! 	NULL for automatic method selection or a suitable low-level method, see details
#! }
#!   \item{return}{
#!      choose the return format, see details
#! }
#!   \item{dnn}{the names to be given to the dimensions in the result (the
#!     \emph{dimnames names}).}
#!   \item{deparse.level}{controls how the default \code{dnn} is
#!     constructed.  See \sQuote{Details}.}
#! }
#! \details{
#!   This function automatically chooses from several low-level functions considering the size of \code{x} and the availability of a cache. 
#!   Suitable methods are \code{\link{hashmaptab}} (simultaneously creating and using a hashmap)
#! , \code{\link{hashtab}} (first creating a hashmap then using it)
#! , \code{\link{sortordertab}} (fast ordering) 
#! and \code{\link{ordertab}} (memory saving ordering).
#! \cr
#!   If the argument \code{dnn} is not supplied, the internal function
#!   \code{list.names} is called to compute the \sQuote{dimname names}.  If the
#!   arguments in \code{\dots} are named, those names are used.  For the
#!   remaining arguments, \code{deparse.level = 0} gives an empty name,
#!   \code{deparse.level = 1} uses the supplied argument if it is a symbol,
#!   and \code{deparse.level = 2} will deparse the argument.
#! 
#!   Arguments \code{exclude}, \code{useNA}, are not supported, i.e. \code{NA}s are always tabulated, and, different from \code{\link{table}} they are sorted first if \code{order="values"}. 
#! }
#! \value{
#!   By default (with \code{return="table"}) \code{\link{table}} returns a \emph{contingency table}, an object of
#!   class \code{"table"}, an array of integer values. Note that unlike S the result is always an array, a 1D array if one factor is given. Note also that for multidimensional arrays this is a \emph{dense} return structure which can dramatically increase RAM requirements (for large arrays with high mutual information, i.e. many possible input combinations of which only few occur) and that \code{\link{table}} is limited to \code{2^31} possible combinations (e.g. two input vectors with 46340 unique values only). Finally note that the tabulated values or value-combinations are represented as \code{dimnames} and that the implied conversion of values to strings can cause \emph{severe} performance problems since each string needs to be integrated into R's global string cache. 
#!   \cr
#!   You can use the other \code{return=} options to cope with these problems, the potential combination limit is increased from \code{2^31} to \code{2^63} with these options, RAM is only rewquired for observed combinations and string conversion is avoided. 
#!   \cr
#!   With \code{return="data.frame"} you get a \emph{dense} representation as a \code{\link{data.frame}} (like that resulting from \code{as.data.frame(table(...))}) where only observed combinations are listed (each as a data.frame row) with the corresponding frequency counts (the latter as component
#!   named by \code{responseName}).  This is the inverse of \code{\link{xtabs}}..
#!   \cr
#!   With \code{return="list"} you also get a \emph{dense} representation as a simple \code{\link{list}} with components 
#!   \item{values }{a integer64 vector of the technically tabulated values, for 1D this is the tabulated values themselves, for kD these are the values representing the potential combinations of input values}
#!   \item{counts}{the frequency counts}
#!   \item{dims}{only for kD: a list with the vectors of the unique values of the input dimensions}
#! }
#! \note{
#!   Note that by using \code{\link{as.integer64.factor}} we can also input 
#!   factors into \code{table.integer64} -- only the \code{\link{levels}} get lost.
#!  \cr
#!   Note that because of the existence of \code{\link{as.factor.integer64}} 
#! the standard \code{\link{table}} function -- within its limits -- can also be used 
#! for \code{\link{integer64}}, and especially for combining \code{\link{integer64}} input 
#! with other data types.
#! }
#! \seealso{
#!   \code{\link{table}} for more info on the standard version coping with Base R's data types, \code{\link{tabulate}} which can faster tabulate \code{\link{integer}s} with a limited range \code{[1L .. nL not too big]}, \code{\link{unique.integer64}} for the unique values without counting them and \code{\link{unipos.integer64}} for the positions of the unique values. 
#! }
#! \examples{
#! message("pure integer64 examples")
#! x <- as.integer64(sample(c(rep(NA, 9), 1:9), 32, TRUE))
#! y <- as.integer64(sample(c(rep(NA, 9), 1:9), 32, TRUE))
#! z <- sample(c(rep(NA, 9), letters), 32, TRUE)
#! table.integer64(x)
#! table.integer64(x, order="counts")
#! table.integer64(x, y)
#! table.integer64(x, y, return="data.frame")
#!
#! message("via as.integer64.factor we can use 'table.integer64' also for factors")
#! table.integer64(x, as.integer64(as.factor(z)))
#! 
#! message("via as.factor.integer64 we can also use 'table' for integer64")
#! table(x)
#! table(x, exclude=NULL)
#! table(x, z, exclude=NULL)
#!
#! \dontshow{
#!  stopifnot(identical(table.integer64(as.integer64(c(1,1,2))), table(c(1,1,2))))
#!  stopifnot(identical(table.integer64(as.integer64(c(1,1,2)),as.integer64(c(3,4,4))), table(c(1,1,2),c(3,4,4))))
#!  message("the following works with three warnings due to coercion")
#!  stopifnot(identical(table.integer64(c(1,1,2)), table(c(1,1,2))))
#!  stopifnot(identical(table.integer64(as.integer64(c(1,1,2)),c(3,4,4)), table(c(1,1,2),c(3,4,4))))
#!  stopifnot(identical(table.integer64(c(1,1,2),as.integer64(c(3,4,4))), table(c(1,1,2),c(3,4,4))))
#!  message("the following works because of as.factor.integer64")
#!  stopifnot(identical(table(as.integer64(c(1,1,2))), table(c(1,1,2))))  
#!  stopifnot(identical(table(as.integer64(c(1,1,2)),as.integer64(c(3,4,4))), table(c(1,1,2),c(3,4,4))))
#!  stopifnot(identical(table(as.integer64(c(1,1,2)),c(3,4,4)), table(c(1,1,2),c(3,4,4))))
#!  stopifnot(identical(table(c(1,1,2),as.integer64(c(3,4,4))), table(c(1,1,2),c(3,4,4))))
#! }
#!
#! }
#! \keyword{category}


table.integer64 <- function(
  ...
, return = c("table","data.frame","list")
, order = c("values","counts")
, nunique = NULL
, method = NULL
, dnn = list.names(...), deparse.level = 1
){
  order <- match.arg(order)
  return <- match.arg(return)
  # this is taken from 'table'
	list.names <- function(...){
        l <- as.list(substitute(list(...)))[-1L]
        nm <- names(l)
        fixup <- if (is.null(nm)) 
            seq_along(l)
        else nm == ""
        dep <- vapply(l[fixup], function(x) switch(deparse.level + 
            1, "", if (is.symbol(x)) as.character(x) else "", 
            deparse(x, nlines = 1)[1L]), "")
        if (is.null(nm)) 
            dep
        else {
            nm[fixup] <- dep
            nm
        }
    }  

	# COPY ON MODIFY is broken for reading from list(...)
	# because list(...) creates a copy of all ... and this invalidates our caches
	# therefore we go this sick workaround
	argsymbols <- as.list(substitute(list(...)))[-1L]
	argframe <- parent.frame()
	A <- function(i)eval(argsymbols[[i]], argframe)
	N <- length(argsymbols)
	if (!N) 
		stop("nothing to tabulate")
	if (N == 1L && is.list(A(1L))){
		args <- A(1L)
		if (length(dnn) != length(args)) 
            dnn <- if (!is.null(argn <- names(args))) argn
				else paste(dnn[1L], seq_along(args), sep = ".")		
		N <- length(args)
		A <- function(i)args[[i]]
	}
	force(dnn)
		
	if (N==1L){
		x <- A(1L)
			if (!is.integer64(x)){
				warning("coercing first argument to integer64")
				x <- as.integer64(x)
			}
	}else{
		a <- A(1L)
		n <- length(a)
		nu <- integer(N)
		d <- integer64(N+1); d[[1]] <- 1L
		dims <- vector("list", N)
		names(dims) <- dnn
		for (i in 1:N){
			a <- A(i)
			if (length(a) != n) 
				stop("all input vectors must have the same length")
			if (!is.integer64(a)){
				warning("coercing argument ", i, " to integer64")
				a <- as.integer64(a)
			}
			c <- cache(a)
			if (is.null(c$order)){
				s <- clone(a)
				o <- seq_along(s)
				ramsortorder(s,o)
				nu[[i]] <- sortnut(s)[["nunique"]]
			}else if (is.null(c$sort)){
				o <- c$order
				s <- a[o]
				nu[[i]] <- c$nunique
			}else{
				o <- c$order
				s <- c$sort
				nu[[i]] <- c$nunique
			}
			d[[i+1]] <- d[[i]] * nu[[i]]
			if (is.na(d[[i+1]]))
				stop("attempt to make a table from more than >= 2^63 hypothetical combinations")
			dims[[i]] <- sortuni(s, nu[[i]])
			if (i==1L)
				x <- sortorderkey(s,o) - 1L
			else
				x <- x + d[[i]] * (sortorderkey(s,o) - 1L)
		}
	}
  c <- cache(x)
  if (is.null(nunique) && !is.null(c))
	nunique <- c$nunique
  if (is.null(method)){
    if (is.null(c)){
		if (order=="values" && (is.null(nunique) || nunique>65536))
			method <- "sorttab"
		else
			method <- "hashmaptab"
	}else{
		if (order=="values"){
			if (exists("sort", envir=c, inherits=FALSE))
				method <- "sorttab"
			else if (exists("hashmap", envir=c, inherits=FALSE) && c$nunique<sqrt(length(x)))
				method <- "hashtab"
			else if (exists("order", envir=c, inherits=FALSE))
				method <- "ordertab"
			else
				method <- "sorttab"
		}else{ # order = "counts"
			if (exists("hashmap", envir=c, inherits=FALSE))
				method <- "hashtab"
			else if (exists("sort", envir=c, inherits=FALSE))
				method <- "sorttab"
			else if (exists("order", envir=c, inherits=FALSE))
				method <- "ordertab"
			else
				method <- "hashmaptab"
		}
	}
  }
  switch(method
  , hashmaptab={
		tmp <- hashmaptab(x, nunique=nunique)
		cnt <- tmp$counts
		val <- tmp$values
		rm(tmp)
    }
  , hashtab={
		if (is.null(c) || !exists("hashmap", envir=c, inherits=FALSE))
			h <- hashmap(x, nunique=nunique)
		else
			h <- c 
		tmp <- hashtab(h, keep.order=FALSE)
		cnt <- tmp$counts
		val <- tmp$values
		rm(tmp)
    }
  , sorttab={
		if (is.null(c) || !exists("sort", c, inherits=FALSE)){
			s <- clone(x)
			ramsort(s, na.last=FALSE)
		}else{
			s <- get("sort", c, inherits=FALSE)
		}
		if (is.null(nunique))
			nunique <- sortnut(s)[1]
		val <- sortuni(s, nunique)
		cnt <- sorttab(s, nunique)
    }
  , ordertab={
		if (is.null(c) || !exists("order", c, inherits=FALSE)){
			o <- seq_along(x)
			ramorder(x, o, na.last=FALSE)
		}else{
			o <- get("order", c, inherits=FALSE)
		}
		if (is.null(nunique))
			nunique <- ordernut(x, o)[1]
		val <- orderuni(x, o, nunique, keep.order=FALSE)
		cnt <- ordertab(x, o, nunique, keep.order=FALSE)
		rm(o)
    }
  , stop("unknown method", method)
  )
	if (order=="values"){
		if (substr(method, 1, 4)=="hash"){
			o <- seq_along(val)
			ramsortorder(val, o, na.last=FALSE)
			cnt <- cnt[o]
		}
	}else{
		# xx workaround until we have implemented ramsort.integer
		o <- sort.list(cnt, na.last=NA, method="quick")
		cnt <- cnt[o]
		# o <- seq_along(cnt)
		# ramsortorder(cnt, o, na.last=FALSE)
		val <- val[o]
	}
  ## attaching names is extremely expensive with many unique values, doing this only for compatibility with 'table' here
  switch(return
  ,  "table" = {
  		if (N==1){
			attr(cnt, "dim") <-length(cnt)
			dn <- list(as.character(val))
			names(dn) <- dnn[1]
			attr(cnt, "dimnames") <- dn
		}else{
			a <- array(0L, dim=nu, dimnames=lapply(dims, as.character))
			a[as.integer(val)+1L] <- as.integer(cnt)
			cnt <- a
		}
		oldClass(cnt) <- "table"
	}
  ,  "data.frame" = {
		if (N==1){
			cnt <- data.frame(values=val, Freq=cnt)
			names(cnt)[[1]] <- dnn[1]
		}else{
			for (i in N:1){
				w <- val %/% d[[i]]
				val <- val - d[[i]]*w
				dims[[i]] <- dims[[i]][as.integer(w)+1L]
			}
			cnt <- data.frame(dims, Freq=cnt)
		}
	}
  , "list" = {
		if (N==1)
			cnt <- list(values=val, counts=cnt)
		else
			cnt <- list(values=val, counts=cnt, dims=dims)
    }
  )
  cnt
}


as.factor.integer64 <- function(x){

	c <- cache(x)
	if (is.null(c$order)){
		s <- clone(x)
		o <- seq_along(s)
		na.count <- ramsortorder(s,o)
		nu <- sortnut(s)[["nunique"]]
	}else if (is.null(c$sort)){
		o <- c$order
		s <- x[o]
		na.count <- c$na.count
		nu <- c$nunique
	}else{
		o <- c$order
		s <- c$sort
		na.count <- c$na.count
		nu <- c$nunique
	}
	dimtab <- sortuni(s, nu)
	dimpos <- sortorderkey(s,o,na.skip.num=na.count) - 1L
	attr(dimpos, "levels") <- dimtab
	oldClass(dimpos) <- "factor"
	dimpos
}

as.ordered.integer64 <- function(x){

	c <- cache(x)
	if (is.null(c$order)){
		s <- clone(x)
		o <- seq_along(s)
		na.count <- ramsortorder(s,o)
		nu <- sortnut(s)[["nunique"]]
	}else if (is.null(c$sort)){
		o <- c$order
		s <- x[o]
		na.count <- c$na.count
		nu <- c$nunique
	}else{
		o <- c$order
		s <- c$sort
		na.count <- c$na.count
		nu <- c$nunique
	}
	dimtab <- sortuni(s, nu)
	dimpos <- sortorderkey(s,o,na.skip.num=na.count) - 1L
	attr(dimpos, "levels") <- dimtab
	oldClass(dimpos) <- c("ordered", "factor")
	dimpos
}

as.integer64.factor <- function(x, ...)as.integer64(unclass(x))



#! \name{keypos}
#! \alias{keypos}
#! \alias{keypos.integer64}
#! \title{Extract Positions in redundant dimension table}
#! \description{
#!   \code{keypos} returns the positions of the (fact table) elements that participate in their sorted unique subset (dimension table)
#! }
#! \usage{
#! keypos(x, \dots)
#! \method{keypos}{integer64}(x, method = NULL, \dots)
#! }
#! \arguments{
#!   \item{x}{a vector or a data frame or an array or \code{NULL}.}
#!   \item{method}{
#! 	NULL for automatic method selection or a suitable low-level method, see details
#! }
#!   \item{\dots}{ignored}
#! }
#! \details{
#!   NAs are sorted first in the dimension table, see \code{\link{ramorder.integer64}}.
#!   \cr
#!   This function automatically chooses from several low-level functions considering the size of \code{x} and the availability of a cache. 
#!   Suitable methods are \code{\link{sortorderkey}} (fast ordering) 
#! and \code{\link{orderkey}} (memory saving ordering).
#! }
#! \value{
#!   an integer vector of the same length as code{x} containing positions relativ to code{sort(unique(x), na.last=FALSE)}
#! }
#! \author{
#! 	Jens Oehlschlägel <Jens.Oehlschlaegel@truecluster.com>
#! }
#! \seealso{
#!   \code{\link{unique.integer64}} for the unique subset and \code{\link{match.integer64}} for finding positions in a different vector.
#! }
#! \examples{
#! x <- as.integer64(sample(c(rep(NA, 9), 1:9), 32, TRUE))
#! keypos(x)
#! 
#! stopifnot(identical(keypos(x),  match.integer64(x, sort(unique(x), na.last=FALSE))))
#! }
#! \keyword{manip}
#! \keyword{univar}



keypos <- function(x, ...)UseMethod("keypos")
keypos.integer64 <- function(x
, method = NULL
, ...
){
  c <- cache(x)
  if (is.null(method)){
    if (is.null(c)){
		method <- "sortorderkey"
	}else{
		if (exists("order", envir=c, inherits=FALSE)){
			if (exists("sort", envir=c, inherits=FALSE))
				method <- "sortorderkey"
			else
				method <- "orderkey"
		}else
			method <- "sortorderkey"
	}
  }
  switch(method
  , sortorderkey={
		if (is.null(c) || !exists("sort", c, inherits=FALSE) || !exists("order", c, inherits=FALSE)){
			s <- clone(x)
			o <- seq_along(x)
			ramsortorder(s, o, na.last=FALSE)
		}else{
			s <- get("sort", c, inherits=FALSE)
			o <- get("order", c, inherits=FALSE)
		}
		p <- sortorderkey(s, o)
    }
  , orderkey={
		if (is.null(c) || !exists("order", c, inherits=FALSE)){
			o <- seq_along(x)
			ramorder(x, o, na.last=FALSE)
		}else{
			o <- get("order", c, inherits=FALSE)
		}
		p <- orderkey(x, o)
    }
  , stop("unknown method", method)
  )
  p
}

#! \name{tiepos}
#! \alias{tiepos}
#! \alias{tiepos.integer64}
#! \title{Extract Positions of Tied Elements}
#! \description{
#!   \code{tiepos} returns the positions of those elements that participate in ties.
#! }
#! \usage{
#! tiepos(x, \dots)
#! \method{tiepos}{integer64}(x, nties = NULL, method = NULL, \dots)
#! }
#! \arguments{
#!   \item{x}{a vector or a data frame or an array or \code{NULL}.}
#!   \item{nties}{
#! 	NULL or the number of tied values (including NA). Providing \code{nties} can speed-up when \code{x} has no cache. Note that a wrong nties can cause undefined behaviour up to a crash.
#! }
#!   \item{method}{
#! 	NULL for automatic method selection or a suitable low-level method, see details
#! }
#!   \item{\dots}{ignored}
#! }
#! \details{
#!   This function automatically chooses from several low-level functions considering the size of \code{x} and the availability of a cache. 
#!   Suitable methods are \code{\link{sortordertie}} (fast ordering) 
#! and \code{\link{ordertie}} (memory saving ordering).
#! }
#! \value{
#!   an integer vector of positions
#! }
#! \author{
#! 	Jens Oehlschlägel <Jens.Oehlschlaegel@truecluster.com>
#! }
#! \seealso{
#!   \code{\link{rank.integer64}} for possibly tied ranks and \code{\link{unipos.integer64}} for positions of unique values.
#! }
#! \examples{
#! x <- as.integer64(sample(c(rep(NA, 9), 1:9), 32, TRUE))
#! tiepos(x)
#! 
#! stopifnot(identical(tiepos(x),  (1:length(x))[duplicated(x) | rev(duplicated(rev(x)))]))
#! }
#! \keyword{manip}
#! \keyword{univar}


tiepos <- function(x, ...)UseMethod("tiepos")
tiepos.integer64 <- function(x
, nties = NULL
, method = NULL
, ...
){
  c <- cache(x)
  if (is.null(nties) && !is.null(c))
	nties <- c$nties
  if (is.null(method)){
    if (is.null(c)){
		method <- "sortordertie"
	}else{
		if (exists("order", envir=c, inherits=FALSE)){
			if (exists("sort", envir=c, inherits=FALSE))
				method <- "sortordertie"
			else
				method <- "ordertie"
		}else
			method <- "sortordertie"
	}
  }
  switch(method
  , sortordertie={
		if (is.null(c) || !exists("sort", c, inherits=FALSE) || !exists("order", c, inherits=FALSE)){
			s <- clone(x)
			o <- seq_along(x)
			ramsortorder(s, o, na.last=FALSE)
		}else{
			s <- get("sort", c, inherits=FALSE)
			o <- get("order", c, inherits=FALSE)
		}
		if (is.null(nties))
			nties <- sortnut(s)[2]
		p <- sortordertie(s, o, nties)
    }
  , ordertie={
		if (is.null(c) || !exists("order", c, inherits=FALSE)){
			o <- seq_along(x)
			ramorder(x, o, na.last=FALSE)
		}else{
			o <- get("order", c, inherits=FALSE)
		}
		if (is.null(nties))
			nties <- ordernut(x, o)[2]
		p <- ordertie(x, o, nties)
    }
  , stop("unknown method", method)
  )
  p
}


#! \name{rank.integer64}
#! \alias{rank.integer64}
#! \title{Sample Ranks from integer64}
#! \description{
#!   Returns the sample ranks of the values in a vector.  Ties (i.e., equal
#!   values) are averaged and missing values propagated.
#! }
#! \usage{
#! 	\method{rank}{integer64}(x, method = NULL, \dots)
#! }
#! \arguments{
#!   \item{x}{a integer64 vector}
#!   \item{method}{
#! 	NULL for automatic method selection or a suitable low-level method, see details
#! }
#!   \item{\dots}{ignored}
#! }
#! \details{
#!   This function automatically chooses from several low-level functions considering the size of \code{x} and the availability of a cache. 
#!   Suitable methods are \code{\link{sortorderrnk}} (fast ordering) 
#! and \code{\link{orderrnk}} (memory saving ordering).
#! }
#! \value{
#!   A numeric vector of the same length as \code{x}.
#! }
#! \author{
#! 	Jens Oehlschlägel <Jens.Oehlschlaegel@truecluster.com>
#! }
#! \seealso{
#!   \code{\link{order.integer64}}, \code{\link{rank}} and \code{\link{prank}} for percent rank.
#! }
#! \examples{
#! x <- as.integer64(sample(c(rep(NA, 9), 1:9), 32, TRUE))
#! rank.integer64(x)
#! 
#! stopifnot(identical(rank.integer64(x),  rank(as.integer(x)
#! , na.last="keep", ties.method = "average")))
#! }
#! \keyword{univar}

rank.integer64 <- function(x
, method = NULL
, ...
){
  c <- cache(x)
  if (is.null(method)){
    if (is.null(c)){
		method <- "sortorderrnk"
	}else{
		if (exists("order", envir=c, inherits=FALSE)){
			if (exists("sort", envir=c, inherits=FALSE))
				method <- "sortorderrnk"
			else
				method <- "orderrnk"
		}else
			method <- "sortorderrnk"
	}
  }
  switch(method
  , sortorderrnk={
		if (is.null(c) || !exists("sort", c, inherits=FALSE) || !exists("order", c, inherits=FALSE)){
			s <- clone(x)
			o <- seq_along(x)
			na.count <- ramsortorder(s, o, na.last=FALSE)
		}else{
			s <- get("sort", c, inherits=FALSE)
			o <- get("order", c, inherits=FALSE)
			na.count <- get("na.count", c, inherits=FALSE)
		}
		p <- sortorderrnk(s, o, na.count)
    }
  , orderrnk={
		if (is.null(c) || !exists("order", c, inherits=FALSE)){
			o <- seq_along(x)
			na.count <- ramorder(x, o, na.last=FALSE)
		}else{
			o <- get("order", c, inherits=FALSE)
			na.count <- get("na.count", c, inherits=FALSE)
		}
		p <- orderrnk(x, o, na.count)
    }
  , stop("unknown method", method)
  )
  p
}

#! \name{prank}
#! \alias{prank}
#! \alias{prank.integer64}
#! \title{(P)ercent (Rank)s}
#! \description{
#! 	Function \code{prank.integer64}  projects the values [min..max] via ranks [1..n] to [0..1]. 
#! 	\code{\link{qtile.integer64}} is the inverse function of 'prank.integer64' and projects [0..1] to [min..max].
#! }
#! \usage{
#! 	prank(x, \dots)
#! 	\method{prank}{integer64}(x, method = NULL, \dots)
#! }
#! \arguments{
#!   \item{x}{a integer64 vector}
#!   \item{method}{
#! 	NULL for automatic method selection or a suitable low-level method, see details
#! }
#!   \item{\dots}{ignored}
#! }
#! \details{
#! 	Function \code{prank.integer64} is based on \code{\link{rank.integer64}}.
#! }
#! \value{
#!   \code{prank} returns a numeric vector of the same length as \code{x}.
#! }
#! \author{
#! 	Jens Oehlschlägel <Jens.Oehlschlaegel@truecluster.com>
#! }
#! \seealso{
#!   \code{\link{rank.integer64}} for simple ranks and \code{\link{qtile}} for the inverse function quantiles.
#! }
#! \examples{
#! x <- as.integer64(sample(c(rep(NA, 9), 1:9), 32, TRUE))
#! prank(x)
#! 
#! x <- x[!is.na(x)]
#! stopifnot(identical(x,  unname(qtile(x, probs=prank(x)))))
#! }
#! \keyword{univar}

prank <- function(x, ...)UseMethod("prank")
prank.integer64 <- function(x
, method = NULL
, ...
)
{	
	n <- nvalid(x)
	if (n<2L)
		return(rep(as.integer64(NA), length(x)))
	(rank.integer64(x, method=method, ...)-1L)/(n-1L)
}

#! \name{qtile}
#! \alias{qtile}
#! \alias{qtile.integer64}
#! \alias{quantile.integer64}
#! \alias{median.integer64}
#! \alias{mean.integer64}
#! \alias{summary.integer64}
#! \title{(Q)uan(Tile)s }
#! \description{
#! 	Function \code{\link{prank.integer64}}  projects the values [min..max] via ranks [1..n] to [0..1]. 
#! 	\code{qtile.ineger64} is the inverse function of 'prank.integer64' and projects [0..1] to [min..max].
#! }
#! \usage{
#! 	qtile(x, probs=seq(0, 1, 0.25), \dots)
#! 	\method{qtile}{integer64}(x, probs = seq(0, 1, 0.25), names = TRUE, method = NULL, \dots)
#! 	\method{quantile}{integer64}(x, probs = seq(0, 1, 0.25), na.rm = FALSE, names = TRUE, type=0L, \dots)
#! 	\method{median}{integer64}(x, na.rm = FALSE)
#!  \method{mean}{integer64}(x, na.rm = FALSE, \dots)
#! 	\method{summary}{integer64}(object, \dots)
#! }
#! \arguments{
#!   \item{x}{a integer64 vector}
#!   \item{object}{a integer64 vector}
#!   \item{probs}{
#! 		numeric vector of probabilities with values in [0,1] - possibly containing \code{NA}s
#! }
#!   \item{names}{
#! 	logical; if \code{TRUE}, the result has a \code{names} attribute. Set to \code{FALSE} for speedup with many probs.
#! }
#!   \item{type}{
#! 	an integer selecting the quantile algorithm, currently only 0 is supported, see details
#! }
#!   \item{method}{
#! 	NULL for automatic method selection or a suitable low-level method, see details
#! }
#!   \item{na.rm}{
#! 	logical; if \code{TRUE}, any \code{NA} and \code{NaN}'s are removed from \code{x} before the quantiles are computed.
#! }
#!   \item{\dots}{ignored}
#! }
#! \details{
#!  Functions \code{quantile.integer64} with \code{type=0} and \code{median.integer64} are convenience wrappers to \code{qtile}.
#!  \cr
#!	Function \code{qtile} behaves very similar to \code{quantile.default} with \code{type=1} 
#!  in that it only returns existing values, it is mostly symetric 
#!  but it is using 'round' rather than 'floor'. 
#!  \cr
#!  Note that this implies that \code{median.integer64} does not interpolate for even number of values 
#! (interpolation would create values that could not be represented as 64-bit integers).
#!  \cr
#!   This function automatically chooses from several low-level functions considering the size of \code{x} and the availability of a cache. 
#!   Suitable methods are \code{\link{sortqtl}} (fast sorting) 
#! and \code{\link{orderqtl}} (memory saving ordering).
#! }
#! \value{
#!   \code{prank} returns a numeric vector of the same length as \code{x}.
#!   \cr
#!   \code{qtile} returns a vector with elements from \code{x} 
#!   at the relative positions specified by \code{probs}.
#! }
#! \author{
#! 	Jens Oehlschlägel <Jens.Oehlschlaegel@truecluster.com>
#! }
#! \seealso{
#!   \code{\link{rank.integer64}} for simple ranks and \code{\link{quantile}} for quantiles.
#! }
#! \examples{
#! x <- as.integer64(sample(c(rep(NA, 9), 1:9), 32, TRUE))
#! qtile(x, probs=seq(0, 1, 0.25))
#! quantile(x, probs=seq(0, 1, 0.25), na.rm=TRUE)
#! median(x, na.rm=TRUE)
#! summary(x)
#! 
#! x <- x[!is.na(x)]
#! stopifnot(identical(x,  unname(qtile(x, probs=prank(x)))))
#! }
#! \keyword{univar}

qtile <- function(x, probs = seq(0, 1, 0.25), ...)UseMethod("qtile")
qtile.integer64 <- function(x, probs = seq(0, 1, 0.25), names = TRUE, method = NULL, ...){
	if (any(is.na(probs) | probs<0 | probs>1))
		stop("p outside [0,1]")
  c <- cache(x)
  if (is.null(method)){
    if (is.null(c)){
		method <- "sortqtl"
	}else{
		if (exists("sort", envir=c, inherits=FALSE))
			method <- "sortqtl"
		else if (exists("order", envir=c, inherits=FALSE))
			method <- "orderqtl"
		else
			method <- "sortqtl"
	}
  }
  switch(method
  , sortqtl={
		if (is.null(c) || !exists("sort", c, inherits=FALSE)){
			s <- clone(x)
			na.count <- ramsort(s, na.last=FALSE)
		}else{
			s <- get("sort", c, inherits=FALSE)
			na.count <- get("na.count", c, inherits=FALSE)
		}
		qs <- sortqtl(s, na.count, probs)
    }
  , orderqtl={
		if (is.null(c) || !exists("order", c, inherits=FALSE)){
			o <- seq_along(x)
			na.count <- ramorder(x, o, na.last=FALSE)
		}else{
			o <- get("order", c, inherits=FALSE)
			na.count <- get("na.count", c, inherits=FALSE)
		}
		qs <- orderqtl(x, o, na.count, probs)
    }
  , stop("unknown method", method)
  )
  if (names){
	np <- length(probs)
	dig <- max(2L, getOption("digits"))
	names(qs) <- paste(if (np < 100) 
		formatC(100 * probs, format = "fg", width = 1, digits = dig)
	else format(100 * probs, trim = TRUE, digits = dig), 
		"%", sep = "")
  }
  qs
}


quantile.integer64 <- function(x, probs = seq(0, 1, 0.25), na.rm = FALSE, names = TRUE, type=0L, ...){
	if (type[[1]]!=0L)
		stop("only type==0 ('qtile') supported")
	if (!na.rm && na.count(x)>0)
		stop("missing values not allowed with 'na.rm='==FALSE")
	qtile.integer64(x, probs = probs, na.rm = na.rm, names = names, ...)
}

median.integer64 <- function(x, na.rm=FALSE){
	if (!na.rm && na.count(x)>0)
		stop("missing values not allowed with 'na.rm='==FALSE")
	qtile.integer64(x, probs = 0.5, na.rm = na.rm, names = FALSE)
}

# mean.integer64 <- function(x, na.rm=FALSE){
	# s <- sum(x, na.rm=na.rm)
	# if (!is.na(s)){
		# if (na.rm)
			# s <- s%/%(length(x)-na.count(x))
		# else
			# s <- s%/%length(x)
	# }
	# s
# }
mean.integer64 <- function(x, na.rm=FALSE, ...){
	ret <- double(1)
	.Call("mean_integer64", x, as.logical(na.rm), ret)
	oldClass(ret) <- "integer64"
	ret
}

summary.integer64 <- function (object, ...){
	nas <- na.count(object)
	qq <- quantile(object, na.rm=TRUE)
	qq <- c(qq[1L:3L], mean(object, na.rm=TRUE), qq[4L:5L])
    names(qq) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.")
	if (any(nas)) 
		c(qq, "NA's" = nas)
	else qq
}
