library(FrF2)
### test programs for fold.design

##  function for independent checks
mult <- function(a1, a2, ...){ 
    l <- list(...)
    if (length(l))
        return(Recall(a1, Recall(a2, ...)))
    if (any(is.na(as.numeric(as.character(a1))))) a1 <- 2*as.numeric(as.factor(a1))-3
    if (any(is.na(as.numeric(as.character(a2))))) a2 <- 2*as.numeric(as.factor(a2))-3
    A <- as.numeric(as.character(a1))*as.numeric(as.character(a2))
    A
}

## check for replications and reptowide designs
plan1 <- FrF2(8,5,randomize=FALSE)
plan2 <- FrF2(8,5,replications=2, randomize=FALSE)
set.seed(98776)
plan3 <- FrF2(8,5,replications=2, repeat.only=TRUE)
fp1 <- fold.design(plan1)
fp1b <- fold.design(plan1, columns=5)
fp2 <- fold.design(plan2)
fp3 <- fold.design(reptowide(plan3))

plan4 <- FrF2(16,6,randomize=FALSE, factor.names=list(eins=c(1,2),zwei=c(43,87),drei=c("alt","neu"),vier="",fuenf="",sechs=""))
plan5 <- FrF2(16,randomize=FALSE, factor.names=list(eins=c(1,2),zwei=c(43,87),drei=c("alt","neu"),vier="",fuenf="",sechs=""))
fp4 <- fold.design(plan4)
fp5 <- fold.design(plan5)
plan6 <- FrF2(16,design="6-2.2",randomize=FALSE)
plan7 <- FrF2(16,generators=catlg[["6-2.2"]]$gen,randomize=FALSE,factor.names=list(eins=c(1,2),zwei=c(43,87),drei=c("alt","neu"),vier="",fuenf="",sechs=""))
fp6 <- fold.design(plan6)
fp7 <- fold.design(plan7)
fp6b <- fold.design(plan6, columns=c(4,6))
fp7b <- fold.design(plan7, columns=c(4,6))

## alias structure for three generators that differ only by sign
plan8 <- FrF2(16,generators=c(7,13,15),randomize=FALSE)
fp8 <- fold.design(plan8) ## seems to be correct, is at least compatible
fp8b <- fold.design(plan8,columns=1:7)
plan9 <- FrF2(16,generators=c(7,-13,15),randomize=FALSE)
fp9 <- fold.design(plan9, columns=5) 
fp9b <- fold.design(plan9, columns=c(5,6,7))
fp9c <- fold.design(plan9)
plan10 <- FrF2(16,generators=c(-7,-13,-15),randomize=FALSE) 
fp10 <- fold.design(plan10, columns=5)  ## wrong sign in design.info somewhere
fp10b <- fold.design(plan10, columns=c(5,6,7))  

## estimable muss noch her

set.seed(98776)
plan11 <- FrF2(estimable=formula("~one+two+three+four+two:three+two:four"), 
       factor.names=c("one","two","three","four"), res3=TRUE)
fp11 <- fold.design(plan11)
  ## clear=FALSE allows to allocate all effects on distinct columns in the 
  ##     8 run MA resolution IV design
set.seed(98776)
plan12 <- FrF2(estimable=formula("~one+two+three+four+two:three+two:four"), 
       factor.names=c("one","two","three","four"), clear=FALSE)
fp12 <- fold.design(plan12)

  ## 7 factors instead of 6, but no requirements for factor G
set.seed(98776)
plan13 <-   FrF2(16, nfactors=7, estimable = formula("~A+B+C+D+E+F+A:(B+C+D+E+F)"), 
       clear=FALSE)
fp13 <- fold.design(plan13)

plan14 <-  FrF2(32,14,WPs=8,nfac.WP=4,randomize=FALSE)
fp14 <- fold.design(plan14)
fp14b <- fold.design(plan14,columns=c(1,2))
fp14c <- fold.design(plan14,columns=5:14)

plan15 <- FrF2(32,14,WPs=8,nfac.WP=2,randomize=FALSE)
fp15 <- fold.design(plan15)
fp15b <- fold.design(plan15, columns=c(1,2))
fp15c <- fold.design(plan15, columns=c(4,2))

liste <- c("fp1","fp1b","fp2","fp3","fp4","fp5","fp6","fp6b","fp7","fp7b","fp8","fp8b",
           "fp9","fp9b","fp9c","fp10","fp10b","fp11","fp12","fp13","fp14","fp14b","fp14c","fp15","fp15b","fp15c")

for (cc in liste){
   print(paste("*********************",cc,"*******************************"))
   print(eval(parse(text=paste("design.info(",cc,")$creator"))))
   print(eval(parse(text=paste("design.info(",cc,")$generators"))))
   print(eval(parse(text=paste("design.info(",cc,")$aliased"))))
   print(eval(parse(text=paste("design.info(",cc,")$res.WP"))))
   print(eval(parse(text=paste("design.info(",cc,")$type"))))
   print(eval(parse(text=cc)))
   }