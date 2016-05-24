remainderReduce <- function(x,keepTruthTable = TRUE)
{ 
    if (!"truthTable"  %in% class(x) ) stop("x is not a truthTable.")
    call <- match.call()
    mydata <- x$truthTable
    conditions <- x$conditions
    nlevels <- x$nlevels
    if (keepTruthTable) {
         x.s <- sort(x, decreasing=FALSE,criter="OUT")
         truthTable <- x.s$truthTable
    }
    else {
        truthTable <- NULL
    }
    dat1 <- mydata[mydata[["OUT"]]=="?",conditions]
    dat0 <- mydata[mydata[["OUT"]]!="?",conditions]
    explained <- dat1
    idExclude <- apply(dat0,1,implicant2Id,nlevels=nlevels)
    if (nrow(explained)==0) stop("No remainder is in the truthTable. Have you passed the full truthTable to me?")
    primesId <- apply(dat1,1,implicant2Id,nlevels=nlevels)
    primesId <- reduce2(primesId,nlevels=nlevels)
    primeImplicants <- id2Implicant(primesId ,nlevels=nlevels,names=conditions)
    PIChart <- PIChart(primeImplicants,explained)
    sl <- solvePIChart(PIChart)
    solutions <- apply(sl,2,function(idx)primeImplicants[idx,])
    commonSolutions <- apply(sl,1,function(idx) {if (length(id <- unique(idx))==1) id })
    ans <- list(solutions=solutions,commonSolutions=commonSolutions,
                solutionsIDX=sl,primeImplicants=primeImplicants,
                truthTable=truthTable,explained=explained,outcome=x$outcome,
                idExclude=idExclude,nlevels=nlevels,
                PIChart=PIChart, call=call)
    class(ans) <- c("remainders","QCA")
    ans
}

