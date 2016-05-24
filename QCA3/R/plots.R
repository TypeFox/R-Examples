## This is part of QCA3 project
## by HUANG Ronggui (2009-2010)

# plot.QCA <- function(x,...){
#     explain <- x$call$explain
#     truthTable <- x$truthTable
#     if (pmatch(explain,"positive",0)==1) Case <- truthTable[rownames(truthTable) %in% rownames(x$explained),]
#     if (pmatch(explain,"negative",0)==1) Case <- truthTable[rownames(truthTable) %in% rownames(x$explained),]
#     conditions <- names(x$explained)
#     idExplained <- apply(Case[,names(x$explained)], 1, implicant2Id,nlevel=x$nlevels)
#     ids <- rep(idExplained,Case$N)
#     Coverage <- apply(x$solutions[[1]],1, function(x) {
#         common <- intersect(esubSet(x),ids)
#         ## subSet or subCombination?
#         ids[ids %in% common]
#     }
#                       )
#     names(Coverage)<- paste("IM",seq_len(length(Coverage)),sep=".")
#     CovList <- lapply(Coverage,as.character)
#     require(venneuler) ## use venneuler instead of Vennerable
#     m <- data.frame(element=unlist(CovList),
#                     set=rep(names(CovList),lapply(CovList,length)))
#     v <- venneuler::venneuler(m)
#     venneuler:::plot.VennDiagram(v)
#     ##require(Vennerable)
#     ##plot(Venn(Sets=CovList),doWeights=TRUE)
# }


fsplot <- function(formula,data,main="fuzzy set plot",xlab=NULL,ylab=NULL,...){
    if (is.null(ylab)) ylab <- deparse(formula[[2]])
    if (is.null(xlab)) xlab <- deparse(formula[[3]])
    formula <- as.formula(paste(deparse(formula),"-1",sep=""))
    ## remove the intercept
    conditions <- model.matrix(formula,data=data)
    conditions <- apply(conditions,1,min)
    out <- model.response(model.frame(formula,data=data))
    plot(conditions,out,xlim=c(0,1),ylim=c(0,1),main=main,xlab=xlab,ylab=ylab,...)
    abline(0,1)
    abline(-0.1,1,lty=3)
    abline(0.1,1,lty=3)
    cos <- consistency(conditions,out,"less")
    ## causal sufficency of the conditions when score is high
    cov <- coverage(conditions,out,"less")
    ## proportion of memership in outcome accounted by conditions
    lab <-sprintf(
                  "Set-theoretic consistency: %.3f  Set-theoretic coverage: %.3f",cos,cov)
    mtext(lab,line=0.3)
}

fsnot <- function(x) {
    1-x
}

fsand <- function(...) pmin(...,na.rm=FALSE)

fsor <- function(...) pmax(...,na.rm=FALSE)
