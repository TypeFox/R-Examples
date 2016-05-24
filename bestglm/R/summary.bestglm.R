summary.bestglm <-
function (object, SubsetsQ=FALSE, ...){
ti<-object$Title
indNL<-(gregexpr(pattern="\n", text=ti))[[1]]
WhichIC <- ifelse(indNL<1, ti, substr(start=1, stop=indNL-1, ti))
MR <- object$ModelReport
FittedModel <- paste(WhichIC, "-", ifelse(MR$LEAPSQ,"leaps",ifelse(MR$gaussianQ,"lm","glm")),sep="")
cat(paste("Fitting algorithm: ", FittedModel), fill=TRUE)
NullModel<-MR$NullModel
BestModel <- object$BestModel
bestDF <- BestModel$df.residual
bestDeviance <- deviance(BestModel)
outM <- matrix(c(bestDF, MR$NullModel["DF"], bestDeviance, MR$NullModel["Deviance"]), ncol=2)
dimnames(outM)<-list(c("Null Model", "Full Model"),c("df", "deviance"))
#
#likelihood-ratio test
outHTest<-vector("list", 3)
X <- MR$NullModel["Deviance"] - bestDeviance
dfX <- MR$NullModel["DF"] - bestDF
outHTest$statistic <- X
names(outHTest$statistic) <- "X"
outHTest$parameter <- dfX
names(outHTest$parameter) <- "df"
outHTest$p.value <- 1 - pchisq(X, dfX)
outHTest$method <- "likelihood-ratio test - GLM"
outHTest$data.name <- paste("H0: Null Model vs. H1: Best Fit", FittedModel )
class(outHTest)<-"htest"
if (dfX==0)
    cat("Best model is null model!", fill=TRUE)
else {
    cat("Best Model:", fill=TRUE)
    print(outM)
    print(outHTest)
    }
#
if (SubsetsQ && !any(is.na(object$Subsets))) {
    cat("Best subsets:", fill=TRUE)
    print.data.frame(object$Subsets)
    }
}

