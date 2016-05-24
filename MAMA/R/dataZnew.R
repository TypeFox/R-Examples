posterior.mean<- function(data, varname, nsamp, permute = 0)
{
N<-length(GEDM(data))
dataZ<-function(data, varname, nsamp)
{
N<-length(GEDM(data))
dum<-GEDM(data)
cl<-clinical(data)
for (i in 1:N) {
group <- grep(varname,names(clinical(data)[[i]]))
lev<-levels(as.factor(as.numeric(as.factor(clinical(data)[[i]][,group]))))
dum[[i]]<-dum[[i]][, c(
sample(which(as.numeric(as.factor(clinical(data)[[i]][,group]))==lev[1]),nsamp),
sample(which(as.numeric(as.factor(clinical(data)[[i]][,group]))==lev[2]),nsamp))]
vars = list("var1", "var2")
tmp = data.frame(cl=c(rep(0,nsamp),rep(1,nsamp)), names=colnames(dum[[i]]))
names(tmp)<-c(varname,"names")
names(vars) = names(tmp)
cl[[i]] = tmp
}
newdata<-new("MetaArray", GEDM=dum, clinical=cl, datanames=datanames(data))

return(newdata)
}

sampl.data<-dataZ(data, varname, nsamp)
sampl.data<-as.list(sampl.data)
arg<-list()
for (i in c(1:N))
	{
	cl<-sampl.data[[i]]$clinical
	row.names(cl)<-cl$names
	ex<-as.matrix(sampl.data[[i]]$GEDM)
	arg[[i]]<- new("ExpressionSet", exprs=ex, phenoData=new("AnnotatedDataFrame", data=cl))
	 }
	 
merged<-mergeExprs2(arg, datanames(data))
z.stat<-Zscore(merged, pheno=rep(1,N), permute = permute)
class(z.stat)<-c(class(z.stat),"posterior.mean")
return(z.stat)

}





mergeExprs2<-function (arg, names)
{
check<-function (x) 
{
    if (!is.element(class(x), c("list", "mergeExpressionSet", 
        "ExpressionSet", "matrix"))) 
        stop("all data must be either a list, a mergeExpressionSet, a matrix, or an ExpressionSet")
}

    #arg <- list(...)
    x <- alist(... = )
    studynames <- alist(... = )
    k <- 0
    for (i in 1:length(arg)) {
        check(arg[[i]])
        if (is.element(class(arg[[i]]), "mergeExpressionSet")) {
            mm <- mergeget(arg[[i]])
            studynames[[i]] <- names(arg[[i]])
            for (j in 1:length(arg[[i]])) {
                k <- k + 1
                x[[k]] <- mm[[j]]
            }
        }
        else {
            k <- k + 1
            x[[k]] <- mergeget(arg[[i]])
            studynames[[i]] <- as.character(as.list(names)[[i]])

        }
    }
    studynames <- unlist(studynames)
    tt <- length(x)
    nnote <- matrix(NA, tt, 2)
    for (i in 1:tt) {
        if (i == 1)
            iid <- as.matrix(featureNames(x[[i]]))
        else iid <- rbind(iid, as.matrix(featureNames(x[[i]])))
        nnote[i, 2] <- ""
    }
    iid <- as.vector(sort(unique(iid)))
    for (i in 1:tt) {
        y <- assayData(x[[i]])[["exprs"]]
        idy <- featureNames(x[[i]])
        y.avg <- AverageDuplicates(y, idy)
        assayData(x[[i]]) <- list(exprs = as.matrix(y.avg$data))
        featureNames(x[[i]]) <- y.avg$acc
    }
    idmatrix <- matrix(0, length(iid), tt)
    index <- as.vector(nnote[, 2])
    for (i in 1:tt) {
        idx <- featureNames(x[[i]])
        cc <- match(iid, idx)
        idmatrix[, i] <- ifelse(is.na(cc), 0, 1)
    }
    colnames(idmatrix) <- studynames
    rownames(idmatrix) <- iid
    names(x) <- studynames
    merged <- new("mergeExpressionSet", data = x, geneStudy = idmatrix,
        notes = "")
    return(merged)
}








