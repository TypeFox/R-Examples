design <-
function(
formula,
data
)

{

# varying intercept at first position
    arg <- strsplit(gsub(" ", "", formula[length(formula)]), "\\+")
    int <- if (grepl("v\\(1",arg)==TRUE) 0 else 1
    if (int == 0) {
        pint <- grep("v\\(1",arg[[1]])
        arg2 <- paste(c(arg[[1]][pint], arg[[1]][(1:length(arg[[1]]))[-pint]]),"+",collapse="")
        formula[length(formula)]<-parse(text=substr(arg2,1,nchar(arg2)-1))
        }

# model.frame
    special <- c("v", "p", "grouped", "grouped.fused", "sp", "SCAD", "elastic", "vspline")
    m <- model.frame(formula=terms(formula, specials=special, data=data), data)
    if (nrow(m) == 0)
        stop("No (non-missing) observations")
        
# model.matrix
    Terms <- attr(m, "terms")
    attr(Terms, "intercept") <- 1
    X <- model.matrix(Terms, m)
    if (int==0) X <- X[,-1]

    namen <- colnames(X)
    namen <- sub(",", ".", namen, fixed=TRUE)
    namen <- gsub(" ", "", namen, fixed=TRUE)
    namen <- sub("vspline(", "", namen, fixed=TRUE)
    namen <- sub("v(", "", namen, fixed=TRUE)
    namen <- sub("sp(", "", namen, fixed=TRUE)
    namen <- sub("p(", "", namen, fixed=TRUE)
    namen <- sub("grouped(", "", namen, fixed=TRUE)
    namen <- sub("grouped.fused(", "", namen, fixed=TRUE)
    namen <- sub("SCAD(", "", namen, fixed=TRUE)
    namen <- sub("elastic(", "", namen, fixed=TRUE)
    namen <- gsub("(", "", namen, fixed=TRUE)
    namen <- gsub(")", "", namen, fixed=TRUE)
    namen <- gsub("n=", "", namen, fixed=TRUE)
    namen <- gsub("TRUE", "", namen, fixed=TRUE)
    namen <- gsub("FALSE", "", namen, fixed=TRUE)
    namen <- gsub(",", "", namen, fixed=TRUE)    
    namen <- gsub("bj=", "", namen, fixed=TRUE)
    namen <- gsub("knots=", "", namen, fixed=TRUE)
    namen <- gsub("\"", "", namen)#, fixed=TRUE)
    colnames(X) <- namen

# formula
    label <- attr(Terms, "term.labels")
    if (int == 1) {
        label <- c("1", label)
        }
    arg3 <- paste(label,"+",collapse="")
    formula[length(formula)]<-parse(text=substr(arg3,1,nchar(arg3)-1))

# return
return(list(X=X, Terms=Terms, m=m, int=int, formula=formula))

}
