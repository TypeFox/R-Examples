"cochran.test" <-
function(object,data,inlying=FALSE)
{
DNAME <- deparse(substitute(object))

if (is.vector(object)) {

by.factor <- as.factor(1:length(data))
vars <- object
names(vars) <- levels(by.factor)

k <- length(data)
df <- mean(data)

}
else
{
 if (missing(data)) 
        data <- environment(object)

bn<-as.character(attr(terms(object),"variables")[-1])
by.factor<-as.factor(data[[bn[2]]])
vars <- tapply(data[[bn[1]]],by.factor,var)
names(vars) <- levels(by.factor)

k <- nlevels(by.factor)
df <- length(data[[bn[1]]])/k
}

if (inlying) {
value <- min(vars)/sum(vars)
group <- levels(by.factor)[which(vars == min(vars))]
method <- "Cochran test for inlying variance"
alt <- paste("Group",group,"has inlying variance")
pval <- pcochran(value,df,k)
}
else {
value <- max(vars)/sum(vars)
group <- levels(by.factor)[which(vars == max(vars))]
method <- "Cochran test for outlying variance"
alt <- paste("Group",group,"has outlying variance")
pval <- 1-pcochran(value,df,k)
}

RVAL <- list(statistic = c(C = value), parameter = c(df = df, k = k),
alternative = alt, p.value = pval, method = method, 
estimate = vars, data.name = DNAME)

class(RVAL) <- "htest"
return(RVAL)


}

