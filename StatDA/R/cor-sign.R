cor.sign <- function(x,method=c("pearson","kendall", "spearman"))
{
# compute correlation matrix of x with method "pearson", "kendall", "spearman"
# print resulting matrix with significance levels
#
R=cor(x,method=method)
# define matrix for p-values
P=matrix(0,nrow(R),ncol(R))
dimnames(P) <- dimnames(R)
# compute p-values:
for (i in 1:(ncol(P)-1)){
  for (j in (i+1):nrow(P)){
    P[i,j] <- cor.test(x[,i],x[,j],method=method)$p.value
    P[j,i] <- P[i,j]
  }
}
out <- list(cor=R,p.value=P)
class(out) <- "corsign"
out
}

print.corsign <-
function (x, symb=c("***","**","*","."," "),cutp=c(0.001,0.01,0.05,0.1),
          digits = 2, quote = FALSE, ...)
{
# x ... resulting object of function cor.sign (class corsign)
# symb ... symbols to be used
# cutp ... cutpoints for symbols
# digits ... number of digits to print
#
# define resulting object to be printed
if (length(symb)!=(length(cutp)+1)) stop("Wrong length of symb or cutp!")
res <- unclass(x$p.value)
res[res<cutp[1]] <- symb[1]
if (length(symb)>2){
  for (i in 2:(length(symb)-1)){
    res[(res>=cutp[i-1]) & (res<cutp[i])] <- symb[i]
  }
}
if (length(symb>1)){
  res[res>=cutp[length(cutp)]] <- symb[length(symb)]
}

res[lower.tri(res,diag=TRUE)] <- format(round(x$cor[lower.tri(x$cor,diag=TRUE)],digits))
print(res,quote=quote, ...)
cat("\n")
invisible(x)
}

