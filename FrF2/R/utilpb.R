## Plackett Burman generating vectors and matrices, according to 1946 article
## are in pb.list (in workspace sysdata.rda)
## most functions in here are used for extending the pb designs with entries 0 in this list
## (double.des, Williamson etc.)

## function for doubling designs
## normalize.row.last is applied later to ensure that last row is "-1"
double.des <- function(des){
   zeil <- nrow(des)
   ### +/- column and reflected columns unconventionally moved to the front !!!
   ### benefit: first columns are resolution IV
   rbind(cbind(rep(1,zeil),des,des),cbind(rep(-1,zeil),-des,des))
}

circ.mat <- function(gen){
   n <- length(gen)+1
   sel <- gen
   for (i in 1:(n-2))
    sel <- c(sel,gen[(n-i):(n-1)],gen[1:(n-i-1)])
    matrix(sel,n-1,n-1,byrow=1)
}

williamson <- function(A,B,C,D){rbind(
  cbind(A, B, C, D),
  cbind(-B, A, D, -C),
  cbind(-C, -D, A, B),
  cbind(-D, C, -B, A))
}

normalize.col.first <- function(mat){
    hilf <- mat[,1]==-1
    mat[hilf,] <- -mat[hilf,]
    mat[,-1]
}

normalize.row.last <- function(mat){
    hilf <- mat[nrow(mat),]==1
    mat[,hilf] <- -mat[,hilf]
    mat
}

des.recode <- function (var, recodes, as.factor.result, char) 
{
    recode.list <- rev(strsplit(recodes, ";")[[1]])
    is.fac <- is.factor(var)
    if (missing(as.factor.result)) 
        as.factor.result <- is.fac
    if (missing(char)) char <- FALSE
    result <- var
    if (is.fac) 
        result <- as.character(result)

    for (term in recode.list) {
            set <- eval(parse(text = strsplit(term, "=")[[1]][1]))
        if (!char)
            target <- eval(parse(text = strsplit(term, "=")[[1]][2]), 
                envir=parent.frame(), enclos=sys.frame(0))
        else 
        target <- strsplit(term, "=")[[1]][2]
        for (val in set) {
                if (is.na(val)) 
                  result[is.na(var)] <- target
                else result[var == val] <- target
        }
    }
    if (as.factor.result) result <- as.factor(result)
    result
}
