memberX <- function(data,minsize){
#membership of observations to bins based on total score X
#this function is needed to compute stable residuals forward search
X.tot <- apply(data,1,sum)
sorted.R <- sort(X.tot)
N <- dim(data)[1]
group <- max(which(sorted.R == sorted.R[minsize]))
repeat {
    if (N - max(group) < minsize)        break
    group <- c(group, max(which(sorted.R == sorted.R[minsize +  max(group)])))
}
group <- group[-length(group)]
group <- c(sorted.R[group], max(sorted.R))
L <- length(group)
lo <- c(min(sorted.R), group[-L] + 1)
member <- apply(1 - outer(X.tot, group, "<="), 1, sum) +    1
return(member)
}
