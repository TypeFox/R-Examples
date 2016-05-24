# other function used in the package

## format cat output
format_cat<-function(x) format(as.character(sprintf("%.7g", x)),width =8,justify="left")


## check range of values
check.range <- function(q,q.set) {q[q<=max(q.set) & q>=min(q.set)]}
