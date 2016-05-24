xval.buffer <-
function(n.cases, n.xval=5, buffer.length=0)
{
    case.indices <- 1:n.cases
    xval <- suppressWarnings(matrix(c(case.indices,
                rep(NA, length(case.indices))), ncol=n.xval*2))
    xval <- xval[,1:(min(which(is.na(xval[1,])))-1),drop=FALSE]
    cases <- list()
    for(i in 1:ncol(xval)) {
        x.cur <- xval[!is.na(xval[,i]),i]
        exclude <- intersect(case.indices, (min(x.cur)-buffer.length):
                             (max(x.cur)+buffer.length))
        cases[[i]] <- list(train=setdiff(case.indices, exclude), valid=x.cur)
    }
    cases
}
