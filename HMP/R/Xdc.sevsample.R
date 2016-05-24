Xdc.sevsample <-
function(group.data, epsilon=10^(-4), est="mom"){
if(missing(group.data))
stop("group.data missing.")

n.groups <- length(group.data)
taxaK <- ncol(group.data[[1]])

if(tolower(est) == "mle"){
Xdc <- Xdc.statistics(group.data, epsilon)
}else if(tolower(est) == "mom"){
Xdc <- Xdc.statistics.MoM(group.data)
}else{
stop(sprintf("Est '%s' not found. Est must be 'mle' or 'mom'.", as.character(est)))
}

p.value <- 1-pchisq(q=Xdc, df=(n.groups-1)*taxaK, ncp=0, lower.tail=TRUE)

xdc.sevsamp.test <- list(Xdc, p.value)
names(xdc.sevsamp.test) <- c("Xdc statistics", "p value")

return(xdc.sevsamp.test)
}
