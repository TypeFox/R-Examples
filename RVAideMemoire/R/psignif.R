.psignif <-
function (p) {
  result <- character(length(p))
  for (i in 1:length(p)) {
    if (p[i]!="NA") {
	if (as.numeric(p[i]) >= 0.1) {result[i] <- " "} else
	if (as.numeric(p[i]) < 0.1 & as.numeric(p[i]) >= 0.05) {result[i] <- "."} else
	if (as.numeric(p[i]) < 0.05 & as.numeric(p[i]) >= 0.01) {result[i] <- "*"} else
	if (as.numeric(p[i]) < 0.01 & as.numeric(p[i]) >= 0.001) {result[i] <- "**"} else
	if (as.numeric(p[i]) < 0.001) {result[i] <- "***"}
    } else {
	result[i] <- " "
    }
  }
  return(result)
}

