updata <-
function(data = data, lastdose, npt, ndlt){#browser()
if (!(lastdose %in% data$dose)) {stop("dose do not treat")}
if (ndlt > npt) {stop("ndlt > npt is impossible")}
if (npt == 0){return(data)}
if (ndlt < 0 | npt < 0) {stop("ndlt & npt must be positive")}
idx<-which(data$dose==lastdose)
data[idx, "npt"] <- data[idx,"npt"] + npt
data[idx, "ndlt"] <- data[idx,"ndlt"] + ndlt
return(data)
}
