var.asdummy <-
function (x) {
x_length <- length(x)   
x_fac <- as.factor(x)   
lev <- levels(x_fac)   
x_count <- nlevels(x_fac)   
dummyvarnames <- list()   
dummyname <- 0   

for (i in 1:x_count) {   
dummyname[i] <- paste (lev[i],"_DUMMY", sep="")
dummyvarnames <- rbind(dummyvarnames, dummyname[i])
}

dummydataset <- data.frame(matrix(0, nrow=x_length, ncol=x_count))
names(dummydataset) <- dummyvarnames

for (i in 1:x_count) {   
positionstrue <- which(x_fac == lev[i])
dummydataset[positionstrue,i] <- 1
}
return(dummydataset)
}
