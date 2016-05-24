dataFiller <-
function(data,NAstring=NA)
{


central.value <- function(x) {
if (is.numeric(x)) median(x,na.rm=T)
else if (is.factor(x)) levels(x)[which.max(table(x))]
else {
f <- as.factor(x)
levels(f)[which.max(table(f))]
}
}


dist.mtx <- as.matrix(daisy(data,stand=T))

ShowMissing=NULL

ShowMissing=data[which(!complete.cases(data)),]

for(r in which(!complete.cases(data)))
data[r,which(is.na(data[r,]))] <-
apply(data.frame(data[c(as.integer(names(sort(dist.mtx[r,])[2:11]))),
which(is.na(data[r,]))]),2,central.value)

cat("the missing case(s) in the orignal dataset ","\n\n")

print(ShowMissing)

cat("\n\n")

return(data)
}
