"fun.diag.ks.g" <-
function(result, data, no.test = 1000, len = floor(0.9*length(data)), param)
{
data <- data[sample(length(data), no.test * len, TRUE)]
sample.fitted <- rgl(len * no.test, result[1], result[2], result[3], result[4], param)
test <- split(data, 1:no.test)
sample.fitted <- split(sample.fitted, 1:no.test)
result.o <- sum(sapply(1:no.test, function(i, test, sample.fitted)
ks.gof(test[[i]], sample.fitted[[i]])$p.value, test, sample.fitted) > 0.05)
return(result.o)
}

