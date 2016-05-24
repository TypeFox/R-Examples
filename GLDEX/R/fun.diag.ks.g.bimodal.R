"fun.diag.ks.g.bimodal" <-
function(result1, result2,prop1,prop2,data, no.test = 1000, len = floor(0.9 * length(data)), param1,param2)
{   

if(length(result1)==9){
prop2<-1-result1[9]
prop1<-result1[9]
result2<-result1[5:8]
result1<-result1[1:4]
}

if(missing(prop2)){
prop2<-1-prop1
}

data <- data[sample(length(data), no.test * len, TRUE)]
no.1<-round(len * no.test*prop1)
sample.fitted1 <- rgl(no.1, result1[1], result1[2], result1[3], result1[4], param1)
sample.fitted2 <- rgl(len * no.test-no.1, result2[1], result2[2], result2[3], result2[4], param2)

test <- split(data, 1:no.test)
sample.fitted <- split(c(sample.fitted1,sample.fitted2)[sample(length(data), no.test * len, TRUE)], 1:no.test)
result.o <- sum(sapply(1:no.test, function(i, test, sample.fitted)
ks.gof(test[[i]], sample.fitted[[i]])$p.value, test, sample.fitted) > 0.05)
return(result.o)
}

