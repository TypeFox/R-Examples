"fun.which.zero" <-
function(data)
{
data.m<-cbind(1:length(data), data)
result <- data.m[data.m[, 2] == 0, 1]
return(result)
}

