normscores <-
function (x) 
(x - mean(x))/sqrt(sum((x - mean(x))^2))

