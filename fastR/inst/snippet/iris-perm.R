data(iris); numSims <-10000;
setosa <- iris[iris$Species == "setosa",];
corStat <- function(x,y) { sum( x * y ) - length(x) * mean(x) * mean(y) };
testStat <- with(setosa, corStat(Sepal.Length,Petal.Length)); testStat;
simResults <- with(setosa, replicate(numSims, 
                corStat(Sepal.Length,sample(Petal.Length))));
# 1-sided p-value
sum( simResults >= testStat ) / numSims;
# a reasonable 2-sided p-value
sum( abs(simResults) >= abs(testStat) ) / numSims;
