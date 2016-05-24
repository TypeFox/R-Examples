numSims <-10000;
versi <- iris[iris$Species == "versicolor",];
corStat <- function(x,y) { sum( x * y ) - length(x) * mean(x) * mean(y) };
testStat <- with(versi, corStat(Sepal.Length,Petal.Length)); testStat;
simResults <- with(versi, replicate(numSims, 
                corStat(Sepal.Length,sample(Petal.Length))));
# 1-sided p-value
sum( simResults >= testStat ) / numSims;
# a reasonable 2-sided p-value
sum( abs(simResults) >= abs(testStat) ) / numSims;
