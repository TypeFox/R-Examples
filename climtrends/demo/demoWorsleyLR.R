# Worsley Likelihood Ratio 
# example based on:
# TREND Tests in Excel Spreadsheet
# eWater toolkit TREND Version 1.0.2
# http://www.toolkit.net.au/Tools/TREND/documentation

test1 <- matrix(c(1940,681,
1941,3661,
1942,8625,
1943,2475,
1944,573,
1945,2794,
1946,10190,
1947,5143,
1948,4139,
1949,8945,
1950,7295,
1951,19883,
1952,12119,
1953,8772,
1954,8848,
1955,16309,
1956,16254,
1957,2303,
1958,7671,
1959,3985,
1960,13742,
1961,5333,
1962,4859,
1963,12381,
1964,12137,
1965,6075,
1966,4669,
1967,378,
1968,7507,
1969,3891,
1970,13046,
1971,12954,
1972,2445,
1973,14759,
1974,20200,
1975,16331,
1976,6922,
1977,6739,
1978,11629,
1979,7351,
1980,2445,
1981,9960,
1982,10,
1983,11786,
1984,10214,
1985,11216,
1986,8393,
1987,10005,
1988,6896,
1989,11632),ncol=2,byrow=TRUE)
colnames(test1) <- c('year','data')

WLR <- WorsleyLikelihoodRatio(test1,returnZk=TRUE)
cat('Worsley Likelihood Ratio\n\nYear of change =',WLR$YearChange,'\nV =',WLR$V,'\nW =',WLR$W)
plot(test1[,1], WLR$Zk ,type='l',xlab='',ylab='Worsley Likelihood Ratio')


