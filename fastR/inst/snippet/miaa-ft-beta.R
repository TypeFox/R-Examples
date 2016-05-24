miaa <- miaa05
length(miaa$FTPct)
beta.mom(miaa$FTPct)
# remove players who took no shots
someshots <- miaa$FTPct[miaa$FTA >=1]
length(someshots)
beta.mom(someshots) -> bmom1; bmom1
plot1 <- qqmath(someshots, 
    dist=function(x)qbeta(x,bmom1["shape1"],bmom1["shape2"]))
# remove players with fewer than 10 shots
tenshots <- miaa$FTPct[miaa$FTA >=10]
length(tenshots)
beta.mom(tenshots) -> bmom2; bmom2
plot2 <- qqmath(tenshots, 
    dist=function(x)qbeta(x,bmom2["shape1"],bmom2["shape2"]))
