require(drfit)
data(antifoul)
rantifoul <- drfit(antifoul)
rantifoul
drplot(rantifoul,antifoul,overlay=TRUE,bw=FALSE)
