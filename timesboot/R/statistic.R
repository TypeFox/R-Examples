statistic <-
function(ts){
cm = acf(ts,plot=FALSE)
return (cm$acf)
}
