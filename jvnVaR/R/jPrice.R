jPrice <-
function(name){
dateList <- NULL 
data('dateList', package='jvnVaR', envir=environment())
dataSelected <- NULL 
data('dataSelected', package='jvnVaR', envir=environment())
t <- dataSelected[dateList,c(name)]
naCommon <- which(is.na(t))
object <- t[-naCommon]
return(object)
}
