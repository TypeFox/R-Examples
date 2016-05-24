`f.score` <-
function (data, sens, measured, parameter, criterion, limit) 
{
	weigth<-sens[parameter,]>limit
	#functions <- array(dim=(max(c.krit)))
	#functions[[c.krit.rmse]] <- function(){return(0)}
	functions <- c(
		#rmse 
		function (x) {
		return(rmse(measured, as.numeric(x), weigth))
		},
		function(x) {
			return(nashS(measured,as.numeric(x), weigth))
		},
		function(x) {
			return(nashS_HF(measured,as.numeric(x), weigth))
		}
	)

	func <- functions[[criterion]]

	return(apply(data, 1, func))

}

