Shift <-
function(parameter,update){
		parameter <- scale(parameter,center=update,scale=FALSE)
		return(parameter)
	}
