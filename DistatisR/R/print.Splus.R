print.Splus <-
function (x,...) {


	res.Splus <- x
	if (!inherits(res.Splus, "Splus")) stop ("no convenient data")
	cat("**Results for S+ Matrix**\n")
	cat("*The results are available in the following objects:\n\n")	
	
	if(res.Splus$compact){
		#list(Splus=Splus)
		res <- array("", c(1, 2), list(1:1, c("name", "description")))
		res[1,] <- c("$Splus","S+ Matrix")
	}else{
		#list(SCP = CP3, F = F, PartialF = PartialF, ProjectionMatrix = Proj, Splus=Splus)
		res <- array("", c(5, 2), list(1:5, c("name", "description")))
		res[1,] <- c("$SCP","Sum of Cross Products matrix")
		res[2,] <- c("$F","Factor scores")
		res[3,] <- c("$PartialF","Partial Factor Scores")
		res[4,] <- c("$ProjectionMatrix", "The Projection Matrix")
		res[5,] <- c("$Splus","S+ Matrix")
	}
	print(res)
}
