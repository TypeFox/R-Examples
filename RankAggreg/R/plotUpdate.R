`plotUpdate` <-
 function(objScores, hist, nb, method, ...){
	mat = matrix(c(1, 2, 1, 3), 2)
	layout(mat)
	if(all(objScores == min(objScores)))
		hist(c(objScores,objScores[1]+.0001), nb, col="blue", border="darkblue", xlab="Objective function scores",
			main=paste("Sample Distribution at Iteration", nrow(hist)), xlim=c(objScores[1]-.5, objScores[1]+.5))
	else		
		hist(objScores, nb, col="blue", border="darkblue", xlab="Objective function scores",
			main=paste("Sample Distribution at Iteration", nrow(hist)))
	plot(1:nrow(hist), hist[,1], col="red", type="l", xlab="Iteration", ylab="Scores",
		main=paste("Fittest Individual in", ifelse(method=="CE","Sample","Population")))
	legend("topright", legend=c("Minimum"), fill=c("red"), box.lty=0)
	plot(1:nrow(hist), hist[,2], col="green", type="l", xlab="Iteration", ylab="Scores",
		main=paste("General", ifelse(method=="CE","Sample","Population"), "Fitness"))
	legend("topright", legend=c("Median"), fill=c("green"), box.lty=0)
}
	