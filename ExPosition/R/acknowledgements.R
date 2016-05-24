acknowledgements <-
function(){
	print("This function serves to keep a list of those who have contributed to ExPosition (and related packages) throughout development.")
	
	num.people <- 11
	
	full.list <- array("", c(num.people, 1), list(1:num.people, c("Person")))
	#rbind(
		i<-1
		full.list[i,] <- c("Francesca Filbey: Data"); 	i<-i+1
		full.list[i,] <- c("Michael Meyners: Feedback and suggestions for prettyGraphs"); 	i<-i+1
		full.list[i,] <- c("Brad Buchsbaum: Feedback, testing, and suggestions"); 	i<-i+1		
		full.list[i,] <- c("Anjali Krishnan: Testing and feedback"); 	i<-i+1
		full.list[i,] <- c("Michael Kriegsman: Testing and feedback"); 	i<-i+1
		full.list[i,] <- c("Jenny Wong: Testing and feedback"); 	i<-i+1
		full.list[i,] <- c("Daniel Faso: Testing and feedback"); 	i<-i+1
		full.list[i,] <- c("Shaikat Hossain: Testing and feedback"); 	i<-i+1
		full.list[i,] <- c("Amy Louise Schwarz: Testing and feedback"); 	i<-i+1
		full.list[i,] <- c("Adam Teed: Testing and feedback"); 	i<-i+1
		full.list[i,] <- c("Rachel Williams: Data entry and creation"); 	i<-i+1
 		full.list[i,] <- c("Students of Research Methods 3 at UT Dallas (2010, 2011, and 2012): Feedback, suggestions, interface testing and quality control, suffering"); 	i<-i+1		
	#)
	print(full.list)
}
