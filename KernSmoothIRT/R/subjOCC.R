 subjOCC <- function(x, stype = c("ObsScore","ExpectedScore","MLScore","Theta","MLTheta")){
 
 if(missing(stype)) stype <- "ObsScore"
 
 Plist <- list(x$nitem)
 
 for (i in 1:x$nitem){
	
	Single <- x$OCC[which(x$OCC[,1] == i),]
	SingleSubj <- matrix(NA,nrow = nrow(Single), ncol =(x$nsubj))
	
	if(stype == "Theta"){
		Thet <- x$evalpoints
		Sub  <- x$subjtheta
	}
	else if(stype == "MLTheta"){
		Thet <- x$evalpoints
		Sub  <- x$subjthetaML
	}
	else if(stype == "MLScore"){
		Thet <- x$expectedscores
		Sub  <- x$subjscoreML
	}
	
	else if(stype == "ExpectedScore"){
		Thet <- x$expectedscores
		Sub  <- subjETS(x)
	}
	
	else if(stype == "ObsScore"){
		Thet <- x$expectedscores
		Sub  <- x$subjscore
	}
	
	for(row in 1:nrow(Single)){
		SingleSubj[row,] <- approx(x=Thet, y=Single[row,-c(1:3)], xout=Sub)$y
	}
	
	Plist[[i]] <- SingleSubj
	
 }
 return(Plist)
 
 }
 
 
 
 
