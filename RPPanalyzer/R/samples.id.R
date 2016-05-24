`samples.id` <-
function(x,id.cols=c("sample","repli")){
	mmnt.lines <- which ( x[[3]][,"sample_type"]=="measurement")
	temp <- c(NULL)
		for (i in id.cols){
			temp <- paste(temp,as.character(x[[3]][mmnt.lines,i]),sep="")
			}
		identifier <- unique(temp)
				replis <- (length(temp)/length(identifier))
				count <- length(identifier)
			print(paste("found",count,"individual samples, spottet in",replis
					,"replicates"))
	 return(identifier)
	}

