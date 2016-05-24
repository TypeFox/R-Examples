readtext.ij <-function(path){
	file.list<-list.files(path)
	file.list <- file.list[grep(".txt$",file.list)]
	if (length(file.list[-grep("macro|bat",file.list)]!=0)) file.list <- file.list[-grep("macro|bat",file.list)]

	size<-length(0)
	file.name<-length(0)

	data <- list()

	temp.slash <- substr(path,nchar(path),nchar(path))
		if(temp.slash!="/" & temp.slash!="\\"){
			path <- paste(path,"/",sep="")
		}		 

	for (i in 1:length(file.list))
	{
	temp <- read.delim(paste(path,file.list[i],sep=""))
	data[[i]] <- data.frame(Area=temp$Area)
	# data[[i]] <- read.delim(paste(path,file.list[i],sep=""))
	names(data)[i]<-paste(file.list[i])
	}
	return(data)
	}

