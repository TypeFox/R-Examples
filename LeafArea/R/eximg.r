eximg <- function(){
os<-.Platform$OS.type
if (os == "windows") {
	temp <- tempdir()
	imagedir <- system.file("images",package="LeafArea")
	imagedir <- gsub("/","\\\\",imagedir)

	bat <- paste("cd ",imagedir,
            "\n copy A1-01.jpeg ",temp,"\\A1-01.jpeg 
            copy A1-02.jpeg ",temp,"\\A1-02.jpeg
            copy A2.jpeg ",temp,"\\A2.jpeg
            copy A123-01.jpeg ",temp,"\\A123-01.jpeg
            copy A123-02.jpeg ",temp,"\\A123-02.jpeg
            copy A300-1.jpeg ",temp,"\\A300-2.jpeg
            copy A300-2.jpeg ",temp,"\\A300-1.jpeg",sep="")

	tempbat <- paste(tempfile('bat'),".bat",sep="")
	write(bat, file=tempbat)
	shell(tempbat,intern=T)
	unlink(tempbat)
	return(paste(tempdir(),"\\",sep="")	)
} else {
	imagedir <- file.path(system.file("images",package="LeafArea"),"")
	temp <- tempdir()
	system(paste("cp -r", imagedir, temp))
	unix.check <-Sys.info()["sysname"]
	if(unix.check=="Linux") return(paste(tempdir(),"/images/",sep="")) else return(paste(tempdir(),"/",sep="")	)
	}
}
