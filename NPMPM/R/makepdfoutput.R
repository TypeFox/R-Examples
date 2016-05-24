makepdfoutput <-
function(text){

textt <<- text

pdfoutput <- system.file("Sweave","pdfoutput.Snw",package="utils")

options(device.ask.default=FALSE)

today <- Sys.Date()


texname <- paste("pdfoutput",today,".tex",sep="",collapse="")

Sweave(file=pdfoutput,output=texname)

tools::texi2dvi(texname,pdf=TRUE,quiet=FALSE)


} # END FUNCTION

