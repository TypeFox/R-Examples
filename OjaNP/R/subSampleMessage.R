`subSampleMessage` <- 
function(p,silent,string){
    if (!silent){  
        printP <- floor(100*p+0.5)
        if (printP==0) printString <- "less than 1"
        else
            if (printP==100) printString <- "more than 99"
            else printString <- paste("about ",printP,sep="")
        message("PLEASE NOTE: The computation of the Oja ",string," is based on a \nrandom sub-sample of ",
              printString,"% of all possible hyperplanes.\n",sep="")
        if (.Platform$OS.type == "windows")
            flush.console()
    }
}
