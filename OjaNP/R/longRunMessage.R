`longRunTimeMessage` <- 
function(p,k,N,silent,string){
    if (!silent){
        message("PLEASE NOTE: You have requested to compute ",floor(p*N),"\nhyperplanes in R^",k,". ",
                    "This may take a while.\n",sep="")
        if (.Platform$OS.type == "windows")
            flush.console()
    }
}
