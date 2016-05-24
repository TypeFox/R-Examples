.versionSuff <- function(name){
    paste(sep="", name, if(getRversion()<"2.16") ".O" else ".N")
}
