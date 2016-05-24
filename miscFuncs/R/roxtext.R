##' roxtext function
##'
##' A function to generate roxygen documentation templates for functions for example,
##'
##'
##' would generate a template for this function. Note that functions with default arguments that include quotes
##' will throw up an error at the moment, just delete these bits from the string, and if shold work.
##'
##' @param s a string enclosed in quotes
##' @return minimal roxygen template
##' @export

roxtext <- function(s){ # s a string
    s <- unlist(strsplit(s,"<-"))
    s <- gsub(" ","",s)
    fname <- s[1]
    s <- s[2]
    left <- gregexpr("\\(",s)[[1]][1]
    right <- rev(gregexpr(")",s)[[1]])[1]
    s <- substr(s,left+1,right-1)
    arglist <- unlist(strsplit(s,","))
    cat("##' ",fname," function\n",sep="")
    cat("##'\n")
    cat("##' A function to \n")
    cat("##'\n")
    for (a in arglist){
        if(length(grep("=",a))==0){
            cat("##' @param ",a," X \n",sep="")
        }
        else{
            x <- unlist(strsplit(a,"="))[1]
            x <- gsub(" ","",x)
            cat("##' @param ",x," X \n",sep="")
        }
    }
    cat("##' @return ...\n")
    cat("##' @export\n")
}
