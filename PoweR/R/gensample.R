gensample <- function(law.index, n, law.pars = NULL, check = TRUE, center = FALSE, scale = FALSE) {

    if (class(law.index) == "function") {

        return(.Call("gensampleRcpp",law.index,n,if (is.null(law.pars)) 0 else law.pars,nbparlaw=length(law.pars),as.character(match.call()[2]),as.integer(center), as.integer(scale),PACKAGE="PoweR"))

    } else {
    
        if(getRversion() < "3.1.0") dontCheck <- identity

        Claw.name <- paste("law",law.index,sep="")
  
        if (check) {  # The following instruction takes time! This is why the check argument exists.
            tmp <- names(getDLLRegisteredRoutines("PoweR")[[".C"]])
            ind.laws <- grep("law",tmp[grep("law",tmp)])
            if (!(law.index %in% ind.laws)) stop(paste("This law (",law.index,") has not been included in the package!",sep=""))
            
            name <- .C(dontCheck(Claw.name),0L,0.0,name=rep(" ",50),1L,rep(0.0,4),0L,1L,PACKAGE="PoweR")$name
            law.name <- gsub('\\','',gsub('$','',sub(' +$', '', paste(name,collapse="")),fixed=TRUE),fixed=TRUE)
            
            if (length(law.pars) > 4) stop("The maximum number of law parameters is 4. Contact the package author to increase this value.")      
        } else {
            law.name <- ""
        }
        
        if (is.null(law.pars)) {law.pars <- rep(0.0,4);nbparlaw <- 0L} else {nbparlaw <- length(law.pars);law.pars <- c(law.pars,rep(0.0,4-nbparlaw))}
        
        out <- .C(dontCheck(Claw.name),as.integer(n),x=rep(0.0,n),rep(" ",50),0L,
                  law.pars=as.double(law.pars),nbparlaw=as.integer(nbparlaw),1L,PACKAGE="PoweR")

        if (center) out$x <- out$x - mean(out$x)
        if (scale) out$x <- out$x / sd(out$x)
        
        return(list(sample=out$x,law=law.name,law.pars=out$law.pars[1:out$nbparlaw]))
    }
}
