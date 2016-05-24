"MCbound" <-
function(type,parms,conf.level=.99){
    if (type=="fixed"){
        if (is.null(names(parms)) & length(parms)==1){
            out<-MCbound.fixed(parms,conf.level) 
        }
        else if (!is.null(names(parms)) & !is.na(parms["Nmax"])){
            out<-MCbound.fixed(parms["Nmax"],conf.level)
        }
        else{ stop("Fixed parms vector must either have length one or have a named value of Nmax") }
    }
    else if (type=="tsprt"){
        if (is.null(names(parms)) & length(parms)==5){
            warning("parms vector not named: assume parms<-c(p0,p1,A,B,Nmax)")
            out<-MCbound.tsprt(parms,conf.level) 
        }
        else if (!is.null(names(parms))){
            out<-MCbound.tsprt(parms,conf.level)
        }
        else{ stop("tsprt parms vector must either have length 5 or named values (see help)") }
    }
    else if (type=="Bvalue"){
         if (is.null(names(parms)) & length(parms)==4){
           warning("parms vector not named: assume parms<-c(Nmax,alpha,e0,e1)")
           out<-MCbound.Bvalue(parms[1],parms[2],parms[3],parms[4],conf.level) 
        }
        else if (!is.null(names(parms))){
            out<-MCbound.Bvalue(parms["Nmax"],parms["alpha"],parms["e0"],parms["e1"],conf.level)
        }
        else{ stop("Bvalue parms must either have length 4 or named values (see help)")}
     }
    else if (type=="BC"){
        if (is.null(names(parms)) & length(parms)==2){
            warning("parms vector not named: assume parms<-c(Nmax,Smax)")
            out<-MCbound.BC(parms[1],parms[2],conf.level) 
        }
        else if (!is.null(names(parms))){
            out<-MCbound.BC(parms["Nmax"],parms["Smax"],conf.level)
        }
        else{ stop("BC parms must either have length 2 or named values (see help)") }
      }
    else { stop("type must equal either: fixed, tsprt, Bvalue, or BC")  }
    out
}

