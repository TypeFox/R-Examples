#' Easy simulation of PK models
#'        
#' Easy simulation of basic PK models
#'        
#' See http://simulx.webpopix.org/mlxr/pkmodel/ for more details.
#' @export
#' @param time a vector
#' @param treatment a list with fields
#' \itemize{
#'   \item \code{time} : a vector of input times,
#'   \item \code{amount} : a scalar or a vector of amounts,
#'   \item \code{rate} : a scalar or a vector of infusion rates (default=\code{Inf}),
#'   \item \code{tinf} : a scalar or a vector of infusion times (default=0),
#' }
#' @param parameter vector of parameters with their names and values
#' 
#' @examples    
#' \dontrun{
#'   adm <- list(time=c(2,14,20), amount=40)
#'   p   <- c(V=8, Cl=0.5,k12=0.3, k21=0.2)
#'   t   <- seq(0, 30, by=0.1)
#'   
#'   res   <- pkmodel(time = t, treatment = adm, parameter = p)
#'   
#'   print(ggplot(data=res, aes(x=time, y=cc)) + geom_line(size=1) +
#'     xlab("time (h)") + ylab("concentration (mg/L)"))
#'   
#'   adm <- list(time = c(1,23,37,45), amount = c(1,0.5,2,0.3))
#'   p <- c(Mtt=5, Ktr=1, ka=0.5, V=10, Vm=1, Km=0.6, p=0.5)
#'   t <- seq(0, 80, by=0.1)
#'   
#'   res <- pkmodel(t,adm,p)
#'   
#'   print(ggplot(data=res, aes(x=time, y=cc)) + geom_line(size=1) +
#'     xlab("time (h)") + ylab("concentration (mg/L)"))
#'   
#'   adm <- list( time = 2, amount = 40)
#'   
#'   p <- inlineDataFrame("
#'   id   ka   V    Cl
#'   1   0.5   4     1
#'   2     1   6     1
#'   3   1.5   6   1.5
#'   ")
#'   
#'   t <- seq(0, 30, by=0.1)
#'   
#'   res <- pkmodel(t,adm,p)
#'   
#'   print(ggplot(data=res, aes(x=time, y=cc, colour=id)) + geom_line(size=1) +
#'     xlab("time (h)") + ylab("concentration (mg/L)"))   
#'   adm <- list(time=seq(2, 100, by=24), amount=40, rate=5)
#'   p <- c(V=8, Cl=0.5, k12=0.3, k21=0.2, ke0=0.2)
#'   t <- seq(0, 50, by=0.1)
#'   
#'   res <- pkmodel(t,adm,p)
#'   
#'   if( require("reshape2") ){
#'     r <- melt(res, id='time', variable.name='c')
#'     print(ggplot(r, aes(time,value)) + geom_line(aes(colour = c),size=1) +
#'             ylab('concentration') + guides(colour=guide_legend(title=NULL)) +
#'             theme(legend.position=c(.9, .8)))
#'   } 
#' }
pkmodel <- function(time,treatment,parameter){
  # ########################################################################################  
  #  pkmodel.R is governed by the CeCILL-B license. 
  #  You can  use, modify and/ or redistribute the software under the terms of 
  #  the CeCILL license as circulated by CEA, CNRS and INRIA at the following URL
  #  http://www.cecill.info/index.en.html
  #
  #  pkmodel.R was developed by Marc Lavielle and Fazia Bellal (Inria) for the DDMoRe project. 
  # ########################################################################################  
  if (!is.list(parameter)){
    parameter <- list(name=names(parameter), value=as.numeric(parameter))
  }
  if (is.data.frame(parameter)){
    np <- names(parameter)
    jid <- which(np=="id")
    parameter <- list(name=np[-jid], colNames=np, value=as.matrix(parameter))
  }
  pn=parameter$name
  
  #   pn=names(parameter)
  if(length(grep("ke0",pn))>0){iop.ke0=1}else{iop.ke0=0}
  if(iop.ke0==0){
    out <- list(name="cc",time=time)
  }
  else{
    out <- list(name=c("cc","ce"),time=time)
  }
  
  
  if ("V2" %in% parameter$name){
    if ("V1" %in% parameter$name){
      iv1 <- which(parameter$name=="V1")
      parameter$name[iv1]="V"
    }
    if ("Q" %in% parameter$name){
      iq <- which(parameter$name=="Q")
      parameter$name[iq]="Q2"
    }
    iv1 <- which(parameter$name=="V")
    iv2 <- which(parameter$name=="V2")
    iq2 <- which(parameter$name=="Q2")
    k12 <- parameter$value[iq2]/parameter$value[iv1]
    k21 <- parameter$value[iq2]/parameter$value[iv2]
    parameter$value[iq2] <- k12
    parameter$value[iv2] <- k21
    parameter$name[iq2]  <- "k12"
    parameter$name[iv2]  <- "k21"
  }
  
  if ("V3" %in% parameter$name){
    iv1 <- which(parameter$name=="V")
    iv3 <- which(parameter$name=="V3")
    iq3 <- which(parameter$name=="Q3")
    k13 <- parameter$value[iq3]/parameter$value[iv1]
    k31 <- parameter$value[iq3]/parameter$value[iv3]
    parameter$value[iq3] <- k13
    parameter$value[iv3] <- k31
    parameter$name[iq3]  <- "k13"
    parameter$name[iv3]  <- "k31"
  }
  
  data <- simulx(model="pkmodel",parameter=parameter,output=out,treatment=treatment)
  if(iop.ke0==0){
    r=data$cc
  }
  else{
    r=merge(data$cc,data$ce)
  }
  
  if (!is.null(r$id)){
    if (length(levels(r[,"id"]))==1)
      r[,"id"]=NULL
  }
  
  return(r)
}
