dmes <- function(x,y) {#Calculation of nonparametric orddom effect sizes
 dom <- dm(x,y)
 RET <- list()
 nam <- c("nx","ny","PSc","Ac","dc","NNTc","PSw","Aw","dw","NNTw","PSb","Ab","db","NNTb")
 nx  <- length(x)
 ny  <- length(y)
 nxny<- length(dom)
 PSc  <- (sum(dom<0)/nxny)
 Ac   <- PSc+(0.5*sum(dom==0)/nxny)
 dc   <- PSc-(sum(dom>0)/nxny)
 if(nx==ny){#assuming paired case
 PSw  <- (sum(diag(dom)<0)/nx)
 Aw   <- PSw+(0.5*sum(diag(dom)==0)/nx)
 dw   <- PSw-(sum(diag(dom)>0)/nx)
 PSb  <- (sum(dom<0)-sum(diag(dom)<0))/(nxny-nx)
 Ab   <- PSb+(0.5*(sum(dom==0)-sum(diag(dom)==0))/(nxny-nx))
 db   <- PSb-((sum(dom>0)-sum(diag(dom)>0))/(nxny-nx))
 RET[nam] <- c(nx,ny,PSc,Ac,dc,(1/dc),PSw,Aw,dw,(1/dw),PSb,Ab,db,(1/db))} else {
 RET[nam] <- c(nx,ny,PSc,Ac,dc,(1/dc),PSc,Ac,dc,(1/dc),PSc,Ac,dc,(1/dc))}
 return(RET)}
