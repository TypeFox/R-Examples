as.data.frame.nPV <-
function(x, ...){

NSETS<-x$NSETS
NPVn<-NULL; PPVn<-NULL
for(i in 1:NSETS){
NPVn <- cbind(NPVn, rep(x$nlist[[i]][["NPV"]][["n"]], length.out=x$nsteps))
PPVn <- cbind(PPVn, rep(x$nlist[[i]][["PPV"]][["n"]], length.out=x$nsteps))
}

SETNAM <- make.names(rownames(x$outDAT))
NAMN <- paste("n.NPV.", SETNAM, sep="")
NAMP <- paste("n.PPV.", SETNAM, sep="")

colnames(NPVn)<-NAMN
colnames(PPVn)<-NAMP

as.data.frame(cbind(propP=x$propP, NPVn, PPVn), ...)

}

