#'## Function to order the nodes in a Network in enaR
#'## Singh P. | July 2014
#'## -----------------------------------------

netOrder <- function(x,order=0) {
    if (class(x) != "network") {
        stop("x is not a network class object")
    }
                                        # Load Initials
    flow <- as.matrix(x, attrname = "flow")
   # input <- x %v% "input"
   # resp <- x %v% "respiration"
   # export <- x %v% "export"
   # output <- x %v% "output"
   # storage <- x %v% "storage"
    living <- x %v% "living"
   # names <- x %v% "vertex.names"
    N <- length(living)
                                        # Determine Order (ordr)

    if(identical(order,0)==TRUE) {
        order<-rep(0,N)
        liv1<-which(living==TRUE)
        liv2<-which(living==FALSE)
        order<-c(liv1,liv2)
        if(identical(order,1:N)==TRUE) {warning('Network meets default conditions, no changes made')}
    }



                                        # Rearrange Network Characteristics
   # living <- living[ordr]
    flow  <- flow[order,order]
   # export <- export[ordr]
   # resp <- resp[ordr]
   # storage <- storage[ordr]
   # output <- output[ordr]
   # input <- input[ordr]
   # names <- names[ordr]


                                        # Modify Network
    x<-permute.vertexIDs(x,order)
    set.edge.attribute(x, 'flow', flow[flow>0])
    #x %v% "input" <- input
    #x %v% "respiration" <- resp
    #x %v% "export" <- export
    #x %v% "output" <- output
    #x %v% "storage" <- storage
    #x %v% "living" <- living
    #x %v% "vertex.names" <- names


                                        # Return the ordered network
    return(x)

}
