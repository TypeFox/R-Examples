ibdhap.summary.calls <-
    function( calls, data.type = c("h", "g","r"), position=NA){
        
        
        ##define parameters based on haplotype/genotype/or reduced data
        
        if(length(data.type)>1){ stop("data.type improperly indicated, please see ?ibdhap.summary")}
        else if( is.element("h", data.type)){ no.ibd.ind = 15}
        else if( is.element("g", data.type)){no.ibd.ind = 9 }
        else if( is.element("r", data.type)){no.ibd.ind = 4}
        else{ stop("data.type improperly indicated, please see ?ibdhap.summary")}

        ##calculate mean proportions 
        mean.any.ibd<-sum((calls<no.ibd.ind)*(calls>0))/((nrow(calls))*ncol(calls))
        mean.no.call<-sum((calls==0))/((nrow(calls))*ncol(calls))
        mean.not.ibd<-( 1- mean.any.ibd - mean.no.call)
        
                                        #mean proportion
        mean.prop = c(any.ibd = mean.any.ibd, not.ibd = mean.not.ibd, no.call = mean.no.call)
        
                                        #calculate mean lengths of segments
        
                                        # to do this, first use ibdhap.seg.lengths to get a list of ibd.states and 
                                        # their respective lengths
        
        seg.lengths<-NULL
        ibd.states<- NULL
        for( icol in 1:ncol(calls)){
            
            temp<-ibdhap.seg.lengths( calls[,icol], position = position )
            ibd.states <-  c(ibd.states, temp[,1]) #ibd.states
            seg.lengths <- c(seg.lengths,temp[,2]) #seg.lengths
            
        }
        
        
        
        states.ibd<-((ibd.states>0)&(ibd.states<no.ibd.ind))
        states.no.call<-(is.element(ibd.states,0))
        states.not.ibd<-(is.element(ibd.states,no.ibd.ind))
        
                                        #mean segment lengths
        mean.length<-   c(
            len.ibd = sum(seg.lengths[states.ibd])/sum(states.ibd),
            len.not.ibd = sum(seg.lengths[states.not.ibd])/sum(states.not.ibd),
            len.nocall = sum(seg.lengths[states.no.call])/sum(states.no.call)
            )
        
                                        #counts
        seg.counts<-c(ibd = sum(states.ibd), no.ibd=sum(states.not.ibd), no.call=sum(states.no.call))	   
        
        return(list( mean.prop = mean.prop, mean.length = mean.length, seg.counts = seg.counts))
        
}

