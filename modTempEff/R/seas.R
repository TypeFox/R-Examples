seas<-function(x,ndx=stop("please, provide `ndx' in seas()")){
        r<-x
        #if(is.null(ndx)) ndx<-min(40,length(x)/4)
        attr(r,"ndx")<-ndx
        return(r)
        }
