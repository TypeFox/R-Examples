incon <-
function(data){ #last col : colclass attribute
            data <- as.data.frame(data)
            p <- ncol(data)
            n <- nrow(data)
            data[,p] <- as.numeric(data[,p])
            Class <- dim(table(data[,p])) #class level
            
            data.ex <- data[,1:(p-1)] # data : except class attribute
            if((p-1)==1){data.ex=as.matrix(data.ex)}    
            mat <- unique(data.ex)   #matching instance
            if(dim(mat)[1]==dim(data)[1]) return(inConRate=0)
            w <- which(duplicated(data.ex)==TRUE)#instance number(ID) of duplicated data 
                                           #-first duplicated data is removed in mat by unique function.
            dup <- unique(data.ex[w,]) #duplicated instance
            if(is.numeric(dup)==TRUE){dup=as.matrix(dup)}
            n.dup=nrow(dup) #number of duplicated instance

           num <- 0
           for(i in 1:n.dup){
                same <- apply(as.numeric(dup[i,1:(p-1)])==t(data.ex),2,prod) #To find same instance for all attribute.
                                                                      #same==1 : It means same instances except class.
                ix <- as.numeric(which(same==1)) # Same instance's instance number
                firs <- ix[1] #ix[1] is representative matching instances
                num <- as.numeric(ifelse(same==0,0,firs))+num  # dupliced data's instance number makes same.
                fnum <- ifelse(num==0,rownames(data),num) #final numbering
                fnum <- as.numeric(fnum)
                freq <- table(fnum,data[,p]) # class frequency for matching instance
           }                              #final table for inconsistency rate, for class
        insMat <- apply(freq,1,sum) # number of matched instance 
        insMax <- apply(freq,1,max) # max class in number of matched instance 
        inConCount <- insMat-insMax # inconsitency count
        inConRate <- sum(inConCount)/n #inconsistency rate of dataset -sum of inconsistency count
        return(inConRate)
}
