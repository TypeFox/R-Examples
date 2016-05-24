sampleDepth = function(rwl, stc=c(5,2,1)) {
                             temp = rbind(colnames(rwl),rwl)
                             NoTREES = function (rw){
                             rw[!is.na(rw)]<-substr(rw[1], 1, stc[1]+stc[2])
                             return(rw)
                             }
                             temp = apply(temp,2,NoTREES)
                             temp=temp[-1,]
                             NoCores = apply(rwl, 1, function(y) sum(!is.na(y)))
                             NoTrees = apply(temp, 1, function(y) length(unique(y[!is.na(y)])))
                             out=data.frame(NoTrees, NoCores)
                             return(out)
                             }

#sampleDepth(rwl, stc=c(5,2,1))