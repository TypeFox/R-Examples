vcov.gcrq<-function(object,...){
        if(is.null(object$boot.coef)) stop(" 'vcov.gcrq' works only with boot")
        coef.boot<-object$boot.coef
        n.boot<-dim(coef.boot)[3]
        n.tau<- dim(coef.boot)[2]
        n.coef<-dim(coef.boot)[1]
#come calcolare la vcov complessiva???
#non saprei. Comunque fissato il percentile, puoi ottenere la vcov
#Una naive puoi ottenerla calcolando per ogni percentile la formula sandwich.. (vedi Koenker bandaid)
        VCOV<-NULL
        for(j in 1:n.tau) {
            m<-var(t(coef.boot[,j,]))
            colnames(m)<-rownames(m)<-rownames(object$coefficients)
            VCOV[[length(VCOV)+1]]<-m
            }
        names(VCOV)<-colnames(object$coefficients)
        return(VCOV)
        }
        
        
        