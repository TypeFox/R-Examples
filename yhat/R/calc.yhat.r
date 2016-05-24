
calc.yhat<-function (lm.out,prec=3) 
{
    ifelse(length(lm.out$call) > 2, new.data <- eval(as.name(lm.out$call[[3]]), 
        parent.frame()), new.data <- model.frame(lm.out$call[[2]]))
    new.scale <- model.matrix(lm.out)
    v <- as.numeric(as.factor(row.names(new.scale)))
    #    new.data <- new.data[v, ]
    IV <- attr(lm.out$terms, "term.labels")
    DV <- dimnames(attr(lm.out$terms, "factors"))[[1]][1]
    IVx <- dimnames(attr(lm.out$terms, "factors"))[[1]][-1]
    new.data<-na.omit(new.data[,c(DV,IVx)])
    new.scale[, 1] <- lm.out$model[, 1]
    new.scale <- data.frame(new.scale)
    colnames(new.scale) <- c(DV, IV)
    betaOut <- coef(lm.out)[-1] * sapply(new.scale[IV], "sd")/sapply(new.scale[DV], "sd")
    structure.coef <- cor(na.omit(fitted.values(lm.out)), new.scale[IV])
    rownames(structure.coef)<-""
    apsOut<-aps(new.data,DV,IV)
    domOut<-dominance(apsOut)
    dombinOut<-dombin(domOut)
    comOut<-commonality(apsOut)
    rlwOut<-t(rlw(new.data,DV,IV))
    corOut<-cor(new.scale)
    corOut<-corOut[,1]
    corOut<-corOut[-1]
    a<-rbind(lm.out$coef[-1],betaOut,corOut,structure.coef,structure.coef^2,u<-comOut[c(1:length(betaOut)),1],corOut^2-u,
             domOut$CD,domOut$GD, corOut*betaOut,rlwOut)
    rownames(a)<-c("b","Beta","r","rs","rs2","Unique","Common",rownames(domOut$CD),"GenDom","Pratt","RLW")
    a<-round(t(a),digits=prec)
    a<-rbind(a,colSums(a))
    l<-nrow(a)
    rownames(a)[l]<-"Total"
    a[l,c(1:4)]<-NA
    apsanal<-cbind(comOut,domOut$DA)
    rownames(apsanal)<-rownames(domOut$DA)
    colnames(apsanal)[1]<-"Commonality" 
    apsanal<-rbind(apsanal,Total=colSums(apsanal))
    apsanal[nrow(apsanal),3:ncol(apsanal)]<-NA
    apsanal<-apsanal[,-3]
    opm<-a[-nrow(a),]
    for (i in 1:ncol(opm)){
      o<-order(abs(opm[,i]),decreasing=TRUE)
      opm[o,i]<-1:nrow(opm)
    }
    return(list(PredictorMetrics= a, OrderedPredictorMetrics=opm,PairedDominanceMetrics=dombinOut, APSRelatedMetrics=round(apsanal,digits=prec)) )
}           