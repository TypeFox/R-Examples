`calcfit2Diffrep` <-
function (C1, C2) 
{
    data.mat <- as.matrix(cbind(C1, C2))
    group <- factor(c(rep(1, dim(C1)[2]), rep(2,  dim(C2)[2])))
    design <- model.matrix(~-1 + group)
    fit = lmFit(data.mat, design)
    contrast.matrix<-makeContrasts("group2-group1",levels=design)
    rownames(contrast.matrix)=colnames(design)
    fit2 <- contrasts.fit(fit, contrast.matrix)
    fit2 <- eBayes(fit2)
    fit2
}

