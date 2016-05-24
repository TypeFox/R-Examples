ConsistencyIdeal <- function(dataset,col.p,col.j,col.lik,id.recogn,type="both",scale.unit=TRUE,ncp=NULL,axes=c(1,2),nbsim=0,replace.na=FALSE,graph=TRUE){
################################################################################
senso.consist <- function(dataset,col.p,col.j,col.lik,id.recogn,consist="both",scale.unit=TRUE,ncp=5,axes=c(1,2),correct=TRUE,replace.na=FALSE,graph=TRUE){
    if (!consist %in% c("consumer","panel","both"))
        stop("Inconvenient 'consist' definition.")
    dataset[,col.p] <- as.factor(dataset[,col.p])
    product <- levels(dataset[,col.p])
    nbprod <- length(product)
    prod.name <- colnames(dataset)[col.p]
    dataset[,col.j] <- as.factor(dataset[,col.j])
    juge <- levels(dataset[,col.j])
    nbjuge <- length(juge)
    juge.name <- colnames(dataset)[col.j]
    descp <- dataset[,c(col.j,col.p)]
    id.seq <- grep(id.recogn,colnames(dataset))
    if (length(id.seq)<2)
        stop("Not convenient 'id.recogn' definition")
    intensity <- dataset[,id.seq-1]
    attribut.int <- colnames(intensity)
    ideal <- dataset[,id.seq]
    attribut.id <- colnames(ideal)
    nbatt <- length(attribut.int)
    if (!is.numeric(col.lik)){
        pos.lik <- NULL
        for (a in 1:ncol(dataset))
            if (colnames(dataset)[a]==col.lik)
                pos.lik=a
        col.lik=pos.lik
    }
    if (is.null(col.lik) || col.lik %in% c(id.seq,(id.seq-1)))
        stop("Inconvenient 'Liking variable' definition")
    lik.name <- colnames(dataset)[col.lik]
    liking <- as.matrix(dataset[,col.lik])
    colnames(liking) <- lik.name
    lik.data <- cbind(descp,liking)
    int.data <- cbind(descp,intensity)
    int.avg.save <- averagetable(int.data,formul=paste("~",prod.name,"+",juge.name,sep=""),firstvar=3)
    id.data <- cbind(descp,ideal)
    id.avg.save <- averagetable(id.data,formul=paste("~",juge.name,"+",prod.name,sep=""),firstvar=3)
    if (consist %in% c("panel","both") || correct==T)
        int.avg.j <- averagetable(int.data,formul=paste("~",juge.name,"+",prod.name,sep=""),firstvar=3)
    res <- vector("list")
    if (consist %in% c("panel","both")){
        int.avg <- int.avg.save
        id.avg <- id.avg.save
        colnames(int.avg.j) <- colnames(int.avg) <- attribut.id
        id.avg <- id.avg-int.avg.j
        if (any(rownames(int.avg) %in% rownames(id.avg)))
            rownames(int.avg) <- paste("P_",rownames(int.avg),sep="")
        lik.info <- as.matrix(product)
        colnames(lik.info) <- prod.name
        for (j in 1:nbjuge){
            lik.j <- lik.data[lik.data[,1]==juge[j],-1]
            lik.info <- merge(lik.info,lik.j,by=prod.name,all=T,sort=F)
            colnames(lik.info)[ncol(lik.info)] <- paste(colnames(dataset)[col.lik],juge[j],sep="_")
        }
        rownames(lik.info) <- lik.info[,1]
        lik.info <- lik.info[,-1]
        lik.info <- t(scale(lik.info,center=T,scale=F))
        rownames(lik.info) <- juge
        data1 <- id.avg
        data2 <- rbind(scale(id.avg,center=T,scale=F),scale(int.avg,scale=F))
        data3 <- merge(id.avg,lik.info,all=T,by=0,sort=F)
        rownames(data3) <- data3[,1]
        data3 <- data3[,-1]
        res.pca1 <- PCA(data1,scale.unit=scale.unit,ncp=ncp,graph=F)
        res.pca2 <- PCA(data2,ind.sup=(nrow(data1)+1):nrow(data2),scale.unit=scale.unit,ncp=ncp,graph=F)
        res.pca3 <- PCA(data3,quanti.sup=(ncol(data1)+1):ncol(data3),scale.unit=scale.unit,ncp=ncp,graph=F)
        dev.new()
        layout(matrix(1:2,1,2))
        plot.PCA(res.pca1,choix="ind",cex=0.8,axes=axes,new.plot=F)
        plot.PCA(res.pca1,choix="var",cex=0.8,axes=axes,new.plot=F)
        dev.new()
        layout(matrix(1:2,1,2))
        plot.PCA(res.pca2,choix="ind",cex=0.8,label="ind.sup",axes=axes,new.plot=F)
        plot.PCA(res.pca3,choix="var",cex=0.8,label="quanti.sup",axes=axes,new.plot=F)
        correl <- cor(res.pca2$ind.sup$coord,res.pca3$quanti.sup$coord)
        rownames(correl) <- paste(rownames(correl),"_ideal.senso",sep="")
        colnames(correl) <- paste(colnames(correl),"_ideal.hedo",sep="")
        res$panel$dataset$ideal <- round(scale(id.avg,center=T,scale=F),2)
        res$panel$dataset$perceived <- round(scale(int.avg,center=T,scale=F),2)
        res$panel$dataset$hedonic <- round(lik.info,2)
        res$panel$PCA.ideal <- res.pca1
        res$panel$PCA.ideal_hedo <- res.pca3
        res$panel$PCA.ideal_senso <- res.pca2
        res$panel$correlation <- correl
    }
    if (consist %in% c("consumer","both")){
        int.avg <- int.avg.save
        id.avg <- id.avg.save
        perclik.cor.res <- matrix(NA,nbjuge,nbatt)
        juge.cor.res <- matrix(NA,nbjuge,2)
        rownames(perclik.cor.res) <- rownames(juge.cor.res) <- juge
        colnames(perclik.cor.res) <- attribut.int
        colnames(juge.cor.res) <- c(juge.name,"p-values")
        colnames(id.avg) <- attribut.int
        if (correct)
            id.avg <- id.avg-int.avg.j
        for (j in 1:nbjuge){
            lik.data.j <- lik.data[lik.data[,1]==juge[j],3]
            if (sd(lik.data.j)==0){
                juge.cor.res[j,1] <- juge.cor.res[j,2] <- NA
            } else {
                for (a in 1:nbatt){
                    int.data.j <- int.data[int.data[,1]==juge[j],(a+2)]
                    perclik.cor.res[j,a] <- cor(lik.data.j,int.data.j)
                    if (replace.na)
                        if (is.na(perclik.cor.res[j,a]))
                            perclik.cor.res[j,a]=0
                }
                res.cor <- cor.test(as.matrix(perclik.cor.res[j,]),t(id.avg[j,]),alternative="greater")
                juge.cor.res[j,1] <- res.cor$estimate
                juge.cor.res[j,2] <- res.cor$p.value
            }
        }
        oo <- order(juge.cor.res[,1],decreasing=T)
        juge.cor.res <- juge.cor.res[oo,]
        if (graph){
            dev.new()
            plot(density(na.omit(juge.cor.res[,1])),main=paste("Relationship between Ideal, Intensity and Liking data","\n","Distribution of cor(ideal,cor(intensity,liking)) for the different ",juge.name,sep=""),xlab="Individual correlation coefficient r(z,r(y,h))")
        }
        res$conso$driver.lik <- magicsort(round(perclik.cor.res,3))
        res$conso$correlations <- round(juge.cor.res,3)
    }
    return(res)
}
################################################################################
hedo.consist <- function(dataset,col.p,col.j,col.lik,id.recogn,family.model="PCR",ncp=NULL,scale.unit=TRUE,nbsim=0,graph=F){
    options(contrasts=c("contr.sum","contr.sum"))
    dataset[,col.p] <- as.factor(dataset[,col.p])
    product <- levels(dataset[,col.p])
    nbprod <- length(product)
    dataset[,col.j] <- as.factor(dataset[,col.j])
    juge <- levels(dataset[,col.j])
    nbjuge <- length(juge)
    descp <- dataset[,c(col.j,col.p)]
    id.pos <- grep(id.recogn,colnames(dataset))
    if (length(id.pos)<2)
        stop("Not convenient 'id.recogn' definition")
    intensity <- dataset[,id.pos-1]
    attribut.int <- colnames(intensity)
    ideal <- dataset[,id.pos]
    attribut.id <- colnames(ideal)
    nbatt <- length(attribut.int)
    if (!is.numeric(col.lik)){
        pos.lik <- NULL
        for (a in 1:ncol(dataset))
            if (colnames(dataset)[a]==col.lik)
                pos.lik=a
        col.lik=pos.lik
    }
    if (is.null(col.lik) || col.lik %in% c(id.pos,(id.pos-1)))
        stop("Inconvenient 'Liking variable' definition")
    lik.name <- colnames(dataset)[col.lik]
    liking <- as.matrix(dataset[,col.lik])
    colnames(liking) <- lik.name
    int.data <- cbind(descp,intensity)
    id.data <- cbind(descp,ideal)
    hedo.data <- cbind(descp,liking)
    nbsim=nbsim+1
#    if (is.null(simulation)){
        simulation <- as.matrix(1:nbprod)
        if (nbsim>1)
            for (sim in 2:nbsim)
                simulation <- cbind(simulation,as.matrix(sample(1:nbprod,nbprod,replace=F)))
#    } else {
#        simulation[,1] <- as.matrix(1:nbprod)
#        nbsim=ncol(simulation)
#    }
    rownames(simulation) <- product
    colnames(simulation) <- paste("Sim.",1:nbsim,sep="")
    res.sim.tf=T
    if (nbsim==1)
        res.sim.tf=F
    analyse.r2 <- matrix(NA,nbjuge,1)
    analyse.r2aj <- matrix(NA,nbjuge,1)
    hedo.idm <- matrix(NA,nbjuge,nbsim)
    rownames(analyse.r2) <- rownames(analyse.r2aj) <- rownames (hedo.idm) <- juge
    colnames(analyse.r2) <- "R2"
    colnames(analyse.r2aj) <- "adjusted R2"
    colnames(hedo.idm) <-  paste("Ideal_Sim.",1:nbsim,sep="")
    hedo.ref <- hedo.id <- matrix(0,nbprod,0)
    rownames(hedo.ref) <- product
    rownames(hedo.id) <- paste("id_",product,sep="")
    if (is.null(ncp))
#        if (family.model=="PLS"){
#            ncp=1
#        } else if (family.model=="Danzart"){
#            ncp=2
#        } else if (family.model=="PCR"){
            ncp=min(5,nbprod-1)
#        }
#    if (family.model=="Danzart"){
#        res.reg <- res.sim <- matrix(0,nbjuge,5)
#        rownames(res.reg) <- rownames(res.sim) <- juge
#        res.sim.sim <- matrix(0,nbsim,5)
#        rownames(res.sim.sim) <- paste("Sim.",1:nbsim,sep="")
#        colnames(res.reg) <- colnames(res.sim) <- colnames(res.sim.sim) <- c("Dim1","Dim2","Dim1_quad","Dim2_quad","Dim1Dim2")
#    } else if (family.model=="PCR"){
        res.reg <- res.sim <- matrix(0,nbjuge,ncp)
        rownames(res.reg) <- rownames(res.sim) <- juge
        res.sim.sim <- matrix(0,nbsim,ncp)
        rownames(res.sim.sim) <- paste("Sim.",1:nbsim,sep="")
        colnames(res.reg) <- colnames(res.sim) <- colnames(res.sim.sim) <- paste("Dim",1:ncp,sep="")
#    }
    juge.remove <- NULL
    for (j in 1:nbjuge){
#        if (nbsim > 1)
#            print(paste(colnames(dataset)[col.j],": ",j,"/",nbjuge,sep=""))
        int.j <- int.data[int.data[,1]==juge[j],]
        rownames(int.j) <- int.j[,2]
        int.j <- int.j[,-c(1,2)]
        id.j <- id.data[id.data[,1]==juge[j],-c(1,2)]
        rownames(id.j) <- paste("id_",id.data[id.data[,1]==juge[j],2],sep="")
        colnames(id.j) <- colnames(int.j)
        id.jm <- t(as.matrix(apply(id.j,2,mean,na.action="na.omit")))
        if (any(is.na(id.jm)))
            for (i in 1:length(id.jm))
                if (is.na(id.jm[1,i]))
                    id.jm[1,i] <- mean(na.omit(id.j[,i]))
        rownames(id.jm) <- "Ideal"
        hedo.j <- hedo.data[hedo.data[,1]==juge[j],]
        hedo.j.rn <- hedo.j[,2]
        hedo.j.cn <- colnames(hedo.j)[3]
        hedo.j <- as.matrix(hedo.j[,-c(1,2)])
        rownames(hedo.j) <- hedo.j.rn
        colnames(hedo.j) <- hedo.j.cn
        hedo.ref <- merge(hedo.ref,hedo.j,all=T,by=0,sort=F)
        hedo.ref.rn <- hedo.ref[,1]
        hedo.ref <- as.matrix(hedo.ref[,-1])
        rownames(hedo.ref) <- hedo.ref.rn
        colnames(hedo.ref)[ncol(hedo.ref)] <- juge[j]
        if (!sd(hedo.j[,1])==0){
#            if (family.model=="PLS"){
#                int.mean <- apply(int.j,2,mean)
#                int.j <- sweep(int.j,2,int.mean,FUN="-")
#                id.j <- sweep(id.j,2,int.mean,FUN="-")
#                id.jm <- sweep(id.jm,2,int.mean,FUN="-")
#                if (scale.unit){
#                    int.sd <- apply(int.j,2,sd)
#                    int.j <- sweep(int.j,2,int.sd,FUN="/")
#                    id.j <- sweep(id.j,2,int.sd,FUN="/")
#                    id.jm <- sweep(id.jm,2,int.sd,FUN="/")
#                }
#                col.rmv <- c(NULL)
#                for (a in 1:ncol(int.j))
#                    if (any(is.na(int.j[,a])))
#                        col.rmv <- c(col.rmv,a)
#                if (!is.null(col.rmv)){
#                    int.j <- int.j[,-col.rmv]
#                    id.jm <- id.jm[,-col.rmv]
#                    id.jm <- t(as.data.frame(id.jm))
#                }
#                for (sim in 1:nbsim){
#                    hedo.j.sim <- scale(hedo.j,center=T,scale=F)
#                    rownames(hedo.j.sim) <- rownames(hedo.j)[as.vector(simulation[,sim])]
#                    data.pls <- merge(int.j,hedo.j.sim,all=T,by=0,sort=F)
#                    rownames(data.pls) <- data.pls[,1]
#                    data.pls <- data.pls[,-1]
#                    res.pls <- plsr(as.formula(paste(colnames(data.pls)[ncol(data.pls)],"~.",sep="")),data=data.pls,ncomp=ncp,method="oscorespls",validation="LOO")
#                    if (sim==1){
#                        data.cor <- merge(hedo.j,mean(hedo.j)+predict(res.pls)[,,ncp],all=T,sort=F,by=0)[,-1]
#                        analyse.r2[j,1] <- cor(data.cor[,1],data.cor[,2])^2
#                        pred.id <- as.matrix(mean(hedo.j)+predict(res.pls,ncomp=ncp,newdata=id.j))
#                        rownames(pred.id) <- rownames(id.j)
#                        hedo.id <- merge(hedo.id,pred.id,all=T,by=0,sort=F)
#                        hedo.id.rn <- hedo.id[,1]
#                        hedo.id <- as.matrix(hedo.id[,-1])
#                        rownames(hedo.id) <- hedo.id.rn
#                        colnames(hedo.id)[ncol(hedo.id)] <- juge[j]
#                    }
#                    hedo.idm[j,sim] <- mean(hedo.j)+predict(res.pls,ncomp=ncp,newdata=id.jm)
#                }
#            } else if (family.model=="Danzart"){
#                data.pca <- rbind(int.j,id.jm)
#                data.pca <- rbind(data.pca,id.j)
#                res.pca <- PCA(data.pca,ind.sup=(nbprod+1):nrow(data.pca),scale.unit=scale.unit,graph=F,ncp=2)
#                score <- cbind(res.pca$ind$coord[,1:2],res.pca$ind$coord[,1]^2,res.pca$ind$coord[,2]^2,res.pca$ind$coord[,1]*res.pca$ind$coord[,2])
#                colnames(score) <- c("Dim1","Dim2","Dim1_quad","Dim2_quad","Dim1Dim2")
#                score.idm <- cbind.data.frame(res.pca$ind.sup$coord[,1:2],res.pca$ind.sup$coord[,1]^2,res.pca$ind.sup$coord[,2]^2,res.pca$ind.sup$coord[,1]*res.pca$ind.sup$coord[,2])
#                colnames(score.idm) <- colnames(score)
#                for (sim in 1:nbsim){
#                    hedo.j.sim <- hedo.j
#                    rownames(hedo.j.sim) <- rownames(hedo.j)[as.vector(simulation[,sim])]
#                    data.danzart <- merge(score,hedo.j.sim,all=T,by=0,sort=F)
#                    rownames(data.danzart) <- data.danzart[,1]
#                    data.danzart <- data.danzart[,-1]
#                    res.regbest <- reg.best(data.danzart[,ncol(data.danzart)],data.danzart[,-ncol(data.danzart)])$best
#                    if (rownames(res.regbest$coefficients)[1]=="(Intercept)"){
#                        effets <- rownames(res.regbest$coefficients)[-1]
#                        coeff <- res.regbest$coefficients[-1,4]
#                    } else {
#                        effets <- rownames(res.regbest$coefficients)
#                        coeff <- res.regbest$coefficients[,4]
#                    }
#                    if (min(coeff)<=0.05){
#                        if (length(effets)>=1){
#                            formul <- paste(colnames(data.danzart)[ncol(data.danzart)],"~",paste(effets,collapse="+"),sep="")
#                        } else if (length(effets)==1){
#                            formul <- paste(colnames(data.danzart)[ncol(data.danzart)],"~",effets,sep="")
#                        }
#                        for (l1 in 1:length(effets)){
#                            if (coeff[l1]<=0.05)
#                                for (l2 in 1:ncol(res.reg))
#                                    if (effets[l1]==colnames(res.reg)[l2])
#                                        if (sim==1){
#                                            res.reg[j,l2]=res.reg[j,l2]+1
#                                        } else {
#                                            res.sim.sim[sim,l2]=res.sim.sim[sim,l2]+1
#                                            res.sim[j,l2]=res.sim[j,l2]+1
#                                        }
#                        }
#                    } else {
#                        formul <- paste(colnames(data.danzart)[ncol(data.danzart)],"~1",sep="")
#                        effets <- 1
#                    }
#                    res.danzart <- aov(as.formula(formul),data=data.danzart)
#                    if (sim==1){
#                        if (effets[1]==1){
#                            analyse.r2[j,1]=NA
#                            analyse.r2aj[j,1]=NA
#                            hedo.idm[j,sim]=NA
#                            pred.id <- matrix(NA,nrow(score.idm[-1,]),1)
#                            rownames(pred.id) <- rownames(id.j)
#                            hedo.id <- merge(hedo.id,pred.id,all=T,by=0,sort=F)
#                            hedo.id.rn <- hedo.id[,1]
#                            hedo.id <- as.matrix(hedo.id[,-1])
#                            rownames(hedo.id) <- hedo.id.rn
#                            colnames(hedo.id)[ncol(hedo.id)] <- juge[j]
#                        } else {
#                            analyse.r2[j,1] <- summary.lm(res.danzart)$r.squared
#                            analyse.r2aj[j,1] <- summary.lm(res.danzart)$adj.r.squared
#                            hedo.idm[j,sim] <- predict(res.danzart,ncomp=ncp,newdata=as.data.frame(t(score.idm[1,])))
#                            pred.id <- as.matrix(predict(res.danzart,ncomp=ncp,newdata=score.idm[-1,]))
#                            rownames(pred.id) <- rownames(id.j)
#                            hedo.id <- merge(hedo.id,pred.id,all=T,by=0,sort=F)
#                            hedo.id.rn <- hedo.id[,1]
#                            hedo.id <- as.matrix(hedo.id[,-1])
#                            rownames(hedo.id) <- hedo.id.rn
#                            colnames(hedo.id)[ncol(hedo.id)] <- juge[j]
#                        }
#                    } else {
#                        hedo.idm[j,sim] <- predict(res.danzart,ncomp=ncp,newdata=score.idm[1,])
#                    }
#                }
#            } else if (family.model=="PCR"){
                data.pca <- rbind(int.j,id.jm)
                data.pca <- rbind(data.pca,id.j)
                res.pca <- PCA(data.pca,ind.sup=(nbprod+1):nrow(data.pca),scale.unit=scale.unit,graph=F,ncp=ncp)
                score <- res.pca$ind$coord[,1:ncp]
                colnames(score) <- paste("Dim",1:ncp,sep="")
                score.idm <- res.pca$ind.sup$coord[,1:ncp]
                colnames(score.idm) <- colnames(score)
                for (sim in 1:nbsim){
                    hedo.j.sim <- hedo.j
                    rownames(hedo.j.sim) <- rownames(hedo.j)[as.vector(simulation[,sim])]
                    data.pcr <- merge(score,hedo.j.sim,all=T,by=0,sort=F)
                    rownames(data.pcr) <- data.pcr[,1]
                    data.pcr <- data.pcr[,-1]
                    res.regbest <- RegBest(data.pcr[,ncol(data.pcr)],data.pcr[,-ncol(data.pcr)])$best
                    if (rownames(res.regbest$coefficients)[1]=="(Intercept)"){
                        effets <- rownames(res.regbest$coefficients)[-1]
                        coeff <- res.regbest$coefficients[-1,4]
                    } else {
                        effets <- rownames(res.regbest$coefficients)
                        coeff <- res.regbest$coefficients[,4]
                    }
                    if (min(coeff)<=0.05){
                        if (length(effets)>=1){
                            formul <- paste(colnames(data.pcr)[ncol(data.pcr)],"~",paste(effets,collapse="+"),sep="")
                        } else if (length(effets)==1){
                            formul <- paste(colnames(data.pcr)[ncol(data.pcr)],"~",effets,sep="")
                        }
                        for (l1 in 1:length(effets)){
                            if (coeff[l1]<=0.05)
                                for (l2 in 1:ncol(res.reg))
                                    if (effets[l1]==colnames(res.reg)[l2])
                                        if (sim==1){
                                            res.reg[j,l2]=res.reg[j,l2]+1
                                        } else {
                                            res.sim.sim[sim,l2]=res.sim.sim[sim,l2]+1
                                            res.sim[j,l2]=res.sim[j,l2]+1
                                        }
                        }
                    } else {
                        effets=1
                        formul <- paste(colnames(data.pcr)[ncol(data.pcr)],"~1",sep="")
                    }
                    res.pcr <- aov(as.formula(formul),data=data.pcr)
                    if (sim==1){
                        if (effets[1]==1){
                            analyse.r2[j,1]=NA
                            analyse.r2aj[j,1]=NA
                            hedo.idm[j,sim]=NA
                            pred.id <- matrix(NA,nrow(score.idm[-1,]),1)
                            rownames(pred.id) <- rownames(id.j)
                            hedo.id <- merge(hedo.id,pred.id,all=T,by=0,sort=F)
                            hedo.id.rn <- hedo.id[,1]
                            hedo.id <- as.matrix(hedo.id[,-1])
                            rownames(hedo.id) <- hedo.id.rn
                            colnames(hedo.id)[ncol(hedo.id)] <- juge[j]
                        } else {
                            analyse.r2[j,1] <- summary.lm(res.pcr)$r.squared
                            analyse.r2aj[j,1] <- summary.lm(res.pcr)$adj.r.squared
                            hedo.idm[j,sim] <- predict(res.pcr,ncomp=ncp,newdata=as.data.frame(t(score.idm[1,])))
                            pred.id <- as.matrix(predict(res.pcr,ncomp=ncp,newdata=as.data.frame(score.idm[-1,])))
                            rownames(pred.id) <- rownames(id.j)
                            hedo.id <- merge(hedo.id,pred.id,all=T,by=0,sort=F)
                            hedo.id.rn <- hedo.id[,1]
                            hedo.id <- as.matrix(hedo.id[,-1])
                            rownames(hedo.id) <- hedo.id.rn
                            colnames(hedo.id)[ncol(hedo.id)] <- juge[j]
                        }
                    } else {
                        hedo.idm[j,sim] <- predict(res.pcr,ncomp=ncp,newdata=as.data.frame(t(score.idm[1,])))
                    }
                }
#            } else {
#                stop("Not convenient 'family.model' definition: it should be 'PLS', 'Danzart' or 'PCR'...")
#            }
        } else {
            hedo.id <- cbind(hedo.id,matrix(NA,nrow(hedo.id),1))
            colnames(hedo.id)[ncol(hedo.id)] <- juge[j]
        }
    }
    data.boxprod <- matrix(0,0,2)
    colnames(data.boxprod) <- c("type","Hedonic consistency (panel)")
    for (j in 1:nbjuge){
        data.boxprod <- rbind(data.boxprod,cbind(as.matrix(rep("Actual product",nbprod)),as.matrix(hedo.ref[,j])))
        data.boxprod <- rbind(data.boxprod,cbind(as.matrix(rep("Ideal product",nbprod)),as.matrix(hedo.id[,j])))
    }
    mean.data.boxprod <- matrix(0,1,2)
    colnames(mean.data.boxprod) <- c("Actual product","Ideal product")
    rownames(mean.data.boxprod) <- "Mean"
    mean.data.boxprod[1,1] <- mean(na.omit(hedo.ref))
    mean.data.boxprod[1,2] <- mean(na.omit(hedo.id))
    data.boxprod <- as.data.frame(data.boxprod)
    data.boxprod[,2] <- as.numeric(as.character(data.boxprod[,2]))
    for (i in 1:nrow(data.boxprod))
        if (!is.na(data.boxprod[i,2])){
            if (data.boxprod[i,2]>=(max(hedo.ref)+3))
                data.boxprod[i,2]=max(hedo.ref)+3
            if (data.boxprod[i,2]<=-1)
                data.boxprod[i,2]=-1
        }
    boxprod(data.boxprod,col.p=1,firstvar=2,numr=1,numc=1)
    hedo.idm.cr <- matrix(NA,nbjuge,1)
    rownames(hedo.idm.cr) <- juge
    colnames(hedo.idm.cr) <- "Ideal"
    for (j in 1:nbjuge){
        hedo.j.mean <- mean(na.omit(hedo.ref[,j]))
        hedo.j.sd <- sd(na.omit(hedo.ref[,j]))
        if (!is.na(hedo.idm[j,1])){
            if ((hedo.idm[j,1]-hedo.j.mean)/hedo.j.sd >= 3){
                hedo.idm.cr[j,1] <- 3
            } else {
                hedo.idm.cr[j,1] <- (hedo.idm[j,1]-hedo.j.mean)/hedo.j.sd
            }
        }
    }
    dev.new()
    plot(analyse.r2,hedo.idm.cr,type="n",xlim=c(0,1),xlab="Rsquare",ylab="Standardized liking scores associated to the ideals",main="Standardized liking potential of the ideal products")
    for (j in 1:nbjuge)
        if (!is.na(hedo.idm.cr[j,1]))
            text(analyse.r2[j,1],hedo.idm.cr[j,1],label=juge[j],cex=0.7)
    abline(h=0,col="red3")
    abline(v=0.5,col="red3")
    abline(h=0.5,col="blue3",lty=3)
    compte.j <- matrix(0,nbjuge,1)
    rownames(compte.j) <- juge
    colnames(compte.j) <- "Count"
    for (j in 1:nbjuge){
        nb.na=0
        for (p in 1:nbprod)
            if (!is.na(hedo.id[p,j])){
                if (hedo.id[p,j]>=hedo.ref[p,j])
                    compte.j[j,1]=compte.j[j,1]+1
            } else {
                nb.na=nb.na+1
            }
        if (nb.na==nbprod)
            compte.j[j,1]=NA
    }
    dev.new()
    hist(compte.j,breaks=seq(0,nbprod,1),xlab="# Potential liking (ideal product) > liking score (actual product)",main=paste("Number of products with","\n","liking potential (ideal) > liking scores (actual product)",sep=""))
    if (res.sim.tf){
        p.val <- matrix(0,nbjuge,1)
        rownames(p.val) <- juge
        colnames(p.val) <- "Pvalue"
        for (j in 1:nbjuge)
            for (sim in 2:nbsim)
                if (!is.na(hedo.idm[j,1])){
                    if (hedo.idm[j,1]<hedo.idm[j,sim])
                        p.val[j,1]=p.val[j,1]+1
                } else {
                    p.val[j,1]=NA
                }
        p.val <- p.val/(nbsim-1)
        if (graph){
            dev.new()
            layout(matrix(1:6,2,3))
            for (j in 1:nbjuge){
                if (!is.na(hedo.idm[j,1])){
                    plot(density(na.omit(hedo.idm[j,])),xlim=c(0,max(liking)+3),xlab="Ideal liking score",main=juge[j])
                    abline(v=hedo.idm[j,1],col="red3")
                } else {
                    plot(c(1:5),c(1:5),type="n",xlim=c(0,max(liking)+3),xlab="Ideal liking score",main=juge[j])
                }
                if (ceiling(j/6)==floor(j/6) && !j==nbjuge){
                    dev.new()
                    layout(matrix(1:6,2,3))
                }
            }
        }
        dev.new()
        layout(matrix(1:2,1,2))
        plot(density(na.omit(p.val)),xlim=c(0,1),xlab="P-value",main=paste("Distribution of the individual p-values (",family.model,")",sep=""))
        abline(v=0.05,col="green3")
        text(0.075,0.1,label="5%",col="green3",cex=0.7)
        abline(v=0.1,col="blue3")
        text(0.13,0.1,label="10%",col="blue3",cex=0.7)
        hist(p.val,breaks=seq(0,1,0.025),xlab="P-value",main=paste("Distribution of the individual p-values (",family.model,")",sep=""))
        abline(v=0.05,col="green3")
        abline(v=0.1,col="blue3")
        dev.new()
        hedo.idm.plot <- as.matrix(hedo.idm[,1])
        for (j in 1:nbjuge)
            if (!is.na(hedo.idm[j,1]))
                if (hedo.idm[j,1]>max(liking)+3)
                    hedo.idm.plot[j,1]=max(liking)+3
        plot(analyse.r2,hedo.idm.plot,type="n",xlim=c(0,1),ylim=c(0,max(liking)+3),xlab="R2",ylab="Ideal liking score",main=paste("Liking potential of the ideal product","\n","in function of the quality of the model",sep=""))
        for (j in 1:nbjuge)
            if (!is.na(p.val[j,1]))
                if (p.val[j,1]<=0.05){
                    text(analyse.r2[j,1],hedo.idm.plot[j,1],label=juge[j],col="green3",cex=0.7)
                } else if (p.val[j,1]<=0.1){
                    text(analyse.r2[j,1],hedo.idm.plot[j,1],label=juge[j],col="blue3",cex=0.7)
                } else {
                    text(analyse.r2[j,1],hedo.idm.plot[j,1],label=juge[j],col="red3",cex=0.7)
                }
        legend("bottom",legend=c("p-val>10%","10%>p-val>5%","5%p-val"),text.col=c("red3","blue3","green3"),bty="n",horiz=T,title="Significance",title.col="black",cex=0.7)
#        if (family.model=="Danzart" || family.model=="PCR"){
            sum.reg <- 100*apply(res.reg,2,sum)/nbjuge
            sum.sim <- 100*apply(res.sim,2,sum)/((nbsim-1)*nbjuge)
            dev.new()
            layout(matrix(1:2,1,2))
            x.pos <- barplot(sum.reg,ylim=c(0,100),ylab="Significanace of the effect",axes=F,axisnames=F,main=paste("Model (",family.model,")",sep=""))
            abline(h=mean(sum.sim),col="red3")
            text(x.pos[1],mean(sum.sim)+2,label=round(mean(sum.sim),2),col="red3",cex=0.7)
            axis(2)
            axis(1,at=x.pos,labels=colnames(res.reg))
#            vps <- baseViewports()
#            pushViewport(vps$inner, vps$figure, vps$plot)
#            grid.text(colnames(res.reg),x=unit(x.pos,"native"),y=unit(-1,"lines"),just="right",rot=45,default.units="native",gp=gpar(fontsize=8))
#            popViewport(3)
            x.pos <- barplot(sum.sim,ylim=c(0,100),ylab="Significance of the effect",axes=F,axisnames=F,main=paste("Simulations (",family.model,")",sep=""))
            abline(h=mean(sum.sim),col="red3")
            text(x.pos[1],mean(sum.sim)+2,label=round(mean(sum.sim),2),col="red3",cex=0.7)
            axis(2)
            axis(1,at=x.pos,labels=colnames(res.sim))
#            vps <- baseViewports()
#            pushViewport(vps$inner, vps$figure, vps$plot)
#            grid.text(colnames(res.sim),x=unit(x.pos,"native"),y=unit(-1,"lines"),just="right",rot=45,default.units="native",gp=gpar(fontsize=8))
#            popViewport(3)
#        }
    }
    res <- list()
#    res$family.mod <- family.model
    res$R2 <- analyse.r2
    res$R2aj <- analyse.r2aj
    res$hedo$product <- hedo.ref
    res$hedo$ideal.product <- hedo.id
    res$hedo$avg.ideal <- hedo.idm
    res$hedo$relativ.ideal <- hedo.idm.cr
    res$hedo$compte <- compte.j
    if (res.sim.tf){
        res$simulation$hedo <- hedo.idm
        res$simulation$pvalue <- p.val
        res$simulation$simul <- simulation
#        if (family.model=="Danzart" || family.model=="PCR"){
            res$model$regression <- res.reg
            res$model$simulation <- res.sim
            res$model$simulation.bysim <- res.sim.sim
#        }
    }
    return(res)
}
################################################################################
    if (type %in% c("sensory","both"))
        res.senso <- senso.consist(dataset=dataset,col.p=col.p,col.j=col.j,col.lik=col.lik,id.recogn=id.recogn,scale.unit=scale.unit,ncp=ncp,axes=axes,graph=graph,correct=TRUE,replace.na=replace.na,consist="both")
    if (type %in% c("hedonic","both")){
        if (is.null(ncp))
            ncp=5
        res.hedo <- hedo.consist(dataset=dataset,col.p=col.p,col.j=col.j,col.lik=col.lik,id.recogn=id.recogn,scale.unit=scale.unit,ncp=ncp,graph=graph,nbsim=nbsim)
    }
    res <- vector("list",2)
    names(res) <- c("Senso","Hedo")
    if (type %in% c("sensory","both"))
        res[[1]] <- res.senso
    if (type %in% c("hedonic","both"))
        res[[2]] <- res.hedo
    return (res)
}