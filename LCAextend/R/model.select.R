model.select <-
function(ped,distribution,trans.const=TRUE,optim.param,optim.probs.indic=c(TRUE,TRUE,TRUE,TRUE),famdep=TRUE,selec="bic",H=5,K.vec=1:7,
                        tol=0.001,x=NULL,var.list=NULL)
{
    fam <- ped[,1]
    id <- ped[,2]
    dad <- ped[,3]
    mom <- ped[,4]
    sex <- ped[,5]
    status <- ped[,6]

	y <- as.matrix(ped[,-(1:6)])

    if(distribution=="normal") init.model <- init.norm
    if(distribution=="multinomial") init.model <- init.ordi

    if((selec=="cross")&length(K.vec)>1)
    {
        if(length(unique(fam))<H) stop("cross-validation model selection cannot be performed with fewer than H families\n")
        per.fam <- NULL
        for(f in unique(fam)) per.fam[f] <- sum(fam==f)
        fam.sorted <- unique(fam)[sort(per.fam,decreasing=TRUE,index.return=TRUE)$ix]
        group.fam <- matrix(NA,ceiling(length(unique(fam))/H),H)
        for(i in 1:floor(length(unique(fam))/H))
        {
            if(i%%2) from.to <- (1+(i-1)*H):(i*H)
            else from.to <- (i*H):(1+(i-1)*H)
            group.fam[i,] <- fam.sorted[from.to]
        }
        if(ceiling(length(unique(fam))/H)!=floor(length(unique(fam))/H))
        {
            remainder <- rep(0,H)
            remainder[1:(length(unique(fam))-floor(length(unique(fam))/H)*H)] <- fam.sorted[(floor(length(unique(fam))/H)*H+1):length(unique(fam))]
            group.fam[ceiling(length(unique(fam))/H),] <- remainder[(H+1)-sort(apply(matrix(per.fam[group.fam[1:floor(length(unique(fam))/H),]],
                                                                                floor(length(unique(fam))/H),H),2,sum),index.return=TRUE,decreasing=TRUE)$ix]
        }
    }

    if(selec=="bic")
    {
        res.lca <- list()
        ll <- bic <- NULL
    }
    if(selec=="cross") ll.valid <- NULL

    for(model in K.vec)
    {
        cat("model ",model,"\n")
	
		if(distribution=="normal"|distribution=="multinomial") param <- init.model(as.matrix(y[status==2,]),model,x,var.list)

		probs <- list()
        probs$p <- rep(1/model,times=model)
        probs$p0 <- 0.5
        if(famdep)
        {
            probs$p0connect <- rep(0.5,times=model)
            probs$p.trans <- init.p.trans(model,trans.const)
            probs$p.found <- probs$p.child <- 0.5
        }
        else probs$p.aff <- 0.5
        if(selec=="bic")
        {
            res.lca[[which(K.vec==model)]] <- lca.model(ped,probs,param,optim.param,fit=TRUE,optim.probs.indic,tol,x,var.list,famdep,modify.init=NULL)
            ll[which(K.vec==model)] <- res.lca[[which(K.vec==model)]]$ll
            bic[which(K.vec==model)] <- -2*ll[which(K.vec==model)]+n.param(y,model,trans.const,optim.param,optim.probs.indic,famdep)*log(nrow(ped))
        }
        if(selec=="cross")
        {
            if(length(K.vec)==1) res <- lca.model(ped,probs,param,optim.param,fit=TRUE,optim.probs.indic,tol,x,var.list,famdep,modify.init=NULL)
            else
            {
                ll.valid.h <- NULL
                for(h in 1:H)
                {
                    cat("model ",model,"test set ",h,"\n")
				  # Jordie: correction d'une erreur: mettre "as.matrix(x[fam %in% group.fam[,h],])" au lieu de 
				  # "x[fam %in% group.fam[,h]]" car x est une matrice.				  
                    res.estim <- try(lca.model(ped[fam%in%group.fam[,-h],],probs,param,optim.param,fit=TRUE,optim.probs.indic,tol,as.matrix(x[fam%in%group.fam[,-h],]),var.list,
                                                famdep,modify.init=NULL))
                    if(inherits(res.estim,"try-error"))
                    {
                        ll.valid.h[h] <- NA
                        cat("there is a problem with likelihood estimation in test set",h,"\n")
                    }
				  # Jordie: correction d'une erreur: mettre "as.matrix(x[fam %in% group.fam[,h],])" au lieu de 
				  # "x[fam %in% group.fam[,h]]" car x est une matrice.				  
                    else ll.valid.h[h] <- lca.model(ped[fam%in%group.fam[,h],],probs=res.estim$probs,param=res.estim$param,optim.param,fit=FALSE,
                                                    optim.probs.indic,tol,as.matrix(x[fam%in%group.fam[,h],]),var.list,famdep,modify.init=NULL)$ll
                }
                ll.valid[which(K.vec==model)] <- H*mean(ll.valid.h,na.rm=TRUE)
            }
        }        
    }
    if(selec=="bic")
    {
        res <- res.lca[[which.min(bic)]]
        res$ll <- ll
        res$bic <- bic
    }
    if((selec=="cross")&(length(K.vec)>1))
    {
        best.model <- K.vec[which.max(ll.valid)]

		if(distribution=="normal"|distribution=="multinomial") param <- init.model(as.matrix(y[status==2,]),best.model,x,var.list)

        probs <- list()
        probs$p <- rep(1/best.model,times=best.model)
        probs$p0 <- 0.5
        if(famdep)
        {
            probs$p0connect <- rep(0.5,times=best.model)
            probs$p.trans <- init.p.trans(best.model,trans.const)
            probs$p.found <- probs$p.child <- 0.5
        }
        else probs$p.aff <- 0.5
        res <- lca.model(ped,probs,param,optim.param,fit=TRUE,optim.probs.indic,tol,x,var.list,famdep,modify.init=NULL)
        res$ll.valid <- ll.valid
    }
    res
}

