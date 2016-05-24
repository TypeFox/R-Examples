`run.analysis` <-
function(form, covariates, FDR = 0.1, norm.post.repl = FALSE, 
        norm.peaks = c("common", "all", "none"), normalization, add.norm = TRUE,  
        repl.method = "max", use.model = "lm", pval.fcn = "default", 
        lrg.only = TRUE, masses = NA, isotope.dist = 7, root.dir = ".", lrg.dir, 
        lrg.file = "lrg_peaks.RData", res.dir, res.file = "analyzed.RData", 
        overwrite = FALSE, use.par.file = FALSE, par.file = "parameters.RData", 
        bhbysubj = TRUE, subs, ...){
    if(missing(res.dir)){res.dir <- paste(root.dir, "/Results", sep="")}
    if(missing(lrg.dir)){lrg.dir <- paste(root.dir, "/Large_Peaks", sep="")}
    if(use.par.file){
        load(paste(root.dir, "/", par.file, sep=""))
        tmp <- match.call(expand.dots = FALSE)
        tmp[[1]] <- as.name("list")
        tmp <- eval(tmp)
        tmp[["..."]] <- NULL
        parameter.list <- extract.pars(par.file, root.dir)
        if(length(tmp) > 0){
            for(i in 1:length(tmp)){
                assign(names(tmp)[i],tmp[[i]])
                parameter.list[[names(tmp)[i]]] <- tmp[[i]]
            }
        }
    } else {
        parameter.list <- NA
    }
    if(!missing(normalization)){
        norm.post.repl <- (normalization == "postrepl")
        if(normalization %in% c("postbase","postrepl")){
            norm.peaks <- "all"
        } else {
            norm.peaks <- normalization
        }
    }
    if(class(use.model) == "character"){
        use.model <- get(use.model)
    }
    if(identical(masses, NA) && !lrg.only){
        stop("Either 'lrg.only' must be TRUE or 'masses' must be defined (or both)")
    }
    if(!file.exists(res.dir)){
        dir.create(res.dir)
    }
    if(!file.exists(paste(res.dir, "/", res.file, sep="")) || overwrite){
        load(paste(lrg.dir, "/", lrg.file, sep=""))
        if(!missing(subs)){
            clust.mat <- clust.mat[,subs,drop=F]
            lrg.mat <- lrg.mat[,subs,drop=F]
            amps <- amps[subs,,drop=F]
            covariates <- covariates[subs,,drop=F]
        }
        if(!norm.post.repl && norm.peaks != "none"){
            norm <- switch(norm.peaks,
                all = apply(clust.mat*lrg.mat,2,sum)/apply(lrg.mat,2,sum),
                common = apply(amps,1,mean)
            )
            if(add.norm){
                clust.mat <- clust.mat - rep(norm, each=dim(clust.mat)[1]) + mean(norm)
            } else {
                clust.mat <- as.matrix(clust.mat) %*% diag(mean(1/norm)/norm)
            }
        }
### to be changed        
        bysubjvar <- covariates$subj
        
        if(bhbysubj){
            num.lrg <- by(t(lrg.mat)+0, bysubjvar, function(x){apply(x,2,max)})
            num.lrg <- do.call(cbind, num.lrg)
            num.lrg <- apply(num.lrg,1,sum)
        } else {
            num.lrg <- apply(lrg.mat,1,sum)
        }
        if(!identical(repl.method,"none")){
            if(is.character(repl.method)){
                repl.method <- get(repl.method)
            }
            clust.mat <- do.call(cbind, by(as.data.frame(t(clust.mat)), bysubjvar, function(x)apply(x,2,repl.method)))
            lrg.mat <- do.call(cbind, by(as.data.frame(t(lrg.mat)), bysubjvar, function(x)apply(x,2,any)))
            amps <- do.call(rbind, by(as.data.frame(amps), bysubjvar, function(x)apply(x,2,repl.method)))
            if(length(cols <- grep("^[^:]*$",attributes(terms(form))$term.labels,value=TRUE))==1){
                covariates <- data.frame(c(by(covariates[,cols], bysubjvar, unique)))
                colnames(covariates) <- cols
            } else {
                covariates <- do.call(rbind, by(covariates[,cols], 
                    bysubjvar, unique))
            }
        }        
        bysubjvar <- covariates$subj
        if(norm.post.repl && norm.peaks != "none"){
            norm <- switch(norm.peaks,
                all = apply(clust.mat*lrg.mat,2,sum)/apply(lrg.mat,2,sum),
                common = apply(amps,1,mean)
            )
            if(add.norm){
                clust.mat <- clust.mat - rep(norm, each=dim(clust.mat)[1]) + mean(norm)
            } else {
                clust.mat <- as.matrix(clust.mat) %*% diag(mean(1/norm)/norm)
            }
        }        
        
        if(!identical(masses, NA)){
            wh <- .get.sp.masses(names(num.lrg), masses, isotope.dist)
            clust.mat <- clust.mat[wh,]
            num.lrg <- num.lrg[wh]
        }
    
        p.value <- rep(0, dim(clust.mat)[1])
        form <- update(form, Y~.)

	    if(identical(pval.fcn,"default")){
            if(identical(use.model, t.test)){
                pval.fcn <- function(x){x$p.value}
            } else if(identical(use.model, lm)){
                pval.fcn <- function(x){
                    args <- c(lapply(summary(x)$fstatistic, c), FALSE)
                    names(args) <- c("q", "df1", "df2", "lower.tail")
                    do.call(pf,args)}
            }
        }
        Delta = rep(NA, dim(clust.mat)[1])
        for(i in 1:length(p.value)){
            tmpdat <- data.frame(Y=t(clust.mat[i,,drop=FALSE]), covariates)
            colnames(tmpdat)[1] <- "Y"
            tmp <- use.model(form, dat=tmpdat, ...)
            p.value[i] <- pval.fcn(tmp)
            if(identical(use.model, t.test)){
                Delta[i] <- diff(tmp$estimate)
            }
        }
        which.sig <- data.frame(Delta, p.value, num.lrg, ord=1:length(p.value))

        which.sig <- which.sig[order(which.sig$p.val),]
        num.lrg.vals <- sort(unique(which.sig$num.lrg))
        tmp <- matrix(NA, ncol=length(num.lrg.vals), nrow=dim(which.sig)[1])
        colnames(tmp) <- paste("S", num.lrg.vals, sep="")
        for(k in num.lrg.vals){
            inds <- which(which.sig$num.lrg >= k)
            tmp[which.sig$ord[inds],match(k,num.lrg.vals)] <- 0
            sp <- .benj.hoch(which.sig$p.val[inds], FDR)
            if(sp > 0){
                tmp[which.sig$ord[inds][1:sp],match(k,num.lrg.vals)] <- 1
            }
            rm(sp,inds)
        }
        if(!lrg.only){
            tmp <- cbind(0, tmp)
            colnames(tmp)[1] <- "S0"
            sp <- .benj.hoch(which.sig$p.val, FDR)
            if(sp > 0){
                tmp[which.sig$ord[1:sp],1] <- 1
            }
        }

        which.sig <- which.sig[order(which.sig$ord),c("Delta","p.value","num.lrg")]
        which.sig <- data.frame(which.sig, tmp)
        sigs <- which.sig[apply(which.sig[,-(1:3),drop=FALSE]==1,1,any, na.rm=TRUE),]
        if(dim(sigs)[1]){
            sigs <- sigs[,c(TRUE,TRUE,TRUE,apply(sigs[,-(1:3),drop=FALSE]==1,2,any,na.rm=TRUE))]
            if(all(is.na(sigs$Delta))){
            	sigs$Delta <- c()
            }
        } else {
            sigs <- sigs[,1:3]
        }
        min.FDR <- sapply(1:max(num.lrg), function(x){
            tmp <- sort(which.sig$p.value[which.sig$num.lrg >= x])
            min(length(tmp)*tmp/(1:length(tmp)))
        })
        names(min.FDR) <- 1:length(min.FDR)

        save(amps,centers,clust.mat,min.FDR,sigs,which.sig,parameter.list,bysubjvar,
             file=paste(res.dir, "/", res.file, sep=""))
    } else {
        warning("Results file exists and overwrite = FALSE; no results file created")
    }
}
