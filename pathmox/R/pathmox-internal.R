.pls.loads <-
function(X, Y.lvs, blocks)
{
    lvs <- length(blocks)
    mvs <- ncol(X)
    blocklist <- as.list(1:lvs)
    for (j in 1:lvs)
         blocklist[[j]] <- rep(j,blocks[j])
    blocklist <- unlist(blocklist)
    loads <- rep(NA, mvs)
    comu <- rep(NA, mvs)
    for (j in 1:lvs)
        loads[blocklist==j] <- cor(X[,blocklist==j], Y.lvs[,j])
    comu <- loads^2
    names(loads) <- colnames(X)  
    names(comu) <- colnames(X)
    res.loads <- list(loads, comu)
    return(res.loads)
}

.pls.paths <-
function(IDM, Y.lvs, plsr)
{
    lvs.names <- colnames(IDM)
    endo = rowSums(IDM)
    endo[endo!=0] <- 1  # vector indicating endogenous LVs
    innmod <- as.list(1:sum(endo))
    Path <- IDM
    residuals <- as.list(1:sum(endo))
    R2 <- rep(0,nrow(IDM))
    for (aux in 1:sum(endo)) 
    {
        k1 <- which(endo==1)[aux]    # index for endo LV
        k2 <- which(IDM[k1,]==1)     # index for indep LVs
        if (length(k2)>1 & plsr) {               
            path.lm <- .plsr1(Y.lvs[,k2], Y.lvs[,k1], nc=2)
            Path[k1,k2] <- path.lm$coeffs
            residuals[[aux]] <- path.lm$resid
            R2[k1] <- path.lm$R2[1]
            inn.val <- c(path.lm$R2[1], path.lm$cte, path.lm$coeffs)
            inn.lab <- c("R2", "Intercept", paste(rep("path_",length(k2)),names(k2),sep=""))
            names(inn.val) <- NULL
            innmod[[aux]] <- data.frame(concept=inn.lab, value=round(inn.val,4))
        }
        if (length(k2)==1 | !plsr) {
            path.lm <- summary(lm(Y.lvs[,k1] ~ Y.lvs[,k2]))
            Path[k1,k2] <- path.lm$coef[-1,1]
            residuals[[aux]] <- path.lm$residuals  
            R2[k1] <- path.lm$r.squared
            inn.val <- c(path.lm$r.squared, path.lm$coef[,1])
            inn.lab <- c("R2", "Intercept", paste(rep("path_",length(k2)),names(k2),sep=""))
            names(inn.val) <- NULL
            innmod[[aux]] <- data.frame(concept=inn.lab, value=round(inn.val,4))
        }
    }
    names(innmod) <- lvs.names[endo!=0]  
    names(R2) <- lvs.names
    res.paths <- list(innmod, Path, R2, residuals)
    return(res.paths)
}

.pls.weights <-
function(X, IDM, blocks, modes, scheme, tol, iter)
{
    lvs <- nrow(IDM)
    mvs <- ncol(X)
    sdv <- sqrt((nrow(X)-1)/nrow(X))   # std.dev factor correction
    blocklist <- as.list(1:lvs)
    for (j in 1:lvs)
         blocklist[[j]] <- rep(j,blocks[j])
    blocklist <- unlist(blocklist)
    # outer design matrix 'ODM' and matrix of outer weights 'W'
    ODM <- matrix(0, mvs, lvs)
    for (j in 1:lvs)
        ODM[which(blocklist==j),j] <- rep(1,blocks[j])
    W <- ODM %*% diag(1/(sd(X %*% ODM)*sdv),lvs,lvs)
    w.old <- rowSums(W)    
    w.dif <- 1
    itermax <- 1
    repeat 
    {            
        Y <- X %*% W  # external estimation of LVs 'Y'
        Y <- scale(Y) * sdv
        # matrix of inner weights 'e' 
        E <- switch(scheme, 
               "centroid" = sign(cor(Y) * (IDM + t(IDM))),
               "factor" = cor(Y) * (IDM + t(IDM)),
               "path" = .path.scheme(IDM, Y))
        Z <- Y %*% E  # internal estimation of LVs 'Z'
        # scaling Z
        Z <- Z %*% diag(1/(sd(Z)*sdv), lvs, lvs)
        # computing outer weights 'w'
        for (j in 1:lvs)
        {
            X.blok = X[,which(blocklist==j)] 
            if (modes[j]=="A")# reflective way
                ODM[which(blocklist==j),j] <- (1/nrow(X)) * Z[,j] %*% X.blok
            if (modes[j]=="B")# formative way
                ODM[which(blocklist==j),j] <- solve.qr(qr(X.blok),Z[,j])
        }
        W <- ODM
        w.new <- rowSums(W)                
        w.dif <- sum((abs(w.old) - abs(w.new))^2)  # difference of out.weights 
        if (w.dif<tol || itermax==iter) break
        w.old <- w.new
        itermax <- itermax + 1
    } # end repeat       
    W <- ODM %*% diag(1/(sd(X %*% ODM)*sdv),lvs,lvs)
    w.new <- rowSums(W)                
    names(w.new) <- colnames(X)
    dimnames(W) <- list(colnames(X),rownames(IDM))       
    res.ws <- list(w.new, W, itermax)
    if (itermax==iter) res.ws=NULL
    return(res.ws)
}

.plsr1 <-
function(x, y, nc=NULL, scaled=TRUE)
{
    # ============ checking arguments ============
    X <- as.matrix(x)
    Y <- as.matrix(y)
    n <- nrow(X)
    p <- ncol(X)
    if (is.null(nc))
        nc <- p
    # ============ setting inputs ==============
    if (scaled) Xx<-scale(X) else Xx<-scale(X,scale=F)
    if (scaled) Yy<-scale(Y) else Yy<-scale(Y,scale=F)
    X.old <- Xx
    Y.old <- Yy
    Th <- matrix(NA, n, nc)# matrix of X-scores
    Ph <- matrix(NA, p, nc)# matrix of X-loadings
    Wh <- matrix(NA, p, nc)# matrix of raw-weights
    Uh <- matrix(NA, n, nc)# matrix of Y-scores
    ch <- rep(NA, nc)# vector of y-loadings
    # ============ pls regression algorithm ==============
    for (h in 1:nc)
    {
        w.old <- t(X.old) %*% Y.old / sum(Y.old^2)
        w.new <- w.old / sqrt(sum(w.old^2)) # normalization
        t.new <- X.old %*% w.new
        p.new <- t(X.old) %*% t.new / sum(t.new^2) 
        c.new <- t(Y.old) %*% t.new / sum(t.new^2)
        u.new <- Y.old / as.vector(c.new)
        Y.old <- Y.old - t.new%*%c.new# deflate y.old
        X.old <- X.old - (t.new %*% t(p.new))# deflate X.old
        Th[,h] <- t.new
        Ph[,h] <- p.new
        Wh[,h] <- w.new
        Uh[,h] <- u.new
        ch[h] <- c.new
    }
    Ws <- Wh %*% solve(t(Ph)%*%Wh)# modified weights
    Bs <- as.vector(Ws %*% ch) # std beta coeffs    
    Br <- Bs * (rep(sd(Y),p)/apply(X,2,sd))   # beta coeffs
    cte <- as.vector(mean(y) - Br%*%apply(X,2,mean))# intercept
    y.hat <- X%*%Br+cte# y predicted
    resid <- as.vector(Y - y.hat)# residuals
    R2 <- as.vector(cor(Th, Yy))^2  # R2 coefficients    
    names(Br) <- colnames(X)
    names(resid) <- rownames(Y)
    names(y.hat) <- rownames(Y)
    names(R2) <- paste(rep("t",nc),1:nc,sep="")
    res <- list(coeffs=Br, coef.std=Bs, cte=cte, R2=R2[1:nc], resid=resid, y.pred=y.hat)    
    return(res)
}

.path.scheme <-
function(IDM, Y)
{
    lvs <- nrow(IDM)
    E <- IDM
    for (k in 1:lvs) 
    {
        if (length(which(IDM[k,]==1)) > 0)
            E[which(IDM[k,]==1),k] <- lm(Y[,k]~Y[,which(IDM[k,]==1)]-1)$coef
        if (length(which(IDM[,k]==1)) > 0)
            E[which(IDM[,k]==1),k] <- cor(Y[,k], Y[,which(IDM[,k]==1)])
    }                 
    return(E)
}

.nominal.split <-
function(v)
{
    # INPUT(S)
    # v: vector containing the nominal categories

    # PARAM(S)
    k = length(v)   # length of vector
    num.parts = (2^(k-1))-1   # number of binary partitions
    # stop limit of order of partitions
    if ((k%%2)==0) stop.parts=k/2 else stop.parts=(k-1)/2    
    part = as.list(1:num.parts)    # list for warehousing the partitions
    antipar = part   # complementary list for part

    # FUNCTION(S)
    # internal function for evaluating candidates splits
    eval.cand = function(x, y, aux)
    {
        # x: list of partitions
        # y: candidate split to be included in part
        # aux: num of parts already included in part
        acum = 0
        for (i in 1:(aux-1))
        {
             true.false = setequal(x[[i]], y)
             if (true.false) acum=acum+1 else acum=acum+0
        }
        return(acum)
    }
       
    # first order partitions
    aux = 1
    for (i in 1:k)
    {
         part[[aux]] = v[i]
         antipar[[aux]] = setdiff(v, part[[i]])
         aux = aux + 1
    }
    
    # partitions of order greater than one
    for (ord in 2:stop.parts)
    {
        if (stop.parts < ord)
            break
        # partitions of the immediate past order
        p1 = aux - choose(k, (ord-1))
        p2 = aux - 1
        part.ant = part[p1:p2]
        antipar.ant = antipar[p1:p2]    
        for (i in 1:length(part.ant))
        { 
             part1 = part.ant[[i]]
             anti1 = antipar.ant[[i]]
             for (j in 1:length(anti1))
             {
                  candidat = c(part1, anti1[j])
                  cand = eval.cand(part, candidat, aux)
                  if (cand==0)
                  {
                       part[[aux]] = candidat
                       antipar[[aux]] = setdiff(v, candidat)
                       aux = aux + 1
                  } else
                  {
                       next
                  }
             }
        }
    }    
    # in case that k is even, only select partitions until num.parts
    if (length(part) > num.parts)
    {
        part = part[1:num.parts]
        antipar = antipar[1:num.parts]
    }   
    parts = list(par1=part, par2=antipar)
    return(parts)
}

.ordinal.split <-
function(cat.exe)
{
    a = length(cat.exe)
    if (a == 2)  # binary variable
    {
        part = cat.exe[1]
        antipar = cat.exe[2]
    }
    if (a > 2)   # more than 2 categories
    {
        part = as.list(1:(a-1))
        antipar = part
        for (i in 1:(a-1))
        {
             part[[i]] = cat.exe[1:i]
             antipar[[i]] = setdiff(cat.exe, part[[i]])
        }
    }
    parts = list(par1=part, par2=antipar)
    return(parts)
}

.pls.basic <- 
function(DT, IDM, blocks, modes, scheme, scaled, tol, iter)
{   
    ##### PARAMETERS
    DM <- DT
    lvs <- nrow(IDM)
    mvs <- sum(blocks)
    blocklist <- as.list(1:lvs)
    for (j in 1:lvs)
         blocklist[[j]] <- rep(j,blocks[j])
    blocklist <- unlist(blocklist)
    endo = rowSums(IDM)
    endo[endo!=0] = 1
    plsr <- FALSE
    ### Variable Names
    lvs.names = colnames(IDM)
    mvs.names = colnames(DM)

    # apply the selected scaling
    if (scaled) {
        sd.X <- sqrt((nrow(DM)-1)/nrow(DM)) * apply(DM, 2, sd)
        X <- scale(DM, scale=sd.X)
    } else {
        X <- scale(DM, scale=FALSE)
    }
    dimnames(X) <- list(rownames(DM), mvs.names)

    # ==================== Stage 1: Iterative procedure ==================
    out.ws <- .pls.weights(X, IDM, blocks, modes, scheme, tol, iter)
    if (is.null(out.ws)) stop("The pls algorithm is non convergent") 
    out.weights <- round(out.ws[[1]], 4)
    cor.XY <- cor(X, X%*%out.ws[[2]])
    w.sig <- rep(NA,lvs)
    for (k in 1:lvs) 
         w.sig[k] <- ifelse(sum(sign(cor.XY[which(blocklist==k),k]))<=0,-1,1)
    Z.lvs <- X %*% out.ws[[2]] %*% diag(w.sig,lvs,lvs)
    Y.lvs <- Z.lvs
    if (!scaled) 
        Y.lvs <- DM %*% out.ws[[2]] %*% diag(w.sig,lvs,lvs)
    dimnames(Y.lvs) <- list(rownames(X), lvs.names)
    dimnames(Z.lvs) <- list(rownames(X), lvs.names)
    loadcomu <- .pls.loads(X, Y.lvs, blocks)    
    loads <- loadcomu[[1]]

    # ============ Stage 2: Path coefficients and total effects ==========
    pathmod <- .pls.paths(IDM, Y.lvs, plsr)
    innmod <- pathmod[[1]]
    Path <- pathmod[[2]]
    R2 <- pathmod[[3]]
    residuals <- pathmod[[4]]

    # ============================= Results ==============================
    model <- list(IDM=IDM, blocks=blocks, scheme=scheme, modes=modes, 
                  scaled=scaled, tol=tol, iter=iter)
    resul = list(out.weights=out.weights, loadings=loads, scores=Y.lvs,   
                 path.coefs=Path, R2=R2, residuals=residuals, model=model)
    resul
}

.sanchez.aluja <- 
function(Y.lvs, IDM, split.a, split.b)
{    
    lvs <- ncol(Y.lvs)	
    endo <- rowSums(IDM)		
    endo[endo!=0] <- 1  
    n.endo <- sum(endo)
    n <- nrow(Y.lvs)
    p.values <- rep(0,n.endo)		
    f.values <- rep(0,n.endo)
    df1.values <- rep(0,n.endo)
    df2.values <- rep(0,n.endo)
    s1s2.values <- rep(0,n.endo)
    who.endo <- which(endo==1)
    SCE0 <- 0
    SCE1 <- 0		
    Sa <- 0
    Sb <- 0
    p.sum <- 0

    for (i in 1:n.endo)
    {
        aux <- who.endo[i]
        indep <- which(IDM[aux,1:aux]==1)
        reg <- cbind(Y.lvs[,aux], Y.lvs[,indep])
        # SS0 calculation
        lm.res <- lm(Y.lvs[,aux] ~ Y.lvs[,indep])$residuals
        SS0 <- sum(lm.res^2)
        SCE0 <- SCE0 + SS0
        # regression matrices "Xa" and "Xb"
        Xa <- reg[split.a,]
        Xb <- reg[split.b,]
        lm.resa <- lm(Xa[,1] ~ Xa[,2:ncol(Xa)])$residuals
        lm.resb <- lm(Xb[,1] ~ Xb[,2:ncol(Xb)])$residuals
        SS1 <- sum(lm.resa^2) + sum(lm.resb^2)
        SCE1 <- SCE1 + SS1
        Sa <- Sa + lm.resa
        Sb <- Sb + lm.resb
        # test
        p <- length(indep)        
        p.sum <- p.sum + p
        df1 <- p
        df2 <- n - 2*p
        f <- ((n-2*p)/p) * ((SS0-SS1)/SS1)
        p.values[i] <- 1 - pf(f, df1, df2)
        f.values[i] <- f
        df1.values[i] <- df1
        df2.values[i] <- df2
        s1s2.values[i] <- round(sd(lm.resa)/sd(lm.resb),4)
    }

    n.sum <- n.endo * (length(split.a) + length(split.b))
    f.sum <- ((n.sum-2*p.sum)/p.sum)  *  ((SCE0 - SCE1)/SCE1)  	# f global
    df1.sum <- p.sum   						# gl num
    df2.sum <- n.sum - 2*p.sum					# gl den
    pval.sum <- 1 - pf(f.sum, df1.sum, df2.sum)	   		# pv
    s1s2.sum <- sd(Sa) / sd(Sb)
    F.global <- c(f.sum, df1.sum, df2.sum, pval.sum, s1s2.sum)   		

    F.partial <- cbind(f.values, df1.values, df2.values, p.values, s1s2.values) 
    rownames(F.partial) <- rownames(IDM)[who.endo]
    colnames(F.partial) <- c("f.stat", "df.num", "df.den", "p.val", "sa/sb")
    res <- list(F.global, F.partial)
    return(res)
}

.xexeloa <-
function(pls, DT, EXEV, type.exev, elemnod, nv, size, mox)
{
     # cuauhmaitl elems
     N <- nrow(EXEV)
     if (nv==0) elems<-1:length(elemnod) else elems<-which(elemnod==nv)
     indelem <- cbind(1:length(elems), elems)# identificador
     E.node <- EXEV[elems,]# elems en cuauhmaitl actual       
     exevs <- apply(E.node, 2, function(x) nlevels(as.factor(x)))
     list.cata <- as.list(1:length(exevs))   # categs a.cuauhmaitl
     list.catb <- as.list(1:length(exevs))   # categs b.cuauhmaitl
     Ftest.global <- matrix(NA, length(exevs), 5)# + chingon global
     Ftest.partial <- as.list(1:length(exevs))# + chingon partial

     # from plspm
     Y.lvs <- pls$scores
     IDM <- pls$model$IDM
     blocks <- pls$model$blocks
     modes <- pls$model$modes
     scheme <- pls$model$scheme
     scaled <- pls$model$scaled
     tol <- pls$model$tol
     iter <- pls$model$iter
     endo <- rowSums(IDM)
     endo[endo!=0] <- 1  
     n.endo <- sum(endo)

     # size limite
     if (size < 1)  size.limit=ceiling(size*N)  else  size.limit=size

     # for each exev
     for (i in 1:length(exevs))    
     {
          if (exevs[i] == 1)  # ya se uso
          {
               Ftest.global[i,4] <- 1   # sta bien chafa
               Ftest.partial[[i]] <- matrix(1,1,5)# sta bien chafa
               list.cata[[i]] <- NA# NAs a.cuauhmaitl
               list.catb[[i]] <- NA# NAs b.cuauhmaitl
               next   
          }  else   # se usa               
          {
               ### moxexeloa
               cat.exe <- unique(E.node[,i])  
               if (type.exev[i] == "ord") { 
                    bin.split <- .ordinal.split(sort(cat.exe))  
               } else {# nom
                    if (exevs[i]==2) {  
                        bin.split <- as.list(cat.exe)
                    } else { 
                        bin.split <- .nominal.split(cat.exe)  
                    }
               }
               ### sanchez-aluja moxexeloa
               Ftest.glo <- matrix(NA, length(bin.split[[1]]), 5)   # global
               Ftest.par <- as.list(1:length(bin.split[[1]]))      # partial
               for (aux in 1:length(bin.split[[1]]))   # p/c moxexeloa
               {
                     split.a <- which(E.node[,i] %in% bin.split[[1]][[aux]])# a.cuauhmaitl
                     split.b <- which(E.node[,i] %in% bin.split[[2]][[aux]])  # b.cuauhmaitl
                     n.a <- length(split.a)
                     n.b <- length(split.b)
                     if (n.a<=size.limit || n.b<=size.limit)# cuahuitl aturat
                     {
                          Ftest.glo[aux,4] <- 1   # sta bien chafa
                          Ftest.par[[aux]] <- matrix(1,1,5)# sta bien chafa
                          next  
                     } else# cuahuitl segueix
                     {
                          idm <- .sanchez.aluja(Y.lvs, IDM, split.a, split.b)
                          Ftest.glo[aux,] <- idm[[1]]   # global
                          Ftest.par[[aux]] <- idm[[2]]  # partial
                     }
               }
               ### el gallo dels exevs
               if (mox=="pathmox")
                   min.p <- which(Ftest.glo[,4]==min(Ftest.glo[,4]))[1]  
               if (mox=="techmox") 
               {
                   pvals.ave <- unlist(lapply(Ftest.par, function(x) exp(mean(log(x[,4]))) ))
                   min.p <- which(pvals.ave==min(pvals.ave))[1]
               }
               Ftest.global[i,] <- Ftest.glo[min.p,]
               Ftest.partial[[i]] <- Ftest.par[[min.p]]
               list.cata[[i]] <- bin.split[[1]][[min.p]]
               list.catb[[i]] <- bin.split[[2]][[min.p]]
          }
     }

     ### els chafillas
     if (mox=="pathmox")
     {
         otras <- order(Ftest.global[,4]) 
         F.otras <- Ftest.global[otras,]
         F.otras[,4] <- round(F.otras[,4],6)
         F.otras[,5] <- round(F.otras[,5],4)
         colnames(F.otras) <- c("f.stat","df.num","df.den","p.val","sa/sb")
         cata.otras <- NULL
         catb.otras <- NULL
         for (i in 1:length(otras))
         {
             cata.otras <- rbind(cata.otras, paste(as.vector(list.cata[[otras[i]]]),sep="",collapse="/"))
             catb.otras <- rbind(catb.otras, paste(as.vector(list.catb[[otras[i]]]),sep="",collapse="/"))
         }
         test.otras <- data.frame(variable=as.vector(colnames(EXEV)[otras]), F.otras, 
                                  categ.a=cata.otras, categ.b=catb.otras)
         test.otras <- test.otras[which(!is.na(F.otras[,1])),]
         test.otras[,5] <- format(test.otras[,5], scientific=FALSE)
         colnames(test.otras)[6] <- "sa/sb"
     }
     if (mox=="techmox")
     {
         pvals.ave <- unlist(lapply(Ftest.partial, function(x) exp(mean(log(x[,4]))) ))
         otras <- order(pvals.ave) 
         cata.otras <- NULL
         catb.otras <- NULL
         for (i in 1:length(otras))
         {
             cata.otras <- rbind(cata.otras, paste(as.vector(list.cata[[otras[i]]]),sep="",collapse="/"))
             catb.otras <- rbind(catb.otras, paste(as.vector(list.catb[[otras[i]]]),sep="",collapse="/"))
         }
         test.otras <- data.frame(variable=as.vector(colnames(EXEV)[otras]), pval.geomean=sort(pvals.ave),
                                  categ.a=cata.otras, categ.b=catb.otras)
         test.otras <- test.otras[which(!is.na(Ftest.global[,1])),]
         test.otras[,2] <- format(test.otras[,2], scientific=FALSE)
     }

     ### el chingon
     pvals.ave <- unlist(lapply(Ftest.partial, function(x) exp(mean(log(x[,4]))) ))
     if (mox=="pathmox")
         optim <- which(Ftest.global[,4]==min(Ftest.global[,4]))[1]
     if (mox=="techmox")
         optim <- which(pvals.ave==min(pvals.ave))[1]
     f.optim <- Ftest.global[optim,]
     cats.optim <- list.cata[[optim]]
     anticat.optim <- list.catb[[optim]]
     aux.a <- which(E.node[,optim] %in% cats.optim)
     aux.b <- which(E.node[,optim] %in% anticat.optim)
     part.a <- indelem[aux.a,2]
     part.b <- indelem[aux.b,2]
     list.elems <- list(part.a, part.b)

     ### resul
     test.global <- data.frame(exev=colnames(EXEV)[optim], f.stat=f.optim[1], df.num=f.optim[2],
                      df.den=f.optim[3], p.val=f.optim[4], s1_s2=round(f.optim[5],4))
     colnames(test.global)[6] <- "sa/sb"
     test.partial <- list(exev=colnames(EXEV)[optim], p.val=pvals.ave[optim], 
                      F.test=Ftest.partial[[optim]])
     categs <- list(as.vector(cats.optim), as.vector(anticat.optim))
     resul <- list(inner.global=test.global, inner.partial=test.partial, categs=categs, 
                   otras=test.otras, list.elems=list.elems)
     return(resul)
}

.fix.xexeloa <-
function(pls, DT, EXEV, type.exev, elemnod, nv, size, mox)
{
     # cuauhmaitl "nv" elems 
     N <- nrow(EXEV)
     if (nv==0) elems<-1:length(elemnod) else elems<-which(elemnod==nv)
     indelem <- cbind(1:length(elems), elems)# identificador
     E.node <- EXEV[elems,]# elems en cuauhmaitl actual      
     exevs <- apply(E.node, 2, function(x) nlevels(as.factor(x)))# vector fo exevs
     list.cata <- as.list(1:length(exevs))   # categs a.cuauhmaitl
     list.catb <- as.list(1:length(exevs))   # categs b.cuauhmaitl
     Ftest.global <- matrix(NA, length(exevs), 5)# + chingon global
     Ftest.partial <- as.list(1:length(exevs))# + chingon partial

     # parameters of plspm
     Y.lvs <- pls$scores# matrix of latent variables
     IDM <- pls$model$IDM# IDM matrix
     blocks <- pls$model$blocks
     modes <- pls$model$modes
     scheme <- pls$model$scheme
     scaled <- pls$model$scaled
     tol <- pls$model$tol
     iter <- pls$model$iter
     endo <- rowSums(IDM)
     endo[endo!=0] <- 1  
     n.endo <- sum(endo)

     # size limite
     if (size < 1)  size.limit=ceiling(size*N)  else  size.limit=size

     # for each exev
     for (i in 1:length(exevs))    
     {
          if (exevs[i] == 1)  # ya se uso
          {
               Ftest.global[i,4] <- 1   # sta bien chafa
               list.cata[[i]] <- NA# NAs a.cuauhmaitl
               list.catb[[i]] <- NA# NAs b.cuauhmaitl
               next   # next exev
          }  else   # se usa               
          {
               ### moxexeloa
               cat.exe <- unique(E.node[,i])  
               if (type.exev[i] == "ord") {
                    bin.split <- .ordinal.split(sort(cat.exe))  
               } else { # nom
                    if (exevs[i]==2) {
                        bin.split <- as.list(cat.exe)
                    } else { 
                        bin.split <- .nominal.split(cat.exe)   # function of nominal splits
                    }
               }
               ### sanchez-aluja moxexeloa
               Ftest.glo <- matrix(NA, length(bin.split[[1]]), 5)   # global
               Ftest.par <- as.list(1:length(bin.split[[1]]))      # partial
               for (aux in 1:length(bin.split[[1]]))   # p/c moxexeloa
               {
                     split.a <- which(E.node[,i] %in% bin.split[[1]][[aux]])# a.cuauhmaitl
                     split.b <- which(E.node[,i] %in% bin.split[[2]][[aux]])  # b.cuauhmaitl
                     n.a <- length(split.a)
                     n.b <- length(split.b)
                     if (n.a<=size.limit || n.b<=size.limit)# cuahuitl aturat
                     {
                          Ftest.glo[aux,4] <- 1   # sta bien chafa
                          Ftest.par[[aux]] <- matrix(1,1,5)# sta bien chafa
                          next  
                     } else# cuahuitl segueix
                     {
                          idm <- .sanchez.aluja(Y.lvs, IDM, split.a, split.b)
                          Ftest.glo[aux,] <- idm[[1]]   # global
                          Ftest.par[[aux]] <- idm[[2]]  # partial
                     }
               }
               ### el gallo dels exev
               if (mox=="pathmox")
                   min.p <- which(Ftest.glo[,4]==min(Ftest.glo[,4]))[1]  
               if (mox=="techmox") 
               {
                   pvals.ave <- unlist(lapply(Ftest.par, function(x) exp(mean(log(x[,4]))) ))
                   min.p <- which(pvals.ave==min(pvals.ave))[1]
               }
               Ftest.global[i,] <- Ftest.glo[min.p,]
               Ftest.partial[[i]] <- Ftest.par[[min.p]]
               list.cata[[i]] <- bin.split[[1]][[min.p]]
               list.catb[[i]] <- bin.split[[2]][[min.p]]
          }
     }
     ### els chafillas
     if (mox=="pathmox")
     {
         otras <- order(Ftest.global[,4]) 
         F.otras <- Ftest.global[otras,]
         F.otras[,4] <- round(F.otras[,4],6)
         F.otras[,5] <- round(F.otras[,5],4)
         colnames(F.otras) <- c("f.stat","df.num","df.den","p.val","sa/sb")
         cata.otras <- NULL
         catb.otras <- NULL
         for (i in 1:length(otras))
         {
             cata.otras <- rbind(cata.otras, paste(as.vector(list.cata[[otras[i]]]),sep="",collapse="/"))
             catb.otras <- rbind(catb.otras, paste(as.vector(list.catb[[otras[i]]]),sep="",collapse="/"))
         }
         test.otras <- data.frame(variable=as.vector(colnames(EXEV)[otras]), F.otras, 
                                  categ.a=cata.otras, categ.b=catb.otras)
         test.otras <- test.otras[which(!is.na(F.otras[,1])),]
         test.otras[,5] <- format(test.otras[,5], scientific=FALSE)
         colnames(test.otras)[6] <- "sa/sb"
         who.otras <- otras[which(!is.na(F.otras[,1]))]
     }
     if (mox=="techmox")
     {
         pvals.ave <- unlist(lapply(Ftest.partial, function(x) exp(mean(log(x[,4]))) ))
         otras <- order(pvals.ave) 
         pvals.otras <- pvals.ave[otras]
         cata.otras <- NULL
         catb.otras <- NULL
         for (i in 1:length(otras))
         {
             cata.otras <- rbind(cata.otras, paste(as.vector(list.cata[[otras[i]]]),sep="",collapse="/"))
             catb.otras <- rbind(catb.otras, paste(as.vector(list.catb[[otras[i]]]),sep="",collapse="/"))
         }
         test.otras <- data.frame(variable=as.vector(colnames(EXEV)[otras]), pval.geomean=sort(pvals.ave),
                                  categ.a=cata.otras, categ.b=catb.otras)
         test.otras <- test.otras[which(!is.na(Ftest.global[,1])),]
         test.otras[,2] <- format(test.otras[,2], scientific=FALSE)
         who.otras <- otras[which(pvals.otras!=1)]
     }

     ### el fijado         

     if (length(who.otras)>0)
     {
         candi.table <- cbind(Num=who.otras, test.otras)
         if (nv==0)
             cat("CANDIDATE SPLITS: ROOT NODE", "\n\n")
         if (nv>0)
             cat("CANDIDATE SPLITS: NODE", nv, "\n\n")
         print(candi.table)
         cat("\n")
         cat("Select a candidate split from column 'Num' and then press Enter:", "\n")
         fijo <- scan(what=integer(1), n=1)
         if (mode(fijo)!="numeric" || length(fijo)!=1 || fijo<1 || (fijo%%1)!=0) {
             cat("Invalid number of candidate split", "\n")
             stop() 
         }
         f.fijo <- Ftest.global[fijo,]
     } else {
         fijo <- which(Ftest.global[,4]==min(Ftest.global[,4]))[1]
         f.fijo <- Ftest.global[fijo,]
     }         
     cats.fijo <- list.cata[[fijo]]
     anticat.fijo <- list.catb[[fijo]]

     aux.a <- which(E.node[,fijo] %in% cats.fijo)
     aux.b <- which(E.node[,fijo] %in% anticat.fijo)
     part.a <- indelem[aux.a,2]
     part.b <- indelem[aux.b,2]
     list.elems <- list(part.a, part.b)

     ### resul
     if (mox=="pathmox")
     {
         test.fijo <- data.frame(exev=colnames(EXEV)[fijo], f.stat=f.fijo[1], df.num=f.fijo[2],
                            df.den=f.fijo[3], p.val=f.fijo[4], sa_sb=round(f.fijo[5],4))
         colnames(test.fijo)[6] <- "sa/sb"
     }
     if (mox=="techmox")
     {
         test.fijo <- list(exev=colnames(EXEV)[fijo], p.val=pvals.ave[fijo], 
                      F.test=Ftest.partial[[fijo]])
     }
     categs <- list(as.vector(cats.fijo), as.vector(anticat.fijo))
     resul <- list(inner.test=test.fijo, categs=categs, list.elems=list.elems, otras=test.otras)
     return(resul)
}

