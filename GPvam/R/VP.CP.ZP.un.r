VP.CP.ZP.un <-                                            
function (Z_mat, fixed_effects, control) 
{   
persistence<-control$persistence
    gr.eta <- function(eta, X, Y, Z, ybetas, R_inv, G.inv) {
        gr.y <- crossprod(Z, R_inv) %*% (Y - X %*% ybetas - Z %*% 
            eta)
        gr.p.eta <- -G.inv %*% eta
        -as.vector(gr.y + gr.p.eta)
    }
    H.eta <- function(G.inv, Z, R_inv) {
        h.eta <- G.inv
        h.y <- crossprod(Z, R_inv) %*% Z
        forceSymmetric(h.eta + h.y)
    }
    ltriangle <- function(x) {
        if (!is.null(dim(x)[2])) {
            resA <- as.vector(x[lower.tri(x, diag = TRUE)])
            resA
        }else {
            nx <- length(x)
            d <- 0.5 * (-1 + sqrt(1 + 8 * nx))
            resB <- .symDiagonal(d)
            resB[lower.tri(resB, diag = TRUE)] <- x
            if (nx > 1) {
                resB <- resB + t(resB) - diag(diag(resB))
            }
            as(resB, "sparseMatrix")
        }
    }
    reduce.G <- function(G, nyear, nteacher) {
        if (!is.null(dim(G)[2])) {
            temp_mat <- G
            index1 <- 0
            resA <- c(NULL)
            for (j in 1:nyear) {
                temp_mat_j <- as.numeric(temp_mat[(index1 + 1):(index1 + 1), (index1 + 1):(index1 + 1)])
                resA <- c(resA, temp_mat_j)
                index1 <- index1 + nteacher[j]
            }
            resA
        }
        else {
            resB <- Matrix(0, 0, 0)
          
            for (j in 1:nyear) {
                resB <- bdiag(resB, suppressMessages(G[j]*Diagonal(nteacher[j])))
               
            }
            rm(j)
            resB
        }
    }
    update.eta <- function(X, Y, Z, R_inv, ybetas, G, nyear, 
        cons.logLik, Ny, nstudent, n_eta) {
        G.chol <- chol(G)
        G.inv <- chol2inv(G.chol)
        H <- H.eta(G.inv, Z, R_inv)
        chol.H <- chol(H)
        var.eta <- as.matrix(chol2inv(chol.H))
        eta<- var.eta%*%as.vector(crossprod(Z, R_inv) %*% (Y - X %*% ybetas))
        log.p.eta <- -(n_eta/2) * log(2 * pi) - sum(log(diag(G.chol))) - 
            0.5 * crossprod(eta, G.inv) %*% eta
        log.p.y <- -(Ny/2) * log(2 * pi) + sum(log(diag(chol(R_inv)))) - 
            0.5 * crossprod(Y - X %*% ybetas - Z %*% eta, R_inv) %*% 
                (Y - X %*% ybetas - Z %*% eta)
        res <- var.eta
        attr(res, "likelihood") <- as.vector(cons.logLik + log.p.eta + 
            log.p.y - sum(log(diag(chol.H))))
        attr(res, "eta") <- eta
        res
    }
    Score <- function(thetas, eta = eta.hat, ybetas, X, Y, Z, 
        year.count, n_ybeta, nyear, n_eta, nstudent, nteacher, 
        Kg, cons.logLik, con = control, mis.list, pattern.parmlist2, 
        pattern.count, pattern.length, pattern.Rtemplate, pattern.diag, 
        pattern.key, Ny, pattern.sum=pattern.sum,persistence,P,alpha.diag,nalpha,alpha) {
        n_Rparm <- nyear * (nyear + 1)/2
        G <- thetas[(n_Rparm + 1):(n_Rparm + nyear)]
        G <- reduce.G(G = G, nyear = nyear, nteacher = nteacher)
        if(persistence=="VP"){
         alpha.parm<-thetas[(n_Rparm + nyear+1):length(thetas)]
         alpha[!((1:nalpha)%in%alpha.diag)]<-alpha.parm
         
                  Z<-Matrix(0,nrow(Z_mat),nteach_effects)
         
       for(i in 1:nalpha){
        comp<-which(tril(ltriangle(1:nalpha))==i,arr.ind=TRUE)
        Z<-Z+alpha[i]*P[[comp[1]]][[comp[2]]]
       }

       if(!huge.flag){Z.dense <- as.matrix(Z)}
           for (p in unique(Z_mat$pat)) {
        if(!huge.flag){
        Z.p[[p]] <- Z.dense[pat[[p]], , drop = FALSE]}else{
        Z.p[[p]] <- Z[pat[[p]], , drop = FALSE]
        }
    }
   # rm(Z.dense)
    }
        R_i <- ltriangle(as.vector(thetas[1:n_Rparm]))
        R_i.parm <- as.vector(thetas[1:n_Rparm])
        if (length(mis.list) > 0) {
            R <- symmpart(suppressMessages(kronecker(suppressMessages(Diagonal(nstudent)), 
                R_i)[-mis.list, -mis.list]))
            R_inv <- solve(R)
        }else {
            R <- symmpart(suppressMessages(kronecker(suppressMessages(Diagonal(nstudent)), 
                R_i)))
            R_inv <- symmpart(suppressMessages(kronecker(suppressMessages(Diagonal(nstudent)), 
                chol2inv(chol(R_i)))))
        }
        new.eta <- update.eta(X = X, Y = Y, Z = Z, 
            R_inv = R_inv, ybetas = ybetas, G = G, nyear = nyear, 
            cons.logLik = cons.logLik, Ny = Ny, nstudent = nstudent, 
            n_eta = n_eta)
        eta.hat <- attr(new.eta, "eta")
        var.eta.hat <- new.eta
        temp_mat <- var.eta.hat + tcrossprod(eta.hat, eta.hat)
                
                pattern.sum <- list()
for (p in unique(Z_mat$pat)) {
pattern.sum[[p]]<-R_mstep2(invsqrtW_=as.matrix(rep(1,Ny)),JYp_=as.matrix(Y.p[[p]]),loopsize_=pattern.count[[p]]/pattern.length[[p]],
  patternlength_=pattern.length[[p]],rownumber_=as.matrix(Y.p.rownumber[[p]]),ybetas_=as.matrix(ybetas),
  etahat_=as.matrix(eta.hat),tempmatR_=as.matrix(temp_mat),
  JXpi_=as.matrix(X.p[[p]]@i),JXpp_=as.matrix(X.p[[p]]@p),JXpx_=as.matrix(X.p[[p]]@x),JXpdim_=as.matrix(X.p[[p]]@Dim),
  JZpi_=as.matrix(Z.p[[p]]@i),JZpp_=as.matrix(Z.p[[p]]@p),JZpx_=as.matrix(Z.p[[p]]@x),JZpdim_=as.matrix(Z.p[[p]]@Dim))
}  
 
        score.R <- -pattern.f.score(R_i.parm, nyear, pattern.parmlist2, 
            pattern.count, pattern.length, pattern.Rtemplate, 
            pattern.diag, pattern.key, pattern.sum)
        temp_mat <- G
        gam_t <- list()
        sv_gam_t <- list()
        index1 <- 0
        for (j in 1:nyear) {
            gam_t[[j]] <- matrix(0, Kg[j], Kg[j])
            temp_mat_j <- temp_mat[(index1 + 1):(index1 + nteacher[j] * 
                Kg[j]), (index1 + 1):(index1 + nteacher[j] * 
                Kg[j])]
            gam_t[[j]] <- temp_mat_j[1:Kg[j], 1:Kg[j]]
            sv_gam_t[[j]] <- chol2inv(chol(gam_t[[j]]))
            index1 <- index1 + nteacher[j] * Kg[j]
        }
        rm(j)
        score_mat <- var.eta.hat + tcrossprod(eta.hat, eta.hat)
        gam_t_sc <- list()
        index1 <- 0
        score.G <- Matrix(0, 0, 0)
        for (j in 1:nyear) {
            gam_t_sc[[j]] <- matrix(0, Kg[j], Kg[j])
            score_mat_j <- score_mat[(index1 + 1):(index1 + nteacher[j] * 
                Kg[j]), (index1 + 1):(index1 + nteacher[j] * 
                Kg[j])]
            index2 <- c(1)
            for (k in 1:nteacher[j]) {
                gam_t_sc[[j]] <- gam_t_sc[[j]] + score_mat_j[(index2):(index2 + 
                  Kg[j] - 1), (index2):(index2 + Kg[j] - 1)]
                index2 <- index2 + Kg[j]
            }
            index1 <- index1 + nteacher[j] * Kg[j]
            der <- -0.5 * (nteacher[j] * sv_gam_t[[j]] - sv_gam_t[[j]] %*% 
                gam_t_sc[[j]] %*% sv_gam_t[[j]])
            if (is.numeric(drop(sv_gam_t[[j]]))) {
                score.eta.t <- der
            }
            else {
                score.eta.t <- 2 * der - diag(diag(der))
            }
            for (k in 1:nteacher[j]) {
                score.G <- bdiag(score.G, score.eta.t)
            }
            
        }
        alpha.parm<-alpha[!((1:nalpha)%in%alpha.diag)]
        score.a<-alpha.score(alpha.parm,alpha,temp_mat=score_mat,nalpha,alpha.diag,P,R_inv,eta.hat,ybetas,X)
        rm(j, k)
       if(persistence=="CP"|persistence=="ZP"){
       -c(score.R, reduce.G(G = score.G, nyear = nyear, nteacher = nteacher))}else if(persistence=="VP"){
       -c(score.R, reduce.G(G = score.G, nyear = nyear, nteacher = nteacher),score.a)
       }
       
    }
    update.ybeta <- function(X, Y, Z, R_inv, eta.hat) {
        A.ybeta <- crossprod(X, R_inv) %*% X
        B.ybeta <- crossprod(X, R_inv) %*% (Y - Z %*% eta.hat)
        as.vector(solve(A.ybeta, B.ybeta))
    }
    bin2dec <- function(s) sum(s * 2^(rev(seq_along(s)) - 1))
    dec2bin <- function(s) {
        L <- length(s)
        maxs <- max(s)
        digits <- floor(logb(maxs, base = 2)) + 1
        res <- array(NA, dim = c(L, digits))
        for (i in 1:digits) {
            res[, digits - i + 1] <- (s%%2)
            s <- (s%/%2)
        }
        if (L == 1) 
            res[1, ]
        else res
    }
    alpha.score<-function(alpha.parm,alpha,temp_mat,nalpha,alpha.diag,P,R_inv,eta.hat,ybetas,X){
        score.a<-c(NULL)
        alpha.s<-alpha
        alpha.s[!((1:nalpha)%in%alpha.diag)]<-alpha.parm
        alpha.set<-(1:nalpha)[!((1:nalpha)%in%alpha.diag)]
        Z.a<-Matrix(0,nrow(Z_mat),nteach_effects)
       for(i in 1:nalpha){
        comp<-which(tril(ltriangle(1:nalpha))==i,arr.ind=TRUE)
        Z.a<-Z.a+alpha.s[i]*P[[comp[1]]][[comp[2]]]    
       }  
        for(i in alpha.set){
        comp<-which(tril(ltriangle(1:nalpha))==i,arr.ind=TRUE)
        score.a<-c(score.a,as.numeric((t(Y)-t(ybetas)%*%t(X))%*%R_inv%*%P[[comp[1]]][[comp[2]]]%*%eta.hat-sum(diag(crossprod(Z.a,R_inv%*%P[[comp[1]]][[comp[2]]]%*%temp_mat)))))
        }
        score.a
        }
    R_mstep2 <- function(invsqrtW_,JYp_,loopsize_, patternlength_,rownumber_,ybetas_,etahat_,tempmatR_,JXpi_,JXpp_,JXpx_,JXpdim_,JZpi_,JZpp_,JZpx_,JZpdim_){
.Call( "R_mstep_cpp",invsqrtW_,JYp_,loopsize_, patternlength_,rownumber_,ybetas_,etahat_,tempmatR_,JXpi_,JXpp_,JXpx_,JXpdim_,JZpi_,JZpp_,JZpx_,JZpdim_)
}

pattern.f.score <- function(R_i.parm, nyear, pattern.parmlist2, pattern.count, pattern.length, pattern.Rtemplate, pattern.diag, pattern.key, pattern.sum) {
     R_i <- as.matrix(ltriangle(as.vector(R_i.parm)))
    pattern.score <- numeric(nyear/2 * (nyear + 1) )
    for (p in nonempty.patterns) {
    pattern.y <- solve(pattern.f.R(R_i, p, nyear, pattern.key))
    YSY<-pattern.y %*% pattern.sum[[p]] %*% pattern.y
    PCL<-pattern.countoverlength[[p]]
        for (r.parm in 1:(nyear/2 * (nyear + 1))) {
            if(is.null(pat.coord <- pat.coord.guide[[r.parm]][[p]])) next
            pattern.yc <- pattern.y[pat.coord]
            pattern.score[r.parm] <- pattern.score[r.parm] - ( PCL* pattern.yc - YSY[pat.coord])

        }

    }
    pattern.score[1:(nyear/2 * (nyear + 1) ) %in% pattern.diag] <- 0.5 * pattern.score[1:(nyear/2 * (nyear + 1) ) %in% pattern.diag]
    -pattern.score
}
  
      pattern.f.R <- function(R, p, nyear, pattern.key) {
        R[pattern.key[p, ] * (1:nyear), pattern.key[p, ] * (1:nyear), 
            drop = FALSE]
    }
    Z_mat$year <- as.numeric(Z_mat$year)
    nyear <- length(unique(Z_mat$year))
    Z_mat$mis <- rep(0, dim(Z_mat)[1])
    student_list <- unique(Z_mat$student)
    Z_mat[is.na(Z_mat$y), ]$mis <- rep(1, dim(Z_mat[is.na(Z_mat$y), 
        ])[1])
    for (g in 1:nyear) {
        mis_stu <- student_list[!(student_list %in% unique(Z_mat[Z_mat$year == 
            g, ]$student))]
        le <- length(mis_stu)
        if (le > 0) {
            temp.exp <- Z_mat[1:le, ]
            temp.exp$year <- rep(g, le)
            temp.exp$mis <- rep(1, le)
            temp.exp$student <- mis_stu
            temp.exp$teacher <- rep(NA, le)
            temp.exp$y <- rep(NA, le)
            Z_mat <- rBind(Z_mat, temp.exp)
        }
    }
    rm(g, le)
    Z_mat.full <- Z_mat
    Z_mat <- Z_mat[!is.na(Z_mat$y), ]
    Z_mat.full <- Z_mat.full[which((Z_mat.full$student %in% Z_mat$student)), 
        ]
    Ny <- sum(Z_mat$mis == 0)
    nstudent <- length(unique(Z_mat$student))
    year.count <- numeric(nyear)
    for (j in 1:nyear) {
        year.count[j] <- sum(Z_mat[Z_mat$year == j, ]$mis == 
            0)
    }
    rm(j)
    RE_s_start_pos <- 1
    Kg <- rep(1,nyear)
    Z_mat <- Z_mat[order(Z_mat$year, Z_mat$teacher), ]
    Z_mat.full <- Z_mat.full[order(Z_mat.full$year, Z_mat.full$teacher), 
        ]
    na_list <- grep("^NA", Z_mat$teacher)
    if (length(na_list) > 0) {
        teachyearcomb <- unique(cBind(Z_mat[-na_list, ]$year, Z_mat[-na_list, ]$teacher))
    }else {
        teachyearcomb <- unique(cBind(Z_mat$year, Z_mat$teacher))
    }
    

    Z_mat <- Z_mat[order(Z_mat$student, Z_mat$year, Z_mat$teacher), 
        , drop = FALSE]
    Z_mat.full <- Z_mat.full[order(Z_mat.full$student, Z_mat.full$year, 
        Z_mat.full$teacher), ]
    nteach_effects <- dim(teachyearcomb)[1]
      teacheffZ_mat <- Matrix(0, nrow = nrow(Z_mat), ncol = nteach_effects)
    t_effects <- rep(NA, nteach_effects)
    indx <- 1
    eblup.tracker <- matrix(0, 0, 3)
    dP<-list()
      for(i in 1:nyear){
       dP[[i]]<-list()
       for(j in 1:i)
       dP[[i]][[j]]<-matrix(0,0,2)
       }
       nalpha<-nyear/2*(nyear+1)
              if(persistence=="ZP"|persistence=="VP"){
       alpha<-ltriangle(diag(nyear))
       }else if(persistence=="CP"){
       alpha<-rep(1,nalpha)
       }
       alpha_key<-tril(ltriangle(1:nalpha))
       alpha.diag<-diag(alpha_key)
       alpha.parm<-alpha[!((1:nalpha)%in%alpha.diag)]
       

      
      for (k in 1:nrow(teachyearcomb)) {
          student_subset <- Z_mat.full$student[Z_mat.full$year == 
              teachyearcomb[k, 1] & Z_mat.full$teacher == teachyearcomb[k, 
              2]]
                            t_effects[k] <- paste(teachyearcomb[k, 1], "_", 
                  teachyearcomb[k, 2], sep = "")
          eblup.tracker <- rBind(eblup.tracker, c(teachyearcomb[k, ], teachyearcomb[k,1]))
          for (yr in as.numeric(teachyearcomb[k, 1]):nyear) {
              if (sum(is.element(Z_mat$student, student_subset) & 
                  Z_mat$year == yr) != 0) {
                  q1<-(1:nrow(Z_mat))[is.element(Z_mat$student, student_subset) & Z_mat$year == yr & !is.na(Z_mat$y)]  
                  q2<-rep(k,length(q1))
                  dP[[yr]][[as.numeric(teachyearcomb[k, 1])]]<-rbind(dP[[yr]][[as.numeric(teachyearcomb[k, 1])]],cbind(q1,q2))
              }     

          }
      }
            P<-list()
      for(i in 1:nyear){
       P[[i]]<-list()
       for(j in 1:i)
       P[[i]][[j]]<-as(sparseMatrix(i=dP[[i]][[j]][,1],j=dP[[i]][[j]][,2],dims=c(nrow(Z_mat),nteach_effects)),"dgCMatrix")
       }
       Z<-Matrix(0,nrow(Z_mat),nteach_effects)
       for(i in 1:nalpha){
        comp<-which(tril(ltriangle(1:nalpha))==i,arr.ind=TRUE)
        Z<-Z+alpha[i]*P[[comp[1]]][[comp[2]]]
       }
 
    mis.list <- which(Z_mat.full$mis == 1)
    nteacher <- as.vector(tapply(teachyearcomb[, 2], teachyearcomb[, 
        1], length))
    colnames(Z) <- t_effects
    X_mat <- sparse.model.matrix(fixed_effects, Z_mat, drop.unused.levels = TRUE)
    X_mat <- X_mat[, !(colSums(abs(X_mat)) == 0), drop = FALSE]
    if (rankMatrix(X_mat,method = 'qrLINPACK')[1] != dim(X_mat)[2]) {
    stop("WARNING: Fixed-effects design matrix not full-rank")

    }

    n_eta <- nteach_effects
    n_ybeta <- dim(X_mat)[2]
    Y <- as.vector(Z_mat$y)
    X <- Matrix(X_mat)
    huge.flag<-TRUE
     if(!huge.flag){
    Z.dense <- as.matrix(Z)
    X.dense <- as.matrix(X)
    }
    Z_mat.full$r <- 1 - Z_mat.full$mis
    pattern.student <- matrix(Z_mat.full$r, nstudent, nyear, 
        byrow = TRUE)
    Z_mat.full$pat <- rep(apply(pattern.student, 1, bin2dec), 
        each = nyear)
    Z_mat$pat <- Z_mat.full[Z_mat.full$r == 1, ]$pat
    pat <- list()
    pattern.count <- list()
    pattern.length <- list()
    X.p <- list()
    Y.p <- list()
    Z.p <- list()
     Y.p.rownumber <- list()
    pattern.countoverlength<-list()
    rownumber <- 1:Ny
    pattern.key <- dec2bin(1:(2^nyear - 1))
    X<-as(X,"sparseMatrix")
    for (p in unique(Z_mat$pat)) {
        pat[[p]] <- which(Z_mat$pat == p)
        if(!huge.flag){
        X.p[[p]] <- X.dense[pat[[p]], , drop = FALSE]     
        Z.p[[p]] <- Z.dense[pat[[p]], , drop = FALSE]
        }else{
        X.p[[p]] <- X[pat[[p]], , drop = FALSE]     
        Z.p[[p]] <- Z[pat[[p]], , drop = FALSE]
        }
        Y.p.rownumber[[p]] <- rownumber[pat[[p]]]
        Y.p[[p]] <- Y[pat[[p]]]
        pattern.count[[p]] <- length(Y.p[[p]])
        pattern.length[[p]] <- sum(pattern.key[p, ])
        pattern.countoverlength[[p]]<-pattern.count[[p]]/pattern.length[[p]]
    }
   # rm(Z.dense, X.dense)
    pattern.yguide <- list()
    for (g in 1:nyear) {
        pattern.yguide[[g]] <- which(pattern.key[, g] == 1)
    }
    pattern.Rtemplate <- ltriangle(1:(nyear/2 * (nyear + 1)))
    pattern.diag <- diag(pattern.Rtemplate)
    pattern.Rtemplate <- ltriangle(1:(nyear/2 * (nyear + 1)))
    pattern.parmlist1 <- list()
    pattern.parmlist2 <- list()
    for (p in unique(Z_mat$pat)) {
        pattern.parmlist1[[p]] <- sort(unique(as.vector(pattern.f.R(pattern.Rtemplate, 
            p, nyear, pattern.key))))
    }
    for (r.parm in 1:(nyear/2 * (nyear + 1))) {
        pattern.parmlist2[[r.parm]] <- which(sapply(pattern.parmlist1, 
            f <- function(x) r.parm %in% x))
    }
    eta.hat <- numeric(n_eta)
    var.eta.hat <- Matrix(0, n_eta, n_eta)
    R.temp.comp <- numeric(nyear)
    for (g in 1:nyear) {
        R.temp.comp[g] <- var(Z_mat[Z_mat$year == g, ]$y)/4
    }
    R_i <- diag(R.temp.comp)
    if (length(mis.list) > 0) {
        R <- symmpart(suppressMessages(kronecker(suppressMessages(Diagonal(nstudent)), 
            R_i)[-mis.list, -mis.list]))
        R_inv <- solve(R)
    }else {
        R <- symmpart(suppressMessages(kronecker(suppressMessages(Diagonal(nstudent)), 
            R_i)))
        R_inv <- symmpart(suppressMessages(kronecker(suppressMessages(Diagonal(nstudent)), 
            chol2inv(chol(R_i)))))
    }
    
        pat.coord.guide<-list()
for (r.parm in 1:(nyear/2 * (nyear + 1) )) {
pat.coord.guide[[r.parm]]<-list()
        for (p in pattern.parmlist2[[r.parm]]) {
            pat.coord.guide[[r.parm]][[p]] <- which(tril(pattern.f.R(pattern.Rtemplate, p, nyear, pattern.key)) == r.parm)

         }
}
    nonempty.patterns<-NULL
for(p in 1:length(pat.coord.guide[[1]])){
if(is.null(retrieve.parm<-pattern.parmlist1[[p]])) next
nonempty.patterns<-c(nonempty.patterns,p)
}

    
    ybetas <- update.ybeta(X, Y, Z, R_inv, eta.hat)
    names(ybetas) <- colnames(X_mat)
    G <- 100*.symDiagonal(n_eta)
    cons.logLik <- 0.5 * n_eta * log(2 * pi)
    iter <- control$max.iter.EM
    Y.mat <- Matrix(0, iter, n_ybeta)
    G.mat <- Matrix(0, iter, length(reduce.G(G = G, nyear = nyear, nteacher = nteacher)))
    R.mat <- Matrix(0, iter, nyear * (nyear + 1)/2)
    lgLik <- numeric(iter)
    conv <- FALSE
   if (control$verbose) cat("Beginning EM algorithm\n")
    flush.console()
    for (it in 1:iter) {
        ptm <- proc.time()[3]
        suppressWarnings(rm(var.eta.hat,temp_mat))
         mresid <- as.numeric(Y - X %*% ybetas)
    cresid <- as.numeric(mresid - Z %*% eta.hat)
    yhat <- as.numeric(X %*% ybetas + Z %*% eta.hat)
    yhat.m <- as.numeric(X %*% ybetas)
        new.eta <- update.eta(X = X, Y = Y, Z = Z, 
            R_inv = R_inv, ybetas = ybetas, G = G, nyear = nyear, 
            cons.logLik = cons.logLik, Ny = Ny, nstudent = nstudent, 
            n_eta = n_eta)
        eta.hat <- attr(new.eta, "eta")
        var.eta.hat <- new.eta
        temp_mat <- var.eta.hat + tcrossprod(eta.hat, eta.hat)
        lgLik[it] <- attr(new.eta, "likelihood")
        rm(new.eta)
        thets1 <- c(Y.mat[it - 1, ], R.mat[it - 1, ], G.mat[it - 
            1, ])
        thets2 <- c(Y.mat[it, ], R.mat[it, ], G.mat[it, ])
        if (it > 5) {
            check.lik <- abs(lgLik[it] - lgLik[it - 1])/abs(lgLik[it] + 
                control$tol1) < control$tol1
            if (check.lik) {
                conv <- TRUE
                if (control$verbose) {
                  cat("\n\n Algorithm converged.\n")
                  cat("\n\niter:", it, "\n")
                  cat("log-likelihood:", sprintf("%.7f", lgLik[it]), 
                    "\n")
                  cat("change in loglik:", sprintf("%.7f", lgLik[it] - 
                    lgLik[it - 1]), "\n")
                  cat("fixed effects:", round(ybetas, 4), "\n")
                  cat("R_i:\n")
                  print(round(as.matrix(R_i), 4))
                  cat("\n")
                  print(round(cov2cor(as.matrix(R_i)), 4))
                  cat("\n")
                  for (j in 1:nyear) {
                    cat("\ngamma_teach_year", j, "\n")
                    print(round(reduce.G(G,nyear,nteacher)[j], 4))
                    cat("\n")
                    flush.console()
                  }
                  rm(j)
                }
                break
            }
            if (it==iter) {
                conv <- TRUE
                if (control$verbose) {
                  cat("\n\n Algorithm converged.\n")
                              cat("\n\niter:", it, "\n")
            cat("log-likelihood:", sprintf("%.7f", lgLik[it]), 
                "\n")
            cat("change in loglik:", sprintf("%.7f", lgLik[it] - 
                lgLik[it - 1]), "\n")
            cat("fixed effects:", round(ybetas, 4), "\n")
            cat("R_i:\n")
            print(round(as.matrix(R_i), 4))
            cat("\n")
            print(round(cov2cor(as.matrix(R_i)), 4))
            cat("\n")
            cat("G:\n")
            print(round(reduce.G(G,nyear,nteacher), 4))
             cat("\n")
            cat("alphas:\n")
            print(round(alpha, 4))
            rm(j)
                }
                break
            }
        }
        if ((control$verbose) & (it > 1)) {
            cat("\n\niter:", it, "\n")
            cat("log-likelihood:", sprintf("%.7f", lgLik[it]), 
                "\n")
            cat("change in loglik:", sprintf("%.7f", lgLik[it] - 
                lgLik[it - 1]), "\n")
            cat("fixed effects:", round(ybetas, 4), "\n")
            cat("R_i:\n")
            print(round(as.matrix(R_i), 4))
            cat("\n")
            print(round(cov2cor(as.matrix(R_i)), 4))
            cat("\n")
            cat("G:\n")
            print(round(reduce.G(G,nyear,nteacher), 4))
             cat("\n")
            cat("alphas:\n")
            print(round(alpha, 4))
            rm(j)
        }  
        
        
        
         
         
         if(persistence=="VP"){
         
alpha.parm<-alpha[!((1:nalpha)%in%alpha.diag)]
alpha.cc <- 1
        hes.count <- 1
       alpha.parm.old <- alpha.parm
        s.prev <- numeric(length(alpha.parm))

 #one step is sufficient because score function
 #is linear in alpha
            s <- alpha.score(alpha.parm,alpha=alpha,temp_mat=temp_mat,nalpha=nalpha,alpha.diag=alpha.diag,P=P,R_inv=R_inv,eta.hat=eta.hat,ybetas=ybetas,X=X)
            j <- jacobian(alpha.score, alpha.parm, method = "simple",alpha=alpha,temp_mat=temp_mat,nalpha=nalpha,alpha.diag=alpha.diag,P=P,R_inv=R_inv,eta.hat=eta.hat,ybetas=ybetas,X=X)
                hesprod <- solve(j, s)
               alpha.cc <- s %*% s
            alpha.parm <- alpha.parm - hesprod
            hes.count <- hes.count + 1
            s.prev <- s

alphan<-alpha
alphan[!((1:nalpha)%in%alpha.diag)]<-alpha.parm  
         
         
         
         }  
        #end update alphas
        d.temp<-diag(temp_mat)
        indx<-0
        Gn<-c(NULL)
        for(i in 1:nyear){
          Gn<-c(Gn,mean(d.temp[(indx+1):(indx+nteacher[i])]))
          indx<-indx+nteacher[i]                 
        }
     rm(indx,i,d.temp)


        ybetasn <- numeric(n_ybeta)
        ybetasn <- update.ybeta(X, Y, Z, R_inv, eta.hat)
        
                pattern.sum <- list()
for (p in unique(Z_mat$pat)) {
pattern.sum[[p]]<-R_mstep2(invsqrtW_=as.matrix(rep(1,Ny)),JYp_=as.matrix(Y.p[[p]]),loopsize_=pattern.count[[p]]/pattern.length[[p]],
  patternlength_=pattern.length[[p]],rownumber_=as.matrix(Y.p.rownumber[[p]]),ybetas_=as.matrix(ybetas),
  etahat_=as.matrix(eta.hat),tempmatR_=as.matrix(temp_mat),
  JXpi_=as.matrix(X.p[[p]]@i),JXpp_=as.matrix(X.p[[p]]@p),JXpx_=as.matrix(X.p[[p]]@x),JXpdim_=as.matrix(X.p[[p]]@Dim),
  JZpi_=as.matrix(Z.p[[p]]@i),JZpp_=as.matrix(Z.p[[p]]@p),JZpx_=as.matrix(Z.p[[p]]@x),JZpdim_=as.matrix(Z.p[[p]]@Dim))
}  
 
        R_i.parm <- ltriangle(as.matrix(R_i))
        R.cc <- 1
        hes.count <- 1
        R_i.parm.old <- R_i.parm
        s.prev <- numeric(length(R_i.parm))
        while (R.cc > 1e-04) {
            s <- pattern.f.score(R_i.parm = R_i.parm, nyear = nyear, 
                pattern.parmlist2 = pattern.parmlist2, pattern.count = pattern.count, 
                pattern.length = pattern.length, pattern.Rtemplate = pattern.Rtemplate, 
                pattern.diag = pattern.diag, pattern.key = pattern.key, 
                pattern.sum = pattern.sum)
            j <- jacobian(pattern.f.score, c(R_i.parm), method = "simple", 
                nyear = nyear, pattern.parmlist2 = pattern.parmlist2, 
                pattern.count = pattern.count, pattern.length = pattern.length, 
                pattern.Rtemplate = pattern.Rtemplate, pattern.diag = pattern.diag, 
                pattern.key = pattern.key, pattern.sum = pattern.sum)
if(it==1& (hes.count < 10)){     
           hesprod <- solve(j + max(c(diag(j), 5)) * diag(length(R_i.parm)), 
                  s)
                R.cc <- s %*% s
if(hes.count==9) R.cc<-0

          }else if ((it <=4) & (hes.count < 30)) {
                hesprod <- solve(j + max(c(diag(j), 5)*((1-hes.count/31)^2)) * diag(length(R_i.parm)), 
                  s)
                R.cc <- s %*% s
if(hes.count==29) R.cc<-0
            } else {
                hesprod <- solve(j, s)
                R.cc <- s %*% s
            }
            R_i.parm <- R_i.parm - hesprod
            hes.count <- hes.count + 1
            s.prev <- s
        }
        R_i <- ltriangle(R_i.parm)
        rm(R_i.parm)
        dimnames(R_i) <- list(NULL, NULL)
        if (length(mis.list) > 0) {
            R <- symmpart(suppressMessages(kronecker(Diagonal(nstudent), 
                R_i))[-mis.list, -mis.list])
            R_inv <- suppressMessages(solve(R))
        }else {
            R <- symmpart(suppressMessages(kronecker(Diagonal(nstudent), 
                R_i)))
            R_inv <- symmpart(suppressMessages(kronecker(Diagonal(nstudent), 
                solve(R_i))))
        }
        
        ybetas <- ybetasn
        G <- reduce.G(Gn,nyear,nteacher)
        
        if(persistence=="VP"){
        
                 alpha<-alphan
         Z<-Matrix(0,nrow(Z_mat),nteach_effects)
         
       for(i in 1:nalpha){
        comp<-which(tril(ltriangle(1:nalpha))==i,arr.ind=TRUE)
        Z<-Z+alpha[i]*P[[comp[1]]][[comp[2]]]
       }

       if(!huge.flag){Z.dense <- as.matrix(Z)}
           for (p in unique(Z_mat$pat)) {
        if(!huge.flag){
        Z.p[[p]] <- Z.dense[pat[[p]], , drop = FALSE]}else{
        Z.p[[p]] <- Z[pat[[p]], , drop = FALSE]
        }
    }
   # rm(Z.dense)
        
        }
        
    if (control$verbose)    cat("Iteration Time: ", proc.time()[3] - ptm, " seconds\n")
        flush.console()
    }
    names(ybetas) <- colnames(X_mat)
           if(persistence=="CP"|persistence=="ZP"){
        thetas <- c(ltriangle(as.matrix(R_i)), reduce.G(G = G, nyear = nyear, nteacher = nteacher))}else if(persistence=="VP"){
        thetas <- c(ltriangle(as.matrix(R_i)), reduce.G(G = G, nyear = nyear, nteacher = nteacher),alpha[!((1:nalpha)%in%alpha.diag)])
        }

    lgLik.hist <- lgLik
    lgLik <- lgLik[it]
    Hessian <- NA
     std_errors <- c(rep(NA, length(thetas)))
    if (control$hessian == TRUE) {
     if (control$verbose)   cat("Calculating Hessian of the variance components...")
        flush.console()
        
                pattern.sum <- list()
for (p in unique(Z_mat$pat)) {
pattern.sum[[p]]<-R_mstep2(invsqrtW_=as.matrix(rep(1,Ny)),JYp_=as.matrix(Y.p[[p]]),loopsize_=pattern.count[[p]]/pattern.length[[p]],
  patternlength_=pattern.length[[p]],rownumber_=as.matrix(Y.p.rownumber[[p]]),ybetas_=as.matrix(ybetas),
  etahat_=as.matrix(eta.hat),tempmatR_=as.matrix(temp_mat),
  JXpi_=as.matrix(X.p[[p]]@i),JXpp_=as.matrix(X.p[[p]]@p),JXpx_=as.matrix(X.p[[p]]@x),JXpdim_=as.matrix(X.p[[p]]@Dim),
  JZpi_=as.matrix(Z.p[[p]]@i),JZpp_=as.matrix(Z.p[[p]]@p),JZpx_=as.matrix(Z.p[[p]]@x),JZpdim_=as.matrix(Z.p[[p]]@Dim))
}  
 
        if (control$hes.method == "richardson") {
            Hessian <- ltriangle(ltriangle(jacobian(Score, thetas, 
                eta = eta.hat, ybetas = ybetas, X = X, Y = Y, 
                Z = Z, pattern.sum = pattern.sum, con = control, 
                year.count = year.count, n_ybeta = n_ybeta, nyear = nyear, 
                n_eta = n_eta, nstudent = nstudent, nteacher = nteacher, 
                Kg = Kg, cons.logLik = cons.logLik, mis.list = mis.list, 
                pattern.parmlist2 = pattern.parmlist2, pattern.count = pattern.count, 
                pattern.length = pattern.length, pattern.Rtemplate = pattern.Rtemplate, 
                pattern.diag = pattern.diag, pattern.key = pattern.key, 
                Ny = Ny,persistence=persistence,P=P,alpha.diag=alpha.diag,nalpha=nalpha,alpha=alpha)))
        }
        else {
            Hessian <- ltriangle(ltriangle(jacobian(Score, thetas, 
                method = "simple", eta = eta.hat, ybetas = ybetas, 
                X = X, Y = Y, Z = Z, pattern.sum = pattern.sum, 
                con = control, year.count = year.count, n_ybeta = n_ybeta, 
                nyear = nyear, n_eta = n_eta, nstudent = nstudent, 
                nteacher = nteacher, Kg = Kg, cons.logLik = cons.logLik, 
                mis.list = mis.list, pattern.parmlist2 = pattern.parmlist2, 
                pattern.count = pattern.count, pattern.length = pattern.length, 
                pattern.Rtemplate = pattern.Rtemplate, pattern.diag = pattern.diag, 
                pattern.key = pattern.key, Ny = Ny,persistence=persistence,P=P,alpha.diag=alpha.diag,nalpha=nalpha,alpha=alpha)))
        }
        std_errors <- try(c(sqrt(diag(solve(Hessian)))), silent = TRUE)
        hes.warn <- FALSE
        Hessian <- round(Hessian, 4)
        if (any(eigen(Hessian)$values <= 0)) {
       if (control$verbose)     cat("Warning: Hessian not PD", "\n")
            std_errors <- c(rep(NA, length(thetas)))
            hes.warn <- TRUE
        }
    }
    c.temp <- crossprod(X, R_inv) %*% Z
    c.1 <- rBind(crossprod(X, R_inv) %*% X, t(c.temp))
    G.inv <- chol2inv(chol(G))
    c.2 <- rBind(c.temp, H.eta(G.inv, Z, R_inv))
    C_inv <- cBind(c.1, c.2)
    C <- solve(C_inv)
    eblup_stderror <- sqrt(diag(C)[-c(1:n_ybeta)])
    ybetas_stderror <- sqrt(diag(C)[1:n_ybeta])
    rm(C,C_inv,c.2,c.1,c.temp)
    eblup <- cBind(eta.hat, eblup_stderror)
    eblup <- round(eblup, 4)
    eblup <- as.data.frame(eblup)
    eblup.tracker <- as.data.frame(eblup.tracker)
    eblup <- as.data.frame(cbind(eblup.tracker, eblup))
    colnames(eblup) <- c("teacher_year", "teacher", "effect_year", 
        "EBLUP", "std_error")
    eblup$teacher <- as.character(eblup$teacher)
    t_lab <- as.vector(NULL)
    r_lab <- as.vector(NULL)
    for (j in 1:nyear) {
        ne <- (Kg[j] * (Kg[j] + 1))/2
        y <- c(NULL)
        x <- c(NULL)
        for (k in 1:Kg[j]) {
            x <- c(x, k:Kg[j])
            y <- c(y, rep(k, (Kg[j] - k + 1)))
        }
        t_lab <- c(t_lab, paste("teacher effect from year", rep(j, 
            ne), sep = ""))
    }
    ne <- nyear * (nyear + 1)/2
    y <- c(NULL)
    x <- c(NULL)
    for (k in 1:nyear) {
        x <- c(x, k:nyear)
        y <- c(y, rep(k, (nyear - k + 1)))
    }
    r_lab <- paste("error covariance", ":[", x, ",", y, "]", 
        sep = "")
    rm(j, ne)
    alpha.label<-c(NULL)
    for(i in 1:(nyear-1)){
    for(j in (i+1):(nyear)){
    alpha.label<-c(alpha.label,paste("alpha_",j,i,sep=""))
    }
    }
           if(persistence=="CP"|persistence=="ZP"){
    effect_la <- c(names(ybetas), r_lab, t_lab)}else if(persistence=="VP"){
    effect_la <- c(names(ybetas), r_lab, t_lab, alpha.label)
    }
    
    if (control$hessian == TRUE) {
        parameters <- round(cBind(c(ybetas, thetas), c(ybetas_stderror, 
            std_errors)), 4)
        colnames(parameters) <- c("Estimate", "Standard Error")
        rownames(parameters) <- as.character(effect_la)
    }
    if (control$hessian == FALSE) {
        parameters <- round(cBind(c(ybetas, thetas), c(ybetas_stderror, 
            rep(NA, length(thetas)))), 4)
        colnames(parameters) <- c("Estimate", "Standard Error")
        rownames(parameters) <- as.character(effect_la)
    }
    R_i <- round(R_i, 4)
 if (control$verbose)   cat("done.\n")
    mresid <- as.numeric(Y - X %*% ybetas)
    cresid <- as.numeric(mresid - Z %*% eta.hat)
    yhat <- as.numeric(X %*% ybetas + Z %*% eta.hat)
    yhat.m <- as.numeric(X %*% ybetas)
    Hessian <- round(Hessian, 5)
    R_i <- round(R_i, 4)
    gam_t<-list()
    for (i in 1:nyear) {
        gam_t[[i]] <- round(as.matrix(reduce.G(G,nyear=nyear,nteacher=nteacher)[i]), 4)
        colnames(gam_t[[i]]) <- paste("year", i, sep = "")
        rownames(gam_t[[i]]) <- colnames(gam_t[[i]])
    }
    rchol <- try(chol(R_inv))
    yhat.s <- try(as.vector(rchol %*% (yhat)))
    sresid <- try(as.vector(rchol %*% Y - yhat.s))
    teach.cov <- lapply(gam_t, function(x) round(x, 4))
   L<- list(loglik = lgLik, teach.effects = eblup, parameters = parameters, 
        Hessian = Hessian, R_i = as.matrix(R_i), teach.cov = gam_t, 
        mresid = mresid, cresid = cresid, y = Y, yhat = yhat, 
        stu.cov = NA, num.obs = Ny, num.student = nstudent, num.year = nyear, 
        num.teach = nteacher, yhat.m = yhat.m, sresid = sresid, 
        yhat.s = yhat.s,iter=it,persistence=control$persistence)
}
