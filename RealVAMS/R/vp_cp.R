vp_cp <-
function(Z_mat,B.mat,control){
fam.binom<-control$outcome.family
persistence <- control$persistence
#gr.eta <- function(eta, X, Y, Z, ybetas, R_inv, G.inv) {
#    gr.y <- crossprod(Z, R_inv) %*% (Y - X %*% ybetas - Z %*% eta)
#    gr.p.eta <- -G.inv %*% eta
#    -as.vector(gr.y + gr.p.eta)
#}
  H.eta <- function(G.inv, Z, R.full.inv) {
      h.eta <- G.inv
      h.y <- crossprod(Z, R.full.inv) %*% Z
      as(symmpart(h.eta + h.y),"sparseMatrix")
  }
ltriangle <- function(x) {
    if (!is.null(dim(x)[2])) {
        resA <- as.vector(suppressMessages(x[lower.tri(x, diag = TRUE)]))
        resA
    } else {
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
reduce.G <- function(G, nyear.score, nteacher) {
    if (!is.null(dim(G)[2])) {
        temp_mat <- G
        index1 <- 0
        resA <- c(NULL)
        for (j in 1:nyear.score) {
            temp_mat_j <- ltriangle(as.matrix(suppressMessages(temp_mat[(index1 + 1):(index1 + 2), (index1 + 1):(index1 + 2)])))
            resA <- c(resA, temp_mat_j)
            index1 <- index1 + 2 * nteacher[j]
        }
        if(control$school.effects){
          temp_mat_j <- ltriangle(as.matrix(suppressMessages(temp_mat[(index1 + 1):(index1 + 2), (index1 + 1):(index1 + 2)])))
            resA <- c(resA, temp_mat_j)
        }
        resA
    } else {
        resB <- Matrix(0, 0, 0)
        for (j in 1:nyear.score) {
            resB <- bdiag(resB, suppressMessages(kronecker(Diagonal(nteacher[j]), ltriangle(G[(3 * j - 2):(3 * j)]))))
        }
        if(control$school.effects){
            resB<-bdiag(resB,suppressMessages(kronecker(Diagonal(nschool_effects), ltriangle(G[(3 * (j+1) - 2):(3 * (j+1))]))))
        }        
        rm(j)
        resB
    }
}
update.eta <- function(X, Y, Z, R.full.inv, ybetas, G, cons.logLik, Ny, nstudent, n_eta) {
    G.chol <- Cholesky(G)
    G.inv <- chol2inv(G.chol)
    H <- H.eta(G.inv = G.inv, Z = Z, R.full.inv = R.full.inv)
    chol.H <- Cholesky(H)
    H.inv <- chol2inv(chol.H)
    if (control$REML) {
        c.temp <- crossprod(X, R.full.inv) %*% Z
        c.1 <- rBind(crossprod(X, R.full.inv) %*% X, t(c.temp))
        c.2 <- rBind(c.temp, H)
        C_inv <- cBind(c.1, c.2)
        chol.C_inv <- Cholesky(forceSymmetric(symmpart(C_inv)))
        cs <- chol2inv(chol.C_inv)
        C.mat <- cs[-c(1:n_ybeta), -c(1:n_ybeta)]
        var.eta <- C.mat
    } else {
        var.eta <- H.inv
    }
    res <- var.eta
    sign.mult<-function(det){
    det$modulus*det$sign}    
    eta <- H.inv %*% as.vector(crossprod(Z, R.full.inv) %*% (Y - X %*% ybetas))
    log.p.eta <- -(n_eta/2) * log(2 * pi) - sign.mult(determinant(G.chol)) - 0.5 * crossprod(eta, G.inv) %*% eta
    log.p.y <- -(Ny/2) * log(2 * pi) +  sign.mult(determinant(Cholesky(R.full.inv)))  - 0.5 * crossprod(Y - X %*% ybetas - Z %*% eta, R.full.inv) %*% (Y - X %*% ybetas - Z %*% eta)
    if (control$REML) {
        attr(res, "likelihood") <- as.vector(cons.logLik + log.p.eta + log.p.y - sign.mult(determinant(chol.C_inv)))
    } else {
        attr(res, "likelihood") <- as.vector(cons.logLik + log.p.eta + log.p.y - sign.mult(determinant(chol.H)))
    }
    attr(res, "eta") <- eta
    attr(res,"h.inv")<-H.inv
    res
}
Score <- function(thetas, eta = eta.hat, ybetas, X, Y, Z, n_ybeta, n_eta, sqrt.W, inv.sqrt.W, nstudent, nteacher, Kg, cons.logLik, con = control, mis.list, pattern.parmlist2, pattern.count, pattern.length,
    pattern.Rtemplate, pattern.diag, pattern.key, Ny, pattern.sum = pattern.sum, persistence, P, alpha.diag, nalpha, alpha) {
    n_Rparm <- nyear.pseudo * (nyear.pseudo + 1)/2 - 1
    G <- thetas[(n_Rparm + 1):(n_Rparm + ncol(G.mat))]
    G <- reduce.G(G = G, nyear.score = nyear.score, nteacher = nteacher)
    if (persistence == "VP") {
        alpha.parm <- thetas[(n_Rparm + ncol(G.mat) + 1):length(thetas)]
        alpha[!((1:nalpha) %in% alpha.diag)] <- alpha.parm
        if(!control$school.effects){
Z <- Matrix(0, nrow(Z_mat), nteach_effects)
}else{
Z <- Matrix(0, nrow(Z_mat), nteach_effects+nschool_effects)
}
        for (i in 1:nalpha) {
            comp <- which(tril(ltriangle(1:nalpha)) == i, arr.ind = TRUE)
            Z <- Z + alpha[i] * P[[comp[1]]][[comp[2]]]
        }
     
     
if(control$school.effects){ 
colnames(Z)<-eta_effects
Z<-Z+Z.school.only

}
        Z<-drop0(Matrix(Z))
        #BZ structure
if(!control$school.effects){
Z.expand <- Matrix(0, nrow(Z_mat), 2 * nteach_effects)
Z.expand[, seq(1, 2 * nteach_effects, by = 2)] <- Z
}else{
Z.expand <- Matrix(0, nrow(Z_mat), 2 * nteach_effects + 2 * nschool_effects)
Z.expand[, seq(1, 2 * nteach_effects + 2*nschool_effects, by = 2)] <- Z
}
        colnames(Z) <- eta_effects
        J.Z <- rBind(Z.expand, B.Z.expand)
        J.Z <- J.Z[order(J.mat.original$student, J.mat.original$year), ]
        colnames(J.Z) <- interleave(colnames(Z), colnames(B.Z))
        for (p in unique(J.mat$pat)) {
            J.Z.p[[p]] <- J.Z[pat[[p]], , drop = FALSE]
        }
        Z <- J.Z
    }
    R_i <- ltriangle(c(as.vector(thetas[1:n_Rparm]), 1))
    R_i.parm <- c(as.vector(thetas[1:n_Rparm]), 1)
    LRI <- length(R_i.parm)
    R_i.parm.constrained <- R_i.parm[-LRI]
    if (length(mis.list) > 0) {
        R <- symmpart(suppressMessages(kronecker(suppressMessages(Diagonal(nstudent)), R_i)[-mis.list, -mis.list]))
        rinv<-Matrix(0,0,0)
R_i.inv<-chol2inv(chol(R_i))
for(i in 1:nstudent){
if(!any(mis.list%in%seq((i-1)*nyear.pseudo+1,i*nyear.pseudo))){
rinv<-bdiag(rinv,R_i.inv)
}else{
inv.indx<-which(!(seq((i-1)*nyear.pseudo+1,i*nyear.pseudo)%in%mis.list))
rinv<-bdiag(rinv,chol2inv(chol(R_i[inv.indx,inv.indx])))
}
}
        R_inv <- rinv
    } else {
        R <- symmpart(suppressMessages(kronecker(suppressMessages(Diagonal(nstudent)), R_i)))
        R_inv <- symmpart(suppressMessages(kronecker(suppressMessages(Diagonal(nstudent)), chol2inv(chol(R_i)))))
    }
    R.full <- drop0(symmpart(sqrt.W %*% R %*% sqrt.W))
    R.full.inv <- drop0(symmpart(inv.sqrt.W %*% R_inv %*% inv.sqrt.W))
    new.eta <- update.eta(X = X, Y = Y, Z = Z, R.full.inv = R.full.inv, ybetas = ybetas, G = G, cons.logLik = cons.logLik, Ny = Ny, nstudent = nstudent, n_eta = n_eta)
    eta.hat <- attr(new.eta, "eta")
    var.eta.hat <- new.eta
    temp_mat <- var.eta.hat + tcrossprod(eta.hat, eta.hat)
    temp_mat_R<-attr(new.eta,"h.inv")+ tcrossprod(eta.hat, eta.hat)
####
                pattern.sum <- list()
for (p in unique(J.mat$pat)) {
pattern.sum[[p]]<-Matrix(R_mstep2(invsqrtW_=as.matrix(diag(inv.sqrt.W)),JYp_=as.matrix(J.Y.p[[p]]),loopsize_=pattern.count[[p]]/pattern.length[[p]],
  patternlength_=pattern.length[[p]],rownumber_=as.matrix(J.Y.p.rownumber[[p]]),ybetas_=as.matrix(ybetas),
  etahat_=as.matrix(eta.hat),tempmatR_=as.matrix(temp_mat_R),
  JXpi_=as.matrix(J.X.p[[p]]@i),JXpp_=as.matrix(J.X.p[[p]]@p),JXpx_=as.matrix(J.X.p[[p]]@x),JXpdim_=as.matrix(J.X.p[[p]]@Dim),
  JZpi_=as.matrix(J.Z.p[[p]]@i),JZpp_=as.matrix(J.Z.p[[p]]@p),JZpx_=as.matrix(J.Z.p[[p]]@x),JZpdim_=as.matrix(J.Z.p[[p]]@Dim)))
}  
####
    score.R <- -pattern.f.score.constrained(R_i.parm.constrained = R_i.parm.constrained, nyear = nyear.pseudo, pattern.parmlist2 = pattern.parmlist2, pattern.count = pattern.count, pattern.length = pattern.length,
        pattern.Rtemplate = pattern.Rtemplate, pattern.diag = pattern.diag, pattern.key = pattern.key, pattern.sum = pattern.sum)
    G.originial <- G
    gam_t <- list()
    sv_gam_t <- list()
    index1 <- 0
    for (j in 1:nyear.score) {
        gam_t[[j]] <- matrix(0, Kg[j], Kg[j])
        G.originial_j <- G.originial[(index1 + 1):(index1 + nteacher[j] * Kg[j]), (index1 + 1):(index1 + nteacher[j] * Kg[j])]
        gam_t[[j]] <- G.originial_j[1:Kg[j], 1:Kg[j]]
        sv_gam_t[[j]] <- chol2inv(chol(gam_t[[j]]))
        index1 <- index1 + nteacher[j] * Kg[j]
    }
        if(control$school.effects){ 
        gam_t_school <- matrix(0, 2, 2)
        G.originial_j <- G.originial[(index1 + 1):(index1 + nschool_effects * 2), (index1 + 1):(index1 + nschool_effects * 2)]
        gam_t_school <- G.originial_j[1:2, 1:2]
        sv_gam_t_school <- chol2inv(chol(gam_t_school))
        index1 <- index1 + nschool_effects * 2
        }
    rm(j)
    gam_t_sc <- list()
    index1 <- 0
    score.G <- Matrix(0, 0, 0)
    for (j in 1:nyear.score) {
        gam_t_sc[[j]] <- matrix(0, Kg[j], Kg[j])
        temp_mat_j <- temp_mat[(index1 + 1):(index1 + nteacher[j] * Kg[j]), (index1 + 1):(index1 + nteacher[j] * Kg[j])]
        index2 <- c(1)
        for (k in 1:nteacher[j]) {
            gam_t_sc[[j]] <- gam_t_sc[[j]] + temp_mat_j[(index2):(index2 + Kg[j] - 1), (index2):(index2 + Kg[j] - 1)]
            index2 <- index2 + Kg[j]
        }
        index1 <- index1 + nteacher[j] * Kg[j]
        der <- -0.5 * (nteacher[j] * sv_gam_t[[j]] - sv_gam_t[[j]] %*% gam_t_sc[[j]] %*% sv_gam_t[[j]])
        if (is.numeric(drop(sv_gam_t[[j]]))) {
            score.eta.t <- der
        } else {
            score.eta.t <- 2 * der - diag(diag(der))
        }
        for (k in 1:nteacher[j]) {
            score.G <- bdiag(score.G, score.eta.t)
        }
    }
    
    if(control$school.effects){ 
       gam_t_school <- matrix(0, 2, 2)
        temp_mat_j <- temp_mat[(index1 + 1):(index1 + nschool_effects * 2), (index1 + 1):(index1 + nschool_effects*2)]
        index2 <- c(1)
        for (k in 1:nschool_effects) {
            gam_t_school <- gam_t_school + temp_mat_j[(index2):(index2 + 2 - 1), (index2):(index2 + 2 - 1)]
            index2 <- index2 + 2
        }
        index1 <- index1 + nschool_effects *2
        der <- -0.5 * (nschool_effects * sv_gam_t_school - sv_gam_t_school %*% gam_t_school %*% sv_gam_t_school)
        if (is.numeric(drop(sv_gam_t_school))) {
            score.eta.t <- der
        } else {
            score.eta.t <- 2 * der - diag(diag(der))
        }
        for (k in 1:nschool_effects) {
            score.G <- bdiag(score.G, score.eta.t)
        }
    }
    rm(j, k)
    if (persistence == "CP" | persistence == "ZP") {
        -c(score.R, reduce.G(G = score.G, nyear.score = nyear.score, nteacher = nteacher))
    } else if (persistence == "VP") {
        alpha.parm <- alpha[!((1:nalpha) %in% alpha.diag)]
        score.a <- alpha.score(alpha.parm, alpha = alpha, temp_mat = temp_mat_R[seq(1, n_eta, 2), seq(1, n_eta, 2)], nalpha = nalpha, alpha.diag = alpha.diag, P = P, R_inv = R.full.inv[J.mat$response ==
            "score", J.mat$response == "score"], eta.hat = eta.hat[seq(1, n_eta, 2)], ybetas = ybetas, X = X[J.mat$response == "score", ], Y = Y[J.mat$response == "score"])
        -c(score.R, reduce.G(G = score.G, nyear.score = nyear.score, nteacher = nteacher), score.a)
    }
}
hessian.f <- function() {
    lgLik.hist <- lgLik
    lgLik <- lgLik[it]
    Hessian <- NA
    std_errors <- c(rep(NA, length(thetas)))
    if (control$hessian == TRUE) {
        if (control$verbose)
            # cat('Calculating the Hessian...\n')
        flush.console()
        if (control$hes.method == "richardson") {
            Hessian <- ltriangle(ltriangle(jacobian(Score, thetas, eta = eta.hat, ybetas = ybetas, X = J.X, Y = J.Y, Z = J.Z, sqrt.W = sqrt.W, inv.sqrt.W = inv.sqrt.W, pattern.sum = pattern.sum, con = control,
                n_ybeta = n_ybeta, n_eta = n_eta, nstudent = nstudent, nteacher = nteacher, Kg = Kg, cons.logLik = cons.logLik, mis.list = mis.list, pattern.parmlist2 = pattern.parmlist2, pattern.count = pattern.count,
                pattern.length = pattern.length, pattern.Rtemplate = pattern.Rtemplate, pattern.diag = pattern.diag, pattern.key = pattern.key, Ny = Ny, persistence = persistence, P = P, alpha.diag = alpha.diag,
                nalpha = nalpha, alpha = alpha)))
        } else {
            Hessian <- ltriangle(ltriangle(jacobian(Score, thetas, method = "simple", eta = eta.hat, ybetas = ybetas, X = J.X, Y = J.Y, Z = J.Z, sqrt.W = sqrt.W, inv.sqrt.W = inv.sqrt.W, pattern.sum = pattern.sum,
                con = control, n_ybeta = n_ybeta, n_eta = n_eta, nstudent = nstudent, nteacher = nteacher, Kg = Kg, cons.logLik = cons.logLik, mis.list = mis.list, pattern.parmlist2 = pattern.parmlist2,
                pattern.count = pattern.count, pattern.length = pattern.length, pattern.Rtemplate = pattern.Rtemplate, pattern.diag = pattern.diag, pattern.key = pattern.key, Ny = Ny, persistence = persistence,
                P = P, alpha.diag = alpha.diag, nalpha = nalpha, alpha = alpha)))
        }
        std_errors <- try(c(sqrt(diag(solve(Hessian)))), silent = TRUE)
        hes.warn <- FALSE
        if (any(eigen(Hessian)$values <= 0)) {
            if (control$verbose)
                cat("Warning: Hessian not PD", "\n")
            std_errors <- c(rep(NA, length(thetas)))
            hes.warn <- TRUE
        }
    }
    chol2inv(chol(Hessian))
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
        res[1, ] else res
}
alpha.score <- function(alpha.parm, alpha, temp_mat_R, nalpha, alpha.diag, P, R_inv, eta.hat, ybetas, X, Y) {
    score.a <- c(NULL)
    alpha.s <- alpha
    alpha.s[!((1:nalpha) %in% alpha.diag)] <- alpha.parm
    alpha.set <- (1:nalpha)[!((1:nalpha) %in% alpha.diag)]
    if(!control$school.effects){
Z.a <- Matrix(0, nrow(Z_mat), nteach_effects)
}else{
Z.a <- Matrix(0, nrow(Z_mat), nteach_effects+nschool_effects)
}
    for (i in 1:nalpha) {
        comp <- which(tril(ltriangle(1:nalpha)) == i, arr.ind = TRUE)
        Z.a <- Z.a + alpha.s[i] * P[[comp[1]]][[comp[2]]]
    }
    if(control$school.effects){ 
colnames(Z.a)<-eta_effects
Z.a<-Z.a+Z.school.only

}
    for (i in alpha.set) {
        comp <- which(tril(ltriangle(1:nalpha)) == i, arr.ind = TRUE)
        score.a <- c(score.a, as.numeric((t(Y) - t(ybetas) %*% t(X)) %*% R_inv %*% P[[comp[1]]][[comp[2]]] %*% eta.hat - sum(diag(crossprod(Z.a, R_inv %*% P[[comp[1]]][[comp[2]]] %*% temp_mat_R)))))
    }
    score.a
}
#Need to make this more efficient like pattern.f.score.constrained
#pattern.f.score <- function(R_i.parm, nyear, pattern.parmlist2, pattern.count, pattern.length, pattern.Rtemplate, pattern.diag, pattern.key, pattern.sum) {
#    R_i <- ltriangle(as.vector(R_i.parm))
#    pattern.score <- numeric(nyear/2 * (nyear + 1))
#    for (r.parm in 1:(nyear/2 * (nyear + 1))) {
#        for (p in pattern.parmlist2[[r.parm]]) {
#            pat.coord <- pat.coord <- pat.coord.guide[[r.parm]][[p]]
#            pattern.y <- solve(pattern.f.R(R_i, p, nyear, pattern.key))
#            pattern.yc <- pattern.y[pat.coord]
#            pattern.score[r.parm] <- pattern.score[r.parm] - (pattern.countoverlength[[p]] * pattern.yc - (pattern.y %*% pattern.sum[[p]] %*% pattern.y)[pat.coord])
#        }
#        if (r.parm %in% pattern.diag)
#            pattern.score[r.parm] <- 0.5 * pattern.score[r.parm]
#    }
#    -pattern.score
#}
pattern.f.score.constrained <- function(R_i.parm.constrained, nyear, pattern.parmlist2, pattern.count, pattern.length, pattern.Rtemplate, pattern.diag, pattern.key, pattern.sum) {
    R_i <- as.matrix(ltriangle(c(as.vector(R_i.parm.constrained), 1)))
    pattern.score <- numeric(nyear/2 * (nyear + 1) - 1)
    for (p in nonempty.patterns) {
    pattern.y <- solve(pattern.f.R(R_i, p, nyear, pattern.key))
    YSY<-pattern.y %*% pattern.sum[[p]] %*% pattern.y
    PCL<-pattern.countoverlength[[p]]
        for (r.parm in 1:(nyear/2 * (nyear + 1) - 1)) {
            if(is.null(pat.coord <- pat.coord.guide[[r.parm]][[p]])) next
            pattern.yc <- pattern.y[pat.coord]
            pattern.score[r.parm] <- pattern.score[r.parm] - ( PCL* pattern.yc - YSY[pat.coord])

        }

    }
pattern.score[1:(nyear/2 * (nyear + 1) - 1) %in% pattern.diag] <- 0.5 * pattern.score[1:(nyear/2 * (nyear + 1) - 1) %in% pattern.diag]
    -pattern.score
}
pattern.f.R <- function(R, p, nyear, pattern.key) {
    R[pattern.key[p, ] * (1:nyear), pattern.key[p, ] * (1:nyear), drop = FALSE]
}
ptm.total <- proc.time()[3]
Z_mat$year <- as.numeric(Z_mat$year)
##
student.delete.list<-c(NULL)
for(s in unique(Z_mat$student)){
if(sum(!is.na(Z_mat[Z_mat$student==s,"y"]))==0) student.delete.list<-c(student.delete.list,s)
}
if(length(student.delete.list)>0){
Z_mat<-Z_mat[!Z_mat$student%in%student.delete.list,]
B.mat<-B.mat[!B.mat$student%in%student.delete.list,]
}
##
nyear.score <- length(unique(Z_mat$year))
B.year <- nyear.score + 1
B.mat$year <- B.year
B.mat$mis <- 0
B.mat$y <- B.mat$r
B.mat$r <- NULL
B.mat <- B.mat[order(B.mat$student, decreasing = TRUE), ]
Z_mat$mis <- rep(0, dim(Z_mat)[1])
student_list <- unique(Z_mat$student)
Z_mat[is.na(Z_mat$y), ]$mis <- rep(1, dim(Z_mat[is.na(Z_mat$y), ])[1])
for (g in 1:nyear.score) {
    mis_stu <- student_list[!(student_list %in% unique(Z_mat[Z_mat$year == g, ]$student))]
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
Z_mat.full <- Z_mat.full[which((Z_mat.full$student %in% Z_mat$student)), ]
B.mat[is.na(B.mat$y), ]$mis <- rep(1, dim(B.mat[is.na(B.mat$y), ])[1])
mis_stu <- student_list[!(student_list %in% unique(B.mat[B.mat$year == B.year, ]$student))]
le <- length(mis_stu)
if (le > 0) {
    temp.exp <- B.mat[1:le, ]
    temp.exp$year <- rep(B.year, le)
    temp.exp$mis <- rep(1, le)
    temp.exp$student <- mis_stu
    temp.exp$y <- rep(NA, le)
    B.mat <- rBind(B.mat, temp.exp)
}
rm(le)
B.mat.full <- B.mat
B.mat <- B.mat[!is.na(B.mat$y), ]
nstudent <- length(unique(Z_mat$student))
Kg <- rep(2, nyear.score)
Z_mat <- Z_mat[order(Z_mat$year, Z_mat$teacher), ]
Z_mat.full <- Z_mat.full[order(Z_mat.full$year, Z_mat.full$teacher), ]
na_list <- grep("^NA", Z_mat$teacher)
if (length(na_list) > 0) {
    teachyearcomb <- unique(cBind(Z_mat[-na_list, ]$year, Z_mat[-na_list, ]$teacher))
} else {
    teachyearcomb <- unique(cBind(Z_mat$year, Z_mat$teacher))
}
Z_mat <- Z_mat[order(Z_mat$student, Z_mat$year, Z_mat$teacher), , drop = FALSE]
Z_mat.full <- Z_mat.full[order(Z_mat.full$student, Z_mat.full$year, Z_mat.full$teacher), ]
nteach_effects <- dim(teachyearcomb)[1]
teacheffZ_mat <- Matrix(0, nrow = nrow(Z_mat), ncol = nteach_effects)
t_effects <- rep(NA, nteach_effects)
indx <- 1
eblup.tracker <- matrix(0, 0, 3)
dP <- list()
for (i in 1:nyear.score) {
    dP[[i]] <- list()
    for (j in 1:i) dP[[i]][[j]] <- matrix(0, 0, 2)
}
nalpha <- nyear.score/2 * (nyear.score + 1)
if (persistence == "ZP") {
    alpha <- ltriangle(diag(nyear.score))
} else if (persistence == "CP" | persistence == "VP") {
    alpha <- rep(1, nalpha)
}
alpha_key <- tril(ltriangle(1:nalpha))
alpha.diag <- diag(alpha_key)
alpha.parm <- alpha[!((1:nalpha) %in% alpha.diag)]
for (k in 1:nrow(teachyearcomb)) {
    student_subset <- Z_mat.full$student[Z_mat.full$year == teachyearcomb[k, 1] & Z_mat.full$teacher == teachyearcomb[k, 2]]
    #t_effects[k] <- paste(teachyearcomb[k, 1], "_", teachyearcomb[k, 2], sep = "")
    t_effects[k] <- paste( teachyearcomb[k, 2], sep = "")
    eblup.tracker <- rBind(eblup.tracker, c(teachyearcomb[k, ], teachyearcomb[k, 1]))
    for (yr in as.numeric(teachyearcomb[k, 1]):nyear.score) {
        if (sum(is.element(Z_mat$student, student_subset) & Z_mat$year == yr) != 0) {
            q1 <- (1:nrow(Z_mat))[is.element(Z_mat$student, student_subset) & Z_mat$year == yr & !is.na(Z_mat$y)]
            q2 <- rep(k, length(q1))
            dP[[yr]][[as.numeric(teachyearcomb[k, 1])]] <- rbind(dP[[yr]][[as.numeric(teachyearcomb[k, 1])]], cbind(q1, q2))
        }
    }
}
if(control$school.effects){ 
school.start<-nteach_effects+1
z.school<-t(fac2sparse(Z_mat$schoolID))
#school.length<-ncol(z.school)
Z.school.only<-cBind(Matrix(0,nrow(Z_mat),nteach_effects),z.school)
nschool_effects<-ncol(z.school)
}
P <- list()
for (i in 1:nyear.score) {
    P[[i]] <- list()
    for (j in 1:i) {
        P[[i]][[j]] <- as(sparseMatrix(i = dP[[i]][[j]][, 1], j = dP[[i]][[j]][, 2], dims = c(nrow(Z_mat), nteach_effects)), "dgCMatrix")
        if(control$school.effects)   P[[i]][[j]] <- cBind(P[[i]][[j]],Matrix(0,nrow(Z_mat),nschool_effects))
    }
}
if(!control$school.effects){
Z <- Matrix(0, nrow(Z_mat), nteach_effects)
}else{
Z <- Matrix(0, nrow(Z_mat), nteach_effects+nschool_effects)
}
for (i in 1:nalpha) {
    comp <- which(tril(ltriangle(1:nalpha)) == i, arr.ind = TRUE)
    Z <- Z + alpha[i] * P[[comp[1]]][[comp[2]]]
}

if(control$school.effects){ 
eta_effects<-colnames(Z)<-colnames(Z.school.only)<- c(t_effects,colnames(z.school))
colnames(Z)<-eta_effects
Z<-Z+Z.school.only
}else{
eta_effects<-colnames(Z)<-t_effects
}
temp_Z <- Z[order(Z_mat$year, Z_mat$student, decreasing = TRUE), ]
temp_Z_mat <- Z_mat[order(Z_mat$year, Z_mat$student, decreasing = TRUE), ]
temp_Z <- temp_Z[!(duplicated(temp_Z_mat$student)), ]
temp_Z_mat <- temp_Z_mat[!(duplicated(temp_Z_mat$student)), ]
temp_Z <- temp_Z[temp_Z_mat$student %in% B.mat$student, ]
temp_Z_mat <- temp_Z_mat[temp_Z_mat$student %in% B.mat$student, ]
#2.2 changed 1 to "student" below
temp_Z_mat <- temp_Z_mat[, "student", drop = FALSE]
B.Z <- temp_Z
colnames(B.Z) <- paste(colnames(temp_Z), "_B", sep = "")
#bznames<-colnames(B.Z)
#if(control$school.effects){ 
#bz.school<-t(fac2sparse(B.mat$schoolID))
#nschool_effects<-ncol(bz.school)
#B.Z<-cBind(B.Z,bz.school)
#colnames(B.Z)<-c(paste(t_effects, "_B", sep = ""),paste(colnames(bz.school),"_B",sep=""))
#}
B.mat <- merge(temp_Z_mat, B.mat, by = "student", sort = FALSE)
B.mat$response <- "B"
Z_mat$response <- "score"
B.mat.full$response <- "B"
Z_mat.full$response <- "score"
J.mat <- J.mat.original <- rbind(Z_mat[, c("student", "year", "mis", "y", "response")], B.mat[, c("student", "year", "mis", "y", "response")])
J.mat.full <- J.mat.full.original <- rbind(Z_mat.full[, c("student", "year", "mis", "y", "response")], B.mat.full[, c("student", "year", "mis", "y", "response")])
nyear.pseudo <- length(unique(J.mat$year))
X_mat <- sparse.model.matrix(control$score.fixed.effects, Z_mat, drop.unused.levels = TRUE)
X_mat <- X_mat[, !(colSums(abs(X_mat)) == 0), drop = FALSE]
if (rankMatrix(X_mat,method = 'qrLINPACK')[1] != dim(X_mat)[2]) {
    cat("WARNING: Fixed-effects design matrix for test scores is not full-rank", "\n")
    flush.console()
   ANSWER <- readline("Continue anyway? (Y/N)")
    if (substr(ANSWER, 1, 1) == "n"|substr(ANSWER, 1, 1) == "N") stop("WARNING: Fixed-effects design matrix not full-rank")
}
B_X_mat <- sparse.model.matrix(control$outcome.fixed.effects, B.mat, drop.unused.levels = TRUE)
B_X_mat <- B_X_mat[, !(colSums(abs(B_X_mat)) == 0), drop = FALSE]
if (rankMatrix(B_X_mat,method = 'qrLINPACK')[1] != dim(B_X_mat)[2]) {
    cat("WARNING: Fixed-effects design matrix for outcomes is not full-rank", "\n")
    flush.console()
   ANSWER <- readline("Continue anyway? (Y/N)")
    if (substr(ANSWER, 1, 1) == "n"|substr(ANSWER, 1, 1) == "N") stop("WARNING: Fixed-effects design matrix not full-rank")
}
J.X <- bdiag(X_mat, B_X_mat)
colnames(J.X) <- c(paste(colnames(X_mat), "_score", sep = ""), paste(colnames(B_X_mat), "_b", sep = ""))
if(!control$school.effects){
Z.expand <- Matrix(0, nrow(Z_mat), 2 * nteach_effects)
Z.expand[, seq(1, 2 * nteach_effects, by = 2)] <- Z
}else{
Z.expand <- Matrix(0, nrow(Z_mat), 2 * nteach_effects + 2 * nschool_effects)
Z.expand[, seq(1, 2 * nteach_effects + 2*nschool_effects, by = 2)] <- Z
}
if(!control$school.effects){
B.Z.expand <- Matrix(0, nrow(B.mat), 2 * nteach_effects)
B.Z.expand[, seq(2, 2 * nteach_effects, by = 2)] <- B.Z
}else{
B.Z.expand <- Matrix(0, nrow(B.mat), 2 * nteach_effects + 2 * nschool_effects)
B.Z.expand[, seq(2, 2 * nteach_effects + 2 * nschool_effects, by = 2)] <- B.Z
}
J.Z <- rBind(Z.expand, B.Z.expand)
J.Z <- J.Z[order(J.mat$student, J.mat$year), ]
J.X <- J.X[order(J.mat$student, J.mat$year), ]
interleave <- function(v1, v2) {
    ord1 <- 2 * (1:length(v1)) - 1
    ord2 <- 2 * (1:length(v2))
    c(v1, v2)[order(c(ord1, ord2))]
}
colnames(J.Z) <- interleave(colnames(Z), colnames(B.Z))
J.mat.full <- J.mat.full[order(J.mat.full$student, J.mat.full$year), ]
J.mat <- J.mat[order(J.mat$student, J.mat$year), ]
mis.list <- which(J.mat.full$mis == 1)
nteacher <- as.vector(tapply(teachyearcomb[, 2], teachyearcomb[, 1], length))
n_eta <- ncol(J.Z)
n_ybeta <- ncol(J.X)
J.Y <- as.vector(J.mat$y)
# huge.flag<-prod(dim(Z))>1e8 if(!huge.flag){ Z.dense <- as.matrix(Z) X.dense <- as.matrix(X) }
J.mat.full$r <- 1 - J.mat.full$mis
pattern.student <- matrix(J.mat.full$r, nstudent, nyear.pseudo, byrow = TRUE)
J.mat.full$pat <- rep(apply(pattern.student, 1, bin2dec), each = nyear.pseudo)
J.mat$pat <- J.mat.full[J.mat.full$r == 1, ]$pat
pat <- list()
pattern.count <- list()
pattern.length <- list()
pattern.countoverlength<-list()
J.Y.p <- list()
J.X.p <- list()
J.Z.p <- list()
J.Y.p.rownumber <- list()
J.mat$rownumber <- 1:nrow(J.mat)
pattern.key <- dec2bin(1:(2^nyear.pseudo - 1))
for (p in unique(J.mat$pat)) {
    pat[[p]] <- which(J.mat$pat == p)
    J.X.p[[p]] <- J.X[pat[[p]], , drop = FALSE] 
    J.Z.p[[p]] <- J.Z[pat[[p]], , drop = FALSE]
    J.Y.p[[p]] <- J.Y[pat[[p]]]
    J.Y.p.rownumber[[p]] <- J.mat$rownumber[pat[[p]]]
    pattern.count[[p]] <- length(J.Y.p[[p]])
    pattern.length[[p]] <- sum(pattern.key[p, ])
    pattern.countoverlength[[p]]<-pattern.count[[p]]/pattern.length[[p]]
}
# rm(Z.dense, X.dense)
pattern.yguide <- list()
for (g in 1:nyear.pseudo) {
    pattern.yguide[[g]] <- which(pattern.key[, g] == 1)
}
pattern.Rtemplate <- ltriangle(1:(nyear.pseudo/2 * (nyear.pseudo + 1)))
pattern.diag <- diag(pattern.Rtemplate)
pattern.parmlist1 <- list()
pattern.parmlist2 <- list()
for (p in unique(J.mat$pat)) {
    pattern.parmlist1[[p]] <- sort(unique(as.vector(pattern.f.R(pattern.Rtemplate, p, nyear.pseudo, pattern.key))))
}
for (r.parm in 1:(nyear.pseudo/2 * (nyear.pseudo + 1))) {
    pattern.parmlist2[[r.parm]] <- which(sapply(pattern.parmlist1, f <- function(x) r.parm %in% x))
}
eta.hat <- numeric(n_eta)
var.eta.hat <- Matrix(0, n_eta, n_eta)
R.temp.comp <- rep(1, nyear.pseudo)
for (g in 1:nyear.score) {
    R.temp.comp[g] <- var(J.mat[J.mat$year == g, ]$y)/6
}
R_i <- diag(R.temp.comp)
if (length(mis.list) > 0) {
    R <- symmpart(suppressMessages(kronecker(suppressMessages(Diagonal(nstudent)), R_i)[-mis.list, -mis.list]))
    rinv<-Matrix(0,0,0)
R_i.inv<-chol2inv(chol(R_i))
for(i in 1:nstudent){
if(!any(mis.list%in%seq((i-1)*nyear.pseudo+1,i*nyear.pseudo))){
rinv<-bdiag(rinv,R_i.inv)
}else{
inv.indx<-which(!(seq((i-1)*nyear.pseudo+1,i*nyear.pseudo)%in%mis.list))
rinv<-bdiag(rinv,chol2inv(chol(R_i[inv.indx,inv.indx])))
}
}

    R_inv <- rinv
} else {
    R <- symmpart(suppressMessages(kronecker(suppressMessages(Diagonal(nstudent)), R_i)))
    R_inv <- symmpart(suppressMessages(kronecker(suppressMessages(Diagonal(nstudent)), chol2inv(chol(R_i)))))
}
ybetas <- update.ybeta(J.X, J.Y, J.Z, R_inv, eta.hat)
ybetas[(ncol(X_mat) + 1):n_ybeta] <- glm.fit(x = B_X_mat, y = B.mat$y, family = binomial("probit"))[[1]]
names(ybetas) <- colnames(J.X)
G <- 100 * .symDiagonal(n_eta)
Ny <- sum(J.mat$mis == 0)
pat.coord.guide<-list()
for (r.parm in 1:(nyear.pseudo/2 * (nyear.pseudo + 1) - 1)) {
pat.coord.guide[[r.parm]]<-list()
        for (p in pattern.parmlist2[[r.parm]]) {
            pat.coord.guide[[r.parm]][[p]] <- which(tril(pattern.f.R(pattern.Rtemplate, p, nyear.pseudo, pattern.key)) == r.parm)

         }
}
nonempty.patterns<-NULL
for(p in 1:length(pat.coord.guide[[1]])){
if(is.null(retrieve.parm<-pattern.parmlist1[[p]])) next
nonempty.patterns<-c(nonempty.patterns,p)
}
         
if (control$REML) {
    cons.logLik <- 0.5 * (n_eta + n_ybeta) * log(2 * pi)
} else {
    cons.logLik <- 0.5 * (n_eta) * log(2 * pi)
}
  J.mat$fit <- J.mat$mu <- J.mat$mu.eta.val <- J.mat$nu <- J.mat$sqrt.w <- NA
    J.mat[J.mat$response == "B", ]$fit <- as.vector(J.X[J.mat$response == "B", ] %*% ybetas + J.Z[J.mat$response == "B", ] %*% eta.hat)
    J.mat[J.mat$response == "B", ]$mu <- as.vector(fam.binom$linkinv(J.mat[J.mat$response == "B", ]$fit))
    J.mat[J.mat$response == "B", ]$mu.eta.val <- fam.binom$mu.eta(J.mat[J.mat$response == "B", ]$fit)
    J.mat[J.mat$response == "B", ]$nu <- as.vector(J.mat[J.mat$response == "B", ]$fit) + (J.mat[J.mat$response == "B", ]$y - J.mat[J.mat$response == "B", ]$mu)/J.mat[J.mat$response == "B", ]$mu.eta.val
    J.mat[J.mat$response == "B", ]$sqrt.w <- sqrt(fam.binom$variance(J.mat[J.mat$response == "B", ]$mu)/J.mat[J.mat$response == "B", ]$mu.eta.val^2)
    J.mat[J.mat$response == "score", ]$sqrt.w <- 1
    J.mat[J.mat$response == "score", ]$nu <- J.mat[J.mat$response == "score", ]$y
    sqrt.W <- Diagonal(x = J.mat$sqrt.w)
    inv.sqrt.W <- Diagonal(x = 1/(J.mat$sqrt.w))
    R.full <- drop0(symmpart(sqrt.W %*% R %*% sqrt.W))
    R.full.inv <- drop0(symmpart(inv.sqrt.W %*% R_inv %*% inv.sqrt.W))
    J.Y <- J.mat$nu
    for (p in unique(J.mat$pat)) {
        J.Y.p[[p]] <- J.Y[pat[[p]]]
    }




for (PQL.it in 1:control$max.PQL.it) {
    if (PQL.it <= control$EM.to.NR) {
        EM.iter <- control$max.iter.EM[1]
    } else {
        EM.iter <- control$max.iter.EM[2]
    }
    Y.mat <- Matrix(0, EM.iter + 1, n_ybeta)
    G.mat <- Matrix(0, EM.iter + 1, length(reduce.G(G = G, nyear = nyear.score, nteacher = nteacher)))
    R.mat <- Matrix(0, EM.iter + 1, nyear.pseudo * (nyear.pseudo + 1)/2)
    lgLik <- numeric(EM.iter + 1)
    conv <- FALSE
    if (control$verbose)
        flush.console()

    if (PQL.it <= control$EM.to.NR)
        {
            cat("Beginning EM algorithm\n")
            for (it in 1:1000) {
                ptm <- proc.time()[3]
                suppressWarnings(rm(var.eta.hat, temp_mat,temp_mat_R))
                mresid <- as.numeric(J.Y - J.X %*% ybetas)
                cresid <- as.numeric(mresid - J.Z %*% eta.hat)
                yhat <- as.numeric(J.X %*% ybetas + J.Z %*% eta.hat)
                yhat.m <- as.numeric(J.X %*% ybetas)
                new.eta <- update.eta(X = J.X, Y = J.Y, Z = J.Z, R.full.inv = R.full.inv, ybetas = ybetas, G = G, cons.logLik = cons.logLik, Ny = Ny, nstudent = nstudent, n_eta = n_eta)
                eta.hat <- attr(new.eta, "eta")
                var.eta.hat <- new.eta
                temp_mat <- var.eta.hat + tcrossprod(eta.hat, eta.hat)
                temp_mat_R<-attr(new.eta,"h.inv")+ tcrossprod(eta.hat, eta.hat)
                loglikelihood <- attr(new.eta, "likelihood")
                lgLik[it] <- attr(new.eta, "likelihood")
                Y.mat[it, ] <- ybetas
                R.mat[it, ] <- suppressMessages(ltriangle(R_i))
                G.mat[it, ] <- suppressMessages(reduce.G(G, nyear.score, nteacher))
                rm(new.eta)
                thets1 <- c(Y.mat[it - 1, ], R.mat[it - 1, ], G.mat[it - 1, ])
                thets2 <- c(Y.mat[it, ], R.mat[it, ], G.mat[it, ])
                if (it >= 2) {
                  check.lik <- abs(lgLik[it] - lgLik[it - 1])/abs(lgLik[it] + control$tol1) < control$tol1
                  if (check.lik | it == (EM.iter)) {
                    conv <- TRUE
                    EM.conv <- TRUE
                    if (control$verbose) {
                      cat("\n\n End EM Algorithm.\n")
                      cat("\n\niter:", it, "\n")
                      cat("log-likelihood:", sprintf("%.7f", lgLik[it]), "\n")
                      cat("change in loglik:", sprintf("%.7f", lgLik[it] - lgLik[it - 1]), "\n")
                      cat("fixed effects:", round(ybetas, 4), "\n")
                      cat("R_i:\n")
                      print(round(as.matrix(R_i), 4))
                      cat("\n")
                      print(round(cov2cor(as.matrix(R_i)), 4))
                      cat("\n")
                      for (j in 1:nyear.score) {
                        cat("\ngamma_teach_year", j, "\n")
                        print(round(ltriangle(reduce.G(G, nyear.score, nteacher)[(3 * j - 2):(3 * j)]), 4))
                        cat("\n")
                        print(round(cov2cor(ltriangle(reduce.G(G, nyear.score, nteacher)[(3 * j - 2):(3 * j)])), 4))
                        cat("\n")
                        flush.console()
                      }
                      if(control$school.effects){
                      j<-j+1
                       cat("\n")
                       cat("\ngamma_school_effects\n")
                              print(round(ltriangle(reduce.G(G, nyear.score, nteacher)[(3 * j - 2):(3 * j)]), 4))
                        cat("\n")
                        print(round(cov2cor(ltriangle(reduce.G(G, nyear.score, nteacher)[(3 * j - 2):(3 * j)])), 4))
                        cat("\n")
                        flush.console()
                      }
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
                  cat("log-likelihood:", sprintf("%.7f", lgLik[it]), "\n")
                  cat("change in loglik:", sprintf("%.7f", lgLik[it] - lgLik[it - 1]), "\n")
                  cat("fixed effects:", round(ybetas, 4), "\n")
                  cat("R_i:\n")
                  print(round(as.matrix(R_i), 4))
                  cat("\n")
                  print(round(cov2cor(as.matrix(R_i)), 4))
                  cat("\n")
                  cat("G:\n")
                  for (j in 1:nyear.score) {
                    cat("\ngamma_teach_year", j, "\n")
                    print(round(ltriangle(reduce.G(G, nyear.score, nteacher)[(3 * j - 2):(3 * j)]), 4))
                    cat("\n")
                    print(round(cov2cor(ltriangle(reduce.G(G, nyear.score, nteacher)[(3 * j - 2):(3 * j)])), 4))
                        cat("\n")
                    flush.console()
                  }
                                        if(control$school.effects){
                      j<-j+1
                       cat("\n")
                       cat("\ngamma_school_effects\n")
                              print(round(ltriangle(reduce.G(G, nyear.score, nteacher)[(3 * j - 2):(3 * j)]), 4))
                        cat("\n")
                        print(round(cov2cor(ltriangle(reduce.G(G, nyear.score, nteacher)[(3 * j - 2):(3 * j)])), 4))
                        cat("\n")
                        flush.console()
                      }
                  cat("\n")
                  cat("alphas:\n")
                  print(round(alpha, 4))
                  rm(j)
                }
                flush.console()
                if (persistence == "VP") {
                  alpha.parm <- alpha[!((1:nalpha) %in% alpha.diag)]
                  alpha.parm.old <- alpha.parm
                  s.prev <- numeric(length(alpha.parm))
                  
                  # one step is sufficient because score function is linear in alpha
                  s <- alpha.score(alpha.parm, alpha = alpha, temp_mat_R = temp_mat_R[seq(1, n_eta, 2), seq(1, n_eta, 2)], nalpha = nalpha, alpha.diag = alpha.diag, P = P, R_inv = R.full.inv[J.mat$response ==
                    "score", J.mat$response == "score"], eta.hat = eta.hat[seq(1, n_eta, 2)], ybetas = ybetas, X = J.X[J.mat$response == "score", ], Y = J.Y[J.mat$response == "score"])
                  j <- jacobian(alpha.score, alpha.parm, method = "simple", alpha = alpha, temp_mat_R = temp_mat_R[seq(1, n_eta, 2), seq(1, n_eta, 2)], nalpha = nalpha, alpha.diag = alpha.diag, P = P, R_inv = R.full.inv[J.mat$response ==
                    "score", J.mat$response == "score"], eta.hat = eta.hat[seq(1, n_eta, 2)], ybetas = ybetas, X = J.X[J.mat$response == "score", ], Y = J.Y[J.mat$response == "score"])
                  hesprod <- solve(j, s)
                  alpha.parm <- alpha.parm - hesprod
                  alphan <- alpha
                  alphan[!((1:nalpha) %in% alpha.diag)] <- alpha.parm
                }
                # end update alphas
                indx <- 1
                Gn <- c(NULL)
                gt <- list()
                for (i in 1:nyear.score) {
                  gt[[i]] <- matrix(0, 2, 2)
                  for (j in 1:nteacher[i]) {
                    gt[[i]] <- gt[[i]] + temp_mat[(2 * (indx - 1) + 1):(2 * indx), (2 * (indx - 1) + 1):(2 * indx)]
                    indx <- indx + 1
                  }
                }
                if(control$school.effects){
                   gt.school <- matrix(0, 2, 2)
                   for (j in 1:nschool_effects) {
                    gt.school <- gt.school + temp_mat[(2 * (indx - 1) + 1):(2 * indx), (2 * (indx - 1) + 1):(2 * indx)]
                    indx <- indx + 1
                   }
                }
                Gn <- Matrix(0, 0, 0)
                for (i in 1:nyear.score) {
                  Gn <- bdiag(Gn, suppressMessages(kronecker(Diagonal(nteacher[i]), gt[[i]]/nteacher[i])))
                }
                if(control$school.effects){
                Gn <- bdiag(Gn,suppressMessages(kronecker(Diagonal(nschool_effects), gt.school/nschool_effects)))
                }

             
                pattern.sum <- list()
for (p in unique(J.mat$pat)) {
pattern.sum[[p]]<-R_mstep2(invsqrtW_=as.matrix(diag(inv.sqrt.W)),JYp_=as.matrix(J.Y.p[[p]]),loopsize_=pattern.count[[p]]/pattern.length[[p]],
  patternlength_=pattern.length[[p]],rownumber_=as.matrix(J.Y.p.rownumber[[p]]),ybetas_=as.matrix(ybetas),
  etahat_=as.matrix(eta.hat),tempmatR_=as.matrix(temp_mat_R),
  JXpi_=as.matrix(J.X.p[[p]]@i),JXpp_=as.matrix(J.X.p[[p]]@p),JXpx_=as.matrix(J.X.p[[p]]@x),JXpdim_=as.matrix(J.X.p[[p]]@Dim),
  JZpi_=as.matrix(J.Z.p[[p]]@i),JZpp_=as.matrix(J.Z.p[[p]]@p),JZpx_=as.matrix(J.Z.p[[p]]@x),JZpdim_=as.matrix(J.Z.p[[p]]@Dim))
}  
  

                R_i.parm <- ltriangle(as.matrix(R_i))
                LRI <- length(R_i.parm)
                R_i.parm.constrained <- R_i.parm[-LRI]
                R.cc <- 1
                hes.count <- 1
                s.prev <- numeric(length(R_i.parm))
                while (R.cc > 1e-04) {
                  s <- pattern.f.score.constrained(R_i.parm.constrained = R_i.parm.constrained, nyear = nyear.pseudo, pattern.parmlist2 = pattern.parmlist2, pattern.count = pattern.count, pattern.length = pattern.length,
                    pattern.Rtemplate = pattern.Rtemplate, pattern.diag = pattern.diag, pattern.key = pattern.key, pattern.sum = pattern.sum)
                  j <- jacobian(pattern.f.score.constrained, c(R_i.parm.constrained), method = "simple", nyear = nyear.pseudo, pattern.parmlist2 = pattern.parmlist2, pattern.count = pattern.count, pattern.length = pattern.length,
                    pattern.Rtemplate = pattern.Rtemplate, pattern.diag = pattern.diag, pattern.key = pattern.key, pattern.sum = pattern.sum)
                  if (it == 1 & (hes.count < 10)& PQL.it <= 3) {
                    hesprod <- solve(j + max(c(diag(j), 5)) * diag(LRI - 1), s)
                    R.cc <- s %*% s
                    if (hes.count == 9)
                      R.cc <- 0
                  } else if ((it <= 4) & (hes.count < 30) & PQL.it == 1) {
                    hesprod <- solve(j + max(c(diag(j), 5) * ((1 - hes.count/31)^2)) * diag(LRI - 1), s)
                    R.cc <- s %*% s
                    if (hes.count == 29)
                      R.cc <- 0
                  } else {
                    hesprod <- solve(j, s)
                    R.cc <- s %*% s
                  }
                  R_i.parm.constrained <- R_i.parm.constrained - hesprod
                  hes.count <- hes.count + 1
                  s.prev <- s
                }
                R_i.parm[-LRI] <- R_i.parm.constrained
                R_i <- ltriangle(R_i.parm)
                rm(R_i.parm)
                dimnames(R_i) <- list(NULL, NULL)
                if (length(mis.list) > 0) {
                  R <- symmpart(suppressMessages(kronecker(Diagonal(nstudent), R_i))[-mis.list, -mis.list])
                  rinv<-Matrix(0,0,0)
R_i.inv<-chol2inv(chol(R_i))
for(i in 1:nstudent){
if(!any(mis.list%in%seq((i-1)*nyear.pseudo+1,i*nyear.pseudo))){
rinv<-bdiag(rinv,R_i.inv)
}else{
inv.indx<-which(!(seq((i-1)*nyear.pseudo+1,i*nyear.pseudo)%in%mis.list))
rinv<-bdiag(rinv,chol2inv(chol(R_i[inv.indx,inv.indx])))
}
}
                  R_inv <- rinv
                } else {
                  R <- symmpart(suppressMessages(kronecker(Diagonal(nstudent), R_i)))
                  R_inv <- symmpart(suppressMessages(kronecker(Diagonal(nstudent), chol2inv(chol(R_i)))))
                }
                R.full <- drop0(symmpart(sqrt.W %*% R %*% sqrt.W))
                R.full.inv <- drop0(symmpart(inv.sqrt.W %*% R_inv %*% inv.sqrt.W))
                names(ybetas) <- colnames(J.X)
                G <- Gn
                if (persistence == "VP") {
                  alpha <- alphan
                  if(!control$school.effects){
Z <- Matrix(0, nrow(Z_mat), nteach_effects)
}else{
Z <- Matrix(0, nrow(Z_mat), nteach_effects+nschool_effects)
}
                  for (i in 1:nalpha) {
                    comp <- which(tril(ltriangle(1:nalpha)) == i, arr.ind = TRUE)
                    Z <- Z + alpha[i] * P[[comp[1]]][[comp[2]]]
                  }
                  if(control$school.effects){ 
colnames(Z)<-eta_effects
Z<-Z+Z.school.only

                  }
if(!control$school.effects){
Z.expand <- Matrix(0, nrow(Z_mat), 2 * nteach_effects)
Z.expand[, seq(1, 2 * nteach_effects, by = 2)] <- Z
}else{
Z.expand <- Matrix(0, nrow(Z_mat), 2 * nteach_effects + 2 * nschool_effects)
Z.expand[, seq(1, 2 * nteach_effects + 2*nschool_effects, by = 2)] <- Z
}
                  colnames(Z) <- eta_effects
                  # B.Z.expand <- Matrix(0, nrow(B.mat), 2 * nteach_effects) B.Z.expand[, seq(2, 2 * nteach_effects, by = 2)] <- B.Z
                  J.Z <- rBind(Z.expand, B.Z.expand)
                  J.Z <- J.Z[order(J.mat.original$student, J.mat.original$year), ]
                  colnames(J.Z) <- interleave(colnames(Z), colnames(B.Z))
                  for (p in unique(J.mat$pat)) {
                    J.Z.p[[p]] <- J.Z[pat[[p]], , drop = FALSE]
                  }
                }
                G.chol <- Cholesky(G)
                G.inv <- chol2inv(G.chol)
                R.inv.Z <- R.full.inv %*% J.Z
                V.1 <- chol2inv(Cholesky(G.inv + t(J.Z) %*% R.inv.Z))
                tX.Rinv.Z <- t(J.X) %*% R.inv.Z
                tX.Rinv.X <- t(J.X) %*% R.full.inv %*% J.X
                ybetas <- as.vector(chol2inv(Cholesky(forceSymmetric(symmpart(tX.Rinv.X - tX.Rinv.Z %*% V.1 %*% t(tX.Rinv.Z))))) %*% (t(J.X) %*% R.full.inv - tX.Rinv.Z %*% V.1 %*% t(R.inv.Z)) %*% J.Y)
                if (control$verbose)
                  cat("Iteration Time: ", proc.time()[3] - ptm, " seconds\n")
                flush.console()
            }
        }  #em loop  # End of EM
    nalpha.free <- sum(!((1:nalpha) %in% alpha.diag))
    if (control$NR == TRUE & PQL.it > control$EM.to.NR) {
            # NR loop
            cat("Beginning NR algorithm\n")
            convh <- abs(loglikelihood)
            nr.counter <- 0
            while (as.numeric(convh/abs(loglikelihood)) > 1e-09) {
                nr.counter <- nr.counter + 1
                new.eta <- update.eta(X = J.X, Y = J.Y, Z = J.Z, R.full.inv = R.full.inv, ybetas = ybetas, G = G, cons.logLik = cons.logLik, Ny = Ny, nstudent = nstudent, n_eta = n_eta)
                eta.hat <- attr(new.eta, "eta")
                var.eta.hat <- new.eta
                temp_mat <- var.eta.hat + tcrossprod(eta.hat, eta.hat)
                temp_mat_R<-attr(new.eta,"h.inv")+ tcrossprod(eta.hat, eta.hat)
                cat("\n", attr(new.eta, "likelihood") - loglikelihood, "\n")
                loglikelihood <- attr(new.eta, "likelihood")
              #  if (nr.counter == 1 | persistence == "VP") {
                 if(TRUE){
                pattern.sum <- list()
for (p in unique(J.mat$pat)) {
pattern.sum[[p]]<-Matrix(R_mstep2(invsqrtW_=as.matrix(diag(inv.sqrt.W)),JYp_=as.matrix(J.Y.p[[p]]),loopsize_=pattern.count[[p]]/pattern.length[[p]],
  patternlength_=pattern.length[[p]],rownumber_=as.matrix(J.Y.p.rownumber[[p]]),ybetas_=as.matrix(ybetas),
  etahat_=as.matrix(eta.hat),tempmatR_=as.matrix(temp_mat_R),
  JXpi_=as.matrix(J.X.p[[p]]@i),JXpp_=as.matrix(J.X.p[[p]]@p),JXpx_=as.matrix(J.X.p[[p]]@x),JXpdim_=as.matrix(J.X.p[[p]]@Dim),
  JZpi_=as.matrix(J.Z.p[[p]]@i),JZpp_=as.matrix(J.Z.p[[p]]@p),JZpx_=as.matrix(J.Z.p[[p]]@x),JZpdim_=as.matrix(J.Z.p[[p]]@Dim)))
}  
  

                }
                names(ybetas) <- colnames(J.X)
                if (persistence == "CP" | persistence == "ZP") {
                  thetas <- c(ltriangle(as.matrix(R_i))[-LRI], reduce.G(G = G, nyear.score = nyear.score, nteacher = nteacher))
                } else if (persistence == "VP") {
                  thetas <- c(ltriangle(as.matrix(R_i))[-LRI], reduce.G(G = G, nyear.score = nyear.score, nteacher = nteacher), alpha[!((1:nalpha) %in% alpha.diag)])
                }
                thetas.old <- thetas
                ybetas.old <- ybetas
                Hessian.inv <- hessian.f()
                score.thetas <- Score(thetas = thetas, eta = eta.hat, ybetas = ybetas, X = J.X, Y = J.Y, Z = J.Z, sqrt.W = sqrt.W, inv.sqrt.W = inv.sqrt.W, pattern.sum = pattern.sum, con = control, n_ybeta = n_ybeta,
                  n_eta = n_eta, nstudent = nstudent, nteacher = nteacher, Kg = Kg, cons.logLik = cons.logLik, mis.list = mis.list, pattern.parmlist2 = pattern.parmlist2, pattern.count = pattern.count,
                  pattern.length = pattern.length, pattern.Rtemplate = pattern.Rtemplate, pattern.diag = pattern.diag, pattern.key = pattern.key, Ny = Ny, persistence = persistence, P = P, alpha.diag = alpha.diag,
                  nalpha = nalpha, alpha = alpha)
                thetas <- thetas - Hessian.inv %*% score.thetas
                convh <- t(score.thetas) %*% Hessian.inv %*% score.thetas
                # score.ybetas<-t(J.X)%*%R.full.inv%*%(J.Y-J.X%*%ybetas-J.Z%*%eta.hat)
                n_Rparm <- nyear.pseudo * (nyear.pseudo + 1)/2 - 1
                G <- thetas[(n_Rparm + 1):(n_Rparm + ncol(G.mat))]
                G <- reduce.G(G = G, nyear.score = nyear.score, nteacher = nteacher)
                if (persistence == "VP") {
                  alpha.parm <- thetas[(n_Rparm + ncol(G.mat) + 1):length(thetas)]
                  alpha[!((1:nalpha) %in% alpha.diag)] <- alpha.parm
                 if(!control$school.effects){
Z <- Matrix(0, nrow(Z_mat), nteach_effects)
}else{
Z <- Matrix(0, nrow(Z_mat), nteach_effects+nschool_effects)
}
                  for (i in 1:nalpha) {
                    comp <- which(tril(ltriangle(1:nalpha)) == i, arr.ind = TRUE)
                    Z <- Z + alpha[i] * P[[comp[1]]][[comp[2]]]
                  }
                  if(control$school.effects){ 
       colnames(Z)<-eta_effects
       Z<-Z+Z.school.only

}
if(!control$school.effects){
Z.expand <- Matrix(0, nrow(Z_mat), 2 * nteach_effects)
Z.expand[, seq(1, 2 * nteach_effects, by = 2)] <- Z
}else{
Z.expand <- Matrix(0, nrow(Z_mat), 2 * nteach_effects + 2 * nschool_effects)
Z.expand[, seq(1, 2 * nteach_effects + 2*nschool_effects, by = 2)] <- Z
}
                  colnames(Z) <- eta_effects
                  J.Z <- rBind(Z.expand, B.Z.expand)
                  J.Z <- J.Z[order(J.mat.original$student, J.mat.original$year), ]
                  colnames(J.Z) <- interleave(colnames(Z), colnames(B.Z))
                  for (p in unique(J.mat$pat)) {
                    J.Z.p[[p]] <- J.Z[pat[[p]], , drop = FALSE]
                  }
                  # if(!huge.flag){Z.dense <- as.matrix(Z)} for (p in unique(Z_mat$pat)) { if(!huge.flag){ Z.p[[p]] <- Z.dense[pat[[p]], , drop = FALSE]}else{ Z.p[[p]] <- Z[pat[[p]], , drop = FALSE] } } rm(Z.dense)
                }
                R_i <- ltriangle(c(as.vector(thetas[1:n_Rparm]), 1))
                R_i.parm <- c(as.vector(thetas[1:n_Rparm]), 1)
                LRI <- length(R_i.parm)
                R_i.parm.constrained <- R_i.parm[-LRI]
                if (length(mis.list) > 0) {
                  R <- symmpart(suppressMessages(kronecker(suppressMessages(Diagonal(nstudent)), R_i)[-mis.list, -mis.list]))
                  rinv<-Matrix(0,0,0)
R_i.inv<-chol2inv(chol(R_i))
for(i in 1:nstudent){
if(!any(mis.list%in%seq((i-1)*nyear.pseudo+1,i*nyear.pseudo))){
rinv<-bdiag(rinv,R_i.inv)
}else{
inv.indx<-which(!(seq((i-1)*nyear.pseudo+1,i*nyear.pseudo)%in%mis.list))
rinv<-bdiag(rinv,chol2inv(chol(R_i[inv.indx,inv.indx])))
}
}

                  R_inv <- rinv
                } else {
                  R <- symmpart(suppressMessages(kronecker(suppressMessages(Diagonal(nstudent)), R_i)))
                  R_inv <- symmpart(suppressMessages(kronecker(suppressMessages(Diagonal(nstudent)), chol2inv(chol(R_i)))))
                }
                R.full <- drop0(symmpart(sqrt.W %*% R %*% sqrt.W))
                R.full.inv <- drop0(symmpart(inv.sqrt.W %*% R_inv %*% inv.sqrt.W))
                G.chol <- Cholesky(G)
                G.inv <- chol2inv(G.chol)
                R.inv.Z <- R.full.inv %*% J.Z
                V.1 <- chol2inv(Cholesky(G.inv + t(J.Z) %*% R.inv.Z))
                tX.Rinv.Z <- t(J.X) %*% R.inv.Z
                tX.Rinv.X <- t(J.X) %*% R.full.inv %*% J.X
                ybetas <- as.vector(chol2inv(Cholesky(symmpart(tX.Rinv.X - tX.Rinv.Z %*% V.1 %*% t(tX.Rinv.Z)))) %*% (t(J.X) %*% R.full.inv - tX.Rinv.Z %*% V.1 %*% t(R.inv.Z)) %*% J.Y)
            }
            new.eta <- update.eta(X = J.X, Y = J.Y, Z = J.Z, R.full.inv = R.full.inv, ybetas = ybetas, G = G, cons.logLik = cons.logLik, Ny = Ny, nstudent = nstudent, n_eta = n_eta)
            eta.hat <- attr(new.eta, "eta")
            var.eta.hat <- new.eta
            temp_mat <- var.eta.hat + tcrossprod(eta.hat, eta.hat)
            temp_mat_R<-attr(new.eta,"h.inv")+ tcrossprod(eta.hat, eta.hat)
            cat("\n", attr(new.eta, "likelihood") - loglikelihood, "\n")
        }  #end NR
    if (persistence == "CP" | persistence == "ZP") {
        thetas.pql <- c(ltriangle(as.matrix(R_i))[-LRI], reduce.G(G = G, nyear.score = nyear.score, nteacher = nteacher))
    } else if (persistence == "VP") {
        thetas.pql <- c(ltriangle(as.matrix(R_i))[-LRI], reduce.G(G = G, nyear.score = nyear.score, nteacher = nteacher), alpha[!((1:nalpha) %in% alpha.diag)])
    }
    parms.new <- c(ybetas, thetas.pql)
    if (PQL.it > 1) {
        max.change <- 2 * max(abs(parms.new - parms.old)/(abs(parms.new) + abs(parms.old) + 1e-06))
        if (control$verbose) {
            cat("\nMaximum relative parameter change between PQL iterations: ", max.change, "\n")
        }
        if (max.change < control$pconv)
            break
    }
    thetas.pql.old <- thetas.pql
    ybetas.pql.old <- ybetas
    parms.old <- c(ybetas.pql.old, thetas.pql.old)
    # we use inverse of Wolfinger's W this part of the code transforms the responses (pseudo-likelihood approach) most of this code was taken from the R function glmmPQL
    J.mat[J.mat$response == "B", ]$fit <- as.vector(J.X[J.mat$response == "B", ] %*% ybetas + J.Z[J.mat$response == "B", ] %*% eta.hat)
    J.mat[J.mat$response == "B", ]$mu <- as.vector(fam.binom$linkinv(J.mat[J.mat$response == "B", ]$fit))
    J.mat[J.mat$response == "B", ]$mu.eta.val <- fam.binom$mu.eta(J.mat[J.mat$response == "B", ]$fit)
    J.mat[J.mat$response == "B", ]$nu <- as.vector(J.mat[J.mat$response == "B", ]$fit) + (J.mat[J.mat$response == "B", ]$y - J.mat[J.mat$response == "B", ]$mu)/J.mat[J.mat$response == "B", ]$mu.eta.val
    J.mat[J.mat$response == "B", ]$sqrt.w <- sqrt(fam.binom$variance(J.mat[J.mat$response == "B", ]$mu)/J.mat[J.mat$response == "B", ]$mu.eta.val^2)
    J.mat[J.mat$response == "score", ]$sqrt.w <- 1
    J.mat[J.mat$response == "score", ]$nu <- J.mat[J.mat$response == "score", ]$y
    sqrt.W <- Diagonal(x = J.mat$sqrt.w)
    inv.sqrt.W <- Diagonal(x = 1/(J.mat$sqrt.w))
    R.full <- drop0(symmpart(sqrt.W %*% R %*% sqrt.W))
    R.full.inv <- drop0(symmpart(inv.sqrt.W %*% R_inv %*% inv.sqrt.W))
    J.Y <- J.mat$nu
    for (p in unique(J.mat$pat)) {
        J.Y.p[[p]] <- J.Y[pat[[p]]]
    }
}  #pql loop

names(ybetas) <- colnames(J.X)
if (persistence == "CP" | persistence == "ZP") {
    thetas <- c(ltriangle(as.matrix(R_i))[-LRI], reduce.G(G = G, nyear.score = nyear.score, nteacher = nteacher))
} else if (persistence == "VP") {
    thetas <- c(ltriangle(as.matrix(R_i))[-LRI], reduce.G(G = G, nyear.score = nyear.score, nteacher = nteacher), alpha[!((1:nalpha) %in% alpha.diag)])
}
lgLik.hist <- lgLik
lgLik <- lgLik[it]
Hessian <- NA
std_errors <- c(rep(NA, length(thetas)))
if (control$hessian == TRUE) {
    if (control$verbose)
        cat("Calculating Hessian of the variance components...")
    flush.console()
                pattern.sum <- list()
for (p in unique(J.mat$pat)) {
pattern.sum[[p]]<-Matrix(R_mstep2(invsqrtW_=as.matrix(diag(inv.sqrt.W)),JYp_=as.matrix(J.Y.p[[p]]),loopsize_=pattern.count[[p]]/pattern.length[[p]],
  patternlength_=pattern.length[[p]],rownumber_=as.matrix(J.Y.p.rownumber[[p]]),ybetas_=as.matrix(ybetas),
  etahat_=as.matrix(eta.hat),tempmatR_=as.matrix(temp_mat_R),
  JXpi_=as.matrix(J.X.p[[p]]@i),JXpp_=as.matrix(J.X.p[[p]]@p),JXpx_=as.matrix(J.X.p[[p]]@x),JXpdim_=as.matrix(J.X.p[[p]]@Dim),
  JZpi_=as.matrix(J.Z.p[[p]]@i),JZpp_=as.matrix(J.Z.p[[p]]@p),JZpx_=as.matrix(J.Z.p[[p]]@x),JZpdim_=as.matrix(J.Z.p[[p]]@Dim)))
}  
  
    if (control$hes.method == "richardson") {
        Hessian <- ltriangle(ltriangle(jacobian(Score, thetas, eta = eta.hat, ybetas = ybetas, X = J.X, Y = J.Y, Z = J.Z, sqrt.W = sqrt.W, inv.sqrt.W = inv.sqrt.W, pattern.sum = pattern.sum, con = control,
            n_ybeta = n_ybeta, n_eta = n_eta, nstudent = nstudent, nteacher = nteacher, Kg = Kg, cons.logLik = cons.logLik, mis.list = mis.list, pattern.parmlist2 = pattern.parmlist2, pattern.count = pattern.count,
            pattern.length = pattern.length, pattern.Rtemplate = pattern.Rtemplate, pattern.diag = pattern.diag, pattern.key = pattern.key, Ny = Ny, persistence = persistence, P = P, alpha.diag = alpha.diag,
            nalpha = nalpha, alpha = alpha)))
    } else {
        Hessian <- ltriangle(ltriangle(jacobian(Score, thetas, method = "simple", eta = eta.hat, ybetas = ybetas, X = J.X, Y = J.Y, Z = J.Z, sqrt.W = sqrt.W, inv.sqrt.W = inv.sqrt.W, pattern.sum = pattern.sum,
            con = control, n_ybeta = n_ybeta, n_eta = n_eta, nstudent = nstudent, nteacher = nteacher, Kg = Kg, cons.logLik = cons.logLik, mis.list = mis.list, pattern.parmlist2 = pattern.parmlist2, pattern.count = pattern.count,
            pattern.length = pattern.length, pattern.Rtemplate = pattern.Rtemplate, pattern.diag = pattern.diag, pattern.key = pattern.key, Ny = Ny, persistence = persistence, P = P, alpha.diag = alpha.diag,
            nalpha = nalpha, alpha = alpha)))
    }
    std_errors <- try(c(sqrt(diag(solve(Hessian)))), silent = TRUE)
    hes.warn <- FALSE
    if (any(eigen(Hessian)$values <= 0)) {
        if (control$verbose)
            cat("Warning: Hessian not PD", "\n")
        std_errors <- c(rep(NA, length(thetas)))
        hes.warn <- TRUE
    }
}
c.temp <- crossprod(J.X, R.full.inv) %*% J.Z
c.1 <- rBind(crossprod(J.X, R.full.inv) %*% J.X, t(c.temp))
G.inv <- chol2inv(Cholesky(G))
c.2 <- rBind(c.temp, H.eta(G.inv, J.Z, R.full.inv))
C_inv <- cBind(c.1, c.2)
C <- solve(C_inv)
eblup_stderror <- sqrt(diag(C)[-c(1:n_ybeta)])
ybetas_stderror <- sqrt(diag(C)[1:n_ybeta])
rm(C, C_inv, c.2, c.1, c.temp)
eblup <- as.matrix(cBind(eta.hat, eblup_stderror))
eblup <- as.data.frame(eblup)
eblup <- as.data.frame(cbind(colnames(J.Z), eblup))
colnames(eblup) <- c("effect", "EBLUP", "std_error")
t_lab <- as.vector(NULL)
r_lab <- as.vector(NULL)
for (j in 1:nyear.score) {
    ne <- (Kg[j] * (Kg[j] + 1))/2
    y <- c(NULL)
    x <- c(NULL)
    for (k in 1:Kg[j]) {
        x <- c(x, k:Kg[j])
        y <- c(y, rep(k, (Kg[j] - k + 1)))
    }
    t_lab <- c(t_lab, paste("teacher effect from year", rep(j, ne), sep = ""))
}
if(control$school.effects){
t_lab<-c(t_lab,rep("School Effect",3))
}
y <- c(NULL)
x <- c(NULL)
for (k in 1:nyear.pseudo) {
    x <- c(x, k:nyear.pseudo)
    y <- c(y, rep(k, (nyear.pseudo - k + 1)))
}
r_lab <- paste("error covariance", ":[", x, ",", y, "]", sep = "")
rm(j, ne)
alpha.label <- c(NULL)
for (i in 1:(nyear.score - 1)) {
    for (j in (i + 1):(nyear.score)) {
        alpha.label <- c(alpha.label, paste("alpha_", j, i, sep = ""))
    }
}
if (persistence == "CP" | persistence == "ZP") {
    effect_la <- c(names(ybetas), r_lab, t_lab)
} else if (persistence == "VP") {
    effect_la <- c(names(ybetas), r_lab, t_lab, alpha.label)
}
constant.element <- which(effect_la == paste("error covariance", ":[", nyear.pseudo, ",", nyear.pseudo, "]", sep = ""))
if (control$hessian == TRUE) {
    parameters <- round(cBind(append(c(ybetas, thetas), 1, after = constant.element - 1), append(c(ybetas_stderror, std_errors), NA, after = constant.element - 1)), 4)
    colnames(parameters) <- c("Estimate", "Standard Error")
    rownames(parameters) <- as.character(effect_la)
}               
if (control$hessian == FALSE) {
    parameters <- round(cBind(append(c(ybetas, thetas), 1, after = constant.element - 1), c(ybetas_stderror, rep(NA, length(thetas) + 1))), 4)
    colnames(parameters) <- c("Estimate", "Standard Error")
    rownames(parameters) <- as.character(effect_la)
}
if (control$verbose) cat("done.\n")
mresid <- as.numeric(J.Y - J.X %*% ybetas)
cresid <- as.numeric(mresid - J.Z %*% eta.hat)
yhat <- as.numeric(J.X %*% ybetas + J.Z %*% eta.hat)
yhat.m <- as.numeric(J.X %*% ybetas)
gam_t <- list()
for (i in 1:nyear.score) {
    gam_t[[i]] <- as.matrix(ltriangle(reduce.G(G, nyear.score = nyear.score, nteacher = nteacher)[(1+3*(i-1)):(3*i)]))
    colnames(gam_t[[i]]) <- c(paste("year", i, sep = ""),"")
}
persistence_parameters = ltriangle(alpha)
persistence_parameters[upper.tri(persistence_parameters)]<-NA
if(!control$school.effects){
school.subset<-NULL
}else{
school.subset<-eblup[(2*nteach_effects+1):n_eta,]
}
teach.cov <- lapply(gam_t, function(x) round(x, 4))
res <- list(loglik = lgLik, teach.effects = eblup[1:(2*nteach_effects),],school.effects=school.subset, parameters = parameters, Hessian = Hessian, R_i = as.matrix(R_i), teach.cov = gam_t, mresid = mresid, cresid = cresid, y = J.Y, yhat = yhat, 
    num.obs = Ny, num.student = nstudent, num.year = nyear.score, num.teach = nteacher,   persistence = control$persistence, persistence_parameters =  persistence_parameters,X=J.X,Z=J.Z,
    G=G,R=R, R.full=R.full,sqrt.W=J.mat$sqrt.w)
# res$teach.effects$teacher_year <- key[match(res$teach.effects$teacher_year, key[, 2]), 1] res$teach.effects$effect_year <- key[match(res$teach.effects$effect_year, key[, 2]), 1]
class(res) <- "RealVAMS"
cat("Total Time: ", proc.time()[3] - ptm.total, " seconds\n")
return(res)
}
