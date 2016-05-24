PB_cre <-
function(Z_mat=Z_mat, first.order=first.order,home.field=home.field, control=control,game.effect=game.effect){
teams <- sort(unique(c(Z_mat$home,Z_mat$away)))
nteams <- length(teams)


if(home.field&!control$OT.flag){
j_fixed_effects <- formula(~Location+0)
}else if(home.field&control$OT.flag){
j_fixed_effects <- formula(~Location+(OT)+0)   
}else if(!home.field&control$OT.flag){
j_fixed_effects <- formula(~(OT)+0) 
}else{
j_fixed_effects <- formula(~1)
}


r_fixed_effects <- formula(~neutral.site)
#iter.EM - maximum number of EM iterations
#tol1 - convergence crieterion for the first order Laplace approximation.
#       tol2 refers to the maximum relative change in parameters between iterations
#       The first order Laplace approximation runs until tol1 signals,
#       at which point the fully exponential corrections for eta begin
#tol2 - The fully exponential iterations run until tol2 signals (maximum relative
#        change in model paramters)
#tolFE - The algorithm runs with the fully exponential corrections only to eta until tolFE
#       signals (maximum relative change in parameters). After this, the fully exponential
#       corrections for both eta and var.eta are calculated
#verbose - controls output


#------------------
#First derivative of complete data log-likelihood (L(eta) in the paper)
gr.eta <- function(eta, J_X, J_Y, J_Z, jbetas, G.inv, n_eta, R_Z,R_Y,R_X,rbetas) {
    gr.j <- as.vector(as.vector( jder1(J_X %*% jbetas + J_Z %*% eta,J_Y))%*%J_Z)
    gr.r <- colSums(R_Z * as.vector((-1)^(1 - R_Y) * rder1((-1)^(1 - R_Y) * (R_X %*% rbetas + R_Z %*% eta))))
gr.p.eta <- -G.inv %*% eta
    -as.vector(gr.j +gr.r+ gr.p.eta)
}
#--------------------
#Second derivative of complete data log-likelihood (\Sigma^w in the paper)
H.eta <- function(eta, jbetas, J_X, J_Y, J_Z,G.inv, n_eta,R_Y,R_X,R_Z,rbetas) {
    h.eta <- G.inv
    temp <- -jder2(J_X %*% jbetas + J_Z %*% eta)
    h.j <- Matrix(crossprod(J_Z, temp * J_Z))
    temp.r <- -rder2((-1)^(1 - R_Y) * (R_X %*% rbetas + R_Z %*% eta))
    h.r <- Matrix(crossprod(R_Z, temp.r * R_Z))

    Matrix(h.eta + h.j+h.r)
}
#---------------------------------------------------
log.factorial<-function(x) {
if(x==0){ res<-1
}else{
res<-0
for(i in 2:x){
res<-res+log(i)}
}
res
}
log.factorial.v<-Vectorize(log.factorial,"x")
#---------------------------------------------------
#ltriangle: extract lower triangle
#and rebuild from elements
#note ltriangle(ltriangle(x))=x
ltriangle <- function(x) {
    if (!is.null(dim(x)[2])) {
        resA <- as.vector(x[lower.tri(x, diag = TRUE)])
        resA
    } else {
        nx <- length(x)
        d <- 0.5 * (-1 + sqrt(1 + 8 * nx))
        resB <- Diagonal(d)
        resB[lower.tri(resB, diag = TRUE)] <- x
        if (nx > 1) {
            resB <- resB + t(resB) - diag(diag(resB))
        }
        resB
    }
}
#--------------
#These functions are derivatives from Section 3.1
jder0 <- function(q,J_Y) {
    q <- as.vector(q)
    -log.factorial.v(J_Y)+J_Y*q-exp(q)
}
jder1 <- function(q,J_Y) {
    q <- as.vector(q)
   J_Y-exp(q)
}
jder2 <- function(q) {
    q <- as.vector(q)
    -exp(q)
}
jder3 <- function(q) {
    q <- as.vector(q)
    -exp(q)
    }
jder4 <- function(q) {
    q <- as.vector(q)
    -exp(q)
}
#These functions are derivatives of the standard normal CDF that are used in the missing data mechanism
rder0 <- function(q) {
    q <- as.vector(q)
    pnorm(q, log.p = TRUE)
}
rder1 <- function(q) {
    q <- as.vector(q)
    dnorm(q)/pmax(pnorm(q), 1e-15)
}
rder2 <- function(q) {
    q <- as.vector(q)
    (-(q/sqrt(2 * pi) * exp(-q^2/2)) * pnorm(q) - dnorm(q)^2)/pmax(pnorm(q)^2, 1e-15)
}
rder3 <- function(q) {
    q <- as.vector(q)
    ((q^2 - 1)/sqrt(2 * pi) * exp(-q^2/2) * pnorm(q)^2 + 3 * (q/sqrt(2 * pi) * exp(-q^2/2)) * pnorm(q) * dnorm(q) + 2 * dnorm(q)^3)/pmax(pnorm(q)^3, 1e-15)
}
rder4 <- function(q) {
    q <- as.vector(q)
    ((-q * (q^2 - 3)/sqrt(2 * pi) * exp(-q^2/2)) * pnorm(q)^3 - 4 * ((q^2 - 1)/sqrt(2 * pi) * exp(-q^2/2)) * pnorm(q)^2 * dnorm(q) - 12 * (q/sqrt(2 * pi) * exp(-q^2/2)) * pnorm(q) * dnorm(q)^2 - 3 * (q/sqrt(2 * pi) * exp(-q^2/2))^2 * pnorm(q)^2 - 6 * 
        dnorm(q)^4)/pmax(pnorm(q)^4, 1e-15)
}
 if(!game.effect) reduce.G<-function(G) ltriangle(as.matrix(G[1:3,1:3]))
 if(game.effect) reduce.G<-function(G) c(ltriangle(as.matrix(G[1:3,1:3])),G[n_eta,n_eta])
#-------------------
#finds the values of the random effects that maximize the
#EM-complete data likelihood. These are the EBLUPS (must be
#corrected if non-normal data used)
update.eta <- function(eta, J_X, J_Y, J_Z, jbetas, G, nyear, n_eta, cons.logLik, R_Z,R_Y,R_X,rbetas) {
    G.chol <- chol(G)
    #For G, chol2inv is faster than solve
    G.inv <- chol2inv(G.chol)
    eta <- as.vector(eta)
    H <- H.eta(eta = eta, jbetas = jbetas,  J_X = J_X, J_Y = J_Y, J_Z = J_Z, G.inv = G.inv, n_eta = n_eta,R_Y=R_Y,R_X=R_X,R_Z=R_Z,rbetas=rbetas)
    #For H, solve is much faster than chol2inv
    var.eta <- as.matrix(solve(H))
    product <- matrix(rep(1, n_eta), n_eta, 1)
    gr <- rep(1, n_eta)
    while (as.numeric(gr %*% product) > 1e-08) {
        gr <- gr.eta(eta = eta, J_X = J_X, J_Y = J_Y, J_Z = J_Z, jbetas = jbetas,  G.inv = G.inv, n_eta = n_eta,R_Z=R_Z,R_Y=R_Y,R_X=R_X,rbetas=rbetas)
        product <- var.eta %*% gr
        eta <- eta - product
        H <- H.eta(eta = eta, jbetas = jbetas, J_X = J_X, J_Y = J_Y, J_Z = J_Z, G.inv = G.inv,  n_eta = n_eta,R_Y=R_Y,R_X=R_X,R_Z=R_Z,rbetas=rbetas)
        var.eta <- as.matrix(solve(H))
    }
    rm(product, gr)
    chol.H <- chol(H)
    log.p.eta <- -(length(eta)/2) * log(2 * pi) - sum(log(diag(G.chol))) - 0.5 * crossprod(eta, as(G.inv,"generalMatrix")) %*% eta
    log.p.j <- sum(jder0(J_X %*% jbetas + J_Z %*% eta,J_Y))
    log.p.r <- sum(pnorm(as.vector((-1)^(1 - R_Y) * (R_X %*% rbetas + R_Z %*% eta)), log.p = TRUE))
    res <- var.eta
    attr(res, "likelihood") <- as.vector(cons.logLik + log.p.eta  + log.p.j +log.p.r - 0.5 * (2 * sum(log(diag(chol.H)))))
    attr(res, "eta") <- eta
    res
}


#-----------------
#M-step update for fixed effects
Sc.jbetas.f <- function(eta, jbetas, J_X, J_Y, J_Z, n_eta, Svar, Svar2) {
    E.jbetas <-  jder1((J_X %*% jbetas + J_Z %*% eta),J_Y)
    R2 <- jder2((J_X %*% jbetas + J_Z %*% eta))
    R3 <- jder3( (J_X %*% jbetas + J_Z %*% eta))
    temp.C <- (R2 * J_Z)
    block <- (Svar %*% t(temp.C)) * R3
    trc.jbetas <- rowSums(Svar2 * R3) + colSums(rowSums(Svar2) * block)
    score.jbetas <- as.vector(t(J_X) %*% (E.jbetas + 0.5 * trc.jbetas))
    rm(temp.C, block)
    score.jbetas
}

Sc.jbetas.f.first.order <- function(eta, jbetas,  J_X, J_Y, J_Z) {
    E.jbetas <- jder1(J_X %*% jbetas + J_Z %*% eta,J_Y)
    score.jbetas <- as.vector(t(J_X) %*% E.jbetas )
    score.jbetas
}


#update jbetas
#update for fixed effects
#requires var.eta, not var.eta.hat
update.jbetas <- function(eta, var.eta, jbetas, J_X, J_Y, J_Z, n_eta) {
    Svar <- as.matrix(J_Z %*% var.eta)
    Svar2 <- Svar * J_Z
    product <- matrix(rep(1, length(jbetas)), length(jbetas), 1)
    scjbetas <- rep(1, length(jbetas))
    while (as.numeric(scjbetas %*% product) > 1e-08) {
        Hjbetas <- as.matrix(jacobian(Sc.jbetas.f, jbetas, method = "simple", eta = eta, J_X = J_X, J_Y = J_Y, J_Z = J_Z, n_eta = n_eta, Svar = Svar, Svar2 = Svar2))
        scjbetas <- Sc.jbetas.f(eta, jbetas, J_X, J_Y, J_Z, n_eta, Svar = Svar, Svar2 = Svar2)
        product <- c(solve(Hjbetas, scjbetas))
        jbetas <- jbetas - product
    }
    rm(Svar, Svar2, scjbetas, product)
    jbetas
}

update.jbetas.first.order <- function(eta, jbetas, J_X, J_Y, J_Z, n_eta) {

    product <- matrix(rep(1, length(jbetas)), length(jbetas), 1)
    scjbetas <- rep(1, length(jbetas))
    while (as.numeric(scjbetas %*% product) > 1e-08) {
        Hjbetas <- as.matrix(jacobian(Sc.jbetas.f.first.order, jbetas, method = "simple", eta = eta, J_X = J_X, J_Y = J_Y, J_Z = J_Z))
        scjbetas <- Sc.jbetas.f.first.order(eta, jbetas, J_X, J_Y, J_Z)
        product <- c(solve(Hjbetas, scjbetas))
        jbetas <- jbetas - product
    }
    rm( scjbetas, product)
    jbetas
    }

 #-----------------
#M-step update for fixed effects in missing data mechanism
Sc.rbetas.f <- function(eta, rbetas, R_X, R_Y, R_Z, n_eta, Svar, Svar2) {
    E.rbetas <- ((-1)^(1 - R_Y) * rder1((-1)^(1 - R_Y) * (R_X %*% rbetas + R_Z %*% eta)))
    R2 <- rder2((-1)^(1 - R_Y) * (R_X %*% rbetas + R_Z %*% eta))
    R3 <- (-1)^(1 - R_Y) * rder3((-1)^(1 - R_Y) * (R_X %*% rbetas + R_Z %*% eta))
    temp.C <- (R2 * R_Z)
    block <- (Svar %*% t(temp.C)) * R3
    trc.rbetas <- rowSums(Svar2 * R3) + colSums(rowSums(Svar2) * block)
    score.rbetas <- as.vector(t(R_X) %*% (E.rbetas + 0.5 * trc.rbetas))
    rm(temp.C, block)
    score.rbetas
}

Sc.rbetas.f.first.order <- function(eta, rbetas,  R_X, R_Y, R_Z) {
    E.rbetas <- ((-1)^(1 - R_Y) * rder1((-1)^(1 - R_Y) * (R_X %*% rbetas + R_Z %*% eta)))
    score.rbetas <- as.vector(t(R_X) %*% E.rbetas )
    score.rbetas
}
#update rbetas
#update for fixed effects in missing data mechansim
#requires var.eta, not var.eta.hat
update.rbetas <- function(eta, var.eta, rbetas, R_X, R_Y, R_Z, n_eta) {
    Svar <- as.matrix(R_Z %*% var.eta)
    Svar2 <- Svar * R_Z
    product <- matrix(rep(1, length(rbetas)), length(rbetas), 1)
    scrbetas <- rep(1, length(rbetas))
    while (as.numeric(scrbetas %*% product) > 1e-08) {
        Hrbetas <- as.matrix(jacobian(Sc.rbetas.f, rbetas, method = "simple", eta = eta, R_X = R_X, R_Y = R_Y, R_Z = R_Z, n_eta = n_eta, Svar = Svar, Svar2 = Svar2))
        scrbetas <- Sc.rbetas.f(eta, rbetas, R_X, R_Y, R_Z, n_eta, Svar = Svar, Svar2 = Svar2)
        product <- c(solve(Hrbetas, scrbetas))
        rbetas <- rbetas - product
    }
    rm(Svar, Svar2, scrbetas, product)
    rbetas
}

update.rbetas.first.order <- function(eta, rbetas, R_X, R_Y, R_Z, n_eta) {
    product <- matrix(rep(1, length(rbetas)), length(rbetas), 1)
    scrbetas <- rep(1, length(rbetas))
    while (as.numeric(scrbetas %*% product) > 1e-08) {
        Hrbetas <- as.matrix(jacobian(Sc.rbetas.f.first.order, rbetas, method = "simple", eta = eta, R_X = R_X, R_Y = R_Y, R_Z = R_Z))
        scrbetas <- Sc.rbetas.f.first.order(eta, rbetas, R_X, R_Y, R_Z)
        product <- c(solve(Hrbetas, scrbetas))
        rbetas <- rbetas - product
    }
    rm( scrbetas, product)
    rbetas
    }
              Trace.calc2 <- function(comp.k,comp.l,dsig.dc.indx5, indx.5) {
                  if (comp.k !=comp.l) {
                    if(comp.k!=indx.5) comp2<-comp.k
                    if(comp.l!=indx.5) comp2<-comp.l
                    q<-var.eta%*%(dsig.dc.indx5%*%var.eta[,comp2])
                    d2.sig<- t_J_Z%*%(J_Z*(as.vector(J_Z%*%q)*temp.trc.C+ Svar[,indx.5]*Svar[,comp2]*temp.trc.D))+
                             t_R_Z%*%(R_Z*(as.vector(R_Z%*%q)*R_temp.trc.C+ R_Svar[,indx.5]*R_Svar[,comp2]*R_temp.trc.D)) 
                    dsig.dc.comp2<- t_J_Z%*%(J_Z*Svar[,comp2]*temp.trc.C)
                    res<-sum(diag(var.eta%*%d2.sig)) + sum((t(dsig.dc.indx5)%*%var.eta)*(var.eta%*%dsig.dc.comp2))
                  }else{

                    q<-var.eta%*%(dsig.dc.indx5%*%var.eta[,indx.5])
                    d2.sig<- t_J_Z%*%(J_Z*(as.vector(J_Z%*%q)*temp.trc.C+ Svar[,indx.5]*Svar[,indx.5]*temp.trc.D)) +
                             t_R_Z%*%(R_Z*(as.vector(R_Z%*%q)*R_temp.trc.C+ R_Svar[,indx.5]*R_Svar[,indx.5]*R_temp.trc.D)) 
                    res<-sum(diag(var.eta%*%d2.sig)) + sum((t(dsig.dc.indx5)%*%var.eta)*(var.eta%*%dsig.dc.indx5))

                  }
                   res

              }

#-------------------------------
#Data format

J_Y<-c(t(cbind(Z_mat$Score.For,Z_mat$Score.Against)))
#number of measured scores
R_mat <- Z_mat
Nj <- length(Z_mat$home_win)
Nr <- length(R_mat$home_win)
R_RE_mat <- Matrix(0,Nr,length(teams))
colnames(R_RE_mat)<-teams
for(i in 1:length(teams)){
R_RE_mat[R_mat$home==teams[i],i]<-rep(1,length(R_RE_mat[R_mat$home==teams[i],i]))
R_RE_mat[R_mat$away==teams[i],i]<-rep(-1,length(R_RE_mat[R_mat$away==teams[i],i]))
}
#populate the X matrix for missing data model
R_X_mat <- sparse.model.matrix(r_fixed_effects, R_mat, drop.unused.levels = TRUE)
#drop 0 columns from X_mat
R_X_mat <- R_X_mat[, !(colSums(abs(R_X_mat)) == 0), drop = FALSE]
R_X_mat <- R_X_mat[, !(colSums(abs(R_X_mat)) == 0), drop = FALSE]
    if (rankMatrix(R_X_mat,method = 'qrLINPACK')[1] != dim(R_X_mat)[2]) {
        stop("WARNING: Fixed-effects design matrix not full-rank")
    }
n_rbeta <- dim(R_X_mat)[2]
R_X <- Matrix(R_X_mat)

R_Y <- as.vector(R_mat$home_win)

if(any(R_Y==2)){
tie.indx<-which(R_Y==2)
X.tie<-R_X[tie.indx,,drop=FALSE]
Z.tie<-R_RE_mat[tie.indx,,drop=FALSE]
R_Y<-R_Y[-tie.indx]
R_X<-R_X[-tie.indx,,drop=FALSE]
R_RE_mat<-R_RE_mat[-tie.indx,,drop=FALSE]
for(i in 1:length(tie.indx)){
R_Y<-c(R_Y,1,0)
R_X<-rBind(R_X,X.tie[i,,drop=FALSE],X.tie[i,,drop=FALSE])
R_RE_mat<-rBind(R_RE_mat,Z.tie[i,,drop=FALSE],Z.tie[i,,drop=FALSE])
}
Nr <- length(R_Y)
}

J_RE_mat <- Matrix(0,2*Nj,2*nteams)
colnames(J_RE_mat)<-rep(teams,each=2)
#offense then defense
  joffense<-c(t(cbind(Z_mat$home,Z_mat$away)))
  jdefense<-c(t(cbind(Z_mat$away,Z_mat$home)))
  J_mat<-cbind(as.numeric(J_Y),joffense,jdefense)
  J_mat<-as.data.frame(J_mat)
  templ<-rep(Z_mat$neutral.site,each=2)
  templ2<-rep("Neutral Site",length(templ))
  for(i in 1:(2*Nj)){
  if(templ[i]==0&i%%2==1)
  templ2[i]<-"Home"
  if(templ[i]==0&i%%2==0)
  templ2[i]<-"Away"
  }
  J_mat<-cbind(J_mat,templ2)
  colnames(J_mat)<-c("J_Y","offense","defense","Location")
     if(control$OT.flag){
   J_mat<-cbind(J_mat,rep(Z_mat$OT,each=2))
   colnames(J_mat)<-c("J_Y","offense","defense","Location","OT")
   }
J_mat<-as.data.frame(J_mat)
J_mat$J_Y<-as.numeric(J_Y)
J_mat$game<-as.factor(rep(1:Nj,each=2))
jreo<- sparse.model.matrix(as.formula(~offense+0),data=J_mat)
jred<- -1*sparse.model.matrix(as.formula(~defense+0),data=J_mat)

J_RE_mat[,seq(1,2*nteams,by=2)]<-jreo
J_RE_mat[,seq(2,2*nteams,by=2)]<-jred
#if(game.effect) J_RE_mat<-cBind(J_RE_mat,sparse.model.matrix(as.formula(~game+0),data=J_mat))
J_X_mat <- sparse.model.matrix(j_fixed_effects, J_mat, drop.unused.levels = TRUE)


if(!game.effect) n_eta <- 3*nteams
if(game.effect) n_eta<- 3*nteams + Nj
n_jbeta <- dim(J_X_mat)[2]


new_jz <- Matrix(0, 2*Nj, 3*nteams)
new_rz <- Matrix(0, Nr, 3*nteams)
j.indicator <- rep(c(1,1,0),length(teams))

if(!game.effect) j_i <- j.indicator * (1:n_eta)
if(game.effect) j_i <- j.indicator * (1:(3*nteams))
j_i <- j_i[j_i != 0]
eta.label <- character(n_eta)
new_jz[, j_i] <- J_RE_mat
eta.label[j_i] <- teams
if(!game.effect) r_i <- (1 - j.indicator) * (1:n_eta)
if(game.effect) r_i<-(1-j.indicator)*(1:(3*nteams))
r_i <- r_i[r_i != 0]
new_rz[, r_i] <- R_RE_mat
#eta.label[r_i] <- c(R_t_effects)


#The results are loaded into the Z matricies
#Z[[]] and J_Z[[]] below -----------------------
FE.count<-0
J_X <- Matrix(J_X_mat)
J_Z <- Matrix(new_jz)
if(game.effect) J_Z<-cBind(J_Z,sparse.model.matrix(as.formula(~game+0),data=J_mat))
t_J_Z <- t(J_Z)
cross_J_Z <- crossprod(J_Z)
R_Z <- Matrix(new_rz)
if(game.effect) R_Z<-cBind(R_Z,Matrix(0,nrow(R_Z),length(unique(J_mat$game))))
t_R_Z <- t(R_Z)
cross_R_Z <- crossprod(R_Z)
#if(game.effect) J_RE_mat<-cBind(J_RE_mat,sparse.model.matrix(as.formula(~game+0),data=J_mat))
#initialize parameters
eta.hat <- trc.y1 <- numeric(n_eta)
var.eta.hat <- Matrix(0, n_eta, n_eta)
G <- Diagonal(n_eta)
cons.logLik <- 0.5 * n_eta * log(2 * pi)
L1.conv <- FALSE
L2.conv <- FALSE
L1.conv.it <- 0
suppressWarnings(jbetas<-log(as.vector(solve(crossprod(J_X))%*%t(J_X)%*%J_Y)))
jbetas[is.na(jbetas)]<-.001
rbetas <- rep(0,ncol(R_X_mat))
iter <- control$iter.EM
j.mat <- Matrix(0, iter, n_jbeta)
r.mat <- Matrix(0, iter, n_rbeta)
time.mat <- Matrix(0, iter, 1)
if(!game.effect) G.mat <- Matrix(0, iter, 6)
if(game.effect) G.mat <- Matrix(0, iter, 7)
lgLik <- numeric(iter)
#Begin EM algorithm
for (it in 1:iter) {
    ptm <- proc.time()
    rm(var.eta.hat)
  
    new.eta <- update.eta(eta = eta.hat, J_X = J_X, J_Y = J_Y, J_Z = J_Z,jbetas = jbetas, G = G, n_eta = n_eta, cons.logLik = cons.logLik, R_Z=R_Z,R_Y=R_Y,R_X=R_X,rbetas=rbetas)
   
    #save parameter values in matrix
    j.mat[it, ] <- c(jbetas)
    r.mat[it, ] <- c(rbetas)
    lgLik[it] <- attr(new.eta, "likelihood")
    trc.y1 <- numeric(n_eta)
    trc.y2 <-Matrix(0, n_eta, n_eta)
    G.mat[it, ] <- reduce.G(G)
    eta <- as.vector(attr(new.eta, "eta"))
    var.eta <- new.eta
    rm(new.eta)
    if(home.field){
     thets1 <- c(r.mat[it - 1,1 ], j.mat[it - 1, ], G.mat[it - 1, ])
    thets2 <- c(r.mat[it,1 ], j.mat[it, ], G.mat[it, ])
    }else{
     thets1 <- c( j.mat[it - 1, ], G.mat[it - 1, ])
    thets2 <- c( j.mat[it, ], G.mat[it, ])    
    }
    # print results if verbose
    if ((control$verbose) & (it > 1)) {
        cat("\n\niter:", it, "\n")
        cat("log-likelihood:", lgLik[it], "\n")
        cat("max % change in parm.", 100 * max(round(abs((thets2 - thets1)/thets1 * as.numeric(abs(thets1) > 0.00001)), 8)), "%\n")
        cat("p.mean:", round(jbetas, 4), "\n")
        cat("b.mean:", round(rbetas, 4), "\n")
        cat("G:", reduce.G(G),"\n")
    }
    #check for convergence of first order Laplace approximation
    if ((it > 10) & (L1.conv == FALSE)) {
        check.parm1 <- max(abs(thets2 - thets1) * as.numeric(abs(thets1) <= 0.00001)) < control$tol1
        check.parm2 <- max(abs((thets2 - thets1)/thets1 * as.numeric(abs(thets1) > 0.00001))) < control$tol1

        if (check.parm1 & check.parm2) {
            L1.conv <- TRUE
            L1.conv.it <- it
            eblup.L1 <- round(cbind(eta, sqrt(diag(var.eta))), 4)
            if (control$verbose) {
                cat("\n First order Laplace has converged \n")
              
                cat("log-likelihood:", lgLik[it], "\n")
                cat("max % change in parm.", 100 * round(max(abs(thets2 - thets1)/abs(thets1) * as.numeric(abs(thets1) > 0.00001)), 5), "%\n")
                cat("p.mean:", round(jbetas, 4), "\n")
                cat("b.mean:", round(rbetas, 4), "\n")
        cat("G:", reduce.G(G),"\n")
            }

#first.order.eblup<-cbind(1:nteams,eblup[order(eblup[,1],decreasing=TRUE),])
            if(first.order) break
        }
    }
    #E-step
    #Fully exponential corrections, calculated after 1st order
    #Laplace approximation converges
           if (L1.conv == TRUE) {
              Svar <- as.matrix(J_Z %*% var.eta)
              R_Svar<- as.matrix(R_Z %*% var.eta)
              temp.trc.D <- as.vector(jder4( (J_X %*% jbetas + J_Z %*% eta)))
          temp.trc.C <-  jder3(J_X %*% jbetas + J_Z %*% eta)
          R_temp.trc.D <- as.vector(rder4((-1)^(1 - R_Y) * (R_X %*% rbetas + R_Z %*% eta)))
          R_temp.trc.C <- (-1)^(1 - R_Y) * rder3((-1)^(1 - R_Y) * (R_X %*% rbetas + R_Z %*% eta))
          trc.y1 <- rep(0,n_eta)
          check.parmFE1 <- max(abs(thets2 - thets1) * as.numeric(abs(thets1) <= 0.0001)) < control$tolFE
          check.parmFE2 <- max(abs((thets2 - thets1)/thets1 * as.numeric(abs(thets1) > 0.0001))) < control$tolFE
          #calculate corrections to var.eta.hat after only calculating corrections to eta.hat for a while
         
              #compute the trace corrections for var.eta

              #The function trace.calc2 is used to calculate the trace corrections

              temp.comp <- comp.list
              length.comp <- dim(comp.list)[1]
              var.corrections <- matrix(0, 0, 3)
              for (indx.5 in n_eta:1) {
                  ptm.cor <- proc.time()[3]
                  guide1 <- unique(which(temp.comp == indx.5, arr.ind = TRUE)[, 1])
                  guide2 <- comp.list[guide1, , drop = FALSE]
                  if (length(guide1) > 0) {
                  if(indx.5==n_eta)cat("Calculating FE corrections for random effects vector...\n")
                    dsig.dc.indx5<- t_J_Z%*%(J_Z*Svar[,indx.5]*temp.trc.C)+t_R_Z%*%(R_Z*R_Svar[,indx.5]*R_temp.trc.C)
                    trc.y1[indx.5]<- sum(diag(var.eta%*%dsig.dc.indx5))
                     if (check.parmFE1 & check.parmFE2| FE.count > 0) {
              if(indx.5==n_eta) cat("Calculating FE corrections for random effects covariance matrix...\n")
              FE.count <- FE.count + 1
                    res.temp <- cBind(guide2, mapply(Trace.calc2, comp.k = guide2[, 1], comp.l = guide2[, 2], MoreArgs = list(dsig.dc.indx5, indx.5 = indx.5)))
                    temp.comp <- temp.comp[-(guide1), , drop = FALSE]
                    var.corrections <- rBind(var.corrections, res.temp)
                    }
   
                  }
                #  cat(1-dim(temp.comp)[1]/length.comp,'%,',indx.5,',',proc.time()[3]-ptm.cor,'\n')
              }
              rm(indx.5)
               if (check.parmFE1 & check.parmFE2  | FE.count > 0) {
           
              trc.y2 <- sparseMatrix(i = var.corrections[, 1], j = var.corrections[, 2], x = var.corrections[, 3], dims = c(n_eta, n_eta))
              trc.y2 <- trc.y2 + t(trc.y2) - diag(diag(trc.y2))
            }
              rm(Svar,dsig.dc.indx5)
            
      }
    
    
    
      eta.hat <- as.vector(eta + 0.5 * trc.y1)
      flush.console()
      ptm.rbet <- proc.time()[3]
 
    if(first.order==FALSE&L1.conv==TRUE){
  jbetasn <- update.jbetas(eta = eta, var.eta = var.eta, jbetas = jbetas, J_X=J_X, J_Y=J_Y, J_Z=J_Z, n_eta=n_eta)
  rbetasn <- update.rbetas(eta = eta, var.eta = var.eta, rbetas = rbetas, R_X=R_X, R_Y=R_Y, R_Z=R_Z, n_eta=n_eta)
}else{
jbetasn <- update.jbetas.first.order(eta = eta, jbetas = jbetas, J_X=J_X, J_Y=J_Y, J_Z=J_Z, n_eta=n_eta)
rbetasn <- update.rbetas.first.order(eta = eta, rbetas = rbetas, R_X=R_X, R_Y=R_Y, R_Z=R_Z, n_eta=n_eta)
}

 if(!home.field) rbetasn<-numeric(length(rbetasn))

    
      rm(eta)
      var.eta <- var.eta + 0.5 * trc.y2
      rm(trc.y2)
      var.eta.hat <- var.eta
      eblup <- as.matrix(cBind(eta.hat, sqrt(diag(var.eta.hat))))

      colnames(eblup) <- c("eblup", "std. error")
      if(!game.effect) rownames(eblup) <- c(rep(teams,each=3))
      if(game.effect) rownames(eblup) <- c(rep(teams,each=3),1:Nj)
      rm(var.eta)
            if(first.order==TRUE&L1.conv==TRUE){

                  if (control$verbose) {
                  cat("\nFirst-order Model has converged!\n")
                  cat("\n\niter:", it, "\n")
                  cat("log-likelihood:", lgLik[it], "\n")
                  cat("max % change in parm.", 100 * round(max(abs((thets2 - thets1)/thets1) * as.numeric(abs(thets1) > 0.00001)), 5), "%\n")
                  cat("p.mean:", round(jbetas, 4), "\n")
                  cat("r.mean:", round(rbetas, 4), "\n")
                  cat("G:", reduce.G(G),"\n")

              }
              break

        }
      # check convergence of Fully exponential model
      if ((L1.conv == TRUE) & (it > L1.conv.it)) {
          check.parm1 <- max((thets2 - thets1) * as.numeric(abs(thets1) <= 0.00001)) < control$tol2
          check.parm2 <- max(abs((thets2 - thets1)/thets1) * as.numeric(abs(thets1) > 0.00001)) < control$tol2
          if (check.parm1 & check.parm2 ) {
              L2.conv <- TRUE
              if (control$verbose) {
                  cat("\nFully Exponential Model has converged!\n")
                  cat("\n\niter:", it, "\n")
                  cat("log-likelihood:", lgLik[it], "\n")
                  cat("max % change in parm.", 100 * round(max(abs((thets2 - thets1)/thets1) * as.numeric(abs(thets1) > 0.00001)), 5), "%\n")
                  cat("p.mean", round(jbetas, 4), "\n")
                  cat("r.mean:", round(rbetas, 4), "\n")
                         cat("G:", reduce.G(G),"\n")
              }
              break
          }
      }
      # M-step
      #The following steps update the G
      # matrix

      temp_mat <- var.eta.hat + tcrossprod(eta.hat, eta.hat)
      if(!game.effect){
      gt1<-gt2<-matrix(0,3,3)
      for(i in 1:nteams){
       gt1<-gt1+temp_mat[(3*(i-1)+1):(3*i),(3*(i-1)+1):(3*i)]
      }
      gt1<-gt1/nteams
      Gn<-suppressMessages(kronecker(Diagonal(nteams),symmpart(gt1)))
      }else{
       gt1<-gt2<-matrix(0,3,3)
      for(i in 1:nteams){
       gt1<-gt1+temp_mat[(3*(i-1)+1):(3*i),(3*(i-1)+1):(3*i)]
      }
      gt1<-gt1/nteams
      Gn<-suppressMessages(kronecker(Diagonal(nteams),symmpart(gt1)))
      Gn<-bdiag(Gn,Diagonal(Nj)*mean(diag(temp_mat)[(3*nteams+1):n_eta]))
      }
      jbetas <- jbetasn
      rbetas <- rbetasn
      rbetas[colnames(R_X)=="neutral.site"]<-0
      names(jbetas)<-colnames(J_X)
      G <- Gn
      if(it==4){
      comp.list <- which(abs(G) > 1e-10, arr.ind = TRUE)
      comp.list <- comp.list[comp.list[,1]<=comp.list[,2],]
      comp.list <- unique(comp.list)
}
      rm(Gn)
        it.time <- (proc.time() - ptm)[3]
      time.mat[it, ] <- c(it.time)
      cat("Iteration", it, "took", it.time, "\n")
      eblup <- cBind(eta.hat, sqrt(diag(var.eta.hat)))

      colnames(eblup) <- c("eblup", "std. error")
      if(!game.effect) rownames(eblup) <- c(rep(teams,each=3))
      if(game.effect) rownames(eblup) <- c(rep(teams,each=3),1:Nj)
      #residual<-as.vector(J_Y-exp(J_X%*%jbetas+J_Z%*%eta.hat))
}  #end EM



if(!game.effect){
G.res<-as.matrix(G[1:3,1:3])
colnames(G.res)<-c("Offense","Defense","Win Propensity")
}
if(game.effect){
G.res<-as.matrix(bdiag(G[1:3,1:3],G[n_eta,n_eta]))
colnames(G.res)<-c("Offense","Defense","Win Propensity","Game Effect")
}
G.res.cor<-cov2cor(G.res)

names(jbetas)<-colnames(J_X_mat)
names(rbetas)[1]<-c("LocationHome")



Score <- function(thetas) {
n_ybeta<-length(jbetas)
Ny<-length(J_Y)
    ybetas <- thetas[1:n_ybeta]
    if(home.field){
    rbetas<-c(thetas[n_ybeta+1])
    rbetas<-c(rbetas,rep(0,ncol(R_X)-1))
    G <- thetas[(n_ybeta+2):length(thetas)]
    }else{
    G <- thetas[(n_ybeta+1):length(thetas)]    
    rep(0,ncol(R_X))
    }
    G<-suppressMessages(kronecker(Diagonal(length(teams)),ltriangle(G)))
 new.eta <- update.eta(eta = eta.hat, J_X = J_X, J_Y = J_Y, J_Z = J_Z,jbetas = ybetas, G = G, n_eta = n_eta, cons.logLik = cons.logLik, R_Z=R_Z,R_Y=R_Y,R_X=R_X,rbetas=rbetas)
   
    eta <- attr(new.eta, "eta")
    eta.hat<-eta
    var.eta <- var.eta.hat <- new.eta
    eta.hat <- as.vector(eta)
    temp_mat <- var.eta.hat + tcrossprod(eta.hat, eta.hat)
   # temp_mat_R <- attr(new.eta, "h.inv") + tcrossprod(eta.hat,
    #        eta.hat)
    rm(new.eta)
    if(home.field){
      score.r <- Sc.rbetas.f.first.order(eta, rbetas, R_X, R_Y, R_Z)
    }
      score.y <- Sc.jbetas.f.first.order(eta, ybetas, J_X, J_Y, J_Z)
    
     gam_t_sc <- list()
        index1 <- 0
        score.G <- Matrix(0, 0, 0)
         gam_t_sc <- matrix(0, 3,3)
         index2 <- c(1)
         for (k in 1:nteams) {
                gam_t_sc <- gam_t_sc + temp_mat[(index2):(index2 + 
                  2), (index2):(index2 + 2)]
                index2 <- index2 + 3
            }
            gam_t <- G[1:3, 1:3]
            sv_gam_t <- chol2inv(chol(gam_t))
        der <- -0.5 * (nteams * sv_gam_t - sv_gam_t %*% 
                gam_t_sc %*% sv_gam_t)
            if (is.numeric(drop(sv_gam_t))) {
                score.eta.t <- der
            }
            else {
                score.eta.t <- 2 * der - diag(diag(der))
            }
            
           # for (k in 1:nteams) {
           #     score.G <- bdiag(score.G, score.eta.t)
           # }
        
      score.G<-ltriangle(score.eta.t)  
   

     if(home.field){      
    -c(score.y, score.r[1], score.G)
    }else{
     -c(score.y, score.G)
    }
}

Score.ge <- function(thetas) {
n_ybeta<-length(jbetas)
Ny<-length(J_Y)
    ybetas <- thetas[1:n_ybeta]
    if(home.field){
    rbetas<-c(thetas[n_ybeta+1])
    rbetas<-c(rbetas,rep(0,ncol(R_X)-1))
    G <- thetas[(n_ybeta+2):length(thetas)]
    }else{
    G <- thetas[(n_ybeta+1):length(thetas)]    
    rep(0,ncol(R_X))
    }
      G<-bdiag(suppressMessages(kronecker(Diagonal(length(teams)),bdiag(ltriangle(G[1:6])))),Diagonal(Nj)*G[7])
 new.eta <- update.eta(eta = eta.hat, J_X = J_X, J_Y = J_Y, J_Z = J_Z,jbetas = ybetas, G = G, n_eta = n_eta, cons.logLik = cons.logLik, R_Z=R_Z,R_Y=R_Y,R_X=R_X,rbetas=rbetas)
   
    eta <- attr(new.eta, "eta")
    eta.hat<-eta
    var.eta <- var.eta.hat <- new.eta
    eta.hat <- as.vector(eta)
    temp_mat <- var.eta.hat + tcrossprod(eta.hat, eta.hat)
   # temp_mat_R <- attr(new.eta, "h.inv") + tcrossprod(eta.hat,
    #        eta.hat)
    rm(new.eta)
    if(home.field){
      score.r <- Sc.rbetas.f.first.order(eta, rbetas, R_X, R_Y, R_Z)
    }
      score.y <- Sc.jbetas.f.first.order(eta, ybetas, J_X, J_Y, J_Z)
    
         gam_t_sc <- list()
        index1 <- 0
        score.G <- Matrix(0, 0, 0)
         gam_t_sc <- matrix(0, 3,3)
         index2 <- c(1)
         for (k in 1:nteams) {
                gam_t_sc <- gam_t_sc + temp_mat[(index2):(index2 + 
                  2), (index2):(index2 + 2)]
                index2 <- index2 + 3
            }
         #gam_t_sc<-bdiag(gam_t_sc,))   
            gam_t <- G[1:3, 1:3]
            sv_gam_t <- chol2inv(chol(gam_t))
        der <- -0.5 * (nteams * sv_gam_t - sv_gam_t %*% 
                gam_t_sc %*% sv_gam_t)
            if (is.numeric(drop(sv_gam_t))) {
                score.eta.t <- der
            }
            else {
                score.eta.t <- 2 * der - diag(diag(der))
            }
          der.g <- as.numeric(-0.5 * (Nj * solve(G[n_eta,n_eta]) - solve(G[n_eta,n_eta]) * 
                sum(diag(temp_mat)[(3*nteams+1):n_eta] * solve(G[n_eta,n_eta]))))
  
        score.G<-c(ltriangle(score.eta.t),der.g)   
   

     if(home.field){      
    -c(score.y, score.r[1], score.G)
    }else{
     -c(score.y, score.G)
    }
}

if(home.field){
thetas <- c(jbetas,rbetas[1], reduce.G(G))
names(thetas)<-c(colnames(J_X),"Binary mean","G[1,1]","G[2,1]","G[3,1]","G[2,2]","G[3,2]","G[3,3]")
if(game.effect)  names(thetas)<-c(colnames(J_X),"Binary mean","G[1,1]","G[2,1]","G[3,1]","G[2,2]","G[3,2]","G[3,3]","G[4,4]")
}else{
thetas <- c(jbetas, reduce.G(G))
names(thetas)<-c("Mean Score","G[1,1]","G[2,1]","G[3,1]","G[2,2]","G[3,2]","G[3,3]")
if(game.effect)  names(thetas)<-c("Mean Score","G[1,1]","G[2,1]","G[3,1]","G[2,2]","G[3,2]","G[3,3]","G[4,4]")
}

Hessian<-NULL
if(control$Hessian){
cat("\nCalculating Hessian with a central difference approximation...\n")
flush.console()
if(game.effect){
Hessian <- symmpart(jacobian(Score.ge, thetas, method="simple"))
rownames(Hessian)<-colnames(Hessian)<-names(thetas)
}else{
Hessian <- symmpart(jacobian(Score, thetas, method="simple"))
rownames(Hessian)<-colnames(Hessian)<-names(thetas)
}
#std_errors <- c(sqrt(diag(solve(Hessian))))
if(class(try(chol(Hessian),silent=TRUE))=="try-error") cat("\nWarning: Hessian not positive-definite\n")
}




   res<-list(n.ratings.mov=NULL,n.ratings.offense=NULL,n.ratings.defense=NULL,p.ratings.offense=eblup[seq(1,3*nteams,by=3),1],p.ratings.defense=eblup[seq(2,3*nteams,by=3),1],b.ratings=eblup[seq(3,3*nteams,by=3),1],n.mean=NULL,p.mean=jbetas,b.mean=rbetas[1],G=G.res,G.cor=G.res.cor,R=NULL,R.cor=NULL,home.field=home.field,Hessian=Hessian,parameters=thetas)
}
