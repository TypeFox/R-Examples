normal_cre <-
function(Z_mat=Z_mat, first.order=first.order,home.field=home.field, control=control){

if(home.field&!control$OT.flag){
y_fixed_effects <- formula(~Location+0)
}else if(home.field&control$OT.flag){
y_fixed_effects <- formula(~Location+(OT)+0)   
}else if(!home.field&control$OT.flag){
y_fixed_effects <- formula(~(OT)+0) 
}else{
y_fixed_effects <- formula(~1)
}


#cutoff<-100
home_field<-home.field


#football_EM_twovar_cont<-function(football_data,home_field=FALSE, fbs_list=fbs_list_o,control=control,cutoff=100){
#Z_mat<-football_data

X<-NULL
      H.eta <- function(sigmas, cross_Z_j, Sig.mat, G.inv, nyear, 
          n_eta,sigmas2, cross_R_Z_j, Sig.mat2,R_R.inv) {
          h.eta <- G.inv
          

        h.r<-crossprod(R_Z, R_R.inv) %*% R_Z
        
          (h.eta + h.r)
      }
      ltriangle <- function(x) {
          if (!is.null(dim(x)[2])) {
              resA <- as.vector(x[lower.tri(x, diag = TRUE)])
              resA
          }
          else {
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

 reduce.G<-function(G) ltriangle(as.matrix(G[1:2,1:2]))
#-------------------
   update.eta <- function(X, Y, Z, cross_Z_j, Sig.mat, 
          ybetas, sigmas, G, nyear, n_eta, cons.logLik,R_X, R_Y, R_Z, 
              cross_R_Z_j, Sig.mat2, ybetas2, 
              sigmas2,R_R.inv) {
          G.chol <- chol(G)
          G.inv <- chol2inv(G.chol)
          H <- H.eta(sigmas = sigmas, cross_Z_j = cross_Z_j, Sig.mat = Sig.mat, 
              G.inv = G.inv, nyear = nyear, n_eta = n_eta, sigmas2=sigmas2,cross_R_Z_j=cross_R_Z_j,Sig.mat2=Sig.mat2,R_R.inv)
          chol.H <- chol(H)
          var.eta <- as.matrix(solve(H))
          rm(H)
          eta<-var.eta%*%t(R_Z)%*%R_R.inv%*%(R_Y-R_X%*%ybetas2)
          log.p.eta <- -(length(eta)/2) * log(2 * pi) - sum(log(diag(G.chol))) - 
              0.5 * crossprod(eta, as(G.inv,"generalMatrix")) %*% eta
          #log.p.y <- sum(dnorm(Y, as.vector(X %*% ybetas + Z %*% 
           #   eta), as.vector(Sig.mat %*% sigmas), log = TRUE))
          log.p.r <- -(Nr/2) * log(2 * pi) + sum(log(diag(chol(R_R.inv)))) - 
            0.5 * crossprod(R_Y - R_X %*% ybetas2 - R_Z %*% eta, R_R.inv) %*% 
                (R_Y - R_X %*% ybetas2 - R_Z %*% eta)
          res <- var.eta
          attr(res, "likelihood") <- as.vector(cons.logLik + log.p.eta + 
              log.p.r - 0.5 * (2 * sum(log(diag(chol.H)))))
          attr(res, "eta") <- eta
          res
      }
  
    update.ybeta <- function(X, Y, Z, R_inv, eta.hat) {
        A.ybeta <- crossprod(X, R_inv) %*% X
        B.ybeta <- crossprod(X, R_inv) %*% (Y - Z %*% eta.hat)
        as.vector(solve(A.ybeta, B.ybeta))
    }
#Data format

Z_mat$home<-as.character(Z_mat$home)
Z_mat$away<-as.character(Z_mat$away)
Z_mat$year<-rep(1,dim(Z_mat)[1])
teams <- sort(unique(c(Z_mat$home,Z_mat$away)))
nteams<-length(teams)
teamsfbs<-teams
nfbs<-length(teamsfbs)


J_Y<-c(t(cbind(Z_mat$Score.For,Z_mat$Score.Against)))
#number of measured scores
Nr <- length(Z_mat$home_win)

R_RE_mat <- Matrix(0,Nr,length(teams))
J_RE_mat <- Matrix(0,2*Nr,2*length(teams))
colnames(R_RE_mat)<-teams
colnames(J_RE_mat)<-rep(teams,each=2)
for(i in 1:length(teams)){
R_RE_mat[Z_mat$home==teams[i],i]<-rep(1,length(R_RE_mat[Z_mat$home==teams[i],i]))
R_RE_mat[Z_mat$away==teams[i],i]<-rep(-1,length(R_RE_mat[Z_mat$away==teams[i],i]))
}

#offense then defense
  joffense<-c(t(cbind(Z_mat$home,Z_mat$away)))
  jdefense<-c(t(cbind(Z_mat$away,Z_mat$home)))
  J_mat<-cbind(as.numeric(J_Y),joffense,jdefense)
  J_mat<-as.data.frame(J_mat)
  templ<-rep(Z_mat$neutral.site,each=2)
  templ2<-rep("Neutral Site",length(templ))
  for(i in 1:(2*Nr)){
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
#J_mat$fcs<-as.factor(J_mat$fcs)
J_mat$J_Y<-as.numeric(J_Y)

jreo<-sparse.model.matrix(as.formula(~offense+0),data=J_mat)
jred<--1*sparse.model.matrix(as.formula(~defense+0),data=J_mat)

J_RE_mat[,seq(1,2*length(teams),by=2)]<-jreo
J_RE_mat[,seq(2,2*length(teams),by=2)]<-jred




J_X_mat <- sparse.model.matrix(y_fixed_effects, J_mat, drop.unused.levels = TRUE)

#populate the X matrix for missing data model
#R_X_mat <- sparse.model.matrix(r_fixed_effects, Z_mat, drop.unused.levels = TRUE)
#drop 0 columns from X_mat
#R_X_mat <- R_X_mat[, !(colSums(abs(R_X_mat)) == 0), drop = FALSE]
#if (rankMatrix(R_X_mat)[1] != dim(R_X_mat)[2]) {
#    cat("WARNING: Fixed-effects design matrix for missing-data model not full-rank", "\n")
#    break
#}
n_eta <- 2*length(teams)
#n_ybeta<-n_rbeta <- dim(R_X_mat)[2]
n_jbeta <- dim(J_X_mat)[2]
#X_mat<-R_X_mat

Sig.mat <- as.matrix(rep(1,Nr))
Sig.mat2 <- as.matrix(rep(1,2*Nr))
nyear<-1
#RE_mat<-R_RE_mat
#The results are loaded into the Z matricies
#Z[[]] and R_Z[[]] below -----------------------
FE.count<-0

R_X <- Matrix(J_X_mat)
R_Y <- as.numeric(as.vector(J_mat$J_Y))
#R_Y<-pmin(abs(R_Y),cutoff)*sign(R_Y)
#Z <- Matrix(new_yz)
Z<-c(NULL)
R_Z <- Matrix(J_RE_mat)
t_R_Z <- t(R_Z)


              #from here on, R means J
#X <- Matrix(R_X_mat)
Y <- as.vector(Z_mat$Score.For-Z_mat$Score.Against)
#Y<-pmin(abs(Y),cutoff)*sign(Y)

#t_Z <- t(Z)
#cross_Z <- crossprod(Z)
cross_R_Z <- crossprod(R_Z)
cross_Z_j <- list()
cross_R_Z_j <- list()
      X_j <- list(NULL)
      R_X_j <- list(NULL)
      cross_X_j <- list(NULL)
      cross_R_X_j <- list(NULL)
      Y_j <- list(NULL)
      R_Y_j <- list(NULL)
      Z_j <- list(NULL)
      R_Z_j <- list(NULL)
      for (j in 1:nyear) {
#          cross_Z_j[[j]] <- crossprod(Matrix(new_yz[Z_mat$year == 
 #             j, ]))
              cross_R_Z_j[[j]] <- crossprod(Matrix(J_RE_mat[Z_mat$year == 
              j, ]))
 #         X_j[[j]] <- X_mat[Z_mat$year == j, ]
          R_X_j[[j]] <- J_X_mat
 #         Y_j[[j]] <- as.vector(Y[Z_mat$year == j ])
          R_Y_j[[j]] <- as.vector(J_Y[Z_mat$year == j ])
 #         Z_j[[j]] <- new_yz[Z_mat$year == j, ]
          R_Z_j[[j]] <- J_RE_mat[Z_mat$year == j, ]
  #        cross_X_j[[j]] <- crossprod(X_j[[j]])
           cross_R_X_j[[j]] <- crossprod(R_X_j[[j]])
      }
#initialize parameters
eta.hat <- numeric(n_eta)
var.eta.hat <- Matrix(0, n_eta, n_eta)
G <- 100*Diagonal(n_eta)
R_R<-R_R.inv<-Diagonal(nrow(J_mat))
cons.logLik <- 0.5 * n_eta * log(2 * pi)
#Partition desigin matrices by year and
#calculate initial parameter values
# from data
         sigmas <- c(rep(0, nyear))
         sigmas2 <- c(rep(0, nyear))
ybetas<-0
ybetas2 <- update.ybeta(X=R_X, Y=R_Y, Z=R_Z, R_inv=R_R.inv, eta.hat=eta.hat)
names(ybetas2)<-colnames(J_X_mat)
#these next few lines are used to populate
#comp.list, a list of components needed for
#the E-step update (not all n_eta^2 are needed)

year.count<-Nr
iter <- control$iter.EM
r.mat <- Matrix(0, iter, length(ybetas2))
time.mat <- Matrix(0, iter, 1)
G.mat <- Matrix(0, iter, 3)
lgLik <- numeric(iter)
L1.conv <- FALSE
L2.conv <- FALSE
L1.conv.it <- 0
#Begin EM algorithm
for (it in 1:iter) {
    ptm <- proc.time()
    rm(var.eta.hat)

    new.eta <- update.eta(X = X, Y = Y, Z = Z, 
              cross_Z_j = cross_Z_j, Sig.mat = Sig.mat, ybetas = ybetas, 
              sigmas = sigmas, G = G, nyear = nyear, n_eta = n_eta, 
              cons.logLik = cons.logLik,R_X = R_X, R_Y = R_Y, R_Z = R_Z, 
              cross_R_Z_j = cross_R_Z_j, Sig.mat2 = Sig.mat2, ybetas2 = ybetas2, 
              sigmas2 = sigmas2,R_R.inv=R_R.inv)

    #save parameter values in matrix
    r.mat[it, ] <- c(ybetas2)
    lgLik[it] <- attr(new.eta, "likelihood")
    trc.y1 <- numeric(n_eta)
    trc.y2 <- Matrix(0, n_eta, n_eta)
    G.mat[it, ] <- reduce.G(G)
    eta <- as.vector(attr(new.eta, "eta"))
    var.eta <- new.eta
        eta.hat <- as.vector(eta)
        var.eta.hat <- var.eta
    rm(new.eta)
    thets1 <- c(r.mat[it - 1, ], G.mat[it - 1, ])
    thets2 <- c(r.mat[it, ], G.mat[it, ])
    # print results if verbose
    if ((control$verbose) & (it > 1)) {
        cat("\n\niter:", it, "\n")
        cat("log-likelihood:", lgLik[it], "\n")
                    cat("change in loglik:", sprintf("%.7f", lgLik[it] - 
                      lgLik[it - 1]), "\n")
        cat("n.mean:", round(ybetas2, 4), "\n")
        cat("G:", reduce.G(G),"\n")
  
    }
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
                    cat("n.mean:", round(ybetas2, 4), "\n")
                  
                             cat("G:", reduce.G(G),"\n")
                      flush.console()
                    }
                    rm(j)
                     break
                  }
                 
              }

          


          

    #Fully exponential corrections, calculated after 1st order
    #Laplace approximation converges
  

    flush.console()
    ptm.rbet <- proc.time()[3]
   # cat("starting update.ybetas", "\n")

    #cat("finished update.ybetas", proc.time()[3] - ptm.rbet, "\n")
    rm(eta)
    
  

    eblup <- as.matrix(cBind(eta.hat, sqrt(diag(var.eta.hat))))
    colnames(eblup) <- c("eblup", "std. error")
    rownames(eblup) <- rep(teams,each=2)
    rm(var.eta)

    # M-step
    #The following steps update the G
    # matrix
    
    temp_mat <- var.eta.hat + tcrossprod(eta.hat, eta.hat)
    gt1<-gt2<-matrix(0,2,2)
    for(i in 1:length(teams)){
     gt1<-gt1+temp_mat[(2*(i-1)+1):(2*i),(2*(i-1)+1):(2*i)]
    }
    gt1<-gt1/nfbs

    Gn<-kronecker(Diagonal(nfbs),symmpart(gt1))
                        sigup<-matrix(0,2,2)
for(i in 1:(nrow(J_mat)/2)){
yb<-R_Y[(2*i-1):(2*i)]
xb<-R_X[(2*i-1):(2*i),,drop=FALSE]
zb<-R_Z[(2*i-1):(2*i),]
if(home.field){
yxb<-yb-xb%*%ybetas2
}else{
yxb<-as.matrix(yb-rep(ybetas2,2))
}
sigup<- suppressWarnings(sigup+suppressMessages(tcrossprod(yxb))-yxb%*%t(zb%*%eta.hat)- (zb%*%eta.hat)%*%t(yxb)+zb%*%temp_mat%*%t(zb))
}
sigup<-symmpart(sigup/(nrow(J_mat)/2))
ybetas2 <- update.ybeta(X=R_X, Y=R_Y, Z=R_Z, R_inv=R_R.inv, eta.hat=eta.hat)

    #if(home_field) ybetas2[colnames(R_X)=="LocationNeutral Site"]<-0

    R_R<-suppressMessages(kronecker(Diagonal(nrow(J_mat)/2),sigup))
    R_R.inv<-suppressMessages(kronecker(Diagonal(nrow(J_mat)/2),solve(sigup)))
    G <- Gn
    rm(Gn)
    it.time <- (proc.time() - ptm)[3]
    time.mat[it, ] <- c(it.time)
    cat("Iteration", it, "took", it.time, "\n")
    eblup <- cBind(eta.hat, sqrt(diag(var.eta.hat)))
    colnames(eblup) <- c("eblup", "std. error")
    rownames(eblup) <- rep(teams,each=2)
}  #end EM




               pattern.f.score <- function(R.i.parm,ybetas,X,Y,Z,Ny) {
        R_i <- ltriangle(as.vector(R.i.parm))
         pattern.Rtemplate <- ltriangle(1:(2/2 * (2 + 1)))
    pattern.diag <- diag(pattern.Rtemplate)
        pattern.score <- numeric(2/2 * (2 + 1))
        

            pattern.sum <- matrix(0, 2, 2)
            for (i in 1:(Ny/2)) {
                X.t <- X[(1 + (i - 1) * 2):(i * 
                  2), , drop = FALSE]
                Y.t <- Y[(1 + (i - 1) * 2):(i * 
                  2)]
                Z.t <- Z[(1 + (i - 1) * 2):(i * 
                  2), , drop = FALSE]      
                temp.t <- Y.t - X.t %*% ybetas
                pattern.sum <- pattern.sum + tcrossprod(temp.t) - 
                  tcrossprod(temp.t, Z.t %*% eta.hat) - tcrossprod(Z.t %*% 
                  eta.hat, temp.t) + as.matrix(Z.t%*%temp_mat%*%t(Z.t))
            }
         pattern.y <- solve(R_i)
         pattern.score<- -ltriangle(((Ny/2) * pattern.y) -(pattern.y %*% pattern.sum %*% pattern.y))
            
         pattern.score<-pattern.score*c(1,.5,1)
        -pattern.score
    }
    

Score <- function(thetas) {
n_ybeta<-length(ybetas2)
Ny<-length(R_Y)
    ybetas <- thetas[1:n_ybeta]
    R.tri<-R.i.parm <-  thetas[(n_ybeta+1):(n_ybeta+3)]
    LRI <- length(R.i.parm)
    #R_i.parm.constrained <- R.i.parm[-LRI]
    R_i<-ltriangle(R.tri)
                R <- symmpart(suppressMessages(kronecker(suppressMessages(Diagonal(Ny/2)),
                R_i)))
            R_inv <- symmpart(suppressMessages(kronecker(suppressMessages(Diagonal(Ny/2)),
                chol2inv(chol(R_i)))))
    G <- thetas[(n_ybeta+4):length(thetas)]
    G<-kronecker(Diagonal(length(teams)),ltriangle(G))
    #update.eta returns new var.eta with eta and likelihood as attr()
      new.eta <- update.eta(X = X, Y = Y, Z = Z, 
              cross_Z_j = cross_Z_j, Sig.mat = Sig.mat, ybetas = ybetas, 
              sigmas = sigmas, G = G, nyear = nyear, n_eta = n_eta, 
              cons.logLik = cons.logLik,R_X = R_X, R_Y = R_Y, R_Z = R_Z, 
              cross_R_Z_j = cross_R_Z_j, Sig.mat2 = Sig.mat2, ybetas2 = ybetas, 
              sigmas2 = sigmas2,R_R.inv=R_R.inv)
    eta <- attr(new.eta, "eta")
    eta.hat<-eta
    var.eta <- var.eta.hat <- new.eta
    eta.hat <- as.vector(eta)
    temp_mat <- var.eta.hat + tcrossprod(eta.hat, eta.hat)
   # temp_mat_R <- attr(new.eta, "h.inv") + tcrossprod(eta.hat,
    #        eta.hat)
    rm(new.eta)
        A.ybeta <- crossprod(R_X, R_inv) %*% R_X
        B.ybeta <- crossprod(R_X, R_inv) %*% (R_Y - R_Z %*% eta.hat)
    score.y <- as.vector(B.ybeta - A.ybeta %*% ybetas)
    score.R <- -pattern.f.score(R.i.parm,ybetas2,R_X,R_Y,R_Z,Ny)
    
     gam_t_sc <- list()
        index1 <- 0
        score.G <- Matrix(0, 0, 0)
         gam_t_sc <- matrix(0, 2,2)
         index2 <- c(1)
         for (k in 1:nteams) {
                gam_t_sc <- gam_t_sc + temp_mat[(index2):(index2 + 
                  1), (index2):(index2 + 1)]
                index2 <- index2 + 2
            }
            gam_t <- G[1:2, 1:2]
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
   

           
    -c(score.y, score.R, score.G)
}




Hessian<-NULL
thetas <- c(ybetas2, ltriangle(sigup), reduce.G(G))
if(control$Hessian){
cat("\nCalculating Hessian with a central difference approximation...\n")
flush.console()

Hessian <- symmpart(jacobian(Score, thetas, method="simple"))
#std_errors <- c(sqrt(diag(solve(Hessian))))
if(class(try(chol(Hessian),silent=TRUE))=="try-error") cat("\nWarning: Hessian not positive-definite\n")
}




G.res<-as.matrix(G[1:2,1:2])
colnames(G.res)<-c("Offense","Defense")
G.res.cor<-cov2cor(G.res)
R.res<-as.matrix(R_R[1:2,1:2])
colnames(R.res)<-c("Home","Away")
R.res.cor<-cov2cor(R.res)
if(!home.field) ybetas2<-ybetas2[1]
names(ybetas2)<-colnames(J_X_mat)

sresid=NULL
cresid=NULL
    mresid <- try(as.numeric(R_Y - R_X %*% ybetas2))
    cresid <- try(as.numeric(mresid - R_Z %*% eta.hat))
    yhat <- try(as.numeric(R_X %*% ybetas2 + R_Z %*% eta.hat))
    rchol <- try(chol(R_R.inv))
    yhat.s <- try(as.vector(rchol %*% (yhat)))
    sresid <- try(as.vector(rchol %*% R_Y - yhat.s))
                                                                                                                                                                                                                                                                                                                           
   res<-list(n.ratings.mov=NULL,n.ratings.offense=eblup[seq(1,2*nteams,by=2),1],n.ratings.defense=eblup[seq(2,2*nteams,by=2),1],p.ratings.offense=NULL,p.ratings.defense=NULL,b.ratings=NULL,n.mean=ybetas2,p.mean=NULL,b.mean=NULL,G=G.res,G.cor=G.res.cor,R=R.res,R.cor=R.res.cor,home.field=home.field,actual=R_Y,pred=R_X%*%ybetas2+R_Z%*%eta.hat,Hessian=Hessian,parameters=thetas,sresid=sresid)
}
