coefH <- function (X, se = TRUE, nice.output = TRUE, group.var = NULL){ # , ACM = FALSE
    X <- check.data(X)
    eps <- 1e-40
    labels <- dimnames(X)[[2]]
    OL <- list()
    if (!is.null(group.var) && nrow(as.matrix(group.var)) != nrow(X)){
       group.var <- NULL
       warning("group.var not the same length/nrow as X: group.var ignored")
    }
    if (!se && is.null(group.var)) {
       S <- var(X)
       Smax <- var(apply(X, 2, sort))
       Hij <- S/Smax
       diag(S) <- 0
       diag(Smax) <- 0
       Hi <- apply(S, 1, sum)/apply(Smax, 1, sum)
       H <- sum(S)/sum(Smax)
       OL <- list(Hij = Hij, Hi = Hi, H = H)
    } else {
       g <- max(X) - min(X) + 1
       J <- ncol(X)
       P <- choose(J, 2)
       N <- nrow(X)
       if (any(apply(X, 2, var) < eps)) stop("One or more variables have zero variance")
       n <- as.matrix(table(apply(X, 1, paste, collapse = "")))
       lab.n <- matrix(names(table(apply(X, 1, paste, collapse = ""))))
       lab.b <- apply(all.patterns(2, g), 2, paste, collapse = "")
       lab.u <- as.character(0:(g - 1))
       Bi <- substr(lab.b, 1, 1)
       Bj <- substr(lab.b, 2, 2)
       R <- t(apply(lab.n, 1, string2integer))
       r <- length(lab.n)
       U <- list()
       for (j in 1:J) U[[j]] <- tabulate(X[,j]+1,g)
       W <- list()
       WA <- list()
       WY <- list()
       WE <- list()
       WF <- list()
       for (i in 1:J) {
          W[[i]] <- list()
          WA[[i]] <- list()
          WY[[i]] <- list()
          WE[[i]] <- list()
          WF[[i]] <- list()
          for (j in i:J) if (j > i) {
             W[[i]][[j]] <- weights(X[, c(i, j)],g-1)
             A1a <- NULL
             for (a in 0:(g - 1)) for (b in 0:(g - 1)) A1a <- rbind(A1a, as.numeric(R[, i] == a & R[, j] == b))
             WA[[i]][[j]] <- W[[i]][[j]] %*% A1a
             Eij <- matrix(U[[i]][as.numeric(Bi) + 1], nrow = g^2, ncol = 1) * matrix(U[[j]][as.numeric(Bj) + 1], nrow = g^2, ncol = 1)/N
             Y22 <- cbind(outer(Bi, lab.u, "=="), outer(Bj, lab.u, "==")) * matrix(Eij, nrow = g^2, ncol = 2 * g)
             Ri <- substr(lab.n, i, i)
             Rj <- substr(lab.n, j, j)
             Z2 <- rbind(outer(lab.u, Ri, "==")[, , 1], outer(lab.u, Rj, "==")[, , 1]) * c(1/U[[i]], 1/U[[j]])
             Z2[is.nan(Z2)] <- 1/eps
             YZ2 <- Y22 %*% Z2 - Eij %*% matrix(1/N, 1, r)
             WY[[i]][[j]] <- W[[i]][[j]] %*% YZ2
             Fij <- complete.observed.frequencies(X[, c(i, j)], 2, g)
             WF[[i]][[j]] <- W[[i]][[j]] %*% Fij
             WE[[i]][[j]] <- W[[i]][[j]] %*% Eij
          }
       }
       g3 <- matrix(c(unlist(WF[[1]][[2]]), unlist(WF), unlist(WE)), nrow = 2 * P + 1, byrow = TRUE)
       A4 <- rbind(matrix(c(1, -1, rep(0, (J * (J - 1)) - 1)), 1, (J * (J - 1)) + 1), cbind(matrix(0, P, 1), diag(P), -1 * diag(P)))
       A5 <- cbind(matrix(1, P, 1), -1 * diag(P))
       g4 <- phi(A4, g3, "log")
       g5 <- phi(A5, g4, "exp")
       G3 <- matrix(c(unlist(WA[[1]][[2]]), unlist(WA), unlist(WY)), nrow = 2 * P + 1, byrow = TRUE)
       G4 <- dphi(A4, g3, G3, "log")
       G5ij <- dphi(A5, g4, G4, "exp")
       Hij <- se.Hij <- matrix(0, J, J)
       Hij[lower.tri(Hij)] <- g5
       Hij <- Hij + t(Hij)
       if (P > 1) G3 <- rbind(apply(G3[2:(P + 1), ], 2, sum), apply(G3[2:(P + 1), ], 2, sum), apply(G3[(P + 2):(2 * P + 1), ], 2, sum))
       g3 <- matrix(c(sum(g3[2:(P + 1), ]), sum(g3[2:(P + 1), ]), sum(g3[(P + 2):(2 * P + 1), ])), ncol = 1)
       A4 <- cbind(matrix(1, 2, 1), -1 * diag(2))
       A5 <- matrix(c(1, -1), 1, 2)
       g4 <- phi(A4, g3, "log")
       g5 <- phi(A5, g4, "exp")
       G4 <- dphi(A4, g3, G3, "log")
       G5 <- dphi(A5, g4, G4, "exp")
       H <- g5
       G3 <- matrix(0, 2 * J + 1, r)
       g3 <- matrix(0, 2 * J + 1, 1)
       # Gaat het hier goed met de haakjes??
       for (j in 1:J) for (k in 1:J) if (k > j) {
          g3[j + 1, ] <- g3[j + 1, ] + WF[[j]][[k]]
          g3[J + j + 1, ] <- g3[J + j + 1, ] + WE[[j]][[k]]
          G3[j + 1, ] <- G3[j + 1, ] + WA[[j]][[k]]
          G3[J + j + 1, ] <- G3[J + j + 1, ] + WY[[j]][[k]]
       } else {
          if (k < j) {
             g3[j + 1, ] <- g3[j + 1, ] + WF[[k]][[j]]
             g3[J + j + 1, ] <- g3[J + j + 1, ] + WE[[k]][[j]]
             G3[j + 1, ] <- G3[j + 1, ] + WA[[k]][[j]]
             G3[J + j + 1, ] <- G3[J + j + 1, ] + WY[[k]][[j]]
          }
       }
       g3[1, ] <- g3[2, ]
       G3[1, ] <- G3[2, ]
       A4 <- rbind(matrix(c(1, -1, rep(0, (J * 2) - 1)), 1, (J * 2) + 1), cbind(matrix(0, J, 1), diag(J), -1 * diag(J)))
       A5 <- cbind(matrix(1, J, 1), -1 * diag(J))
       g4 <- phi(A4, g3, "log")
       g5 <- phi(A5, g4, "exp")
       G4 <- dphi(A4, g3, G3, "log")
       G5i <- dphi(A5, g4, G4, "exp")
       Hi <- matrix(g5)
       if (se) {
          ACM.Hij = G5ij %*% (as.numeric(n) * t(G5ij))
          se.Hij[lower.tri(se.Hij)] <- sqrt(diag(ACM.Hij))
          se.Hij <- se.Hij + t(se.Hij)
          dimnames(se.Hij) <- dimnames(Hij) <- list(labels, labels)
          ACM.Hi = G5i %*% (as.numeric(n) * t(G5i))
          se.Hi <- matrix(sqrt(diag(ACM.Hi)))
          dimnames(se.Hi)[[1]] <- dimnames(Hi)[[1]] <- labels
          ACM.H <- G5 %*% (as.numeric(n) * t(G5))
          se.H <- sqrt(diag(ACM.H))
       }
       # OUTPUT
       if(nice.output && se) {
          output.matrix.Hij <- matrix(NA, J, J * 2)
          for (j in 2 * (1:J)) {
             output.matrix.Hij[, j - 1] <- format(paste(" ", formatC(round(Hij[, j/2], 3), digits = 3, format = "f"), " ", sep = ""), width = 7, justify = "right")
             output.matrix.Hij[, j] <- format(paste("(", formatC(round(se.Hij[, j/2], 3), digits = 3, format = "f"), ")", sep = ""), width = 7, justify = "right")
          }
          new.labels <- rep(labels, each = 2)
          new.labels[2 * (1:J)] <- "se"
          dimnames(output.matrix.Hij)[[1]] <- labels
          dimnames(output.matrix.Hij)[[2]] <- new.labels
          output.matrix.Hij[row(output.matrix.Hij) == 0.5 * col(output.matrix.Hij)] <- format("", width = 7, justify = "right")
          output.matrix.Hij[(row(output.matrix.Hij)) == (0.5 * col(output.matrix.Hij) + 0.5)] <- format("", width = 7, justify = "right")
          output.matrix.Hij <- noquote(output.matrix.Hij)
          output.matrix.Hi <- matrix(NA, J, 2)
          output.matrix.Hi[, 1] <- format(formatC(round(Hi, 3), digits = 3, format = "f"), width = 7, justify = "right")
          output.matrix.Hi[, 2] <- format(paste("(", formatC(round(se.Hi, 3), digits = 3, format = "f"), ")", sep = ""), width = 7, justify = "right")
          dimnames(output.matrix.Hi) <- list(labels, c("Item H", "se"))
          output.matrix.Hi <- noquote(output.matrix.Hi)
          output.matrix.H <- matrix(NA, 1, 2)
          output.matrix.H[, 1] <- format(formatC(round(H, 3), digits = 3, format = "f"), width = 7, justify = "right")
          output.matrix.H[, 2] <- format(paste("(", formatC(round(se.H, 3), digits = 3, format = "f"), ")", sep = ""), width = 7, justify = "right")
          dimnames(output.matrix.H) <- list("", c("Scale H", "se"))
          output.matrix.H <- noquote(output.matrix.H)
          OL <- list(Hij = output.matrix.Hij, Hi = output.matrix.Hi, H = output.matrix.H)
       } else {
          if (se) OL <- list(Hij = Hij, se.Hij = se.Hij, Hi = Hi, se.Hi = se.Hi, H = H, se.H = se.H)
          if (!se) OL <- list(Hij = Hij, Hi = Hi, H = H)
       }
        
       if (!is.null(group.var)){
          group.item <- length(OL) + 1
          OL[[group.item]] <- list() 
          names(OL)[[group.item]] <- "Groups"
          group.var <- apply(as.matrix(group.var),1,paste, sep="", collapse="/")
          group.names <- sort(unique(group.var))
          K <- length(group.names)
          for (group in 1:K){
             X. <- X[group.var == group.names[group],]
             if(length(X.)==ncol(X)){
                warning(paste("No scalability coefficients computed for group",group.names[group],". Group contains less than two cases."))
             } else { 
                X. <- check.data(X.)
                OL[[group.item]][[group]] <- list()
                names(OL[[group.item]])[[group]] <- group.names[group]
                N. <- nrow(X.)
                if (any(apply(X., 2, var) < eps)) warning(paste("In group",group.names[group],", some variables have zero variance"))
                n. <- as.matrix(table(apply(X., 1, paste, collapse = "")))
                lab.n. <- matrix(names(table(apply(X., 1, paste, collapse = ""))))
                lab.b. <- apply(all.patterns(2, g), 2, paste, collapse = "")
                Bi. <- substr(lab.b., 1, 1)
                Bj. <- substr(lab.b., 2, 2)
                R. <- t(apply(lab.n., 1, string2integer))
                r. <- length(lab.n.)
                U. <- list()
                for (j in 1:J) U.[[j]] <- tabulate(X.[,j]+1, g)
                WA. <- list()
                WY. <- list()
                WE. <- list()
                WF. <- list()
                for (i in 1:J) {
                   WA.[[i]] <- list()
                   WY.[[i]] <- list()
                   WE.[[i]] <- list()
                   WF.[[i]] <- list()
                   for (j in i:J) if (j > i) {
                      A1a. <- NULL
                      for (a in 0:(g - 1)) for (b in 0:(g - 1)) A1a. <- rbind(A1a., as.numeric(R.[, i] == a & R.[, j] == b))
                      WA.[[i]][[j]] <- W[[i]][[j]] %*% A1a.
                      Eij. <- matrix(U.[[i]][as.numeric(Bi.) + 1], nrow = g^2, ncol = 1) * matrix(U.[[j]][as.numeric(Bj.) +  1], nrow = g^2, ncol = 1)/N.
                      Y22. <- cbind(outer(Bi., lab.u, "=="), outer(Bj., lab.u, "==")) * matrix(Eij., nrow = g^2, ncol = 2 * g)
                      Ri. <- substr(lab.n., i, i)
                      Rj. <- substr(lab.n., j, j)
                      Z2. <- rbind(outer(lab.u, Ri., "==")[, , 1], outer(lab.u, Rj., "==")[, , 1]) * c(1/U.[[i]], 1/U.[[j]])
                      Z2.[is.nan(Z2.)] <- 1/eps
                      YZ2. <- Y22. %*% Z2. - Eij. %*% matrix(1/N., 1, r.)
                      WY.[[i]][[j]] <- W[[i]][[j]] %*% YZ2.
                      Fij. <- complete.observed.frequencies(X.[, c(i, j)], 2, g)
                      WF.[[i]][[j]] <- W[[i]][[j]] %*% Fij.
                      WE.[[i]][[j]] <- W[[i]][[j]] %*% Eij.
                   }  
                }
                g3. <- matrix(c(unlist(WF.[[1]][[2]]), unlist(WF.), unlist(WE.)), nrow = 2 * P + 1, byrow = TRUE)
                A4 <- rbind(matrix(c(1, -1, rep(0, (J * (J - 1)) - 1)), 1, (J * (J - 1)) + 1), cbind(matrix(0, P, 1), diag(P), -1 * diag(P)))
                A5 <- cbind(matrix(1, P, 1), -1 * diag(P))
                g4. <- phi(A4, g3., "log")
                g5. <- phi(A5, g4., "exp")
                G3. <- matrix(c(unlist(WA.[[1]][[2]]), unlist(WA.), unlist(WY.)), nrow = 2 * P + 1, byrow = TRUE)
                G4. <- dphi(A4, g3., G3., "log")
                G5ij. <- dphi(A5, g4., G4., "exp")
                Hij. <- matrix(0, J, J)
                Hij.[lower.tri(Hij.)] <- g5.
                Hij. <- Hij. + t(Hij.)
                if (P > 1) G3. <- rbind(apply(G3.[2:(P + 1), ], 2, sum), apply(G3.[2:(P + 1), ], 2, sum), apply(G3.[(P + 2):(2 * P + 1),], 2, sum))
                g3. <- matrix(c(sum(g3.[2:(P + 1), ]), sum(g3.[2:(P + 1), ]), sum(g3.[(P + 2):(2 * P + 1), ])), ncol = 1)
                A4 <- cbind(matrix(1, 2, 1), -1 * diag(2))
                A5 <- matrix(c(1, -1), 1, 2)
                g4. <- phi(A4, g3., "log")
                g5. <- phi(A5, g4., "exp")
                G4. <- dphi(A4, g3., G3., "log")
                G5. <- dphi(A5, g4., G4., "exp")
                H. <- g5.
                G3. <- matrix(0, 2 * J + 1, r.)
                g3. <- matrix(0, 2 * J + 1, 1)
                for (j in 1:J) for (k in 1:J) if (k > j) {
                   g3.[j + 1, ] <- g3.[j + 1, ] + WF.[[j]][[k]]
                   g3.[J + j + 1, ] <- g3.[J + j + 1, ] + WE.[[j]][[k]]
                   G3.[j + 1, ] <- G3.[j + 1, ] + WA.[[j]][[k]]
                   G3.[J + j + 1, ] <- G3.[J + j + 1, ] + WY.[[j]][[k]]
                } else { 
                   if (k < j) {
                      g3.[j + 1, ] <- g3.[j + 1, ] + WF.[[k]][[j]]
                      g3.[J + j + 1, ] <- g3.[J + j + 1, ] + WE.[[k]][[j]]
                      G3.[j + 1, ] <- G3.[j + 1, ] + WA.[[k]][[j]]
                      G3.[J + j + 1, ] <- G3.[J + j + 1, ] + WY.[[k]][[j]]
                   }
                }
                g3.[1, ] <- g3.[2, ]
                G3.[1, ] <- G3.[2, ]
                A4 <- rbind(matrix(c(1, -1, rep(0, (J * 2) - 1)), 1, (J * 2) + 1), cbind(matrix(0, J, 1), diag(J), -1 * diag(J)))
                A5 <- cbind(matrix(1, J, 1), -1 * diag(J))
                g4. <- phi(A4, g3., "log")
                g5. <- phi(A5, g4., "exp")
                G4. <- dphi(A4, g3., G3., "log")
                G5i. <- dphi(A5, g4., G4., "exp")
                Hi. <- matrix(g5.)
                if (se) { # if (se || ACM)
                   se.Hij. <- matrix(0, J, J)
                   ACM.Hij = G5ij. %*% (as.numeric(n.) * t(G5ij.))
                   se.Hij.[lower.tri(se.Hij.)] <- sqrt(diag(ACM.Hij))
                   se.Hij. <- se.Hij. + t(se.Hij.)
                   dimnames(se.Hij.) <- dimnames(Hij.) <- list(labels, labels)
                   ACM.Hi = G5i. %*% (as.numeric(n.) * t(G5i.))
                   se.Hi. <- matrix(sqrt(diag(ACM.Hi)))
                   dimnames(se.Hi)[[1]] <- dimnames(Hi.)[[1]] <- labels
                   ACM.H = G5. %*% (as.numeric(n.) * t(G5.))
                   se.H. <- sqrt(diag(ACM.H))
                } 
                # OUTPUT  
                if(nice.output && se) {
                   output.matrix.Hij. <- matrix(NA, J, J * 2)
                   for (j in 2 * (1:J)) {
                      output.matrix.Hij.[, j - 1] <- format(paste(" ", formatC(round(Hij.[, j/2], 3), digits = 3, format = "f"), " ", sep = ""), width = 7, justify = "right")
                      output.matrix.Hij.[, j] <- format(paste("(", formatC(round(se.Hij.[, j/2], 3), digits = 3, format = "f"), ")", sep = ""), width = 7, justify = "right")
                   }
                   dimnames(output.matrix.Hij.)[[1]] <- labels
                   dimnames(output.matrix.Hij.)[[2]] <- new.labels
                   output.matrix.Hij.[row(output.matrix.Hij.) == 0.5 * col(output.matrix.Hij.)] <- format("", width = 7, justify = "right")
                   output.matrix.Hij.[(row(output.matrix.Hij.)) == (0.5 * col(output.matrix.Hij.) + 0.5)] <- format("", width = 7, justify = "right")
                   output.matrix.Hij. <- noquote(output.matrix.Hij.)
                   output.matrix.Hi. <- matrix(NA, J, 2)
                   output.matrix.Hi.[, 1] <- format(formatC(round(Hi., 3), digits = 3, format = "f"), width = 7, justify = "right")
                   output.matrix.Hi.[, 2] <- format(paste("(", formatC(round(se.Hi., 3), digits = 3, format = "f"), ")", sep = ""), width = 7, justify = "right")
                   dimnames(output.matrix.Hi.) <- list(labels, c("Item H", "se"))
                   output.matrix.Hi. <- noquote(output.matrix.Hi.)
                   output.matrix.H. <- matrix(NA, 1, 2)
                   output.matrix.H.[, 1] <- format(formatC(round(H., 3), digits = 3, format = "f"), width = 7, justify = "right")
                   output.matrix.H.[, 2] <- format(paste("(", formatC(round(se.H., 3), digits = 3, format = "f"), ")", sep = ""), width = 7, justify = "right")
                   dimnames(output.matrix.H.) <- list("", c("Scale H", "se"))
                   output.matrix.H. <- noquote(output.matrix.H.)
                   OL[[group.item]][[group]] <- list(Hij = output.matrix.Hij., Hi = output.matrix.Hi., H = output.matrix.H.)
                } else {
                   if (se) OL[[group.item]][[group]] <- list(Hij = Hij., se.Hij = se.Hij., Hi = Hi., se.Hi = se.Hi., H = H., se.H = se.H.)
                   if (!se) OL[[group.item]][[group]] <- list(Hij = Hij., Hi = Hi., H = H.)            
                }
             }# end else (warning)
          }# i-loop through group.names
       }#end if (group.var!=NULL)
    }
    return(OL)    
}#end function coefH
 
"string2integer" <- function(s) as.numeric(unlist(strsplit(s,NULL)))

"all.patterns" <- function(J,m){
  grid <- list()
  j <- 0;
  p <- m^J
  for (j in 1:J){
    grid <- c(grid, j)
    grid[[j]] <- 0:(m-1)
  }
  X <- t(expand.grid(grid))
  dimnames(X) <- NULL
  return(X[J:1,])
}

"weights" <-
# X: Data matrix N x 2 of integer scores [0,1, ..., maxx]
# w: Guttman weights 1 x g^2
# depends on "all.patterns"
function(X, maxx=max.x, minx=0){
 max.x <- max(X)
 g <- maxx + 1
 N <- nrow(X)
 if (ncol(X) != 2){
   warning('X contains more than two columns. Only first two columns will be used')
   X <- X[,1:2]
 }
# Compute order of the ISRFs
 if (maxx == 1) tmp.1 <- matrix(apply(X,2,tabulate, maxx), nrow=1) else tmp.1 <- apply(X,2,tabulate, maxx)
 tmp.2 <- apply(tmp.1,2,function(x) rev(cumsum(rev(x))))+runif(2*maxx,0,1e-3)

 # runif is added to avoid equal ranks
 order.of.ISRFs <- matrix(rank(-tmp.2),1,maxx*2)
# Compute
 Y <- matrix(all.patterns(2,g),nrow=1)
 Z <- matrix(rep(Y, maxx), nrow = maxx, byrow = TRUE)
 Z <- ifelse(Z < row(Z),0,1)
 Z <- matrix(as.vector(Z), ncol = maxx*2, byrow = T)
# COMPUTE WEIGHTS
 Z <- Z[,order(order.of.ISRFs)]
 w <- matrix(apply(Z,1,function(x){sum(x*cumsum(abs(x-1)))}),nrow=1)
 return(w)
}


"phi" <- function(A,f, action){
# Numerical values are translations h(A %*% f) = A %*% f -
  eps = 1E-80;
  switch(action,
    "identity" = A %*% f,
    "exp"      = A %*% exp(f),
    "log"      = A %*% log(abs(f)+eps),
    "sqrt"     = A %*% sqrt(f),
    "xlogx"    = A %*% (-f*log(f+eps)),
    "xbarx"    = A %*% (f*(1-f))  # x(1-x)
  )
}

"dphi" <- function(A,f,df, action){
  eps=1E-80;
  switch(action,
    "identity" = A %*% df,
    "exp"      = A %*% (as.numeric(exp(f)) * df),
    "log"      = A %*% (as.numeric(1/(f+eps)) * df),
    "sqrt"     = A %*% (as.numeric(1/(2*sqrt(f))) * df),
    "xlogx"    = A %*% (as.numeric(-1-log(f+eps)) * df),
    "xbarx"    = A %*% (as.numeric(1-2*f) * df),  # x(1-x)
  )
}

"complete.observed.frequencies" <- function(data,J,m, order.items=FALSE){
  if(order.items) order <- rev(order(apply(data,2,mean))) else order <- 1:J
  data <- as.matrix(data[,order])
  t.R <- cbind(t(all.patterns(J,m)),0)
  p <- m^J
  N <- nrow(data)
  for (i in 1:p){
    size <- abs(data - matrix(1,N,1) %*% t.R[i,1:J]) %*% matrix(1,J,1) == 0
    t.R[i,J+1] <- length(size[size==TRUE])
  }
  return(matrix(t.R[,J+1]))
}
