GP.csh <-
function (Z_mat, fixed_effects, control) 
  {
      gr.eta <- function(eta, X, Y, Z, Sig.mat, ybetas, sigmas, 
          G.inv, nyear, n_eta) {
          gr.y <- colSums(Z * (as.vector((Y - X %*% ybetas - Z %*% 
              eta)/as.vector(Sig.mat %*% sigmas^2))))
          gr.p.eta <- -G.inv %*% eta
          -as.vector(gr.y + gr.p.eta)
      }
      H.eta <- function(sigmas, cross_Z_j, Sig.mat, G.inv, nyear, 
          n_eta) {
          h.eta <- G.inv
          h.y <- Matrix(0, n_eta, n_eta)
          for (j in 1:nyear) {
              h.y <- h.y + (1/sigmas[j]^2) * cross_Z_j[[j]]
          }
          rm(j)
          forceSymmetric(h.eta + h.y)
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
      reduce.G <- function(G, nstudent, nyear, nteacher, Kg) {
          if (!is.null(dim(G)[2])) {
              temp_mat <- G
              gam_stu <- c(temp_mat[1, 1])
              resA <- gam_stu
              index1 <- nstudent
              for (j in 1:nyear) {
                  temp_mat_j <- temp_mat[(index1 + 1):(index1 + 
                    Kg[j]), (index1 + 1):(index1 + Kg[j])]
                  resA <- c(resA, ltriangle(as.matrix(temp_mat_j)))
                  index1 <- index1 + nteacher[j] * Kg[j]
              }
              resA
          }
          else {
              resB <- as.vector(G[1]) * .symDiagonal(nstudent)
              index <- c(2)
              for (j in 1:nyear) {
                  ne <- (Kg[j] * (Kg[j] + 1))/2
                  resB <- suppressMessages(bdiag(resB, suppressMessages(kronecker(suppressMessages(.symDiagonal(nteacher[j])), 
                    ltriangle(G[index:(index + ne - 1)])))))
                  index <- index + ne
              }
              rm(j)
              resB
          }
      }
      update.eta <- function(X, Y, Z, cross_Z_j, Sig.mat, 
          ybetas, sigmas, G, nyear, n_eta, cons.logLik) {
          G.chol <- chol(G)
          G.inv <- chol2inv(G.chol)
          H <- H.eta(sigmas = sigmas, cross_Z_j = cross_Z_j, Sig.mat = Sig.mat, 
              G.inv = G.inv, nyear = nyear, n_eta = n_eta)
          chol.H <- chol(H)
          var.eta <- as.matrix(solve(H))
          rm(H)
          eta<-var.eta%*%colSums(Z *(as.vector((Y - X %*% ybetas)/as.vector(Sig.mat %*% sigmas^2))))
          log.p.eta <- -(length(eta)/2) * log(2 * pi) - sum(log(diag(G.chol))) - 
              0.5 * crossprod(eta, G.inv) %*% eta
          log.p.y <- sum(dnorm(Y, as.vector(X %*% ybetas + Z %*% 
              eta), as.vector(Sig.mat %*% sigmas), log = TRUE))
          res <- var.eta
          attr(res, "likelihood") <- as.vector(cons.logLik + log.p.eta + 
              log.p.y - 0.5 * (2 * sum(log(diag(chol.H)))))
          attr(res, "eta") <- eta
          res
      }
      Score <- function(thetas, eta, ybetas, X, Y, Z, cross_Z_j, 
          X_j, cross_X_j, Y_j, Z_j, Sig.mat, year.count, n_ybeta, 
          nyear, n_eta, nstudent, nteacher, Kg, cons.logLik, con) {
          sigmas_sq <- pmax(thetas[(1):(nyear)], 1e-08)
          G <- thetas[seq(nyear + 1, length(thetas))]
          G <- reduce.G(G = G, nstudent = nstudent, nyear = nyear, 
              nteacher = nteacher, Kg = Kg)
          new.eta <- update.eta(X = X, Y = Y, Z = Z, 
              cross_Z_j = cross_Z_j, Sig.mat = Sig.mat, ybetas = ybetas, 
              sigmas = sqrt(sigmas_sq), G = G, nyear = nyear, n_eta = n_eta, 
              cons.logLik = cons.logLik)
          eta <- attr(new.eta, "eta")
          eta.hat <- eta
          var.eta.hat <- new.eta
          rm(new.eta)
          score.sigmas_sq <- as.vector(-1/(2 * sigmas_sq) * year.count + 
              1/(2 * sigmas_sq^2) * ((t(Sig.mat) %*% (as.vector(Y - 
                  X %*% ybetas) * as.vector(Y - X %*% ybetas - 
                  2 * Z %*% eta.hat))) + sapply(cross_Z_j, function(x) sum(diag(x %*% 
                  var.eta.hat))) + sapply(cross_Z_j, function(x) as.numeric(crossprod(eta.hat, 
                  x) %*% eta.hat))))
          temp_mat <- G
          gam_stu <- temp_mat[1, 1]
          sv_gam_stu <- 1/gam_stu
          gam_t <- list()
          sv_gam_t <- list()
          index1 <- nstudent
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
          gam_stu_sc <- sum(diag(score_mat)[1:nstudent])
          score.eta.stu <- drop(-0.5 * (nstudent * sv_gam_stu - 
              sv_gam_stu %*% gam_stu_sc %*% sv_gam_stu))
          score.G <- score.eta.stu * diag(nstudent)
          gam_t_sc <- list()
          index1 <- nstudent
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
          rm(j, k)
          -c(score.sigmas_sq, reduce.G(G = score.G, nstudent = nstudent, 
              nyear = nyear, nteacher = nteacher, Kg = Kg))
      }
      update.ybeta <- function(X_j, cross_X_j, Y_j, Z_j, sigmas, 
          eta.hat, n_ybeta, nyear) {
          A.ybeta <- matrix(0, n_ybeta, n_ybeta)
          for (g in 1:nyear) {
              A.ybeta <- A.ybeta + 1/(sigmas[g]^2) * cross_X_j[[g]]
          }
          rm(g)
          B.ybeta <- numeric(n_ybeta)
          for (g in 1:nyear) {
              B.ybeta = B.ybeta + 1/(sigmas[g]^2) * (crossprod(X_j[[g]], 
                  (Y_j[[g]] - Z_j[[g]] %*% eta.hat)))
          }
          rm(g)
          as.vector(solve(A.ybeta, B.ybeta))
      }
      Z_mat$year <- as.factor(Z_mat$year)
      nyear <- nlevels(Z_mat$year)
      Z_mat$teacher <- as.character(Z_mat$teacher)
      Z_mat.full <- Z_mat
      Z_mat <- Z_mat[!is.na(Z_mat$y), ]
      Z_mat.full <- Z_mat.full[which((Z_mat.full$student %in% Z_mat$student)), 
          ]
      Ny <- length(Z_mat$y)
      nstudent <- length(unique(Z_mat$student))
      year.count <- numeric(nyear)
      for (j in 1:nyear) {
          year.count[j] <- length(Z_mat[Z_mat$year == j, ]$y)
      }
      rm(j)
      RE_s_start_pos <- 1
      Z_mat <- Z_mat[order(Z_mat$student), ]
      s_effects <- paste("ygam.", unique(Z_mat$student), sep = "")
      RE_mat <- sparse.model.matrix(~as.factor(student) + 0, Z_mat)
      colnames(RE_mat) <- s_effects
      Kg <- nyear - 1:nyear + 1
      RE_mat <- RE_mat[order(Z_mat$year, Z_mat$teacher), ]
      Z_mat <- Z_mat[order(Z_mat$year, Z_mat$teacher), ]
      na_list <- grep("^NA", Z_mat$teacher)
      if (length(na_list) > 0) {
          teachyearcomb <- unique(cBind(Z_mat[-na_list, ]$year, 
              Z_mat[-na_list, ]$teacher))
      }
      else {
          teachyearcomb <- unique(cBind(Z_mat$year, Z_mat$teacher))
      }
      nteach_effects <- sum(nyear - as.numeric(teachyearcomb[, 
          1]) + 1)
      teacheffZ_mat <- Matrix(0, nrow = nrow(Z_mat), ncol = nteach_effects)
      t_effects <- rep(NA, nteach_effects)
      indx <- 1
      eblup.tracker <- matrix(0, 0, 3)
      for (k in 1:nrow(teachyearcomb)) {
          student_subset <- Z_mat.full$student[Z_mat.full$year == 
              teachyearcomb[k, 1] & Z_mat.full$teacher == teachyearcomb[k, 
              2]]
          for (yr in teachyearcomb[k, 1]:nyear) {
              if (sum(is.element(Z_mat$student, student_subset) & 
                  Z_mat$year == yr) != 0) {
                  teacheffZ_mat[is.element(Z_mat$student, student_subset) & 
                    Z_mat$year == yr & !is.na(Z_mat$y), indx] <- 1
              }
              t_effects[indx] <- paste(teachyearcomb[k, 1], "_", 
                  teachyearcomb[k, 2], "_", yr, sep = "")
              indx <- indx + 1
              eblup.tracker <- rBind(eblup.tracker, c(teachyearcomb[k, 
                  ], yr))
          }
      }
      eblup.tracker <- as.data.frame(eblup.tracker)
      colnames(eblup.tracker) <- c("year", "teacher", "effect_year")
      rm(k, yr, indx)
      nteacher <- as.vector(tapply(teachyearcomb[, 2], teachyearcomb[, 
          1], length))
      colnames(teacheffZ_mat) <- t_effects
      RE_mat <- cBind(RE_mat, teacheffZ_mat)
      Sig.mat <- sparse.model.matrix(~as.factor(Z_mat$year) + 0)
      X_mat <- sparse.model.matrix(fixed_effects, Z_mat, drop.unused.levels = TRUE)
      X_mat <- X_mat[, !(colSums(abs(X_mat)) == 0)]
      if (rankMatrix(X_mat,method = 'qrLINPACK')[1] != dim(X_mat)[2]) {
          stop("WARNING: Fixed-effects design matrix not full-rank")

      }
      n_eta <- nstudent + nteach_effects
      n_ybeta <- dim(X_mat)[2]
      Z <- Matrix(RE_mat)
      cross_Z <- crossprod(Z)
      Y <- as.vector(Z_mat$y)
      X <- Matrix(X_mat)
      cross_Z_j <- list()
      X_j <- list(NULL)
      cross_X_j <- list(NULL)
      Y_j <- list(NULL)
      Z_j <- list(NULL)
      for (j in 1:nyear) {
          cross_Z_j[[j]] <- crossprod(Matrix(RE_mat[Z_mat$year == 
              j, ]))
          X_j[[j]] <- X_mat[Z_mat$year == j, ]
          Y_j[[j]] <- as.vector(Z_mat[Z_mat$year == j, ]$y)
          Z_j[[j]] <- RE_mat[Z_mat$year == j, ]
          cross_X_j[[j]] <- crossprod(X_j[[j]])
      }
      rm(j)
      eta.hat <- numeric(n_eta)
      var.eta.hat <- Matrix(0, n_eta, n_eta)
      sigmas <- c(rep(1, nyear))
      ybetas <- update.ybeta(X_j = X_j, cross_X_j = cross_X_j, 
          Y_j = Y_j, Z_j = Z_j, sigmas = sigmas, eta.hat = eta.hat, 
          n_ybeta = n_ybeta, nyear = nyear)
      names(ybetas) <- colnames(X_mat)
      G <- 100 * .symDiagonal(nstudent + nteach_effects)
      cons.logLik <- 0.5 * n_eta * log(2 * pi)
      iter <- control$max.iter.EM
      Y.mat <- Matrix(0, iter, n_ybeta)
      G.mat <- Matrix(0, iter, length(reduce.G(G = G, nstudent = nstudent, 
          nyear = nyear, nteacher = nteacher, Kg = Kg)))
      sig.mat <- Matrix(0, iter, nyear)
      lgLik <- numeric(iter)
      conv <- FALSE
      if (control$verbose) cat("Beginning EM algorithm\n")
      flush.console()
      for (it in 1:iter) {
          ptm <- proc.time()[3]
          new.eta <- update.eta(X = X, Y = Y, Z = Z, 
              cross_Z_j = cross_Z_j, Sig.mat = Sig.mat, ybetas = ybetas, 
              sigmas = sigmas, G = G, nyear = nyear, n_eta = n_eta, 
              cons.logLik = cons.logLik)
          Y.mat[it, ] <- c(ybetas)
          sig.mat[it, ] <- c(sigmas)
          G.mat[it, ] <- reduce.G(G = G, nstudent = nstudent, nyear = nyear, 
              nteacher = nteacher, Kg = Kg)
          eta.hat <- attr(new.eta, "eta")
          var.eta.hat <- new.eta
          lgLik[it] <- attr(new.eta, "likelihood")
          rm(new.eta)
          thets1 <- c(Y.mat[it - 1, ], sig.mat[it - 1, ], G.mat[it - 
              1, ])
          thets2 <- c(Y.mat[it, ], sig.mat[it, ], G.mat[it, ])
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
                    cat("yearly error variances:", round(sigmas^2, 
                      4), "\n")
                    cat("student variance\n")
                    print(round(as.matrix(gam_stu), 4))
                    cat("\n")
                    print(try(round(cov2cor(as.matrix(gam_stu)), 
                      4), silent = TRUE))
                    for (j in 1:nyear) {
                      cat("\ngamma_teach_year", j, "\n")
                      print(round(as.matrix(gam_t[[j]]), 4))
                      cat("\n")
                      print(try(round(cov2cor(as.matrix(gam_t[[j]])), 
                        4), silent = TRUE))
                      flush.console()
                    }
                    rm(j)
                  }
                  break
              }
                            if (check.lik) {
                  conv <- TRUE
                  if (it==iter) {
                    cat("\n\n Maximum number of iterations reached.\n")
                    cat("\n\niter:", it, "\n")
                    cat("log-likelihood:", sprintf("%.7f", lgLik[it]), 
                      "\n")
                    cat("change in loglik:", sprintf("%.7f", lgLik[it] - 
                      lgLik[it - 1]), "\n")
                    cat("fixed effects:", round(ybetas, 4), "\n")
                    cat("yearly error variances:", round(sigmas^2, 
                      4), "\n")
                    cat("student variance\n")
                    print(round(as.matrix(gam_stu), 4))
                    cat("\n")
                    print(try(round(cov2cor(as.matrix(gam_stu)), 
                      4), silent = TRUE))
                    for (j in 1:nyear) {
                      cat("\ngamma_teach_year", j, "\n")
                      print(round(as.matrix(gam_t[[j]]), 4))
                      cat("\n")
                      print(try(round(cov2cor(as.matrix(gam_t[[j]])), 
                        4), silent = TRUE))
                      flush.console()
                    }
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
              cat("yearly error variances:", round(sigmas^2, 4), 
                  "\n")
              cat("student variance\n")
              print(round(as.matrix(gam_stu), 4))
              cat("\n")
              print(try(round(cov2cor(as.matrix(gam_stu)), 4), 
                  silent = TRUE))
              for (j in 1:nyear) {
                  cat("\ngamma_teach_year", j, "\n")
                  print(round(as.matrix(gam_t[[j]]), 4))
                  cat("\n")
                  print(try(round(cov2cor(as.matrix(gam_t[[j]])), 
                    4), silent = TRUE))
                  flush.console()
              }
              rm(j)
          }
          sigmas2n <- as.vector(t(Sig.mat) %*% (as.vector(1/(Sig.mat %*% 
              year.count)) * ((Y - X %*% ybetas) * (Y - X %*% ybetas - 
              2 * Z %*% eta.hat))) + 1/(year.count) * sapply(cross_Z_j, 
              function(x) as.numeric(sum(diag(x %*% var.eta.hat)))) + 
              1/(year.count) * sapply(cross_Z_j, function(x) as.numeric(crossprod(eta.hat, 
                  x %*% eta.hat))))
          sigmasn <- sqrt(sigmas2n)
          eblup <- cBind(eta.hat, sqrt(diag(var.eta.hat)))
          temp_mat <- var.eta.hat
          rm(var.eta.hat)
          temp_mat <- temp_mat + tcrossprod(eta.hat, eta.hat)
          gam_stu <- (1/nstudent) * sum(diag(temp_mat)[1:nstudent])
          gam_t <- list()
          index1 <- nstudent
          for (j in 1:nyear) {
              gam_t[[j]] <- Matrix(0, Kg[j], Kg[j])
              temp_mat_j <- temp_mat[(index1 + 1):(index1 + nteacher[j] * 
                  Kg[j]), (index1 + 1):(index1 + nteacher[j] * 
                  Kg[j])]
              index2 <- c(1)
              for (k in 1:nteacher[j]) {
                  gam_t[[j]] <- gam_t[[j]] + temp_mat_j[(index2):(index2 + 
                    Kg[j] - 1), (index2):(index2 + Kg[j] - 1)]
                  index2 <- index2 + Kg[j]
              }
              index1 <- index1 + nteacher[j] * Kg[j]
              gam_t[[j]] <- suppressMessages(as(symmpart(suppressMessages(gam_t[[j]]/nteacher[j])), "sparseMatrix"))
          }
          rm(j, k, index2, index1)
          rm(temp_mat)
          ybetasn <- numeric(n_ybeta)
          Gn <- gam_stu * Diagonal(nstudent)
          for (j in 1:nyear) {
              Gn <- suppressMessages(bdiag(Gn, suppressMessages(kronecker(suppressMessages(.symDiagonal(nteacher[j])), 
                  gam_t[[j]]))))
          }
          rm(j)
          ybetasn <- update.ybeta(X_j = X_j, cross_X_j = cross_X_j, 
              Y_j = Y_j, Z_j = Z_j, sigmas = sigmas, eta.hat = eta.hat, 
              n_ybeta = n_ybeta, nyear = nyear)
          ybetas <- ybetasn
          sigmas <- sigmasn
          G <- Gn
         if (control$verbose) cat("Iteration Time: ", proc.time()[3] - ptm, "\n")
          flush.console()
      }
      names(ybetas) <- colnames(X_mat)
      thetas <- c(sigmas^2, reduce.G(G = G, nstudent = nstudent, 
          nyear = nyear, nteacher = nteacher, Kg = Kg))
      lgLik.hist <- lgLik
      lgLik <- lgLik[it]
      if (!control$hessian) Hessian<-NA
      if (control$hessian) {
          if (control$verbose) cat("Calculating Hessian of the variance components...")
          flush.console()
          if (control$hes.method == "richardson") {
              Hessian <- jacobian(Score, thetas, eta = eta.hat, 
                  ybetas = ybetas, X = X, Y = Y, Z = Z, cross_Z_j = cross_Z_j, 
                  X_j = X_j, cross_X_j = cross_X_j, Y_j = Y_j, 
                  Z_j = Z_j, Sig.mat = Sig.mat, year.count = year.count, 
                  n_ybeta = n_ybeta, nyear = nyear, n_eta = n_eta, 
                  nstudent = nstudent, nteacher = nteacher, Kg = Kg, 
                  cons.logLik = cons.logLik, con = control)
          }
          else {
              Hessian <- jacobian(Score, thetas, method = "simple", 
                  eta = eta.hat, ybetas = ybetas, X = X, Y = Y, 
                  Z = Z, cross_Z_j = cross_Z_j, X_j = X_j, cross_X_j = cross_X_j, 
                  Y_j = Y_j, Z_j = Z_j, Sig.mat = Sig.mat, year.count = year.count, 
                  n_ybeta = n_ybeta, nyear = nyear, n_eta = n_eta, 
                  nstudent = nstudent, nteacher = nteacher, Kg = Kg, 
                  cons.logLik = cons.logLik, con = control)
          }
          std_errors <- try(c(sqrt(diag(solve(Hessian)))), silent = TRUE)
          Hessian <- round(Hessian, 5)
          hes.warn <- FALSE
          if (any(eigen(Hessian)$values <= 0)) {
             if (control$verbose) cat("Warning: Hessian not PD", "\n")
              flush.console()
              std_errors <- c(rep(NA, length(thetas)))
              hes.warn <- TRUE
          }
      }
      R_inv <- Diagonal(x = as.vector(Sig.mat %*% (1/(sigmas^2))))
      c.temp <- as(crossprod(X, R_inv) %*% Z, "TsparseMatrix")
      c.1 <- rBind(as(crossprod(X, R_inv) %*% X, "TsparseMatrix"), 
          t(c.temp))
      G.inv <- chol2inv(chol(G))
      c.2 <- rBind(c.temp, H.eta(sigmas, cross_Z_j, Sig.mat, G.inv, 
          nyear, n_eta))
      C_inv <- cBind(c.1, c.2)
      C <- solve(C_inv)
      eblup_stderror <- sqrt(diag(C)[-c(1:n_ybeta)])
      ybetas_stderror <- sqrt(diag(C)[1:n_ybeta])
      eblup <- cBind(eta.hat, eblup_stderror)
      eblup <- round(eblup, 4)
      eblup <- eblup[-(1:nstudent), ]
      eblup <- as.data.frame(eblup)
      eblup.tracker <- as.data.frame(eblup.tracker)
      eblup <- as.data.frame(cbind(eblup.tracker, eblup))
      colnames(eblup) <- c("teacher_year", "teacher", "effect_year", 
          "EBLUP", "std_error")
      eblup$teacher <- as.character(eblup$teacher)
      t_lab <- as.vector(NULL)
      for (j in 1:nyear) {
          ne <- (Kg[j] * (Kg[j] + 1))/2
          y <- c(NULL)
          x <- c(NULL)
          for (k in 1:Kg[j]) {
              x <- c(x, k:Kg[j])
              y <- c(y, rep(k, (Kg[j] - k + 1)))
          }
          t_lab <- c(t_lab, paste("teacher effect from year", rep(j, 
              ne), ":[", x, ",", y, "]", sep = ""))
      }
      rm(j)
      effect_la <- c(names(ybetas), paste("sigma^2_", control$key[1:nyear,1], 
          sep = ""), "student", t_lab)
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
      if (control$verbose) cat("done.\n")
      mresid <- as.numeric(Y - X %*% ybetas)
      cresid <- as.numeric(mresid - Z %*% eta.hat)
      yhat <- as.numeric(X %*% ybetas + Z %*% eta.hat)
      yhat.m <- as.numeric(X %*% ybetas)
      for (i in 1:nyear) {
          gam_t[[i]] <- round(as.matrix(gam_t[[i]]), 4)
          colnames(gam_t[[i]]) <- paste("year", control$key[i:nyear,1], sep = "")
          rownames(gam_t[[i]]) <- colnames(gam_t[[i]])
      }
      Hessian <- round(Hessian, 5)
      stu.cov <- round(as.matrix(gam_stu), 4)
      rchol <- try(chol(Diagonal(x = as.vector(Sig.mat %*% (1/sigmas^2)))))
      yhat.s <- try(as.vector(rchol %*% (yhat)))
      sresid <- try(as.vector(rchol %*% Y - yhat.s))
    list(loglik = lgLik, teach.effects = eblup, parameters = parameters, 
        Hessian = Hessian, R_i = NA, teach.cov = gam_t, mresid = mresid, 
        cresid = cresid, y = Y, yhat = yhat, stu.cov = gam_stu, 
        num.obs = Ny, num.student = nstudent, num.year = nyear, 
        num.teach = nteacher, yhat.m = yhat.m, sresid = sresid, 
        yhat.s = yhat.s,iter=it,persistence=control$persistence)
}
