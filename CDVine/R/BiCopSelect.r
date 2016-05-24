BiCopSelect <- function(u1, u2, familyset = NA, selectioncrit = "AIC", indeptest = FALSE, level = 0.05) {
  if (is.null(u1) == TRUE || is.null(u2) == TRUE) 
    stop("u1 and/or u2 are not set or have length zero.")
  if (length(u1) != length(u2)) 
    stop("Lengths of 'u1' and 'u2' do not match.")
  if (length(u1) < 2) 
    stop("Number of observations has to be at least 2.")
  if (any(u1 > 1) || any(u1 < 0)) 
    stop("Data has be in the interval [0,1].")
  if (any(u2 > 1) || any(u2 < 0)) 
    stop("Data has be in the interval [0,1].")
  if (!is.na(familyset[1])) 
    for (i in 1:length(familyset)) if (!(familyset[i] %in% c(0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 13, 14, 
      16, 17, 18, 19, 20, 23, 24, 26, 27, 28, 29, 30, 33, 34, 36, 37, 38, 39, 40))) 
      stop("Copula family not implemented.")
  if (selectioncrit != "AIC" && selectioncrit != "BIC") 
    stop("Selection criterion not implemented.")
  if (level < 0 & level > 1) 
    stop("Significance level has to be between 0 and 1.")
  
  out <- list()
  
  data1 <- u1
  data2 <- u2
  
  if (!is.na(familyset[1]) & any(familyset == 0)) {
    
    out$p.value.indeptest <- NA
    out$family <- 0
    out$par <- c(0, 0)
    
  } else {
    
    if (!is.na(familyset[1]) && (!any(c(1, 2, 5, 23, 24, 26:30, 33, 34, 36:40) %in% familyset) || !any(c(1:10, 
      13, 14, 16:20) %in% familyset))) 
      stop("'familyset' has to include at least one bivariate copula family for positive and one for negative dependence.")
    
    emp_tau <- fasttau(data1, data2)
    
    if (indeptest == TRUE) {
      out$p.value.indeptest <- BiCopIndTest(data1, data2)$p.value
    } else {
      out$p.value.indeptest <- NA
    }
    
    if (!is.na(out$p.value.indeptest) & out$p.value.indeptest >= level) {
      
      out$family <- 0
      out$par <- c(0, 0)
      
    } else {
      
      start <- list()
      start[[1]] <- sin(pi * emp_tau/2)
      start[[2]] <- c(sin(emp_tau * pi/2), 10)
      start[[3]] <- start[[13]] <- 2 * abs(emp_tau)/(1 - abs(emp_tau))
      start[[4]] <- start[[14]] <- 1/(1 - abs(emp_tau))
      start[[5]] <- Frank.itau.JJ(emp_tau)
      start[[6]] <- start[[16]] <- Joe.itau.JJ(abs(emp_tau))
      start[[7]] <- start[[17]] <- c(0.5, 1.5)
      start[[8]] <- start[[18]] <- c(1.5, 1.5)
      start[[9]] <- start[[19]] <- c(1.5, 0.5)
      start[[10]] <- start[[20]] <- c(1.5, 0.5)
      start[[23]] <- start[[33]] <- -2 * abs(emp_tau)/(1 - abs(emp_tau))
      start[[24]] <- start[[34]] <- -1/(1 - abs(emp_tau))
      start[[26]] <- start[[36]] <- -Joe.itau.JJ(abs(emp_tau))
      start[[27]] <- start[[37]] <- c(-0.5, -1.5)
      start[[28]] <- start[[38]] <- c(-1.5, -1.5)
      start[[29]] <- start[[39]] <- c(-1.5, -0.5)
      start[[30]] <- start[[40]] <- c(-1.5, -0.5)
      
      if (emp_tau < 0) {
        todo <- c(1, 2, 5, 23, 24, 26:30, 33, 34, 36:40)
      } else {
        todo <- c(1:10, 13, 14, 16:20)
      }
      if (!is.na(familyset[1])) 
        todo <- todo[which(todo %in% familyset)]
      
      optiout <- list()
      
      if (any(todo == 2)) {
        optiout[[2]] <- suppressWarnings(BiCopEst(data1, data2, family = 2, max.df = 30))
        optiout[[2]]$par <- c(optiout[[2]]$par, optiout[[2]]$par2)
        if (optiout[[2]]$par[2] >= 30 - 0.01) {
          todo[todo == 2] <- 1
          todo <- unique(todo)
          optiout[[2]] <- list()
        }
      }
      
      if (any(todo == 7)) {
        optiout[[7]] <- MLE_intern(cbind(data1, data2), start[[7]], 7)
        if (optiout[[7]]$par[1] <= 0.1 | optiout[[7]]$par[2] <= 1.1) {
          if (optiout[[7]]$par[1] <= 0.1) {
          todo[todo == 7] <- 4
          todo <- unique(todo)
          } else if (optiout[[7]]$par[2] <= 1.1) {
          todo[todo == 7] <- 3
          todo <- unique(todo)
          }
          optiout[[7]] <- list()
        }
      }
      
      if (any(todo == 8)) {
        optiout[[8]] <- MLE_intern(cbind(data1, data2), start[[8]], 8)
        if (optiout[[8]]$par[1] <= 1.1 | optiout[[8]]$par[2] <= 1.1) {
          if (optiout[[8]]$par[1] <= 1.1) {
          todo[todo == 8] <- 4
          todo <- unique(todo)
          } else if (optiout[[8]]$par[2] <= 1.1) {
          todo[todo == 8] <- 6
          todo <- unique(todo)
          }
          optiout[[8]] <- list()
        }
      }
      
      if (any(todo == 9)) {
        optiout[[9]] <- MLE_intern(cbind(data1, data2), start[[9]], 9)
        if (optiout[[9]]$par[1] <= 1.1 | optiout[[9]]$par[2] <= 0.1) {
          if (optiout[[9]]$par[1] <= 1.1) {
          todo[todo == 9] <- 3
          todo <- unique(todo)
          } else if (optiout[[9]]$par[2] <= 0.1) {
          todo[todo == 9] <- 6
          todo <- unique(todo)
          }
          optiout[[9]] <- list()
        }
      }
      
      if (any(todo == 10)) {
        optiout[[10]] <- MLE_intern(cbind(data1, data2), start[[10]], 10)
        if (optiout[[10]]$par[2] >= 0.99) {
          todo[todo == 10] <- 6
          todo <- unique(todo)
          optiout[[10]] <- list()
        }
      }
      
      if (any(todo == 17)) {
        optiout[[17]] <- MLE_intern(cbind(data1, data2), start[[17]], 17)
        if (optiout[[17]]$par[1] <= 0.1 | optiout[[17]]$par[2] <= 1.1) {
          if (optiout[[17]]$par[1] <= 0.1) {
          todo[todo == 17] <- 14
          todo <- unique(todo)
          } else if (optiout[[17]]$par[2] <= 1.1) {
          todo[todo == 17] <- 13
          todo <- unique(todo)
          }
          optiout[[17]] <- list()
        }
      }
      
      if (any(todo == 18)) {
        optiout[[18]] <- MLE_intern(cbind(data1, data2), start[[18]], 18)
        if (optiout[[18]]$par[1] <= 1.1 | optiout[[18]]$par[2] <= 1.1) {
          if (optiout[[18]]$par[1] <= 1.1) {
          todo[todo == 18] <- 14
          todo <- unique(todo)
          } else if (optiout[[18]]$par[2] <= 1.1) {
          todo[todo == 18] <- 16
          todo <- unique(todo)
          }
          optiout[[18]] <- list()
        }
      }
      
      if (any(todo == 19)) {
        optiout[[19]] <- MLE_intern(cbind(data1, data2), start[[19]], 19)
        if (optiout[[19]]$par[1] <= 1.1 | optiout[[19]]$par[2] <= 0.1) {
          if (optiout[[19]]$par[1] <= 1.1) {
          todo[todo == 19] <- 13
          todo <- unique(todo)
          } else if (optiout[[19]]$par[2] <= 0.1) {
          todo[todo == 19] <- 16
          todo <- unique(todo)
          }
          optiout[[19]] <- list()
        }
      }
      
      if (any(todo == 20)) {
        optiout[[20]] <- MLE_intern(cbind(data1, data2), start[[20]], 20)
        if (optiout[[20]]$par[2] >= 0.99) {
          todo[todo == 20] <- 16
          todo <- unique(todo)
          optiout[[20]] <- list()
        }
      }
      
      if (any(todo == 27)) {
        optiout[[27]] <- MLE_intern(cbind(data1, data2), start[[27]], 27)
        if (optiout[[27]]$par[1] >= -0.1 | optiout[[27]]$par[2] >= -1.1) {
          if (optiout[[27]]$par[1] >= -0.1) {
          todo[todo == 27] <- 24
          todo <- unique(todo)
          } else if (optiout[[27]]$par[2] >= -1.1) {
          todo[todo == 27] <- 23
          todo <- unique(todo)
          }
          optiout[[27]] <- list()
        }
      }
      
      if (any(todo == 28)) {
        optiout[[28]] <- MLE_intern(cbind(data1, data2), start[[28]], 28)
        if (optiout[[28]]$par[1] >= -1.1 | optiout[[28]]$par[2] >= -1.1) {
          if (optiout[[28]]$par[1] >= -1.1) {
          todo[todo == 28] <- 24
          todo <- unique(todo)
          } else if (optiout[[28]]$par[2] >= -1.1) {
          todo[todo == 28] <- 26
          todo <- unique(todo)
          }
          optiout[[28]] <- list()
        }
      }
      
      if (any(todo == 29)) {
        optiout[[29]] <- MLE_intern(cbind(data1, data2), start[[29]], 29)
        if (optiout[[29]]$par[1] >= -1.1 | optiout[[29]]$par[2] >= -0.1) {
          if (optiout[[29]]$par[1] >= -1.1) {
          todo[todo == 29] <- 23
          todo <- unique(todo)
          } else if (optiout[[29]]$par[2] >= -0.1) {
          todo[todo == 29] <- 26
          todo <- unique(todo)
          }
          optiout[[29]] <- list()
        }
      }
      
      if (any(todo == 30)) {
        optiout[[30]] <- MLE_intern(cbind(data1, data2), start[[30]], 30)
        if (optiout[[30]]$par[2] <= -0.99) {
          todo[todo == 30] <- 26
          todo <- unique(todo)
          optiout[[30]] <- list()
        }
      }
      
      if (any(todo == 37)) {
        optiout[[37]] <- MLE_intern(cbind(data1, data2), start[[37]], 37)
        if (optiout[[37]]$par[1] >= -0.1 | optiout[[37]]$par[2] >= -1.1) {
          if (optiout[[37]]$par[1] >= -0.1) {
          todo[todo == 37] <- 34
          todo <- unique(todo)
          } else if (optiout[[37]]$par[2] >= -1.1) {
          todo[todo == 37] <- 33
          todo <- unique(todo)
          }
          optiout[[37]] <- list()
        }
      }
      
      if (any(todo == 38)) {
        optiout[[38]] <- MLE_intern(cbind(data1, data2), start[[38]], 38)
        if (optiout[[38]]$par[1] >= -1.1 | optiout[[38]]$par[2] >= -1.1) {
          if (optiout[[38]]$par[1] >= -1.1) {
          todo[todo == 38] <- 34
          todo <- unique(todo)
          } else if (optiout[[38]]$par[2] >= -1.1) {
          todo[todo == 38] <- 36
          todo <- unique(todo)
          }
          optiout[[38]] <- list()
        }
      }
      
      if (any(todo == 39)) {
        optiout[[39]] <- MLE_intern(cbind(data1, data2), start[[39]], 39)
        if (optiout[[39]]$par[1] >= -1.1 | optiout[[39]]$par[2] >= -0.1) {
          if (optiout[[39]]$par[1] >= -1.1) {
          todo[todo == 39] <- 33
          todo <- unique(todo)
          } else if (optiout[[39]]$par[2] >= -0.1) {
          todo[todo == 39] <- 36
          todo <- unique(todo)
          }
          optiout[[39]] <- list()
        }
      }
      
      if (any(todo == 40)) {
        optiout[[40]] <- MLE_intern(cbind(data1, data2), start[[40]], 40)
        if (optiout[[40]]$par[2] <= -0.99) {
          todo[todo == 40] <- 36
          todo <- unique(todo)
          optiout[[40]] <- list()
        }
      }
      
      for (i in todo[!(todo %in% c(2, 7:10, 17:20, 27:30, 37:40))]) {
        
        optiout[[i]] <- MLE_intern(cbind(data1, data2), start[[i]], i)
      }
      
      if (selectioncrit == "AIC") {
        
        AICs <- rep(Inf, 40)
        
        for (i in todo) {
          if (i %in% c(2, 7:10, 17:20, 27:30, 37:40)) {
          ll <- sum(log(BiCopPDF(data1, data2, i, optiout[[i]]$par, optiout[[i]]$par[2])))
          AICs[i] <- -2 * ll + 4
          } else {
          ll <- sum(log(BiCopPDF(data1, data2, i, optiout[[i]]$par)))
          AICs[i] <- -2 * ll + 2
          }
        }
        
        out$family <- todo[which.min(AICs[todo])]
        
      } else if (selectioncrit == "BIC") {
        
        BICs <- rep(Inf, 40)
        
        for (i in todo) {
          if (i %in% c(2, 7:10, 17:20, 27:30, 37:40)) {
          ll <- sum(log(BiCopPDF(data1, data2, i, optiout[[i]]$par, optiout[[i]]$par[2])))
          BICs[i] <- -2 * ll + 2 * log(length(data1))
          } else {
          ll <- sum(log(BiCopPDF(data1, data2, i, optiout[[i]]$par)))
          BICs[i] <- -2 * ll + log(length(data1))
          }
        }
        
        out$family <- todo[which.min(BICs[todo])]
        
      }
      
      out$par <- optiout[[out$family]]$par
      if (!(out$family %in% c(2, 7:10, 17:20, 27:30, 37:40))) 
        out$par[2] <- 0
      
      # out$emp_tau = emp_tau
      
      # out$CondOn.1 =
      # .C('Hfunc2',as.integer(out$family),as.integer(length(data1)),as.double(data1),as.double(data2),as.double(out$par[1]),as.double(out$par[2]),as.double(rep(0,length(data1))),PACKAGE='VineCopula')[[7]]
      # out$CondOn.2 =
      # .C('Hfunc1',as.integer(out$family),as.integer(length(data1)),as.double(data2),as.double(data1),as.double(out$par[1]),as.double(out$par[2]),as.double(rep(0,length(data1))),PACKAGE='VineCopula')[[7]]
      
    }
  }
  
  out$par2 <- out$par[2]
  out$par <- out$par[1]
  
  return(out)
  
} 
