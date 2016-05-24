.Rn = function(x, copula, dims) {
  switch(class(copula),
         claytonCopula = {
           Rn.ac.c = -sum(.V.clay.12(copula@parameters, x[,1], x[,2])^2) / sum(.S.clay.12(copula@parameters, x[,1], x[,2])) - 1
           Rn.ac.c
         },
         gumbelCopula = {
           Rn.ac.g = -sum(.V.gumb.12(copula@parameters, x[,1], x[,2])^2) / sum(.S.gumb.12(copula@parameters, x[,1], x[,2])) - 1
           Rn.ac.g
         },
         frankCopula = {
           Rn.ac.f = -sum(.V.fran.12(copula@parameters, x[,1], x[,2])^2) / sum(.S.fran.12(copula@parameters, x[,1], x[,2])) - 1
           Rn.ac.f
         },
         tCopula = {
           psn.sample = qt(x, df = copula@parameters[2])
           Rn.t = -sum(apply(psn.sample, 1, FUN = .V.t, rho = copula@parameters[1], nu = copula@parameters[2])^2)/sum(apply(psn.sample, 1, FUN = .S.t, rho = copula@parameters[1], nu = copula@parameters[2])) - 1
           Rn.t
         },
         normalCopula = {
           sig = matrix(c(1, copula@parameters, copula@parameters, 1), ncol = 2, byrow = T)
           sig.inv = rbind(c(sig[1,1], -sig[1,2]), c(-sig[2,1], sig[2,2]))/det(sig)
           psn.sample = qnorm(x)
           Rn.g = -sum(apply(psn.sample, 1, FUN = .V.ga, sig.inv = sig.inv, dims = dims)[2,]^2)/sum(apply(psn.sample, 1, FUN = .S.ga.12, sig =sig)) - 1
           Rn.g
         })
  
}

.Tn = function(x, copula, B, m, dims, param.est) {
  switch(class(copula),
         claytonCopula = {
           l.ac.c = log(.clay.12.density(copula@parameters, x))
           l.ac.c.b = 0
           if (param.est == T){
             for(b in 1:B){
               l.ac.c.b = c(l.ac.c.b, .clay.12.density(fitCopula(claytonCopula(dim = dims), data = x[-(((b-1)*m+1):(b*m)),], method = "mpl")@estimate, x[((b-1)*m+1):(b*m),]))
             }
           } else {
             for(b in 1:B){
               l.ac.c.b = c(l.ac.c.b, .clay.12.density(copula@parameters, x[((b-1)*m+1):(b*m),]))
             }
           }
           l.ac.c.b = l.ac.c.b[-1]
           stat = sum(l.ac.c - log(l.ac.c.b)) - 1
           stat
         },
         gumbelCopula = {
           l.ac.g = log(.gumb.12.density(copula@parameters, x))
           l.ac.g.b = 0
           if (param.est == T){
             for(b in 1:B){
               l.ac.g.b = c(l.ac.g.b, .gumb.12.density(fitCopula(gumbelCopula(dim = dims), data = x[-(((b-1)*m+1):(b*m)),], method = "mpl")@estimate, x[((b-1)*m+1):(b*m),]))
             }
           } else {
             for(b in 1:B){
               l.ac.g.b = c(l.ac.g.b, .gumb.12.density(copula@parameters, x[((b-1)*m+1):(b*m),]))
             }
           }
           
           l.ac.g.b = l.ac.g.b[-1]
           stat = sum(l.ac.g - log(l.ac.g.b)) - 1
           stat
         },
         frankCopula = {
           l.ac.f = log(.fran.12.density(copula@parameters, x))
           l.ac.f.b = 0
           if (param.est == T){
             for(b in 1:B){
               l.ac.f.b = c(l.ac.f.b, .fran.12.density(fitCopula(frankCopula(dim = dims), data = x[-(((b-1)*m+1):(b*m)),], method = "mpl")@estimate, x[((b-1)*m+1):(b*m),]))
             }
           } else {
             for(b in 1:B){
               l.ac.f.b = c(l.ac.f.b, .fran.12.density(copula@parameters, x[((b-1)*m+1):(b*m),]))
             }
           }
           l.ac.f.b = l.ac.f.b[-1]
           stat = sum(l.ac.f - log(l.ac.f.b)) - 1
           stat
         },
         tCopula = {
           psn.sample = qt(x, df = copula@parameters[2])
           l.t = log(.t.12.dens(copula@parameters[1], psn.sample, nu = copula@parameters[2]))
           l.t.b = 0
           if (param.est == T){
             for(b in 1:B){
               l.t.b = c(l.t.b, .t.12.dens(fitCopula(tCopula(dim = dims, df = copula@parameters[2], df.fixed = T), data = psn.sample[-(((b-1)*m+1):(b*m)),], method = "mpl")@estimate, psn.sample[((b-1)*m+1):(b*m),], nu = copula@parameters[2]))
             }
           } else {
             for(b in 1:B){
               l.t.b = c(l.t.b, .t.12.dens(copula@parameters[1], psn.sample[((b-1)*m+1):(b*m),], nu = copula@parameters[2]))
             }
           }
           l.t.b = l.t.b[-1]
           stat = sum(l.t - log(l.t.b)) - 1
           stat
         },
         normalCopula = {
           psn.sample = qnorm(x)
           l.ac.ga.b = as.vector(sapply(1:B, FUN = .opt.ga, psn.sample = psn.sample, m = m, dims = dims))
           stat = sum(log(dCopula(x, normalCopula(copula@parameters, dims))) - log(l.ac.ga.b)) - 1 # change here an estimation of the copula!!!!!
           stat
         })
  
}

.Kernel = function(x, copula, dims, n, nodes.Integration, MJ, delta.J) {
  switch(class(copula),
         claytonCopula = {
           Int_Grid = createIntegrationGrid("GQU", dimension = dims, k = nodes.Integration)
           bootsample = rCopula(MJ, claytonCopula(copula@parameters, dim = dims))
           h = as.vector((diag(2.6073*n^(-1/6)*chol(cov(x))) * delta.J))
           Jn_c = c()
           for(i in 1:dim(Int_Grid$nodes)[1]){
             Jn_c = c(Jn_c, Int_Grid$weights[i] * .integrand(Int_Grid$nodes[i,], x, Lbootsample = bootsample, h))
           }
           stat = sum(Jn_c)
           stat
         },
         gumbelCopula = {
           Int_Grid = createIntegrationGrid("GQU", dimension = dims, k = nodes.Integration)
           bootsample = rCopula(MJ, gumbelCopula(copula@parameters, dim = dims))
           h = as.vector((diag(2.6073*n^(-1/6)*chol(cov(x))) * delta.J))
           Jn_c = c()
           for(i in 1:dim(Int_Grid$nodes)[1]){
             Jn_c = c(Jn_c, Int_Grid$weights[i] * .integrand(Int_Grid$nodes[i,], x, Lbootsample = bootsample, h))
           }
           stat = sum(Jn_c)
           stat
         },
         frankCopula = {
           Int_Grid = createIntegrationGrid("GQU", dimension = dims, k = nodes.Integration)
           bootsample = rCopula(MJ, frankCopula(copula@parameters, dim = dims))
           h = as.vector((diag(2.6073*n^(-1/6)*chol(cov(x))) * delta.J))
           Jn_f = c()
           for(i in 1:dim(Int_Grid$nodes)[1]){
             Jn_f = c(Jn_f, Int_Grid$weights[i] * .integrand(Int_Grid$nodes[i,], x, Lbootsample = bootsample, h))
           }
           stat = sum(Jn_f)
           stat
         },
         tCopula = {
           psn.sample = qt(x, df = copula@parameters[2])
           Int_Grid = createIntegrationGrid("GQU", dimension = dims, k = nodes.Integration)
           bootsample = rCopula(MJ, tCopula(copula@parameters[1], df = copula@parameters[2], dim = dims))
           h = as.vector((diag(2.6073*n^(-1/6)*chol(cov(x))) * delta.J))
           Jn_c = c()
           for(i in 1:dim(Int_Grid$nodes)[1]){
             Jn_c = c(Jn_c, Int_Grid$weights[i] * .integrand(Int_Grid$nodes[i,], x, Lbootsample = bootsample, h))
           }
           stat = sum(Jn_c)
           stat
         },
         normalCopula = {
           Int_Grid = createIntegrationGrid("GQU", dimension = dims, k = nodes.Integration)
           bootsample = rCopula(MJ, normalCopula(copula@parameters, dim = dims))
           h = as.vector((diag(2.6073*n^(-1/6)*chol(cov(x))) * delta.J))
           Jn_c = c()
           for(i in 1:dim(Int_Grid$nodes)[1]){
             Jn_c = c(Jn_c, Int_Grid$weights[i] * .integrand(Int_Grid$nodes[i,], x, Lbootsample = bootsample, h))
           }
           stat = sum(Jn_c)
           stat
         })
  
}