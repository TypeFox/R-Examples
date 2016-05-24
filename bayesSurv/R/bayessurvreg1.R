###########################################
#### AUTHOR:     Arnost Komarek        ####
####             (2004)                ####
####                                   ####
#### FILE:       bayessurvreg1.R       ####
####                                   ####
#### FUNCTIONS:  bayessurvreg1         ####
###########################################

### ======================================
### bayessurvreg1
### ======================================
bayessurvreg1 <- function(
     formula,
     random,
     data = parent.frame(),
     subset,
     na.action = na.fail,
     x = FALSE,
     y = FALSE,                          
     onlyX = FALSE,
     nsimul = list(niter = 10, nthin = 1, nburn = 0, nnoadapt = 0, nwrite = 10),
     prior = list(kmax = 5, k.prior = "poisson", poisson.k = 3,
                  dirichlet.w = 1,
                  mean.mu = NULL, var.mu = NULL,
                  shape.invsig2 = 1.5, shape.hyper.invsig2 = 0.8, rate.hyper.invsig2 = NULL,
                  pi.split = NULL, pi.birth = NULL,
                  Eb0.depend.mix = FALSE),
     prior.beta,
     prior.b,
     prop.revjump,
     init = list(iter = 0, mixture = NULL, beta = NULL, b = NULL, D = NULL,
                 y = NULL, r = NULL, otherp = NULL, u = NULL),
     store = list(y = TRUE, r = TRUE, b = TRUE, u = TRUE, MHb = FALSE, regresres = FALSE),
     dir = getwd(),
     toler.chol = 1e-10,
     toler.qr = 1e-10,
     ...)
{  
   thispackage = "bayesSurv"
   #thispackage = NULL
  
   transform = function(t){log(t)}
   dtransform = function(t){1/t}
  
   sim.to.R <- FALSE                  ## this is here for a compatibility with an older code
                                      ## If you try to change sim.to.R into TRUE, you also need to redefine 'row.need' variable below 
   store <- bayessurvreg1.checkStore(store)
   nsimul <- bayessurvreg.checknsimul(nsimul)
   
   ## Give a function call to be recorded in a resulting object.
   call <- match.call(expand.dots = TRUE)

   ## Extract all the design information from the function call
   m <- match.call(expand.dots = FALSE)
   des <- bayessurvreg.design(m, formula, random, data, transform, dtransform)
   if (onlyX) return (des$X)   

   
   ## =========================================================
   ## Manipulate with initial values and the prior information
   ## =========================================================
   priordi <- bayessurvreg1.priorInit(prior, init, des$Yinit, des$Xinit, des$n, des$nX, des$nrandom, des$ncluster, des$indb, des$randomInt, toler.chol)
   prior <- attr(priordi, "prior")
   init <- attr(priordi, "init")
   
   if (missing(prop.revjump)) prop.revjump <- list()
   revjumpdi <- bayessurvreg1.revjump(prop.revjump)
   prop.revjump <- attr(revjumpdi, "prop.revjump")
   
   if (missing(prior.beta)) prior.beta <- list()
   betadi <- bayessurvreg1.priorBeta(prior.beta, des$nX, des$indb, des$factors, des$n.factors, des$n.in.factors)
   prior.beta <- attr(betadi, "prior.beta")   
   prior.beta.noadapt <- attr(betadi, "prior.beta")
   if (des$nX){     
     adapts <- pmatch(prior.beta$type.upd, table = c("adaptive.metropolis"), nomatch = 0, duplicates.ok = TRUE)       
     n.adapt <- sum(adapts)
     prior.beta.noadapt$type.upd[adapts == 1] <- "random.walk.metropolis"
   }
   else{
     n.adapt <- 0
   }       
   betadi.noadapt <- bayessurvreg1.priorBeta(prior.beta.noadapt, des$nX, des$indb, des$factors, des$n.factors, des$n.in.factors) 
   
   if (missing(prior.b)) prior.b <- list()
   bdi <- bayessurvreg1.priorb(prior.b, des$nrandom, des$ncluster, toler.chol)
   prior.b <- attr(bdi, "prior.b")
   

   ## ===================================================================
   ## Compute quantities to determine the space needed to be allocated
   ##   and numbers of iterations in different phases
   ## ===================================================================
   if (nsimul$nburn >= nsimul$niter) nsimul$nburn <- nsimul$niter - 1
   if (nsimul$nburn < 0) nsimul$nburn <- 0

   if (n.adapt == 0) nsimul$nnoadapt <- nsimul$nburn
   if (nsimul$nnoadapt > nsimul$nburn) nsimul$nnoadapt <- nsimul$nburn
   if (nsimul$nnoadapt < 0) nsimul$nnoadapt <- 0
   
   if (nsimul$nburn == 0) nruns <- 1
   else                   if (nsimul$nnoadapt == nsimul$nburn | nsimul$nnoadapt == 0) nruns <- 2
                          else                                                        nruns <- 3

   nrun <- numeric(3)
   nrun[3] <- nsimul$niter - nsimul$nburn
   nrun[2] <- nsimul$nburn - nsimul$nnoadapt
   nrun[1] <- nsimul$nnoadapt

   nwrite.run <- nrun
   nwrite.run[nsimul$nwrite <= nrun] <- nsimul$nwrite   
   max.nwrite <- max(nwrite.run)

   if (!des$nrandom){ store$b <- FALSE;  store$MHb <- FALSE}
   #row.need <- ifelse(sim.to.R, max(nsimul$nburn, nafterburn), max.nwrite)
   # I do not know any more what 'nafterburn' should be, so change the above row into row.need <- max.nwrite
   row.need <- max.nwrite
   
   ## =====================================================================================
   ## Write headers to files with stored values
   ## =====================================================================================
   bayessurvreg1.writeHeaders(dir, prior, store, des$nX, des$X, des$names.random, des$ncluster, des$nrandom, des$rnamesX,
                              unique(des$cluster), betadi$integer[1], bdi$integer[4])

   ## =========================================================
   ## Combine similar parameters into one vector
   ## =========================================================
   dims <- c(des$n, des$ncluster, des$nwithin, des$nY, des$nX, des$nfixed, des$nrandom, 1*des$randomInt, row.need)
   storeV <- c(store$y, store$r, store$b, store$u, store$MHb, store$regresres)
   nsimul.run1 <- c(nrun[1], nsimul$nthin, nwrite.run[1])
   nsimul.run2 <- c(nrun[2], nsimul$nthin, nwrite.run[2])
   nsimul.run3 <- c(nrun[3], nsimul$nthin, nwrite.run[3])   
   tolers <- c(toler.chol, toler.qr)

   ## =====================================
   ## Keep some parameters to be returned
   ## =====================================
   keep.init <- init
   
   cat("Simulation started on                       ", date(), "\n", sep = "")      
     ## Run without adaptation of a proposal covariance matrices
     ## Either the whole burn up or first part of burn up
   if (nruns == 3 | (nruns == 2 & nsimul$nnoadapt == nsimul$nburn)){
     fit <- .C("bayessurvreg1", as.character(dir),
                                dims = as.integer(dims),               
                                Y = as.double(des$Y),
                                X = as.double(des$X),
                                indb = as.integer(des$indb),               
                                iter = as.integer(init$iter),
                                loglik = as.double(c(0, 0)),
                                mixture = as.double(init$mixture),
                                mixmoment = as.double(c(0, 0)),
                                beta = as.double(init$beta),
                                b = as.double(init$b),
                                D = as.double(init$D),
                                r = as.integer(init$r),
                                Ys = as.double(init$y),       
                                otherp = as.double(init$otherp),
                                u = as.double(init$u),
                                prior.pari = as.integer(priordi$integer),
                                prior.pard = as.double(priordi$double),
                                revJump.pari = as.integer(revjumpdi$integer),
                                revJump.pard = as.double(revjumpdi$double),
                                prior.betai = as.integer(betadi.noadapt$integer),
                                prior.betad = as.double(betadi.noadapt$double),
                                prior.bi = as.integer(bdi$integer),
                                prior.bd = as.double(bdi$double),
                                nsimul = as.integer(nsimul.run1),
                                store = as.integer(storeV),
                                tolers = as.double(tolers),
                                err = integer(1),
               PACKAGE = thispackage)     
     if (fit$err != 0) stop ("Something went wrong during the simulation.")
     cat("Simulation without adaptation finished on   ", date(), "   (iteration ", fit$iter, ")", "\n", sep = "")
   
       ## Recalculate intial proposal covariance matrices
       ##  * do not change it if sample covariance matrix is not positive definite
     if (des$nX){
       betasam <- read.table(paste(dir, "/beta.sim", sep = ""), header = TRUE)
       betaaver <- apply(betasam, 2, mean)
       prior.beta$mean.sampled <- betaaver
       to.adapt <- (1:betadi$integer[1])[adapts == 1]
       if (length(to.adapt) > 0){
         for (i in 1:length(to.adapt)){
           colsX <- prior.beta$blocks$ind.block[[to.adapt[i]]]
           sample <- as.data.frame(betasam[, colsX])
           colnames(sample) <- colnames(betasam)[colsX]
           covmat <- sampleCovMat(sample)
           covmat <- prior.beta$sd.AM[length(colsX)] * (covmat + prior.beta$eps.AM[to.adapt[i]]*diag(length(colsX)))
           vind <- 0:(dim(covmat)[1] - 1)
           diagI <- (vind * (2*dim(covmat)[1] - vind + 1)) / 2
           covmatt <- covmat[lower.tri(covmat, diag = TRUE)]
           chol <- .C("cholesky", A = as.double(covmatt), rank = integer(1), as.integer(dim(covmat)[1]),
                                    as.integer(diagI), as.double(toler.chol),
                      PACKAGE = thispackage)
           if (chol$rank == dim(covmat)[1]){
             prior.beta$blocks$cov.prop[[to.adapt[i]]] <- as.numeric(covmatt)
           } 
           else{
             warning("Sample covariance matrix after no adapt period was not positive definite.")
           }
         }
       }
       betadi <- bayessurvreg1.priorBeta(prior.beta, des$nX, des$indb, des$factors, des$n.factors, des$n.in.factors)
     }       
     
       ## Give new initials
     init$iter <- fit$iter;   init$mixture <- fit$mixture;   init$beta <- fit$beta;
     init$b <- fit$b;         init$D <- fit$D;               init$r <- fit$r;
     init$y <- fit$Ys;        init$otherp <- fit$otherp;     init$u <- fit$u;

       ## Rewrite sampled values by new files
     bayessurvreg1.writeHeaders(dir, prior, store, des$nX, des$X, des$names.random, des$ncluster, des$nrandom, des$rnamesX,
                                unique(des$cluster), betadi$integer[1], bdi$integer[4])
   }     

     ## Burn up with adaptation
     ## Either the whole burn up or the second part of burn up
   if (nruns == 3 | (nruns == 2 & nsimul$nnoadapt == 0)){
     fit <- .C("bayessurvreg1", as.character(dir),
                                dims = as.integer(dims),               
                                Y = as.double(des$Y),
                                X = as.double(des$X),
                                indb = as.integer(des$indb),               
                                iter = as.integer(init$iter),
                                loglik = as.double(c(0, 0)),
                                mixture = as.double(init$mixture),
                                mixmoment = as.double(c(0, 0)),               
                                beta = as.double(init$beta),
                                b = as.double(init$b),
                                D = as.double(init$D),
                                r = as.integer(init$r),
                                Ys = as.double(init$y),       
                                otherp = as.double(init$otherp),
                                u = as.double(init$u),
                                prior.pari = as.integer(priordi$integer),
                                prior.pard = as.double(priordi$double),
                                revJump.pari = as.integer(revjumpdi$integer),
                                revJump.pard = as.double(revjumpdi$double),
                                prior.betai = as.integer(betadi$integer),
                                prior.betad = as.double(betadi$double),
                                prior.bi = as.integer(bdi$integer),
                                prior.bd = as.double(bdi$double),
                                nsimul = as.integer(nsimul.run2),
                                store = as.integer(storeV),
                                tolers = as.double(tolers),
                                err = integer(1),
              PACKAGE = thispackage)
     if (fit$err != 0) stop ("Something went wrong during the simulation.")
     cat("Burn-up finished on                         ", date(), "   (iteration ", fit$iter, ")", "\n", sep = "")
     
       ## Give new initials
     init$iter <- fit$iter;   init$mixture <- fit$mixture;   init$beta <- fit$beta;
     init$b <- fit$b;         init$D <- fit$D;               init$r <- fit$r;
     init$y <- fit$Ys;        init$otherp <- fit$otherp;     init$u <- fit$u;     

       ## Rewrite sampled values by new files
     bayessurvreg1.writeHeaders(dir, prior, store, des$nX, des$X, des$names.random, des$ncluster, des$nrandom, des$rnamesX,
                                unique(des$cluster), betadi$integer[1], bdi$integer[4])
   }     
   
     ## Main simulation
   fit <- .C("bayessurvreg1", as.character(dir),
                              dims = as.integer(dims),               
                              Y = as.double(des$Y),
                              X = as.double(des$X),
                              indb = as.integer(des$indb),               
                              iter = as.integer(init$iter),
                              loglik = as.double(c(0, 0)),
                              mixture = as.double(init$mixture),
                              mixmoment = as.double(c(0, 0)),             
                              beta = as.double(init$beta),
                              b = as.double(init$b),
                              D = as.double(init$D),
                              r = as.integer(init$r),
                              Ys = as.double(init$y),       
                              otherp = as.double(init$otherp),
                              u = as.double(init$u),
                              prior.pari = as.integer(priordi$integer),
                              prior.pard = as.double(priordi$double),
                              revJump.pari = as.integer(revjumpdi$integer),
                              revJump.pard = as.double(revjumpdi$double),
                              prior.betai = as.integer(betadi$integer),
                              prior.betad = as.double(betadi$double),
                              prior.bi = as.integer(bdi$integer),
                              prior.bd = as.double(bdi$double),
                              nsimul = as.integer(nsimul.run3),
                              store = as.integer(storeV),
                              tolers = as.double(tolers),
                              err = integer(1),
            PACKAGE = thispackage)
   if (fit$err != 0) stop ("Something went wrong during the simulation.")
   cat("Simulation finished on                      ", date(), "   (iteration ", fit$iter, ")", "\n", sep = "")   

   toreturn <- fit$iter
   attr(toreturn, "call") <- call
   attr(toreturn, "prior") <- attr(priordi, "prior")
   attr(toreturn, "init") <- attr(priordi, "init")
   attr(toreturn, "prop.revjump") <- attr(revjumpdi, "prop.revjump")
   attr(toreturn, "prior.beta") <- attr(betadi, "prior.beta")
   attr(toreturn, "prior.b") <- attr(bdi, "prior.b")
   if (x) attr(toreturn, "x") <- des$X
   if (y) attr(toreturn, "y") <- des$Y
   class(toreturn) <- "bayessurvreg1"
   
   return(toreturn)
}


