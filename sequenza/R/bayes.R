mufreq.dbinom <- function(mufreq, mufreq.model, depth.t, seq.errors = 0.01, ...) {
   mufreq.model[mufreq.model == 0] <- seq.errors
   n.success       <- round(mufreq * depth.t, 0)
   dbinom( x = n.success, size = depth.t, prob = mufreq.model, ...)
}

mufreq.dpois <- function(mufreq, mufreq.model, depth.t, seq.errors = 0.01, ...) {
   mufreq.model[mufreq.model == 0] <- seq.errors
   n.success       <- round(mufreq * depth.t, 0)
   dpois( x = n.success, lambda = mufreq.model * depth.t, ...)
}

baf.dbinom <- function(baf, baf.model, depth.t, ...) {
   n.success       <- round(baf * depth.t, 0)
   dbinom( x = n.success, size = depth.t, prob = baf.model, ...)
}

baf.dpois <- function(baf, baf.model, depth.t, ...) {
   n.success       <- round(baf * depth.t, 0)
   dpois( x = n.success, lambda = baf.model * depth.t, ...)
}

depth.ratio.dbinom <- function(size, depth.ratio, depth.ratio.model, ...) {
   #n.success        <- round(depth.n * depth.ratio, 0)
   n.success        <- round(size * (depth.ratio/(1 + depth.ratio)), 0)
   prob             <- depth.ratio.model / (1 + depth.ratio.model)
   dbinom( x = n.success, size = size, prob = prob, ...)
}

depth.ratio.dpois <- function(size, depth.ratio, depth.ratio.model, ...) {
   x        <- round(size * depth.ratio, 0)
   dpois( x = x, lambda = depth.ratio.model * size, ...)
}

mufreq.bayes <- function(mufreq, depth.ratio, cellularity, ploidy, avg.depth.ratio,
                         weight.mufreq = 100, weight.ratio = 100, CNt.min = 1, CNt.max = 7,
                         CNn = 2, priors.table = data.frame(CN = CNt.min:CNt.max, value = 1)) {

   mufreq.tab <- data.frame(F = mufreq, ratio = depth.ratio,
                            weight.mufreq = weight.mufreq, weight.ratio = weight.ratio)
   types <- types.matrix(CNt.min = CNt.min, CNt.max = CNt.max, CNn = CNn)
   mufreq.depth.ratio <- cbind(types, model.points(cellularity = cellularity, ploidy = ploidy,
                                                   types = types, avg.depth.ratio = avg.depth.ratio))
   rows.x             <- 1:nrow(mufreq.tab)

   priors <- rep(1, nrow(mufreq.depth.ratio))
   for (i in 1:nrow(priors.table)) {
      priors[mufreq.depth.ratio$CNt == priors.table$CN[i]] <- priors.table$value[i]
   }
   priors <- priors / sum(priors)

   bayes.fit <- function (x, mat, model.pts, priors) {
      test.ratio <- model.pts$depth.ratio
      test.mufrq <- model.pts$mufreqs
      min.offset <- 1e-323
      score.r    <- depth.ratio.dbinom(size = mat[x,]$weight.ratio, depth.ratio = mat[x,]$ratio, test.ratio)
      score.m    <- mufreq.dbinom(mufreq = mat[x,]$F, depth.t = mat[x,]$weight.mufreq, test.mufrq)

      score.r    <- score.r * priors
      score.m    <- score.m

      post.model <- score.r * score.m

      post.model[post.model == 0] <- min.offset
      # max.lik <-  which.max(post.model)
      # max.post <- c(as.numeric(model.pts[max.lik,1:3]), log2(post.model[max.lik]))
      # max.post

      res.cn     <- model.pts$CNt[which.max(score.r)]
      idx.pts    <- model.pts$CNt == res.cn
      model.lik  <- cbind(model.pts[idx.pts, 1:3], log(post.model[idx.pts]))
      if (is.null(dim(model.lik))) {
         max.post <- model.lik
      } else {
         max.post   <- model.lik[which.max(model.lik[,4]),]
      }

      max.post
   }
   types.L           <- mapply(FUN = bayes.fit, rows.x,
                         MoreArgs = list(mat = mufreq.tab,
                                         model.pts = mufreq.depth.ratio,
                                         priors = priors),
                                         SIMPLIFY = FALSE)
   types.L           <- do.call(rbind, types.L)
   colnames(types.L) <- c("CNn","CNt","Mt", "LPP")
   types.L
}

baf.bayes <- function(Bf, depth.ratio, cellularity, ploidy, avg.depth.ratio,
                      sd.Bf = 0.1, sd.ratio = 0.5, weight.Bf = 1, weight.ratio = 1, CNt.min = 0,
                      CNt.max = 7, CNn = 2, priors.table = data.frame(CN = CNt.min:CNt.max,
                      value = 1), ratio.priority = FALSE) {

   mufreq.tab <- data.frame(Bf = Bf, ratio = log(depth.ratio),
                            sd.Bf = sd.Bf, sd.ratio = sd.ratio,
                            weight.Bf = weight.Bf,
                            weight.ratio = weight.ratio)
   mufreq.depth.ratio <- model.points(cellularity = cellularity, ploidy = ploidy,
                                      types = cbind(CNn = CNn, CNt = CNt.min:CNt.max, Mt = 0),
                                      avg.depth.ratio = avg.depth.ratio)
   model.d.ratio      <- cbind(CNt = CNt.min:CNt.max, depth.ratio = log(mufreq.depth.ratio[, 2]))
   model.baf          <- expected.baf(sd = mean(sd.Bf[Bf > 0], na.rm = TRUE), CNn = CNn, CNt = CNt.max, cellularity = cellularity)
   #model.baf          <- theoretical.baf(CNn = CNn, CNt = CNt.max, cellularity = cellularity)
   # B-allele freq are never 0.5, always smaller. just a work around on this.. to be better fixed!
   #model.baf$BAF[model.baf$A==model.baf$B] <- quantile(rep(mufreq.tab$Bf, times = mufreq.tab$N.Bf),
   #                                                    na.rm = TRUE, probs = 0.95)
   if(CNt.min == 0) {
     model.baf          <- as.data.frame(rbind(c(0, 0, max(model.baf$BAF), 0), model.baf))
   }
   #                                                na.rm = TRUE, probs = skew.baf)
   model.pts          <- merge(model.baf, model.d.ratio)
   # model.pts          <- cbind(baf.type = apply(model.pts[, 1:3], 1, FUN = function(x) paste(x, collapse = "_")),
   #                             model.pts[, 4:5])
   rows.x             <- 1:nrow(mufreq.tab)

   priors <- rep(1, nrow(model.pts))
   for (i in 1:nrow(priors.table)) {
      priors[model.pts$CNt == priors.table$CN[i]] <- priors.table$value[i]
   }
   priors <- priors / sum(priors)

   bayes.fit <- function (x, mat, model.pts, priors, ratio.priority) {
      test.ratio <- model.pts$depth.ratio
      test.baf   <- model.pts$BAF
      min.offset <- 1e-323
      #score.r    <- depth.ratio.dbinom(size = mat[x,]$sd.ratio, depth.ratio = mat[x,]$ratio, test.ratio)
      score.r    <- dt2(sd = mat[x,]$sd.ratio/sqrt(mat[x,]$weight.ratio), mean = mat[x,]$ratio, x = test.ratio, df = 5, log = TRUE)
      score.r    <- score.r + log(priors)
      if (!is.na(mat[x,]$Bf) & !is.na(mat[x,]$sd.Bf/sqrt(mat[x,]$weight.Bf))) {
         #score.b    <- baf.dpois(baf = mat[x,]$Bf, depth.t = mat[x,]$N.Bf, test.baf, log = TRUE)
         score.b    <- dt2(mean = mat[x,]$Bf, sd = mat[x,]$sd.Bf/sqrt(mat[x,]$weight.Bf), x = test.baf, df = 5, log = TRUE)
         post.model <- score.r + score.b
      } else {
         post.model <- score.r         
      }

      #post.model[post.model == 0] <- min.offset
      if (ratio.priority == FALSE) {
         max.lik <-  which.max(post.model)
         max.post <- c(as.numeric(model.pts[max.lik,1:3]), post.model[max.lik])
      } else {
         res.cn     <- model.pts$CNt[which.max(score.r)]
         idx.pts    <- model.pts$CNt == res.cn
         model.lik  <- cbind(model.pts[idx.pts, 1:3], post.model[idx.pts])
         if (is.null(dim(model.lik))) {
            max.post <- model.lik
         } else {
            max.post   <- model.lik[which.max(model.lik[,4]),]
         }
      }
      if (is.na(mat[x,]$Bf)) {
         max.post[2:3] <- NA
      }
      max.post
   }
   bafs.L           <- mapply(FUN = bayes.fit, rows.x,
                         MoreArgs = list(mat = mufreq.tab,
                                         model.pts = model.pts,
                                         priors = priors,
                                         ratio.priority = ratio.priority),
                                         SIMPLIFY = FALSE)
   bafs.L           <- do.call(rbind, bafs.L)
   colnames(bafs.L) <- c("CNt", "A", "B", "LPP")
   bafs.L
}

mufreq.model.fit <- function(cellularity = seq(0.3, 1, by = 0.01),
                          ploidy = seq(1, 7, by = 0.1),
                          mc.cores = getOption("mc.cores", 2L), ...) {

   result <- expand.grid(ploidy = ploidy, cellularity = cellularity,
                         KEEP.OUT.ATTRS = FALSE)

   fit.cp <- function(ii) {
      L.model <- mufreq.bayes(cellularity = result$cellularity[ii],
                           ploidy = result$ploidy[ii], ...)
      sum(L.model[,4])
   }
   bayes.res <- mclapplyPb(X = 1:nrow(result), FUN = fit.cp, mc.cores = mc.cores)
   result$LPP <- unlist(bayes.res)
   z <- tapply(result$LPP, list(result$ploidy, result$cellularity), mean)
   x <- as.numeric(rownames(z))
   y <- as.numeric(colnames(z))
   max.lik <- max(result$LPP, na.rm = TRUE)
   LogSumLik <- log(sum(exp(result$LPP - max.lik))) + max.lik
   znorm <- exp(z - LogSumLik)
   list(ploidy = x, cellularity = y, lpp = znorm)
}

baf.model.fit <- function(cellularity = seq(0.3, 1, by = 0.01),
                          ploidy = seq(1, 7, by = 0.1),
                          mc.cores = getOption("mc.cores", 2L), ...) {

   result <- expand.grid(ploidy = ploidy, cellularity = cellularity,
                         KEEP.OUT.ATTRS = FALSE)

   fit.cp <- function(ii) {
      L.model <- baf.bayes(cellularity = result$cellularity[ii],
                           ploidy = result$ploidy[ii], ...)
      sum(L.model[,4])
   }
   bayes.res <- mclapplyPb(X = 1:nrow(result), FUN = fit.cp, mc.cores = mc.cores)
   result$LPP <- unlist(bayes.res)
   z <- tapply(result$LPP, list(result$ploidy, result$cellularity), mean)
   x <- as.numeric(rownames(z))
   y <- as.numeric(colnames(z))
   max.lik <- max(result$LPP, na.rm = TRUE)
   LogSumLik <- log(sum(exp(result$LPP - max.lik))) + max.lik
   znorm <- exp(z - LogSumLik)
   list(ploidy = x, cellularity = y, lpp = znorm)
}
