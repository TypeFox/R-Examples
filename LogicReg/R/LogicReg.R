# =====================================================================
# Copyright (C) 1999-2008  Ingo Ruczinski and Charles Kooperberg

# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# The text of the GNU General Public License, version 2, is available
# as http://www.gnu.org/copyleft or by writing to the Free Software
# Foundation, 59 Temple Place - Suite 330, Boston, MA 02111-1307, USA.

# The main reference for Logic Regression is
# Ruczinski I, Kooperberg C, LeBlanc M (2003). Logic Regression,
# Journal of Computational and Graphical Statistics, 12, 475-511.
# Other references can be found on our homepages
# http://www.biostat.jhsph.edu/~iruczins/
# http://kooperberg.fhcrc.org/~clk
# You can contact us at ingo@jhu.edu and clk@fredhutch.org
# =======================================================================

eval.logreg <- function(ltree, data)
{
   if(class(ltree) == "logregtree") {
      ltree <- ltree$trees
      nn <- matrix(1, ncol = length(ltree[, 1]), nrow = length(data[ , 1]))
      for(i in 1:length(ltree[, 1])) {
         if(ltree[i, 2] == 3) {
            nn[, i] <- data[, ltree[i, 3]]
            if(ltree[i, 4] == 1)
               nn[, i] <- 1 - nn[, i]
         }
      }
      ox <- length(ltree[, 1])
      ox <- floor(ox/2)
      for(i in ox:1) {
         if(ltree[i, 2] == 1)
            nn[, i] <- nn[, 2 * i] * nn[, 2 * i + 1]
         if(ltree[i, 2] == 2)
            nn[, i] <- nn[, 2 * i] + nn[, 2 * i + 1]
      }
      nn[nn > 1] <- 1
      if(sum(ltree != 0) == 0)
         nn[, 1] <- 0
      nn[, 1]
   }
   else {
      if(class(ltree) == "logregmodel") {
         n1 <- eval.logreg(ltree$trees[[1]], data)
         if(ltree$ntrees[1] > 1)
            for(i in 2:ltree$ntrees[1])
               n1 <- cbind(n1, eval.logreg(ltree$trees[[i]], data))
         n1
      }
      else {
         stop("ltree not of class logregtree of logregmodel")
      }
   }
}
frame.logreg <- function(fit, msz, ntr, newbin, newresp, newsep, newcens, newweight)
{
   if(class(fit) != "logreg")
      stop("fit not of class logreg")
   if(missing(newbin)) {
      outframe <- data.frame(y = fit$response)
      outframe <- data.frame(outframe, wgt = fit$weight)
      if(fit$type == "proportional.hazards" || fit$type == "exponential.survival")
         outframe <- data.frame(outframe, cens = fit$censor)
      if(fit$nsep > 0)
         outframe <- data.frame(outframe, fit$separate)
      binhere <- fit$binary
   }
   else {
      binhere <- newbin
      lbinhere <- length(binhere)
      if(is.data.frame(binhere))
         binhere <- as.matrix(binhere)
      if(is.matrix(binhere) == FALSE)
         binhere <- matrix(binhere, ncol = fit$nbinary)
      n1 <- length(binhere[, 1])
      n2 <- length(binhere[1,  ])
      if(n2 != fit$nbinary)
         stop("new number of binary predictors doesn't match fit")
      if(!missing(newweight)) {
         if(length(newweight) != n1)
            stop("length(newweight) != length(newbin[,1])")
      }
      else {
         newweight <- rep(1, n1)
      }
      if(!missing(newresp)) {
         if(length(newresp) != n1)
            stop("length(newresp) != length(newbin[,1])")
         outframe <- data.frame(y = newresp, wgt = newweight)
      }
      else {
         outframe <- data.frame(wgt = newweight)
      }
      if(fit$type == "proportional.hazards" || fit$type == "exponential.survival"){
         if(missing(newcens)) {
            warning("newcens missing, taking all censoring indicators to be 1")
            outframe <- data.frame(outframe, cens = rep(1, n1))
         }
         else {
            if(length(newcens) != n1)
               stop("length(newcens) != length(newbin[,1])")
            outframe <- data.frame(outframe, cens = newcens)
         }
      }
      if(fit$nsep > 0) {
         if(missing(newsep))
            stop("you need to specify newsep")
         if(is.matrix(newsep) == FALSE)
            newsep <- matrix(newsep, ncol = fit$nsep)
         if(is.data.frame(newsep))
            newsep <- as.matrix(newsep)
         if(length(newsep[, 1]) != n1)
            stop("length(newsep[,1]) != length(newbin[,1])")
         if(length(newsep[1,  ]) != fit$nsep)
            stop("new number of separate predictors doesn't match fit")
         outframe <- data.frame(outframe, newsep)
      }
   }
   if(fit$choice!=1 && fit$choice!=2 && fit$choice!=6)
      stop("fit$choice needs to be 1, 2, or 6")
   if(fit$choice == 1) {
      ntr <- fit$ntrees[1]
      msz <- fit$nleaves[1]
      for(j in 1:ntr) {
         mtree <- fit$model$trees[[j]]$trees
         if(msz < 0)
            msz <- 0.5 * (length(mtree[, 1]) + 1)
         if(mtree[1, 2] == 0)
            outframe <- data.frame(outframe, tmp = 0)
         else {
            tmp <- eval.logreg(fit$model$trees[[j]], binhere)
            outframe <- data.frame(outframe, tmp = tmp)
         }
         l1 <- length(outframe[1,  ])
         names(outframe)[l1] <- paste("tree", msz, ".", ntr, ".", j, sep = "")
      }
   }
   else {
      for(i in 1:fit$nmodels) {
         xfit <- fit$alltrees[[i]]
         ntrx <- xfit$ntrees[1]
         mszx <- xfit$nleaves[1]
         l1 <- 1
         if(!missing(ntr) && fit$choice == 2)
            l1 <- sum(ntrx == ntr)
         if(!missing(msz))
            l1 <- l1 * sum(mszx == msz)
         if(l1 > 0)
            for(j in 1:ntrx) {
               if(xfit$trees[[j]]$trees[1, 2] == 0)
                  outframe <- data.frame(outframe, tmp = 0)
               else {
                  tmp <- eval.logreg(xfit$trees[[j]], binhere)
                  outframe <- data.frame(outframe, tmp = tmp)
               }
               l1 <- length(outframe[1,  ])
               names(outframe)[l1] <- paste("tree", mszx, ".", ntrx, ".", j, sep = "")
            }
      }
   }
   outframe
}
logreg.anneal.control <- function(start = 0, end = 0, iter = 0, earlyout = 0, update = 0)
{
   if(is.list(start)) {
      if(length(start$end) > 0)
         end <- start$end
      if(length(start$iter) > 0)
         iter <- start$iter
      if(length(start$earlyout) > 0)
         earlyout <- start$earlyout
      if(length(start$update) > 0)
         update <- start$update
      if(length(start$start) > 0)
         start <- start$start
      else start <- 0
   }
   if(iter < 100 && (start != 0 || end != 0))
      stop(" not enough repetitions ")
   if(start < end)
      stop("starting temperature below ending temperature")
   list(start = start, end = end, iter = iter, earlyout = earlyout, update = update)
}
logreg.myown <- function()
{
   help(logreg.myown)
}
logreg.tree.control <- function(treesize = 8, opers = 1, minmass = 0, n1)
{
   if(is.list(treesize)) {
      if(length(treesize$opers) > 0)
         opers <- treesize$opers
      if(length(treesize$minmass) > 0)
         minmass <- treesize$minmass
      if(length(treesize$treesize) > 0)
         treesize <- treesize$treesize
      else treesize <- 16
   }
   treesize <- 2^floor(logb(abs(treesize + 0.0001), 2))
   if(opers != 2 && opers != 3 && opers!=4)
      opers <- 1
   if(minmass < 0)
      minmass <- 0
   if(!missing(n1))
      if(minmass > n1/4)
         stop("minmass should be at most samplesize/4")
   list(treesize = treesize, opers = opers, minmass = minmass)
}
logregmodel <- function()
{
   help(logregmodel)
}
logregtree <- function()
{
   help(logregtree)
}
plot.logreg <- function(x, pscript = FALSE, title = TRUE, ...)
{
   if(class(x) != "logreg")
      stop("x not of class logreg")
   logregtmp <- hist(1:10, plot = FALSE)
   if(length(logregtmp$mids) > 0)
      logregwhite <- "white"
   if(length(logregtmp$mids) == 0)
      logregwhite <- 0
   nmsx <- 0
   if(length(x$binnames) > 0) {
      nms <- x$binnames
      if(nms[1] != "V1" && nms[1] != "X1")
         nmsx <- 1
   }
   if(x$choice == 7){
      barplot(x$size[,2]/(sum(x$size[,2]))*100,names=as.character(x$size[,1]),
              main="size distribution of visited models (%)",...)
      if(length(x$single)>0){
         tmp <- c(x$single)
         if(length(x$binnames)==0)
            nms <- 1:length(tmp)
         else
            nms <- x$binnames
         barplot(tmp/(sum(x$size[,2]))*100,names=as.character(nms),
              main="marginal frequency of being in the model (%)",...)
         if(length(x$double)>0){
            tmp2 <- x$double
            tmp3 <- 1:length(tmp2[,1])
            tmp4 <- outer(tmp3,tmp3,"-")
            tmp2[tmp4<=0] <- NA
            image(tmp3,tmp3,tmp2,main="frequency of jointly being in model",...)
            tmp5 <- outer(tmp/sum(tmp),tmp)
            image(tmp3,tmp3,tmp2/tmp5,main="observed/expected ratio of jointly being in model",...)

            
          }
       }
    }
   if(x$choice == 1) {
      if(nmsx == 1)
         plot.logregmodel(x$model, pscript = pscript, title = title, nms = nms, ...)
      else plot.logregmodel(x$model, pscript = pscript, title = title, ...)
   }
   if(x$choice == 2 || x$choice == 6) {
      for(i in 1:x$nmodels) {
         xfit <- x$alltrees[[i]]
         if(nmsx == 1)
            plot.logregmodel(xfit, pscript = pscript, title = title, nms = nms, ...)
         else plot.logregmodel(xfit, pscript = pscript, title = title, ...)
      }
      data <- x$allscores
      rng.tr <- range(data[, 3])
      rng.sz <- range(data[, 2])
      rng <- range(data[, 1])
      if(pscript)
         postscript("scoreplot.tr.ps", print.it = FALSE, horizontal = TRUE, ...)
      plot(rng.sz, rng, type = "n", xlab = "model size", ylab = "scores")
      for(j in rng.tr[1]:rng.tr[2]) {
         which <- (1:dim(data)[1])[data[, 3] == j]
         lines(data[which, 2], data[which, 1], lty = 2)
         points(data[which, 2], data[which, 1], pch = 15, cex = 2.5)
         text(data[which, 2], data[which, 1], as.character(j), col = logregwhite, cex = 0.8)
      }
      if(pscript)
         dev.off()
   }
   if(x$choice == 3) {
      data <- x$cvscores
      kf <- data[1, 4]
      rng <- range(data[, 2])
      rng.tr <- range(data[seq(kf, dim(data)[1], kf), 6])
      rng.ts <- range(data[seq(kf, dim(data)[1], kf), 8])
      rng.ntr <- range(data[, 1])
      if(pscript)
         postscript("cvplot.tr.ps", print.it = FALSE, horizontal = TRUE, ...)
      plot(rng, rng.tr, type = "n", xlab = "model size", ylab = "training scores")
      for(j in rng.ntr[1]:rng.ntr[2]) {
         which <- (1:dim(data)[1])[data[, 3] == kf & data[, 1] == j]
         lines(data[which, 2], data[which, 6], lty = 2)
         points(data[which, 2], data[which, 6], pch = 15, cex = 2.5)
         text(data[which, 2], data[which, 6], as.character( j), col = logregwhite, cex = 0.8)
      }
      if(pscript)
         dev.off()
      if(pscript)
         postscript("cvplot.ts.ps", print.it = FALSE, horizontal = TRUE, ...)
      plot(rng, rng.ts, type = "n", xlab = "model size", ylab = "test scores")
      for(j in rng.ntr[1]:rng.ntr[2]) {
         which <- (1:dim(data)[1])[data[, 3] == kf & data[, 1] == j]
         lines(data[which, 2], data[which, 8], lty = 2)
         points(data[which, 2], data[which, 8], pch = 15, cex = 2.5)
         text(data[which, 2], data[which, 8], as.character(j), col = logregwhite, cex = 0.8)
      }
      if(pscript)
         dev.off()
      if(pscript)
         postscript("cvplot.ts2.ps", print.it = FALSE, horizontal = TRUE, ...)
      hmse <- 1
      data2 <- data[seq(kf, dim(data)[1], kf), c(1, 2, 8)]
      nn <- dim(data2)[1]
      vrb <- rep(0, nn)
      for(j in 1:nn)
         vrb[j] <- sqrt(var(data[(j - 1) * kf + (1:kf), 7])/ (kf - 1))
      data2 <- cbind(data2, vrb)
      rng.ts2 <- range(c(data2[, 3] + hmse * data2[, 4], data2[,
         3] - hmse * data2[, 4]))
      plot(rng, rng.ts2, type = "n", xlab = "model size", ylab = 
         "test scores +/- one standard error")
      for(j in rng.ntr[1]:rng.ntr[2]) {
         which <- (1:dim(data2)[1])[data2[, 1] == j]
         lines(data2[which, 2], data2[which, 3], lty = 2)
         for(k in which)
            lines(rep(data2[k, 2], 2), c(data2[k, 3] - hmse * data2[k, 4], data2[k, 3] + hmse * data2[ k, 4]))
         points(data2[which, 2], data2[which, 3], pch = 15, cex = 2.5)
         text(data2[which, 2], data2[which, 3], as.character(j), col = logregwhite, cex = 0.8)
      }
      if(pscript)
         dev.off()
   }
   if(x$choice == 4) {
      data2 <- c(x$nullscore, x$bestscore)
      data <- x$randscore
      rng <- range(c(data, data2))
      rng[1] <- rng[1] - sign(rng[1]) * 0.02 * rng[1]
      rng[2] <- rng[2] + sign(rng[2]) * 0.02 * rng[2]
      if(pscript)
         postscript("nmplot.ps", print.it = FALSE, horizontal = TRUE, ...)
      test <- hist(data, plot = FALSE)
      hist(data, xlim = rng, xlab = "score", ylab = "counts", main = "")
      lines(rep(data2[1], 2), 0.98 * c(0, max(test$counts)))
      lines(rep(data2[2], 2), 0.98 * c(0, max(test$counts)))
      text(data2[1], max(test$counts) * (1.02), "null model", cex = 0.7)
      text(data2[2], max(test$counts) * (1.02), "best model", cex = 0.7)
      if(title)
         title(main = "Null model randomization test")
      if(pscript)
         dev.off()
   }
   if(x$choice == 5) {
      modulo <- function(a, b)
      {
         fac <- floor(a/b)
         rem <- a - b * fac
         fac
      }
      hmint <- 25
      scb <- x$bestscore
      scb2 <- x$nullscore
      scb3 <- x$allscores
      data <- x$randscores
      data <- data[, data[1,  ] > 0]
      whm <- data[1,  ] + 1000 * data[2,  ]
      hm <- length(whm)
      rng <- range(c(unlist(data[ - (1:2),  ]), scb))
      dff <- (rng[2] - rng[1])
      sc <- floor(log10(dff))
      rng2 <- 10^( - sc) * rng
      rng2[1] <- floor(rng2[1])
      rng2[2] <- ceiling(rng2[2])
      dv <- max(modulo(hmint, diff(rng2)), round(diff(rng2)/hmint))
      brks <- (seq(rng2[1], rng2[2], 1/dv)) * (10^sc)
      if(pscript)
         postscript("plot.rand.ps", print.it = FALSE, horizontal = FALSE, width = 6, height = 10)
      op.logreg <- par(no.readonly = TRUE)
      par(mfrow = c(hm, 1), mar = c(1., 0, 0, 0), oma = c(3, 3, 2, 2))
      for(j in 1:hm) {
         if(j == hm)
            par(mar = c(1.9, 0, 0, 0))
         dd <- data[ - c(1:2), j]
         test <- hist(dd, breaks = brks, plot = FALSE)
         test <- c(test$density, test$counts)
         hist(dd, probability = TRUE, xlab = "", xaxt = "n",
            yaxt = "n", main = "", breaks = brks)
         abline(v = scb)
         abline(v = scb2)
         if(length(scb3) == 3) abline(v = scb3[3], lty=2)
         else{
            ww <- scb3[scb3[,3]==data[1,j] & scb3[,2]==data[2,j],1]
            if(length(ww)==1)
               abline(v=ww, lty=2)
         }
         pp <- paste("\ntrees:", data[1, j], "leaves:", data[2, j])
         text(brks[1], 0.95 * max(test), pp, adj = 0)
         if(j == 1)
            text(scb, 0.95 * max(test), " <== best model for unrandomized data", adj = 0)
         if(j == 1)
            text(scb2, 0.95 * max(test), " <== null model", adj = 0)
         box()
      }
      axis(1, brks)
      par(op.logreg)
      if(pscript)
         dev.off()
   }
   invisible()
}
plot.logregmodel <- function(x, pscript = FALSE, title = TRUE, nms, ...)
{
   if(class(x) != "logregmodel")
      stop("x not of class logregmodel")
   logregtmp <- hist(1:10, plot = FALSE)
   if(length(logregtmp$mids) > 0)
      logregwhite <- "white"
   if(length(logregtmp$mids) == 0)
      logregwhite <- 0
   ntr <- x$ntrees[1]
   msz <- x$nleaves[1]
   if(msz == -1)
      msz <- 0
   scores <- x$coef
   lscores <- length(scores)
   scores <- scores[(lscores - ntr + 1):lscores]
   for(j in 1:ntr) {
      if(pscript)
         postscript(paste("tree", ntr, ".", msz, ".", j, ".ps", sep = ""), print.it = FALSE, horizontal = TRUE, ...)
      if(missing(nms))
         plot.logregtree(x$trees[[j]], indents = c(0.2, 0, 0, 0), coef = scores[j], ...)
      else plot.logregtree(x$trees[[j]], indents = c(0.2, 0, 0, 0), coef = scores[j], nms = nms, ...)
      if(title && msz > 0)
         title(main = paste("tree", j, "out of", ntr, "total size is", msz))
      if(title && msz <= 0)
         title(main = paste("tree", j, "out of", ntr))
      if(pscript)
         dev.off()
   }
}

plot.logregtree <- function(x, nms, full = TRUE, and.or.cx = 1.0, leaf.sz = 1.0, 
   leaf.txt.cx = 1.0, coef.cx = 1.0, indents = rep(0, 4), coef = TRUE, coef.rd = 4, ...)
{
   ltree <- x
   if(class(ltree) != "logregtree")
      stop("ltree is not of class logregtree")
   ntchx1 <- indents[2]
   ntchy1 <- indents[1]
   ntchx2 <- indents[4]
   ntchy2 <- indents[3]
   logregtmp <- hist(1:10, plot = FALSE)
   if(length(logregtmp$mids) > 0)
      logregwhite <- "white"
   if(length(logregtmp$mids) == 0)
      logregwhite <- 0
   data <- ltree$trees
   cf <- ltree$coef
   if(is.na(cf))
      coef <- FALSE
   names(data) <- c("number", "conc", "knot", "neg", "pick")
   if(data$pick[1] == 0) {
      cat("empty tree", "\n")
      plot(c(0, 1), c(0, 1), type = "n", axes = FALSE, xlab = "", ylab = "")
      text(0.5, 0.5, "Empty tree", cex = 1.5*coef.cx)
   }
   else {
      if(!missing(nms)) {
         which <- (1:dim(data)[1])[data$conc == 3]
         data$knot[which] <- nms[data$knot[which]]
      }
      level <- ceiling(log(1 + max(data$number[data$pick >= 1]))/ log(2))
      max.p <- 2^level - 1
      data <- data[1:max.p,  ]
      data$x <- rep(0, max.p)
      data$y <- rep(0, max.p)
      left.b <- 2^(0:(level - 1))
      right.b <- 2^(1:level) - 1
      l.max <- max(data$number[left.b][data$pick[left.b] >= 1])
      r.max <- max(data$number[right.b][data$pick[right.b] >= 1])
      if(full) {
         l.shift <- 1/(l.max + 1)
         r.shift <- 1/(r.max + 1)
         w.shift <- min(l.shift, r.shift)
         plot(c(0.85 * (l.shift - ntchx1), 0.95 * (1 - r.shift + ntchx2)), c(-1 + ntchy1, 1.05 * ( - level -
            ntchy2)), type = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n", col = 0)
         if(coef) {
            text(0.85 * (l.shift - ntchx1), -1 + ntchy1, paste("Parameter =", round(cf, coef.rd)), cex = 1.5*coef.cx, adj = 0)
         }
      }
      else {
         plot(c(0, 1), c(-1,  - level), type = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n", col = 0)
         text(0.02, -1.2, paste("Parameter =", round(cf, coef.rd)), cex = 1.5*coef.cx, adj = 0)
      }
      for(k in level:1) {
         y.val <- ( - k)
         min.val <- 2^( - y.val - 1)
         max.val <- 2^( - y.val) - 1
         diff.x <- 1/(min.val + 1)
         x.val <- seq(diff.x, 1 - diff.x, diff.x)
         if(level == 1)
            x.val <- x.val * (0.95 * (1 - r.shift + ntchx2) - 0.85 * (l.shift - ntchx1)) + 0.85 * (l.shift - ntchx1)
         for(j in 1:length(x.val)) {
            which <- 2^(k - 1) + j - 1
            if(data$pick[which] == 1) {
               data$x[which] <- x.val[j]
               data$y[which] <- y.val
               if(data$conc[which] == 1) {
                  text(x.val[j], y.val, "and", cex = 2.0*and.or.cx)
               }
               if(data$conc[which] == 2) {
                  text(x.val[j], y.val, "or", cex = 2.0*and.or.cx)
               }
               if(data$conc[which] == 3) {
                  if(data$neg[which] == 0) {
                     points(x.val[j], y.val, pch = 0, cex = 7.0*leaf.sz)
                     text(x.val[j], y.val, as.character( data$knot[ which]), cex = 1.0*leaf.txt.cx)
                  }
                  if(data$neg[which] == 1) {
                     points(x.val[j], y.val, pch = 15, cex = 7.0*leaf.sz)
                     text(x.val[j], y.val, as.character( data$knot[ which]), col = logregwhite, cex = 1.0*leaf.txt.cx)
                  }
               }
               if(data$conc[which] == 4) {
                  if(data$neg[which] == 0) {
                     points(x.val[j], y.val, pch = 1, cex = 7.0*leaf.sz)
                     text(x.val[j], y.val, as.character( data$knot[ which]), cex = 1)
                  }
                  if(data$neg[which] == 1) {
                     points(x.val[j], y.val, pch = 16, cex = 7.0*leaf.sz)
                     text(x.val[j], y.val, as.character( data$knot[ which]), cex = 1, col = logregwhite)
                  }
               }
            }
         }
      }
      for(k in 1:max.p) {
         if((ceiling(log(k)/log(2)) < level) & (2 * k < max.p)) 
            {
            if((data$pick[k] == 1) & (data$pick[2 * k] == 1)) {
               lines(data$x[c(k, 2 * k)], c(data$ y[k] - 0.2, data$y[2 * k] + 0.2))
               lines(data$x[c(k, 2 * k + 1)], c(data$ y[k] - 0.2, data$y[2 * k + 1] + 0.2))
            }
         }
      }
      invisible()
   }
}


predict.logreg <- function(object, msz, ntr, newbin, newsep, ...) 
{
    if(class(object) != "logreg") 
        stop("object not of class logreg")
    if(object$choice > 2 && object$choice !=6) 
        stop("object$choice needs to be 1, 2, or 6")
    if(!missing(newbin)) {
        if(missing(msz) && missing(ntr) && missing(newsep)) 
            y <- frame.logreg(fit = object, newbin = newbin)
        if(missing(msz) && missing(ntr) && (missing(newsep) == FALSE)) 
            y <- frame.logreg(fit = object, newbin = newbin, newsep = newsep)
        if(missing(msz) && (missing(ntr) == FALSE) && missing(newsep)) 
            y <- frame.logreg(fit = object, newbin = newbin, ntr = ntr)
        if(missing(msz) && (missing(ntr) == FALSE) && (missing(newsep) == FALSE)) 
            y <- frame.logreg(fit = object, newbin = newbin, newsep = newsep, ntr = ntr)
        if((missing(msz) == FALSE) && missing(ntr) && missing(newsep)) 
            y <- frame.logreg(fit = object, newbin = newbin, msz = msz)
        if((missing(msz) == FALSE) && missing(ntr) && (missing(newsep) == FALSE)) 
            y <- frame.logreg(fit = object, newbin = newbin, newsep = newsep, msz = msz)
        if((missing(msz) == FALSE) && (missing(ntr) == FALSE) && missing(newsep)) 
            y <- frame.logreg(fit = object, newbin = newbin, ntr = ntr, msz = msz)
        if((missing(msz) == FALSE) && (missing(ntr) == FALSE) && (missing(newsep) == FALSE)) 
            y <- frame.logreg(fit = object, newbin = newbin, newsep = newsep, ntr = ntr, msz = msz)
    }
    else {
        if(missing(msz) && missing(ntr)) 
            y <- frame.logreg(fit = object)[, -1]
        if(missing(msz) && (missing(ntr) == FALSE)) 
            y <- frame.logreg(fit = object, ntr = ntr)[, -1]
        if((missing(msz) == FALSE) && missing(ntr)) 
            y <- frame.logreg(fit = object, msz = msz)[, -1]
        if((missing(msz) == FALSE) && (missing(ntr) == FALSE)) 
            y <- frame.logreg(fit = object, ntr = ntr, msz = msz)[, -1]
    }
    iik <- 0
    if(missing(msz) && missing(ntr) && object$select == "greedy") iik <- 1
    unstrip <- function(x) 
    {
        dd <- dim(x)
        y <- x
        if(length(dd) == 2) {
            dd2 <- dd[2]
            if(dd2 == 1) 
                y <- c(x[, 1])
            if(dd2 == 2) 
                y <- cbind(c(x[, 1]), c(x[, 2]))
            if(dd2 > 2) 
                y <- cbind(c(x[, 1]), c(x[, 2]), c(x[, 3]))
            if(dd2 > 3) 
                for(i in 4:dd2)
                   y <- cbind(y, c(x[, i]))
            y
        }
        if(length(dd) == 1 || length(dd) == 0) {
            y <- c(unlist(c(unlist(x))))
            names(y) <- NULL
        }
        y
    }
    ly <- length(y[, 1])
    if(is.matrix(y) == FALSE) 
        y <- matrix(unstrip(y), nrow = ly)
    lly <- length(y[,1])
    y <- y[, -1]
    if(lly==1) y <- matrix(y, nrow=1)
    if(object$type == "proportional.hazards") 
        y <- y[, -1]
    z <- NULL
    if(length(y) == ly) 
        y <- matrix(y, ncol = 1)
    y <- cbind(y, rep(1, ly))
    if(object$select == "single.model") {
        z <- rep(object$model$coef[1], length(y[, 1]))
        for(i in 2:length(object$model$coef)) 
            z <- z + object$model$coef[i] * y[, i - 1]
    }
    if(object$select == "multiple.models" |object$select == "greedy") {
        if(missing(msz)) 
            msz <- min(object$nleaves):max(object$nleaves)
        if(missing(ntr)) 
            ntr <- min(object$ntrees):max(object$ntrees)
        z <- NULL
        if(object$nsep > 0) {
            y1 <- y[, 1:object$nsep]
            y <- y[, - (1:object$nsep)]
        }
        jk <- 0
        for(i in 1:object$nmodels) {
            nt1 <- object$alltrees[[i]]$ntrees
            ms1 <- object$alltrees[[i]]$nleaves
            if(sum(nt1 == ntr) * sum(ms1 == msz) > 0| iik == 1) {
                jk <- jk + 1
                if(ly==1)
                  y <- matrix(y, nrow = 1)
                else{
                   if(length(y) == ly) 
                     y <- matrix(y, ncol = 1)
                }
                y3 <- y[, (1:nt1)]
                if(i != object$nmodels) 
                  y <- y[, - (1:nt1)]
                if(nt1 == 1) 
                  y3 <- matrix(y3, ncol = 1)
                if(object$nsep > 0) 
                  y3 <- cbind(y1, y3)
                cc <- object$alltrees[[i]]$coef
                z2 <- cc[1]
                if(length(cc) > 1) 
                  for(ii in 2:length(cc))
                     z2 <- z2 + cc[ii] * y3[, ii - 1]
                str <- paste("tr", nt1, ".lf", ms1, sep = "")
                if(length(z) > 0) {
                  z <- data.frame(z, z2)
                }
                else {
                  z <- data.frame(z2)
                }
                names(z)[jk] <- str
            }
        }
    }
    if(object$type == "classification") 
        z[z != 0] <- 1
    if(object$type == "logistic") 
        z <- exp(z)/(1 + exp(z))
    z
}
print.logreg <- function(x, nms, notnms, pstyle = 1, ...)
{
   if(class(x) != "logreg")
      stop("x not of class logreg")
   inn <- 2
   if(missing(notnms))
      inn <- 1
   if(missing(nms))
      inn <- 0
   if(inn == 0 & length(x$binnames) > 0) {
      nms <- x$binnames
      inn <- 1
   }
   if(x$choice == 1) {
      cat("score", round(x$model$score, 3),"\n")
      if(inn == 2)
         print.logregmodel(x$model, nms, notnms, pstyle = pstyle)
      if(inn == 1)
         print.logregmodel(x$model, nms, pstyle = pstyle)
      if(inn == 0)
         print.logregmodel(x$model, pstyle = pstyle)
   }
   if(x$choice == 4) {
       bb <- x$bestscore
       nn <- x$nullscore
       rr <- x$randscore
       cat("Null Score",round(nn,3),"; best score",round(bb,3),"\n")
       cat("Summary",length(rr),"Randomized scores\n")
       print(summary(rr))
       ii <- sum(rr<bb)
       cat(ii," randomized scores (",round(100*ii/length(rr),2),"%) are better than the best score\n")
   }
   if(x$choice == 5){
     uu <- x$allscores[,c(3,2,1)]
     uu[uu[,2]==0,1] <- 0
     uu <- cbind(uu[,1:2],x$nullscore,uu[,3],x$bestscore,0,0,0,0,0,0,0)
     for(i in 1:length(uu[,1])){
        vv <- x$randscores[-(1:2),i]
        uu[i,6:11] <- summary(vv)
        uu[i,12] <- sum(vv<x$bestscore)*100/length(vv)
     }
     uu <- data.frame(uu)
     names(uu) <- c("trees","leaves","null","start","best","rand: min","1st Qu","median","mean","3rd Qu","max","% < best")
     uu[uu[,2]==0,1] <- 0
     cat(length(vv),"randomizations\n")
     print(uu)
   }
   if(x$choice == 3){
      uu <- x$cvscores
      uu1 <- uu[,1]*1000+uu[,2]
      uu2 <- unique(uu1)
      vv <- uu[uu[,3]==uu[,4],c(1,2,6,5,8,7)]
      for(i in 1:length(uu2)){
         vv[i,4] <- sqrt(var(uu[uu1==uu2[i],5]))
         vv[i,6] <- sqrt(var(uu[uu1==uu2[i],7]))
      }
      names(vv)[4:6] <- c("train.sd","cv/test","cv/test.sd")
      print(vv)
   }
   if(x$choice==7)
      stop("works not when x$choice 7")
   if(x$choice == 2 || x$choice==6) {
      j <- length(x$allscores[, 1])
      for(i in 1:j) {
         cat(x$allscores[i, 3], "trees with", x$allscores[i, 2], "leaves: score is", round(x$allscores[i, 1], 3), "\n")
         if(inn == 2)
            print.logregmodel(x$alltrees[[i]], nms, notnms, pstyle = pstyle)
         if(inn == 1)
            print.logregmodel(x$alltrees[[i]], nms, pstyle = pstyle)
         if(inn == 0)
            print.logregmodel(x$alltrees[[i]], pstyle = pstyle)
      }
   }
   invisible()
}
print.logregmodel <- function(x, nms, notnms, pstyle = 1, ...)
{
   if(class(x) != "logregmodel")
      stop("x not of class logregmodel")
   inn <- 2
   if(missing(notnms))
      inn <- 1
   if(missing(nms))
      inn <- 0
   ntr <- x$ntrees[1]
   scores <- x$coef
   lscores <- length(scores)
   if(scores[1]!=0)cat(signif(scores[1],3))
    if(lscores>(ntr+1)){
       sscores <- scores[2:(lscores-ntr)]
       xl <- length(sscores)
       for(i in 1:xl){
          if(sscores[i]>0)
              cat(" +",signif(sscores[i], 3), " * sep[,",i,"]",sep="")
          else
              cat(" ",signif(sscores[i], 3), " * sep[,",i,"]",sep="")
      }
    }
    tscores <- scores[(lscores - ntr + 1):lscores]
    for(j in 1:ntr) {
        if(x$trees[[j]]$trees[1, 5] != 0) {
            tmp <- (tscores[j] >= 0)
            if (is.na(tmp) == FALSE) {
                if (tscores[j] >= 0)
                  cat(" +", signif(tscores[j], 3), " * ", sep = "")
                else cat(" ", signif(tscores[j], 3), " * ", sep = "")
                if (inn == 2)
                  print.logregtree(x$trees[[j]], nms, notnms,
                    pstyle = pstyle)
                if (inn == 1)
                  print.logregtree(x$trees[[j]], nms, pstyle = pstyle)
                if (inn == 0)
                  print.logregtree(x$trees[[j]], pstyle = pstyle)
            }
        }
        else {
            if (j == 1) {
                if (tscores[j] >= 0)
                  cat(" +", signif(tscores[j], 3), " * 1", sep = "")
                else cat(" ", signif(tscores[j], 3), " *  1",
                  sep = "")
            }
        }
   }
   cat("\n")
   invisible()
}
print.logregtree <- function(x, nms, notnms, pstyle = 1, ...)
{
   if(class(x) == "logregtree")
      x <- x$trees
   tree <- x
   tree <- matrix(c(unlist(tree)), ncol = 5)
   if(missing(nms)) {
      i <- max(tree[, 3]) + 1
      nms <- rep("X", i)
      for(j in 1:i)
         nms[j] <- paste("X", j, sep = "")
   }
   if(missing(notnms)) {
      notnms <- nms
      for(i in 1:length(notnms)) {
         if(pstyle != 1)
            notnms[i] <- paste("not", notnms[i])
         else notnms[i] <- paste("(not ", notnms[i], ")", sep = "")
      }
   }
   l1 <- length(tree[, 1])
   if(tree[1, 2] == 1 || tree[1, 2] == 2) {
      l3 <- rep(l1, 0)
      l3[1:3] <- 1:3
      if(l1 > 3) {
         l2 <- floor(l1/2)
         for(i in 2:l2) {
            j <- c(2 * i, 2 * i + 1)
            j <- j[j <= l1]
            l3[j] <- l3[i]
         }
      }
      oper <- "and"
      if(tree[1, 2] == 2)
         oper <- "or"
      cat("(")
      print.logregtree(tree[l3 == 2,  ], nms, notnms, pstyle)
      if(pstyle == 1) 
         cat(" ",oper," ",sep="")
      else 
         cat(" ) ",oper," ( ",sep="")
      print.logregtree(tree[l3 == 3,  ], nms, notnms, pstyle)
      cat(")")
   }
   else {
      if(tree[1, 2] == 0)
         cat("(0) ")
      else {
         if(tree[1, 4] == 1)
            cat(notnms[tree[1, 3]])
         else cat(nms[tree[1, 3]])
      }
   }
        invisible()
}
logreg <- function(resp, bin, sep, wgt, cens, type, select, ntrees, nleaves, penalty,
   seed, kfold, nrep, oldfit, anneal.control, tree.control, mc.control)
{ 
   call <- match.call()
   logreg.binary <- function(x)
   {
      l1 <- length(x)
      l2 <- sum(x == 0)
      l3 <- sum(x == 1)
      (l1 == (l2 + l3))
   }
   logreg.storetree <- function(x)
   {
      i1 <- matrix(x[ - (1:3)], ncol = 4, byrow = TRUE)
      i2 <- length(i1[, 1])
      i3 <- data.frame(1:i2, i1)
      names(i3) <- c("number", "conc", "knot", "neg", "pick")
      i3
   }
   n1 <- 0
   miss.select <- 1 * (missing(select))
   miss.sep <- 1 * (missing(sep))
   miss.wgt <- 1 * (missing(wgt))
   miss.cens <- 1 * (missing(cens))
   miss.type <- 1 * (missing(type))
   miss.ntrees <- 1 * (missing(ntrees))
   miss.nleaves <- 1 * (missing(nleaves))
   miss.penalty <- 1 * (missing(penalty))
   miss.seed <- 1 * (missing(seed))
   miss.kfold <- 1 * (missing(kfold))
   miss.nrep <- 1 * (missing(nrep))
   miss.tree <- 1 * (missing(tree.control))
   miss.anneal <- 1 * (missing(anneal.control))
   miss.mc <- 1 * (missing(mc.control))
   miss.binnames <- 1
   if(!missing(oldfit)) {
      if(class(oldfit) != "logreg")
         stop("oldfit is not of class logreg")
      if(missing(resp) && length(oldfit$response) > 0)
         resp <- oldfit$response
      if(length(oldfit$binnames) > 0 && missing(bin)) {
         binnames <- oldfit$binnames
         miss.binnames <- 0
      }
      if(missing(bin) && length(oldfit$binary) > 0)
         bin <- oldfit$binary
      if(missing(sep) && oldfit$nsep > 0) {
         miss.sep <- 0
         sep <- oldfit$separate
      }
      if(missing(cens) && length(oldfit$censor) > 0) {
         miss.cens <- 0
         cens <- oldfit$censor
      }
      if(missing(wgt) && length(oldfit$weight) > 0) {
         miss.wgt <- 0
         wgt <- oldfit$weight
      }
      if(missing(nrep) && length(oldfit$nrep) > 0) {
         miss.nrep <- 0
         nrep <- oldfit$nrep
      }
      if(missing(kfold) && length(oldfit$kfold) > 0) {
         miss.kfold <- 0
         kfold <- oldfit$kfold
      }
      if(missing(penalty) && length(oldfit$penalty) > 0) {
         miss.penalty <- 0
         penalty <- oldfit$penalty
      }
      if(missing(anneal.control) && length(oldfit$anneal.control) > 0) {
         miss.anneal <- 0
         anneal.control <- oldfit$anneal.control
      }
      if(missing(mc.control) && length(oldfit$mc.control) > 0) {
         miss.mc <- 0
         mc.control <- oldfit$mc.control
      }
      if(missing(tree.control) && length(oldfit$tree.control) > 0) {
         miss.tree <- 0
         tree.control <- oldfit$tree.control
      }
      if(missing(ntrees) && length(oldfit$ntrees) > 0) {
         miss.ntrees <- 0
         ntrees <- oldfit$ntrees
         if(length(ntrees) == 1)
            ntrees <- c(ntrees, ntrees)
      }
      if(missing(nleaves) && length(oldfit$nleaves) > 0) {
         miss.nleaves <- 0
         nleaves <- oldfit$nleaves
         if(length(nleaves) == 1)
            nleaves <- c(nleaves, nleaves)
      }
      if(missing(seed) && length(oldfit$seed) > 0) {
         miss.seed <- 0
         seed <- oldfit$seed
      }
      if(missing(type) && length(oldfit$type) > 0) {
         miss.type <- 0
         if(oldfit$type == "own.scoring")
            type <- 0
         if(oldfit$type == "classification")
            type <- 1
         if(oldfit$type == "regression")
            type <- 2
         if(oldfit$type == "logistic")
            type <- 3
         if(oldfit$type == "proportional.hazards")
            type <- 4
         if(oldfit$type == "exponential.survival")
            type <- 5
      }
      if(missing(select) && length(oldfit$select) > 0) {
         miss.select <- 0
         if(oldfit$select == "single.model")
            select <- 1
         if(oldfit$select == "multiple.models")
            select <- 2
         if(oldfit$select == "cross.validation")
            select <- 3
         if(oldfit$select == "null.model.test")
            select <- 4
         if(oldfit$select == "randomization.test")
            select <- 5
         if(oldfit$select == "bayesian")
            select <- 7
         if(oldfit$select == "greedy")
            select <- 6
      }
   }
   if(miss.type == 1)
      type <- 2
   if(miss.penalty == 1)
      penalty <- 0
   if(miss.kfold == 1)
      kfold <- 10
   if(miss.nrep == 1)
      nrep <- 25
   if(miss.seed == 1)
      seed <- 0
   if(miss.binnames == 1) {
      if(is.data.frame(bin))
         binnames <- names(bin)
      else binnames <- NULL
   }
   if(is.numeric(type))
      mdl <- type
   else {
      if(type == "c" | type == "C")
         mdl <- 1
      if(type == "r" | type == "R")
         mdl <- 2
      if(type == "l" | type == "L")
         mdl <- 3
      if(type == "s" | type == "S")
         mdl <- 4
      if(type == "e" | type == "E")
         mdl <- 5
      if(type == "o" | type == "O")
         mdl <- 0
   }
   if(mdl < 0 || mdl > 5)
      stop("not an implemented model")
   choice <- 99
   if(is.numeric(select))
      choice <- select
   else {
      if(select == "s" | select == "S")
         choice <- 1
      if(select == "m" | select == "M")
         choice <- 2
      if(select == "c" | select == "C")
         choice <- 3
      if(select == "n" | select == "N")
         choice <- 4
      if(select == "r" | select == "R")
         choice <- 5
      if(select == "b" | select == "B")
         choice <- 7
      if(select == "g" | select == "G")
         choice <- 6
   }
   if(choice < 0 || choice > 7)
      stop("not a valid selection approach")
   n1 <- length(resp)
   if((mdl == 1 || mdl == 3) && (logreg.binary(resp) == FALSE))
      stop("some non binary response data, and a binary fitting function")
   if((mdl == 4 || mdl == 5) && (min(resp) < 0))
      stop("survival data needs to be positive")
   bin <- as.matrix(bin)
   n2 <- length(bin)/n1
   if(n2 != floor(n2))
      stop("length of binary data not a multiple of response")
   if(n2 < 1)
      stop("binary predictors mandatory")
   if(logreg.binary(bin) == FALSE)
      stop("some non binary data among binary predictors")
   nsep <- 0
   if(miss.sep == 0) {
      sep <- as.matrix(sep)
      nsep <- length(sep)/n1
      if(nsep != floor(nsep))
         stop("length of separate predictors not a multiple of response")
   }
   else {
      sep <- 0
   }
  if(nsep > 0 && (mdl == 2 || mdl == 5 || mdl == 4)){
      tep <- sep
      k <- 0
      for(i in 1:nsep){
         if((min(tep[,i])!=0 || max(tep[,i])!=1) && k!=2){
            k <- 1
            if(max(tep[,i])!=min(tep[,i])) {
               tep[,i] <- (tep[,i]-min(tep[,i]))/(max(tep[,i])-min(tep[,i]))
               if(logreg.binary(tep[,i])==FALSE) k <- 2
            }
            else
               tep[,i] <- 0
         }
      }
      if(k==1 && mdl!=4){
         sep <- tep
         print("separate covariates that were binary recoded to 0/1")
      }
      if(k==2 && mdl!=4){
         if(mdl==5)
            stop("Exponential survival models are only implemented if\nall separate predictors are binary")
         if(mdl==3)
         cat("Logistic regression runs much fast if\nall separate predictors are binary\n")
      }
      if(mdl==4 && k!=2){
         cat(paste("The code will be run much faster if you use exponential regression\n",
                                    "instead of proportional hazards models after a cumulative hazards\n",
                                    "transformation, as all your separate predictors are binary\n",
                                    "type help(cumhaz) for more information\n"))
      }
   }
   if(nsep==0 && mdl==4){
      cat(paste("The code will be run much faster if you use exponential regression\n",
                            "instead of proportional hazards models after a cumulative hazards\n",
                            "transformation, as you have no separate continuous predictors\n",
                            "type help(cumhaz) for more information\n"))
   }
   if(miss.tree == 1)
      tree.control <- logreg.tree.control(n1 = n1)
   else tree.control <- logreg.tree.control(tree.control, n1 = n1)
   if(miss.wgt == 1)
      wgt <- rep(1, n1)
   else {
      if(length(wgt) != n1)
         stop("weight does not have the right length")
      if(min(wgt) < 0 && mdl != 0)
         stop("weight has to be larger than 0")
   }
   if(miss.cens == 1)
      cens <- rep(1, n1)
   else {
      if(length(cens) != n1)
         stop("cens doesn't have right length")
      if((mdl == 4 || mdl == 5) && (logreg.binary(cens) == FALSE))
         stop("censoring data has to be 0/1")
   }
   if(choice == 1 || choice == 4 || choice == 7) {
      if(miss.ntrees)
         ntr <- c(1, 1)
      else ntr <- rep(ntrees[1], 2)
      if(miss.nleaves)
         msz <- c(-1, -1)
      else msz <- rep(nleaves[1], 2)
   }
   else {
      if(miss.ntrees)
         ntr <- c(1, 1)
      else ntr <- ntrees[1:2]
      if(is.na(ntr[2]))
         ntr[2] <- ntr[1]
            if(choice == 6 && miss.nleaves == TRUE)
                nleaves <- rep(tree.control$treesize*ntr[2],2)
      msz <- nleaves[1:2]
      if(is.na(msz[2]))
         msz[2] <- msz[1]
      if(ntr[1] > ntr[2])
         stop(" upper limit of number of trees smaller than lower limit")
      if(msz[1] > msz[2])
         stop(" upper limit of number of leaves smaller than lower limit")
      if(msz[1] < 0)
         stop("number of leaves needs to be at least 0")
   }
   if(ntr[1] < 1)
      stop("number of trees needs to be at least 1")
   if(mdl == 1 && max(ntr) > 1)
      stop("for classification only 1 tree is possible")
   if(mdl == 1 && nsep > 0)
      stop("for classification no separate predictors are possible")
   if(miss.anneal == 1)
      anneal.control <- logreg.anneal.control()
   else anneal.control <- logreg.anneal.control(anneal.control)
   if(miss.mc == 1)
      mc.control <- logreg.mc.control()
   else mc.control <- logreg.mc.control(mc.control)
   if(miss.mc != 0 && select == 7)
      anneal.control$update <- mc.control$update
   if(penalty < 0 && (choice == 1 || choice == 6))
      stop("penalty should be at least 0")
   if(choice == 3) {
      nrep <- kfold
      if(nrep < 2)
         stop("kfold needs to be at least 2")
   }
   if(choice > 3 && nrep < 1 && choice < 6)
      stop("nrep needs to be at least 1")
   xseed <- seed
   if(xseed == 0)
      xseed <- floor(runif(1) * 1000000) + 1
   ipars <- c(mdl, msz[1:2], tree.control$treesize, ntr[1:2], tree.control$
      opers, anneal.control$update, xseed, nrep, choice, 
      anneal.control$earlyout, tree.control$minmass, mc.control$
      nburn, mc.control$niter, mc.control$output)
   rpars <- c(anneal.control$start, anneal.control$end, anneal.control$
      iter, penalty, mc.control$hyperpars)
   orders <- order(rank(resp) + runif(n1)/1000000)
   nkn <- tree.control$treesize * 2 - 1
   nxx <- 2
   
   if(choice == 1 || choice == 7) {
      na <- ntr[1] * (nkn * 4 + 3)
      nb <- nsep + ntr[1] + 1
      nc <- 2
   }
   if(choice == 2) {
      nd <- (ntr[2] - ntr[1] + 1) * (msz[2] - msz[1] + 1)
      na <- nd * ntr[2] * (nkn * 4 + 3)
      nb <- nd * (nsep + ntr[2] + 1)
      nc <- nd
   }
   if(choice == 6) {
      nd <- msz[2] + 2
      na <- nd * ntr[2] * (nkn * 4 + 3)
      nb <- nd * (nsep + ntr[2] + 1)
      nc <- nd
   }
   if(choice == 3) {
      na <- nb <- 2
      nd <- (ntr[2] - ntr[1] + 1) * (msz[2] - msz[1] + 1)
      nc <- nd * 8 * nrep
   }
   if(choice == 4) {
      na <- nb <- 2
      nc <- nrep + 2
   }
   n100 <- -100
   if(choice == 7) {
      na <- 256
      nb <- n2
      nc <- n2 * n2
      if(abs(mc.control$output) < 2)
         nc <- 1
      nxx <- nc * n2
      if(abs(mc.control$output) < 3)
         nxx <- 1
      n100 <- 0
   }
   if(choice != 5)
      xtree <- rep(n100, na)
   if(choice == 5) {
      nb <- 2
      nd <- (ntr[2] - ntr[1] + 1) * (msz[2] - msz[1] + 1)
      nc <- (nrep + 2) * nd + 2
      xtree <- NULL
      if(tree.control$treesize != oldfit$tree.control$treesize)
         stop("treesize should match treesize in oldfit")
      if(oldfit$choice != 2)
         stop("oldfit not an object from a multiple model fit")
      if(oldfit$choice == 2) {
         for(i in 1:oldfit$nmodels) {
            tmp <- oldfit$alltrees[[i]]
            if(i == 1)
               bestscore <- tmp$score
            else {
               xscore <- tmp$score
               if(xscore < bestscore)
                  bestscore <- xscore
            }
            for(j in 1:tmp$ntrees[1]) {
               xtree <- c(xtree, tmp$nleaves[1], tmp$
                  ntrees[1], j)
               xtree <- c(xtree, c(t(as.matrix(tmp$
                  trees[[j]]$trees[, -1]))))
            }
         }
         allscores <- oldfit$allscores
      }
      xtree <- c(xtree, rep(n100, 4 * nkn * 10))
   }
   ip4 <- 2*ipars[4]+1
   fit <- .Fortran("slogreg",
      as.integer(n1),
      as.integer(n2),
      as.integer(nsep),
      ip = as.integer(ipars),
      as.single(rpars),
      as.single(t(sep)),
      as.integer(cens),
      as.integer(orders),
      as.single(resp),
      as.single(wgt),
      as.integer(t(bin)),
      trees = as.integer(xtree),
      coef = as.single(rep(n100, nb)),
      scores = as.single(rep(n100, nc)),
      as.integer(ipars[6]),
      as.integer(ip4),
      as.integer(rep(0,2*ipars[6]*ip4*n1)),
      as.integer(rep(0,7*ipars[6]*(ip4+1)*n2*4)),
      as.single(rep(0,7*ipars[6]*(ip4+1)*n2*4)),
      as.integer(t(bin)),
      rd4 = as.integer(rep(0, nxx)),
      PACKAGE="LogicReg")
   if(fit$ip[1]<(-900))stop("fatal declaration error - reduce problem or recompile package")
   ipars <- (choice)
   rpars <- c(penalty)
   if(mdl == 0)
      type <- "own.scoring"
   if(mdl == 1)
      type <- "classification"
   if(mdl == 2)
      type <- "regression"
   if(mdl == 3)
      type <- "logistic"
   if(mdl == 4)
      type <- "proportional.hazards"
   if(mdl == 5)
      type <- "exponential.survival"
   if(choice == 1)
      chs <- "single.model"
   if(choice == 2)
      chs <- "multiple.models"
   if(choice == 3)
      chs <- "cross.validation"
   if(choice == 4)
      chs <- "null.model.test"
   if(choice == 5)
      chs <- "randomization.test"
   if(choice == 7)
      chs <- "bayesian"
   if(choice == 6)
      chs <- "greedy"
   if(tree.control$opers != 2 && tree.control$opers != 3)
      tree.control$operators <- "both"
   if(tree.control$opers == 2)
      tree.control$operators <- "and"
   if(tree.control$opers == 3)
      tree.control$operators <- "or"
   m1 <- list(nsample = n1, nbinary = n2, nseparate = nsep, type = type,
      select = chs, anneal.control = anneal.control, tree.control = 
      tree.control, seed = seed, choice = choice)
   if(choice == 7) {
      m1$anneal.control <- NULL
      v1 <- fit$trees
      v3 <- 1:length(v1)
      v3 <- max(v3[v1 > 0])
      v1 <- v1[1:v3]
      v1 <- data.frame(0:(v3 - 1), v1)
      names(v1) <- c("size", "number")
      m1$size <- v1
      m1$single <- fit$coef
      m1$single[m1$single < 0] <- 0
      m1$mc.control <- mc.control
      if(abs(mc.control$output) > 1) {
         m1$double <- matrix(fit$scores, ncol = length(fit$coef))
         m1$double[m1$double < 0] <- 0
      }
      if(abs(mc.control$output) > 2) {
         m1$triple <- array(fit$rd4, dim = c(length(fit$coef), length(fit$coef), length(fit$coef)))
      }
   }
   if(choice == 1 || choice == 4) {
      m1$nleaves <- msz[1]
      m1$ntrees <- ntr[1]
   }
   else {
      m1$nleaves <- msz
      m1$ntrees <- ntr
   }
   if(choice == 4 || choice == 5)
      m1$nrep <- nrep
   if(choice == 3)
      m1$kfold <- kfold
   if(choice == 1)
      m1$penalty <- penalty
   m1$response <- resp
   m1$binary <- bin
   m1$separate <- sep
   m1$censor <- cens
   m1$weight <- wgt
   if(choice == 3) {
      cvscores <- matrix(fit$scores, ncol = 8, byrow = TRUE)
      cvscores <- cvscores[cvscores[, 1] > -99,  ]
      cvscores <- data.frame(cvscores)
      names(cvscores) <- c("ntree", "nleaf", "k", "kfold", "train", "train.ave", "test", "test.ave")
      m1$cvscores <- cvscores
   }
   if(choice == 4) {
      m1$nullscore <- fit$scores[1]
      m1$bestscore <- fit$scores[2]
      m1$randscore <- fit$scores[3:(2 + nrep)]
   }
   if(choice == 5) {
      m1$nullscore <- fit$scores[1]
      m1$bestscore <- fit$scores[2]
      m1$allscores <- allscores
      randscores <- matrix(fit$scores[-(1:2)], nrow = nrep + 2, byrow = FALSE)
      randscores <- randscores[, randscores[1,  ] > -99]
      randscores <- data.frame(randscores)
      row.names(randscores) <- c("ntree", "nleaf", as.character(1:nrep))
      m1$randscores <- randscores
   }
   if(choice == 1) {
      m2 <- list()
      m2$ntrees <- ntr
      m2$nleaves <- msz
      m2$score <- fit$scores[1]
      m2$coef <- fit$coef[1:(nsep + ntr[1] + 1)]
      if(mdl == 4)
         m2$coef[1] <- 0
      if(mdl == 1)
         m2$coef[1:2] <- c(0, 1)
      class(m2) <- "logregmodel"
      lx <- 3 + 4 * nkn
      m2$trees <- list()
      for(i in 1:ntr[1]) {
         m3 <- list()
         m3$whichtree <- i
         m3$coef <- fit$coef[1 + nsep + i]
         if(mdl == 1)
            m3$coef <- NA
         m3$trees <- logreg.storetree(fit$trees[(i - 1) * lx + (1:lx)])
         class(m3) <- "logregtree"
         m2$trees[[i]] <- m3
      }
      m1$model <- m2
   }
   if(choice == 6) {
      m1$nmodels <- fit$ip[1]
      m1$alltrees <- list()
      m1$allscores <- matrix(0, nrow = m1$nmodels, ncol = 3)
      j <- 0
      k <- 0
      for(i in 1:m1$nmodels) {
         m2 <- list()
         m2$score <- fit$scores[i]
         m2$nleaves <- fit$trees[j + 1]
         m2$ntrees <- fit$trees[j + 2]
         m1$allscores[i,  ] <- c(fit$scores[i], fit$trees[j + (1:2)])
         m2$coef <- fit$coef[k + (1:(m2$ntrees + nsep + 1))]
         if(mdl == 4)
            m2$coef[1] <- 0
         if(mdl == 1)
            m2$coef[1:2] <- c(0, 1)
         class(m2) <- "logregmodel"
         m2$trees <- list()
         lx <- 3 + 4 * nkn
         for(l in 1:m2$ntrees) {
            m3 <- list()
            m3$whichtree <- l
            m3$coef <- m2$coef[1 + nsep + l]
            if(mdl == 1)
               m3$coef <- NA
            m3$trees <- logreg.storetree(fit$trees[j + ( l - 1) * lx + (1:lx)])
            class(m3) <- "logregtree"
            m2$trees[[l]] <- m3
         }
         j <- j + lx * m2$ntrees
         k <- k + (m2$ntrees + nsep + 1)
         m1$alltrees[[i]] <- m2
      }
   }
   if(choice == 2) {
      m1$nmodels <- fit$ip[1]
      m1$alltrees <- list()
      m1$allscores <- matrix(0, nrow = m1$nmodels, ncol = 3)
      j <- 0
      k <- 0
      for(i in 1:m1$nmodels) {
         m2 <- list()
         m2$score <- fit$scores[i]
         m2$nleaves <- fit$trees[j + 1]
         m2$ntrees <- fit$trees[j + 2]
         m1$allscores[i,  ] <- c(fit$scores[i], fit$trees[j + (1:2)])
         m2$coef <- fit$coef[k + (1:(m2$ntrees + nsep + 1))]
         if(mdl == 4)
            m2$coef[1] <- 0
         if(mdl == 1)
            m2$coef[1:2] <- c(0, 1)
         class(m2) <- "logregmodel"
         m2$trees <- list()
         lx <- 3 + 4 * nkn
         for(l in 1:m2$ntrees) {
            m3 <- list()
            m3$whichtree <- l
            m3$coef <- m2$coef[1 + nsep + l]
            if(mdl == 1)
               m3$coef <- NA
            m3$trees <- logreg.storetree(fit$trees[j + (l - 1) * lx + (1:lx)])
            class(m3) <- "logregtree"
            m2$trees[[l]] <- m3
         }
         j <- j + lx * m2$ntrees
         k <- k + (m2$ntrees + nsep + 1)
         m1$alltrees[[i]] <- m2
      }
   }
   m1$call <- call
   m1$binnames <- binnames
   class(m1) <- "logreg"
   m1
}
logreg.mc.control <- function(nburn = 1000, niter = 25000, hyperpars = 0, update = 0, output = 4)
{
   if(is.list(nburn)) {
      if(length(nburn$niter) > 0)
         niter <- nburn$niter
      if(length(nburn$hyperpars) > 0)
         hyperpars <- nburn$hyperpars
      if(length(nburn$output) > 0)
         output <- nburn$output
      if(length(nburn$update) > 0)
         update <- nburn$update
      if(length(nburn$nburn) > 0)
         nburn <- nburn$nburn
      else nburn <- 1000
   }
   hyperpars <- c(hyperpars, rep(0, 4))[1:4]
   list(nburn = nburn, niter = niter, hyperpars = hyperpars, update = 
      update, output = output)
}
cumhaz <- function(y,d)
{
   logregtmp <- hist(1:10, plot = FALSE)
   #if(length(logregtmp$mids) > 0)
   #   library(survival)
   if(missing(d)) d <- rep(1,length(y))
   d - coxph(Surv(y,d)~1,iter.max=0)$residuals
}
