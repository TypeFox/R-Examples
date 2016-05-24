theoretical.depth.ratio <- function(CNt, CNn = 2, cellularity, ploidy, normal.ploidy = 2, avg.depth.ratio = 1) {
   cellu.copy.term   <- (1 - cellularity) + (CNt / CNn * cellularity)
   ploidy.cellu.term <- (ploidy / normal.ploidy * cellularity) + 1 - cellularity
   avg.depth.ratio * cellu.copy.term / ploidy.cellu.term
}

theoretical.mufreq <- function(Mt, CNt, CNn = 2, cellularity) {
   normal.alleles <- (CNt - Mt) * cellularity + CNn * (1 - cellularity)
   all.alleles    <- (CNt * cellularity) + CNn * (1 - cellularity)
   1 - (normal.alleles / all.alleles)
}

types.matrix <- function(CNt.min, CNt.max, CNn = 2) {
   cn.ratio.vect <- seq(from = CNt.min / CNn, to =  CNt.max / CNn, by = 1 / CNn)
   CNt           <- cn.ratio.vect * CNn
   mut.comb      <- lapply(CNt, FUN = function(x) seq(from = 0, to = x))
   times.muts    <- sapply(mut.comb, length)
   data.frame(CNn = CNn, CNt = rep(CNt, times = times.muts),
              Mt = unlist(mut.comb))
}

model.points <- function(cellularity, ploidy,
                         types, avg.depth.ratio) {
   mufreqs     <-  theoretical.mufreq(cellularity = cellularity , CNn = types[, 1], CNt = types[, 2], Mt = types[, 3])
   depth.ratio <-  theoretical.depth.ratio(cellularity = cellularity, ploidy = ploidy,
                                       CNt = types[, 2], CNn = types[, 1],
                                       avg.depth.ratio = avg.depth.ratio)
   cbind(mufreqs, depth.ratio)
}

# theoretical.baf <- function(cellularity = 0.5, CNt = 2, B = 1, CNn = 2){
#    B.tot <- ((B * cellularity)  + (1 - cellularity)) /
#             ((CNt * cellularity) + CNn*(1 - cellularity))
#    B.tot
# }

theoretical.baf <- function(CNt, CNn = 2, cellularity) {
   alleles       <- seq(from = 1, to = CNt, by = 1)
   max.b <- function(CNt) {
      max.b.alleles <- CNt / 2
      if (CNt %% 2 != 0 ) {
         max.b.alleles <- trunc(max.b.alleles)
      }
      max.b.alleles
   }
   fract.normal.alleles <- (1 - cellularity)
   res                  <- list()
   for (i in 1:length(alleles)) {
      max.b.alleles <- max.b(alleles[i])
      max.a.alleles <- alleles[i] - max.b.alleles
      decrements.b  <- seq(from = max.b.alleles, to = 0, by = -1)
      res[[i]]     <- list()
      for (n in 1:length(decrements.b)) {
         A.i <- (max.a.alleles + decrements.b[n])
         B.i <- (max.b.alleles - decrements.b[n])
         BAF <- (fract.normal.alleles + (cellularity * B.i)) / ((alleles[i] * cellularity) + (CNn * fract.normal.alleles))
         res[[i]][[n]] <- cbind(A = A.i, B = B.i, BAF = BAF, CNt = alleles[i])
      }
   }
   for (i in 1:length(res)) {
      res[[i]] <- do.call(rbind,res[[i]])
   }
   as.data.frame(do.call(rbind,res))
}

expected.baf <- function(sd, ...) {
   baf      <- theoretical.baf(...)
   baf.t2 <- function(BAF, sd, by = 0.001){
      bafs   <- seq(0,1,0.001)
      b.b    <- dt2(x = bafs, mean = BAF, sd = sd, df = 5)
      b.a    <- dt2(x = bafs, mean = 1-BAF, sd = sd, df = 5)
      half.b <- bafs[bafs <= 0.5]
      b <- (b.b+b.a)[bafs <= 0.5]
      weighted.mean(half.b,b)
   }
   BAF <- mapply(FUN = baf.t2, baf$BAF,
                 MoreArgs = list(sd = sd))
   wgh <- dt2(x = baf$BAF, mean = 0.5, sd = sd, df = 5)
   wgh <- wgh/max(wgh)
   mean.bf <- function(x) {
      weighted.mean(x=c(x["BAF"], x["eBAF"]), w = c((1-x["wgh"]), x["wgh"]))
   }
   baf$BAF <- apply(cbind(BAF = baf$BAF, eBAF = BAF, wgh = wgh), 1, FUN = mean.bf)
   baf
}

baf.model.points <- function (CNt.min, CNt.max, CNn = 2, cellularity,
                              ploidy, avg.depth.ratio) {
   mufreq.depth.ratio <- model.points(cellularity = cellularity, ploidy = ploidy,
                                      types = cbind(CNn = CNn, CNt = CNt.min:CNt.max, Mt = 0),
                                      avg.depth.ratio = avg.depth.ratio)
   model.d.ratio      <- cbind(CNt = CNt.min:CNt.max, depth.ratio = mufreq.depth.ratio[, 2])
   model.baf          <- theoretical.baf(CNn = CNn, CNt = CNt.max, cellularity = cellularity)
   if (CNt.min == 0) {
      model.baf <- as.data.frame(rbind(c(0, 0, 0.5, 0), model.baf))
   }
   model.pts          <- merge(model.baf, model.d.ratio)
   model.pts
}
