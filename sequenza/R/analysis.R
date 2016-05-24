read.seqz <- function (file, nrows = -1, fast = FALSE, gz = TRUE, header = TRUE,
    colClasses = c('character', 'integer', 'character', 'integer',
      'integer', 'numeric', 'numeric', 'numeric', 'character',
      'numeric', 'numeric', "character", "character", "character"),
    chr.name = NULL, n.lines = NULL, ...) {
   if (!is.null(n.lines) & is.null(chr.name)) fast <-  FALSE
   if(fast && nrows == -1) {
    if(gz) {
       if (!is.null(chr.name)) {
          wc <- system(paste('gzip -d -c ',file,' | grep -c "^', chr.name, '\t"', sep = ''), intern = TRUE)
       } else {
          wc <- system(paste('gzip -d -c', file, '| wc'), intern = TRUE)
       }
    } else {
       if (!is.null(chr.name)) {
          wc <- system(paste(paste('grep -c "^', chr.name, '\t"', sep = ''), file, sep = ' '), intern = TRUE)
       } else {
          wc <- system(paste('wc', file), intern = TRUE)
       }
    }
    if (is.null(chr.name)) {
       wc <- sub("^ +", "", wc)
       wc <- strsplit(wc, ' ')[[1]][1]
    }
    nrows <- max(as.integer(wc), 1)
    message('Reading ', nrows, ' lines...')
  }
   if (!is.null(chr.name)) {
      if (gz) {
         grep.part <- paste("gzip -d -c ", file," | grep '^", chr.name, "\t'", sep = "")
      } else {
         grep.part <- paste("grep '^", chr.name, "\t'", sep = "")
      }
      seqz.data   <- read.delim(pipe(grep.part), nrows = nrows, colClasses = colClasses, header = FALSE, ...)
      if (header == TRUE) {
         head       <- colnames(read.table(file, header = TRUE, nrows = 1 ))
         colnames(seqz.data) <- head
      }
      seqz.data
   } else {
      if (!is.null(n.lines)){
         if (!is.numeric(n.lines) | length(n.lines) != 2) stop("n.lines must be a vector of 2 integers")
         n.lines <- round(sort(n.lines), 0)
         if (header == TRUE) {
            n.lines <- n.lines + 1
         }
         if(gz) {
            seqz.data <- read.delim(pipe(paste("gzip -d -c", file,"| sed -n '", paste(n.lines[1], n.lines[2], sep = ","),"p'")),
                       colClasses = colClasses, nrows = 1 + n.lines[2] - n.lines[1], header = FALSE,...)
         }  else{
            seqz.data <- read.delim(pipe(paste("sed -n '", paste(n.lines[1], n.lines[2], sep = ","),"p'", file)),
                       colClasses = colClasses, nrows = 1 + n.lines[2] - n.lines[1], header = FALSE, ...)
         }
         if (header == TRUE) {
            head  <- colnames(read.table(file, header = TRUE, nrows = 1 ))
            colnames(seqz.data) <- head
         }
         seqz.data
      } else {
         read.delim(file, nrows = nrows, colClasses = colClasses, header = header, ...)
      }
   }
}

read.acgt <- function (file, colClasses = c('character', 'integer', 'character', 'integer',
                                            'integer', 'integer', 'integer', 'integer', "character"),
                       ...) {
   read.seqz(file = file , colClasses = colClasses, ...)
}

gc.norm <- function (x, gc) {
   dr.by.gc <- split(x = x, f = gc)
   raw <- t(sapply(dr.by.gc, quantile, probs = c(0.25, 0.5, 0.75), na.rm = TRUE))
   dr.by.gc.mean <- sapply(dr.by.gc, FUN = mean, na.rm = TRUE)
   dr.by.gc.median <- sapply(dr.by.gc, FUN = median, na.rm = TRUE)
   adj <- sweep(raw, 1, dr.by.gc.median, '/')
   list(raw = raw, adj = adj, gc.values = as.numeric(names(dr.by.gc)),
        raw.mean = dr.by.gc.mean, raw.median = dr.by.gc.median)
}

gc.sample.stats <- function (file, gz = TRUE) {
   colClasses = c('character', 'numeric', 'numeric')
   if (gz) {
      seqz.data <- read.delim(pipe(paste('gzip -d -c', file, '| cut -f 1,6,10')), colClasses = colClasses)
   } else {
      seqz.data <- read.delim(pipe(paste('cut -f 1,6,10', file)), colClasses = colClasses)
   }
   gc.stats <- gc.norm(x  = seqz.data$depth.ratio,
                       gc = seqz.data$GC.percent)
   chr.ord  <- unique(seqz.data$chromosome)
   chr.dim  <- lapply(X = split(seqz.data$chromosome, seqz.data$chromosome), FUN = length)
   chr.dim  <- data.frame(chr = chr.ord, n.lines = do.call(rbind,chr.dim[chr.ord]))
   chr.dim$start <- cumsum(c(1, chr.dim$n.lines[-length(chr.dim$n.lines)]))
   chr.dim$end   <- chr.dim$start + chr.dim$n.lines - 1
   gc.stats$file.metrics <- chr.dim
   gc.stats
}

windowValues <- function(x, positions, chromosomes, window = 1e6, overlap = 0,
                         weight = rep.int( x = 1, times = length(x)), start.coord = 1) {
  weight   <- sqrt(weight)
  overlap  <- as.integer(overlap)
  window.offset <- as.integer(window - round(window * (overlap / (overlap + 1))))
  chr.ordered <- unique(chromosomes)
  data.splitByChr    <- split(data.frame(pos = positions, x = x, weight = weight), 
                              f = factor(chromosomes, levels = chr.ordered))
  lapply(data.splitByChr, function(data.oneChr) {
    range.pos <- range(data.oneChr$pos, na.rm = TRUE)
    if (!is.null(start.coord)) {
      range.pos[1] <- as.integer(start.coord)
    }
    beam.coords <- seq.int(range.pos[1], range.pos[2], by = window.offset)
    if (max(beam.coords) != range.pos[2] ) {
      beam.coords <- c(beam.coords, range.pos[2])
    }
    nWindows <- length(beam.coords) - overlap - 1
    pos.cut <- cut(data.oneChr$pos, breaks = beam.coords)
    x.split <- split(data.oneChr$x, f = pos.cut)
    weight.split <- split(data.oneChr$weight, f = pos.cut)
    window.starts <- beam.coords[1:nWindows]
    window.ends <- beam.coords[(1:nWindows) + 1 + overlap]
    idx.list <- lapply(1:nWindows, function(ii) ii + (0:overlap))
    x.window <- lapply(idx.list, function(idx) unlist(x.split[idx], use.names = FALSE))
    weight.window <- lapply(idx.list, function(idx) unlist(weight.split[idx], use.names = FALSE))
    window.means <- mapply(weighted.mean, x = x.window, w = weight.window)
    window.quantiles <- sapply(x.window, quantile, probs = c(0.25, 0.75), na.rm = TRUE, names = FALSE)
    window.counts <- sapply(x.window, length)
    data.frame(start = window.starts, end = window.ends, mean = window.means, 
               q0 = window.quantiles[1,], q1 = window.quantiles[2,], N = window.counts)
  })
}

get.ci <- function(cp.table, level = 0.95) {
  znormsort <- sort(cp.table$lpp, decreasing = TRUE)
  znormcumLik <- cumsum(znormsort)
  n <- sapply(level, function(x) sum(znormcumLik < x) + 1)
  LikThresh <- znormsort[n]
  values.x <- data.frame(x = cp.table$ploidy, y = apply(cp.table$lpp, 1, max))
  values.y <- data.frame(x = apply(cp.table$lpp, 2, max), y = cp.table$cellularity)
  up.x  <- max(values.x$x[values.x$y >= LikThresh])
  low.x <- min(values.x$x[values.x$y >= LikThresh])
  max.x <- values.x$x[which.max(values.x$y)]
  up.y  <- max(values.y$y[values.y$x >= LikThresh])
  low.y <- min(values.y$y[values.y$x >= LikThresh])
  max.y <- values.y$y[which.max(values.y$x)]
  results <- list()
  values.x$y <- values.x$y/sum(values.x$y)
  values.y$x <- values.y$x/sum(values.y$x)
  results$values.ploidy <- values.x
  results$confint.ploidy <- c(low.x, up.x)
  results$max.ploidy <- max.x
  results$values.cellularity <- values.y
  results$confint.cellularity <- c(low.y, up.y)
  results$max.cellularity <- max.y
  results
}

mut.fractions <- function(AB.tumor, Af, tumor.strand) {
  F = 1 - Af
   base.mut <- lapply(X = AB.tumor, FUN = function(x) unlist(strsplit(as.character(x), split = '[:]')))
   base.fw  <- lapply(X = tumor.strand, FUN = function(x) unlist(strsplit(as.character(x), split = '[:]')))
   frequencify <- function (x) {
      base.name <- substr(unlist(x), 1, 1)
      base.val  <- as.numeric(substr(unlist(x), 2, nchar(x)))
      setNames(base.val, base.name)
   }
   base.freqs <- lapply(X = base.mut, FUN = frequencify)
   fw.freqs   <- lapply(X = base.fw, FUN = frequencify)
   n.base.mut <- do.call(c, lapply(X = base.mut, FUN = length))
   max.fq <- function (x) {
      freq.rel <- base.freqs[[x]] / F[x]
      f.max    <- which.max(freq.rel)
      c(freq.rel[f.max], names(base.freqs[[x]])[f.max], base.freqs[[x]][f.max], fw.freqs[[x]][f.max])
   }
   max.freqs  <- do.call(rbind, lapply(1:length(F), max.fq))
   data.frame(base.count = as.integer(n.base.mut), maj.base.freq = as.numeric(max.freqs[, 1]),
              base = as.character(max.freqs[,2]), freq = as.numeric(max.freqs[,3]),
              fw.freq = as.numeric(max.freqs[,4]))
}

mutation.table <- function(seqz.tab, mufreq.treshold = 0.15, min.reads = 40, min.reads.normal = 10,
                           max.mut.types = 3, min.type.freq = 0.9, min.fw.freq = 0,
                           segments = NULL) {
   chroms      <- unique(seqz.tab$chromosome)
   hom.filt    <- seqz.tab$zygosity.normal == 'hom'
   seqz.tab     <- seqz.tab[hom.filt, ]
   reads.filt  <- seqz.tab$good.reads >= min.reads & seqz.tab$depth.normal >= min.reads.normal
   seqz.tab     <- seqz.tab[reads.filt, ]
   mufreq.filt <- seqz.tab$Af <= (1 - mufreq.treshold)
   seqz.tab     <- seqz.tab[mufreq.filt, ]
   if (!is.null(segments)) {
      for (i in 1:nrow(segments)) {
         pos.filt <- seqz.tab$chromosome == segments$chromosome[i] & seqz.tab$position >= segments$start.pos[i] & seqz.tab$position <= segments$end.pos[i]
         seqz.tab$adjusted.ratio[pos.filt] <- segments$depth.ratio[i]
      }
   }
   seqz.dummy   <- data.frame(chromosome = chroms, position = 1, GC.percent = NA, good.reads = NA,
                             adjusted.ratio = NA, F = 0, mutation = 'NA', stringsAsFactors= FALSE)
   if (nrow(seqz.tab) >= 1) {
      mu.fracts   <- mut.fractions(AB.tumor = seqz.tab$AB.tumor, Af = seqz.tab$Af,
                                   tumor.strand = seqz.tab$tumor.strand)
      mufreq.filt <- mu.fracts$freq >= mufreq.treshold
      type.filt   <- mu.fracts$base.count <= max.mut.types
      prop.filt   <- mu.fracts$maj.base.freq >= min.type.freq
      if (!is.na(min.fw.freq)) {
         fw.2 = 1 - min.fw.freq
         fw.2 <- sort(c(fw.2, min.fw.freq))
         fw.filt     <- mu.fracts$fw.freq > fw.2[1] & mu.fracts$fw.freq < fw.2[2]
         mufreq.filt <- mufreq.filt & type.filt  & prop.filt & fw.filt
      } else {
         mufreq.filt <- mufreq.filt & type.filt  & prop.filt
      }
      mut.type    <- paste(seqz.tab$AB.normal, mu.fracts$base, sep = '>')
      seqz.tab     <- seqz.tab[,c('chromosome', 'position', 'GC.percent', 'good.reads', 'adjusted.ratio')]
      seqz.tab     <- cbind(seqz.tab, F = mu.fracts$freq, mutation = mut.type)
      rbind(seqz.tab[mufreq.filt, ], seqz.dummy)
   } else {
      seqz.dummy
   }
}

find.breaks <- function(seqz.baf, gamma = 80, kmin = 10, baf.thres = c(0, 0.5),
                        verbose = FALSE, seg.algo = "aspcf", ...) {
   chromosome <- gsub(x = seqz.baf$chromosome, pattern = "chr", replacement = "")
   logR = data.frame(chrom = chromosome,
                     pos = seqz.baf$position,
                     s1 = log2(seqz.baf$adjusted.ratio))
   logR.wins <- copynumber::winsorize(logR, verbose = verbose)
   if (seg.algo == "aspcf"){
      BAF = data.frame(chrom = chromosome,
                       pos = seqz.baf$position,
                       s1 = seqz.baf$Bf)
      allele.seg <- copynumber::aspcf(logR = logR.wins, BAF = BAF, baf.thres = baf.thres,
                          verbose = verbose, gamma = gamma, kmin = kmin, ...)
   } else if (seg.algo == "pcf") {
      allele.seg <- copynumber::pcf(data = logR.wins, verbose = verbose,
                                    gamma = gamma, kmin = kmin, ...)
   } else {
      stop("Supported segmentation algorithms are only \'aspcf\' or \'pcf\' from the copynumber package.")
   }
    if (length(grep("chr", seqz.baf$chromosome)) > 0) {
        allele.seg$chrom <- paste("chr", allele.seg$chrom, sep = "")
    }
   breaks   <- allele.seg[, c("chrom", "start.pos", "end.pos", "arm")]
   not.uniq <- which(breaks$end.pos == c(breaks$start.pos[-1],0))
   breaks$end.pos[not.uniq]<- breaks$end.pos[not.uniq] - 1
   breaks
}

segment.breaks <- function(seqz.tab, breaks, min.reads.baf = 1,
                           weighted.mean = TRUE) {
   if (weighted.mean){
      w.r     <- sqrt(seqz.tab$depth.normal)
      rw      <- seqz.tab$adjusted.ratio * w.r
      w.b     <- sqrt(seqz.tab$good.reads)
      bw      <- seqz.tab$Bf * w.b
      seqz.tab <- cbind(seqz.tab[, c("chromosome", "position", "zygosity.normal", "good.reads")],
                    rw = rw, w.r = w.r, bw = bw, w.b = w.b)
   }
   chr.order <- unique(seqz.tab$chromosome)
   seqz.tab <- split(seqz.tab, f = seqz.tab$chromosome)
   segments <- list()
   for (i in 1:length(seqz.tab)) {
      seqz.b.i    <- seqz.tab[[i]][seqz.tab[[i]]$zygosity.normal == 'het', ]
      seqz.b.i    <- seqz.b.i[seqz.b.i$good.reads >= min.reads.baf, ]
      breaks.i    <- breaks[breaks$chrom == names(seqz.tab)[i], ]
      nb          <- nrow(breaks.i)
      breaks.vect <- do.call(cbind, split.data.frame(breaks.i[,c("start.pos", "end.pos")], f = 1:nb))
      unique.breaks <- function(b, offset = 1) {
         while(any(diff(b) == 0)) {
            b[which(diff(b) == 0) + 1] <- b[diff(b) == 0] + offset
         }
         b
      }
      breaks.vect <- unique.breaks(b = as.numeric(breaks.vect), offset = 1)
      fact.r.i    <- cut(seqz.tab[[i]]$position, breaks.vect)
      fact.b.i    <- cut(seqz.b.i$position, breaks.vect)
      seg.i.s.r   <- sapply(X = split(seqz.tab[[i]]$chromosome, f = fact.r.i), FUN = length)
      seg.i.s.b   <- sapply(X = split(seqz.b.i$chromosome, f = fact.b.i), FUN = length)
     
      if (weighted.mean){
         seg.i.rw    <- sapply(X = split(seqz.tab[[i]]$rw, f = fact.r.i), FUN = function(a) sum(a, na.rm = TRUE))
         seg.i.w.r   <- sapply(X = split(seqz.tab[[i]]$w.r, f = fact.r.i), FUN = function(a) sum(a, na.rm = TRUE))
         seg.i.r.sd  <- sapply(X = split(seqz.tab[[i]]$rw/seqz.tab[[i]]$w.r, f = fact.r.i), FUN = function(a) sd(a, na.rm = TRUE))
         seg.i.b.sd  <- sapply(X = split(seqz.b.i$bw/seqz.b.i$w.b, f = fact.b.i), FUN = function(a) sd(a, na.rm = TRUE))         
         seg.i.bw    <- sapply(X = split(seqz.b.i$bw, f = fact.b.i), FUN = function(a) sum(a, na.rm = TRUE))
         seg.i.w.b   <- sapply(X = split(seqz.b.i$w.b, f = fact.b.i), FUN = function(a) sum(a, na.rm = TRUE))
         segments.i <- data.frame(chromosome  = names(seqz.tab)[i], start.pos = as.numeric(breaks.vect[-length(breaks.vect)]),
                               end.pos = as.numeric(breaks.vect[-1]), Bf = seg.i.bw/seg.i.w.b, N.BAF = seg.i.s.b, sd.BAF = seg.i.b.sd,
                               depth.ratio = seg.i.rw/seg.i.w.r, N.ratio = seg.i.s.r, sd.ratio = seg.i.r.sd, stringsAsFactors = FALSE)
      } else {
        seg.i.r    <- sapply(X = split(seqz.tab[[i]]$adjusted.ratio, f = fact.r.i), FUN = function(a) mean(a, na.rm = TRUE))
        seg.i.b    <- sapply(X = split(seqz.b.i$Bf, f = fact.b.i), FUN = function(a) mean(a, na.rm = TRUE))
        seg.i.r.sd <- sapply(X = split(seqz.tab[[i]]$adjusted.ratio, f = fact.r.i), FUN = function(a) sd(a, na.rm = TRUE))
        seg.i.b.sd <- sapply(X = split(seqz.b.i$Bf, f = fact.b.i), FUN = function(a) sd(a, na.rm = TRUE)) 
        segments.i <- data.frame(chromosome  = names(seqz.tab)[i], start.pos = as.numeric(breaks.vect[-length(breaks.vect)]),
                                 end.pos = as.numeric(breaks.vect[-1]), Bf = seg.i.b, N.BAF = seg.i.s.b, sd.BAF = seg.i.b.sd,
                                 depth.ratio = seg.i.r, N.ratio = seg.i.s.r, sd.ratio = seg.i.r.sd, stringsAsFactors = FALSE)
      }
      segments[[i]] <- segments.i[seq(from = 1, to = nrow(segments.i), by = 2),]
   }
   segments <- do.call(rbind, segments[as.factor(chr.order)])
   row.names(segments) <- 1:nrow(segments)
   len.seg <- (segments$end.pos - segments$start.pos) / 1e6
   segments[(segments$N.ratio/len.seg) >= 2, ]
}

alternative.cp.solutions <- function(cp.table) {
   ci <- get.ci(cp.table)
   p.alt <- which(diff(sign(diff(ci$values.ploidy$y))) == -2) + 1
   get.alt <- function(idx.p, cp.table) {
      idx.c <- which.max(cp.table$lpp[idx.p,])
      c(cellularity = cp.table$cellularity[idx.c],
        ploidy = cp.table$ploidy[idx.p],
        SLPP = cp.table$lpp[idx.p, idx.c])
   }
   res <- lapply(p.alt, FUN = function (x) get.alt(x, cp.table))
   res <- as.data.frame(do.call(rbind, res))
   if (nrow(res) > 0 ){
      res[order(res$SLPP, decreasing = TRUE), ]
   } else {
      data.frame(cellularity = ci$max.cellularity, 
                 ploidy = ci$max.ploidy,
                 SLPP =  cp.table$lpp[which(cp.table$ploidy == ci$max.ploidy),
                                      which(cp.table$cellularity == ci$max.cellularity)])
   }
}
