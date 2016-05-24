subclonal.matrix <- function(mut.tab, segments = NULL, cellularity = seq(0.1, 1, 0.05), mc.cores = getOption("mc.cores", 2L)){
  if (!is.null(segments)) {
    mut.tab$CNt <- 2
    for (i in 1:nrow(segments)) {
      pos.filt <- mut.tab$chromosome == segments$chrom[i] & mut.tab$position >= segments$start.pos[i] & mut.tab$position <= segments$end.pos[i]
      mut.tab$CNt[pos.filt] <- segments$CNt[i]
    }
  }
  mut.types.list <- lapply(X = 1:nrow(mut.tab),
                           FUN = function(x) {
                             types.matrix(CNn = mut.tab[x, 'CNn'],
                                          CNt.min = mut.tab[x, 'CNt'],
                                          CNt.max = mut.tab[x, 'CNt'])})
  mut.clonality <- function(F, depth, types, cellularity) {
    theoretical.F <- theoretical.mufreq(cellularity = cellularity, 
                                  CNn = types[, 1], CNt = types[, 2], Mt = types[, 3])
    max(mufreq.dpois(mufreq = F, mufreq.model = theoretical.F, depth.t = depth), na.rm = TRUE)
  }
  res <- mclapplyPb (mc.cores = mc.cores, X = 1:nrow(mut.tab),
                     FUN = function (i) {
                             sapply(X = cellularity, FUN = function(x) {
                                mut.clonality(F = mut.tab$F[i],
                                           depth = mut.tab$good.reads[i],
                                           types = mut.types.list[[i]],
                                           cellularity = x)
                           })
                         })
  res <- do.call(rbind, res)
  res <- res / rowSums(res)
  colnames(res) <- cellularity
  res
}

sequenza2PyClone <- function(mut.tab, seg.cn, sample.id, norm.cn = 2) {
  mut.tab <- cbind(mut.tab[,c('chromosome', 'position', 'good.reads','F', 'mutation')], CNt = NA, A = NA, B = NA)
  for (i in 1:nrow(seg.cn)) {
     pos.filt <- mut.tab$chromosome == seg.cn$chromosome[i] & mut.tab$position >= seg.cn$start.pos[i] & mut.tab$position <= seg.cn$end.pos[i]
     mut.tab[pos.filt, c("CNt", "A", "B")] <- seg.cn[i, c("CNt", "A", "B")]
  }
  id          <- paste(sample.id, mut.tab$chromosome, mut.tab$position, sep = ':')
  var.counts  <- round(mut.tab$good.reads * mut.tab$F, 0)
  nor.counts  <- mut.tab$good.reads - var.counts
  pyclone.tsv <- data.frame(mutation_id = id, ref_counts = nor.counts, var_counts = var.counts,
                            normal_cn = norm.cn,	minor_cn = mut.tab$B, major_cn = mut.tab$A,
                            variant_case = sample.id,	variant_freq = mut.tab$F,	genotype = mut.tab$mutation)
  na.exclude(pyclone.tsv)
}

VarScan2seqz <- function(varscan.somatic, varscan.copynumber = NULL, normal_var_freq = 0.25) {
   
   iupac.nucs     <- setNames(c('A', 'C', 'G', 'GT', 'AC', 'AG', 'CG', 'T', 'AT', 'CT'),
                              c('A', 'C', 'G', 'K', 'M', 'R', 'S', 'T', 'W', 'Y'))
   zygosity.vect  <- setNames(c('hom', 'hom', 'hom', 'het', 'het', 'het', 'het', 'hom', 'het', 'het'),
                              c('A', 'C', 'G', 'K', 'M', 'R', 'S', 'T', 'W', 'Y'))
   varscan.somatic  <- varscan.somatic[varscan.somatic$somatic_status != 'Unknown', ]
   varscan.somatic$normal_var_freq <- as.numeric(sub('%', '', varscan.somatic$normal_var_freq)) / 100
   varscan.somatic$tumor_var_freq  <- as.numeric(sub('%', '', varscan.somatic$tumor_var_freq)) / 100
   zygosity.normal <- zygosity.vect[varscan.somatic$normal_gt]
   if (normal_var_freq > 0.5){
      normal_var_freq <- 1 - normal_var_freq
   }
   idx             <- which(zygosity.normal == "het" &
                   (normal_var_freq >= varscan.somatic$normal_var_freq |
                   varscan.somatic$normal_var_freq >= (1-normal_var_freq)))
   varscan.somatic <- varscan.somatic[-idx, ]
   zygosity.normal <- zygosity.vect[varscan.somatic$normal_gt]
   AB.normal    <- iupac.nucs[varscan.somatic$normal_gt]
   AB.tumor     <- rep('.', length(AB.normal))
   strand       <- AB.tumor
   depth.normal <- varscan.somatic$normal_reads1 + varscan.somatic$normal_reads2
   depth.tumor  <- varscan.somatic$tumor_reads1 + varscan.somatic$tumor_reads1
   depth.ratio  <- depth.tumor/depth.normal
   Af <- 1 - varscan.somatic$tumor_var_freq
   Bf <- rep(0, length(Af))
   idx <- zygosity.normal == 'het' & Af  < 0.5
   Af[idx] <- 1 - Af[idx]
   idx <- zygosity.normal == 'het'
   Bf[idx]  <- 1 - Af[idx]
   ## Overwrite the Genotype base orders to
   ## support pseudo haplotype analysis.
   genotype_1     <- varscan.somatic[idx, c("ref", "var")]
   genotype_1     <- apply(genotype_1, 1, paste, collapse = "")
   AB.normal[idx] <- genotype_1
   idx            <- idx & varscan.somatic$tumor_var_freq > 0.5
   genotype_2     <- varscan.somatic[idx, c("var", "ref")]
   genotype_2     <- apply(genotype_2, 1, paste, collapse = "")
   AB.normal[idx] <- genotype_2
   ## Not the best way, but I'm running out of ideas
   idx <- zygosity.normal == 'hom' & varscan.somatic$somatic_status == 'Somatic'
   if (sum(idx) > 0) {
      mut.b <- cbind(as.character(iupac.nucs[varscan.somatic$tumor_gt[idx]]),
                     as.character(varscan.somatic$normal_gt[idx]))
      mut.b <- sapply(X = 1:sum(idx),
                      FUN = function(x) gsub(x = mut.b[x, 1],
                                        pattern = mut.b[x, 2],
                                        replacement = ''))
      mut   <- paste0(mut.b, varscan.somatic$tumor_var_freq[idx])
      strand[idx]   <- paste0(mut.b,
                              varscan.somatic$tumor_reads2_plus[idx]/(varscan.somatic$tumor_reads2_plus[idx]+varscan.somatic$tumor_reads2_minus[idx]))
      AB.tumor[idx] <- mut
   }
   res <- data.frame(chromosome = as.character(varscan.somatic$chrom), position = varscan.somatic$position,
                     base.ref = as.character(varscan.somatic$ref), depth.normal = depth.normal,
                     depth.tumor = depth.tumor, depth.ratio = depth.ratio,
                     Af = round(Af, 3), Bf = round(Bf, 3), zygosity.normal = zygosity.normal, GC.percent = 50,
                     good.reads = round(depth.tumor, 2), AB.normal = AB.normal,
                     AB.tumor = AB.tumor, tumor.strand = strand, stringsAsFactors = FALSE)
   normal.pos <- res$zygosity.normal == 'hom' & res$AB.tumor == '.'
   res <- res[res$depth.ratio > 0 & !is.infinite(res$depth.ratio) & !normal.pos, ]
   if (!is.null(varscan.copynumber)){
      smart.id <- order(c(1:nrow(varscan.copynumber), 1:nrow(varscan.copynumber) + 0.5))
      varscan.copynumber$log2_ratio   <- 2 ^ (varscan.copynumber$log2_ratio)
      varscan.copynumber$normal_depth <- round(varscan.copynumber$normal_depth, 0)
      varscan.copynumber$tumor_depth <- round(varscan.copynumber$tumor_depth, 0)
      varscan.copynumber$chrom <- as.character(varscan.copynumber$chrom)
      mat.t <- data.frame(chromosome = c(varscan.copynumber$chrom, varscan.copynumber$chrom)[smart.id],
                          position = c(varscan.copynumber$chr_start, varscan.copynumber$chr_stop)[smart.id], base.ref = 'N',
                          depth.normal = c(varscan.copynumber$normal_depth, varscan.copynumber$normal_depth)[smart.id],
                          depth.tumor = c(varscan.copynumber$tumor_depth, varscan.copynumber$tumor_depth)[smart.id],
                          depth.ratio  = c(varscan.copynumber$log2_ratio, varscan.copynumber$log2_ratio)[smart.id],
                          Af = 1, Bf = 0, zygosity.normal = 'hom',
                          GC.percent = c(varscan.copynumber$gc_content, varscan.copynumber$gc_content)[smart.id],
                          stringsAsFactors = FALSE)
      mat.t <- cbind(mat.t, good.reads = mat.t$depth.tumor, AB.normal = 'N', AB.tumor = '.', tumor.strand = '.')
      chrom.order <- unique(mat.t$chromosome)
      l.cnv <- split(mat.t, mat.t$chromosome)
      l.snp <- split(res, res$chromosome)     
      for (i in names(l.snp)){
         tab.i <- rbind(l.cnv[[i]], l.snp[[i]])
         tab.i <- tab.i[order(tab.i$position), ]
         #idx.i <- diff(tab.i$position)
         #dups.i  <- which(idx.i == 0)
         #snp.i <- which(tab.i$GC.percent == -100)
         l.cnv[[i]] <- tab.i
      }
      res <- do.call(rbind, l.cnv[chrom.order])
   }
   res
}