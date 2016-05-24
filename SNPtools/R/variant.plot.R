################################################################################
# Plot the variant alleles.
################################################################################
# Arguments: variants: data.frame, with ID, CHROM, POS, REF and ALT in the first
#                  five columns and the numeric allele calls for the strains in 
#                  the remaining columns. Variants in rows.
#            col: character vector, with three colors for missing, ALT and REF
#                 alleles.
snp.plot = function(variants, col = c("black", "grey50", "white"), cluster = T, 
           ref, highlight, pattern.snps, mgi, qtl) {

  old.par = par(no.readonly = T)

  type = attr(variants, "type")

  # Verify that the reference strain is in the list of strains that we have.
  if(!missing(ref)) {
    if(!ref %in% colnames(variants)) {
      stop(paste("The reference strain (", ref, 
           ") is not in the colnames of variants."), sep = "")
    } # if(!ref %in% colnames(variants))
  } # if(!missing(ref))

  # If there are quality values, strip them out.
  variants = strip.quality.columns(variants)

  # Cluster the strains, if needed.
  if(cluster) {    
    variants = cluster.strains(variants)
  } # if(cluster)

  # Separate the header from the variants.
  hdr = variants[,1:5]
  variants = as.matrix(variants[,-1:-5])

  # If we have a reference strain, change the matrix to reflect this.
  if(!missing(ref)) {
    variants = matrix(as.numeric(variants != variants[,ref]), nrow(variants),
               ncol(variants), dimnames = dimnames(variants))
  } # if(!missing(ref))

  # Convert variant positions to Mb, if needed.
  if(type %in% c("snp", "indel")) {
    if(max(hdr$POS) > 200) {
      hdr$POS = hdr$POS * 1e-6
    } # if(max(hdr$POS) > 200)
  } else if(type == "sv") {
    if(max(hdr$START) > 200) {
      hdr$START = hdr$START * 1e-6
    } # if(max(hdr$START) > 200)

    if(max(hdr$END) > 200) {
      hdr$END = hdr$END * 1e-6
    } # if(max(hdr$END) > 200)
  } # else if(type == "sv")

  # If we have a QTL data.frame, convert the positions to Mb, if needed.
  if(!missing(qtl)) {
    if(ncol(qtl) != 3) {
      stop(paste("The QTL data.frame must contain three columns with chr,",
           "position and score."))
    } # if(ncol(qtl) != 3)

    if(max(qtl[,2])> 200) {
      qtl[,2] = qtl[,2] * 1e-6
    } # if(max(qtl[,2])> 200)
  } # if(!missing(qtl))

  # Flip the SNPs so that they plot with the first column at the top.
  variants = variants[,ncol(variants):1]

  # Make NA = -1 for plotting.
  variants[is.na(variants)] = -1

  # If we have pattern SNPs and genes, then look for their intersection
  # and recored which genes to color orange.
  genes.to.mark = NULL
  if(!missing(pattern.snps) && !missing(mgi)) {
    if(max(pattern.snps$POS) > 200) {
      pattern.snps$POS = pattern.snps$POS * 1e-6
    } # if(max(pattern.snps$POS) > 200)
    genes.to.mark = rep(F, nrow(mgi))
    for(i in 1:nrow(mgi)) {
      variants.in.gene = which(pattern.snps$POS >= mgi$start[i] * 1e-6 & 
                           pattern.snps$POS <= mgi$stop[i] * 1e-6)
      genes.to.mark[i] = length(variants.in.gene) > 0
    } # for(i)
    genes.to.mark = mgi$Name[genes.to.mark]
  } # if(!missing(snp.genes))

  # Plot the SNPs.
  # Make a two panel plot if MGI genes are provided.
  rplt = 0.15
  lplt = 0.95
  if(!missing(mgi)) {
    if(!missing(qtl)) {
      if(!missing(pattern.snps)) {
        layout(matrix(1:3, nrow = 3), heights = c(0.3, 0.1, 0.6))
        par(plt = c(rplt, lplt, 0.05, 0.75))
      } else {
        layout(matrix(1:3, nrow = 3), heights = c(0.3, 0.1, 0.6))
        par(plt = c(rplt, lplt, 0.01, 0.75))
      } # else
    } else {
      if(!missing(pattern.snps)) {
        layout(matrix(1:2, nrow = 2), heights = c(0.3, 0.7))
        par(plt = c(rplt, lplt, 0.07, 0.88))
      } else {
        layout(matrix(1:2, nrow = 2), heights = c(0.3, 0.7))
        par(plt = c(rplt, lplt, 0.01, 0.75))
      } # else
    } # else
  } else {
    if(!missing(qtl)) {
      if(!missing(pattern.snps)) {
        layout(matrix(1:2, nrow = 2), heights = c(0.9, 0.1))
        par(plt = c(rplt, lplt, 0.05, 0.88))
      } else {
        layout(matrix(1:2, nrow = 2), heights = c(0.9, 0.1))
        par(plt = c(rplt, lplt, 0.01, 0.88))
      } # else
    } else {
      if(!missing(pattern.snps)) {
        par(plt = c(rplt, lplt, 0.05, 0.88))
      } else {
        par(plt = c(rplt, lplt, 0.05, 0.88))
      } # else
    } # else
  } # else 

  x = NULL
  if(type %in% c("snp", "indel")) {
    x = hdr$POS
  } else if(type == "sv") {
    x = hdr$START
  } #  else if(type == "sv")

  image(x, 1:ncol(variants), variants, breaks = c(-1.5, -0.5, 0.5, 1.5),
        col = col, ann = F, axes = F)
  par(las = 1)
  saved.xlim = par("usr")[1:2]
  mtext(text = paste("Chr", hdr[1, grep("CHR", colnames(hdr))], "(Mb)"), 
        side = 3, line = 2)
  abline(h = 1:ncol(variants) - 0.5, col = "grey80")
#  par(cex.lab = 1.5, cex.axis = 1.5)
#  abline(v = axis(3, padj = 0.4), lty = 2)
  axis(3, padj = 0.4)
  usr = par("usr")

  # Plot the strain names.
  par(las = 2)
  mtext(text = colnames(variants), side = 2, line = 0.2, at = 1:ncol(variants))
  rect(usr[1], usr[3], usr[2], usr[4], lwd = 2)

  # If the user has given strains to highlight, draw a transparent green
  # box over them.
  if(!missing(highlight)) {
    hcols = match(highlight, colnames(variants))
    if(any(is.na(hcols))) {
      stop(paste("All of the highlight strains are not in the SNP matrix.",
           paste(highlight[is.na(hcols)], collapse = ",")))
    } # if(any(is.na(hcols)))

    for(i in 1:length(hcols)) {
      rect(usr[1], hcols[i] - 0.5, usr[2], hcols[i] + 0.5,
           col = rgb(1, 0, 0, 0.1))
    } # for(i)
  } # if(!missing(highlight))

  # If the user provided pattern SNPs, plot them below the SNP tile plot.
  if(!missing(pattern.snps)) {
    par(xpd = NA)
    if(max(pattern.snps$POS) > 200) {
      pattern.snps$POS = pattern.snps$POS * 1e-6
    } # if(max(pattern.snps$POS) > 200)

    coord = par(c("usr", "fin", "plt"))
    ex = (coord$usr[4] - coord$usr[3]) / (coord$plt[4] - coord$plt[3])
    top = coord$usr[3]
    bottom = coord$usr[3] - ex * coord$plt[3]
    x = matrix(rep(pattern.snps$POS, each = 2), nrow = 2)
    x = rbind(x, rep(NA, ncol(x)))
    y = rbind(rep(top, ncol(x)), rep(bottom, ncol(x)), rep(NA, ncol(x)))
    lines(as.vector(x), as.vector(y), col = rgb(1, 0.6, 0))
    par(xpd = NA)
    text(coord$usr[1], (top + bottom) * 0.5, "SNP Pattern ", adj = 1)
    par(xpd = F)
  } # if(!missing(pattern.snps))

  # If the user has given us QTL values, plot them below the SNPs.
  if(!missing(qtl)) {
    qtl = qtl[qtl[,1] == hdr$CHROM[1] & qtl[,2] >= hdr$POS[1] - 1 &
              qtl[,2] <= hdr$POS[nrow(hdr)] + 1,]

    xout = seq(hdr$POS[1], hdr$POS[nrow(hdr)], length.out = 1000)
    val = approx(qtl[,2], qtl[,3], xout = xout)
    val$y[is.na(val$y)] = min(val$y, na.rm = T)
    br = quantile(val$y, 0:100/100)
    cl = colorRampPalette(c(rgb(0,0,0), rgb(1,0,0)))(length(br) - 1)
    par(plt = c(rplt, lplt, 0.1, 0.9))
    image(1:length(val$x), 1, as.matrix(val$y), breaks = br, col = cl,
          ann = F, axes = F)
    mtext(text = "QTL", side = 2, line = 0.2)

    # Plot the SNP locations where mapping occured.
#    par(xpd = NA)
#    coord = par(c("usr", "fin", "plt"))
#    ex = (coord$usr[4] - coord$usr[3]) / (coord$plt[4] - coord$plt[3])
#    top    = coord$usr[4] + ex * (1.0 - coord$plt[4])
#    bottom = coord$usr[4] 
#    x = qtl[qtl[,2] >= hdr$POS[1] & qtl[,2] <= hdr$POS[nrow(hdr)],2]
#    x = matrix(rep(x, each = 2), nrow = 2)
#    x = rbind(x, rep(NA, ncol(x)))
#    y = rbind(rep(bottom, ncol(x)), rep(top, ncol(x)), rep(NA, ncol(x)))
#    lines(as.vector(x), as.vector(y), lwd = 3, lend = 2)
#    text(coord$usr[1], (top + bottom) * 0.5, "Mapped ", adj = 1)
#    par(xpd = F)
  } # if(!missing(qtl))

  # Draw MGI features, if requested (by including mgi).
  if(!missing(mgi) && nrow(mgi) > 0) {
    gene.cols = rep(rgb(0.4, 0.4, 0.9), nrow(mgi))
    if(!missing(genes.to.mark)) {
      gene.cols[mgi$Name %in% genes.to.mark] = rgb(1, 0.6, 0)
    } # if(!missing(genes.to.mark))
    
    par(plt = c(rplt, lplt, 0.15, 1.0), las = 1)
# , cex.lab = 1.5, cex.axis = 1.5
    gene.plot(mgi, col = gene.cols, xlim = saved.xlim, xaxs = "i")
  } # if(!missing(mgi))

  par(old.par)
} # snp.plot()


################################################################################
# High level function to allow plotting of variant alleles, genes, variant pattern
# and QTL in one command.
# Arguments: var.file: character, location of variant file.
#            mgi.file: character, location of variantfile.
#            chr: character or numeric, the Chr to plot on.
#            start: numeric, the bp or Mb location to start plotting.
#            end: numeric, the bp or Mb location to stop plotting.
#            strains: character vector, the strains to plot. Default = CC founders.
#            ref: character, one of the strains above that should be the 
#                 references strain for SNP plotting. Default = C57BL/6J
#            pattern: character vector, the strains to use in marking the SNP
#                     pattern.
#            qtl: data.frame with three columns containing chr, position and
#                 mapping statistic in that order.
# Returns: data.frame with SNPs that match that allele effect pattern
################################################################################
variant.plot = function(var.file = "http://cgd.jax.org/tools/SNPtools/Build38/sanger.snps.NCBI38.txt.gz",
               mgi.file = "http://cgd.jax.org/tools/SNPtools/MGI/MGI.20130305.sorted.txt.gz",
               chr, start, end, type = c("snp", "indel", "sv"), strains = 
               c("A/J", "C57BL/6J", "129S1/SvImJ", "CAST/EiJ", "NOD/ShiLtJ",
               "NZO/HlLtJ", "PWK/PhJ", "WSB/EiJ"), ref = "C57BL/6J", pattern, qtl)
{
  print("Checking arguments...")

  if(missing(chr)) {
    stop(paste("The chr argument is NULL. Please enter the chromosome and",
         "start and end locations that you would like to plot."))
  } # if(missing(chr))

  if(missing(start)) {
    stop(paste("The start argument is NULL. Please enter the chromosome and",
         "start and end locations that you would like to plot."))
  } # if(missing(start))

  if(missing(end)) {
    stop(paste("The end argument is NULL. Please enter the chromosome and",
         "start and end locations that you would like to plot."))
  } # if(missing(end))

  if(!missing(ref)) {
    if(!ref %in% strains) {
      stop(paste("The refrence strain (", ref, ")is not in the strains.",
      "Please select a reference strain from the strains you are plotting."))
    } # if(!ref %in% strains)
  } # if(!missing(ref))

  if(!missing(qtl)) {
    if(ncol(qtl) != 3) {
      stop(paste("The QTL data.frame must contain three columns with chr,",
           "position and mapping statistic, in that order."))
    } # if(ncol(qtl) != 3)
  } # if(!missing(qtl))

  type = match.arg(type)

  print("Retreiving Variants...")
  variants = get.variants(file = var.file, chr = chr, start = start, 
         end = end, type = type, strains = strains, polymorphic = F)
  print(paste(nrow(variants), "Variants in region."))
  nvariants = convert.variants.to.numeric(variants)
  rm(variants)

  cat.variants = NULL
  if(!missing(pattern)) {
    print("Getting SNPs that match the pattern...")
    pattern.snps = get.pattern.variants(nvariants, pattern)
    print(paste(nrow(pattern.snps), "variants fit the pattern."))
    print("Categorizing pattern variants...")
    cat.variants = categorize.variants(pattern.snps, mgi.file)
    if(!is.null(cat.variants)) {
      gene.models = grep("^Gm", cat.variants$symbol)
      if(length(gene.models) > 0) {
        cat.variants = cat.variants[-gene.models,]
      } # if(length(gene.models) > 0)
      print(paste(length(unique(cat.variants$symbol[!is.na(cat.variants$symbol)])),
            "genes have variants that fit the pattern",
            "between the 3'UTR and the 5'UTR (including introns)."))
    } else {
     print("No variants match the pattern.")
    } # else
  } # if(!missing(pattern))

  print("Getting genes in region...")
  mgi = get.mgi.features(mgi.file, chr = chr, start = start, end = end,
                       source = "MGI", type = "gene")
  gene.models = grep("^Gm", mgi$Name)
  if(length(gene.models) > 0) {
    mgi = mgi[-gene.models,]
  } # if(length(gene.models) > 0)
  print(paste(nrow(mgi), "genes in region."))

  print("Drawing plot...")
  snp.plot(variants = nvariants, pattern.snps = pattern.snps, mgi = mgi, qtl = qtl)

  if(!missing(cat.variants)) {
    return(cat.variants)
  } # if(!missing(cat.variants))

} # variant.plot()

