################################################################################
# Grobs
################################################################################
# create gene grobs
gene_grob <- function(gene, head_len=200, i=0, ...){
  if (!is.dna_seg(gene)) stop("A dna_seg object is required")
  if (nrow(gene) > 1) stop ("gene must be single-row")
  mid <- (gene$start + gene$end)/2
  name <- paste("seg.", i, ".", gene$name, sep="")
  # arrows
  if (gene$gene_type == "arrows" || gene$gene_type == "headless_arrows"){
    if (gene$gene_type == "arrows"){
      arrow <- arrow_coord(x1=gene$start, x2=gene$end,
                           y=0.5, strand=gene$strand, head_len=head_len)
    }
    else {
      arrow <- block_coord(gene$start, gene$end, strand=1, y=0.25)
    }
    grob <- polygonGrob(arrow$x, arrow$y, name=name,
                        gp=gpar(fill=gene$col, lty=gene$lty, lwd=gene$lwd),
                        default.units="native")
  }
  # blocks
  else if (gene$gene_type == "blocks" || gene$gene_type == "side_blocks"){
    if (gene$gene_type == "side_blocks"){
      block <- block_coord(gene$start, gene$end, strand=gene$strand, y=0.5)
    }
    else {
      block <- block_coord(gene$start, gene$end, strand=2, y=0)
    }
    grob <- polygonGrob(block$x, block$y, name=name,
                        gp=gpar(fill=gene$col, col=gene$col, lty=gene$lty,
                          lwd=gene$lwd),
                        default.units="native")  
  }
  # lines
  else if (gene$gene_type == "lines" || gene$gene_type == "side_lines"){
    x <- c(gene$start, gene$end)
    if (gene$gene_type == "side_lines") {
      y <- gene$strand/4 + 0.5
    }
    else {
      y <- 0.5
    }
    grob <- segmentsGrob(x0=gene$start, y0=y, x1=gene$end, y1=y, name=name,
                         gp=gpar(col=gene$col, lwd=gene$lwd, lty=gene$lty),
                         default.units="native")
  }
  # exons
  else if (gene$gene_type == "exons" || gene$gene_type == "side_exons"){
    if (gene$gene_type == "side_exons"){
      block <- exon_coord(gene$start, gene$end, gene$strand)
    }
    else {
      block <- exon_coord(gene$start, gene$end, 0)
    }
    grob <- polygonGrob(block$x, block$y, name=name,
                        gp=gpar(fill=gene$col, col=gene$col, lty=gene$lty,
                          lwd=gene$lwd),
                        default.units="native")  
  }
  
  # bars
  else if (gene$gene_type == "bars" || gene$gene_type == "side_bars") {
    if (gene$gene_type == "side_bars"){
      y0 <- 0.5; y1 <- 0.5+gene$strand/2
    }
    else {
      y0 <- 0; y1 <- 1
    }
    grob <- segmentsGrob(x0=mid, y0=y0, x1=mid, y1=y1, name=name,
                         gp=gpar(col=gene$col, lwd=gene$lwd, lty=gene$lty),
                         default.units="native")
  }

  # introns
  else if (gene$gene_type == "introns") {
    intron <- list(x0=c(gene$start, mid), x1=c(mid, gene$end),
                   y0=c(0.3, 0.5)*gene$strand + 0.5,
                   y1=c(0.5, 0.3)*gene$strand + 0.5)
    grob <- segmentsGrob(x0=intron$x0, y0=intron$y0,
                         x1=intron$x1, y1=intron$y1, name=name,
                         gp=gpar(lty=gene$lty, lwd=gene$lwd),
                         default.units="native")
  }
  # points
  else if (gene$gene_type == "points" || gene$gene_type == "side_points") {
    if (gene$gene_type == "side_points"){
      y <- 0.5+gene$strand/4
    }
    else {
      y <- 0.5
    }
    grob <- pointsGrob(x=mid, y=y, name=name,
                       pch=gene$pch, size=unit(gene$cex/2, "char"),
                       gp=gpar(col=gene$col), default.units="native")
  }
  # text
  else if (gene$gene_type == "text" || gene$gene_type == "side_text") {
    if (gene$gene_type == "side_text"){
      just <- c("centre", c("top", "bottom")[gene$strand/2 + 1.5])
    }
    else {
      just <- c("centre", "centre")
    }
    grob <- textGrob(label=gene$name, x=mid, y=0.5, name=name,
                     just=just, gp=gpar(col=gene$col, cex=gene$cex),
                     default.units="native")
  }
  else {
    
    grob <- try(do.call(gene$gene_type, list(gene, ...)), silent=FALSE)
    #browser()
    if (!(is.grob(grob) || all(sapply(grob, is.grob))))
      stop(paste(gene$gene_type, "is an invalid gene_type or",
                 "does not return a grob"))
  }
  grob
}
# create dna_seg grobs
dna_seg_grob <- function(dna_seg, ...){
  if(!is.dna_seg(dna_seg)) stop("A dna_seg object is required")
  grob_list <- gList()
  if (nrow(dna_seg) < 1) return(grob_list)
  for (i in 1:nrow(dna_seg)){
    gene <- as.dna_seg(dna_seg[i,])
    grob_list[[i]] <- gene_grob(gene, ...)
  }
  grob_list
}
# create similarity grobs
similarity_grob <- function(similarity, i){
  if (!is.comparison(similarity)) stop("A comparison object is required")  
  if (nrow(similarity) > 1) stop ("gene must be single-row")
  if (is.null(similarity$col)) similarity$col <- grey(0.5)
  x1 <- c(similarity$start1, similarity$end1)
  x2 <- c(similarity$end2, similarity$start2)
  polygonGrob(x=c(x1, x2), y=c(1, 1, 0, 0),
              name=paste("comp.", i, ".",
                similarity$start1, "-", similarity$end1, "_",
                similarity$start2, "-", similarity$end2, sep=""),
              gp=gpar(fill=similarity$col, col=similarity$col, lwd=0.1),
              default.units="native")
}
# create comparisons grobs
comparison_grob <- function(comparison, ...){
  if (!is.comparison(comparison)) stop("A comparison object is required")
  grob_list <- gList()
  if (nrow(comparison) < 1) return(grob_list)
  # go in the reverse order to plot strongest comparison last
  for (i in 1:nrow(comparison)){
    grob_list[[i]] <- similarity_grob(comparison[i,], ...)
  }
  grob_list
}
# create annot grob
label_grob <- function(label, cex=0.8){
  y <- 0
  w <- 0.1
  if (!is.annotation(label)) stop("An annotation object is required")
  if (nrow(label) > 1) stop("A single-line annotation is required")
  grob_list <- gList()
  # range
  if (!is.na(label$x2)){
    bracket_coord <- bracket_coord(label$x1, label$x2, y=y, w=w)
    grob_list[[2]] <- linesGrob(x=bracket_coord$x, y=bracket_coord$y,
                                name=paste("annot", "line",
                                  gsub(" ", "_", label$text), sep="."),
                                default.units="native",
                                gp=gpar(col=label$col))
    x <- mean(c(label$x1, label$x2))
    w <- w*2
  } else {
    x <- label$x1
  }
  if (label$rot == 0){
    just <- c(0.5, 0)
  } else {
    just <- c(-0.1, 0.5)
  }
  grob_list[[1]] <- textGrob(label$text, x=x, y=y+w, just=just,
                             name=paste("annot", "label",
                               gsub(" ", "_", label$text), sep="."),
                             rot=label$rot,
                             default.units="native",
                             gp=gpar(col=label$col, cex=cex))
  grob_list
}
# create annotation grob
annotation_grob <- function(annotation, ...){
  if (!is.annotation(annotation)) stop("An annotation object is required")
  grob_list <- gList()
  if (nrow(annotation) < 1) return(grob_list)
  for (i in 1:nrow(annotation)){
    label <- as.annotation(annotation[i,])
    grob_list[[i]] <- label_grob(label, ...)
  }
  grob_list
}
# create tree grob
dna_seg_label_grob <- function(labels, cex, col){
  n_label <- length(labels)
  y <- seq(1, 0, len=n_label)
  labelGrobs <- gList()
  for (i in 1:n_label) {
    labelGrobs[[i]] <-
      textGrob(x=0, y=y[i], name=paste("label", i, sep="."),
               label=labels[i], just="left", gp=gpar(cex=cex, col=col[i]),
               default.units="native")
  }
  width <- unit(1, "grobwidth", labelGrobs[[which.max(nchar(labels))]])
  labelTree <- gTree(children=labelGrobs,
                     vp=viewport(xscale=c(0, 1), yscale=c(0, 1),
                       width=width, name="labels"), name="labelsTree")
  list(grob=labelTree, width=width)
}
# create scale grob
scale_grob <- function(max_length){
  rng <- diff(pretty(c(0, max_length), n=8))[1]
  gList(segmentsGrob(x0=max_length, y0=0, x1=max_length-rng, y1=0,
                     name="scale.lines", default.units="native"),
        textGrob(label=human_nt(rng)$text, x=max_length-rng/2, y=0.5,
                 name="scale.text", default.units="native")
        )
}
# create dna_seg scale grob
dna_seg_scale_grob <- function(range, cex=0.6, unit, i, j){
  range <- as.numeric(range)
  if (length(range) != 2 && !is.numeric(range))
    stop("range must be numeric and length 2")
  x0 <- ceiling(range[1]/unit)*unit
  x1 <- floor(range[2]/unit)*unit
  if (x1 < x0) x1 <- x0
  ticks <- seq(x0, x1, by=unit)
  labels <- human_nt(ticks)
  gList(segmentsGrob(x0=ticks, x1=ticks, y0=0, y1=1,
                     gp=gpar(col=grey(0.3)),
                     name=paste("dna_seg_scale", i, j, "lines", sep="."),
                     default.units="native"),
        textGrob(labels$text, x=ticks, y=0.5, hjust=-0.05,
                 gp=gpar(col=grey(0.3), cex=cex),
                 name=paste("dna_seg_scale", i, j, "labels", sep="."),
                 default.units="native")
        )
}
# create gap grob
gap_grob <- function(w, m, i, j){
  segmentsGrob(x0=c(m-w/4, m-w/8),
               x1=c(m+w/8, m+w/4),
               y0=c(0.2,0.2),
               y1=c(0.8,0.8),
               gp=gpar(col=grey(0.3)),
               default.units="native",
               name=paste("gap", i, j, sep="."))
}
## yaxis grob
yaxis_grob <- function(ylim=c(0, 1), cex=0.6, n=3, i){
  at <- pretty(ylim, n=n)
  at <- at[at >= ylim[1] & at <= ylim[2]]
  coords <- yaxis_coords(at, x0=0, x1=0.5*cex)
  gList(segmentsGrob(x0=unit(coords$x0, "lines"), x1=unit(coords$x1, "lines"),
                     y0=unit(coords$y0, "native"), y1=unit(coords$y1, "native"),
                     name=paste("yaxis.segments", i, sep="."),
                     default.units="native"),
        textGrob(at, x=unit(1, "lines"), y=unit(at, "native"),
                 gp=gpar(cex=cex), just=c("left", "centre"),
                 name=paste("yaxis.labels", i, sep="."),
                 default.units="native")
        )
}
