# Various plot functions for xesq
#
# Date: 
#   Revised: February 15, 2015
#
# Author: Jiarui Ding <jiaruid@cs.ubc.ca>
#   Department of Computer Science, UBC
#   Department of Molecular Oncology, BC Cancer Agency 
#

#==============================================================================
#' Heatmap showing the connected genes' dysregulation probabilities 
#' 
#' @export
#' @param gene A character vector of gene names
#' @param posterior The xseq posteriors, output of \code{InferXseqPosterior} 
#' or \code{LearnXseqParameter}
#' @param mut A data.frame of mutations. The data.frame should have at least
#' three columns of characters: sample, hgnc_symbol, and variant_type. 
#' The variant_type column cat be either "HOMD", "HLAMP", 
#'   "MISSENSE", "NONSENSE", "FRAMESHIFT", "INFRAME", "SPLICE", "NONSTOP", 
#'   "STARTGAINED", "SYNONYMOUS", "OTHER", "FUSION", "COMPLEX".
#' @param subtype A vector representing a character of each patient, 
#' e.g., subtype
#' @param main The heatmap title
#' @param ... Other parameters passed to heatmap.2
#' 
#' @importFrom grDevices colorRampPalette
#' @importFrom stats hclust
#' 
PlotRegulationHeatmap = function(gene, posterior, mut, 
                          subtype=NULL, main="in_Cancer", ...) {
  gene = intersect(gene, names(posterior$posterior.g))
  prob.g = lapply(posterior$posterior.g[gene],
                  function(z) sapply(z, function(z1)
                    sign(z1[,3] - z1[,1]) * pmax(z1[,3],  z1[,1])))
  id = sapply(prob.g, function(z) is.matrix(z))
  gene   = gene[id]
  prob.g = prob.g[id] 
  
  id  = mut[, "hgnc_symbol"] %in% gene
  mut = mut[id, ]
  mut = ConvertMutType(mut)
  col.mut.legend  = GenerateColourSymbol("mutation")
  col.prob.legend = GenerateColourSymbol("probability")
  
  col.lab = lapply(gene, function(z) {
    id = mut[, "hgnc_symbol"] %in% z
    col.mut = col.mut.legend[mut[id, "variant_type"]]
    names(col.mut) = mut[id, "sample"]
    
    prob = posterior$posterior.f[[z]][, 2]
    prob[prob == 0] = 0.001
    col.prob = col.prob.legend[match(ceiling(prob*10)/10, 
                                     names(col.prob.legend))]
    
    names(col.prob) = names(prob)
    pat = intersect(names(col.mut), names(col.prob))
    
    tmp = data.frame(col.mut[pat], col.prob[pat])
    tmp = tmp[colnames(prob.g[[z]]), ]
    colnames(tmp) = c("mutation", "P(F)")
    
    return(tmp)
  })
  names(col.lab) = gene
  
  cols = GenerateColourSymbol(type="subtype")
  sample.all = unique(unlist(lapply(subtype, names)))
  col.subtype = lapply(subtype, function(z) {
    z = z[sample.all]
    uniq.subtype = unique(z)
    col.subtype = colorRampPalette(cols)(length(uniq.subtype))
    names(col.subtype) = uniq.subtype
    col.subtype.z = col.subtype[z]
    
    return(list(col.subtype.z, col.subtype))
  })
  
  legend = lapply(col.subtype, function(z) z[[2]])
  legend.lab = unlist(lapply(legend, names))
  col.legend = unlist(legend)
  
  col.subtype = sapply(col.subtype, function(z) z[[1]])
  rownames(col.subtype) = sample.all
  colnames(col.subtype) = names(subtype)
  
  if (length(legend.lab) > 0) {
    col.lab = lapply(col.lab, function(z)
      cbind(z, col.subtype[rownames(z), ]))
  }
  
  breaks = seq(-1, 1, by=1/32)
  col.heat = rev(colorRampPalette(RColorBrewer::brewer.pal(10, "RdBu"))(length(breaks)-1))
  tmp = sapply(gene, function(z) {
    heatmap.3(as.matrix(prob.g[[z]]),
              scale  = "none", 
              col    = col.heat, 
              breaks = breaks,
              main   =paste(z, main, sep="_"), 
              trace  = "none",
              srtCol = 45, 
              NumColSideColors = dim(col.lab[[z]])[2],
              ColSideColors = as.matrix(col.lab[[z]]),
              hclustfun = function(d) hclust(d, method="ward.D"),
              ...
    )
    
    if (length(legend.lab) > 0) {
      sample.mut = colnames(prob.g[[z]])
      
      legend.gene = lapply(names(legend), function(z) {
        id = legend[[z]] %in% col.subtype[sample.mut, z]
        legend[[z]][id] 
      })
      legend.lab = unlist(lapply(legend.gene, names))
      col.legend = unlist(legend.gene)
      
      legend("topleft", legend=legend.lab,
             fill=col.legend, border=FALSE,
             bty="n", y.intersp = 0.9, cex=0.7, ncol=ceiling(length(legend.lab)/8),
             x.intersp=0.5, title="All")
    } 
    
    mut.legend.col = unique(as.character(col.lab[[z]][,1]))
    mut.legend = sapply(mut.legend.col, 
                        function(z) names(which(col.mut.legend == z)))
    legend("topleft", legend=tolower(mut.legend), 
           fill=mut.legend.col, border=FALSE, 
           bty="n", y.intersp = 0.9, cex=0.7, ncol=3, 
           x.intersp=0.5, title="Mutation type")
    
    
    names(col.prob.legend) = c("0.0 - 0.1", "0.1 - 0.2", "0.2 - 0.3", 
                               "0.3 - 0.4", "0.4 - 0.5", "0.5 - 0.6", 
                               "0.6 - 0.7", "0.7 - 0.8", "0.8 - 0.9", 
                               "0.9 - 1.0")
    legend("topright", legend=names(col.prob.legend), 
           fill=col.prob.legend, border=FALSE, 
           bty="n", y.intersp = 0.9, cex=0.7, ncol=5, 
           x.intersp=0.5, title="P(F)")
  })
  
  return(list(prob.g=prob.g, legend=legend))
}


#==============================================================================
#' @importFrom graphics pie plot rect
#' @importFrom grDevices colorRampPalette rgb
GenerateColourSymbol = function(type="probability", show.plot=FALSE)  {
  # Probability in the range of 0--2.5 are given the same colour
  if (type == "probability") {
    col = colorRampPalette(RColorBrewer::brewer.pal(9, "YlOrRd"))(10)
    name = as.character(1:10/10)
    names(col) = name
  } else if (type == "mutation") {
    variant.type = c("HLAMP", "HOMD", "MISSENSE", "NONSENSE",
                     "FRAMESHIFT", "INFRAME", "SPLICE", "NONSTOP", 
                     "SYNONYMOUS", "STARTGAINED", "OTHER", "FUSION", 
                     "COMPLEX")
    col = c(RColorBrewer::brewer.pal(9, "Set1"))
    col = c(col, RColorBrewer::brewer.pal(12, "Set3"))
    col = col[-c(11, 13, 14)]
    names(col) = variant.type
    col = col[c("HOMD", "HLAMP", "MISSENSE", "NONSENSE", "FRAMESHIFT", 
                "INFRAME", "SPLICE", "NONSTOP", "STARTGAINED", "SYNONYMOUS", 
                "OTHER", "FUSION", "COMPLEX")]
    nomut = rgb(0.95, 0.95, 0.95, 0.8)
    names(nomut) = "NOMUTATION"
    col = c(col, nomut)
  } else if (type == "subtype") {
    col = RColorBrewer::brewer.pal(8, "Accent")
    col = c(col, RColorBrewer::brewer.pal(9, "Pastel1"))
    #col = col[-c(17)]
  } else if (type == "symbol") {
    mut.type = c("MISSENSE", "NONSENSE", "FRAMESHIFT", 
                 "INFRAME", "SPLICE", "NONSTOP", "STARTGAINED", "SYNONYMOUS", 
                 "OTHER", "FUSION", "COMPLEX", "HOMD", "HLAMP")
    
    #pch = c(21, 25, 25, 21, 25, 22, 23, 21, 21, 22, 23, 25, 24)
    pch = c(21, 25, 23, 24, 22, 21, 21, 24, 21, 21, 21, 21, 21)
    names(pch) = mut.type
    
    return(pch)
  }
  
  if (show.plot == TRUE) {
    #par(mfrow = c(2, 1))
    
    pie(rep(1,length(col)), col=col)
    plot(0, type="n", ylab="", xlab="",
         axes=FALSE, ylim=c(1,0), xlim=c(1,length(col)))
    
    for (i in 1:length(col)) {
      rect(i-0.5, -0.5, i+0.5, 0.5, border="black", col=col[i])
    }    
  }
  
  return(col)
}


# ==============================================================================
#' @importFrom graphics boxplot axis mtext points legend
boxjitter = function(x, y, mut=NULL, mut.col=NA, mut.pch=NA, 
                     xlab=NA, ylab=NA, border="gray", horizontal=TRUE, 
                     factor=0.5, ...) {
  # R function that plots boxplot and adds data points (jittered vertically)
  #
  # Date: 
  #   updated date: 2013-05-23
  #   Revised: February 15, 2015
  #
  # Author: Jiarui Ding <jiaruid@cs.ubc.ca>
  #   Department of Computer Science, UBC
  #   Department of Molecular Oncology, BC Cancer Agency 
  #
  
  box = boxplot(x~y, horizontal=horizontal, border=border, 
                outline=FALSE, axes=FALSE, xlab=NA, ylab=NA, ...)
  
  ## Change the axis marks
  box()
  axis(side=1, tck=-0.015, mgp=c(3,0.5,0))
  axis(side=2, tck=-0.015, mgp=c(3,0.5,0), 
       las=1, labels=box$names, at=1:length(box$names))
  
  mtext(side=1, xlab, line=1.5)
  mtext(side=2, ylab, line=1.5)
  
  ## Change the yaxis
  ycoord = y
  i = 1
  for (name in box$names) {
    id = (y == name)
    ycoord[id] = i
    i = i + 1
  } 
  ycoord = jitter(as.numeric(ycoord), factor=factor) 
  
  ## Change the mutations
  col = rep("dodgerblue", length(x))
  pch = rep(4, length(x))
  cex = rep(1, length(x))
  lwd = rep(1, length(x))
  
  names(col) = names(x)
  names(pch) = names(x)
  names(cex) = names(x)
  names(lwd) = names(x)
  
  sample = unique(mut[, "sample"])
  variant = rep("MISSENSE", length(sample))
  names(variant) = sample
  for (sample.i in sample) {
    id = mut[, "sample"] %in% sample.i
    if(sum(id) > 1) {
      variant[sample.i] = "COMPLEX"
    } else {
      variant[sample.i] = mut[id, "variant_type"]
    }
  }
  
  pch[sample] = mut.pch[variant]
  cex[sample] = 1.0
  lwd[sample] = 1.0
  col[sample] = mut.col[variant]
  
  id = which(names(x) %in% sample)
  id.nomut = setdiff(1:length(x), id)
  points(x[id.nomut], ycoord[id.nomut], col=AddTrans(col[id.nomut], 0.2), 
         pch=pch[id.nomut], cex=cex[id.nomut], lwd=lwd[id.nomut])
  points(x[id], ycoord[id], col="ivory4", pch=pch[id], cex=cex[id], 
         lwd=lwd[id], bg=col[id])
  
  id = which(!duplicated(variant[sample]))
  mut.type = 1:13
  names(mut.type) = c("SYNONYMOUS", "MISSENSE", "INFRAME",  "COMPLEX", 
                      "FRAMESHIFT", "NONSENSE", "SPLICE", "NONSTOP", 
                      "STARTGAINED", "FUSION", "OTHER", "HOMD", "HLAMP")
  id.c = order(mut.type[variant[id]])
  id = id[id.c]
  
  ncol = 1
  if (length(id) > 3) {
    ncol = 2
  }
  if (sum(id) > 0) {
    legend("topleft", legend=tolower(variant[id]), col="ivory4", bty="n", 
           pch=pch[sample][id], pt.bg=col[sample][id], ncol=ncol)
  }
  
  return()
}


# ==============================================================================
#' @importFrom grDevices col2rgb 
AddTrans = function(col, trans) {
  # Add transcripancy to colours. Works with either color and trans a vector 
  # of equal length, or one of the two of length 1.
  
  # col  :
  #      colour vector
  # trans:
  #      0 being fully visable and 1 being fully transpant
  # 
  if (any(trans > 1) | any(trans < 0)) {
    stop("Error: trans should in the range of [0, 1]") 
  }
  if (length(col) != length(trans) & !any(c(length(col), length(trans))==1)) 
    stop("Error: vector lengths doesn't match")
  if (length(col)==1 & length(trans)>1)
    col = rep(col, length(trans))
  if (length(trans)==1 & length(col)>1) 
    trans = rep(trans, length(col))
  
  num2hex = function(x) {
    hex = unlist(strsplit("0123456789ABCDEF", split=""))
    return (paste(hex[(x-x%%16)/16 + 1], hex[x%%16 + 1], sep=""))
  }
  rgb = rbind(col2rgb(col), round(255 - trans*255))
  res = paste("#", 
              apply(apply(rgb, 2, num2hex), 2, paste, collapse=""), 
              sep="")
  
  return(res)
}


#==============================================================================
#' @importFrom graphics plot axis mtext 
NeatPlot = function(..., xlab, ylab, xaxt="s", yaxt="s", 
                    xtck.label=TRUE, ytck.label=TRUE, cex.axis=0.8) {
  
  plot(xaxt="n", yaxt="n", xlab=NA, ylab=NA, ...) 
  
  if (xaxt == "s") {
    axis(side=1, tck=-0.015, mgp=c(3,0.5,0), 
         labels=xtck.label, cex.axis=cex.axis)
  }
  if (yaxt == "s") {
    axis(side=2, tck=-0.015, mgp=c(3,0.5,0), 
         labels=ytck.label, cex.axis=cex.axis)
  }
  
  if (!missing(xlab)) {
    mtext(side=1, text=xlab, line=1.5, cex=cex.axis)
  }
  if (!missing(ylab)) {
    mtext(side=2, text=ylab, line=1.5, cex=cex.axis)
  }
}


