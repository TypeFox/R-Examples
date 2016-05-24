#' @title Plot PATHMOX and TECHMOX trees
#' 
#' @description
#' The function \code{plot.treemox} allows to display binary trees of PATHMOX
#' and TECHMOX analyses. If \code{shadow.size=0}, no shadows are drawn.
#' 
#' @param x An object of class \code{"treemox"} returned by
#' \code{\link{pathmox}} or \code{\link{techmox}}.
#' @param root.col Fill color of root node.
#' @param root.bor Border color of root node.
#' @param root.txt Text color of root node.
#' @param root.cex magnification to be used for text in root node.
#' @param root.lwd Line width of border in the root node.
#' @param root.shadow Color of shadow of root node.
#' @param node.col Fill color of child nodes.
#' @param node.bor Border color of child nodes.
#' @param node.txt Text color of child nodes.
#' @param node.cex magnification to be used for text in child nodes.
#' @param node.lwd Line width of border in child nodes.
#' @param node.shadow Color of shadow of child nodes.
#' @param leaf.col Fill color of leaf nodes.
#' @param leaf.bor Border color of leaf nodes.
#' @param leaf.txt Text color of leaf nodes.
#' @param leaf.cex magnification to be used for text in leaf nodes.
#' @param leaf.lwd Line width of border in leaf nodes.
#' @param leaf.shadow Color of shadow of leaf nodes.
#' @param shadow.size Relative size of shadows.
#' @param arr.lwd Line width of the tree branches.
#' @param lcol color of lines
#' @param arr.col color of arrows
#' @param seg.cex A numerical value indicating the magnification to be used for
#' plotting text.
#' @param seg.col The color to be used for the labels of the segmentation
#' variables.
#' @param cat.cex magnification to be used for the categories
#' @param cat.col The color to be used for the labels of the categories
#' @param show.pval Logical value indicating whether the p-values should be
#' plotted.
#' @param pval.col The color to be used for the labels of the p-values.
#' @param main A main title for the plot.
#' @param cex.main The magnification to be used for the main title.
#' @param col.main Color to be used for the main title
#' @param \dots Further arguments are ignored.
#' @method plot treemox
#' @S3method plot treemox
#' @examples
#'
#'  \dontrun{
#'  ## example of PLS-PM in customer satisfaction analysis
#'  ## model with seven LVs and reflective indicators
#'  data(csimobile)
#'  
#'  # select manifest variables
#'  data_mobile = csimobile[,8:33]
#'  
#'  # define path matrix (inner model)
#'  IMAG = c(0, 0, 0, 0, 0, 0, 0)
#'  EXPE = c(1, 0, 0, 0, 0, 0, 0)
#'  QUAL = c(0, 1, 0, 0, 0, 0, 0)
#'  VAL = c(0, 1, 1, 0, 0, 0, 0)
#'  SAT = c(1, 1, 1, 1, 0, 0, 0)
#'  COM = c(0, 0, 0, 0, 1, 0, 0)
#'  LOY = c(1, 0, 0, 0, 1, 1, 0)
#'  mob_path = rbind(IMAG, EXPE, QUAL, VAL, SAT, COM, LOY)
#'  
#'  # blocks of indicators (outer model)
#'  mob_outer = list(1:5, 6:9, 10:15, 16:18, 19:21, 22:24, 25:26)
#'  mob_modes = rep("A", 7)
#'  
#'  # apply plspm
#'  mob_pls = plspm(data_mobile, mob_path, mob_blocks, modes = mob_modes, 
#'                  scheme = "factor", scaled = FALSE)
#'
#'  # re-ordering those segmentation variables with ordinal scale 
#'  # (Age and Education)
#'  csimobile$Education = factor(csimobile$Education, 
#'      levels=c("basic","highschool","university"),
#'      ordered=TRUE)
#'  
#'  # select the segmentation variables
#'  seg_vars = csimobile[,1:7]
#'  
#'  # Pathmox Analysis
#'  mob_pathmox = pathmox(mob_pls, seg_vars, signif=.10, size=.10, deep=2)
#'  
#'  # plot pathmox tree
#'  plot(mob_pathmox, root.col="lightblue", node.col="turquoise", leaf.col="skyblue3", 
#'       shadow.size=0, seg.col="blue2", pval.col="magenta")
#'  }
#'
plot.treemox <-
function(x, root.col = "#eeeeee", root.bor = "#cccccc", root.txt = "#757575", 
         root.cex = 0.8, root.lwd = 3, root.shadow = "gray40",
         node.col = "#feb769", node.bor = "#FE9929", node.txt = "#555555", 
         node.cex = 0.7, node.lwd = 3, node.shadow = "gray30",
         leaf.col = "#93c4e5", leaf.bor = "#5a99c5", leaf.txt ="#555555", 
         leaf.cex = 0.7, leaf.lwd = 3, leaf.shadow = "gray30",
         shadow.size = 0, arr.lwd=3, lcol = "#ddddddbb",  arr.col = "gray95",
         seg.cex = 0.7, seg.col = "#2cb0a7", 
         cat.cex = 0.8, cat.col = "#555555",
         show.pval = TRUE, pval.col = "#2cb0a7", 
         main=NULL, cex.main=1, col.main="gray50", ...)
{
  # =======================================================
  # inputs setting
  # =======================================================  
  MOX <- x$MOX
  last <- nrow(MOX)
  last.level <- MOX$Depth[last]
  num.levels <- rep(1,last.level+1)
  for (i in 1:last.level)
    num.levels[i+1] <- 2^i

  # =======================================================
  # plot parameters
  # =======================================================  
  #dev.new()
  op = par(mar=c(.4,.4,1,1.5))
  openplotmat()
  # estimates coordinates of elements, arranged on a plot
  elpos <- coordinates(num.levels)
  fromto <- cbind(MOX[-1,2], MOX[-1,1])
  nr <- nrow(fromto)
  arrpos <- matrix(ncol=2, nrow=nr)
  # drawing the binary tree structure
  for (i in 1:nr) {
    arrpos[i,] = straightarrow(to=elpos[fromto[i,2],], from=elpos[fromto[i,1],],
                      lwd=arr.lwd, lcol=lcol, arr.col = arr.col, arr.pos=0.6, arr.length=0)
  }
  # drawing the root node (i.e. parent node)
  textellipse(elpos[1,], 0.045, 0.045, lab=c("Root",MOX[1,6]), box.col=root.col, 
              col=root.txt, lcol=root.bor, lwd=root.lwd,
              shadow.size=shadow.size, shadow.col=root.shadow, cex=root.cex)
  ## drawing the child nodes
  for (i in 2:last)
  {
    posi <- MOX$Node[i]
    nodlab <- c(paste("Node",posi), MOX$Size[i])
    if (MOX$Type[i]=="node") {
      textellipse(elpos[posi,], 0.05, 0.03, lab=nodlab, box.col=node.col, 
                  col=node.txt, lcol=node.bor, lwd=node.lwd,
                  shadow.col=node.shadow, shadow.size=shadow.size, cex=node.cex)
    } else { # "leaf"
      textrect(elpos[posi,], 0.045, 0.025, lab=nodlab, box.col=leaf.col, 
               col=leaf.txt, lcol=leaf.bor, lwd=leaf.lwd,
               shadow.col=leaf.shadow, shadow.size=shadow.size, cex=leaf.cex)
    }
  }
  # adding segmentation variables and p.values
  aux <- 1
  for (i in seq(1,nr,by=2))
  {   
    if (i==1) k<-1 else k<-1.15     
    x1 <- (arrpos[i,1] + arrpos[i+1,1]) / 2
    text(x1, k*arrpos[i,2], MOX$Variable[i+1], cex=seg.cex, col=seg.col)
    if (show.pval)
    {
      if (x$model$mox=="pathmox")
        text(x1, k*arrpos[i,2], paste("p.val=",round(x$FT$p.val[aux],4),sep=""), 
             cex=seg.cex, col=pval.col, pos=1) 
      if (x$model$mox=="techmox")
      {
        pvals <- unlist(lapply(x$FT, function(x) exp(mean(log(x[,4]))) ))
        text(x1, k*arrpos[i,2], paste("pv.gm=",round(pvals[aux],4),sep=""), 
             cex=seg.cex, col=pval.col, pos=1) 
      }
    }
    aux <- aux + 1
  }     
  ## adding segmentation categories
  for (i in 1:nr)
  {
    posi <- MOX$Node[i+1]        
    seg.cat <- as.character(MOX$Category[i+1])
    seg.cat <- unlist(strsplit(seg.cat,"/"))
    if (posi%%2 == 0) {
      for (h in 1:length(seg.cat))
        text(arrpos[i,1]-0.03, arrpos[i,2]+h/50, seg.cat[h], cex=cat.cex, col=cat.col)
    } else {
      for (h in 1:length(seg.cat))
        text(arrpos[i,1]+0.03, arrpos[i,2]+h/50, seg.cat[h], cex=cat.cex, col=cat.col)
    }
  }
  ## adding main title
  if (is.null(main)) {
    if (x$model$mox=="pathmox")
      text(.5, .95, c("PATHMOX Tree"), cex=cex.main, col=col.main)
    if (x$model$mox=="techmox")
      text(.5, .95, c("TECHMOX Tree"), cex=cex.main, col=col.main)
  } else {
    text(.5, .95, main, cex=cex.main, col=col.main)
  }
  # reset par
  par(op)
}
