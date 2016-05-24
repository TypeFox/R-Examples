# Copyright (C) 2015
# Benoit Gautier, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
# Ignacio Gonzalez, Genopole Toulouse Midi-Pyrenees, France
# Francois Bartolo, Institut National des Sciences Appliquees et Institut de Mathematiques, Universite de Toulouse et CNRS (UMR 5219), France
# Florian Rohart, Australian Institute for Bioengineering and Nanotechnology, University of Queensland, Brisbane, QLD.
# Kim-Anh Le Cao, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD

# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.


#----------------------------------------------------------------------------------------------------------#
#-- Includes plotIndiv for PLS, sPLS, PLS-DA, SPLS-DA, rCC, PCA, sPCA, IPCA, sIPCA, rGCCA, sGCCA, sGCCDA --#
#----------------------------------------------------------------------------------------------------------#

plotIndiv <-
  function(object,
           comp = NULL,
           rep.space = NULL,
           blocks = NULL, # to choose which block data to plot, when using GCCA module
           ind.names = TRUE,
           group,  # factor indicating the group membership for each sample, useful for ellipse plots. Coded as default for the -da methods, but needs to be input for the unsupervised methods (PCA, IPCA...)           
           col.per.group,
           style="ggplot2", # can choose between graphics,3d, lattice or ggplot2
           plot.ellipse = FALSE,
           ellipse.level = 0.95,
           plot.centroid=FALSE,
           plot.star=FALSE,
           main=NULL,
           add.legend=FALSE,
           X.label = NULL,
           Y.label = NULL,
           Z.label = NULL,
           abline.line = FALSE,
           xlim = NULL,
           ylim = NULL,
           col,
           cex,
           pch,
           alpha=0.2,
           axes.box = "box",
           ...
  )
{
    
    class.object = class(object)
    object.pls=c("pls","spls","splsda","plsda","mlspls","mlsplsda","rcc")
    object.pca=c("ipca","sipca","pca","spca","prcomp")
    object.blocks=c("sgcca","rgcca", "sgccda")
    
    ### Start: Validation of arguments
    ncomp = object$ncomp
    
    
    #-- rep.space
    if(is.null(rep.space) && class.object[1] %in% c(object.blocks,object.pca,"splsda","plsda","mlsplsda"))
      rep.space="X-variate"
    else if(is.null(rep.space))
      rep.space="multi"
    rep.space = match.arg(rep.space, c("XY-variate", "X-variate", "Y-variate","multi"))
    
    
    
    if (class.object[1] %in% object.blocks) {
      
      if (is.null(blocks)){
        blocks = object$names$blocks
        
        if (class.object[1] == "sgccda")
          blocks = blocks[-object$indY]
      } else if (is.numeric(blocks) & min(blocks) > 0 &  max(blocks) <= length(object$names$blocks)) {
        blocks = object$names$blocks[blocks]
      } else if (is.character(blocks)) {
        if (!any(blocks %in% object$names$blocks))
          stop("One element of 'blocks' does not match with the names of the blocks")
      } else {
        stop("Incorrect value for 'blocks", call. = FALSE)
      }
      object$variates = object$variates[names(object$variates) %in% blocks]
      
      if (any(object$ncomp[blocks] == 1)) {
        stop(paste("The number of components for one selected block '", paste(blocks, collapse = " - "),"' is 1. The number of components must be superior or equal to 2."), call. = FALSE)
      }
      
      if (is.null(object$indY) & rep.space %in% c("Y-variate", "XY-variate")) {
        stop("For an object of class 'blocks', 'rep.space' must be 'X-variate", call. = FALSE)
      }
      ncomp = object$ncomp[blocks]
    }
    
    
    
    #-- style
    if (!style %in% c("ggplot2", "lattice", "graphics","3d"))
      stop("'style' must be one of 'ggplot2', '3d' , lattice' or 'graphics'.", call. = FALSE)
    
    #-- axes.box
    choices = c("box", "bbox", "both")
    axes.box = choices[pmatch(axes.box, choices)] 
    
    if (is.na(axes.box)) 
      stop("'axes.box' should be one of 'box', 'bbox' or 'both'.", 
           call. = FALSE)
    
    
    #-- ellipse.level
    if ((ellipse.level > 1) | (ellipse.level < 0))
      stop("The value taken by 'ellipse.level' must be between 0 and 1")
    
    #-- add.legend
    if (length(add.legend) != 1 || !is.logical(add.legend))
      stop("'add.legend' must be a logical value.", call. = FALSE)
    
    
    #-- alpha correlation
    if (!is.numeric(alpha) | (alpha > 1) | (alpha < 0))
      stop("The value taken by 'alpha' must be between 0 and 1", call. = FALSE)
    
    
    #-- comp
    if(is.null(comp))
    {if (style=="3d")
      comp=c(1:3)
     else
       comp=c(1:2)}
    if (length(comp) != 2 && !(style=="3d"))
      stop("'comp' must be a numeric vector of length 2.", call. = FALSE) 
    else if(length(comp) != 3 && (style=="3d"))
      stop("'comp' must be a numeric vector of length 3.", call. = FALSE)
    
    if (!is.numeric(comp))
      stop("Invalid vector for 'comp'.")
    
    if (any(ncomp < max(comp)))
      stop("Each element of 'comp' must be smaller or equal than ", max(object$ncomp), ".", call. = FALSE)
    
    comp1 = round(comp[1]); comp2 = round(comp[2])
    if (style=="3d") comp3 = round(comp[3])
    
    #-- Specific to pls object (choice between X, Y and XY)
    if (class.object[1] %in% object.pls){
      if (rep.space == "multi")
      {blocks=c("X","Y");object$variates = object$variates[names(object$variates) %in% blocks]}
      if (rep.space == "X-variate")
      {object$variates = object$variates["X"]; blocks = "X"}
      if (rep.space == "Y-variate")
      {object$variates = object$variates["Y"]; blocks = "Y"}
      if (rep.space == "XY-variate"){
        object$variates$XYvariates = (object$variates$X + object$variates$Y)/2
        object$variates = object$variates["XYvariates"]; blocks = "XY combined"
      }
    }
    
    if (class.object[1] %in% object.pca)
      blocks = "X"
    
    
    #-- xlim,ylim
    lim.X = FALSE
    if (!is.null(xlim)){
      if (class.object[1] %in% object.blocks ||!class.object[1] %in% c(object.pca,"splsda","plsda","mlsplsda")){
        if (!is.list(xlim) || length(xlim) != length(blocks) || length(unlist(xlim)) != 2 * length(blocks))
          stop("'xlim' must be a list of ",length(blocks)," vectors of length 2.",call. = FALSE)
      } else {
        if (!is.numeric(xlim) || length(xlim) != 2)
          stop("'xlim' must be a vector of length 2.",call. = FALSE)
        xlim=list(xlim)
        
      }
      lim.X = TRUE
    }
    
    lim.Y = FALSE
    if (!is.null(ylim)){
      if (class.object[1] %in% object.blocks ||!class.object[1] %in% c(object.pca,"splsda","plsda","mlsplsda")){
        if (!is.list(ylim) || length(ylim) != length(blocks) || length(unlist(ylim)) != 2 * length(blocks))
          stop("'ylim' must be a list of ",length(blocks)," vectors of length 2.",call. = FALSE)
      } else {
        if (!is.numeric(ylim) || length(ylim) != 2)
          stop("'ylim' must be a vector  of length 2.",call. = FALSE)
        ylim=list(ylim)
      }
      lim.Y = TRUE
    }
    
    #-- Start: Retrieve variates from object
    x = y = z=list()
    if (class.object[1] %in%  c(object.pls, object.blocks)) {
      
      if (is.logical(ind.names) & isTRUE(ind.names)) {
        ind.names = object$names$indiv
      }
      if (length(ind.names) > 1) {
        if (length(ind.names) != length(object$names$indiv))
          stop("'ind.names' must be a character vector of length ", nrow(object$X), " or a boolean atomic vector.")
      }
      
      x = lapply(object$variates, function(x){x[, comp1, drop = FALSE]})
      y = lapply(object$variates, function(x){x[, comp2, drop = FALSE]})
      if(style=="3d") z = lapply(object$variates, function(x){x[, comp3, drop = FALSE]})
      
      if (is.null(X.label)) {
        if (rep.space == "multi") {X.label = paste("variate", comp1)}
        if (rep.space == "X-variate") {X.label = paste("X-variate", comp1)}
        if (rep.space == "Y-variate") {X.label = paste("Y-variate", comp1)}
        if (rep.space == "XY-variate") {X.label = paste("XY-variate", comp1)}
      }
      if (is.null(Y.label)) {
        if (rep.space == "multi") {Y.label = paste("variate", comp2)}
        if (rep.space == "X-variate") {Y.label = paste("X-variate", comp2)}
        if (rep.space == "Y-variate") {Y.label = paste("Y-variate", comp2)}
        if (rep.space == "XY-variate") {Y.label = paste("XY-variate", comp2)}
      }
      if (is.null(Z.label)&&style=="3d") {
        if (rep.space == "multi") {Z.label = paste("variate", comp3)}
        if (rep.space == "X-variate") {Z.label = paste("X-variate", comp3)}
        if (rep.space == "Y-variate") {Z.label = paste("Y-variate", comp3)}
        if (rep.space == "XY-variate") {Z.label = paste("XY-variate", comp3)}
      }
      
    } else if (class.object[1] %in%  object.pca) {
      
      if (is.logical(ind.names)) {
        if (isTRUE(ind.names))
          ind.names = rownames(object$x)
      }
      if (length(ind.names) > 1) {
        if (length(ind.names) != nrow(object$x))
          stop("'ind.names' must be a character vector of length ", nrow(object$x), " or a boolean atomic vector.")
      }
      if (rep.space != "X-variate")
        stop("'rep.space' must be equal to 'X-variate' for an object belonging to 'pca' or 'ipca'")
      
      x[[1]] = object$x[, comp[1]]
      y[[1]] = object$x[, comp[2]]
      if(style=="3d") z[[1]] = object$x[, comp[3]]
      
      
      #-- Variance explained
      if (class.object[1] %in% "pca"){
        if (style == "3d")
        {
          inf=c((object$sdev)[comp1]/sum(object$sdev),(object$sdev)[comp2]/sum(object$sdev),(object$sdev)[comp3]/sum(object$sdev))
          inf = round(inf,2)
        } else {
          inf = c((object$sdev)[comp1]/sum(object$sdev),(object$sdev)[comp2]/sum(object$sdev))
          inf = round(inf,2)}
        
      }   
      if (is.null(X.label)) {
        X.label = paste("PC", comp1,sep='')
        if (class.object[1] %in% "pca") {percentage=paste(inf[1]*100,"% expl. var")
                                         X.label = paste(X.label,percentage,sep=": ")}
      }
      if (is.null(Y.label)){
        Y.label = paste("PC", comp2,sep='')
        if (class.object[1] %in% "pca") {percentage=paste(inf[2]*100,"% expl. var")
                                         
                                         Y.label = paste(Y.label,percentage,sep=": ")}
      }
      if (is.null(Z.label)&&style=="3d"){ 
        Z.label = paste("PC", comp3,sep='')
        if (class.object[1] %in% "pca") {percentage=paste(inf[3]*100,"% expl. var")
                                         Z.label = paste(Z.label,percentage,sep=": ")}
      }
    }
    
    
    #-- End: Retrieve variates from object
    
    #-- ind.names
    display.names = FALSE
    if (length(ind.names) == length(x[[1]])){
      display.names = TRUE
    }
    
    
    
    
    #-- Define group
    missing.group = FALSE
    if (missing(group) & any(class.object %in% c("plsda","splsda","mlsplsda"))){
      group = factor(map(object$ind.mat), labels = object$names$Y)
    } else if (missing(group) & any(class.object %in% c("sgccda"))){
      group = factor(map(object$ind.mat), labels = object$names$colnames$Y)
    } else if (!missing(group)) {
      missing.group = TRUE
      if (!is.factor(group)){
        group = as.factor(group)
      }
      object$ind.mat = unmap(group)
      
      if (length(group) != length(x[[1]]))
        stop("Length of 'group' should be of length ", length(x[[1]]), ", the sample size of your data")
    } else {
      if(plot.star || plot.centroid || plot.ellipse)
        warning('plot.star , plot.ellipse and plot.centroid work only if !group==NULL')
      plot.star=plot.centroid=plot.ellipse=FALSE
      group = factor(rep("No group", length(x[[1]])))
      object$ind.mat = unmap(group)
    }
    
    #-- col.per.group argument
    if (missing(col.per.group)){
      if (nlevels(group) < 10) {
        #only 10 colors in color.mixo
        col.per.group = color.mixo(1:nlevels(group))
      } else {
        #use color.jet
        col.per.group = color.jet(nlevels(group))
      }
    } else {
      if (length(col.per.group) == 1) {
        col.per.group = rep(col.per.group, nlevels(group))
      } else if (length(col.per.group) != length(x[[1]]) & length(col.per.group) != nlevels(group)) {
        stop("Length of 'col.per.group' should be either of length 1 or of length ", nlevels(group), " (the number of groups) or of length ", length(x[[1]]), " (the sample size or your data).
             Alternatively, use the argument 'col' to give one color per sample")
      }
      missing.group = TRUE
      }
    
    levels.color = vector(, length(x[[1]]))
    if (length(col.per.group) != length(x[[1]])) {
      for (i in 1 : nlevels(group)){
        levels.color[group == levels(group)[i]] = col.per.group[i]
      }
    } else {
      levels.color = col.per.group
    }
    
    #-- col argument   
    missing.col = FALSE
    if (!missing(col)){
      if (length(col) > length(x[[1]]))
        stop("Length of 'col' should be of length inferior or equal to ", length(x[[1]]),".")
      
      col = factor(rep(col, ceiling(length(x[[1]])/length(col)))[1 : length(x[[1]])])
      if (!missing.group) {
        group = col
        levels.color = col
        col.per.group = levels(col)
        object$ind.mat = unmap(group)
      }
      missing.col = TRUE
    } else {
      col = levels.color
    }
    
    #-- cex argument
    if (missing(cex)){
      if (style == "ggplot2"){
        cex = rep(5, length(x[[1]]))
        cex = cex[as.factor(group)]
      } else {
        cex = rep(1, length(x[[1]]))
        cex = cex[as.factor(group)]
      }
    } else {
      if (length(cex) == 1){
        cex = rep(cex, length(x[[1]]))
        cex = cex[as.factor(group)]
      } else if (length(cex) > length(x[[1]])) {
        stop("Length of 'cex' should be of length inferior or equal to ", length(x[[1]]),".")
      } else if (length(cex) == length(unique(group))){
        cex = cex[as.factor(group)]
      }else {
        cex = rep(cex, ceiling(length(x[[1]])/length(cex)))[1 : length(x[[1]])]
      }
    }
    
    #-- pch argument
    if (missing(pch)){
      if (missing.col){
        if(style=="3d"){
          pch = unlist(lapply(1 : length(length(levels(col))), function(x){rep(c("sphere", "tetra", "cube", "octa", "icosa", "dodeca")[x], length(col==x))}))}
        else
          pch = as.numeric(col)
      } else {
        if(style=="3d"){
          pch = unlist(lapply(1 : length(length(levels(group))), function(x){rep(c("sphere", "tetra", "cube", "octa", "icosa", "dodeca")[x], length(group==x))}))}
        else
          pch = as.numeric(group)
      }
    } else {
      if (style=="3d"){
        if (!all(unlist(pch) %in% c("sphere", "tetra", "cube", "octa", "icosa", "dodeca")))
          stop("pch' must be a simple character or character vector from {'sphere', 'tetra', 'cube', 'octa', 'icosa', 'dodeca'}.", 
               call. = FALSE)
      }
      if (length(pch) == 1){
        pch = rep(pch, length(x[[1]]))
      } else if (length(pch) > length(x[[1]])){
        stop("Length of 'pch' should be of length inferior or equal to ", length(group),".")
      } else if (length(pch) == length(unique(group))){
        pch = pch[as.factor(group)]
      } else {
        pch = rep(pch, ceiling(length(x[[1]])/length(pch)))[1 : length(x[[1]])]
      }
      display.names = FALSE
    }
    
    
    
    if (plot.ellipse) {
      #-- Start: Computation ellipse
      min.ellipse = max.ellipse = xlim.min = xlim.max = ylim.min = ylim.max = list()
      ind.gp = matrice = cdg = variance = list()
      ind.gp = lapply(1 : ncol(object$ind.mat), function(x){which(object$ind.mat[, x]==1)})
      matrice = lapply(1 : length(x), function(z1) {lapply(ind.gp, function(z2){matrix(c(x[[z1]][z2], y[[z1]][z2]), ncol = 2)})})
      cdg = lapply(1 : length(x), function(z){ lapply(matrice[[z]], colMeans)})
      variance = lapply(1 : length(x), function(z){lapply(matrice[[z]], var)})
      coord.ellipse = lapply(1 : length(x), function(z1){ lapply(1 : ncol(object$ind.mat), function(z2){ellipse(variance[[z1]][[z2]],
                                                                                                                centre = cdg[[z1]][[z2]],
                                                                                                                level = ellipse.level)})})
      max.ellipse = lapply(1 : length(x), function(z1) {sapply(coord.ellipse[[z1]], function(z2){apply(z2, 2, max)})})
      min.ellipse = lapply(1 : length(x), function(z1) {sapply(coord.ellipse[[z1]], function(z2){apply(z2, 2, min)})})
      #-- End: Computation ellipse
      if (is.null(xlim))
        xlim = lapply(1 : length(x), function(z) {c(min(x[[z]], min.ellipse[[z]][1, ]), max(x[[z]], max.ellipse[[z]][1, ]))})
      if (is.null(ylim))
        ylim = lapply(1 : length(x), function(z) {c(min(y[[z]], min.ellipse[[z]][2, ]), max(y[[z]], max.ellipse[[z]][2, ]))})
    } else {
      if (is.null(xlim))
        xlim = lapply(1 : length(x), function(z) {c(min(x[[z]]), max(x[[z]]))})
      if (is.null(ylim))
        ylim = lapply(1 : length(x), function(z) {c(min(y[[z]]), max(y[[z]]))})
    }
    
    
    #-- Start: data set
    df = list()
    if(style=="3d")
    {for (i in 1 : length(x)) {
      df[[i]] = data.frame(x = x[[i]], y = y[[i]],z=z[[i]],group = group)
    }}
    else
    {
      for (i in 1 : length(x)) {
        df[[i]] = data.frame(x = x[[i]], y = y[[i]], group = group)
      }
    }
    
    if(class.object[1] %in% c("ipca","sipca","pca","spca","prcomp", "splsda","plsda","mlsplsda"))
    { 
      if(is.null(main)){
        df = data.frame(do.call(rbind, df), "Block" = "PlotIndiv")}
      else
      {
        df = data.frame(do.call(rbind, df), "Block" = main) 
      }}
    
    else
      df = data.frame(do.call(rbind, df), "Block" = paste0("Block: ", unlist(lapply(1 : length(df), function(z){rep(blocks[z], nrow(df[[z]]))}))))
    
    
    
    
    if(style=="3d")
      names(df)[1:3] = c("x", "y","z")
    else
      names(df)[1:2] = c("x", "y")
    
    if (display.names)
      df$names = rep(ind.names, length(x))
    
    df$pch = pch; df$cex = cex; df$col.per.group = levels.color[group]; df$col = as.character(col)
    
    if (plot.centroid == TRUE || plot.star==TRUE){
      df=cbind(df,rep(0,nrow(df)))
      n=ncol(df)
      df=cbind(df,rep(0,nrow(df)))
      for (i in 1:nlevels(group)){
        if(length(x)>1){
          for (k in 1 : length(x)){
            x0=mean(df[df$group == levels(group)[i] & df$Block %in% paste0("Block: ", blocks[k]), "x"])
            y0=mean(df[df$group == levels(group)[i] & df$Block %in% paste0("Block: ", blocks[k]) , "y"])
            df[df$group == levels(group)[i] & df$Block %in% paste0("Block: ", blocks[k]), n]=x0
            df[df$group == levels(group)[i] & df$Block %in% paste0("Block: ", blocks[k]), n+1]=y0
            names(df)[c(ncol(df)-1,ncol(df))]=c("x0","y0")
          }}
        else
        {
          x0=mean(df[df$group == levels(group)[i] , "x"])
          y0=mean(df[df$group == levels(group)[i]  , "y"])
          df[df$group == levels(group)[i] , n]=x0
          df[df$group == levels(group)[i] , n+1]=y0
          names(df)[c(ncol(df)-1,ncol(df))]=c("x0","y0")
        }
        
      }
    }
    
    if (plot.ellipse == TRUE){
      df.ellipse = data.frame(do.call("rbind", lapply(1 : length(x), function(k){do.call("cbind", coord.ellipse[[k]])})), "Block" = paste0("Block: ", rep(blocks, each = 100)))
      
      names(df.ellipse)[1 : (2*nlevels(group))] = paste0("Col", 1 : (2*nlevels(group)))
    }
    
    if(plot.ellipse == TRUE && class.object[1] %in% c("ipca","sipca","pca","spca","prcomp", "splsda","plsda","mlsplsda"))
      if(is.null(main)){
        df.ellipse$Block="PlotIndiv"}
    else
    {
      df.ellipse$Block=main 
    }
    
    pch.legend=NULL
    if(missing.col)
    {
      
      for (i in 1:nlevels(factor(col))){
        pch.legend=c(pch.legend,df[df$col == levels(factor(col))[i], ]$pch)}
    }
    else{
      for (i in 1:nlevels(group)){
        pch.legend=c(pch.legend,df[df$group == levels(group)[i], ]$pch)}}
    df$pch.legend=pch.legend
    
    
    #-- End: data set
    
    #-- Start: ggplot2
    if (style == "ggplot2"){
      
      #-- Initialise ggplot2
      p = ggplot(df, aes(x = x, y = y, color = group),
                 main = main, xlab = X.label, ylab = Y.label) + theme_bw()
      
      #-- Display sample or row.names
      for (i in levels(group)){
        if (display.names) {
          p = p +geom_point(data = subset(df, group == i),size = 0, shape = 0)+ geom_text(data = subset(df, group == i), aes(label = names), size = 0,show_guide  = F)
        } else {
          p = p + geom_point(data = subset(df, group == i), size = 0, shape = 0)
        }
        if(plot.centroid==TRUE)
        {
          p = p + geom_point(data = subset(df[,c("col","x0","y0","Block","cex","pch","group")], group == i),aes(x=x0,y=y0), size = 0, shape = 0)
        }
      }
      
      #-- Modify scale colour - Change X/Ylabel - split plots into Blocks  
      p = p + scale_colour_manual(values = unique(col.per.group)[match(levels(factor(as.character(group))), levels(group))], name = "Legend", breaks = levels(group))
      p = p + labs(list(title = main, x = X.label, y = Y.label)) + facet_wrap(~ Block, ncol = 2, scales = "free", as.table = FALSE)
      
      #-- xlim, ylim
      if (lim.X) p = p + coord_cartesian(xlim=xlim)
      if (lim.Y) p = p + coord_cartesian(ylim=ylim)
      
      #-- color samples according to col
      if(class.object[1] %in% c("ipca","sipca","pca","spca","prcomp", "splsda","plsda","mlsplsda"))
      {for (i in unique(col)){
        if (display.names) {
          p = p + geom_point(data = subset(df, col == i),size = 0, shape = 0,
                             color = df[df$col == i , ]$col)+geom_text(data = subset(df, col == i), 
                                                                       aes(label = names), 
                                                                       color = df[df$col == i , ]$col,
                                                                       size = df[df$col == i , ]$cex,show_guide  = F)
        } else {
          p = p + geom_point(data = subset(df, col == i), 
                             color = df[df$col == i , ]$col,
                             size = df[df$col == i , ]$cex, 
                             shape = df[df$col == i , ]$pch)
        }
        
        if(plot.centroid==TRUE)
        {
          p = p + geom_point(data = subset(df[,c("col","x0","y0","Block","cex","pch","group")], col == i),aes(x=x0,y=y0), 
                             color = df[df$col == i , ]$col,
                             size = df[df$col == i , ]$cex, 
                             shape = 8)
        }
      }
      }    
      else
      {for (i in unique(col)){
        if (display.names) {
          p = p +geom_point(data = subset(df, col == i),size = 0, shape = 0,
                            color = unique(df[df$col == i & df$Block == paste0("Block: ", blocks[1]), ]$col))+ 
            geom_text(data = subset(df, col == i), 
                      aes(label = names), 
                      color = unique(df[df$col == i & df$Block == paste0("Block: ", blocks[1]), ]$col),
                      size = unique(df[df$col == i & df$Block == paste0("Block: ", blocks[1]), ]$cex),show_guide  = F)
        } else {
          p = p + geom_point(data = subset(df, col == i), 
                             color = unique(df[df$col == i & df$Block == paste0("Block: ", blocks[1]), ]$col),
                             size = unique(df[df$col == i & df$Block == paste0("Block: ", blocks[1]), ]$cex), 
                             shape = unique(df[df$col == i & df$Block == paste0("Block: ", blocks[1]), ]$pch))
        }
        if(plot.centroid==TRUE)
        {
          p = p + geom_point(data = subset(df[,c("col","x0","y0","Block","cex","pch","group")], col == i),aes(x=x0,y=y0), 
                             color = unique(df[df$col == i & df$Block == paste0("Block: ", blocks[1]), ]$col),
                             size = unique(df[df$col == i & df$Block == paste0("Block: ", blocks[1]), ]$cex), 
                             shape = 8)
        }
        
        
      }
      
      
      }
      
      
      #-- Legend
      if (!add.legend) {
        p = p + theme(legend.position="none")
      } else {
        p = p + guides(colour = guide_legend(override.aes = list(shape = if(display.names) {19} else unique(pch.legend), size = unique(cex))))
      }
      
      #-- abline
      if (abline.line)
        p = p + geom_vline(aes(xintercept = 0), linetype = 2, colour = "darkgrey") + geom_hline(aes(yintercept = 0),linetype = 2,colour = "darkgrey")
      
      #-- star
      if (plot.star == TRUE) {
        for (i in 1 : nlevels(group)){
          p = p + geom_segment(data = subset(df, group == levels(group)[i]),
                               aes(x = x0, y = y0,xend=x,yend=y,
                                   label = "Block"),color = unique(col.per.group)[i])
        }
      }
      
      #-- ellipse
      if (plot.ellipse == TRUE) {
        for (i in 1 : nlevels(group)){
          p = p + geom_path(data = df.ellipse,
                            aes_string(x = paste0("Col", 2*(i - 1) + 1), y = paste0("Col", 2 * i),
                                       label = "Block", group = NULL), color = unique(col.per.group)[i])
        }
      }     
      plot(p)
    }
    #-- End: ggplot2
    
    #-- Start: Lattice
    if(style=="lattice") {
      p = xyplot(y ~ x | Block, data = df, xlab = X.label, ylab = Y.label, main = main,
                 group = if (display.names) {names} else {group},
                 scales= list(x = list(relation = "free", limits = xlim),
                              y = list(relation = "free", limits = ylim)),
                 
                 #-- Legend
                 key = if(add.legend == TRUE) {list(space = "right", title = "Legend", cex.title = 1.5,
                                                    point = list(col =  col.per.group),cex=1.3, pch = if(display.names) {16} else unique(pch.legend),text = list(levels(group)))}
                 
                 else {NULL},
                 
                 panel = function(x, y, subscripts, groups, display = display.names,...) {
                   #-- Abline
                   if (abline.line) { panel.abline(v = 0, lty = 2, col = "darkgrey")
                                      panel.abline(h = 0, lty = 2, col = "darkgrey")}
                   
                   #-- Display sample or row.names
                   for (i in 1 : nlevels(group)){
                     
                     if (display){
                       ltext(x = x[group == levels(group)[i]], y = y[group == levels(group)[i]],
                             labels = groups[subscripts & group == levels(group)[i]], col = "white", cex = 0) 
                     } else {
                       lpoints(x = x[group == levels(group)[i]], y = y[group == levels(group)[i]], col = "white", cex = 0, pch = 0)
                     }
                   }
                   
                   
                   
                   
                   
                   
                   #-- color samples according to col
                   for (i in unique(col)){
                     if (display) {
                       ltext(x = x[col == i], y = y[col == i], labels =  groups[subscripts & col == i],
                             col = df[df$col == i, ]$col, cex = df[df$col == i, ]$cex)
                     } else {
                       lpoints(x = x[col == i],  y = y[col == i],
                               col = df[df$col == i, ]$col, cex = df[df$col == i, ]$cex, pch = df[df$col == i, ]$pch)
                     }
                   }
                 })
      print(p) #-- the lattice plot needs to be printed in order to display the ellipse(s)
      
      #-- centroid
      if (plot.centroid) {
        panels = trellis.currentLayout(which = "panel")
        for (k in 1 : length(x)){
          if(class.object[1] %in% c("ipca","sipca","pca","spca","prcomp", "splsda","plsda","mlsplsda"))
          {other=TRUE}
          else
          {
            other=df$Block %in% paste0("Block: ", blocks[k])}
          ind = which(panels == k, arr.ind = TRUE)
          trellis.focus("panel",ind[2], ind[1],highlight = FALSE)
          
          for (i in 1 : nlevels(group)) {
            x0=mean(df[other & group == levels(group)[i] , ]$x)
            y0=mean(df[other & group == levels(group)[i] , ]$y)
            panel.points(x = x0,
                         y = y0,
                         col = unique(col.per.group)[i],pch=8,cex=df[other & group == levels(group)[i] , ]$cex)
          }
        }
        trellis.unfocus()
      }
      
      
      
      #-- star
      if (plot.star) {
        panels = trellis.currentLayout(which = "panel")
        for (k in 1 : length(x)){
          if(class.object[1] %in% c("ipca","sipca","pca","spca","prcomp", "splsda","plsda","mlsplsda"))
          {other=TRUE}
          else
          {
            other=df$Block %in% paste0("Block: ", blocks[k])}
          ind = which(panels == k, arr.ind = TRUE)
          trellis.focus("panel",ind[2], ind[1],highlight = FALSE)
          
          for (i in 1 : nlevels(group)) {
            for (q in 1: length(df[other & group == levels(group)[i]  , "x"]))
            {
              x0=mean(df[other & group == levels(group)[i] , ]$x)
              y0=mean(df[other & group == levels(group)[i] , ]$y)
              panel.segments(x0,y0,df[other & group == levels(group)[i],]$x[q],df[other & group == levels(group)[i],]$y[q], col = unique(col.per.group)[i], cex = df[other & group == levels(group)[i], ]$cex, pch = df[other & group == levels(group)[i], ]$pch)
            }
          }
        }
        trellis.unfocus()
      }
      
      
      
      
      
      #-- ellipse
      if (plot.ellipse) {
        panels = trellis.currentLayout(which = "panel")
        for (k in 1 : length(x)){
          if(class.object[1] %in% c("ipca","sipca","pca","spca","prcomp", "splsda","plsda","mlsplsda"))
          {other.ellipse=TRUE}
          else
          {
            other.ellipse=df.ellipse$Block %in% paste0("Block: ", blocks[k])}
          ind = which(panels == k, arr.ind = TRUE)
          trellis.focus("panel",ind[2], ind[1],highlight = FALSE)
          
          for (i in 1 : nlevels(group)) {
            panel.lines(x = df.ellipse[other.ellipse, paste0("Col", 2*(i - 1) + 1)],
                        y = df.ellipse[other.ellipse, paste0("Col", 2 * i)],
                        col = unique(col.per.group)[i])
          }
        }
        trellis.unfocus()
      }
    }
    #-- End: Lattice
    
    #-- Start: graphics
    if(style=="graphics") {
      
      opar <- par()[! names(par()) %in% c("cin", "cra", "csi", "cxy", "din", "page")]
      par(mfrow=c(1,length(blocks)))
      #-- Define layout
      if (add.legend) 
        par(mai=c( 1.360000, 1.093333, 1.093333,(max(strwidth(levels(group),"inches")))+0.6),xpd=TRUE)
      else 
        par(mar=c(5,4,4,2))
      
      for (k in 1 : length(x)){
        #-- initialise plot
        if(class.object[1] %in% c("ipca","sipca","pca","spca","prcomp", "splsda","plsda","mlsplsda"))
        {titlemain=NULL
         other=TRUE
         if(plot.ellipse)
           other.ellipse=TRUE}
        else
        {titlemain=paste0("Block: ", blocks[k])
         other=df$Block %in% paste0("Block: ", blocks[k])
         if(plot.ellipse)
           other.ellipse=df.ellipse$Block %in% paste0("Block: ", blocks[k])}
        
        plot(df[other, "x" ],
             df[other, "y" ],
             type = "n", xlab = X.label, ylab = Y.label, main = titlemain,
             xlim = c(xlim[[k]][1], xlim[[k]][2]), ylim = c(ylim[[k]][1], ylim[[k]][2]),...)
        
        #-- Display sample or row.names
        for (i in 1 : nlevels(group)){
          if (display.names) {
            text(x = df[group == levels(group)[i] & other, "x"],
                 y = df[group == levels(group)[i] & other, "y"],
                 labels = df[group == levels(group)[i] & other, "names"],
                 col = "white", cex = 0,...)
          } else {
            points(x = df[group == levels(group)[i] & other, "x"],
                   y = df[group == levels(group)[i] & other, "y"],
                   col = "white", cex = 0, pch = 0,...)
          }
        }  
        
        #-- color samples according to col
        for (i in unique(col)){
          if (display.names) {
            text(x = df[df$col == i & other, "x"],
                 y = df[df$col == i & other, "y"],
                 labels = df[df$col == i & other, "names"],
                 col = df[df$col == i, ]$col, cex = df[df$col == i, ]$cex,...)
          } else {
            points(x = df[df$col == i & other, "x"],
                   y = df[df$col == i & other, "y"],
                   col = df[df$col == i, ]$col, cex = df[df$col == i, ]$cex, pch = df[df$col == i, ]$pch,...)
          }
        }
        
        pch.legend=NULL
        if(missing.col)
        {
          
          for (i in 1:nlevels(factor(col))){
            pch.legend=c(pch.legend,df[df$col == levels(factor(col))[i], ]$pch)}
        }
        else{
          for (i in 1:nlevels(group)){
            pch.legend=c(pch.legend,df[df$group == levels(group)[i], ]$pch)}}
        
        if (add.legend) {              
          legend(par()$usr[2]+0.1,par()$usr[4] - (par()$usr[4]-par()$usr[3])/2, col = col.per.group, legend = levels(group), pch = if(display.names) {16} else unique(pch.legend), title = 'Legend', cex = 0.8)
          
        }
        
        #-- Abline
        if (abline.line)
          abline(v = 0, h = 0, lty = 2,...)
        
        #-- Star
        if (plot.star == TRUE) {
          for (i in 1 : nlevels(group)){
            x0=mean(df[group == levels(group)[i] & other, "x"])
            y0=mean(df[group == levels(group)[i] & other, "y"])
            for (q in 1: length(df[group == levels(group)[i] & other, "x"]))
            {
              segments(x0,y0,df[group == levels(group)[i] & other, "x"][q],df[group == levels(group)[i] & other, "y"][q],
                       cex=df$df[group == levels(group)[i] & other, "cex"],col=df[group == levels(group)[i] & other, "col"],...)
            }
          }
        }
        
        #-- Centroid
        if (plot.centroid == TRUE) {
          for (i in 1 : nlevels(group)){
            x0=mean(df[group == levels(group)[i] & other, "x"])
            y0=mean(df[group == levels(group)[i] & other, "y"])
            points(cbind(x0,y0),pch=8,
                   cex=df$df[group == levels(group)[i] & other, "cex"],col=unique(col.per.group)[i],...)
          }
        }
        
        
        #-- Ellipse
        if (plot.ellipse == TRUE) {
          for (i in 1 : nlevels(group)){
            lines(x = df.ellipse[other.ellipse, paste0("Col", 2*(i - 1) + 1)],
                  y = df.ellipse[other.ellipse, paste0("Col", 2 * i)],
                  col = unique(col.per.group)[i],...)
          }
        }
        
        title(main, outer = TRUE, line = -1,...)
      }
      
      opar["usr"]=par()["usr"]
      
      par(opar)
    }
    #-- End: graphics
    
    
    #-- Start: 3d
    if(style=="3d") {
      
      for (k in 1 : length(x)){
        open3d()
        par3d(windowRect = c(500, 30, 1100, 630))
        Sys.sleep(0.1)
        
        if (!is.null(main)) {
          mat = matrix(1:2, 2)
          layout3d(mat, heights = c(1, 10), model = "inherit")
          next3d()
          text3d(0, 0, 0, main)  
          next3d()
        }
        
        if(class.object[1] %in% c("ipca","sipca","pca","spca","prcomp", "splsda","plsda","mlsplsda"))
        {
          other=TRUE
          if(plot.ellipse)
            other.ellipse=TRUE}
        else
        {other=df$Block %in% paste0("Block: ", blocks[k])
         if(plot.ellipse)
           other.ellipse=df.ellipse$Block %in% paste0("Block: ", blocks[k])}
        
        par3d(userMatrix = rotationMatrix(pi/80, 1, -1/(100*pi), 0))
        
        if (add.legend) {
          legend3d(x="right",
                   legend = levels(group), 
                   col = unique(col),
                   pch = rep(16,length(unique(pch))), 
                   pt.cex = unique(cex),
                   bty="n")
        }
        
        #-- Display sample or row.names
        
        for (i in unique(col)){
          if (display.names) {
            text3d(x = df[df$col == i &  other, "x"],
                   y = df[df$col == i &  other, "y"],
                   z = df[df$col == i &  other, "z"],
                   texts = df[df$col == i &  other, "names"],
                   color = df[df$col == i, ]$col, size = unique(df[df$col == i, ]$cex))
          } 
          else {
            cex=unique(df[df$col == i, ]$cex)*20
            switch(unique(df[df$col == i, ]$pch), 
                   sphere = plot3d(x = df[df$col == i & other, "x"],
                                   y = df[df$col == i & other, "y"],
                                   z = df[df$col == i & other, "z"],
                                   col = df[df$col == i, ]$col, size = cex, radius = cex/20, add = TRUE),
                   tetra = shapelist3d(tetrahedron3d(), x = df[df$col == i &other, "x"],
                                       y = df[df$col == i & other, "y"],
                                       z = df[df$col == i & other, "z"],
                                       col = df[df$col == i, ]$col, size = cex/25),
                   cube = shapelist3d(cube3d(),x = df[df$col == i & other, "x"],
                                      y = df[df$col == i & other, "y"],
                                      z = df[df$col == i & other, "z"],
                                      col = df[df$col == i, ]$col, size = cex/30),
                   octa = shapelist3d(octahedron3d(), x = df[df$col == i & other, "x"],
                                      y = df[df$col == i & other, "y"],
                                      z = df[df$col == i & other, "z"],
                                      col = df[df$col == i, ]$col, size = cex/17),
                   icosa = shapelist3d(icosahedron3d(), x = df[df$col == i & other, "x"],
                                       y = df[df$col == i & other, "y"],
                                       z = df[df$col == i &other, "z"],
                                       col = df[df$col == i, ]$col, size = cex/20),
                   dodeca = shapelist3d(dodecahedron3d(), x = df[df$col == i &other, "x"],
                                        y = df[df$col == i & other, "y"],
                                        z = df[df$col == i & other, "z"],
                                        col = df[df$col == i, ]$col, size = cex/20))
          }
        }
        
        #-- Ellipse
        if(plot.ellipse)
        {
          coords=matrix(cbind(df[other, "x"],
                              df[other, "y"],
                              df[other,"z"]),ncol=3)
          centr.coords <- apply(coords, 2, function(x) tapply(x, group, mean))
          if(length(unique(group)) == 1)
            centr.coords <- matrix(centr.coords, nrow=1)
          
          
          rownames(centr.coords) <- levels(group)
          lg <- levels(group)
          for(i in 1:length(lg)) {
            
            g   <- lg[i]
            
            sel <- group == g
            
            s   <- cov(coords[sel,,drop=FALSE])
            cc  <- centr.coords[i,]
            
            
            # lines(ellipse(s, centre=cc), col=unique(col.per.group)[i])
            shade3d(ellipse3d(s, centre=cc, level=ellipse.level), col=unique(col.per.group)[i], alpha=alpha)
            
          }
        }
        
        #-- Centroid
        if (plot.centroid == TRUE) {
          for (i in 1 : nlevels(group)){
            x0=mean(df[group == levels(group)[i] & other, "x"])
            y0=mean(df[group == levels(group)[i] & other, "y"])
            z0=mean(df[group == levels(group)[i] & other, "z"])
            points3d(x=x0,y=y0,z=z0,
                     cex=df$df[group == levels(group)[i] & other, "cex"],col=unique(col.per.group)[i])
          }
        }
        
        #-- Star
        if (plot.star == TRUE) {
          for (i in 1 : nlevels(group)){
            x0=mean(df[group == levels(group)[i] & other, "x"])
            y0=mean(df[group == levels(group)[i] & other, "y"])
            z0=mean(df[group == levels(group)[i] & other, "z"])
            for (q in 1: length(df[group == levels(group)[i] & other, "x"]))
            {
              segments3d(x=c(x0,df[group == levels(group)[i] & other, "x"][q]),y=c(y0,df[group == levels(group)[i] & other, "y"][q]),
                         z=c(z0,df[group == levels(group)[i] & other, "z"][q]),
                         cex=df$df[group == levels(group)[i] & other, "cex"],col=df[group == levels(group)[i] & other, "col"])
            }
          }
        }
        
        
        #-- draws axes/box --#
        if (axes.box == "box") {
          axes3d(marklen = 25)
          box3d()
        }
        if (axes.box == "bbox") {
          bbox3d(color = c("#333377", "black"), emission = gray(0.5), 
                 specular = gray(0.1), shininess = 5, alpha = 0.8, marklen = 25)
        }
        if (axes.box == "both")  {
          axes3d(marklen = 25); box3d()
          bbox3d(color = c("#333377", "black"), emission = gray(0.5), 
                 specular = gray(0.1), shininess = 5, alpha = 0.8, marklen = 25)         
        }
        #-- add axes labels --#
        mtext3d(X.label, "x-+", line = 1)
        mtext3d(Y.label, "y-+", line = 1.5)
        mtext3d(Z.label, "z+-", line = 1)
        if( ! class.object[1] %in% c("ipca","sipca","pca","spca","prcomp", "splsda","plsda","mlsplsda"))
          title3d(main=paste0("Block: ", blocks[k]))
      }
      #-- output --#
      return(invisible(cbind(x, y, z)))
      
    }
    #-- End: 3d
    
  }