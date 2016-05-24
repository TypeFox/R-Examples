plotEig <- function(eigvalue, nf, xax = 1, yax = 2, col.plot = "black", col.kept = "grey", col = "white", facets = NULL, plot = TRUE, storeData = FALSE, pos = -1, ...) {
  ## prepare
  col <- rep(col, length(eigvalue))
  col[nf] <- col.kept
  col[c(xax, yax)] <- col.plot
 
  ## parameters management
  sortparameters <- sortparamADEg(...)
  params <- list()
  params$adepar <- list(ppolygons = list(col = col), p1d = list(horizontal = FALSE), psub = list(position = "topright"), pgrid = list(draw = FALSE), pbackground = list(box = FALSE))
  sortparameters$adepar <- modifyList(params$adepar, sortparameters$adepar, keep.null = TRUE)
  
  if(is.null(facets) || isTRUE(sortparameters$g.args$samelimits)) {
    lim <- c(0, ifelse(is.null(facets), length(eigvalue), max(table(facets)))) + 0.5
    if(isTRUE(sortparameters$adepar$p1d$horizontal))
    	params$g.args <- list(ylim = lim)
    else
      params$g.args <- list(xlim = lim)
    
  	lim.val <- range(eigvalue)
  	if(lim.val[1] >= 0) {
	    lim.val <- c(0, lim.val[2] + diff(c(lim.val[1], lim.val[2])) / 10)
    	if(isTRUE(sortparameters$adepar$p1d$horizontal))
	  	  params$g.args <- list(xlim = lim.val, ylim = params$g.args$ylim)
    	else
	      params$g.args <- list(xlim = params$g.args$xlim, ylim = lim.val)
	  }
  } else {
    params$g.args <- list(xlim = NULL, ylim = NULL)
  }
  
  sortparameters$g.args <- modifyList(params$g.args, sortparameters$g.args, keep.null = TRUE)
  do.call("s1d.barchart", c(list(score = substitute(eigvalue), pos = pos - 2, plot = plot, facets = facets, storeData = storeData), sortparameters$adepar, sortparameters$trellis, sortparameters$g.args, sortparameters$stats, sortparameters$s.misc, sortparameters$rest))
}


## si ADEgS contenu dans un plus petit espace;
## oldposition: matrice de position: nrow:number of graphs, col: x0, y0, x1, y1
## newposition: vector, length 4: x0, y0, x1, y1
## Calcul: toute les oldpositions: dans newposition.
## renvoie d'une matrice, 4col, nrow(oldposition) rows.
## cette indique les nouvelles positions des graphiques dans le reférentiel de refposition
## test:
## oldpos <- t(rbind(rep(c(0, 1 / 3, 2 / 3), 2), c(rep(0.5, 3), rep(0, 3)), rep(c(1 / 3, 2 / 3, 1), 2), c(rep(1, 3), rep(0.5, 3))))
## newpos <- c(0.5, 0.5, 1, 1)
## .updateadegsposition(oldpos, refpositions)
.updateADEgSposition <- function(oldposition, refposition) {
  ## test arguments
  if(NCOL(oldposition) > 4)
    stop("wrong position, only 4columns expected")
  if(any(oldposition[, 1] >= oldposition[, 3]))
    stop("wrong position given, x0>=x1 cannot work")
  if(any(oldposition[, 2] >= oldposition[, 4]))
    stop("wrong position given, y0>=y1 cannot work")
  if(NCOL(refposition) != 1)
    stop("error in .updateADEgSposition, several containing graphs given, only one possible")  ## ne devrait jamais jamais arriver!
  
  ## formula:
  ## xnewi <- xoldi * wnew + x0new
  ## ynewi <- yoldi * hnew + y0new
  x0o <- oldposition[, 1]
  x1o <- oldposition[, 3]
  y0o <- oldposition[, 2]
  y1o <- oldposition[, 4]
  wnew <- refposition[3] - refposition[1]
  hnew <- refposition[4] - refposition[2]
  ## peut mieux faire (optimisation)
  calcNew <- function(old, new, wh) {return(old * wh + new)}
  return(cbind(calcNew(x0o, refposition[1], wnew),
               calcNew(y0o, refposition[2], hnew),
               calcNew(x1o, refposition[1], wnew),
               calcNew(y1o, refposition[2], hnew)))
}


## .getposition: mainly for placing eigen plot.
## gives coordinates according to string position and width, height wanted
.getposition <- function(position, w = 0.25, h = 0.25) {
  if(is.numeric(position) && length(position) == 4)	
    posnum <- position
  else if(is.numeric(position) && length(position) == 2)
    posnum <- c(position[1], position[2], position[1] +  w, position[2] + h)
  else if(is.character(position)) {
    position <- match.arg(position[1], choices = c("bottomleft", "bottomright", "topleft", "topright", "none"), several.ok = FALSE)
    if(position == "bottomleft")
      posnum <- c(0.0, 0.0, w, h)
    else if(position == "bottomright")
      posnum <- c(1 - w, 0.0, 1, h)
    else if(position == "topleft")
      posnum <- c(0.0, 1 - h, w, 1)
    else if(position == "topright")
      posnum <- c(1 - w, 1 - h, 1, 1)
    else if(position == "none")
      posnum <- NULL
    else
      stop("Wrong indication of position")
  }
  else 
    stop("Wrong indication of position")
  return(posnum)
}


## pour adeGs, on doit etre capable de separer facilement les parametres pour pouvoir avoir un adressage specifique pour chaque graphique (ie pas la meme chose poru le sarrow et le slabel dans un scatterdudi)
## selon les graphiques adeGs nous aurons des pattern differents:
## ex pour scatter.dudi, nous pouvons imager 'col', 'row', 'eigen' pour distinguer les paramètres spécifiques au graph
.partoadeg <- function(..., pattern = NULL) {
  if(is.null(pattern))
    stop("error in .partoadeg, pattern should be filled")
  if(try(is.list(...), silent = TRUE) == TRUE)
    dots <- as.list(...)
  else dots <- list(...)
  result <- vector("list", length = length(pattern))
  result <- lapply(result, as.list)
  names(result) <- pattern
	## si deja indique en list
  if(length(dots)) {
    whichG <- c()
    then <- c()
    ## pour ceux indiquer avec des .
    splitgrp <- sapply(names(dots), FUN = function(x) {strsplit(x, ".", fixed=TRUE)})
    for(i in 1:length(splitgrp)) {
      ## premier niveaux quel graph
      whichG <- c(whichG, splitgrp[[i]][1])
      ## deuxieme niveau si il y a le nom suivant (qui etait colle avec un .)
      if(length(splitgrp[[i]]) > 1) { ## un second element
        then <- c(then, paste(splitgrp[[i]][2:length(splitgrp[[i]])], collapse = "."))
      }
      else
        then <- c(then, NA)
    }
    indix <- pmatch(whichG, pattern, duplicates.ok = TRUE)
    notna <- which(!is.na(indix)) ## ne garder que les non na 
    arena <- which(is.na(indix))  ## position dans indix des NA ie: ceux qui n'ont pas de match
    for(i in 1:length(result)) {
      sublist <- result[[i]]    ## sous list deja trouve... a priori list
      if(any(indix[notna] == i)) {  ## si un indix vaut i=> a mettre dans result
                                    ## soit dire une liste soiton a dans then
        toselect <- which(indix == i)
        for(have2 in 1:length(toselect))
          if(!is.na(then[toselect[have2]])) { ## a ete renseigne avec un point ensuite
            newlist <- c(list(), dots[toselect[have2]])
            names(newlist) <- then[toselect[have2]]
            sublist <- c(sublist, newlist)
          }
          else  ## c un na na donc ensuite on avait une liste
            sublist <- c(sublist, dots[[toselect[have2]]])
      }
      if(length(arena)) ## on a en plus des na, donc des parameteres pour tous
        selectNa <- indix[arena]
      sublist <- c(sublist, dots[arena])
      if(!is.null(sublist))
        result[[i]] <- sublist
    }}
  return(result)
}

.n2mfrow <- function(nr.plots) {
  ## inspired by n2mfrow but we change the default when the number of graphs is <6
  if (nr.plots <= 3) 
    c(1, nr.plots)
  else if (nr.plots <= 6) 
    c(2, (nr.plots + 1) %/% 2)
  else if (nr.plots <= 9)
    c((nr.plots + 2) %/% 3, 3)
  else if (nr.plots <= 12)
    c((nr.plots + 3) %/% 4, 4)
  else c(nrow <- ceiling(sqrt(nr.plots)), ceiling(nr.plots / nrow))
}

## Get positions matrix for ADEgs according  a given layout
## strongly inspired by the layout function
## ng: number of positions to get
layout2position <- function(mat, widths = rep(1, NCOL(mat)), heights = rep(1, NROW(mat)), ng, square = FALSE) {
  if(is.vector(mat)) {
    if(missing(ng)) ng <- mat[1] * mat[2]
    mat <- matrix(c(1:ng, rep(0, length.out = ((mat[1] * mat[2]) - ng))), nrow = mat[1], byrow = TRUE)
    if(missing(widths))
      widths <- rep(1, ncol(mat))
    if(missing(heights))
      heights <- rep(1, nrow(mat))
  }
  if(NROW(mat) != length(heights)) stop("wrong number of heigths given", call. = FALSE)
  if(NCOL(mat) != length(widths)) stop("wrong number of widths given", call. = FALSE) 
  nbgraph <- max(mat)
  ## get xi position and yi position
  xi <- c(0)
  yi <- c(0)
  ## here, width given such as proportional colums.
  ## so the sum(width)/length(widths) == 1
  ## more units to take in account"
  if(square == TRUE) {
    wi <- widths / max(length(widths), length(heights))
  	hi <- heights / max(length(widths), length(heights))
  } else {    
  	wi <- widths / sum(widths)
  	hi <- heights / sum(heights)
  }
  
  ## layout from left to right, up to bottom 
  for(i in 1:length(wi))
    xi <- c(xi, xi[i] + wi[i])
  for(i in 1:length(hi))
    yi <- c(yi, yi[i] + hi[i])
  
  yi <- rev(yi)
  pos <- c()
  for(i in 1:nbgraph) { ## for each graph, get the positions as x0, y0, x1, y1
    indx <- which(mat == i, arr.ind = TRUE)
    if(length(indx) == 0) { ## just in case
      warning(paste("in layout2position, a graph position is missing, no graph", i, "defined", sep = " "), call. = FALSE)
      pos <- rbind(pos, rep(0, 4))
    }
    else
      pos <- rbind(pos, c(xi[min(indx[, 2])], yi[(max(indx[, 1]) + 1)], xi[(max(indx[, 2]) + 1)], yi[min(indx[, 1])]))
  }
  return(pos)
}


## For analysis plot (ADEgS creation)
sortparamADEgS <- function(..., graphsnames, nbsubgraphs = rep(1, length(graphsnames))) {
  seppara <- .partoadeg(..., pattern = graphsnames)
  sortparameters <- lapply(seppara, FUN = sortparamADEg)
  alist <- function(x) {
    aa <- list()
    for(i in 1:length(x))
      aa <- c(aa, x[[i]])
    aa
  }
  tomerge <- lapply(sortparameters, alist)
  oki <- lapply(tomerge, .mergingList)
  if(!all(nbsubgraphs == rep(1, length(graphsnames))))
    for (i in 1:length(nbsubgraphs))
      oki[[i]] <- repList(oki[[i]], nbsubgraphs[i])
  return(oki)
}

