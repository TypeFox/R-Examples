## TODO
# table ne passe pas: is.numeric mais pas is.vector. Passe avec table(x, y) : donne une matrice (TRUE avec is.matrix())
# quelque chose pour les attributs d'un vector : sinon regexp ne passe pas
# mettre une première ligne qui donne le code invoqué.
# : matrice : largeur des colonnes quand les indices sont à deux chiffres
# subset: argument as formula

# TODO : data frame and matrix share too much code

debuging <- FALSE;

drawExpression <- function (expr, draw.index=FALSE, draw.names=FALSE, filename=NULL) {
  if (mode(expr) != "character") {
    stop("expr must be a characters string");
  }

  if (length(expr) != 1) {
    stop("expr must be a characters string of length 1");
  }

  ## Build the syntax tree
  e <- parse(text=expr);
  if (mode(e) != "expression") {
    stop(paste(expr, "must be a parsable expression"));
  }

  ## Create an intermediary representation
  drawable <- .drawableTree(e[[1]], 1);

  plot.new();

  ## draw this representation with grid function
  .drawTree(drawable, filename=filename);

}

########################################################
########################################################

# ---------------------------------
# Walk recursively through the syntax tree and populate the "drawable tree"
# structure accordingly.
#
# param
# -----
# call: the syntax tree.
# level: the level in the tree (1-based)
#
# value
# -----
# A list which is the root node of the tree. Each node in the tree is a nammed
# list; children of a node are listed in the "$children" entry of the node.
# ---------------------------------
.drawableTree <- function(call, level) {
  #mode(call) may be name, call, or primitive (numeric, etc.).
  l <- list();
  l$eval = eval(call, envir=parent.frame(2));
  l$type = "";
  l$level = level;
  lengthCall <- length(call);
  if (lengthCall > 1) {
    children = list();
    if (mode(eval(call[[1]])) == "function") {
      if (as.character(call[[1]]) == "[") {
        children[[1]] <- .drawableTree(call[[2]], level+1);
        children[[2]] <- makeOpenningBracket(level + 1);
        for (z in 3:length(call)) {
          if (z > 3) {
            children[[length(children) + 1]] <- makeComma(level+1);
          }
          if (as.character(call[[z]]) == "") {
          } else {
            children[[length(children) + 1]] <- .drawableTree(call[[z]], level+1);
          }
        }
        children[[length(children) + 1]] <- makeClosingBracket(level + 1);
      } else if (as.character(call[[1]]) == "$") {
        children[[1]] <- .drawableTree(call[[2]], level+1);
        children[[2]] <- makeDollar(level + 1);
        children[[3]] <- makeFunction(deparse(call[[3]]), level+1);
      } else if (as.character(call[[1]]) == "[[") {
        children[[1]] <- .drawableTree(call[[2]], level+1);
        children[[2]] <- makeOpenningDoubleBracket(level + 1);
        children[[3]] <- .drawableTree(call[[3]], level+1);
        children[[4]] <- makeClosingDoubleBracket(level + 1);
      } else if (
          (as.character(call[[1]]) == "<")
          || (as.character(call[[1]]) == ":")
          || (as.character(call[[1]]) == ">")
          || (as.character(call[[1]]) == "+")
          || (as.character(call[[1]]) == "-")
          || (as.character(call[[1]]) == "*")
          || (as.character(call[[1]]) == "/")
          || (as.character(call[[1]]) == "==")
          || (as.character(call[[1]]) == ">=")
          || (as.character(call[[1]]) == "<=")
          || (as.character(call[[1]]) == "&")
          || (as.character(call[[1]]) == "|")
          || (as.character(call[[1]]) == "!=")
          || (as.character(call[[1]]) == "$")
          ) {
        children[[1]] <- .drawableTree(call[[2]], level+1);
        children[[2]] <- makeOperator(as.character(call[[1]]), level + 1);
        children[[3]] <- .drawableTree(call[[3]], level+1);
      } else if (grepl("<-$", as.character(call[[1]]))) {
        children[[1]] <- makeFunction(deparse(call[[2]]), level+1);
        children[[2]] <- makeFunction(as.character(call[[1]]), level + 1);
        for (i in 3:lengthCall) {
          children[[i]] <- .drawableTree(call[[i]], level+1);
        }
      } else {
        #print(as.character(call[[1]]));
        children[[1]] <- makeFunction(as.character(call[[1]]), level + 1);
        children[[2]] <- makeOpenningParenthesis(level + 1);
        offset <- 1;
        for (i in 2:lengthCall) {

        #  print(as.character(call[[i]]));

# test if it is the name of a parameter (such as "from" in "seq(from=...)")
          paramName <- names(call[i]);
        #  print(paramName);
          if (! is.null(paramName)) {
            children[[i + offset]] <- makeParam(names(call[i]), level + 1);
            offset = offset + 1;
          }

# This argument is a function, handle it nicely.
        if (mode(eval(call[[i]])) == "function") {
          children[[i + offset]] <- makeFunction(as.character(call[[i]]), level + 1);
        } else {
          children[[i + offset]] <- .drawableTree(call[[i]], level+1);
        }
        #  children[[i + offset]] <- .drawableTree(call[[i]], level + 1);

# add a comma between argument
          if (i < lengthCall) {
            offset = offset + 1;
            children[[i + offset]] <- makeComma(level + 1);
          }
        }
        children[[length(children) + 1]] <- makeClosingParenthesis(level + 1);
      }
    } else {
      for (i in 1:lengthCall) {
        children[[i]] <- .drawableTree(call[[i]], level+1);
      }
    }
    l$children <- children;
  } else {
    l$children = NULL;
  }
  l;
}

########################################################
########################################################

.drawTree <- function(drawable, filename) {
  depth <- getMaxDepth(drawable);
  linesGrob <- vector(length=depth, mode="list");
  for (i in depth:1) {
    line <- list();
    line <- .getLineRec(line, drawable, i);
    lineGrob <- lineBoxGrob(line, i, depth);
    linesGrob[[i]] <- lineGrob;
  }

  heights <- numeric(length(linesGrob));
  widths <- numeric(length(linesGrob));
  for (i in 1:length(linesGrob)) {
    if (i == 1) {
      heights <- grobHeight(linesGrob[[1]]);
      widths <- grobWidth(linesGrob[[1]]);
    } else {
      heights <- unit.c(heights, grobHeight(linesGrob[[i]]));
      widths <- unit.c(widths, grobWidth(linesGrob[[i]]));
    }
  }
  totalHeight <- sum(heights);
  maxWidth <- max(widths);

  if(!is.null(filename)) {
    height_inches <- as.numeric(convertUnit(totalHeight, "inches"));
    width_inches <- as.numeric(convertUnit(maxWidth, "inches"));
    pdf(filename, height=height_inches + .2, width=width_inches + .2, onefile=TRUE);
  }

  listvp <- viewport(x=.5, y=.5, width=maxWidth, height=totalHeight);
  pushViewport(listvp);
  if (debuging) grid.rect(gp=gpar(lty="longdash"));
  for (i in length(linesGrob):1) {
    line <- linesGrob[[i]];
    line$y <- sum(heights[1:i]) - heights[i] * 0.5;
    grid.draw(line);
  }
  popViewport();

  if(!is.null(filename)) {
    dev.off();
  }
}

# A line is a list of representations of R syntax components: R object or
# pieces of syntax (function name, coma, operator...).
.getLineRec <- function(line, drawable, level) {
  if ((drawable$level < level) & (! is.null(drawable$children))) {
    for (i in 1:length(drawable$children)) {
      line <- .getLineRec(line, drawable$children[[i]], level);
    }
  } else {
    line <- c(line, list(drawLineComponent(drawable)));
  }
  return(line);
}

drawLineComponent <- function(drawable) {
  if (drawable$type == "special") {
    return(functionTextBoxGrob(drawable$eval));
  } else {
    if (is.atomic(drawable$eval)) {
      if (is.vector(drawable$eval)) {
        return(vectorBoxGrob(drawable$eval));
      } else if (is.matrix(drawable$eval)) {
        return(matrixBoxGrob(drawable$eval));
      } else if (is.array(drawable$eval)) {
        ## TODO good for tapply
        #return(listBoxGrob(as.list(drawable$eval)));
        ## TODO good for margin.table and table.
        ## not so good for table(): lost its name
        return(matrixBoxGrob(as.matrix(drawable$eval)));
      } else if (is.factor(drawable$eval)) {
        return(vectorBoxGrob(drawable$eval, a.factor=TRUE));
      } else {
        stop(paste("unknown case", drawable$eval));
      }
    } else if (is.data.frame(drawable$eval)) {
      return(dataframeBoxGrob(drawable$eval));
    } else if (is.list(drawable$eval)) {
      return(listBoxGrob(drawable$eval));
    } else if (is.function(drawable$eval)) {
      return(functionTextBoxGrob(drawable$eval));
      # TODO
    } else {
      stop(paste("Unknown type:", drawable$eval));
    }
  }
}

makeDefaultSpecial <- function(text, level) {
  l <- list();
  l$eval = text;
  l$type = "special";
  l$level = level;
  l$children = NULL;
  l;
}

makeOpenningDoubleBracket <- function(level) {
  makeDefaultSpecial("[[", level);
}

makeClosingDoubleBracket <- function(level) {
  makeDefaultSpecial("]]", level);
}

makeDollar <- function(level) {
  makeDefaultSpecial("$", level);
}

makeOpenningBracket <- function(level) {
  makeDefaultSpecial("[", level);
}

makeClosingBracket <- function(level) {
  makeDefaultSpecial("]", level);
}

makeOperator <- function(op, level) {
  makeDefaultSpecial(op, level);
}

makeOpenningParenthesis <- function(level) {
  makeDefaultSpecial("(", level);
}

makeClosingParenthesis <- function(level) {
  makeDefaultSpecial(")", level);
}

makeParam <- function (param, level) {
  makeDefaultSpecial(paste(param, " = ", sep=""), level);
}

makeFunction <- function (functionName, level) {
  makeDefaultSpecial(functionName, level);
}

makeComma <- function (level) {
  makeDefaultSpecial(",", level);
}

getMaxDepth <- function(drawable) {
  depth <- drawable$level;
  if (! is.null(drawable$children)) {
    for (i in 1:length(drawable$children)) {
      depth <- max(depth, getMaxDepth(drawable$children[[i]]));
    }
  }
  depth;
}

getMaxHeightForRaw <- function(drawable, level, height) {
  if ((drawable$level < level) & (! is.null(drawable$children))) {
    for (i in 1:length(drawable$children)) {
      height <- max(height, getMaxHeightForRaw(drawable$children[[i]], level, height));
    }
  } else {
    height <- max(height, drawable$height);
  }
  print(paste(drawable$eval, ":", height));
# add margin
  height <- height + 0.005;
}

############## ############## ############## ############## ##############
############## ############## ############## ############## ##############
############## ############## ############## ############## ##############

############## ############## ############## ############## ##############
# grid object (Grob)
############## ############## ############## ############## ##############

objectGrob <- function(obj) {
 if (is.data.frame(obj)) {
   return(dataframeBoxGrob(obj));
 } else if (is.list(obj)) {
   return(listBoxGrob(obj));
 } else if (is.factor(obj)) {
   return(vectorBoxGrob(as.character(obj)));
 } else if (is.vector(obj)) {
   return(vectorBoxGrob(obj));
 } else if (is.matrix(obj)) {
   return(matrixBoxGrob(obj));
 }
 stop(paste("Object not known", mode(obj)));
}

############## ############## ############## ############## ##############
# Draw a line in the graphic
############## ############## ############## ############## ##############

draw.lineBox <- function(l, x=.5, y=.5, height, width, components, comp.height, comp.width, draw.index=FALSE, draw.names=FALSE, margin) {
  linevp <- viewport(x=x, y=y, width=width, height=height);
  pushViewport(linevp);

  for (i in 1:length(components)) {
    obj <- components[[i]];
    obj$x <- unit(margin, "mm") * i - unit(margin, "mm") + sum(comp.width[1:i]) - comp.width[i] * 0.5;
    obj$y <- comp.height[i] * 0.5 + unit(margin, "mm");

    # if (i < length(components)) { 
    #   grid.lines(
    #   x=unit(4, "mm") * i + sum(comp.width[1:i]),
    #   gp=gpar(lty="dashed")
    #   );
    # }
    grid.draw(obj);
  }

  popViewport();
}

lineBoxGrob <- function(l, draw.index=FALSE, draw.names=FALSE, margin=2) { 
  components <- list();
  comp.height <- vector();
  comp.width <- vector();
  for (i in 1:length(l)) {
    components[[i]] <- l[[i]];
    if (i == 1) {
      comp.height <- grobHeight(components[[i]]);
      comp.width <- grobWidth(components[[i]]);
    } else {
      comp.height <- unit.c(comp.height, grobHeight(components[[i]]));
      comp.width <- unit.c(comp.width, grobWidth(components[[i]]));
    }
  }
  height <- unit(margin, "mm") + max(comp.height)
  width <- unit(margin, "mm") * length(components) + sum(comp.width);

  grob(labels=l, components=components, height=height, width=width, comp.height=comp.height, comp.width=comp.width, x=.5, y=.5, margin=margin, cl="lineBox");
}

# level is the current level in the syntax tree, nlevel the total number of levels
# , level, nlevel, ...
drawDetails.lineBox <- function(x, recording) {
  draw.lineBox(x$labels, x$x, x$y, x$height, x$width, x$components, x$comp.height, x$comp.width, margin=x$margin);
}

xDetails.lineBox <- function(x, theta) {
  grobX(roundrectGrob(x=x$x, y=x$y, width=x$width, height=x$height), theta);
}

yDetails.lineBox <- function(x, theta) {
  grobY(rectGrob(x=x$x, y=x$y, width=x$width, height=x$height), theta);
}

grobWidth.lineBox <- function(x) {
  x$width
}

grobHeight.lineBox <- function(x) {
  x$height
}

############## ############## ############## ############## ##############
# Draw a function name (and parenthesis, comma, operator)
############## ############## ############## ############## ##############

draw.functionText <- function(l, x=.5, y=.5, height, width) {
  functionTextVP <- viewport(x=x, y=y, width=width, height=height);
  pushViewport(functionTextVP);
  #grid.rect(x=x, y=y, width=width, height=height, gp=gpar(lty="longdash"));
  #grid.rect(gp=gpar(lty="longdash"));
  grid.text(l);
  popViewport();
}

functionTextBoxGrob <- function(l, x=.5, y=.5) { 
  height <- unit(1, "line");
  width <- stringWidth(l) + unit(1, "mm");
  grob(labels=l, height=height, width=width, x=x, y=y, cl="functionText");
}

drawDetails.functionText <- function(x, ...) {
  draw.functionText(x$labels, x$x, x$y, x$height, x$width);
}

xDetails.functionText <- function(x, theta) {
  grobX(roundrectGrob(x=x$x, y=x$y, width=x$width, height=x$height), theta);
}

yDetails.functionText <- function(x, theta) {
  grobY(rectGrob(x=x$x, y=x$y, width=x$width, height=x$height), theta);
}

grobWidth.functionText <- function(x) {
  x$width
}

grobHeight.functionText <- function(x) {
  x$height
}

############## ############## ############## ############## ##############
# Draw a list
############## ############## ############## ############## ##############

draw.listBox <- function(l, x=.5, y=.5, height, width, components, comp.height, comp.width, draw.index=FALSE, draw.names=FALSE, marginheight) {
  listvp <- viewport(x=x, y=y, width=width, height=height);
  pushViewport(listvp);

  # content
  content.height <- height - marginheight;
  content.vp <- viewport(
      x=width * .5,
      y=content.height * .5,
      width=width,
      height=content.height
      );
  pushViewport(content.vp);
  grid.rect(gp=gpar(lty="dashed"));
  for (i in 1:length(components)) {
    obj <- components[[i]];
    obj$x <- unit(4, "mm") * i - unit(2, "mm") + sum(comp.width[1:i]) - comp.width[i] * 0.5;
    obj$y <- comp.height[i] * 0.5 + unit(2, "mm");

    if (i < length(components)) { 
      grid.lines(
      x=unit(4, "mm") * i + sum(comp.width[1:i]),
      gp=gpar(lty="dashed")
      );
    }
    grid.draw(obj);
  }
  popViewport();

# names
  if (draw.names | draw.index) {
    margin.vp <- viewport(
        x=width * .5,
        y=content.height + marginheight * .5,
        width=width,
        height=marginheight
        );
    pushViewport(margin.vp);
    grid.rect(gp=gpar(fill="lightgray", lwd=0))
      if (draw.names) {
#y <- marginheight * .5;
        y <- marginheight ;
# if (draw.index) {
#   y <- (marginheight - unit(1, "lines")) * .5;
# }
        for (j in 1:length(l)) {
          x <- unit(4, "mm") * j - unit(2, "mm") + sum(comp.width[1:j]) - comp.width[j] * 0.5;
          grid.text(names(l)[j],
              y=y,
              x=x,
#hjust="center",
              vjust="bottom",
              just="right",
              rot=60
              );
        }
      }
    if (draw.index) {
      for (j in 1:length(l)) {
          x <- unit(4, "mm") * j - unit(2, "mm") + sum(comp.width[1:j]) - comp.width[j] * 0.5;
        grid.text(j,
            y=unit(1, "lines") * .5,
            x=x,
            just="right"
            );
      }
    }
    popViewport();
  }


  popViewport();
}

listBoxGrob <- function(l, x=.5, y=.5, draw.index=FALSE, draw.names=FALSE) { 

  if (!is.null(names(l))) {
    draw.index=TRUE;
    draw.names=TRUE;
  }

  components <- list();
  comp.height <- vector();
  comp.width <- vector();

  if (length(l) == 0) {
    stop("empty list");
  }
  for (i in 1:length(l)) {
    components[[i]] <- objectGrob(l[[i]]);
    if (i == 1) {
      comp.height <- grobHeight(components[[i]]);
      comp.width <- grobWidth(components[[i]]);
    } else {
      comp.height <- unit.c(comp.height, grobHeight(components[[i]]));
      comp.width <- unit.c(comp.width, grobWidth(components[[i]]));
    }
  }
  height <- unit(4, "mm") + max(comp.height)
  width <- unit(4, "mm") * length(components) + sum(comp.width);

  marginheight <- unit(0, "mm");
  if (draw.index) {
    marginheight <- marginheight + unit(1, "lines");
  }
  if (draw.names) {
    if (is.null(names(l))) {
      stop("Cannot draw names if the list has no column names.");
    }
    imax <- which.max(nchar(names(l)));
    max.name.length <- stringWidth(names(l)[imax]);
    marginheight <- marginheight + max.name.length;
  }
  height <- height + marginheight;

  grob(labels=l, components=components, height=height, width=width, comp.height=comp.height, comp.width=comp.width, x=x, y=y, draw.index=draw.index, draw.names=draw.names, cl="listBox", marginheight=marginheight);
}

drawDetails.listBox <- function(x, ...) {
  draw.listBox(x$labels, x$x, x$y, x$height, x$width, x$components, x$comp.height, x$comp.width, x$draw.index, x$draw.names, marginheight=x$marginheight);
}

xDetails.listBox <- function(x, theta) {
  grobX(roundrectGrob(x=x$x, y=x$y, width=x$width, height=x$height), theta);
}

yDetails.listBox <- function(x, theta) {
  grobY(rectGrob(x=x$x, y=x$y, width=x$width, height=x$height), theta);
}

grobWidth.listBox <- function(x) {
  x$width
}

grobHeight.listBox <- function(x) {
  x$height
}

############## ############## ############## ############## ##############
# Draw a vector
############## ############## ############## ############## ##############

draw.vectorBox <- function(vect, x=.5, y=.5, draw.index=FALSE, draw.names=FALSE, width, height) {
  len.vect <- length(vect);

  if (draw.names) {
    v.names <- names(vect);
    max.name.length <- which.max(nchar(v.names));
  }
  if (draw.index) {
    indexes <- 1:len.vect;
  }

  tablevp <- viewport(x=x, y=y, width=width, height=height);
  pushViewport(tablevp);
  #grid.rect(gp=gpar(lty="dashed"));
  grid.lines(y=unit(0, "lines"));
  grid.lines(y=unit(1, "lines"));
  #grid.lines(y=unit(2, "lines"), gp=gpar(lty="dashed"));

  for (i in 1:(len.vect)) {
    if (i == 1) { 
      grid.lines(x=unit(0, "npc"), y=unit(c(0,1), "lines"));
    } else { 
      grid.lines(x=sum(stringWidth(vect[1:(i-1)])) + unit(2, "mm") * (i-1), y=unit(c(0,1), "lines"));
    }
  }
  grid.lines(x=unit(1, "npc"), y=unit(c(0,1), "lines"));

  for (i in 1:(len.vect)) {
    grid.text(vect[i],
        y=unit(2.5, "mm"),
        x=unit(2, "mm") * (i-1) + unit(1, "mm") + sum(stringWidth(vect[1:i])),
        just="right"
        );
  }


#grid.rect(x = unit(0.5, "npc"), y = unit(0.5, "npc"),
#          width = unit(1, "npc"), height = unit(1, "npc"),

#  row.margin.vp <- viewport(
#      x=marginwidth * .5,
#      y=content.height * .5,
#      width=marginwidth,
#      height=content.height
#      );
#  pushViewport(row.margin.vp);

    if (draw.index) {
    grid.rect(x = unit(0.5, "npc"), y = unit(2.5, "mm") + unit(1, "lines"),
      width = width, height = unit(.9, "lines"), gp=gpar(fill="lightgray", lwd=0));
    }
    if (draw.names) {
      y <- unit(1, "lines") + 0.5 * stringWidth(v.names[max.name.length]);
      # unit(2.5, "mm") + 
      if (draw.index) {
        y <- y + unit(1, "lines");
      }    
    grid.rect(x = unit(0.5, "npc"), y = y,
      width = width, height = unit(.9, "lines"), gp=gpar(fill="lightgray", lwd=0));
    }

  for (i in 1:(len.vect)) {
    if (draw.index) {
      if (i == 1) { 
        grid.text(indexes[i],
            y=unit(2.5, "mm") + unit(1, "lines"),
            x=unit(1, "mm") + stringWidth(vect[i]) * .5,
            just=c("centre", "center")
            , gp=gpar(fill="lightgray")
            )
      } else {
        grid.text(indexes[i],
            y=unit(2.5, "mm") + unit(1, "lines"),
            x=unit(2, "mm") * (i-1) + unit(1, "mm") + sum(stringWidth(vect[1:(i-1)])) + stringWidth(vect[i]) * .5,
            just=c("centre", "center")
            , gp=gpar(fill="lightgray")
            )
      }
    }
    if (draw.names) {
      y <- unit(1, "lines") + 0.5 * stringWidth(v.names[max.name.length]);
      # unit(2.5, "mm") + 
      if (draw.index) {
        y <- y + unit(1, "lines");
      }    
      if (i == 1) { 
        grid.text(v.names[i],
            y=y,
            x=unit(1, "mm") + stringWidth(vect[i]) * .5,
            hjust="top",
            vjust="top",
           # just="right",
            rot=60
            );
      } else {
        grid.text(v.names[i],
            y=y,
            x=unit(2, "mm") * (i-1) + unit(1, "mm") + sum(stringWidth(vect[1:(i-1)])) + stringWidth(vect[i]) * .5,
            hjust="top",
            vjust="top",
           # just="right",
            rot=60
            );
      }
    }
  }
  popViewport();
  return(unit.c(height, width));
}

vectorBoxGrob <- function(v, x=.5, y=.5, draw.index=FALSE, draw.names=FALSE, a.factor=FALSE) { 

  if (!is.null(names(v))) {
    draw.index=TRUE;
    draw.names=TRUE;
  }

  height <- unit(1, "lines");
  if (draw.index) height <- height + unit(1, "lines");
  if (draw.names) {
    if (is.null(names(v))) stop("Cannot draw names if the vector has no names.");
    #max.name.length <- which.max(nchar(names(v)));
    #height <- height + stringWidth(names(v)[max.name.length]);
    height <- height + max(stringWidth(names(v)));
  }

  if (!a.factor) {
    if (is.character(v)) {
      n <- names(v);
      v <- paste("\"", v, "\"", sep="");
      names(v) <- n;
    }
  }
  width <- sum(stringWidth(v)) + unit(2, "mm") * length(v);

  grob(labels=v, x=x, y=y, width=width, height=height, draw.index=draw.index, draw.names=draw.names, cl="vectorBox");
}

drawDetails.vectorBox <- function(x, ...) {
  draw.vectorBox(x$labels, x$x, x$y, x$draw.index, x$draw.names, x$width, x$height);
}

xDetails.vectorBox <- function(x, theta) {
  grobX(roundrectGrob(x=x$x, y=x$y, width=x$width, height=x$height), theta);
}

yDetails.vectorBox <- function(x, theta) {
  grobY(rectGrob(x=x$x, y=x$y, width=x$width, height=x$height), theta);
}

grobWidth.vectorBox <- function(x) {
  return(x$width);
}

grobHeight.vectorBox <- function(x) {
  return(x$height);
}

############## ############## ############## ############## ##############
# Draw a matrix
############## ############## ############## ############## ##############

charactermatrix <- function(m) {
  if (is.character(m)) {
    d <- dim(m);
    dn <- dimnames(m);
    m <- paste("\"", m, "\"", sep="");
    dim(m) <- d;
    dimnames(m) <- dn;
  }
  return(m);
}

draw.matrixBox <- function(matrice, x=.5, y=.5, draw.index, draw.names, width, height, colwidth, marginwidth, marginheight) {

  #matrice <- charactermatrix(matrice);
  nr <- nrow(matrice);

  tablevp <- viewport(x=x, y=y, width=width, height=height);
  pushViewport(tablevp);

  # ------ the content of the matrix ---------

  content.width <- width - marginwidth;
  content.height <- height - marginheight;
  content.vp <- viewport(
      x=marginwidth + content.width * .5,
      y=content.height * .5,
      width=content.width,
      height=content.height
      );
  pushViewport(content.vp);

  grid.lines(x=unit(0, "npc"), y=unit(c(0, nr), "lines"));
  for (i in 1:ncol(matrice)) {
    grid.lines(x=sum(colwidth[1:i]) + unit(2, "mm") * i, y=unit(c(0, nr), "lines"));
  }

  grid.lines(y=unit(0, "lines"));
  for (i in 1:nrow(matrice)) {
    grid.lines(y=unit(i, "lines"));
  }

  for (i in 1:nrow(matrice)) {
    for (j in 1:ncol(matrice)) {
      grid.text(matrice[nrow(matrice) - (i-1), j],
          y=unit(i-1, "lines") + unit(2.5, "mm"),
          x=sum(colwidth[1:j]) + unit(2, "mm") * (j-1) + unit(2, "mm") * .5,
          just="right"
          );
    }
  }
  popViewport();

  # ------ the row margin: row names and row index

  row.margin.vp <- viewport(
      x=marginwidth * .5,
      y=content.height * .5,
      width=marginwidth,
      height=content.height
      );
  pushViewport(row.margin.vp);
  grid.rect(gp=gpar(fill="lightgray", lwd=0))
  if (draw.names) {
    for (i in 1:nrow(matrice)) {
      grid.text(rownames(matrice)[nrow(matrice) - (i-1)],
          y=unit(i-1, "lines") + unit(2.5, "mm"),
          x=unit(0, "npc"),
          just="left"
          );
    }
  }
  if (draw.index) {
    for (i in 1:nrow(matrice)) {
      grid.text(nrow(matrice) - i+1,
          y=unit(i-1, "lines") + unit(2.5, "mm"),
          x=marginwidth - unit(1, "mm"),
          just="right"
          );
    }
  }
  popViewport();

  # ------

  # ------ the col margin: col names and col index

  col.margin.vp <- viewport(
      x=marginwidth + content.width * .5,
      y=content.height + marginheight * .5,
      width=content.width,
      height=marginheight
      );
  pushViewport(col.margin.vp);
  grid.rect(gp=gpar(fill="lightgray", lwd=0))
  if (draw.names) {
    #y <- marginheight * .5;
    y <- marginheight ;
   # if (draw.index) {
   #   y <- (marginheight - unit(1, "lines")) * .5;
   # }
    for (j in 1:ncol(matrice)) {
      grid.text(colnames(matrice)[j],
          y=y,
          x=sum(colwidth[1:j]) + unit(2, "mm") * (j-1) + unit(2, "mm") * .5,
          #hjust="center",
          vjust="bottom",
          just="right",
          rot=60
          );
    }
  }
  if (draw.index) {
    for (j in 1:ncol(matrice)) {
      grid.text(j,
          y=unit(1, "lines") * .5,
          x=sum(colwidth[1:j]) + unit(2, "mm") * (j-1) + unit(2, "mm") * .5,
          just="right"
          );
    }
  }
  popViewport();

  # ------

  popViewport();
}

matrixBoxGrob <- function(m, x=.5, y=.5, draw.index=FALSE, draw.names=FALSE) { 

  m <- charactermatrix(m);
  if (!is.null(rownames(m)) & !is.null(colnames(m))) {
    draw.index=TRUE;
    draw.names=TRUE;
  }

  lchar <- nchar(m);
  ilonguest <- apply(nchar(m), 2, which.max);
  ncell <- nrow(m) * ncol(m);
  longest <- m[seq(0, ncell-1, nrow(m)) + ilonguest]
  colwidth <- stringWidth(longest);

  width <- sum(colwidth) + unit(2, "mm") * ncol(m);
  height <- unit(1, "lines") * nrow(m);

  marginwidth <- unit(0, "mm");
  if (draw.index) {
    marginwidth <- marginwidth + stringWidth(nrow(m)) + unit(1, "mm");
  }
  if (draw.names) {
    if (is.null(colnames(m)) | is.null(rownames(m))) {
      stop("Cannot draw names if the matrix has no row names.");
    }
    imax <- which.max(nchar(rownames(m)));
    marginwidth <- marginwidth + stringWidth(rownames(m)[imax]) + unit(1, "mm");
  }
  width <- width + marginwidth

  marginheight <- unit(0, "mm");
  if (draw.index) {
    marginheight <- marginheight + unit(1, "lines");
  }
  if (draw.names) {
    if (is.null(colnames(m))) {
      stop("Cannot draw names if the matrix has no column names.");
    }
    imax <- which.max(nchar(colnames(m)));
    max.name.length <- stringWidth(colnames(m)[imax]);
    marginheight <- marginheight + max.name.length;
  }
  height <- height + marginheight;

  grob(labels=m, x=x, y=y, draw.index=draw.index, draw.names=draw.names, width=width, height=height, cl="matrixBox", colwidth=colwidth, marginwidth=marginwidth, marginheight=marginheight);
}

drawDetails.matrixBox <- function(x, ...) {
  draw.matrixBox(x$labels, x$x, x$y, draw.index=x$draw.index, draw.names=x$draw.names, width=x$width, height=x$height, colwidth=x$colwidth, marginwidth=x$marginwidth, marginheight=x$marginheight);
}

xDetails.matrixBox <- function(x, theta) {
  grobX(roundrectGrob(x=x$x, y=x$y, width=x$width, height=x$height), theta);
}

yDetails.matrixBox <- function(x, theta) {
  grobY(rectGrob(x=x$x, y=x$y, width=x$width, height=x$height), theta);
}

grobWidth.matrixBox <- function(x) {
  x$width
}

grobHeight.matrixBox <- function(x) {
  x$height;
}

############## ############## ############## ############## ##############
# Draw a data frame
############## ############## ############## ############## ##############

draw.dataframeBox <- function(matrice, x=.5, y=.5, draw.index, draw.names, width, height, colwidth, marginwidth, marginheight) {

  nr <- nrow(matrice);

  tablevp <- viewport(x=x, y=y, width=width, height=height);
  pushViewport(tablevp);

  # ------ the content of the matrix ---------

  content.width <- width - marginwidth;
  content.height <- height - marginheight;
  content.vp <- viewport(
      x=marginwidth + content.width * .5,
      y=content.height * .5,
      width=content.width,
      height=content.height
      );
  pushViewport(content.vp);

  # The dashed lines list-like box

  grid.lines(x=unit(0, "npc"), y=unit.c(unit(c(0), "lines"), unit(c(nr), "lines") + unit(2, "mm")), gp=gpar(lty=2));
  for (i in 1:ncol(matrice)) {
    grid.lines(x=sum(colwidth[1:i]) + unit(4, "mm") * i, y=unit.c(unit(c(0), "lines"), unit(c(nr), "lines") + unit(2, "mm")), gp=gpar(lty=2));
  }

  grid.lines(y=unit(0, "lines"), gp=gpar(lty=2));
  grid.lines(y=unit(nr, "lines") + unit(2, "mm"), gp=gpar(lty=2));

  # the inner solid line box

  for (i in 1:ncol(matrice)) {
    if (i == 1) {
## lines on both sides
      x1 <- unit(0, "npc") + unit(1, "mm");
      x2 <- sum(colwidth[1:i]) + unit(3, "mm"); 
      grid.lines(x=x1, y=unit.c(unit(c(0), "lines") + unit(1, "mm"), unit(c(nr), "lines") + unit(1, "mm")), gp=gpar(lty=1));
      grid.lines(x=x2, y=unit.c(unit(c(0), "lines") + unit(1, "mm"), unit(c(nr), "lines") + unit(1, "mm")), gp=gpar(lty=1));
## lines on all row

      grid.lines(y=unit(0, "lines") + unit(1, "mm"), x=unit.c(x1, x2));
      for (i in 1:nrow(matrice)) {
        grid.lines(
            y=unit(i, "lines") + unit(1, "mm"),
            x=unit.c(x1, x2),
            gp=gpar(lty=1)
            );
      }
    } else {
## lines on both sides
      x1 <- sum(colwidth[1:i]) - colwidth[i] + unit(4, "mm") * (i-1) + unit(1, "mm"); 
      x2 <- sum(colwidth[1:i]) + unit(4, "mm") * (i-1) + unit(3, "mm");
      grid.lines(x=x1, y=unit.c(unit(c(0), "lines") + unit(1, "mm"), unit(c(nr), "lines") + unit(1, "mm")), gp=gpar(lty=1));
      grid.lines(x=x2 , y=unit.c(unit(c(0), "lines") + unit(1, "mm"), unit(c(nr), "lines") + unit(1, "mm")), gp=gpar(lty=1));
## lines on all row

      grid.lines(y=unit(0, "lines") + unit(1, "mm"), x=unit.c(x1, x2));
      for (i in 1:nrow(matrice)) {
        grid.lines(
            y=unit(i, "lines") + unit(1, "mm"),
            x=unit.c(x1, x2),
            gp=gpar(lty=1)
            );
      }
    }
  }

  for (i in 1:nrow(matrice)) {
    for (j in 1:ncol(matrice)) {
      grid.text(matrice[nrow(matrice) - (i-1), j],
          y=unit(i-1, "lines") + unit(3.5, "mm"),
          x=sum(colwidth[1:j]) + unit(4, "mm") * (j-1) + unit(4, "mm") * .5,
          just="right"
          );
    }
  }
  popViewport();

  # ------ the row margin: row names and row index

  row.margin.vp <- viewport(
      x=marginwidth * .5,
      y=content.height * .5,
      width=marginwidth,
      height=content.height
      );
  pushViewport(row.margin.vp);
  grid.rect(gp=gpar(fill="lightgray", lwd=0))
  if (draw.names) {
    for (i in 1:nrow(matrice)) {
      grid.text(rownames(matrice)[nrow(matrice) - (i-1)],
          y=unit(i-1, "lines") + unit(2.5, "mm"),
          x=unit(0, "npc"),
          just="left"
          );
    }
  }
  if (draw.index) {
    for (i in 1:nrow(matrice)) {
      grid.text(nrow(matrice) - i+1,
          y=unit(i-1, "lines") + unit(2.5, "mm"),
          x=marginwidth - unit(1, "mm"),
          just="right"
          );
    }
  }
  popViewport();

  # ------

  # ------ the col margin: col names and col index

  col.margin.vp <- viewport(
      x=marginwidth + content.width * .5,
      y=content.height + marginheight * .5,
      width=content.width,
      height=marginheight
      );
  pushViewport(col.margin.vp);
  grid.rect(gp=gpar(fill="lightgray", lwd=0))
  if (draw.names) {
    #y <- marginheight * .5;
    y <- marginheight ;
   # if (draw.index) {
   #   y <- (marginheight - unit(1, "lines")) * .5;
   # }
    for (j in 1:ncol(matrice)) {
      grid.text(colnames(matrice)[j],
          y=y,
          x=sum(colwidth[1:j]) + unit(2, "mm") * (j-1) + unit(2, "mm") * .5,
          #hjust="center",
          vjust="bottom",
          just="right",
          rot=60
          );
    }
  }
  if (draw.index) {
    for (j in 1:ncol(matrice)) {
      grid.text(j,
          y=unit(1, "lines") * .5,
          x=sum(colwidth[1:j]) + unit(2, "mm") * (j-1) + unit(2, "mm") * .5,
          just="right"
          );
    }
  }
  popViewport();

  # ------

  popViewport();
}

dataframeBoxGrob <- function(m, x=.5, y=.5, draw.index=FALSE, draw.names=FALSE) { 

  m <- charactermatrix(m);
  if (!is.null(rownames(m)) & !is.null(colnames(m))) {
    draw.index=TRUE;
    draw.names=TRUE;
  }

  asm <- as.matrix(m); # as.matrix since it will not work with a data.frame
  lchar <- nchar(m);
  ilonguest <- apply(nchar(asm), 2, which.max);
  ncell <- nrow(m) * ncol(m);
  longest <- asm[seq(0, ncell-1, nrow(m)) + ilonguest]
  colwidth <- stringWidth(longest);

  width <- sum(colwidth) + unit(4, "mm") * ncol(m); # 2mm added for the outer list    # DATA FRAME
  height <- unit(1, "lines") * nrow(m) + unit(2, "mm"); # added for the outer list    # DATA FRAME

  marginwidth <- unit(0, "mm");
  if (draw.index) {
    marginwidth <- marginwidth + stringWidth(nrow(m)) + unit(1, "mm");
  }
  if (draw.names) {
    if (is.null(colnames(m)) | is.null(rownames(m))) {
      stop("Cannot draw names if the matrix has no row names.");
    }
    imax <- which.max(nchar(rownames(m)));
    marginwidth <- marginwidth + stringWidth(rownames(m)[imax]) + unit(1, "mm");
  }
  width <- width + marginwidth

  marginheight <- unit(0, "mm");
  if (draw.index) {
    marginheight <- marginheight + unit(1, "lines");
  }
  if (draw.names) {
    if (is.null(colnames(m))) {
      stop("Cannot draw names if the vector has no column names.");
    }
    imax <- which.max(nchar(colnames(m)));
    max.name.length <- stringWidth(colnames(m)[imax]);
    marginheight <- marginheight + max.name.length;
  }
  height <- height + marginheight;

  grob(labels=m, x=x, y=y, draw.index=draw.index, draw.names=draw.names, width=width, height=height, cl="dataframeBox", colwidth=colwidth, marginwidth=marginwidth, marginheight=marginheight);
}

drawDetails.dataframeBox <- function(x, ...) {
  draw.dataframeBox(x$labels, x$x, x$y, draw.index=x$draw.index, draw.names=x$draw.names, width=x$width, height=x$height, colwidth=x$colwidth, marginwidth=x$marginwidth, marginheight=x$marginheight);
}

xDetails.dataframeBox <- function(x, theta) {
  grobX(roundrectGrob(x=x$x, y=x$y, width=x$width, height=x$height), theta);
}

yDetails.dataframeBox <- function(x, theta) {
  grobY(rectGrob(x=x$x, y=x$y, width=x$width, height=x$height), theta);
}

grobWidth.dataframeBox <- function(x) {
  x$width
}

grobHeight.dataframeBox <- function(x) {
  x$height;
}

