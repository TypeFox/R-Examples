# Copyright (c) 2015 Santiago Barreda
# All rights reserved.


errorbars = function(x, y, top, bottom = top, length = .2, add = TRUE, ...){
  if (add) arrows(x, y+top, x, y-bottom, angle=90, code=3, length=length, ...)
  if (!add){
    plot (x,y,pch=16, ylim = range(y) + c(-top, bottom))
    arrows(x, y+top, x, y-bottom, angle=90, code=3, length=length, ...)
  }
}


errorbar = function(x, y, top, bottom = top, length = .2, add = TRUE, ...){
  cl = match.call()
  args = sapply (2:length(cl), function(x) cl[[x]])
  names(args) = names(cl)[-1]
  do.call (errorbars, args)
}

