
# Create a list straight from the RGML file
readPicture <- function(rgmlFile, tidy=TRUE, ...) {
    funGetX = function(x) { as.numeric(xmlAttrs(x)["x"]) }
    funGetY = function(x) { as.numeric(xmlAttrs(x)["y"]) }
    funGetGC = function(x) { pars <- xmlApply(x, xmlAttrs) }
    funPathType = function(x) { xmlAttrs(x)["type"] }
    funToLineJoin = function(x) {
        x <- as.numeric(x)
        if (x == 0)
            1 # PS 0 is miter, R is round
        else if (x == 1)
            0 # PS 1 is round, R is miter
        else
            x # should be 2, bevel
    }
    funToLineEnd = function(x) {
        x <- as.numeric(x)
        if (x == 0)
            1 # PS 0 is butt, R is round
        else if (x == 1)
            0 # PS 1 is round, R is butt
        else
            x # should be 2, square
    }

    funPath = function(x) {
        switch(xmlName(x),
               path={ xval = unlist(xmlApply(x, funGetX))
                      yval = unlist(xmlApply(x, funGetY))

                      # get context
                      gc = funGetGC(xmlElementsByTagName(x, "context")[[1]])

                      # convert to a list of S4 objects
                      switch(funPathType(x),
                             stroke=new("PictureStroke", x=xval, y=yval,
                               lwd=as.numeric(gc$style["lwd"]),
                               lty=readLTY(gc$style["lty"]),
                               lineend=funToLineEnd(gc$style["lineend"]),
                               linejoin=funToLineJoin(gc$style["linejoin"]),
                               linemitre=as.numeric(gc$style["linemitre"]),
                               rgb=rgb(as.numeric(gc$rgb["r"]),
                                 as.numeric(gc$rgb["g"]),
                                 as.numeric(gc$rgb["b"]))),
                             fill=new("PictureFill", x=xval, y=yval,
                               rule="nonzero",
                               lwd=as.numeric(gc$style["lwd"]),
                               lty=readLTY(gc$style["lty"]),
                               lineend=funToLineEnd(gc$style["lineend"]),
                               linejoin=funToLineJoin(gc$style["linejoin"]),
                               linemitre=as.numeric(gc$style["linemitre"]),
                               rgb=rgb(as.numeric(gc$rgb["r"]),
                                 as.numeric(gc$rgb["g"]),
                                 as.numeric(gc$rgb["b"]))),
                             eofill=new("PictureFill", x=xval, y=yval,
                               rule="evenodd",
                               lwd=as.numeric(gc$style["lwd"]),
                               lty=readLTY(gc$style["lty"]),
                               lineend=funToLineEnd(gc$style["lineend"]),
                               linejoin=funToLineJoin(gc$style["linejoin"]),
                               linemitre=as.numeric(gc$style["linemitre"]),
                               rgb=rgb(as.numeric(gc$rgb["r"]),
                                 as.numeric(gc$rgb["g"]),
                                 as.numeric(gc$rgb["b"]))),
                             char=new("PictureChar", x=xval, y=yval,
                               char=xmlAttrs(x)["char"],
                               lwd=as.numeric(gc$style["lwd"]),
                               lty=readLTY(gc$style["lty"]),
                               lineend=funToLineEnd(gc$style["lineend"]),
                               linejoin=funToLineJoin(gc$style["linejoin"]),
                               linemitre=as.numeric(gc$style["linemitre"]),
                               rgb=rgb(as.numeric(gc$rgb["r"]),
                                 as.numeric(gc$rgb["g"]),
                                 as.numeric(gc$rgb["b"]))))
                  },
               text={ xa <- xmlAttrs(x)
                   
                      # get context
                      gc = funGetGC(xmlElementsByTagName(x, "context")[[1]])

                      type = xa["type"]

                      if (type == "text") {
                          letters = list()
                      } else {
                          # list of PictureChar or PictureText
                          letters = xmlApply(x, funPath)
                          cntxt = which(names(letters) == "context")
                          letters = letters[-cntxt]
                      }
                      
                      new("PictureText",
                          string=
                            if (tidy)
                              tidyString(xa["string"])
                            else
                              xa["string"],
                          x=as.numeric(xa["x"]),
                          y=as.numeric(xa["y"]),
                          w=as.numeric(xa["width"]),
                          h=as.numeric(xa["height"]),
                          bbox=as.numeric(strsplit(xa["bbox"], " ")[[1]]),
                          angle=as.numeric(xa["angle"]),
                          letters=letters,
                          lwd=as.numeric(gc$style["lwd"]),
                          rgb=rgb(as.numeric(gc$rgb["r"]),
                            as.numeric(gc$rgb["g"]),
                            as.numeric(gc$rgb["b"])))
                  },
               summary={ attrs <- xmlAttrs(x)
                         numattrs <- as.numeric(xmlAttrs(x)) 
                         names(numattrs) <- names(attrs) 
                         numattrs })
    }

    xmlDoc = xmlTreeParse(rgmlFile, ...)
    version <- as.numeric(xmlAttrs(xmlRoot(xmlDoc))["version"])
    if (version != 3)
        stop(paste("Version mismatch:",
                   "RGML file needs to be recreated with PostScriptTrace()"))
    RGMLlist <- xmlApply(xmlRoot(xmlDoc), funPath)
    new("Picture",
        paths=RGMLlist[-length(RGMLlist)],
        summary=new("PictureSummary",
          numPaths=length(RGMLlist) - 1, # RGMLlist$summary["count"],
          xscale=RGMLlist$summary[c("xmin", "xmax")],
          yscale=RGMLlist$summary[c("ymin", "ymax")]))
}

# Given a list of paths, determine the bounding box
pathBounds <- function(paths) {
    pathXmin <- function(path) {
        if (is(path, "PictureText"))
            path@bbox[1]
        else
            min(path@x)
    }
    pathYmin <- function(path) {
        if (is(path, "PictureText"))
            path@bbox[2]
        else 
            min(path@y)
    }
    pathXmax <- function(path) {
        if (is(path, "PictureText"))
            path@bbox[3]
        else
            max(path@x)
    }
    pathYmax <- function(path) {
        if (is(path, "PictureText"))
            path@bbox[4]
        else
            max(path@y)
    }
    list(xscale=c(min(sapply(paths, pathXmin)),
           max(sapply(paths, pathXmax))),
         yscale=c(min(sapply(paths, pathYmin)), 
           max(sapply(paths, pathYmax))))
}

setMethod("[", "Picture",
          function(x, i, j, drop) {
              paths <- x@paths[i]
              scales <- pathBounds(paths)
              new("Picture",
                  paths=paths,
                  summary=new("PictureSummary",
                    numPaths=length(paths),
                    xscale=scales$xscale,
                    yscale=scales$yscale))
          })

setMethod("[[", "Picture",
          function(x, i, j, drop) {
              if (length(i) > 1)
                  stop("Index must be length 1")
              path <- x@paths[[i]]
              if (is(path, "PictureText") &&
                  length(path@letters) > 0) {
                  paths <- path@letters
                  scales <- pathBounds(paths)
                  new("Picture",
                      paths=paths,
                      summary=new("PictureSummary",
                        numPaths=length(paths),
                        xscale=scales$xscale,
                        yscale=scales$yscale))
              } else {
                  x[i]
              }
          })



