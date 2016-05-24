rysgran.ternary <-
  function (x = NULL, method = "shepard", lang= "en-US", main = NULL,
            z = NULL, show.labels = FALSE, label.points = FALSE, labels = NULL, axis.labels = NULL, 
            show.names = TRUE, show.lines = TRUE, show.legend = TRUE, show.grid = FALSE, 
            z.cex.range = NULL, cex.labels = 1, cex.points = 1, cex.axis = 1, cex.names = 0.8,
            col.names = "gray2", col = "black", col.labels = "black", col.axis = "black", 
            col.lines = "black", col.grid = "gray", pos = 1, pch = 19,
            lty.grid = 3, ...)
    
  {
    if (is.null(x))
      stop ("x is missing")
    if (!is.matrix(x) && !is.data.frame(x)) 
      stop("x must be a matrix or data frame with at least 3 columns and one row.")
    if (method!="shepard" && method!="pejrup" && method!="flemming")
      stop("methods supported are shepard, pejrup and flemming.")
    if (any(x > 1) || any(x < 0))
    {
      if (any(x < 0)) 
        stop("All proportions must be between zero and one.")
      if (any(x > 100)) 
        stop("All percentages must be between zero and 100.")
      x <- x/100
    }
    if (any(abs(rowSums(x) - 1) > 0.01)) 
      warning("At least one set of proportions does not equal 100%.")
    
    par(xpd = TRUE)
    
    if (is.null(main))
    {
      if (method=="shepard")
      {
        if (lang=="en-US" | lang=="en-GR"| lang=="eng"| lang=="e")
          main <- "Shepard Diagram"
        if (lang=="pt-BR" | lang=="pt-PT"| lang=="port"| lang=="p")
          main <- "Diagrama de Shepard"
      }
      if (method=="pejrup")
      {
        if (lang=="en-US" | lang=="en-GR"| lang=="eng"| lang=="e")
          main <- "Pejrup Diagram"
        if (lang=="pt-BR" | lang=="pt-PT"| lang=="port"| lang=="p")
          main <- "Diagrama de Pejrup"
      }
      if (method=="flemming")
      {
        if (lang=="en-US" | lang=="en-GR"| lang=="eng"| lang=="e")
          main <- "Flemming Diagram"
        if (lang=="pt-BR" | lang=="pt-PT"| lang=="port"| lang=="p")
          main <- "Diagrama de Flemming"
      }
    }
    
    if (label.points | show.labels)
    {
      if (is.null(labels)) labels<-row.names(x)
    }
    
    plot(0.5, type = "n", axes = FALSE, xlim = c(0, 1), ylim = c(0,1), main = main, xlab = "", ylab = "")
    
    if (method=="shepard")
    {
      shepard.plot(x = x, lang = lang, axis.labels = axis.labels, show.names = show.names, 
                   show.lines = show.lines, show.legend = show.legend,show.grid = show.grid, 
                   cex.axis = cex.axis, cex.names = cex.names, col.labels= col.labels,
                   col.axis = col.axis, col.names = col.names, col.lines = col.lines, col.grid = col.grid, 
                   lty.grid = par("lty"))
    }
    
    if (method=="pejrup")
    { 
      pejrup.plot(x = x, lang = lang, axis.labels = axis.labels, show.names = show.names, 
                  show.lines = show.lines, show.legend = show.legend,show.grid = show.grid, 
                  cex.axis = cex.axis, cex.names = cex.names, col.labels= col.labels,
                  col.axis = col.axis, col.names = col.names, col.lines = col.lines, col.grid = col.grid, 
                  lty.grid = par("lty"))
    }
    
    if (method=="flemming")
    {
      flemming.plot(x = x, lang = lang, axis.labels = axis.labels, show.names = show.names, 
                    show.lines = show.lines, show.legend = show.legend, show.grid = show.grid, 
                    cex.axis = cex.axis, cex.names = cex.names, col.labels= col.labels,
                    col.axis = col.axis, col.names = col.names, col.lines = col.lines, col.grid = col.grid, 
                    lty.grid = par("lty"))
      if (show.legend==TRUE){
        if(show.names==FALSE){
          warning ("if is TRUE show.legend show.names should also be TRUE")
        }
        tab<- flemming.table(lang=lang)
        print(tab)
        }
      }
    
    par(xpd = FALSE)
        
    ypos <- x[, 3] * (sin(pi/3)) 
    xpos <- 1 - (x[, 1] + x[, 3] * 0.5) 
        
    if (show.labels)
    {
      label.points=FALSE
      z=NULL
      points(xpos, ypos, col = col, pch = NA, cex= cex.points,type = "n",...)
      text(xpos, ypos, labels=labels, pos=NULL, col=col.labels, cex = cex.labels)
    }
    
    if (label.points & is.null(z))
    {
      points(xpos, ypos, col = col, pch = pch, cex= cex.points,type = "p",...)
      text(xpos, ypos, labels=labels, pos=pos, col=col.labels, cex = cex.labels)
    }
    
    if(!is.null(z))
    {
      if(nrow(x) != length(z)) stop("z must be a vector with the same length of x")
      if(is.null(z.cex.range)) z.cex.range<-c(1,3)
      cex.bubbles <- TT.str(z, z.cex.range[1], z.cex.range[2])
      points(xpos, ypos, pch = pch, col = col, type = "p", cex = cex.bubbles,...)
      text(xpos, ypos, labels=labels, pos=pos, col=col.labels, cex = cex.labels)
    }
    
    if(label.points==FALSE & show.labels==FALSE & is.null(z))
      points(xpos, ypos, col = col, pch = pch, cex= cex.points,type = "p",...)
    }
