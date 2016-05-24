terminal_quint <-
structure(function(obj,
                      col = "black",
                      fill = "lightgray",
                      width = 0.5,
                      yscale = NULL,
                      ylines = 3,
                      cex = 0.5,
                      id = TRUE,
                      gp = gpar())
  {
    ni <- obj$ni
    ni.temp <- ni
    d.CI <- data.frame(
                       lower = ni[,8] - 1.96 * ni[,9],
                       mean = ni[,8],
                       upper = ni[,8] + 1.96 * ni[,9]
                       )
    yscale <- c(0-max(abs(d.CI)),max(abs(d.CI))) + c(-0.1, 0.1) * max(abs(d.CI))
    rval <- function(node) { # core plotting function
      nid <- id_node(node)
      top_vp <- viewport(layout = grid.layout(nrow = 3, ncol = 3, # define viewport
                           widths = unit(c(ylines, 1, 1),
                             c("lines", "null", "lines")),  
                           heights = unit(c(1, 1, 2), c("lines","null","lines"))),
                         width = unit(1, "npc"), 
                         height = unit(1, "npc") - unit(2, "lines"),
                         name = paste("node_quint", nid, sep = ""),
                         gp = gp)

      pushViewport(top_vp)
      grid.rect(gp = gpar(fill = "transparent", col =0))

      ind2 <- nid == nodeids(obj,terminal=TRUE)
      
      ## main title
      bottom <- viewport(layout.pos.col=2, layout.pos.row=3)
      pushViewport(bottom)
      if(id){
        grid.text(sprintf("Leaf %s\nP%i",
                          which(nid==nodeids(obj,terminal=TRUE)), ni[ind2,10])
                  )
      }else{grid.text("")}
      popViewport()
      
      plot <- viewport(layout.pos.col = 2, layout.pos.row = 2,
                       xscale = c(0, 1), yscale = yscale,
                       name = paste("node_quint", nid, "plot", 
                         sep = ""))
      pushViewport(plot)
      
      xl <- 0.5 - width/4
      xr <- 0.5 + width/4

      ## box & whiskers
      grid.rect(gp = gpar(fill = c("#92D050","#FF572F","#DDD8C2")[ni[ind2,10]]))
      ##refline
      grid.lines(unit(c(0, 1), "npc"), 
                 unit(0, "native"), gp = gpar(col = col,lwd=unit(width/4,"npc"),lty="dashed"))
      grid.lines(unit(c(xl, xr), "npc"), 
                 unit(d.CI$lower[ind2], "native"), gp = gpar(col = col,lwd=unit(width,"npc")))
      grid.lines(unit(0.5, "npc"), 
                 unit(c(d.CI$lower[ind2],d.CI$upper[ind2]), "native"), gp = gpar(col = col,lwd=unit(width/2,"npc")))
      meanline=FALSE ## if FALSE mean point
      if(meanline){
      grid.lines(unit(c(0.5 - width/2, 0.5+width/2), "npc"), 
                 unit(d.CI$mean[ind2], "native"), gp = gpar(col = col, lwd = 2))
    }else{      grid.points(unit(0.5, "npc"), 
                 unit(d.CI$mean[ind2], "native"), gp = gpar(col = col,lwd=1),size=unit(1,"lines"),pch=20)}
      grid.lines(unit(c(xl, xr), "npc"), unit(d.CI$upper[ind2], "native"), 
                 gp = gpar(col = col,lwd=unit(width,"npc")))
      
      grid.yaxis(label=TRUE) ## TO DO only TRUE for left terminal
      grid.rect(gp = gpar(fill = "transparent"))
      upViewport(2)
    }
    return(rval)
  }, class = "grapcon_generator")
