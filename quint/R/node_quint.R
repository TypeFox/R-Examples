node_quint <-
structure(function(obj, id = TRUE, abbreviate = FALSE, fill = "white", gp = gpar())
{
  meta <- obj$data
  nam <- names(obj)

  extract_label <- function(node) {
    if(is.terminal(node)) return(rep.int("", 2))

    varlab <- character_split(split_node(node), meta)$name
    if(abbreviate > 0) varlab <- abbreviate(varlab, as.numeric(abbreviate))

    plab <- ""
    return(c(varlab, plab))
  }

  maxstr <- function(node) {
      lab <- extract_label(node)
      klab <- if(is.terminal(node)) "" else unlist(lapply(kids_node(node), maxstr))
      lab <- c(lab, klab)
      lab <- unlist(lapply(lab, function(x) strsplit(x, "\n")))
      return(lab[which.max(nchar(lab))])
  }

  nstr <- maxstr(node_party(obj))
  if(nchar(nstr) < 6) nstr <- "aAAAAa"

  ### panel function for the inner nodes
  rval <- function(node) {  
    node_vp <- viewport(
      x = unit(0.5, "npc"),
      y = unit(0.5, "npc"),
      width = unit(1, "strwidth", nstr) * 1.3, 
      height = unit(3, "lines"),
      name = paste("node_inner", id_node(node), sep = ""),
      gp = gp
    )
    pushViewport(node_vp)

    xell <- c(seq(0, 0.2, by = 0.01),
  	      seq(0.2, 0.8, by = 0.05),
  	      seq(0.8, 1, by = 0.01))
    yell <- sqrt(xell * (1-xell))

    lab <- extract_label(node)
    fill <- rep(fill, length.out = 2)

    grid.polygon(x = unit(c(xell, rev(xell)), "npc"),
        	 y = unit(c(yell, -yell)+0.5, "npc"),
        	 gp = gpar(fill = fill[1]))

    grid.text(lab[1], y = unit(1.5 + 0.5 * FALSE, "lines"))

    if(id) {
      nodeIDvp <- viewport(x = unit(0.5, "npc"), y = unit(1, "npc"),
        width = max(unit(1, "lines"), unit(1.3, "strwidth", nam[id_node(node)])),
        height = max(unit(1, "lines"), unit(1.3, "strheight", nam[id_node(node)])))
      pushViewport(nodeIDvp)
      popViewport()
    }
    upViewport()
  }
  return(rval)
}, class = "grapcon_generator")
