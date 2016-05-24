editLattice <-
function (formLatticeOutput) 
{
  #
  #  This function is a highly modified version of 'edit.nb' in 
  #  package 'spdep' originally by Roger Bivand.
  #
    latt <- formLatticeOutput$latt
    nodes <- formLatticeOutput$nodes
    polys <- formLatticeOutput$poly
    hole.list <- formLatticeOutput$hole.list
    clatt <- card(latt)
    if (!inherits(latt, "nb")) 
        stop("not a neighbours list")
    cl <- class(latt)
    if (length(cl) > 1) 
        icl <- cl[-match("nb", cl)]
    else icl <- NULL
    x <- nodes[, 1]
    y <- nodes[, 2]
    n <- length(latt)
    row.names <- attr(latt, "region.id")
    if (is.null(row.names)) 
        row.names <- as.character(1:n)
    xlim <- range(x)
    ylim <- range(y)
    plot.new()
    plot.window(xlim = xlim, ylim = ylim, "", asp = 1)
    lines(polys,col=6,lwd=0.5)
    n.holes = length(hole.list)
    for(i in 1:n.holes){
      lines(rbind(hole.list[[i]],hole.list[[i]][1,]),col=6,lwd=0.5)
      }
    for (i in 1:n) {
        if (clatt[i] > 0) 
            segments(x[i], y[i], x[latt[[i]]], y[latt[[i]]])
    }
    if(min(card(latt))==0){points(nodes[card(latt)==0,], col=2, cex=1,pch=19)}
    finished <- "n"
    deletions <- NULL
    additions <- NULL
   erase.col <- par()$bg
    while (finished == "n") {
        cat("Identifying contiguity for deletion ...\n")
        cand <- identify(x, y, n = 2,plot=FALSE)
        lines(x[cand], y[cand], col = "red",lwd=2)
        if (.Platform$OS.type == "windows") 
            bringToTop(-1)
        if ((cand[2] %in% latt[[cand[1]]]) && (cand[1] %in%     
          latt[[cand[2]]])) {
            delete <- readline("Delete this line (y/n) ")
            if (delete != "y") 
                delete <- "n"
            else {
           #deletions <- c(deletions, paste(cand, collapse = "-"))
                          latt[[cand[1]]] <- latt[[cand[1]]][latt[[cand[1]]] != cand[2]]
                          if (length(latt[[cand[1]]]) == 0) {
                            latt[[cand[1]]] <- as.integer(0)
                            cat(cand[1], "is now an island\n")
                          }
                          latt[[cand[2]]] <- latt[[cand[2]]][latt[[cand[2]]] != 
                            cand[1]]
                          if (length(latt[[cand[2]]]) == 0) {
                            latt[[cand[2]]] <- as.integer(0)
                            cat(cand[2], "is now an island\n")
                          }
                          lines(x[cand], y[cand], col = erase.col)
                cat("deleted contiguity between point", cand[1], 
                  "and", cand[2], "\n")
    plot.new()
    plot.window(xlim = xlim, ylim = ylim, "", asp = 1)
    lines(polys,col=6,lwd=0.5)
    for (i in 1:n) {
        if (clatt[i] > 0) 
            segments(x[i], y[i], x[latt[[i]]], y[latt[[i]]])
    }
    if(min(card(latt))==0){points(nodes[card(latt)==0,], col=2, cex=1,pch=19)}

            }
        }
        else {
            if (length(cand == 2)) {
                cat("No contiguity between chosen points\n")
                addcont <- readline("Add contiguity? (y/n) ")
                if (addcont != "y") 
                  addcont <- "n"
                if (addcont == "y") {
                  latt[[cand[1]]] <- sort(unique(c(latt[[cand[1]]], 
                    cand[2])))
                  latt[[cand[2]]] <- sort(unique(c(latt[[cand[2]]], 
                    cand[1])))
                  cat("added contiguity between point", cand[1], 
                    "and", cand[2], "\n")
                  #additions <- c(additions, paste(cand, collapse = "-"))
                  lines(x[cand], y[cand], col = "black")
                  latt[[cand[2]]] = latt[[cand[2]]][ latt[[cand[2]]]>0]
                  latt[[cand[1]]] = latt[[cand[1]]][ latt[[cand[1]]]>0]
                }
            }
        }
        finished <- readline("Options: quit[q] continue[c] ")
        if (finished != "q") 
            finished <- "n"
    }
    attr(latt, "region.id") <- row.names
    if (is.null(icl)) 
        class(latt) <- "nb"
    else class(latt) <- c("nb", icl)
    formLatticeOutput$latt <- latt
    return(formLatticeOutput)
}

