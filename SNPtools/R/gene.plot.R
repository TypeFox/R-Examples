################################################################################
# Given a genomic range, plot the genes in the interval.
# Daniel Gatti
# Dan.Gatti@jax.org
# Feb. 13, 2013
################################################################################
# Arguments: mgi.file: data.frame, containing genes (or other features) as 
#                      returned by get.mgi.features().
#            col: color, the color of the gene rectangles.
#            ...: other arguments to be passed to plot.
# Returns: data.frame with gene symbol, gene start, gene end, text start, 
#          text end and row.
gene.plot = function(mgi, col = "black", ...) {

  previous.cex = par("cex")

  call = match.call(expand.dots = T)

  # Get the Chr, start and end.
  chr   = mgi$seqid[1]
  start = min(mgi$start)
  end   = max(mgi$stop)

  # Line the genes up sequentially in columns.
  # Locs holds the gene symbol, the gene start and end, the text start and end,
  # as well as the row to plot on.
  locs = data.frame(name = mgi$Name, gstart = mgi$start * 1e-6,
         gend = mgi$stop * 1e-6, tstart = rep(0, nrow(mgi)),
         tend = rep(0, nrow(mgi)), row = 1:nrow(mgi))
  par(lend = 2)

  m = which(names(call) == "xlim")
  if(length(m) > 0) {
    plot(0, 0, col = 0, xlab = paste("Chr", chr, "(Mb)"), ylab = "",
         yaxt = "n", ...)
  } else {
    plot(0, 0, col = 0,   xlim = round(c(min(mgi$start), max(mgi$stop)) * 1e-6),
         xlab = paste("Chr", chr, "(Mb)"), ylab = "", yaxt = "n", ...)
  } # else

  old.cex = 100
  iter = 0
  while(abs(par("cex") - old.cex) > 0.1 & iter < 9) {

    # Offset between gene end and text start.
    offset = strwidth("i")
    locs$tstart = locs$gend   + offset
    locs$tend   = locs$tstart + strwidth(mgi$Name)

    locs$row = 1:nrow(locs)
    locs = line.up.genes(locs)
    locs = resolve.collisions(locs)

    # Determine the number of rows and the row height.
    max.row = max(locs$row) + 1
    usr = par("usr") 
    rowht = (usr[4] - usr[3]) / max.row
    old.cex = par("cex")
    par(cex = old.cex * rowht / strheight("W"))
    iter = iter + 1
  } # for(i)

  par(cex = 0.9 * par("cex"))

  # Plot the genes.
  rect(xleft = locs$gstart, ybottom = usr[4] - rowht * (locs$row + 1) + 0.05 * rowht, 
       xright = locs$gend, ytop = usr[4] - rowht * locs$row - 0.05 * rowht,
       density = -1, col = col, border = col)
  text(x = locs$tstart, y = usr[4] - rowht * locs$row, 
       labels = locs$name, adj = c(0,1), col = col)

  par(cex = previous.cex)

  return(locs)

} # gene.plot()


line.up.genes = function(locs) {

  topleft = locs$tend[1]
  row = 1
  for(i in 1:nrow(locs)) {
    if(locs$gstart[i] > topleft) {
      topleft = locs$tend[i]
      row = 1
    } # if(locs$gstart[i] > topleft)
    locs$row[i] = row
    row = row + 1
    i = i + 1
  } # for(i)
  locs$row[nrow(locs)] = row

  return(locs)
} # line.up.genes


resolve.collisions = function(locs) {

  # Look for collisions and move the gene down until it doesn't collide.
  for(i in 1:nrow(locs)) {
    # Get all of the genes in this row.
    genes = which(locs$row == i)
    # See if any overlap each other.
    overlap = genes[which(locs$gstart[genes[-1]] < locs$tend[genes[-length(genes)]])]
    if(length(overlap) > 0) {
      # Go through each overlapping gene and move it down one row until it
      # doesn't collide with any genes.
      for(j in overlap) {
        done = F
        while(!done) {
          locs$row[j] = locs$row[j] + 1
          rowgenes = which(locs$row == locs$row[j])
          curr.gene = which(rowgenes == j)
          rng = rowgenes[curr.gene + c(-1,1)]
          rng = rng[!is.na(rng)]
          if(length(rng) == 0) {
            # This occurs when we have gone past the maximum number of rows.
            done = T
          } else if(length(rng) == 1) {
            # This occurs when the gene being moved is at the beginning or 
            # end of the plot.
            if(rng <= curr.gene) {
              if(locs$gstart[j] > locs$tend[rng]) {
                 done = T
              } #
            } else {
              if(locs$tend[j] < locs$gstart[rng]) {
                 done = T
              } #            
            } # else
          } else if(length(rng) == 2) {
            if(locs$gstart[j] > locs$tend[rng[1]] & 
               locs$tend[j] < locs$gstart[rng[2]]) {
               done = T
            } #
          } # else
        } # while(!done)
      } # for(j)
    } # if(length(overlap) > 0)
  } # for(i)

  return(locs)
} # resolve.collisions()
