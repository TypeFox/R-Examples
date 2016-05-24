## TODO: figure out a better approach for alignment of dendrogram / image axis labels

# f: SPC with diagnostic boolean variables
# v: named variables
# k: number of groups to highlight
# id: id to print next to dendrogram
diagnosticPropertyPlot <- function(f, v, k, grid.label='pedon_id', dend.label='pedon_id') {
  
  # get internal, unique ID
  id <- idname(f)
  
  # extract site data
  s <- site(f)
  
  # keep only those variables that exist
  v <- names(s)[na.omit(match(v, names(s)))]
  
  ## TODO: why would there be NA in here?
  # filter NA
  no.na.idx <- which(complete.cases(s[, v]))
  s <- s[no.na.idx, ]
  
  # save diagnostic properties
  m <- s[, v]
  
  # optionally check for any vars that are all FALSE and kick them out
  vars.not.missing <- apply(m, 2, any)
  
  # if any are all FALSE, then remove from m and v
  if(any(!vars.not.missing)) {
    not.missing <- which(vars.not.missing)
    m <- m[, not.missing]
    v <- v[not.missing]
  }
  
  # convert to factors, we have to specify the levels as there are cases with all TRUE or FALSE
  m <- as.data.frame(lapply(m, factor, levels=c('FALSE', 'TRUE')))
  
  # make a copy of the matrix for plotting, as numerical data and transpose
  m.plot <- t(as.matrix(as.data.frame(lapply(m, as.numeric))))
  
  # compute dissimilarity between profiles
  d <- daisy(m, metric='gower')
  h.profiles <- as.hclust(diana(d))
  # store text labels for dendrogram
  h.profiles$labels <- as.character(s[[dend.label]]) # factors will break tiplabels()
  p <- as.phylo(h.profiles)
  
  # cut tree at user-specified number of groups
  h.cut <- cutree(h.profiles, k=k)
  
  # setup plot layout
  layout(matrix(c(1,2), nrow=1, ncol=2), widths=c(1,1))
  
  # get number of vars + number of profiles
  n.vars <- ncol(m)
  n.profiles <- nrow(m)
    
  # plot profile dendrogram
  par(mar=c(1,1,6,1))
  plot(p, cex=0.75, label.offset=0.05, y.lim=c(1.125, n.profiles))
  tiplabels(pch=15, col=h.cut, cex=1.125, adj=0.52)
  
  ## note: transpose converts logical -> character, must re-init factors
  # compute dissimilarity between variables
  d.vars <- daisy(data.frame(t(m), stringsAsFactors=TRUE), metric='gower')
  h.vars <- as.hclust(diana(d.vars))
  
  # order of profiles in dendrogram
  o.profiles <- h.profiles$order
  
  # vector of variable names as plotted in dendrogram
  o.vars <- h.vars$order
  
  # plot image matrix, with rows re-ordered according to dendrogram
  par(mar=c(1,6,6,1))
  image(x=1:n.vars, y=1:n.profiles, z=m.plot[o.vars, o.profiles], axes=FALSE, col=c(grey(0.9), 'RoyalBlue'), xlab='', ylab='', ylim=c(0.5, n.profiles+0.5))
  axis(side=2, at=1:n.profiles, labels=s[[grid.label]][o.profiles], las=1, cex.axis=0.75)
  axis(side=3, at=1:n.vars, labels=v[o.vars], las=2, cex.axis=0.75)
  abline(h=1:(n.profiles+1)-0.5)
  abline(v=1:(n.vars+1)-0.5)
  
  # return values
  rd <- cbind(s[, c(id, grid.label)], g=h.cut)
  return(invisible(list(rd=rd, profile.order=o.profiles, var.order=o.vars)))
}



diagnosticPropertyPlot2 <- function(f, v, k, grid.label='pedon_id') {
  
  # get internal, unique ID
  id <- idname(f)
  
  # extract site data
  s <- site(f)
  
  # keep only those variables that exist
  v <- names(s)[na.omit(match(v, names(s)))]
  
  ## TODO: why would there be NA in here?
  # filter NA
  no.na.idx <- which(complete.cases(s[, v]))
  s <- s[no.na.idx, ]
  
  # save grid labels
  s.gl <- as.character(s[[grid.label]])
  
  # save diagnostic properties
  m <- s[, v]
  
  # optionally check for any vars that are all FALSE and kick them out
  vars.not.missing <- apply(m, 2, any)
  
  # if any are all FALSE, then remove from m and v
  if(any(!vars.not.missing)) {
    not.missing <- which(vars.not.missing)
    m <- m[, not.missing]
    v <- v[not.missing]
  }
  
  # convert to factors, we have to specify the levels as there are cases with all TRUE or FALSE
  m <- as.data.frame(lapply(m, factor, levels=c('FALSE', 'TRUE')))
  
  # get number of vars + number of profiles
  n.vars <- ncol(m)
  n.profiles <- nrow(m)
  
  # compute dissimilarity between profiles
  d <- daisy(m, metric='gower')
  h.profiles <- as.hclust(diana(d))
  
  ## note: transpose converts logical -> character, must re-init factors
  # compute dissimilarity between variables
  d.vars <- daisy(data.frame(t(m), stringsAsFactors=TRUE), metric='gower')
  h.vars <- as.hclust(diana(d.vars))
  
  # cut tree at user-specified number of groups
  h.cut <- cutree(h.profiles, k=k)
  
  # format for plotting
  m.plot <- data.frame(id=s[[id]], m, stringsAsFactors=FALSE)
  m.plot.long <- melt(m.plot, id.vars='id')
  # convert TRUE/FALSE into factor
  m.plot.long$value <- factor(m.plot.long$value, levels=c('FALSE', 'TRUE'))
  
  # order of profiles in dendrogram
  o.profiles <- h.profiles$order
  
  # vector of variable names as plotted in dendrogram
  o.vars <- h.vars$order
  
  # set factor levels for ordering of level plot
  m.plot.long$id <- factor(m.plot.long$id, levels=m.plot$id[o.profiles])
  m.plot.long$variable <- factor(m.plot.long$variable, levels=v[o.vars])
  
  # lattice plot
  p <- levelplot(value ~ variable * id, data=m.plot.long,
  col.regions=c(grey(0.9), 'RoyalBlue'), cuts=1, xlab='', ylab='', 
  colorkey = FALSE, 
  scales=list(tck=0, x=list(rot=90), y=list(at=1:length(o.profiles), labels=s.gl[o.profiles])),
  legend=list(
      right=list(fun=dendrogramGrob, args=list(x = as.dendrogram(h.profiles), side="right", size=15, add=list(
        rect=list(fill=h.cut, cex=0.5)))),
      top=list(fun=dendrogramGrob, args=list(x=as.dendrogram(h.vars), side="top", size=4))
      ),
  panel=function(...) {
    panel.levelplot(...)
    # horizontal lines
    panel.segments(x0=0.5, y0=1:(n.profiles+1)-0.5, x1=n.vars+0.5, y1=1:(n.profiles+1)-0.5)
    # vertical lines
    panel.segments(x0=1:(n.vars+1)-0.5, y0=0.5, x1=1:(n.vars+1)-0.5, y1=n.profiles + 0.5)
  }
  )
  
  # print to graphics device
  print(p)
  
  # return values
  rd <- cbind(s[, c(id, grid.label)], g=h.cut)
  return(invisible(list(rd=rd, profile.order=o.profiles, var.order=o.vars)))
}


## failed attempt to include multi-nominal variables
## not going to work with current implementation, mostly due to how colors are mapped to values in image()

# diagnosticPropertyPlot3 <- function(f, v, k, grid.label='pedon_id', dend.label='pedon_id') {
#   
#   # get internal, unique ID
#   id <- idname(f)
#   
#   # extract site data
#   s <- site(f)
#   
#   # keep only those variables that exist
#   v <- names(s)[na.omit(match(v, names(s)))]
#   
#   ## TODO: why would there be NA in here?
#   # filter NA
#   no.na.idx <- which(complete.cases(s[, v]))
#   s <- s[no.na.idx, ]
#   
#   # save diagnostic properties
#   m <- s[, v]
#   
#   # keep track of binary / multinominal variables
#   binary.vars <- which(sapply(m, class) == 'logical')
#   multinom.vars <- which(sapply(m, class) == 'factor')
#   
#   # setup colors:
#   # binary colors, then 'NA' repeated for number of levels in any multinominal data
#   binary.cols <- c(grey(0.9), 'RoyalBlue')
#   multinom.cols <- rep(NA, times=length(levels(m[, multinom.vars])))
#   cols <- c(binary.cols, multinom.cols)
#   
#   # split: multinominal data are assumed to be a factor
#   m.binary <- m[, binary.vars, drop=FALSE]
#   if(length(multinom.vars) > 0)
#     m.multinom <- m[, multinom.vars, drop=FALSE]
#   
#   # optionally check for any vars that are all FALSE and kick them out
#   vars.not.missing <- sapply(m.binary, any)
#   
#   # if any are all FALSE, then remove from m and v
#   if(any(!vars.not.missing)) {
#     not.missing <- which(vars.not.missing)
#     m.binary <- m.binary[, not.missing]
#   }
#   
#   # convert binary data into factors, we have to specify the levels as there are cases with all TRUE or FALSE
#   m.binary <- as.data.frame(lapply(m.binary, factor, levels=c('FALSE', 'TRUE')))
#   
#   # merge binary + multinominal data
#   if(length(multinom.vars) > 0)
#     m <- cbind(m.binary, m.multinom)
#   else
#     m <- m.binary
#   
#   # update variable names, in case any were removed due to missingness
#   v <- names(m)
#   
#   # make a copy of the matrix for plotting, as numerical data and transpose
#   m.plot <- t(as.matrix(as.data.frame(lapply(m, as.numeric))))
#   
#   # compute dissimilarity between profiles
#   d <- daisy(m, metric='gower')
#   h.profiles <- as.hclust(diana(d))
#   # store text labels for dendrogram
#   h.profiles$labels <- as.character(s[[dend.label]]) # factors will break tiplabels()
#   p <- as.phylo(h.profiles)
#   
#   # cut tree at user-specified number of groups
#   h.cut <- cutree(h.profiles, k=k)
#   
#   # setup plot layout
#   layout(matrix(c(1,2), nrow=1, ncol=2), widths=c(1,1))
#   
#   # get number of vars + number of profiles
#   n.vars <- ncol(m)
#   n.profiles <- nrow(m)
#   
#   # plot profile dendrogram
#   par(mar=c(1,1,6,1))
#   plot(p, cex=0.75, label.offset=0.05, y.lim=c(1.125, n.profiles))
#   tiplabels(pch=15, col=h.cut, cex=1.125, adj=0.52)
#   
#   ## note: transpose converts logical -> character, must re-init factors
#   # compute dissimilarity between variables
#   d.vars <- daisy(data.frame(t(m), stringsAsFactors=TRUE), metric='gower')
#   h.vars <- as.hclust(diana(d.vars))
#   
#   # order of profiles in dendrogram
#   o.profiles <- h.profiles$order
#   
#   # vector of variable names as plotted in dendrogram
#   o.vars <- h.vars$order
#   
#   # plot image matrix, with rows re-ordered according to dendrogram
#   par(mar=c(1,6,6,1))
#   image(x=1:n.vars, y=1:n.profiles, z=m.plot[o.vars, o.profiles], axes=FALSE, col=cols, xlab='', ylab='', ylim=c(0.5, n.profiles+0.5))
#   axis(side=2, at=1:n.profiles, labels=s[[grid.label]][o.profiles], las=1, cex.axis=0.75)
#   axis(side=3, at=1:n.vars, labels=v[o.vars], las=2, cex.axis=0.75)
#   abline(h=1:(n.profiles+1)-0.5)
#   abline(v=1:(n.vars+1)-0.5)
#   
#   # return values
#   rd <- cbind(s[, c(id, grid.label)], g=h.cut)
#   return(invisible(list(rd=rd, profile.order=o.profiles, var.order=o.vars)))
# }

