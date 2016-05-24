"plot.strat.column" <-
function(x = NULL,
         counts = NULL,
         depths = NULL,
         at.depths = NULL,
         sample.labels = NULL,
         taxa = NULL,
         short.names = NULL,
         higher.grp = NULL,
         tax.cat = NULL,
         metadata = NULL,
         prop.cutoff = NULL,
         cols = c('depths', 'tick', 'samples', 'tick', 'blank',
                  'data',
                  'tick', 'depths', 'totals'),
         colwidths = NULL,
         output.size = NULL,
         output.file = NULL,
         style = 'mountains',
         reorder = NULL,
         outer.margins = list(bottom = unit(1, "lines"),
                              left = unit(1, "lines"),
                              top = unit(1, "lines"),
                              right = unit(1, "lines")),
         fontsize = 9,
         width.multiplier = 2,
         height.increment = 0.25, # inches per stratigraphic level
         fill.cols = NULL,
         line.cols = NULL,
         box.pars = NULL,
         proportional = FALSE,
         use.color = TRUE,
         bottom.axis = TRUE,
         count.label = TRUE,
         sample.totals = TRUE,
         sample.category.totals = FALSE,
         sample.nos = TRUE,
         show.taxon.count.totals = TRUE,
         plot.depths.increasing.down = TRUE,
         debug = FALSE,
         ...){

## TO DO
# plot.log = (plotting an arbitrary data vector in blank column)
# absolute.dates = NULL
# clust = 'samples', 'taxa', '2way'

if(is.strat.column(x)){
  counts <- x$counts
  if(is.null(depths)) depths <- x$depths
  if(is.null(taxa)) taxa <- x$taxa
  if(is.null(short.names)) short.names <- x$short.names
  if(is.null(higher.grp)) higher.grp <- x$higher.grp
  if(is.null(tax.cat)) tax.cat <- x$tax.cat
}

# Read in the count data:
if(is.data.frame(counts) || is.matrix(counts)){
  counts <- apply(counts, 2, 'as.numeric')
}else if(is.character(counts) && length(counts) == 1){
  counts <- read.csv(counts, header = TRUE, skip = 0, colClasses = '')
  counts <- apply(counts, 2, 'as.numeric')
}else{
  stop('argument to counts not understood')
}

# Check that depths is the right length
if(is.null(depths)){
  warning('no depths provided; plotting samples at regular intervals')
  depths <- 1:nrow(counts)
}
depths <- as.numeric(depths)
if(length(depths) != nrow(counts)){
  stop(paste(length(depths), ' depths, and ', nrow(counts),
             ' rows in the count matrix.', sep = ''))
}
if(is.null(at.depths)) at.depths <- depths

# Check that sample.labels is the right length
if(is.null(sample.labels)){
  warning('no sample.labels provided; assigning row numbers as labels')
  sample.labels <- 1:nrow(counts)
}
if(length(sample.labels) != nrow(counts)){
  stop(paste(length(sample.labels), ' sample labels, and ', nrow(counts),
             ' rows in the count matrix.', sep = ''))
}

# Replace NA in counts with zeros
if(sum(is.na(counts)) > 0){
  warning(paste(sum(is.na(counts)), 'missing values in count matrix replaced with zeros'))
  counts[is.na(counts)] <- 0
}

# Remove empty columns
emptycols <- !(colSums(counts, na.rm = TRUE) > 0)
if(any(emptycols)){
  counts <- counts[,!emptycols]
  tax.cat <- tax.cat[!emptycols]
  taxa <- taxa[!emptycols]
  warning(paste(sum(emptycols),
                'columns with zero counts at all levels removed'))
}

# Deal with prop.cutoff
cutoffcols <- colSums(counts)/sum(counts) < prop.cutoff
if(any(cutoffcols)){
  counts <- counts[,!cutoffcols]
  tax.cat <- tax.cat[!cutoffcols]
  taxa <- taxa[!cutoffcols]
  warning(paste(sum(cutoffcols),
     'columns with total counts below prop.cutoff removed'))
}

# Convert to proportions/relative abundances (if flag set)
if(proportional){
  prop.funct <- function(x){
    x / sum(x, na.rm = TRUE)
  }
  props <- apply(counts, 1, prop.funct)
  props <- t(props)
  show.taxon.count.totals <- FALSE
  count.label <- FALSE
  bottom.axis <- FALSE
}

# Deal with tax.cat
if(!is.null(tax.cat)){
  if(length(tax.cat) != ncol(counts)){
    warning('taxon category labels seem to be the wrong length')
    tax.cat <- NULL
  }
}
if(is.null(tax.cat)){
  tax.cat <- rep('', ncol(counts))
}

# Deal with reorder
if(!is.null(reorder)){
  if(length(reorder) == 1 && is.character(reorder)){
    funny <- function(x) return((1:length(x))[x > 0])
    if(pmatch(reorder, 'fad.by.category', nomatch = FALSE)){
      fads <- depths[as.numeric(lapply(apply(counts, 2, funny), max))]
      reorder.vect <- sort(fads, decreasing = TRUE,
                             index.return = TRUE)$ix
      counts <- counts[,reorder.vect]
      tax.cat <- tax.cat[reorder.vect]
      reorder.vect <- sort(as.character(tax.cat),
                             index.return = TRUE)$ix
      counts <- counts[,reorder.vect]
      tax.cat <- tax.cat[reorder.vect]
    }else if(pmatch(reorder, 'lad.by.category', nomatch = FALSE)){
      lads <- depths[as.numeric(lapply(apply(counts, 2, funny), min))]
      reorder.vect <- sort(lads, decreasing = TRUE,
                      index.return = TRUE)$ix
      counts <- counts[,reorder.vect]
      tax.cat <- tax.cat[reorder.vect]
      reorder.vect <- sort(as.character(tax.cat),
                      index.return = TRUE)$ix
      counts <- counts[,reorder.vect]
      tax.cat <- tax.cat[reorder.vect]
    }else if(pmatch(reorder, 'by.count', nomatch = FALSE)){
      reorder.vect <- sort(colSums(counts), decreasing = TRUE,
                           index.return = TRUE)$ix
      counts <- counts[,reorder.vect]
      tax.cat <- tax.cat[reorder.vect]
    }
  }else if(length(reorder) == ncol(counts)){
    reorder.vect <- as.numeric(reorder)
    counts <- counts[,reorder.vect]
    tax.cat <- tax.cat[reorder.vect]
  }else{
    stop('argument to reorder not understood')
  }
}

# Make up a vector of colours (fill.cols and line.cols) for the actual spindles

if(length(fill.cols) != ncol(counts)){
  if(is.null(fill.cols)){
    fill.cols <- rep('green', ncol(counts))
  }else if(length(fill.cols) == 1){
    fill.cols <- rep(fill.cols, ncol(counts))
  }else{
    fill.cols <- rep('black', ncol(counts))  
    warning('argument to fill.cols not understood; using black')
  }
}

if(length(line.cols) != ncol(counts)){
  if(is.null(line.cols)){
    line.cols <- rep('black', ncol(counts))
  }else if(length(line.cols) == 1){
    line.cols <- rep(line.cols, ncol(counts))
  }else{
    line.cols <- rep('black', ncol(counts))  
    warning('argument to line.cols not understood; using black')
  }
}

if(!use.color){
  fill.cols <- rep('black', ncol(counts))
  line.cols <- rep('black', ncol(counts))
  box.pars <- 'bw'
}

if(debug) box.pars <- 'debug.colors'

# Box styles
if(is.null(box.pars)){
  b1 <- gpar(lwd = 0.5, fill = NULL, col = 'white') #around data
  b2 <- gpar(lwd = 0.5, fill = NULL, col = 'black')
  #b3
  #b3A
  b4 <- gpar(lwd = 1, fill = NULL, col = NULL) #left buffer
  b5 <- gpar(lwd = 1, fill = NULL, col = NULL) #taxon counts
  b6 <- gpar(lwd = 1, fill = NULL, col = 'black') #taxon category
  b7 <- gpar(lwd = 1, fill = NULL, col = 'black') #taxon names
  b8 <- gpar(lwd = 1, fill = NULL, col = 'black') #headings
  b9 <- gpar(lwd = 1, fill = NULL, col = NULL) #left margin
  b10 <- gpar(lwd = 1, fill = NULL, col = NULL)	#both margins
  b11 <- gpar(lwd = 1.5, fill = NULL, col = NULL) #NA
  b12 <- gpar(lwd = 1.5, fill = NULL, col = NULL) #right margin
  b13 <- gpar(lwd = 1.5, fill = NULL, col = 'black') #metadata
  b14 <- gpar(lwd = 2, fill = NULL, col = 'black') #border
}

if(length(box.pars) == 1 && is.character(box.pars)){
  if(box.pars == 'bw'){
    b1 <- gpar(lwd = 0.5, fill = NULL, col = 'white')
    b2 <- gpar(lwd = 0.5, fill = NULL, col = 'black')
    #b3
    #b3A
    b4 <- gpar(lwd = 1, fill = NULL, col = 'black')
    b5 <- gpar(lwd = 1, fill = NULL, col = NULL)
    b6 <- gpar(lwd = 1, fill = NULL, col = NULL)
    b7 <- gpar(lwd = 1, fill = NULL, col = NULL)
    b8 <- gpar(lwd = 1, fill = NULL, col = NULL)
    b9 <- gpar(lwd = 1, fill = NULL, col = NULL)
    b10 <- gpar(lwd = 1, fill = NULL, col = NULL)
    b11 <- gpar(lwd = 1.5, fill = NULL, col = 'black')
    b12 <- gpar(lwd = 1.5, fill = NULL, col = 'black')
    b13 <- gpar(lwd = 1.5, fill = NULL, col = 'black')
    b14 <- gpar(lwd = 2, fill = NULL, col = 'black')
  }else if(box.pars == 'debug.colors'){
    b1 <- gpar(lwd = 0.5, fill = 'yellow', col = NULL)
    b2 <- gpar(lwd = 0.5, fill = NULL, col = 'black')
    #b3
    #b3A
    b4 <- gpar(lwd = 1, fill = grey(0.6), col = 'black')
    b5 <- gpar(lwd = 1, fill = 'orange', col = 'black')
    b6 <- gpar(lwd = 1, fill = 'lightblue', col = 'black')
    b7 <- gpar(lwd = 1, fill = 'wheat', col = 'black')
    b8 <- gpar(lwd = 1, fill = 'lightpink', col = 'black')
    b9 <- gpar(lwd = 1, fill = 'lightgreen', col = 'black')
    b10 <- gpar(lwd = 1, fill = 'tomato1', col = 'black')
    b11 <- gpar(lwd = 1.5, fill = 'plum3', col = 'black')
    b12 <- gpar(lwd = 1.5, fill = 'steelblue1', col = 'black')
    b13 <- gpar(lwd = 1.5, fill = grey(0.9), col = 'black')
    b14 <- gpar(lwd = 2, fill = NULL, col = 'black')   
  }
}

# Fill cols argument
cols <- rep(cols, times =  1 + ((cols == 'data') * (ncol(counts) - 1)))
offset <- min(grep('data', cols)) - 1
max.counts <- as.numeric(apply(counts, 2, max))
width.mult.vector <- as.numeric(cols == 'data') + 1
width.mult.vector[width.mult.vector == 2] <- 
   (max.counts/max(max.counts) * (width.multiplier - 1)) + 1

lineh <- fontsize * get.gpar('cex')$cex *
            get.gpar('lineheight')$lineheight
dev.off() # gets rid of the graphics device started by get.gpar()
# lineh is minimum height (width) of a data column in points

offset.from.max <- 5 # in points
#if(offset.from.max > (0.5 * lineh)){
#  stop('is this a problem? see variable offset.from.max')
#}
# offset.from.max is min. necessary extra space for
#  each data column in points

# Calculate the size of the plot if not specified
if(is.null(output.size)){ #NO OVERALL PLOT SIZE SPECIFIED

  dat.col.width <- unit(lineh, 'points') + unit(5, 'points')

  # Get widths for columns...
  longest.depth <- unit(1, 'strwidth', as.character(max(depths[nchar(depths)
                           == max(nchar(depths))])))
  longest.samp <- unit(1, 'strwidth',
                       as.character(max(depths[nchar(sample.labels)
                          == max(nchar(sample.labels))])))
  if(is.null(colwidths)){
    for(i in seq_along(cols)){
      if(i == 1){
        if(cols[i] == 'data') colwidths <- dat.col.width
        if(cols[i] == 'depths') colwidths <- longest.depth
        if(cols[i] == 'tick') colwidths <- unit(1/16, 'inches')
        if(cols[i] == 'samples') colwidths <- longest.samp
        if(cols[i] == 'totals') colwidths <- unit(1, 'inches')
        if(cols[i] == 'blank') colwidths <- unit(0.5, 'inches')
      }else{
  	    if(cols[i] == 'data') colwidths <- unit.c(colwidths, dat.col.width)
  	    if(cols[i] == 'depths') colwidths <- unit.c(colwidths,
  	       longest.depth)
        if(cols[i] == 'tick') colwidths <- unit.c(colwidths,
           unit(1/16, 'inches'))
        if(cols[i] == 'samples') colwidths <- unit.c(colwidths,
           longest.samp)
        if(cols[i] == 'totals') colwidths <- unit.c(colwidths,
           unit(1, 'inches'))
        if(cols[i] == 'blank') colwidths <- unit.c(colwidths,
           unit(0.5, 'inches'))
      }
    }
  }
  colwidths <- colwidths * width.mult.vector
  inner.width <- convertX(sum(colwidths), 'inches')
  dev.off() # gets rid of the graphics device started by convertX()

if(debug) print(inner.width)

  width <- as.numeric(inner.width) + 3 #in inches
  height <- (height.increment * nrow(counts)) + 5 # in inches
}else if(length(output.size == 2)){ #OVERALL PLOT SIZE SPECIFIED
  width <- output.size[1]
  height <- output.size[2]
  inner.width <- unit(width - 3, 'inches')
  if(is.null(colwidths)){
    colwidths <- rep(inner.width * (1/length(cols)), length(cols))
  }
}else{
  stop('arguments to output.size not understood')
}

if(debug){
  print(cols)
  print(colwidths)
}

if(debug){
  cat('Specified output size is:', output.size, 'in inches\n')
  cat('  Height:', height, 'inches\n')
  cat('  Width:', width, 'inches\n')
}

metadata <- x$metadata
if(is.null(metadata)){
  warning('typically required metadata includes: Site Name, Lat./Long., Country, Principle investigator with contact information, and information on absolute dates of section')
  metadata <- list(SiteName = 'Site Name Unspecified',
                   Latitude = '0 N',
                   Longitude = '0 E',
                   Country = 'Country Unspecified',
                   Location = 'No additional information',
                   SiteType = 'Unknown',
                   Interval = '',
                   Scale = '',
                   StartDate = '',
                   EndDate = '',
                   PlotDate = date(),
                   Contact = 'Unknown',
                   Publications = '',
                   Notes = '')
}

################################################################################
# Subroutines

## NOTE: Viewports are numbered from inside (starting with 1) out

spplot <- function(count.vector, depths, cat.label,
                   fill.col, line.col,
                   taxon, proportional, count.label = TRUE,
                   left.axis = FALSE,
                   bottom.axis = TRUE, top.axis = FALSE){

  # Rectangle/background for whole taxon share of plot width (vp2)
  grid.rect(gp = b2)

  # Print taxon name aligned with center of taxon share
  #  of plot width
  #grid.text(label = taxon, x = 0.5,
  #          y = unit(1, 'npc') + unit(1, 'char'),
  #          rot = 90, just = 'left')

  vp2 <- viewport(x = 0, y = 0,
                  width = unit(1, 'npc') -
                      unit(offset.from.max, 'points'),
                  height = 1,
                  just = c('left', 'bottom'), clip = 'off')
  pushViewport(vp2)

  # Rectangle/background for just around spindle (vp1)
  #  probably only for debugging use
  grid.rect(gp = b1)

  # Print taxon name aligned with left side of taxon share
  #  of plot width
  grid.text(label = taxon, x = 0,
            y = unit(1, 'npc')
                + unit(2, 'points')
                + unit(2, 'lines'),
            rot = 90, just = c('left', 'top'))

  # Print tax.cat (or whatever label is passed to 'tax.cat')
  if(!is.null(cat.label)){
    if(show.taxon.count.totals){
      grid.text(label = as.character(cat.label), x = 0,
                y = unit(1, 'npc') +
                    unit((lineh / 2) - 2, 'points') +
                    unit(lineh, 'points'),
                rot = 0, just = c('left', 'bottom'),
                gp = gpar(cex = 0.8))
    }else{
      grid.text(label = as.character(cat.label), x = 0,
                y = unit(1, 'npc') +
                    unit((lineh / 2) - 2, 'points'),
                rot = 0, just = c('left', 'bottom'),
                gp = gpar(cex = 0.8))
    }
  }

  # Print taxon sums
  if(show.taxon.count.totals){
    grid.text(label = sum(count.vector), x = 0,
              y = unit(1, 'npc') +
              unit((lineh / 2) - 2, 'points'),
              rot = 0, just = c('left', 'bottom'),
              gp = gpar(cex = 0.8))
  }
  
  vp1 <- viewport(xscale = c(0, max(count.vector)),
                 yscale = c(max(depths), min(depths)),
                 default.units = 'native',
                 clip = 'off')
  pushViewport(vp1)

  if(pmatch(style,'mountains', nomatch = FALSE)){
    grid.polygon(x = c(0, count.vector, 0),
                 y = c(depths[1], depths, depths[length(depths)]),
                 default.units = 'native',
                 gp = gpar(fill = fill.col,
                           col = line.col, lwd = 0.75))
    grid.segments(x0 = rep(0, length(depths)),
                  y0 = depths,
                  x1 = count.vector,
                  y1 = depths,
                  default.units = 'native',
                  gp = gpar(col = 'white', lwd = 0.5))
  }else if(pmatch(style,'bars', nomatch = FALSE)){
    grid.rect(x = unit(0, 'npc'), y = depths,
              width = count.vector,
              height = unit(0.15 * height.increment, 'inches'),
              default.units = 'native',
              gp = gpar(fill = fill.col, col = line.col, lwd = 0.75),
              just = c('left','centre'))
  }else if(pmatch(style,'ranges', nomatch = FALSE)){

    occur <- as.logical(count.vector)
    s.range <- range(depths[occur])
    grid.rect(x = unit(0.5, 'npc'), y = s.range[1], width = unit(0.7, 'char'),
              height = s.range[2] - s.range[1],
              default.units = 'native',
              gp = gpar(fill = 'black', col = 'black', lwd = 0.25),
              just = c('center', 'bottom'))
    if(FALSE){ # to be replaced by a condition later....
      grid.rect(x = rep(unit(0.5, 'npc'), length(depths))[occur],
                y = depths[occur], width = unit(0.7, 'char'),
                height = unit(1, 'char'),
                default.units = 'native',
                gp = gpar(fill = 'white', col = 'black', lwd = 0.25))
    }else if(TRUE){ # to be replaced by a condition later....
      grid.points(x = rep(unit(0.5, 'npc'), length(depths))[occur],
                  y = depths[occur],
                  default.units = 'native',
                  pch = 21, gp = gpar(fill = 'white', col = NULL, cex = 0.7))
    }else{ # currently not operational....
      na.counts <- as.character(count.vector)
      na.counts[na.counts == 0] <- ''
      grid.text(label = na.counts[occur],
                  x = rep(unit(0.5, 'npc'), length(depths))[occur],
                  y = depths[occur],
                  default.units = 'native',
                  gp = gpar(cex = 0.6))
    }
    count.label <- FALSE
    bottom.axis <- FALSE
  }else{
    stop('argument to style not understood')
  }
  
  if(count.label){
    na.counts <- as.character(count.vector)
    na.counts[na.counts == 0] <- ''
    grid.text(label = na.counts, x = count.vector, y = depths,
              default.units = 'native',
              hjust = -0.25, vjust = -0.5,
              gp = gpar(cex = 0.6))
  }
  if(bottom.axis) grid.xaxis()
  if(is.null(sample.labels)) sample.labels <- rep('', length(depths))
  if(left.axis == 'depths') grid.yaxis()
  else if(left.axis == 'labels') grid.yaxis(at = depths, label = sample.labels)
  popViewport() #pop vp1
  popViewport() #pop vp2
}

################################################################################
# Main

if(!is.null(output.file)){
  pdf(file = output.file, height = height,
      width = width, family = 'Times', pointsize = fontsize)
}else{
  #warning('displaying a quartz window only works on a macintosh; specify an output.file name on all other platforms')
  dev.new(height = height, width = width, family = 'Times', pointsize = fontsize)
}

# The largest thing drawn
vp14 <- viewport(x = outer.margins$left,
                 y = outer.margins$bottom, 
                 width = unit(1, "npc") -
                     outer.margins$right -
                     outer.margins$left, 
                 height = unit(1, "npc") -
                     outer.margins$top -
                     outer.margins$bottom, 
                 just = c("left", "bottom"), clip = "off")
pushViewport(vp14)

grid.rect(gp = b14) #largest thing drawn

# Metadata viewport
vp13 <- viewport(x = 0.5, y = 0.5, 
        width = unit(1, "npc") - unit(2, 'lines'),
        height = unit(1, "npc") - unit(2, 'lines'), 
        just = c("centre", "centre"), clip = "off")
pushViewport(vp13)

# Probably only for debugging
grid.rect(gp = b13)

## PRINT METADATA HERE

# Calculate some statistics
total.counts <- sum(counts, na.rm = TRUE)
counts.per.sample <- round(total.counts/length(depths))

grid.text(label = metadata$SiteName,
          x = unit(0, 'npc') + unit(0.2, 'inches'),
          y = unit(1, 'npc') - unit(0.2, 'inches'),
          gp = gpar(fontsize = 18), just = c('left', 'center'))

grid.text(label = paste('Lat./Long.: ', metadata$Latitude,'/',
                        metadata$Longitude, sep = ''),
          x = unit(0, 'npc') + unit(0.2, 'inches'),
          y = unit(1, 'npc') - unit(0.2, 'inches') - unit(20, 'points'),
          gp = gpar(fontsize = 10), just = c('left', 'center'))

grid.text(label = paste('Location: ', metadata$Country,', ',
                        metadata$Location, sep = ''),
          x = unit(0, 'npc') + unit(0.2, 'inches'),
          y = unit(1, 'npc') - unit(0.2, 'inches') - unit(32, 'points'),
          gp = gpar(fontsize = 10), just = c('left', 'center'))

grid.text(label = paste('Site Type: ', metadata$SiteType,' Interval: ',
                        metadata$Interval, ' Scale: ', metadata$Scale,
                        ' StartDate: ', metadata$StartDate,
                        ' EndDate: ', metadata$EndDate, sep = ''),
          x = unit(0, 'npc') + unit(0.2, 'inches'),
          y = unit(1, 'npc') - unit(0.2, 'inches') - unit(44, 'points'),
          gp = gpar(fontsize = 10), just = c('left', 'center'))

grid.text(label = paste('Contact: ', metadata$Contact,', or see: ',
                        metadata$Publications, sep = ''),
          x = unit(0, 'npc') + unit(0.2, 'inches'),
          y = unit(1, 'npc') - unit(0.2, 'inches') - unit(56, 'points'),
          gp = gpar(fontsize = 10), just = c('left', 'center'))

grid.text(label = paste('Notes: ', metadata$Notes, sep = ''),
          x = unit(0, 'npc') + unit(0.2, 'inches'),
          y = unit(1, 'npc') - unit(0.2, 'inches') - unit(68, 'points'),
          gp = gpar(fontsize = 10), just = c('left', 'center'))

grid.text(label = paste('Plotted: ', metadata$PlotDate, sep = ''),
          x = unit(0, 'npc') + unit(0.2, 'inches'),
          y = unit(1, 'npc') - unit(0.2, 'inches') - unit(80, 'points'),
          gp = gpar(fontsize = 10), just = c('left', 'center'))

grid.text(label = paste('Total of ', total.counts,
                        ' items counted; approximately ',
                        counts.per.sample,
                        ' per sample', sep = ''),
          x = unit(0, 'npc') + unit(0.2, 'inches'),
          y = unit(1, 'npc') - unit(0.2, 'inches') - unit(92, 'points'),
          gp = gpar(fontsize = 10), just = c('left', 'center'))

metadata.space <- unit(2, 'inches')

# Taxon category heading viewport
vp12 <- viewport(x = 0, y = 0, 
        width = unit(1, "npc"), height = unit(1, "npc") - metadata.space, 
        just = c("left", "bottom"), clip = "off")
pushViewport(vp12)
grid.rect(gp = b12)

sample.totals.space <- unit(1.5, 'inches')

# Space to plot total abundances for each sample (depth)
if(sample.totals){
  vp11 <- viewport(x = 0, y = 0, 
                   width = unit(1, "npc") -
                       sample.totals.space,
                   height = unit(1, "npc"), 
                   just = c("left", "bottom"), clip = "off")
  pushViewport(vp11)
}
grid.rect(gp = b11)

sample.category.totals.space <- unit(0.5, 'inches')
#Space to plot taxon category total abundances for each sample
if(sample.category.totals){
  vp10 <- viewport(x = 0, y = 0, 
          width = unit(1, "npc") - sample.category.totals.space,
          height = unit(1, "npc"), 
          just = c("left", "bottom"), clip = "off")
  pushViewport(vp10)
}
grid.rect(gp = b10)

sample.labels.space <- unit(2, 'lines')
#Space to plot sample labels
if(sample.nos){
  vp9 <- viewport(x = unit(0, 'npc') + unit(1, 'lines'),
                  y = 0, 
                  width = unit(1, "npc") - 2 * sample.labels.space,
                  height = unit(1, "npc"), 
          just = c("left", "bottom"), clip = "off")
  pushViewport(vp9)
}
grid.rect(gp = b9)

strat.space = unit(6, 'lines')
grid.text(label = 'Depth',
          x = unit(0, 'npc') +
              (strat.space * (1 / 1.5))
              - unit(0.5, 'lines'),
          y = unit(1, 'npc') - unit(1, 'inches'),
          just = c('right','center'),
          gp = gpar(fontface = 'bold', cex = 1.5))
if(debug){
  grid.text(label = '<-strat.space->',
            x = unit(0, 'npc'),
            y = unit(1, 'npc') - unit(2, 'lines'),
            just = c('left','center'))
}

# Leaving space for the depth labels and strat column (if any)
vp8 <- viewport(x = 1, y = 0,
        width = unit(1, "npc") - strat.space, height = unit(1, "npc"), 
        just = c("right", "bottom"), clip = "off")
pushViewport(vp8)
grid.rect(gp = b8)

# Leaving one line for the taxon category heading labels
vp7 <- viewport(x = 0, y = 0,
                width = unit(1, "npc"),
                height = unit(1, "npc") - unit(1, 'lines'), 
                just = c("left", "bottom"), clip = "off")
pushViewport(vp7)
grid.rect(gp = b7)

namewidths <- strwidth(colnames(counts), units = 'inches')
longestnw <- max(strwidth(colnames(counts), units = 'inches'))
longest.taxon.name <- (1:ncol(counts))[namewidths == longestnw]

# Leaving room for species names
vp6 <- viewport(x = 0, y = unit(3, 'lines'),
        width = unit(1, "npc"),
        height = unit(1, "npc") - unit(4, 'lines') -
            unit(1, 'strwidth',
                 data = colnames(counts)[longest.taxon.name]), 
        just = c("left", "bottom"), clip = "off")
pushViewport(vp6)
grid.rect(gp = b6)

# Leaving one line for taxon count totals
if(show.taxon.count.totals){
  vp5 <- viewport(x = 0, y = 0,
          width = unit(1, "npc"),
          height = unit(1, "npc") - unit(1, 'lines'), 
          just = c("left", "bottom"), clip = "off")
  pushViewport(vp5)
  grid.rect(gp = b5)
}

# Finally the plotting area, leaving one line for taxon category labels
vp4 <- viewport(x = 0, y = 0,
        width = unit(1, "npc"),
        height = unit(1, "npc") - unit(1, 'lines'), 
        just = c("left", "bottom"), clip = "off")
pushViewport(vp4)
grid.rect(gp = b4)

# Preparatory calculations for plotting actual data columns

xloc <- unit(0, 'lines')
for(i in seq_along(cols)){ #loop through columns
  vpCol <- viewport(x = xloc,
                  y = 0, height = 1,
                  width = colwidths[i],
                  just = c('left', 'bottom'))
  xloc <- xloc + colwidths[i]
  pushViewport(vpCol)
  left.axis <- FALSE
  if(i == 1){
    vpL <- viewport(xscale = c(0,1),
                   yscale = c(max(depths), min(depths)),
                   default.units = 'native',
                   clip = 'off')
    pushViewport(vpL)
    grid.yaxis()
    popViewport() #pop vpL
  }  
  if(cols[i] == 'data'){
  	i <- i - offset
  	if(proportional) cv <- props[,i]
  	else cv <- counts[,i]
    spplot(count.vector = cv, proportional = proportional,
           depths = depths, cat.label = tax.cat[i],
           fill.col = fill.cols[i], line.col = line.cols[i],
           taxon = colnames(counts)[i],
           left.axis = left.axis,
           count.label = count.label,
           bottom.axis = bottom.axis)
  }else if(cols[i] == 'depths'){
    vpDep <- viewport(xscale = c(0,1),
                   yscale = c(max(depths), min(depths)),
                   default.units = 'native',
                   clip = 'off')
    pushViewport(vpDep)
    grid.text(label = at.depths, x = unit(0.5, 'npc'),
              y = unit(at.depths, 'native'))
    grid.rect(x = 0, y = 0, gp = gpar(col = 'black', lwd = 0.5),
              just = c('left', 'bottom'))
    popViewport() #pop vpDep
  }else if(cols[i] == 'tick'){
    vpDT <- viewport(xscale = c(0,1),
                   yscale = c(max(depths), min(depths)),
                   default.units = 'native',
                   clip = 'off')
    pushViewport(vpDT)
    grid.segments(unit(rep(0, length(at.depths)), 'npc'),
                  unit(at.depths, 'native'),
                  unit(rep(1, length(at.depths)), 'npc'),
                  unit(at.depths, 'native'),
                  gp = gpar(col = 'black', lwd = 0.5))
    popViewport() #pop vpDT
  }else if(cols[i] == 'samples'){
    if(is.null(at.depths)) at.depths <- depths
    vpSamp <- viewport(xscale = c(0,1),
                   yscale = c(max(depths), min(depths)),
                   default.units = 'native',
                   clip = 'off')
    pushViewport(vpSamp)
    grid.text(label = sample.labels, x = unit(0.5, 'npc'),
              y = unit(at.depths, 'native'))
    grid.rect(x = 0, y = 0, gp = gpar(col = 'black', lwd = 0.5),
              just = c('left', 'bottom'))
    popViewport() #pop vpSamp

  }else if(cols[i] == 'blank'){
    grid.rect()
  }else if(cols[i] == 'totals'){
  
  }else{
  	warning('argument to cols not understood')
    grid.rect(x = 0, y = 0, gp = gpar(col = 'red', lwd = 2),
              just = c('left', 'bottom'))
  }
  popViewport() #pop vpCol
}






if(sample.totals && !sample.category.totals){

  if(use.color){
    totals.col <- 'brown'
  }else{
    totals.col <- 'black'
  }
  width.of.total.plot <- unit(2, 'inches') - 
                         sample.category.totals.space

  vp3A <- viewport(x = unit(1, 'npc') - unit(1, 'lines') + 
                       2 * sample.labels.space, 
                   y = 0, width = width.of.total.plot,
                   just = c('left','bottom'))
  pushViewport(vp3A)

  spplot(count.vector = rowSums(counts), depths = depths,
         proportional = FALSE,
         fill.col = totals.col, line.col = 'black',
         taxon = 'TOTAL', cat.label = NULL,
         left.axis = 'labels')
  popViewport() #pop vp3A
}

popViewport() #pop vp4
if(show.taxon.count.totals) popViewport() #pop vp5
popViewport() #pop vp6
popViewport() #pop vp7
popViewport() #pop vp8
if(sample.nos) popViewport() #pop vp9
if(sample.category.totals) popViewport() #pop vp10
if(sample.totals) popViewport() #pop vp11
popViewport() #pop vp12
popViewport() #pop vp13
popViewport() #pop vp14

if(!is.null(output.file)) dev.off()

} # End of function

#roving debugger
#if(debug){
#  grid.rect(gp = gpar(col = 'black', lwd = 2))
#  grid.segments(gp = gpar(col = 'black', lwd = 2))
#  grid.text(label = '10',
#            x = unit(0, 'npc'),
#            y = unit(1, 'npc'),
#            gp = gpar(col = 'black', cex = 0.8),
#            just = c('left', 'top'))
#}