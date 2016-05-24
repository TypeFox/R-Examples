#' 
#' @name GOCluster
#' @title Circular dendrogram.
#' @description GOCluster generates a circular dendrogram of the \code{data} 
#'   clustering using by default euclidean distance and average linkage.The 
#'   inner ring displays the color coded logFC while the outside one encodes the
#'   assigned terms to each gene.
#' @param data A data frame which should be the result of 
#'   \code{\link{circle_dat}} in case the data contains only one logFC column. 
#'   Otherwise \code{data} is a data frame whereas the first column contains the
#'   genes, the second the term and the following columns the logFCs of the 
#'   different contrasts.
#' @param process A character vector of selected processes (ID or term
#'   description)
#' @param metric A character vector specifying the distance measure to be used 
#'   (default='euclidean'), see \code{dist}
#' @param clust A character vector specifying the agglomeration method to be 
#'   used (default='average'), see \code{hclust}
#' @param clust.by A character vector specifying if the clustering should be 
#'   done for gene expression pattern or functional categories. By default the 
#'   clustering is done based on the functional categories.
#' @param nlfc If TRUE \code{data} contains multiple logFC columns (default= 
#'   FALSE)
#' @param lfc.col Character vector to define the color scale for the logFC of 
#'   the form c(high, midpoint,low)
#' @param lfc.min Specifies the minimium value of the logFC scale (default = -3)
#' @param lfc.max Specifies the maximum value of the logFC scale (default = 3)
#' @param lfc.space The space between the leafs of the dendrogram and the ring 
#'   for the logFC
#' @param lfc.width The width of the logFC ring
#' @param term.col A character vector specifying the colors of the term bands
#' @param term.space The space between the logFC ring and the term ring
#' @param term.width The width of the term ring
#' @details The inner ring can be split into smaller rings to display multiply
#'   logFC values resulting from various comparisons.
#' @import ggplot2
#' @import ggdendro
#' @import RColorBrewer
#' @import stats
#' @examples
#' \dontrun{
#' #Load the included dataset
#' data(EC)
#' 
#' #Generating the circ object
#' circ<-circular_dat(EC$david, EC$genelist)
#' 
#' #Creating the cluster plot
#' GOCluster(circ, EC$process)
#' 
#' #Cluster the data according to gene expression and assigning a different color scale for the logFC
#' GOCluster(circ,EC$process,clust.by='logFC',lfc.col=c('darkgoldenrod1','black','cyan1'))
#' }
#' @export
#' 

GOCluster<-function(data, process, metric, clust, clust.by, nlfc, lfc.col, lfc.min, lfc.max, lfc.space, lfc.width, term.col, term.space, term.width){
  x <- y <- xend <- yend <- width <- space <- logFC <- NULL
  if (missing(metric)) metric<-'euclidean'
  if (missing(clust)) clust<-'average'
  if (missing(clust.by)) clust.by<-'term'
  if (missing(nlfc)) nlfc <- 0
  if (missing(lfc.col)) lfc.col<-c('firebrick1','white','dodgerblue')
  if (missing(lfc.min)) lfc.min <- -3
  if (missing(lfc.max)) lfc.max <- 3
  if (missing(lfc.space)) lfc.space<- (-0.5) else lfc.space<-lfc.space*(-1)
  if (missing(lfc.width)) lfc.width<- (-1.6) else lfc.width<-lfc.space-lfc.width-0.1
  if (missing(term.col)) term.col<-brewer.pal(length(process), 'Set3')
  if (missing(term.space)) term.space<- lfc.space+lfc.width else term.space<-term.space*(-1)+lfc.width
  if (missing(term.width)) term.width<- 2*lfc.width+term.space else term.width<-term.width*(-1)+term.space
  

  if (clust.by=='logFC') distance <- stats::dist(chord[,dim(chord)[2]], method=metric)
  if (clust.by=='term') distance <- stats::dist(chord, method=metric)
  cluster <- stats::hclust(distance, method=clust)
  dendr <- dendro_data(cluster)
  y_range <- range(dendr$segments$y)
  x_pos <- data.frame(x=dendr$label$x, label=as.character(dendr$label$label))
  chord <- as.data.frame(chord)
  chord$label <- as.character(rownames(chord))
  all <- merge(x_pos, chord, by='label')
  all$label <- as.character(all$label)
  if (nlfc){
    lfc_rect <- all[,c(2, dim(all)[2])]
    for (l in 4:dim(data)[2]) lfc_rect <- cbind(lfc_rect, sapply(all$label, function(x) data[match(x, data$genes), l]))
    num <- dim(data)[2]-1
    tmp <- seq(lfc.space, lfc.width, length = num)
    lfc<-data.frame(x=numeric(),width=numeric(),space=numeric(),logFC=numeric())
    for (l in 1:(length(tmp)-1)){
      tmp_df<-data.frame(x=lfc_rect[,1],width=tmp[l+1],space=tmp[l],logFC=lfc_rect[,l+1])
      lfc<-rbind(lfc,tmp_df)
    }
  }else{
    lfc <- all[,c(2, dim(all)[2])]  
    lfc$space <- lfc.space
    lfc$width <- lfc.width
  }
  term <- all[,c(2:(length(process)+2))]
  color<-NULL;termx<-NULL;tspace<-NULL;twidth<-NULL
  for (row in 1:dim(term)[1]){
    idx <- which(term[row,-1] != 0)
    if(length(idx) != 0){
      termx<-c(termx,rep(term[row,1],length(idx)))
      color<-c(color,term.col[idx])
      tmp<-seq(term.space,term.width,length=length(idx)+1)
      tspace<-c(tspace,tmp[1:(length(tmp)-1)])
      twidth<-c(twidth,tmp[2:length(tmp)])
    }
  }
  tmp <- sapply(lfc$logFC, function(x) ifelse(x > lfc.max, lfc.max, x))
  logFC <- sapply(tmp, function(x) ifelse(x < lfc.min, lfc.min, x))
  lfc$logFC <- logFC
  term_rect <- data.frame(x = termx, width = twidth, space = tspace, col = color)
  legend <- data.frame(x = 1:length(process),label = process)

  ggplot()+
    geom_segment(data=segment(dendr), aes(x=x, y=y, xend=xend, yend=yend))+
    geom_rect(data=lfc,aes(xmin=x-0.5,xmax=x+0.5,ymin=width,ymax=space,fill=logFC))+
    scale_fill_gradient2('logFC', space = 'Lab', low=lfc.col[3],mid=lfc.col[2],high=lfc.col[1],guide=guide_colorbar(title.position='top',title.hjust=0.5),breaks=c(min(lfc$logFC),max(lfc$logFC)),labels=c(round(min(lfc$logFC)),round(max(lfc$logFC))))+
    geom_rect(data=term_rect,aes(xmin=x-0.5,xmax=x+0.5,ymin=width,ymax=space),fill=term_rect$col)+
    geom_point(data=legend,aes(x=x,y=0.1,size=factor(label,levels=label),shape=NA))+
    guides(size=guide_legend("GO Terms",ncol=4,byrow=T,override.aes=list(shape=22,fill=term.col,size = 8)))+
    coord_polar()+
    scale_y_reverse()+
    theme(legend.position='bottom',legend.background = element_rect(fill='transparent'),legend.box='horizontal',legend.direction='horizontal')+
    theme_blank  
    
}

#' 
#' @name GOChord
#' @title Displays the relationship between genes and terms.
#' @description The GOChord function generates a circularly composited overview 
#'   of selected/specific genes and their assigned processes or terms. More 
#'   generally, it joins genes and processes via ribbons in an intersection-like
#'   graph. The input can be generated with the \code{\link{chord_dat}} 
#'   function.
#' @param data The matrix represents the binary relation (1= is related to, 0= 
#'   is not related to) between a set of genes (rows) and processes (columns); a
#'   column for the logFC of the genes is optional
#' @param title The title (on top) of the plot
#' @param space The space between the chord segments of the plot
#' @param gene.order A character vector defining the order of the displayed gene
#'   labels
#' @param gene.size The size of the gene labels
#' @param gene.space The space between the gene labels and the segement of the 
#'   logFC
#' @param nlfc Defines the number of logFC columns (default=1)
#' @param lfc.col The fill color for the logFC specified in the following form: 
#'   c(color for low values, color for the mid point, color for the high values)
#' @param lfc.min Specifies the minimium value of the logFC scale (default = -3)
#' @param lfc.max Specifies the maximum value of the logFC scale (default = 3)
#' @param ribbon.col The background color of the ribbons
#' @param border.size Defines the size of the ribbon borders
#' @param process.label The size of the legend entries
#' @param limit A vector with two cutoff values (default= c(0,0)). The first 
#' value defines the minimum number of terms a gene has to be assigned to. The 
#' second the minimum number of genes assigned to a selected term.
#' @details The \code{gene.order} argument has three possible options: "logFC", 
#'   "alphabetical", "none", which are quite self- explanatory.
#'   
#'   Maybe the most important argument of the function is \code{nlfc}.If your 
#'   \code{data} does not contain a column of logFC values you have to set
#'   \code{nlfc = 0}. Differential expression analysis can be performed for
#'   multiple conditions and/or batches. Therefore, the data frame might contain
#'   more than one logFC value per gene. To adjust to this situation the
#'   \code{nlfc} argument is used as well. It is a numeric value and it defines
#'   the number of logFC columns of your \code{data}. The default is "1"
#'   assuming that most of the time only one contrast is considered.
#'   
#'   To represent the data more useful it might be necessary to reduce the 
#'   dimension of \code{data}. This can be achieved with \code{limit}. The first
#'   value of the vector defines the threshold for the minimum number of terms a
#'   gene has to be assigned to in order to be represented in the plot. Most of
#'   the time it is more meaningful to represent genes with various functions. A
#'   value of 3 excludes all genes with less than three term assignments. 
#'   Whereas the second value of the parameter restricts the number of terms 
#'   according to the number of assigned genes. All terms with a count smaller 
#'   or equal to the threshold are excluded.
#' @seealso \code{\link{chord_dat}}
#' @import ggplot2
#' @import grDevices
#' @examples
#' \dontrun{
#' # Load the included dataset
#' data(EC)
#' 
#' # Generating the binary matrix
#' chord<-chord_dat(circ,EC$genes,EC$process)
#' 
#' # Creating the chord plot
#' GOChord(chord)
#' 
#' # Excluding process with less than 5 assigned genes
#' GOChord(chord, limit = c(0,5))
#' 
#' # Creating the chord plot genes ordered by logFC and a different logFC color scale
#' GOChord(chord,space=0.02,gene.order='logFC',lfc.col=c('red','black','cyan'))
#' }
#' @export

GOChord <- function(data, title, space, gene.order, gene.size, gene.space, nlfc = 1, lfc.col, lfc.min, lfc.max, ribbon.col, border.size, process.label, limit){
  y <- id <- xpro <- ypro <- xgen <- ygen <- lx <- ly <- ID <- logFC <- NULL
  Ncol <- dim(data)[2]
  
  if (missing(title)) title <- ''
  if (missing(space)) space = 0
  if (missing(gene.order)) gene.order <- 'none'
  if (missing(gene.size)) gene.size <- 3
  if (missing(gene.space)) gene.space <- 0.2
  if (missing(lfc.col)) lfc.col <- c('brown1', 'azure', 'cornflowerblue')
  if (missing(lfc.min)) lfc.min <- -3
  if (missing(lfc.max)) lfc.max <- 3
  if (missing(border.size)) border.size <- 0.5
  if (missing (process.label)) process.label <- 11
  if (missing(limit)) limit <- c(0, 0)
  
  if (gene.order == 'logFC') data <- data[order(data[, Ncol], decreasing = T), ]
  if (gene.order == 'alphabetical') data <- data[order(rownames(data)), ]
  if (sum(!is.na(match(colnames(data), 'logFC'))) > 0){
    if (nlfc == 1){
      cdata <- check_chord(data[, 1:(Ncol - 1)], limit)
      lfc <- sapply(rownames(cdata), function(x) data[match(x,rownames(data)), Ncol])
    }else{
      cdata <- check_chord(data[, 1:(Ncol - nlfc)], limit)
      lfc <- sapply(rownames(cdata), function(x) data[, (Ncol - nlfc + 1)])
    }
  }else{
    cdata <- check_chord(data, limit)
    lfc <- 0
  }
  if (missing(ribbon.col)) colRib <- grDevices::rainbow(dim(cdata)[2]) else colRib <- ribbon.col
  nrib <- colSums(cdata)
  ngen <- rowSums(cdata)
  Ncol <- dim(cdata)[2]
  Nrow <- dim(cdata)[1]
  colRibb <- c()
  for (b in 1:length(nrib)) colRibb <- c(colRibb, rep(colRib[b], 202 * nrib[b]))
  r1 <- 1; r2 <- r1 + 0.1
  xmax <- c(); x <- 0
  for (r in 1:length(nrib)){
    perc <- nrib[r] / sum(nrib)
    xmax <- c(xmax, (pi * perc) - space)
    if (length(x) <= Ncol - 1) x <- c(x, x[r] + pi * perc)
  }
  xp <- c(); yp <- c()
  l <- 50
  for (s in 1:Ncol){
    xh <- seq(x[s], x[s] + xmax[s], length = l)
    xp <- c(xp, r1 * sin(x[s]), r1 * sin(xh), r1 * sin(x[s] + xmax[s]), r2 * sin(x[s] + xmax[s]), r2 * sin(rev(xh)), r2 * sin(x[s]))
    yp <- c(yp, r1 * cos(x[s]), r1 * cos(xh), r1 * cos(x[s] + xmax[s]), r2 * cos(x[s] + xmax[s]), r2 * cos(rev(xh)), r2 * cos(x[s]))
  }
  df_process <- data.frame(x = xp, y = yp, id = rep(c(1:Ncol), each = 4 + 2 * l))
  xp <- c(); yp <- c(); logs <- NULL
  x2 <- seq(0 - space, -pi - (-pi / Nrow) - space, length = Nrow)
  xmax2 <- rep(-pi / Nrow + space, length = Nrow)
  for (s in 1:Nrow){
    xh <- seq(x2[s], x2[s] + xmax2[s], length = l)
    if (nlfc <= 1){
      xp <- c(xp, (r1 + 0.05) * sin(x2[s]), (r1 + 0.05) * sin(xh), (r1 + 0.05) * sin(x2[s] + xmax2[s]), r2 * sin(x2[s] + xmax2[s]), r2 * sin(rev(xh)), r2 * sin(x2[s]))
      yp <- c(yp, (r1 + 0.05) * cos(x2[s]), (r1 + 0.05) * cos(xh), (r1 + 0.05) * cos(x2[s] + xmax2[s]), r2 * cos(x2[s] + xmax2[s]), r2 * cos(rev(xh)), r2 * cos(x2[s]))
    }else{
      tmp <- seq(r1, r2, length = nlfc + 1)
      for (t in 1:nlfc){
        logs <- c(logs, data[s, (dim(data)[2] + 1 - t)])
        xp <- c(xp, (tmp[t]) * sin(x2[s]), (tmp[t]) * sin(xh), (tmp[t]) * sin(x2[s] + xmax2[s]), tmp[t + 1] * sin(x2[s] + xmax2[s]), tmp[t + 1] * sin(rev(xh)), tmp[t + 1] * sin(x2[s]))
        yp <- c(yp, (tmp[t]) * cos(x2[s]), (tmp[t]) * cos(xh), (tmp[t]) * cos(x2[s] + xmax2[s]), tmp[t + 1] * cos(x2[s] + xmax2[s]), tmp[t + 1] * cos(rev(xh)), tmp[t + 1] * cos(x2[s]))
      }}}
  if(lfc[1] != 0){
    if (nlfc == 1){
      df_genes <- data.frame(x = xp, y = yp, id = rep(c(1:Nrow), each = 4 + 2 * l), logFC = rep(lfc, each = 4 + 2 * l))
    }else{
      df_genes <- data.frame(x = xp, y = yp, id = rep(c(1:(nlfc*Nrow)), each = 4 + 2 * l), logFC = rep(logs, each = 4 + 2 * l))  
    }
  }else{
    df_genes <- data.frame(x = xp, y = yp, id = rep(c(1:Nrow), each = 4 + 2 * l))
  }
  aseq <- seq(0, 180, length = length(x2)); angle <- c()
  for (o in aseq) if((o + 270) <= 360) angle <- c(angle, o + 270) else angle <- c(angle, o - 90)
  df_texg <- data.frame(xgen = (r1 + gene.space) * sin(x2 + xmax2/2),ygen = (r1 + gene.space) * cos(x2 + xmax2 / 2),labels = rownames(cdata), angle = angle)
  df_texp <- data.frame(xpro = (r1 + 0.15) * sin(x + xmax / 2),ypro = (r1 + 0.15) * cos(x + xmax / 2), labels = colnames(cdata), stringsAsFactors = FALSE)
  cols <- rep(colRib, each = 4 + 2 * l)
  x.end <- c(); y.end <- c(); processID <- c()
  for (gs in 1:length(x2)){
    val <- seq(x2[gs], x2[gs] + xmax2[gs], length = ngen[gs] + 1)
    pros <- which((cdata[gs, ] != 0) == T)
    for (v in 1:(length(val) - 1)){
      x.end <- c(x.end, sin(val[v]), sin(val[v + 1]))
      y.end <- c(y.end, cos(val[v]), cos(val[v + 1]))
      processID <- c(processID, rep(pros[v], 2))
    }
  }
  df_bezier <- data.frame(x.end = x.end, y.end = y.end, processID = processID)
  df_bezier <- df_bezier[order(df_bezier$processID,-df_bezier$y.end),]
  x.start <- c(); y.start <- c()
  for (rs in 1:length(x)){
    val<-seq(x[rs], x[rs] + xmax[rs], length = nrib[rs] + 1)
    for (v in 1:(length(val) - 1)){
      x.start <- c(x.start, sin(val[v]), sin(val[v + 1]))
      y.start <- c(y.start, cos(val[v]), cos(val[v + 1]))
    }
  }	
  df_bezier$x.start <- x.start
  df_bezier$y.start <- y.start
  df_path <- bezier(df_bezier, colRib)
  if(length(df_genes$logFC) != 0){
    tmp <- sapply(df_genes$logFC, function(x) ifelse(x > lfc.max, lfc.max, x))
    logFC <- sapply(tmp, function(x) ifelse(x < lfc.min, lfc.min, x))
    df_genes$logFC <- logFC
  }
  
  g<- ggplot() +
    geom_polygon(data = df_process, aes(x, y, group=id), fill='gray70', inherit.aes = F,color='black') +
    geom_polygon(data = df_process, aes(x, y, group=id), fill=cols, inherit.aes = F,alpha=0.6,color='black') +	
    geom_point(aes(x = xpro, y = ypro, size = factor(labels, levels = labels), shape = NA), data = df_texp) +
    guides(size = guide_legend("GO Terms", ncol = 4, byrow = T, override.aes = list(shape = 22, fill = unique(cols), size = 8))) +
    theme(legend.text = element_text(size = process.label)) +
    geom_text(aes(xgen, ygen, label = labels, angle = angle), data = df_texg, size = gene.size) +
    geom_polygon(aes(x = lx, y = ly, group = ID), data = df_path, fill = colRibb, color = 'black', size = border.size, inherit.aes = F) +		
    labs(title = title) +
    theme_blank
  
  if (nlfc >= 1){
    g + geom_polygon(data = df_genes, aes(x, y, group = id, fill = logFC), inherit.aes = F, color = 'black') +
      scale_fill_gradient2('logFC', space = 'Lab', low = lfc.col[3], mid = lfc.col[2], high = lfc.col[1], guide = guide_colorbar(title.position = "top", title.hjust = 0.5), 
                           breaks = c(min(df_genes$logFC), max(df_genes$logFC)), labels = c(round(min(df_genes$logFC)), round(max(df_genes$logFC)))) +
      theme(legend.position = 'bottom', legend.background = element_rect(fill = 'transparent'), legend.box = 'horizontal', legend.direction = 'horizontal')
  }else{
    g + geom_polygon(data = df_genes, aes(x, y, group = id), fill = 'gray50', inherit.aes = F, color = 'black')+
      theme(legend.position = 'bottom', legend.background = element_rect(fill = 'transparent'), legend.box = 'horizontal', legend.direction = 'horizontal')
  }
}