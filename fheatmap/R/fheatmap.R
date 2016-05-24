
# Graphics Settings
## Theme Blank
theme_blank <- function(){theme_blank<-theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank(), 
                                             panel.background = element_blank(), axis.line = element_blank(), axis.ticks = element_blank(), 
                                             axis.text= element_blank(),panel.grid=element_blank(),panel.border=element_blank(),
                                             plot.margin=rep(unit(0,"null"),4),panel.margin = unit(0,"null"),plot.margin=rep(unit(0,"null"),4),
                                             legend.position = "none")
                          return(theme_blank)}

##Get default colors or get "n" distinct ggplot colors
ggplotColours <- function(n=6, h=c(0, 360) +15){
  if ((diff(h)%%360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

##Scale / Z transform
Zscore <- function(x){
  transformedx <- (x-median(x))/(sd(x))
  return(transformedx)
}

## Read user given Filename 
read_file <- function(filename, header ,rowname, ...)
{ 
  r = regexpr("\\.[a-zA-Z]*$", filename)
  ending = substr(filename, r + 1, r + attr(r, "match.length"))
  f = switch(ending,
             xls  = function(x, ...) read.xls(x, header= header, ...),
             xlsx = function(x, ...) read.xls(x, header= header, ...),
             csv  = function(x, ...) read.csv(x, header= header, ...)
  )
  if(is.null(f)) { f = function(x, ...) read.delim(x, header= header, sep="\t", ...) }
  
  df<-f(filename)
  if (rowname ) { df_edit <- as.data.frame(df[,2:ncol(df)]) ; rownames(df_edit) <- df[,1] }
  else { df_edit <- df }
  
  return(df_edit)
}

#Generate breaks for the matrix
generate_breaks = function(x, n, center = F){
  if(center){
    m = max(abs(c(min(x, na.rm = T), max(x, na.rm = T))))
    res = seq(-m, m, length.out = n + 1)
  }
  else{ res = seq(min(x, na.rm = T), max(x, na.rm = T), length.out = n + 1) }
  
  return(res)
}

#Get color for matrix
mat_vec_colours = function(x, color = NULL, breaks = NULL){
  list_color<-color[as.numeric(cut(x, breaks = breaks, include.lowest = T))]
  names(list_color)<-x
  
  return(list_color)
}

##Calculate width and height
length_npc <- function(list_name=NULL,dim=NULL,font_size=5,fontface="plain"){
  gp = list(fontsize = font_size/0.352777778,fontface=fontface,lwd=font_size/0.352777778)
  convertheight <- function(x){
    width <- convertHeight(unit(1, "grobheight", textGrob(x, rot = 90,gp = do.call(gpar, gp))), "npc",valueOnly = FALSE )
    return(width)
  }
  convertwidth <- function(x){
    width <-  convertWidth( unit(1, "grobwidth", textGrob(x, gp = do.call(gpar, gp))), "npc",valueOnly = FALSE )
    return(width)
  }
  if (dim == "width"){ label_Width <-unlist(lapply(list_name, convertwidth) ) }
  else {  label_Width <- unlist(lapply(list_name,convertheight)) }
  
  return(label_Width)
}

## Layout pos
vplayout = function(x, y){
  return(viewport(layout.pos.row = x, layout.pos.col = y))
}
## layout grid
vpgrid = function(x, y){
  return(viewport(layout = grid.layout(nrow=x, ncol=y)))
}

#Clustering Rows and columns
cluster_data = function(data, distance, method){
  if(!(method %in% c("ward.D","ward.D2", "single", "complete", "average", "mcquitty", "median", "centroid"))){
    stop("clustering method has to one form the list: 'ward', 'single', 'complete', 'average', 'mcquitty', 'median' or 'centroid'.")
  }
  if(!(distance[1] %in% c("correlation", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")) & class(distance) != "dist"){
    print(!(distance[1] %in% c("correlation", "euclidean", "maximum", "manhattan", "canberra", "binary", "minkowski")) | class(distance) != "dist")
    stop("distance has to be a dissimilarity structure as produced by dist or one measure  form the list: 'correlation', 'euclidean', 'maximum', 'manhattan', 'canberra', 'binary', 'minkowski'")
  }
  if(distance[1] == "correlation"){
    d = as.dist(1 - cor(t(data)))
  }
  else{
    if(class(distance) == "dist"){ 
      d = distance
    }
    else{
      d = dist(data, method = distance)
    }
  }
  
  return(hclust(d, method = method))
}

#Draw title
draw_title <- function(title=title,font_size=10,colour="black",family="",fontface="plain",dim =NULL){
  pushViewport(vplayout(dim$title[1]:dim$title[2],dim$title[3]:dim$title[4]))
  draw_title<-ggplot()+annotate(geom="text", x=2, y=1, label=title,size=font_size,colour=colour,family=family,fontface=fontface)+theme_blank()+ xlab(NULL) + ylab(NULL)
  print(draw_title,newpage=FALSE)
  upViewport()
}

#Draw Dendogram
draw_tree = function(hc=NULL,dim=NULL,rows=T,cut=0){
  #Height relative to max height & #Convert height into "npc"
  height=hc$height / (max(hc$height))
  
  #merge data frame and height
  merge <- na.omit(data.frame(hc$merge,height)[hc$order,])
  merge_melt <- melt(merge,id.vars=c("height"))
  
  #Get single and group units
  singletons <- hc$merge[hc$merge < 0]
  singletons <- singletons[order(match(singletons,-(hc$order)))]
  groups <-  sort(hc$merge[hc$merge > 0])
  
  #Determine min and max coordinates of each singletone line
  if (rows){
    pushViewport(vplayout(dim$row_tree[1]:dim$row_tree[2],dim$row_tree[3]:dim$row_tree[4]))
    single_min <- matrix(nrow=length(singletons),ncol=3,c(singletons,rep(max(merge_melt$height),length(singletons)),rev(seq(1.5,length(singletons)+1,1))),byrow=F)   
    merge_singleton <- merge_melt[merge_melt$value %in% singletons,]
    single_max <- matrix(nrow=nrow(single_min),ncol=3,c(single_min[,1],max(merge_melt$height)-merge_singleton[order(match(merge_singleton$value,singletons)),]$height,single_min[,3]),byrow=F)
    xlimits = c(limits=c(-0.05,1))
    ylimits = c(1,length(hc$order)+1)
  }
  else {
    pushViewport(vplayout(dim$col_tree[1]:dim$col_tree[2],dim$col_tree[3]:dim$col_tree[4]))  
    single_min <- matrix(nrow=length(singletons),ncol=3,c(singletons,seq(1.5,length(singletons)+1,1),rep(0,length(singletons))),byrow=F)
    merge_singleton <- merge_melt[merge_melt$value %in% singletons,]
    single_max <- matrix(nrow=nrow(single_min),ncol=3,c(single_min[,1],single_min[,2],merge_singleton[order(match(merge_singleton$value,singletons)),]$height),byrow=F)
    xlimits = c(1,length(hc$order)+1)
    ylimits  =c(0,1.05)
  }
  
  #Populate the dimensions of group unit
  get_xy_value <- function(x){ 
    ymax <- single_max[single_max[,1]== x][3]
    xmax <- single_max[single_max[,1]== x][2]
    return(c(xmax,ymax))
  }
  get_y_recursive <- function(x,rows=T){  
    unit1 <- hc$merge[x,][1]
    unit2 <- hc$merge[x,][2]
    if ( !(unit1 %in% single_max[,1])) { dim1<-get_y_recursive(unit1) } else { dim1 <- get_xy_value(unit1) }
    if ( !(unit2 %in% single_max[,1])) { dim2<-get_y_recursive(unit2) } else { dim2 <- get_xy_value(unit2) }
    
    if (rows) {
      single_min <- c(x,dim1[1],(dim1[2]+dim2[2])/2)
      single_max <- c(x,max(merge_melt$height)-merge_melt[which(merge_melt$value == x),]$height,(dim1[2]+dim2[2])/2 )
    }
    else {
      single_min <- c(x,(dim1[1]+dim2[1])/2,dim1[2])
      single_max <- c(x,(dim1[1]+dim2[1])/2,merge_melt[which(merge_melt$value == x),]$height)
    }
    
    return(list(single_min,single_max))
  }
  for(i in 1:length(groups)){
    group_dim <- get_y_recursive(groups[i],rows=rows)
    single_min <- rbind(single_min,group_dim[[1]] )
    single_max <- rbind(single_max,group_dim[[2]] )
  }
  start <- max(single_max[,1])
  for(i in 1:nrow(hc$merge)){
    start=start+1
    dim1 <- get_xy_value(hc$merge[i,][1]) ; single_min <- rbind(single_min,c(start,dim1))
    dim2 <- get_xy_value(hc$merge[i,][2]) ; single_max <- rbind(single_max,c(start,dim2))
  }
  
  #making final coordinate data frame
  cond <- "cond" ;  xval<- "xval" ; yval <- "yval"
  coordinates <- as.data.frame(rbind(single_min,single_max))
  coordinates <-  coordinates[order(coordinates[,1]),]
  colnames(coordinates) <- c(cond,xval,yval)
  
  ##Draw the dendogram
  dendogram_plot <- ggplot(coordinates, aes(x=xval, y=yval, group = cond,color=as.factor(cond)))+
    geom_line(size=1) + 
    #geom_point() +
    scale_y_continuous(limits=ylimits,expand=c(0,0) ) + 
    scale_x_continuous(limits=xlimits,expand=c(0,0) ) +
    theme_blank() + xlab(NULL) + ylab(NULL) +
    theme(plot.margin=unit(c(0.05,0.05,0.05,0.05),"cm"))
  
  get_color_recursive <- function(x,color_vec=NULL){ 
    lines_color_dm <-na.omit( matrix(ncol=2,nrow=1) )
    
    unit1 <- hc$merge[x,][1] ; unit2 <- hc$merge[x,][2]
    
    if ( unit1 > 0 ){ 
      lines_color_dm <- rbind(lines_color_dm,c(unit1,color_vec[1]))
      x<-get_color_recursive(unit1,rep(color_vec[1],2)) 
      lines_color_dm <- rbind(x,lines_color_dm) 
    }  
    else {  lines_color_dm <- rbind(lines_color_dm,c(unit1,color_vec[1])) }
    
    if ( unit2 > 0 ){ 
      lines_color_dm <- rbind(lines_color_dm,c(unit2,color_vec[2]))
      y<-get_color_recursive(unit2,rep(color_vec[2],2)) 
      lines_color_dm <- rbind(y,lines_color_dm) 
    }  
    else {  lines_color_dm <- rbind(lines_color_dm,c(unit2,color_vec[2])) }  
    
    return(lines_color_dm)
  }
  
  make_color_dm <- function(hc,input_colors){
    dm <-get_color_recursive((nrow(hc$merge)+1)-j,color_vec=input_colors)
    return(dm)
  }
  
  #Check if row/col tree to color,if cut is mentioned
  if ( cut > 0 ) {
    colors <- ggplotColours(n=cut)
    row    <- nrow(hc$merge)
    dm     <- na.omit(matrix(ncol=2,nrow=1)) 
    for( j in 1:(cut-1)){ 
      input_colors <- colors[j:(j+1)]
      dm_ret <-make_color_dm(hc,input_colors)
      dm<-rbind(dm_ret,matrix(dm[!(dm[, 1] %in% dm_ret[,1])], ncol = ncol(dm)))
      dm[dm[,1]==c(hc$merge[row,][1])][2] <- "#000000"
      dm[dm[,1]==c(hc$merge[row,][2])][2] <- "#000000"
      row = row-1  
    }
    # Add connecting lines
    start <- max(as.numeric(dm[,1]))
    for(i in 1:nrow(hc$merge)){
      start=start+1
      dm<- rbind(dm,c(start,dm[dm[,1]==hc$merge[i][1]][2]))  
    }
    #Make color named vector
    color_vec<- dm[,2]
    names(color_vec)<-as.numeric(dm[,1])   
    dendogram_plot <- dendogram_plot + scale_color_manual(values = color_vec)
  }
  else { color_vec <- rep("black",nrow(coordinates)) 
         dendogram_plot <- dendogram_plot + scale_color_manual(values = color_vec)
  }
  
  print(dendogram_plot,newpage=FALSE)
  upViewport()
  
}

##Draw Dimensions
grid_dim_function <- function(display_tree_row=F,display_rownames=F,annotation_row=annotation_row,annotation_col=annotation_col,
                              display_tree_col=F,display_colnames=F){
  #Title
  ih_title=1;jh_title=1;iw_title=1;jw_title=6
  #matrix
  ih_mat=4 ; jh_mat=4 ; iw_mat=3 ; jw_mat=3 
  #Matrix legend
  ih_mlegend=4 ; jh_mlegend=4 ; iw_mlegend=5 ; jw_mlegend=5
  #Row_tree
  ih_rtree=4;jh_rtree=4;iw_rtree=1;jw_rtree=1
  #Col_tree
  ih_ctree=2;jh_ctree=2;iw_ctree=3;jw_ctree=3
  #row_names
  ih_rname=4;jh_rname=4;iw_rname=4;jw_rname=4 
  #col names
  ih_cname=5;jh_cname=5;iw_cname=3;jw_cname=3
  #row annotation
  ih_rannot=4;jh_rannot=4;iw_rannot=2;jw_rannot=2
  #colannotation
  ih_cannot=3;jh_cannot=3;iw_cannot=3;jw_cannot=3  
  #row & col annotation legend
  ih_annot_legend=4;jh_annot_legend=4;iw_annot_legend=6;jw_annot_legend=6
  
  #Matrix & Treerow & Annrow
  if (!display_tree_row){ 
    if(is.null(annotation_row)){ 
      iw_mat=iw_mat-2; iw_cannot=iw_mat; iw_cname=iw_mat;iw_ctree=iw_mat
    } 
    else{ iw_mat=iw_mat-1; iw_cannot=iw_mat; iw_cname=iw_mat;iw_ctree=iw_mat; iw_rannot = iw_rannot -1 ; jw_rannot=jw_rannot-1 } 
  }
  else { if(is.null(annotation_row)){
    iw_mat=iw_mat-1 ;iw_cannot=iw_mat ; iw_cname=iw_mat;iw_ctree=iw_mat
  } 
  }
  
  #matrix & Treecol & Anncol
  if (!display_tree_col){ 
    if(is.null(annotation_col)){ 
      ih_mat=ih_mat-2;ih_rannot=ih_mat;ih_rname=ih_mat;ih_rtree=ih_mat;ih_mlegend=ih_mat;
      ih_annot_legend=ih_mat
    } 
    else{ ih_mat=ih_mat-1;ih_rannot=ih_mat;ih_rname=ih_mat;ih_rtree=ih_mat;ih_mlegend=ih_mat;
          ih_annot_legend=ih_mat
          ih_cannot=ih_cannot-1;jh_cannot=jh_cannot-1
    } 
  }
  else { if(is.null(annotation_col)){ 
    ih_mat=ih_mat-1 ;ih_rannot=ih_mat ;ih_rname=ih_mat;ih_rtree=ih_mat;ih_mlegend=ih_mat
    ih_annot_legend=ih_mat
  } 
  }
  
  #Colname
  if (!display_colnames){ jh_mat =jh_mat+1;jh_rannot=jh_mat;jh_rtree=jh_mat;jh_rname=jh_mat}
  
  #rownames
  if (!display_rownames){ jw_mat =jw_mat+1;jw_cannot=jw_mat;jw_ctree=jw_mat;jw_cname=jw_mat }
  
  #Legend
  if (is.null(annotation_col) & is.null(annotation_row)){
    iw_mlegend = iw_mlegend+1 ;jw_mlegend = jw_mlegend+1
    iw_rname=iw_rname+1;jw_rname=jw_rname+1
    jw_mat=jw_mat+1 ; jw_cname=jw_mat ; jw_ctree = jw_mat
  }
  
  #dim_list
  dim_list =list(title=c(ih_title,jh_title,iw_title,jw_title),
                 matrix=c(ih_mat,jh_mat,iw_mat,jw_mat),
                 mat_legend=c(ih_mlegend, jh_mlegend, iw_mlegend, jw_mlegend),
                 row_tree=c(ih_rtree,jh_rtree,iw_rtree,jw_rtree),
                 col_tree=c(ih_ctree,jh_ctree,iw_ctree,jw_ctree),
                 row_name=c(ih_rname,jh_rname,iw_rname,jw_rname),
                 col_name=c(ih_cname,jh_cname,iw_cname,jw_cname),
                 row_annot=c(ih_rannot,jh_rannot,iw_rannot,jw_rannot),
                 col_annot=c(ih_cannot,jh_cannot,iw_cannot,jw_cannot),
                 annot_legend=c(ih_annot_legend,jh_annot_legend,iw_annot_legend,jw_annot_legend)             
  )
  return(dim_list)
}
#Draw Annotation plots                         
##Push Viewports
make_viewports_grid <- function(annotation_row=NULL,annot_row_color=NULL,annotation_col=NULL,annot_col_color=NULL,annotation_palette=NULL, 
                                npalette_col=5,col_str_length=NULL,row_str_length=NULL,dim=NULL,legend_fontsize=2,col_names=NULL,row_names=NULL,
                                family="",legend_fontface="",legend_color="black", ...){
  
  row_legend_list <-list() ;col_legend_list <-list() ; legend_list <- list();
  legend_row = 0 ; legend_col =0 ; total_row = 0
  
  if(!(is.null(annotation_row ))){
    legend_row=length(annotation_row)
    pushViewport(vplayout(dim$row_annot[1]:dim$row_annot[2],dim$row_annot[3]:dim$row_annot[4]))
    pushViewport(vpgrid(1,ncol(annotation_row)))
    row_legend_list<- generate_annot_plot(annot_object=annotation_row,row=T, annot_color=annot_row_color,names_order=row_names,annot_palette=annotation_palette,npalette_col=npalette_col, ...)
    upViewport(2) 
    
    legend_list <- row_legend_list
  }
  
  if(!(is.null(annotation_col))){
    legend_col=length(annotation_col)
    pushViewport(vplayout(dim$col_annot[1]:dim$col_annot[2],dim$col_annot[3]:dim$col_annot[4]))
    pushViewport(vpgrid(ncol(annotation_col),1))
    col_legend_list<- generate_annot_plot(annot_object=annotation_col,row=F, annot_color=annot_col_color,names_order=col_names,annot_palette=annotation_palette,npalette_col=npalette_col, ...)
    upViewport(2)  
    
    legend_list <- c(legend_list,col_legend_list)
  }
  
  total_row <- sum(c(legend_row,legend_col))
  
  if (total_row > 0){
    generate_legend_plot(dim=dim,fill_legend=legend_list,legend_fontsize=legend_fontsize,family=family,legend_fontface=legend_fontface,legend_color=legend_color)
  }
}
## Generate Annotation object
generate_annot_plot <- function(annot_object=NULL,annot_color=NULL,npalette_col=5,annot_palette=NULL,names_order=NULL,row=F,seed=13,...){
  
  set.seed(13)
  annot_object<- as.data.frame(annot_object[rev(names_order),])
  fill_legend_list <- list()
  for (j in 1:ncol(annot_object) ) {
    color_count = length(levels(annot_object[,j]))+50
    levels_count  = length(levels(annot_object[,j]))
    
    if (row){ 
      pushViewport(vplayout(1,j)) 
      #df_row_annot <- data.frame(
      x_row_annot<-rep(1,nrow(annot_object))
      y_row_annot<-seq(1,nrow(annot_object),1)
      #)
    }
    else {
      pushViewport(vplayout(j,1)) 
      #df_row_annot <- data.frame(
      y_row_annot<-rep(1,nrow(annot_object))
      x_row_annot<-seq(1,nrow(annot_object),1)
      #)
    }
    
    if(!is.integer(annot_object[,j])){
      if(!(is.null(annot_color))){
        fill= as.character(as.data.frame(annot_color[rev(names_order),])[,j]) ;
        names(fill) <-  as.character(annot_object[,j])
      }
      else {
        if(!is.null(annot_palette)){ 
          factor_color=sample( colorRampPalette(brewer.pal(npalette_col,annot_palette))(levels_count),levels_count,replace=F)
          names(factor_color)=levels(annot_object[,j]);fill= factor_color[annot_object[,j]]
        }
        else {
          factor_color = sample( ggplotColours(color_count),levels_count,replace=F )
          names(factor_color) <- levels(annot_object[,j]); fill= factor_color[annot_object[,j]]  
        }
      }
      fill_legend_list[paste("string",j,row,sep="")] <- list(fill[unique(names(fill))]) 
    }
    else{
      factor_color<- colorRampPalette(brewer.pal(9,"Blues")[3:9])(length(unique(annot_object[,j]))+2)
      names(factor_color) <- sort(unique(annot_object[,j])) ; fill= factor_color[as.character(annot_object[,j])] 
      
      #Make list for Legend
      mat_color <- colorRampPalette(brewer.pal(9,"Blues")[3:9])(100)
      pretty_range <- grid.pretty(range(as.matrix(annot_object[,j],ncol=1), na.rm = TRUE))
      labels <-rep("",length(mat_color))
      
      y<- (round((length(labels)-1)/length(pretty_range)))
      labels[seq(y,length(labels)-y,length.out= length(pretty_range)-1)] <- pretty_range[-length(pretty_range)]
      labels[length(labels)] <- pretty_range[length(pretty_range)]
      names(mat_color) <- labels
      fill_legend_list[paste("numeric",j,row,sep="")] <- list(mat_color)
    }
    
    df_row_annot <- data.frame(x_row_annot,y_row_annot)
    print(ggplot(df_row_annot, mapping=aes(xmin=x_row_annot, xmax=x_row_annot+1, ymin=y_row_annot,ymax=y_row_annot+1 ) )+ geom_rect(fill=fill) + theme_blank() +
            xlab(NULL) + ylab(NULL) +scale_x_discrete(expand=c(0,0)) + 
            scale_y_discrete(expand=c(0,0)) +
            theme(plot.margin=unit(c(0.05,0.05,0.05,0.05),"cm")),newpage = FALSE)
    upViewport()
  }
  
  return(fill_legend_list)
}

#Generate Legend Plot
generate_legend_plot <- function(dim=NULL,fill_legend=NULL,legend_fontsize=4,family="",legend_fontface="",legend_color="black", ...){
  
  pushViewport(vplayout(dim$annot_legend[1]:dim$annot_legend[2],dim$annot_legend[3]:dim$annot_legend[4]))
  
  index_num_ele=NULL ; index_str_ele=NULL ; num_height = NULL
  width_rect  = max(unlist(length_npc(list_name=c("XYZ"),dim="width",font_size=legend_fontsize)))
  width_text = max(unlist(length_npc(list_name=fill_legend,dim="width",font_size=legend_fontsize)))
  height_text = max(unlist(length_npc(list_name=c("X"),dim="height",font_size=legend_fontsize)))
  padding = 0.05
  
  nplots <- length(fill_legend)
  height_list = c(1:nplots)
  index_str_ele = grep( "string\\d+",names(fill_legend),perl=TRUE, value=FALSE )
  index_num_ele = grep( "numeric\\d+",names(fill_legend),perl=TRUE, value=FALSE )
  
  str_height <- unname(unlist(lapply(fill_legend[index_str_ele],length)))*(2*height_text)
  if(length(index_num_ele) > 0){ 
    num_height <- rep((1-sum(str_height))/length(index_num_ele),length(index_num_ele)) 
  }
  
  height_list[c(index_str_ele,index_num_ele)] <- as.numeric(c(str_height,num_height))
  unused_ht <-  1-sum(height_list)
  
  pushViewport(viewport(layout = grid.layout( nplots+1, 1 ,heights= unit.c(unit(c(height_list,unused_ht),"npc")))))
  
  for (i in 1:length(fill_legend)){
    pushViewport(vplayout(i,1))
    ylimit  <- length(fill_legend[[i]])
    fill=unname(fill_legend[[i]])
    labels=names(fill_legend[[i]])
    
    xmin    <- rep(padding+0.05,length(fill_legend[[i]]))
    xmax    <- xmin + width_rect
    text_x <- xmax + padding
    scalex <- scale_x_continuous(limits=c(padding+0.01,1.5+width_text),expand=c(0,0))
    ymin = seq(0.1,0.9,length.out =ylimit )
    ymax = ymin + (ymin[2]-ymin[1]-0.05)
    text_y <- (ymin+ymax)/2
    scaley <- scale_y_continuous(limits=c(0,max(ymax)+0.05),expand=c(0,0))
    
    df <- data.frame(xmin,xmax,ymin,ymax,text_x,text_y)
    p<-ggplot(df,mapping=aes(xmin = xmin, xmax  = xmax, ymin = ymin , ymax = ymax ))+
      geom_rect(fill=unname(fill_legend[[i]])) +
      geom_text(x=text_x, y=text_y, label= labels,size=legend_fontsize,family = family,fontface=legend_fontface,colour=legend_color,hjust=0,vjust=0.5)+
      theme_blank() + scalex + scaley 
    print(p,newpage=F)
    upViewport()
  } 
  upViewport(2)
}  

## Generate Colnames and Rownames
draw_names <- function(col_names=NULL,row_names=NULL,total_names=0,font_size=0.5,family="",names_fontface="",names_color="black",dim=NULL){
  
  if(!(is.null(col_names))){
    
    pushViewport(vplayout(dim$col_name[1]:dim$col_name[2],dim$col_name[3]:dim$col_name[4]))
    x_coltext <- seq(1.5,total_names+1,1)
    y_coltext <- rep(max(nchar(col_names))-0.5,total_names)
    df_test <- data.frame(x_coltext,y_coltext,col_names)
    ## Edited the missing column name
    names_plot <- ggplot(df_test)  +
      scale_x_continuous(limits=c(1, total_names + 1),expand=c(0,0) ) + scale_y_continuous( limits=c(1, max(nchar(col_names))) ,expand=c(0,0)) +
      geom_text(x= x_coltext ,y=y_coltext,label=col_names,angle = 90,size=font_size,family=family,fontface= names_fontface, colour= names_color,hjust=1) +
      theme_blank()  + theme(plot.margin=unit(c(0.05,0.05,0.05,0.05),"cm")) +
      xlab(NULL) + ylab(NULL)
    print(names_plot,newpage=FALSE)
    upViewport()
  }
  
  if(!(is.null(row_names))){
    pushViewport(vplayout(dim$row_name[1]:dim$row_name[2],dim$row_name[3]:dim$row_name[4]))
    y_coltext <- rev(seq(1.5,total_names+1,1))
    x_coltext <- rep(1.5,total_names)
    df_test <- data.frame(x_coltext,y_coltext,row_names)
    names_plot <- ggplot(df_test)  +
      scale_y_continuous(limits=c(1, total_names+1),expand=c(0,0) ) + scale_x_continuous( limits=c(1, max(nchar(row_names))) ,expand=c(0,0)) +
      geom_text(x= x_coltext ,y=y_coltext,label=row_names,size=font_size,family=family, fontface= names_fontface, colour= names_color,hjust=0) +
      theme_blank()  + theme(plot.margin=unit(c(0.05,0.05,0.05,0.05),"cm")) +
      xlab(NULL) + ylab(NULL)
    print(names_plot,newpage=FALSE)
    upViewport()
  }
}

#Draw heatmap matrix 
draw_mat <- function(data,dim=NULL,breaks=NULL,mat_color=NULL,cell_border=NULL,cell_border_col=NULL,display_number=F){
  
  #Make space around the matrix
  xmin <- rep(1:ncol(data),nrow(data))
  xmax <- xmin+1
  ymin <- rep(((nrow(data):1)),each=ncol(data))
  ymax <- ymin+1
  df   <- data.frame(xmin,xmax,ymin,ymax) 
  text <- unlist(unname(as.list((t(data)))))
  x_axis_breaks <- ((xmin+xmax)/2)
  y_axis_breaks <- ((ymin+ymax)/2)
  
  if (is.null(breaks)) {
    breaks = generate_breaks(as.vector(as.matrix(data)),length(mat_color))
  }
  else { breaks = breaks }
  
  mat_select_color = mat_vec_colours(as.numeric(as.matrix(data)), color = mat_color, breaks = breaks)
  ##Edited : Added line to each rectangle
  p <-ggplot(df, mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin,ymax=ymax ))
  
  if(cell_border== T){ p <- p+geom_rect(color=cell_border_col,fill=mat_select_color[as.character(t(data))]) } 
  else { p <- p+geom_rect(fill=mat_select_color[as.character(t(data))]) }
  
  if(display_number == T){  p <- p+geom_text(x=x_axis_breaks, y=y_axis_breaks, label=as.character(text),size=3) }
  
  p <- p+theme_blank()  + xlab(NULL) + ylab(NULL) + scale_x_discrete(expand=c(0,0)) + scale_y_discrete(expand=c(0,0)) +
        theme(plot.margin=unit(c(0.05,0.05,0.05,0.05),"cm"))
  
  pushViewport(vplayout(dim$matrix[1]:dim$matrix[2],dim$matrix[3]:dim$matrix[4]))
  print(p,newpage=FALSE)
  upViewport()
  return(mat_select_color)
  
}

#Draw heatmap legend
draw_mat_legend <- function(data,breaks=NULL,mat_color,mat_legend_size=5,dim=NULL){
  
  pretty_range<-grid.pretty(range(as.matrix(data), na.rm = TRUE))
  pretty_range[1] <- min(as.matrix(data))
  xmin <- rep(0.1,length(mat_color))
  ##pretty_range[1] <- min(as.matrix(data))
  xmin <- rep(0.5,length(mat_color))
  xmax <- xmin+0.5
  ymin = seq(5,10,length.out =length(mat_color) )
  ymax = ymin + (ymin[2]-ymin[1])
  df_mat_legend <- data.frame(xmin,xmax,ymin,ymax)
  x_axis_breaks <- rep(1.3,length(pretty_range))
  y_axis_breaks <- seq(5.5,max(ymax)-0.5,length.out =length(pretty_range))
  
  p<-ggplot(df_mat_legend, mapping=aes(xmin=xmin, xmax=xmax, ymin=ymin,ymax=ymax ))+
    geom_rect(fill=mat_color)+annotate(geom="text",x=x_axis_breaks,y=y_axis_breaks,label=pretty_range,size=mat_legend_size)+ 
    xlab(NULL) + ylab(NULL)+
    scale_x_continuous(limits=c(0, 2),expand=c(0,0)) +scale_y_continuous(limits=c(0, max(ymax)+0.6),expand=c(0,0))+
    theme_blank()
  
  pushViewport(vplayout(dim$mat_legend[1]:dim$mat_legend[2],dim$mat_legend[3]:dim$mat_legend[4]))
  print(p,newpage=FALSE)
  
  upViewport()
}

#' A function to draw clustered heatmaps.
#' 
#' A function to draw clustered heatmaps where one has better control over some graphical 
#' parameters such as color of heatmap, dendograms and annotations of rows and columns etc. 
#' 
#' @param data Data, file name or numeric dataframe/matrix to plot a heatmap.For example "\code{\link{fheatmap_data}}".
#' @param header Logical, determining if the data file contains the names of the variables as its first line. (Default: TRUE)
#' @param scale To standardize the given dataframe/matrix. (Default: FALSE)
#' @param title Title of the heatmap. (Default: NA)
#' @param title_fontsize Fontsize of the title. (Default: 6)
#' @param title_color Color of the title. (Default: black)
#' @param title_font_style Fontstyle of the title. 
#' @param title_fontface Fontface of the title. (Default: plain)
#' @param rowname Logical determining if the file contains row names  as its first column. (Default: TRUE)
#' @param breaks A sequence of numbers that covers the range of values in matrix and is one 
#' element longer than mat_color vector. It is useful to assign different colors to different set of values.
#' If value is NULL then the breaks are calculated automatically.
#' @param mat_color Vector of colors used in heatmap. If value is NULL colors are calculated automatically. 
#' (Default: "green", "yellow","red")
#' @param cluster_rows Logical determining if clustering in rows to be done or not. (Default: FALSE)
#' @param cluster_cols Logical determining if clustering in columns to be done or not. (Default: TRUE)
#' @param cut_rowtree  Numeric value(N) to color N clusters in row dendogram. (Default: 0)
#' @param cut_coltree  Numeric value(N) to color N clusters in column dendogram. (Default: 0)
#' @param display_tree_col Logical determining if column dendogram to be shown. (Default: FALSE)
#' @param display_tree_row Logical determining if row dendogram to be shown. (Default: FALSE)
#' @param cluster_distance_rows Distance measure used in clustering rows. Possible values are \code{"correlation"}
#' for Pearson correlation and all the distances supported by \code{\link{dist}}, such as \code{"euclidean"}, etc. (Default: "Euclidean")
#' @param cluster_distance_cols Distance measure used in clustering columns. Possible values the same as for 
#' clustering_distance_rows. (Default: "Euclidean")
#' @param clustering_method clustering method used. For Clustering methods refer the same values as \code{\link{hclust}}. (Default: "ward.D")
#' @param annotation_row Matrix or file with annotations of rows,for example "\code{\link{annotation_row}}" .Each column defines the features 
#' for a specific row. The rows in the data and columns in the annotation are matched using corresponding row and column names.
#'  Note that color schemes takes into account if variable is continuous or discrete. If NULL then annotation for row is not drawn.(Default: NULL)
#' @param annotation_col Matrix or file with annotations of columns,for example "\code{\link{annotation_col}}" . Each row defines the features for a
#' specific column. The columns in the data and rows in the annotation are matched using corresponding row and column names. Note that color schemes takes into 
#' account if variable is continuous or discrete. If NULL then annotation for columns is not drawn.(Default: NULL)
#' @param annot_row_color Matrix or file with colors for annotation in row, for example "\code{\link{annotation_row_color}}" . If a file is given 
#' then row names  and column names  of annotation in row should be given in the first row & column respectively. If NULL then default colors are printed.
#' (Default: NULL)
#' @param annot_col_color Matrix or file with colors for annotation in column,for example "\code{\link{annotation_col_color}}". If a file is given then 
#' row names  and column names  of annotation in column should be given in the first row & column respectively. If NULL then default colors are printed. 
#' (Default: NULL)
#' @param annotation_palette RColorBrewer palette for colors of rows & column annotations.(Default: NULL)
#' Distinct color is selected for each feature in anotation. Use "seed" parameter to fix random selection of colors. If NULL then color
#'  matrix or file is expected to color the annotation. Either  annotation_palette or file can be given at one time (Default: NULL)
#' @param npalette_col Numeric(N) If "annotation_palette" given ,then N distinct colors are selected from annotation_palette.
#' Minimum 3, maximum depending on palette. (Default: 5)
#' @param seed  Numeric(N) If "annotation_palette" given ,then random number N to fix random selection of colors from the given palette.
#' @param display_colnames Logical determining if column names  to be shown or not.  (Default: TRUE)
#' @param display_rownames Logical determining if row names  to be shown or not.  (Default: TRUE)
#' @param fontsize Base Fontsize for the heatmap. (Default: 4)
#' @param row_fontsize Fontsize for row names . (Default: fontsize)
#' @param col_fontsize Fontsize for column names . (Default: fontsize)
#' @param names_font_style Fontstyle for row names  and column names .
#' @param names_fontface Fontface for row names  and column names . (Default: plain)
#' @param names_color Color for row names  and column names .
#' @param legend_fontsize Fontsize for legend annotation. (Default: fontsize)
#' @param legend_font_style Fontstyle for legend annotations.
#' @param legend_fontface Fontface for legend annotations.
#' @param legend_color Color for legend annotations.
#' @param mat_legend_size Fontsize for matrix legend. (Default: fontsize)
#' @param display_number Logical determining if the numeric values are also printed to the cells.
#' 
#' @author  Vaishali Tumulu and Sivasish Sindiri
#' 
#' @examples
#' 
#' # Draw heatmaps
#' fheatmap(fheatmap_data, title="Example Heatmap", title_fontsize=15, title_color="red", title_font_style="mono", 
#'          title_fontface="italic")
#' fheatmap(fheatmap_data, scale =T )
#' fheatmap(fheatmap_data, display_number = T)
#' fheatmap(fheatmap_data, mat_color=c("lightblue","orange","maroon"))
#' fheatmap(fheatmap_data, mat_legend_size=5)
#' fheatmap(fheatmap_data, names_font_style="Courier", names_fontface="italic", names_color="brown", fontsize=6)
#' 
#' 
#' #Draw Dendograms
#' fheatmap(fheatmap_data, cluster_rows=T)
#' fheatmap(fheatmap_data, cluster_rows=T, cut_rowtree=2, cut_coltree=2, display_tree_row=T)
#' fheatmap(fheatmap_data, cluster_rows=T, cluster_distance_rows="maximum", cluster_distance_cols="maximum",
#'          clustering_method= "complete")
#' 
#' # Generate annotations
#' fheatmap(fheatmap_data, annotation_row=annotation_row, annotation_col=annotation_col)
#' fheatmap(fheatmap_data, annotation_row=annotation_row, annotation_col=annotation_col, annot_row_color=annotation_row_color, annot_col_color=annotation_col_color)
#' fheatmap(fheatmap_data, annotation_row=annotation_row, annotation_col=annotation_col, annotation_palette = "Dark2", npalette_col=5, seed=3)
#' fheatmap(fheatmap_data, annotation_row=annotation_row, annotation_col=annotation_col, legend_fontsize=5, legend_font_style="Courier", legend_fontface="bold", legend_color="red")
#' 
#' 
#' 
#' @name fheatmap-package
#' @docType package
#' 
#' @export
#' @import grid
#' @import RColorBrewer
#' @import gdata
#' @import ggplot2
#' @import reshape2
#' @import gplots

#Main Function
fheatmap <- function(data,header=T,scale=F,title=NA,title_fontsize=6,title_color="black",title_font_style="",title_fontface="plain",
                     rowname=T,breaks=NULL,mat_color=NULL,cell_border=T,cell_border_col="slategrey",cluster_rows=F,cluster_cols=T,cut_rowtree=0,cut_coltree=0,
                     display_tree_col=T,display_tree_row=F, cluster_distance_rows="euclidean",cluster_distance_cols="euclidean",
                     clustering_method="ward.D",annotation_palette=NULL,npalette_col=5,annotation_row=NULL,annotation_col=NULL,annot_row_color=NULL,
                     annot_col_color=NULL,display_colnames=T,display_rownames=T,fontsize=4,row_fontsize=fontsize,col_fontsize=fontsize,
                     names_font_style="",names_fontface="plain",names_color="black",legend_fontsize=fontsize,legend_font_style="",legend_fontface="plain",
                     legend_color="black",mat_legend_size=fontsize,display_number=F,seed=13, ...){
  
  #Make new page
  grid.newpage()
  
  #Checking conditions
  if(!(is.null(annot_col_color)) || !(is.null(annot_row_color)) ){
    if(!(is.null(annotation_palette))){ 
      stop(sprintf("Select annotation_palette or annot_row_color/annot_col_color, but not both")) }
  }
  
  #Check if data is a matrix, data frame or string
  if (class(data) == "character"){ 
    tryCatch({ 
      data <- read_file(data,header= header,rowname=rowname) },
      error = function(err) { print(err) ;  stop(sprintf("Improper data file provided !!"))  } )
  }
  else { data = as.matrix(data)  } 
  if (scale){ data <- t(apply(data,1,Zscore)) }
  
  #Setting Variables
  #Get col and row names
  col_names<- colnames(data)
  row_names<- rownames(data)
  
  #Set the default parameters
  padding_w <- max(unlist(length_npc(list_name=c("XY"), dim="width",font_size=row_fontsize))) 
  padding_h <- max(unlist(length_npc(list_name=c("XY"), dim="HEIGHT",font_size=col_fontsize))) 
  max_row_annot_str <- c()
  max_col_annot_str <- c()
  
  #Determine Title Dimensions
  if (!is.na(title) ) { title_height = convertHeight(unit(1, "grobheight", textGrob(title,gp = do.call(gpar,list(fontsize = title_fontsize/0.352777778)))), "npc",valueOnly = TRUE )+padding_h }
  else { title_height = 0.01 }
  
  #Set default Widths
  legend_row_width     = 0.01
  legend_col_width     = 0.01
  mat_legend_width = padding_w*3 
  row_name_width   = max(unlist(length_npc(list_name=row_names, dim="width",font_size=row_fontsize))) +padding_w
  matrix_width     = (1-(legend_row_width+mat_legend_width+row_name_width))*0.98
  ann_row_width    = (1-(legend_row_width+mat_legend_width+row_name_width))*0.01
  tree_row_width   = (1-(row_name_width+legend_row_width+mat_legend_width))*0.01
  
  #Set default Heights
  col_name_height  = max(unlist(length_npc(list_name=col_names, dim="height",font_size=col_fontsize))) +padding_h
  matrix_height    = (1-(col_name_height+title_height))*0.98
  ann_col_height   = (1-(col_name_height+title_height))*0.01
  tree_col_height  = (1-(col_name_height+title_height))*0.01        
  
  #Read Annotation Block
  #Checking annotaion_row
  if(!(is.null(annotation_row))){
    if (class(annotation_row) == "character"){ 
      tryCatch({ 
        annotation_row <- read_file(annotation_row,header= T,rowname=T) 
      },
      error = function(err) { print(err) ;  stop(sprintf("Improper annot_row_color file provided !!"))  } )
    }
    else { annotation_row = annotation_row  }
    
    for( i in 1:ncol(annotation_row)){ max_row_annot_str[i]<-as.character(annotation_row[,i][which.max(nchar(as.character(annotation_row[,i])))]) }
    ncols_rowannot   = ncol(annotation_row)
    ann_row_width    = ncols_rowannot*(padding_w*1.5)
    legend_row_width= max(unlist(length_npc(list_name=max_row_annot_str,dim="width",font_size=legend_fontsize))) +padding_w*2.5  
  } 
  # Checking annotation_row colors
  if(!(is.null(annot_row_color))){
    if (class(annot_row_color) == "character"){  
      tryCatch({ 
        annot_row_color <- read_file(annot_row_color,header= T,rowname=T) },
        error = function(err) { print(err) ;  stop(sprintf("Improper annot_row_color file provided !!"))  } )
    }
    else { annot_row_color = annot_row_color  }  
  }
  
  #Checking annotaion_col
  if(!(is.null(annotation_col))){
    if (class(annotation_col) == "character"){ 
      tryCatch({ 
        annotation_col <- read_file(annotation_col,header= T,rowname=T) 
      },
      error = function(err) { print(err) ;  stop(sprintf("Improper annotaion_col file provided !!"))  } )
    }
    else { annotation_col = annotation_col  } 
    
    for( i in 1:ncol(annotation_col)){ max_col_annot_str[i]<-as.character(annotation_col[,i][which.max(nchar(as.character(annotation_col[,i])))]) }
    ncols_colannot = ncol(annotation_col)
    ann_col_height = ncols_colannot*(padding_h*1.5)
    legend_col_width= max(unlist(length_npc(list_name=max_col_annot_str,dim="width",font_size=legend_fontsize))) +padding_w*2.5
  }
  
  # Checking annotation_col colors
  if(!(is.null(annot_col_color))){
    if (class(annot_col_color) == "character"){  
      tryCatch({
        annot_col_color <- read_file(annot_col_color,header= T,rowname=T) },
        error = function(err) { print(err) ;  stop(sprintf("Improper annot_col_color file provided !!"))  } )
    }
    else { annot_col_color = annot_col_color  }  
  }
  
  #Check if clustering of row or column given
  # Do clustering
  if(cluster_rows){
    tree_row = cluster_data(data, distance = cluster_distance_rows, method = clustering_method)
    data = data[tree_row$order, , drop = FALSE]
    row_names<- rownames(data)
  }
  else { display_tree_row = F }
  if(cluster_cols){
    tree_col = cluster_data(t(data), distance = cluster_distance_cols, method = clustering_method)
    data = data[, tree_col$order, drop = FALSE]
    col_names<- colnames(data)
  }
  else { display_tree_col = F }
  
  #ReSet Parameters
  #ReSet Legend width
  if( !is.null(annotation_row) || !is.null(annotation_col) ){
    legend_width = c(legend_row_width,legend_col_width)[which.max(c(legend_row_width,legend_col_width))]
  }
  else{
    legend_width     = mat_legend_width
    mat_legend_width = row_name_width
    row_name_width   = 0.01
  }
  
  #Reset Width
  if (display_tree_row){ 
    tree_row_width = (1-(legend_width+mat_legend_width+row_name_width+ann_row_width))*0.2
    matrix_width   = (1-(legend_width+mat_legend_width+row_name_width+ann_row_width))*0.8
  }
  else{
    tree_row_width = ann_row_width
    matrix_width  = 1-(legend_width+mat_legend_width+row_name_width+ann_row_width+tree_row_width)
  }
  #Reset Height
  if (display_tree_col){
    tree_col_height  = (1-(title_height+ann_col_height+col_name_height))*0.2
    matrix_height   =  (1-(title_height+ann_col_height+col_name_height))*0.8
  }
  else{
    tree_col_height = ann_col_height
    matrix_height   = (1-(title_height+tree_col_height+ann_col_height+col_name_height))
  }
  
  #Make plot Grid
  pushViewport(plotViewport(c(0.1,1,0.1,0.1),layout=grid.layout(5,6, widths=unit.c(unit(c(tree_row_width,ann_row_width,matrix_width,row_name_width,mat_legend_width,legend_width),"npc")),
                                                                heights= unit.c(unit(c(title_height,tree_col_height,ann_col_height,matrix_height,col_name_height),"npc"))  )))
  
  
#   ##Draw grid for reference
#   for( i in 1:5)
#   {for (j in 1:6){
#     pushViewport(vplayout(i,j))
#     grid.rect()
#     upViewport()
#   }}
  
  #make Grid Dimensions
  grid_dim_list<-grid_dim_function(display_rownames=display_rownames,display_colnames=display_colnames,annotation_row=annotation_row,
                                   annotation_col=annotation_col,display_tree_col=display_tree_col,display_tree_row=display_tree_row)
  
  #   #Draw Title
  if( !is.na(title)){ 
    draw_title(title=title, font_size=title_fontsize,colour=title_color,family=title_font_style,fontface=title_fontface,dim=grid_dim_list)    
  }
  
  #Draw tree for the rows
  if( display_tree_row ){
    if(!is.na(tree_row[[1]][1])){ 
      draw_tree(hc=tree_row,dim=grid_dim_list,rows=T,cut=cut_rowtree)
    }
  }
  #Draw tree for the cols
  if( display_tree_col ){
    if(!is.na(tree_col[[1]][1])){ 
      draw_tree(hc=tree_col,dim=grid_dim_list,rows=F,cut=cut_coltree)
    }
  }
  
  ##Draw annotations
  make_viewports_grid(annotation_row=annotation_row, annot_row_color=annot_row_color,annotation_col=annotation_col, 
                      annot_col_color=annot_col_color,annotation_palette=annotation_palette,npalette_col=npalette_col,col_names=col_names,row_names=row_names, 
                      dim=grid_dim_list,row_str_length=max_row_annot_str,col_str_length=max_col_annot_str,legend_fontsize=legend_fontsize,
                      family=legend_font_style, legend_fontface=legend_fontface, legend_color=legend_color, ...)
  
  #Draw matrix  
  if(is.null(mat_color)){ 
    pairs.breaks <- seq(min(data),max(data), by=0.05)
    brk <- length(pairs.breaks)-1
    #mat_color <- colorpanel(n=brk,low="green",mid="yellow",high="red")
    mat_color = colorRampPalette(rev(brewer.pal(n = 9, name = "RdYlGn")))(nrow(data)*ncol(data))
  }
  mat_select_color <- draw_mat(data,dim=grid_dim_list,breaks=breaks,mat_color= mat_color,cell_border=cell_border,cell_border_col=cell_border_col,
                               display_number=display_number)
  
  ## Draw matrix Legend
  draw_mat_legend(data,mat_color,breaks=breaks,mat_legend_size=mat_legend_size,dim=grid_dim_list)
  
  #Draw names block
  if (display_colnames){
    draw_names(col_names=col_names,total_names=length(col_names),font_size=col_fontsize,family= names_font_style, names_fontface=names_fontface,
               names_color=names_color,dim=grid_dim_list)
  }
  if (display_rownames){
    draw_names(row_names=row_names,total_names=length(row_names),font_size=row_fontsize,family= names_font_style,names_fontface=names_fontface,
               names_color=names_color,dim=grid_dim_list)
  }
  
}
