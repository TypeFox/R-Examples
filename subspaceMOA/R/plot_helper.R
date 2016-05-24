make_plot_matrix <- function(list_of_plots,ncol) {
  if(is.null(list_of_plots)) {
    return()
  }
  gridExtra::arrangeGrob(grobs=list_of_plots,ncol=ncol,padding=grid::unit(0,"cm"),widths=rep(100,ncol),heights=rep(100,ncol))
}

#Helper function that reads the micro and macro cluster positions and weights
#from the dsc object and then puts it into a data frame
#'@import stats 
format_data_from_dsc <- function(dsc,points=NULL) {
  
  microclusters <- get_centers(dsc,type="micro")
  if(is.null(microclusters)) return(NULL)
  microclusters$weight <- get_weights(dsc,type="micro")
  microclusters$type <- "Microcluster"
  
  macroclusters <- get_centers(dsc,type="macro")
  if(!is.null(macroclusters)) {
    macroclusters$weight <- get_weights(dsc,type="macro")
    macroclusters$type <- "Macrocluster"
    dataframe <- rbind(microclusters,macroclusters)
  } else if(!is.null(microclusters)) {
    dataframe <- microclusters
  }
  #Add the points if there are any
  if(!is.null(points)) {
    points$weight <- rep(1,nrow(points))
    if(!is.null(points$class)) {
      points$type <- sapply(points$class,function(class_num){return(paste0("Data Point in Class ",class_num))})
      points$class <- NULL
    } else {
      points$type <- "Data Point"
    }
    dataframe <- rbind(points,dataframe)
  }
  #Turn the type attribute into a factor so that ggplot can handle it better
  dataframe$type <- factor(dataframe$type)
  #Set the levels correctly so that microclusters and macroclusters come first
  if(!is.null(macroclusters))dataframe$type <- relevel(dataframe$type,"Macrocluster")
  dataframe$type <- relevel(dataframe$type,"Microcluster")
  return(dataframe)
}
create_plot_matrix <- function(dataframe) {
  if(is.null(dataframe)) {
    return(NULL)
  }
  #Make a list and fill it with the plots that are going to fill the plot matrix(row-wise)
  number_of_dimensions <- ncol(dataframe) - 2
  list_of_plots <- list()
  for(i in 1:number_of_dimensions) {
    j <- 1
    for(j in 1:number_of_dimensions) {
      axis_label_below <- F
      axis_label_left <- F
      if(i==number_of_dimensions){
        axis_label_below <- T
      }
      if(j==1) {
        axis_label_left <- T
      }
      list_of_plots[[(i-1)*number_of_dimensions + j]] <- 
        dataframe %>%
        basic_plot_from_dataframe(dims=c(j,i)) %>%
        style_plot_for_plotmatrix(axis_label_left=axis_label_left,
                                  axis_label_below=axis_label_below)
    }
  }
  return(list_of_plots)
}

basic_plot_from_dataframe <- function(dataframe,dims) {
  #Create a basic ggplot plot based on a data frame as
  #produced by format_data_from_dsc.
  first_dimension_string <- names(dataframe)[dims[1]]
  second_dimension_string <- names(dataframe)[dims[2]]
  p <- ggplot(dataframe,aes_string(x=first_dimension_string,y=second_dimension_string)) + 
    geom_point(aes_string(shape="type",color="type",size="weight"))
  
  #Now we define plot symbols for Microclusters, Macroclusters and points.
  #This is done implicitly based on the ordering of the types of data
  #in the data frame. Shapes used are 1 (circles), 3 (crosses). 
  #Some data might not contain macroclusters so in these cases, no color
  #and shape should be defined for them.
  macroclusters_present <- any(dataframe$type=="Macrocluster")
  if(macroclusters_present) {
    p <- p +
    scale_color_manual(values=c("red","blue",rep("grey",20))) +
    scale_shape_manual(values=c(1,3,4:20))
  } else {
    p <- p +
    scale_color_manual(values=c("red",rep("grey",20))) +
    scale_shape_manual(values=c(1,4:20))
  }
  
  return(p)
}
#Styles an existing ggplot in such a way that it can be drawn in a scatterplot matrix
style_plot_for_plotmatrix <- function(plot,axis_label_left=F,axis_label_below=F) {
  if(is.null(plot)) {
    return(NULL)
  }
  plot <- plot + 
    #Set the minimum and maximum sizes of the micro and macrocluster marks
    scale_size(range=c(1,10)) +
    #Show no legend and set the margins to a size of 1mm
    theme(legend.position="none",
          #Make margins in order top right bottom left
          plot.margin=grid::unit(c(0,0,0,0),"npc"),
          axis.text.y=element_blank(),
          axis.ticks.y=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x=element_blank(),
          panel.border=element_rect(colour="black",fill=NA,size=1))
  
  if(!axis_label_left) {
    plot <- plot + labs(y=NULL)
  }
  if(!axis_label_below) {
    plot <- plot + labs(x=NULL)
  }
  return(plot)
}
#Styles a plot in such a way that it can be shown on its own
style_plot_for_detail <- function(plot) {
  plot <- plot +
    #Set the minimum and maximum sizes of the micro and macrocluster marks
    scale_size(range=c(3,30)) + 
    #Do not show a legend for sizes
    guides(size=F) 
  
}
r_logical_to_js_boolean_string <- function(x) {
  if(x) return("true")
  else return("false")
}
