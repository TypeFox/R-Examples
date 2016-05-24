plotTimeCourseII<-function(x,plotgroup="",filename="timeseries_multiplot.pdf",numpage=4,cols=2,xname="time",yname="signal",legpos="top",legrow=2,legtitle="treatment",legtitlepos="top",legtextsize=10,legtextcolor="black",legtitlesize=10,legtitlecolor="black",legtitleface="bold",legitemsize=1,plottitlesize=12,plottitleface="bold",xaxissize=10,yaxissize=10,xaxisface="bold",yaxisface="bold",xaxistextsize=8,xaxistextangle=0,yaxistextsize=8,linecolor="Set1")
{

dataset<-x

dataset$time <- as.numeric(as.character(dataset$time))

# a list of unique antibodies from the dataset
ablist<-unique(dataset$ab)
celllinelist<-unique(dataset$cell_line)

# defines plotlist
plotlist<-c("")
plotlist<-as.list(plotlist)


# define standard color scheme if linecolor is not defined
if(linecolor[1]=="Set1"){linecolor<-c("#E41A1C","#377EB8","#4DAF4A","#984EA3","#FF7F00","#FFFF33","#A65628","#F781BF","#999999")}
if(linecolor[1]=="Dark2"){linecolor<-c("#1B9E77","#D95F02","#7570B3","#E7298A","#66A61E","#E6AB02","#A6761D","#666666")}
if(linecolor[1]=="Paired"){linecolor<-c("#A6CEE3","#1F78B4","#B2DF8A","#33A02C","#FB9A99","#E31A1C","#FDBF6F","#FF7F00","#CAB2D6","#6A3D9A","#FFFF99","#B15928")}

# plots single plot for each antibody in plotlist
inum=1
signal<-NULL
sigma<-NULL
treatment<-NULL
for(icell in celllinelist){
datacellline<-dataset[dataset$cell_line == icell,]

for(i in ablist){
# creates a temporary dataset for a single antibody 
tempdata<-datacellline[datacellline$ab == i,]
# renames a handed name of the cell treatment (e.g. perturbation) and renames it "treatment" since aes does not accept variables as names
colnames(tempdata)[colnames(tempdata)==plotgroup]<-c("treatment")


P <- ggplot(tempdata, aes(x=time, y=signal, ymin=signal-sigma, ymax=signal+sigma, colour=treatment, group=paste(treatment,"connection"))) +
      # General Theme: backgroundcolor ...
      theme_bw() + 
      # Sets the color scheme                                                                                                                                           
      scale_color_manual(values=linecolor) +  
      # Set lines between datapoints                                          
      geom_line() +
      # Sets errorbars to each datapoint
      geom_point() + geom_errorbar(width=0) +
      # Creates Title composed of name of the cell line and antibody
      ggtitle(paste(tempdata$cell_line,i,sep=" "))  +
      # Define size and style of the plot title                       
      theme(plot.title = element_text(size=plottitlesize, face=plottitleface)) +
      # Define size and style of the X-axis title
      theme(axis.title.x = element_text(face=xaxisface, size=xaxissize)) +
      # Define size and style of the Y-axis title
      theme(axis.title.y = element_text(face=yaxisface, size=yaxissize)) +
      # Define size and angle of the X-axis text
      theme(axis.text.x = element_text(angle=xaxistextangle, size=xaxistextsize)) +
      # Define font of the Y-axis text 
      theme(axis.text.y = element_text(size=yaxistextsize)) +
      # Define legend position
      theme(legend.position=legpos) + 
      # Define size and color of the legend text
      theme(legend.text = element_text(colour=legtextcolor, size=legtextsize)) + 
      # Define font color, size and style of the legend title                                    
      theme(legend.title = element_text(colour=legtitlecolor, size=legtitlesize, face=legtitleface)) +  
      # Define title of x and y axis
      xlab(xname) + ylab(yname) +                                                    
      #scale_fill_continuous(guide = "legend") +
      # Define item size, number of item rows, title, and title position of the legend
      guides(colour = guide_legend(keywidth = legitemsize, keyheight = legitemsize, nrow=legrow ,title=paste(legtitle),title.position=legtitlepos)) +   
      # Define Y-axis minimum as 0
      expand_limits(y=0) +
      # Define Y-axis maximum as signal + sigma + 0.5
      expand_limits(y=max(tempdata$signal)+max(tempdata$sigma)+0.5)
plotlist[[inum]]<- P
inum<-inum+1
}
}

# Define function whole number
is.wholenumber <- function(x, tol = .Machine$double.eps^0.5)  abs(x - round(x)) < tol 

# Insert dummy if number of plots divided by the number of plots per page is not a whole number
i=length(plotlist)
while(is.wholenumber((length(plotlist))/numpage) != TRUE) {
  i<-i+1
  plotlist[[i]] <-0
  }  
i=0

# Open new pdf document
pdf(file=filename)

# Modified multiplot function for automatic plotting
# Loop which prints each page of the pdf document

ipage=1
while(ipage<=length(plotlist)) {
    # Make a list from the ... arguments and plotlist
    plots<-c("")
    plots<-as.list(plots)
    plots <- c(plotlist[ipage:((ipage+numpage)-1)])

    numPlots = length(plots)

    # Make the panel
    # ncol: Number of columns of plots
    # nrow: Number of rows needed, calculated from # of cols
    layout <- matrix(seq(1, cols * ceiling(numPlots/cols)),
      ncol = cols, nrow = ceiling(numPlots/cols))


      if (numPlots==1) {
        print(plots[[1]])

      } else {
        # Set up the page
        grid.newpage()
        pushViewport(viewport(layout = grid.layout(nrow(layout), ncol(layout))))

        # Make each plot, in the correct location
        for (i in 1:numPlots) {
        # Get the i,j matrix positions of the regions that contain this subplot
        matchidx <- as.data.frame(which(layout == i, arr.ind = TRUE))

        print(plots[[i]], vp = viewport(layout.pos.row = matchidx$row,
                                      layout.pos.col = matchidx$col))
    }
  }

ipage<-ipage+numpage
}

dev.off()
}

