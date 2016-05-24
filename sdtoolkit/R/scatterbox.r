
#A function to combine scatter plots of color-coded points with boxes
#boxseq: a box sequence object, such as that output by sdprim,
#or, if "manboxlist"=TRUE, a list object of appropriate form for making your
#own boxes
#boxnum: Which box in the box sequence would you like to plot?
#data: object of data
#y=output - will default to last column of data
#xdim, ydim=which dimensions to make x and y in plot? Safer to use dimnmame
#filterdims=Should the scatter plot remove points that are filtered by 
#the dimension not shown?
#filterothers= Should the scatter plot remove points that are captured by 
#other boxes, and thus only show the effect of this box?  Can be a logical (FALSE for no removal), or TRUE for remove all other boxes, "PREV" or a vector
#PREV removes those boxes with lower box indices
#vector specifies which points in box seq to remove.
#nobox=just plot the points, don't plot any box
#manboxlist: A list, each component of which represents a box as a nother list, with the first component being a vector indicating which dimensions are restricted and the second component d by 2 matrix giving the d dimensional restrictions. 
#addbox means only add a boxplot to an existing plot
#... gets passed on to colptplot and box plot       

        
scatterbox <- function(boxseq, boxnum=1, xdim, ydim, filterdims=FALSE, filterothers=FALSE, nobox=FALSE,addbox=FALSE, data=attr(boxseq,"data"), y=NULL, manualbox=FALSE,...){

  if(is.null(y)){  #assumes is of form output by sdprim
    ycol <- which(colnames(data)=="output")
    filtered <- data[,1:ycol]
  }else{
    filtered <- cbind(data,y)
  }
  
#Check if the dimensions to be plotted are given as numeric, or as dimnames,
#and make them numeric as appropriate

  if(is.character(xdim)){
    nxdim = which(colnames(filtered)==xdim)
  }else{nxdim=xdim}
  
  if(is.character(ydim)){
    nydim = which(colnames(filtered)==ydim)
  }else{nydim=ydim}

  xlims <- range(filtered[,nxdim])
  ylims <- range(filtered[,nydim])

#FIND which boxes to filter out, if any:

if(filterothers==FALSE){
  filterbs <- c()
} else if (filterothers==TRUE){
  filterbs <- c(1:length(boxseq))[-boxnum]
} else if (filterothers=="PREV"){
  filterbs <- c(1:(boxnum-1))
} else if (is.numeric(filterothers)){#do nothing
  filterbs <- filterothers
}

#NEXT: Make sure box sequence is in appropriate format

  lboxseq <- list()

  for(i in 1:length(boxseq)){
    
    if(manualbox){
    
      A <- is.vector(boxseq[[i]][[1]])  
      B <- is.matrix(box[[i]][[2]])
      C <- dim(boxseq[[i]][[2]])==c(length(box[[i]][[1]]),2)
  
      if(!(A && B && C)){stop("Error - box sequence is not in the appropriate format")}
  
      lboxseq <- boxseq
  
    } else{
      
      A <- is.vector(boxseq[[i]]$lbox[[1]])  
      B <- is.matrix(boxseq[[i]]$lbox[[2]])
      C <- dim(boxseq[[i]]$lbox[[2]])==c(length(boxseq[[i]]$lbox[[1]]),2)
  
      lboxseq[[i]] <- boxseq[[i]]$lbox
  
      if(!(A && B && C)){stop("Error - box sequence is not in the appropriate format")}
    }
  }
 
  #Filter out pts as appropriate:

  for (i in filterbs){   #filter out other boxes
      filtered <- filterpts(filtered, box=lboxseq[[i]], keepin=FALSE)
  }
  
  
  if(filterdims  && (length(lboxseq[[boxnum]][[1]]>2))){   #filter out other dimensions
  
    #make a subbox lacking the two dims
    ttbox <- lboxseq[[boxnum]]
    tbox  <- list()
    tbox[[1]] <- ttbox[[1]][ttbox[[1]]!=nxdim & ttbox[[1]]!=nydim]
    tbox[[2]] <- matrix(ttbox[[2]][ttbox[[1]]!=nxdim & ttbox[[1]]!=nydim,], ncol=2)
    
    filtered <- filterpts(filtered, box=tbox, keepin=TRUE)
  }

  if(!addbox){

  #actually plot the points:
  colptplot(filtered,nxdim,nydim,outdim=ncol(filtered),lowcol="transparent",
     hicol="black",xname=colnames(filtered)[nxdim],yname=colnames(filtered)[nydim],xlims=xlims,ylims=ylims,...)
  
  }
  
  #actually plot the box[es]

  pbox(sdobj=lboxseq[[boxnum]],xdim=nxdim, ydim=nydim, boxnum=NA, fromtype="oldbox",lwd=3,gborder="blue",mdborder="red",col=NA)
  
}
 
#
#        
#        
#          #CONVERT BOX to 'old' form so box can plot it
#          boxy <- list(box=trajinf$box[[boxind]],dimlist=trajinf$dimlist[[boxind]])
#          ofbox <- boxconverter(boxy)
#          
#        
#          #EXTRACT TWO MOST IMPORTANT DIMS
#          d1name <- as.character(curmats[[2]][1,1])
#          d1 <- curmats[[2]][1,4]
#          
#          if(curmats[[2]][1,1]==curmats[[2]][2,1]){
#            d2name <- as.character(curmats[[2]][3,1])
#            d2 <- curmats[[2]][3,4]
#          }else{
#            d2name <- as.character(curmats[[2]][2,1])
#            d2 <- curmats[[2]][2,4]
#          }
#          
#   
#       
#        
#          cat("Enter the row numbers for two dimensions (in the left most column, 
#  not actual name) you would like to use for the x and y axis separated","\n",
#          "by a comma, or enter 'd' to accept the default, which is the two highest","\n","ranked variables.")
#          xandy <- readline(cat("\n","  (In this case",d1name,"and",d2name,".)","\n"))
#          
#          if(xandy!="d"){
#            
#            #EXTRACT x and y
#            
#            ds <- strsplit(xandy,",")[[1]]
#            
#            d1t <- as.numeric(ds[1])
#            d2t <- as.numeric(ds[2])
#            
#            d1 <- curmats[[2]][d1t,4]
#            d2 <- curmats[[2]][d2t,4]
#            
#            d1name <- as.character(curmats[[2]][d1t,1])
#            d2name <- as.character(curmats[[2]][d2t,1])
#            
#            #OLD and I think bad way:
#            #d1name <- curmats[[2]][,1][which(curmats[[2]][,1]==d1)]
#            #d2name <- curmats[[2]][,1][which(curmats[[2]][,1]==d2)]
#            
#          }
#          
#          #OPEN NEW DEVICE SO TRAJ STAYS OPEN  - but keep first box
#          #IF, stuff...
#      
#          #PLOT POINTS
#      
#          options("device")$device()
#          colptplot(cbind(x,y),xdim=d1,ydim=d2,outdim=(ncol(x)+1),lowcol="transparent",hicol="black",xname=d1name,yname=d2name)
#      
#      
#          #PLOT BOX
#          pbox(sdobj=ofbox,xdim=d1,ydim=d2,boxnum=NA,fromtype="oldbox",lwd=3,gborder="blue",mdborder="red",col=NA)
#
#        
#          bringToTop(-1)
#          #PROMPT FOR PLOT OF DIFFERENT VARIABLES - loop back to TOP
#          redo <- readline(cat("Would you like to plot different variables? ('y' or 'n')","\n"))
#          if (redo == "y"){
#            psatisfied <- FALSE
#          } else{
#            psatisfied <- TRUE
#          }
#          
#          #IF NOT FIRST COVERING - offer to plot an earlier box in the sequence?
#         
#         
#         