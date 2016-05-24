#' Draw graph function
#' 
#' A function to draw the target inhibition network.
#' 
#' @param draw_data a data frame combining drug-target interaction data with drug sensitivity. The column names
#' must be upper case.
#' @return An image in both pdf and nnf format of the estimated target inhibition network.
#' @author Jing Tang \email{jing.tang@@helsinki.fi}
#' @references Tang J, Karhinen L, Xu T, Szwajda A, Yadav B, Wennerberg K, Aittokallio T. 
#' Target inhibition networks: predicting selective combinations of druggable targets to block cancer 
#' survival pathways. PLOS Computational Biology 2013; 9: e1003226.
#' @examples 
#' \dontrun{
#' data(tyner_interaction_binary)
#' data(tyner_sensitivity)
#' y<-tyner_sensitivity[,1]
#' k_selected<-sffs(tyner_interaction_binary, y)$k_sel
#' x<-data.frame(tyner_interaction_binary[, k_selected])
#' #binarize the sensitivity data
#' one<-which(y>0.5)
#' zero<-which(y<=0.5)
#' SENS<-y
#' SENS[one]<-1
#' SENS[zero]<-0
#' draw_data<-cbind(x, SENS)
#' drawGraph(draw_data)
#' }

drawGraph <- function(draw_data) {
    
    # column names must be upper case and contains only letters!
    column_names_actual <- dimnames(draw_data)[[2]]
    column_names <- unlist(lapply(seq(1:length(column_names_actual)),function(i) paste(LETTERS[i],LETTERS[i],sep="")))
    colnames(draw_data) = column_names
    dimnames(draw_data)[[2]][length(column_names)] <- "SENS"
    
    boolean <- eqmcc(draw_data, outcome = "SENS",row.dom=F,omit=1)
    #boolean$essential
    num_components <- length(boolean$essential)
    expressions <- sapply(boolean$essential, function(x) strsplit(x, "\\*"))
    
    expressions_compact <- list()
    # only treat the true conditions
    is.upper <- "[A-Z]"
    num <- 1
    for (i in 1:num_components) {
        
        result <- grepl(pattern = is.upper, x = expressions[[i]])
        if(any(result)){
          expressions_compact[[num]] <- unlist(lapply(expressions[[i]][which(result == TRUE)],
                                                    function(x) unlist(strsplit(x, ".", fixed=TRUE))[1]))
          num <- num + 1
        }
        
    }
    
    # remove duplicated compact expressions
    expressions_compact = expressions_compact[which(duplicated(expressions_compact)==F)]
    
    # change the labels back to the target names
    expressions_compact = lapply(expressions_compact,function(i) column_names_actual[grep(paste(i,collapse="|"),column_names)])
    
    num_components <- length(expressions_compact)
    height_components <- lapply(expressions_compact, length)
    
    
    seg <- 1 # the width for one component
    seg_in <- 0.2 # the width of the small line within the components
    seg_out <-0.2 # the width between the components
    scale = mean(unlist(lapply(expressions_compact,function(x) max(nchar(x))/3)))
    
    height_figure <- max(unlist(height_components))
    width_figure <- num_components * (seg*scale + seg_out)
    
    # drawing
    dummy <- 0
    pdf(file="targetInhibitionNetwork.pdf", width=width_figure, height=height_figure,pointsize=12)
    par(mar=c(0,0,0,0))
    plot(dummy, dummy, type = "n", axes = FALSE, ann = FALSE, xlim = c(0, width_figure), ylim = c(0, height_figure))  # no axes, no labels
    start_line <- seg_in
    
    for (i in 1:num_components) {
        leg <- height_components[[i]]
        len_target <- c()
        for(j in 1:leg){
          len_target <- c(len_target, nchar(expressions_compact[[i]][j]))
        }
        max_len <- max(len_target) # the length of the longest target names in the component

        y0 <- height_figure/2 - (leg - 1)/2 # y for the lowest target
        y0 <- seq(y0, y0 + leg - 1) # y for the other targets

        x0 <- rep(start_line, leg)
        x1 <- rep(start_line+seg_in, leg)
        y1 <- y0
        x2 <- x1 + (seg - 2*seg_in)*max_len/3
        x3 <- x2 + seg_in
        
        segments(x0, y0, x1, y1)
        segments(x2, y0, x3, y1)
        lines(x0[c(1, 1)], y0[c(1, leg)], type = "l")
        lines(x3[c(1, 1)], y0[c(1, leg)], type = "l")
        lines(c(x3[1], x3[1] + seg_out), c(height_figure/2, height_figure/2), type = "l")
        
        start_line <- x3[1] + seg_out
        
        # write names
        text(x1, y0, labels = expressions_compact[[i]], pos = 4)
    }
    dev.off()
    
    # two terminal version, output to cytoscape nnf file
    cat(c(),file="targetInhibitionNetwork.nnf")
    for (i in 1:num_components){
      write(paste(c('TargetInhibitionNetwork',paste('M',i,sep="")),collapse="\t"),file="targetInhibitionNetwork.nnf",sep="\n",append=T)
    }
    
    Terminals <- paste("T",1:2, sep="") # the number of terminals for the SIF format
    
    for (i in 1:num_components){
      for (j in expressions_compact[[i]]){
      lines1 <-paste(c(paste('M',i,sep=""),Terminals[1],"pp",j), collapse="\t")
      lines2 <-paste(c(paste('M',i,sep=""),j,"pp",Terminals[2]), collapse="\t")
      write(lines1,file="targetInhibitionNetwork.nnf",sep="\n",append=T)
      write(lines2,file="targetInhibitionNetwork.nnf",sep="\n",append=T)
      }
    } 
    cat()
} 
