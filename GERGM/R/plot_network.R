#' @title Plots of value-edged networks.
#' @description Generates a visualization of a value-edged network.
#'
#' @param sociomatrix A square numeric matrix (socimatrix) with real valued edges
#'(no NA's).
#' @param threshold The threshold for removing edges from the network in order to
#' calculate the positions for the nodes using the futcherman reingold algorithm.
#' The value is multiplied against max(abs(sociomatrix)) to determine the
#' threshold. Defaults to 0.5.
#' @param save_pdf Logical indicating whether the plot should be saved to a PDF.
#' @param pdf_name The name we would like to give to the output file. Be sure to
#' include a ".pdf" extension.
#' @param output_directory The directory where the user would like to output the
#' PDF if save_pdf == TRUE.
#' @examples
#' set.seed(12345)
#' sociomatrix <- matrix(rnorm(400,0,20),20,20)
#' colnames(sociomatrix) <- rownames(sociomatrix) <- letters[1:20]
#' plot_network(sociomatrix)
#' @export
plot_network <- function(sociomatrix,
                         threshold = 0.5,
                         save_pdf = FALSE,
                         pdf_name = "Test.pdf",
                         output_directory = "./"
                         ){

  # check input
  if(class(sociomatrix) != "matrix" & class(sociomatrix) != "data.frame"){
    stop("You must provide the network as a numeric matrix.")
  }

  if(nrow(sociomatrix) != ncol(sociomatrix)){
    stop("You must provide a square matrix.")
  }

  diag(sociomatrix) <- 0

  # create temporary matrices that can be altered
  temp <- temp2 <- matrix(sociomatrix[,],nrow(sociomatrix),ncol(sociomatrix))

  # determine the threshold for removing edges
  cutoff <- max(abs(temp))*threshold

  # remove edges
  temp[which(abs(temp ) < cutoff)] <- 0

  # create a network object using adjacency matrix with edges removed
  net <- igraph::graph.adjacency(temp ,mode="directed",
                                 weighted=TRUE,diag=FALSE)

  #create layout with Fuchterman Reingold
  layout <- igraph::layout_with_fr(net, weights = igraph::E(net)$weight)

  # create a second network object with the un-truncated network
  net2 <- igraph::graph.adjacency(temp2,mode="directed",
                                 weighted=TRUE,diag=FALSE)

  # get an edgelist
  edgelist <- igraph::get.edgelist(net2)
  # get the edge weights
  weights <- igraph::E(net2)$weight

  # order edgeweights from smallest absolute value to largest
  ordering <-order(abs(weights), decreasing = F)
  edgelist <- edgelist[ordering,]
  weights <- weights[ordering]

  # generate edge colors
  negcolors <- colorRampPalette(c('red','black'))
  poscolors <- colorRampPalette(c('black','blue'))
  negcolors <- negcolors(25)
  poscolors <- poscolors(25)

  # generate edge widths
  negbreaks <- seq(min(weights), 0, length.out = 26)
  posbreaks <- seq(0, max(weights), length.out = 26)
  widbreaks <- seq(0,max(abs(weights)),length.out = 50)
  widths <- seq(0,5,length.out = 50)


  ##### If we are saving a PDF
  if(save_pdf){
    #get current working directory
    cur_directory <- getwd()
    setwd(output_directory)

    pdf(file = pdf_name, width = 12, height = 12)
    #start plot
    par(bg = "black", mar = c(2,2,2,2),xpd=TRUE)
    plot(layout,pch = 20, cex = 1, col = "black", axes = F, xlab = "", ylab = "",
         xlim = c((min(layout[,1])-2), (max(layout[,1])+2)),
         ylim = c((min(layout[,2])-2), (max(layout[,2])+2)))

    # add in edges
    for(i in 1:length(weights)){
      cur1 <- layout[edgelist[i,1],]
      cur2 <- layout[edgelist[i,2],]
      curweight <- weights[i]

      # find edge color
      nf <- TRUE
      counter <- 1
      bin <- 1
      while(nf){
        if(curweight > 0){
          if(posbreaks[counter] >= curweight){
            bin <- counter
            nf <- FALSE
          }
        }else{
          if(negbreaks[counter] >= curweight){
            bin <- counter
            nf <- FALSE
          }
        }
        counter <- counter +1
      }

      # find edge width
      nf <- TRUE
      counter <- 1
      wid <- 1
      while(nf){
        if(widbreaks[counter] >= abs(curweight)){
          wid <- counter
          nf <- FALSE
        }
        counter <- counter +1
      }
      if(curweight > 0){
        lines(c(cur1[1],cur2[1]) , c(cur1[2],cur2[2]),
              col = poscolors[bin], lwd = widths[wid])
      }else{
        lines(c(cur1[1],cur2[1]) , c(cur1[2],cur2[2]),
              col = negcolors[bin], lwd = widths[wid])
      }
    }
    text(layout,labels = rownames(sociomatrix), col = "white")
    legend("bottom", inset=0, title = "Edge Values",title.col = "white",
           legend =c(round(min(sociomatrix),2), round(max(sociomatrix),2)),
           fill=c("red","blue"), horiz=T, bg = "black",text.col = "white")
    dev.off()


    #reset working directory
    setwd(cur_directory)
  }else{
    #start plot
    par(bg = "black", mar = c(2,2,2,2),xpd=TRUE)
    plot(layout,pch = 20, cex = 1, col = "black", axes = F, xlab = "", ylab = "",
         xlim = c((min(layout[,1])-2), (max(layout[,1])+2)),
         ylim = c((min(layout[,2])-2), (max(layout[,2])+2)))

    # add in edges
    for(i in 1:length(weights)){
      cur1 <- layout[edgelist[i,1],]
      cur2 <- layout[edgelist[i,2],]
      curweight <- weights[i]

      # find edge color
      nf <- TRUE
      counter <- 1
      bin <- 1
      while(nf){
        if(curweight > 0){
          if(posbreaks[counter] >= curweight){
            bin <- counter
            nf <- FALSE
          }
        }else{
          if(negbreaks[counter] >= curweight){
            bin <- counter
            nf <- FALSE
          }
        }
        counter <- counter +1
      }

      # find edge width
      nf <- TRUE
      counter <- 1
      wid <- 1
      while(nf){
        if(widbreaks[counter] >= abs(curweight)){
          wid <- counter
          nf <- FALSE
        }
        counter <- counter +1
      }
      if(curweight > 0){
        lines(c(cur1[1],cur2[1]) , c(cur1[2],cur2[2]),
              col = poscolors[bin], lwd = widths[wid])
      }else{
        lines(c(cur1[1],cur2[1]) , c(cur1[2],cur2[2]),
              col = negcolors[bin], lwd = widths[wid])
      }
    }
    text(layout,labels = rownames(sociomatrix), col = "white")
    legend("bottom", inset=0, title = "Edge Values",title.col = "white",
           legend =c(round(min(sociomatrix),2), round(max(sociomatrix),2)),
           fill=c("red","blue"), horiz=T, bg = "black",text.col = "white")
  }

  par(bg = "white")
  # do not return anything.
}

