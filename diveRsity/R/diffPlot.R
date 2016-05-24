################################################################################
# diffPlot, plot all pairwise population pairs                                  #
################################################################################
#' @export
diffPlot <- function (x, outfile= NULL, interactive = FALSE) {
  #x <- diveRsity::fastDivPart("test.gen", pairwise = T, fst = T)
  #y <- diveRsity::diffCalc("test.gen", pairwise = T, fst = T)
  on <- outfile
  inta <- interactive
  #output from divPart
  if(!is.element("pairwise", names(x))){
    stop(cat(paste("To plot pairwise differentiation the argument ",
                   "'pairwise' in\neither 'diffCalc' or fastDivPart'",
                   "must be TRUE!", sep = "")))
  }
  if(is.element("gstEst", names(x$pairwise))){
    idx <- c(1, 4, 2, 3)
    fun <- "fastDivPart" 
  } else {
    idx <- c(1, 5, 2, 4)
    fun <- "fastDivPart"
  }
  
  if(is.null(on) == TRUE && inta == TRUE){
    of = paste(getwd(),"/", sep = "")
  } else {
    suppressWarnings(dir.create(paste(getwd(), "/", 
                                      on, "-[diffPlot]", "/", sep="")))
    of=paste(getwd(),"/",on,"-[diffPlot]","/",sep="")
  }
  
  if(!exists("inta", -1)){
    inta <- FALSE
  }
  if(inta == TRUE) {
    sp.header<-list()
    colleer<-list()
    colleer <<- colorRampPalette(c("blue","white"))
    requireNamespace("sendplot")
    direct<-of
    pwc<-combn(ncol(x[[3]][[1]]),2)
    pwNames<-paste(colnames(x[[3]][[1]])[pwc[1,]],
                   colnames(x[[3]][[1]])[pwc[2,]],
                   sep=' vs ')
    
    gst_lab <- round(as.vector(x[[3]][[idx[1]]]), 4)
    gst_lab <- na.omit(gst_lab)
    collab111<-list()
    #
    if(length(x[[3]]) > 3 && fun == "fastDivPart"){
      fst_lab <- round(as.vector(x[[3]][[idx[2]]]), 4)
      fst_lab<-na.omit(fst_lab)
    }
    #
    gpst_lab <- round(as.vector(x[[3]][[idx[3]]]), 4)
    gpst_lab<-na.omit(gpst_lab)
    #
    Dest_lab <- round(as.vector(x[[3]][[idx[4]]]), 4)
    Dest_lab<-na.omit(Dest_lab)
    #
    
    fl_ext<-c(".tif","Dot.png","Dot.tif")
    if(length(x[[3]]) > 3 && fun == "fastDivPart"){
      xy.labels <-  data.frame(pops = pwNames,
                               Nei_Gst = gst_lab,
                               Weir_Theta = fst_lab,
                               Hedrick_Gst = gpst_lab,
                               Jost_D = Dest_lab)
    } else {
      xy.labels <-  data.frame(pop = pwNames,
                               Nei_Gst = gst_lab,
                               Hedrick_Gst = gpst_lab,
                               Jost_D = Dest_lab)
    }
    #Nei Gst
    abx<-list()
    abx<<-x
    collab111 <<- c(round(min(gst_lab),3),
                    round(mean(gst_lab),3),
                    round(max(gst_lab),3))
    
    plot.call <- paste("image(1:nrow(abx[[3]][[", 
                       idx[1],
                       "]]), 1:nrow(abx[[3]][[", 
                       idx[1], 
                       "]]), abx[[3]][[",
                       idx[1], 
                       "]], ylab='', xlab='', main='Pairwise Gst', ",
                       "xaxt='n', yaxt='n', col = colleer(50), ",
                       "las=1, cex.main=3)", sep = "")
    ##
    plot.extras <- paste("plotrix::color.legend(nrow(abx[[3]][[",
                         idx[1], 
                         "]])/5, nrow(abx[[3]][[",
                         idx[1],
                         "]])/3, nrow(abx[[3]][[",
                         idx[1],
                         "]])/4, nrow(abx[[3]][[",
                         idx[1],
                         "]])/1.2, collab111, rect.col=colleer(50), ",
                         "gradient='y', cex=3)", sep = "")
    ##
    suppressWarnings(sendplot::imagesend(plot.call=plot.call,
                                         x.pos=pwc[2,],
                                         y.pos=pwc[1,],
                                         xy.type="points",
                                         xy.labels = xy.labels,
                                         plot.extras=plot.extras,
                                         fname.root="Gst_matrix",
                                         dir=of,
                                         image.size="1050x800",
                                         font.size=18,
                                         spot.radius = 10,
                                         font.type = "Arial",
                                         window.size="1100x800"))
    #clean up
    unlink(paste(of,"Gst_matrix",fl_ext,sep=""))
    #
    #
    #
    #
    #
    #Fst
    if(length(x[[3]]) > 3 && fun == "fastDivPart"){
      collab111 <<- c(round(min(fst_lab),3),
                      round(mean(fst_lab),3),
                      round(max(fst_lab),3))
      plot.call <- paste("image(1:nrow(abx[[3]][[",
                         idx[2],
                         "]]),1:nrow(abx[[3]][[",
                         idx[2],
                         "]]), abx[[3]][[",
                         idx[2],
                         "]],ylab = '',xlab = '',xaxt = 'n',yaxt = 'n', ",
                         "main = 'Pairwise Fst', col = colleer(50), ",
                         "las = 1,cex.main = 3)", sep = "")
      ##
      plot.extras <- paste("plotrix::color.legend(nrow(abx[[3]][[",
                           idx[2],
                           "]])/5, nrow(abx[[3]][[",
                           idx[2],
                           "]])/3, nrow(abx[[3]][[",
                           idx[2],
                           "]])/4, nrow(abx[[3]][[",
                           idx[2],
                           "]])/1.2, collab111, rect.col=colleer(50), ",
                           "gradient='y', cex=3)", sep = "")
      #
      suppressWarnings(sendplot::imagesend(plot.call=plot.call,
                                           x.pos=pwc[2,],
                                           y.pos=pwc[1,],
                                           xy.type="points",
                                           xy.labels = xy.labels,
                                           plot.extras=plot.extras,
                                           fname.root="Fst_matrix",
                                           dir=of,
                                           image.size="1050x800",
                                           font.size=18,
                                           spot.radius = 10,
                                           font.type ="Arial",
                                           window.size="1100x800"))
      #clean up
      unlink(paste(of,"Fst_matrix",fl_ext,sep=""))
    }
    #
    #
    #
    #
    #
    #G'st
    collab111 <<- c(round(min(gpst_lab),3),
                    round(mean(gpst_lab),3),
                    round(max(gpst_lab),3))
    plot.call <- paste("image(1:nrow(abx[[3]][[",
                       idx[3],
                       "]]),1:nrow(abx[[3]][[",
                       idx[3],
                       "]]), abx[[3]][[",
                       idx[3],
                       "]], ylab='', xlab='', xaxt='n', yaxt='n', ",
                       "main='Pairwise Gst (Hedrick)',col = colleer(50), ",
                       "las=1,cex.main=3)", sep = "")
    
    plot.extras <- paste("plotrix::color.legend(nrow(abx[[3]][[",
                         idx[3],
                         "]])/5,nrow(abx[[3]][[",
                         idx[3],
                         "]])/3, nrow(abx[[3]][[",
                         idx[3],
                         "]])/4,nrow(abx[[3]][[",
                         idx[3],
                         "]])/1.2,collab111, rect.col=colleer(50), ",
                         "gradient='y',cex=3)", sep = "")
    ##
    suppressWarnings(sendplot::imagesend(plot.call=plot.call,
                                         x.pos=pwc[2,],
                                         y.pos=pwc[1,],
                                         xy.type="points",
                                         xy.labels = xy.labels,
                                         plot.extras=plot.extras,
                                         fname.root="G_prime_st_matrix",
                                         dir=of,
                                         image.size="1050x800",
                                         font.size=18,
                                         spot.radius = 10,
                                         font.type = "Arial",
                                         window.size="1100x800"))
    #clean up
    unlink(paste(of,"G_prime_st_matrix",fl_ext,sep=""))
    #
    #
    #
    #
    #
    #
    #Dest
    collab111 <<- c(round(min(Dest_lab),3),
                    round(mean(Dest_lab),3),
                    round(max(Dest_lab),3))
    plot.call <- paste("image(1:nrow(abx[[3]][[",
                       idx[4],
                       "]]),1:nrow(abx[[3]][[",
                       idx[4],
                       "]]), abx[[3]][[", 
                       idx[4],
                       "]], ylab = '',xlab = '',xaxt = 'n', ",
                       "yaxt = 'n',main = 'Pairwise D (Jost)', ",
                       "col = colleer(50),las=1,cex.main=3)", sep = "")
    
    plot.extras <- paste("plotrix::color.legend(nrow(abx[[3]][[",
                         idx[4],
                         "]])/5,nrow(abx[[3]][[",
                         idx[4],
                         "]])/3, nrow(abx[[3]][[",
                         idx[4],
                         "]])/4,nrow(abx[[3]][[",
                         idx[4],
                         "]])/1.2,collab111, rect.col = colleer(50), ",
                         "gradient='y',cex=3)", sep = "")
    ##
    suppressWarnings(sendplot::imagesend(plot.call=plot.call,
                                         x.pos=pwc[2,],
                                         y.pos=pwc[1,],
                                         xy.type="points",
                                         xy.labels = xy.labels,
                                         plot.extras=plot.extras,
                                         fname.root="D_matrix_",
                                         dir=of,
                                         image.size="1050x800",
                                         font.size=18,
                                         spot.radius = 10,
                                         font.type = "Arial",
                                         window.size="1100x800"))
    #lean up
    
    unlink(paste(of,"D_matrix_",fl_ext,sep=""))
    
  } else {
    
    #
    if(length(x[[3]]) > 3 && fun == "fastDivPart"){
      par(mfrow=c(2,2))
    } else {
      par(mfrow=c(3,1))
    }
    colleer<-colorRampPalette(c("blue","white"))
    cols<-colleer(50)
    #Gst
    image(1:nrow(x[[3]][[idx[1]]]),
          1:nrow(x[[3]][[idx[1]]]),
          x[[3]][[idx[1]]],
          ylab="population",
          xlab="population",
          main="Pairwise Gst",
          col=cols,
          las=1)
    gst<-as.vector(x[[3]][[idx[1]]])
    gst<-as.vector(na.omit(gst))
    collab111<-c(round(min(gst),3),
                 round(mean(gst),3),
                 round(max(gst),3))
    
    plotrix::color.legend(nrow(x[[3]][[idx[1]]])/5,
                          nrow(x[[3]][[idx[1]]])/3,
                          nrow(x[[3]][[idx[1]]])/4,
                          nrow(x[[3]][[idx[1]]])/1.2,
                          collab111,
                          cols,
                          gradient="y")
    if(length(x[[3]]) > 3 && fun == "fastDivPart"){
      #Fst
      image(1:nrow(x[[3]][[idx[2]]]),
            1:nrow(x[[3]][[idx[2]]]),
            x[[3]][[idx[2]]],
            ylab="population",
            xlab="population",
            main="Pairwise Theta",
            col = cols,
            las=1)
      fst<-as.vector(x[[3]][[idx[2]]])
      fst<-as.vector(na.omit(fst))
      collab111<-c(round(min(fst),3),round(mean(fst),3),round(max(fst),3))
      
      plotrix::color.legend(nrow(x[[3]][[idx[2]]])/5,
                            nrow(x[[3]][[idx[2]]])/3,
                            nrow(x[[3]][[idx[2]]])/4,
                            nrow(x[[3]][[idx[2]]])/1.2,
                            collab111,
                            cols,
                            gradient="y")
    }
    #Hedrick's Gst
    image(1:nrow(x[[3]][[idx[3]]]),
          1:nrow(x[[3]][[idx[3]]]),
          x[[3]][[idx[3]]],
          ylab="population",
          xlab="population",
          main="Pairwise G'st",
          col = cols)
    gprimest<-as.vector(x[[3]][[idx[3]]])
    gprimest<-as.vector(na.omit(gprimest))
    collab111<-c(round(min(gprimest),3),
                 round(mean(gprimest),3),
                 round(max(gprimest),3))
    
    plotrix::color.legend(nrow(x[[3]][[idx[3]]])/5,
                          nrow(x[[3]][[idx[3]]])/3,
                          nrow(x[[3]][[idx[3]]])/4,
                          nrow(x[[3]][[idx[3]]])/1.2,
                          collab111,
                          cols,
                          gradient="y")
    #Jost's D
    image(1:nrow(x[[3]][[idx[4]]]),
          1:nrow(x[[3]][[idx[4]]]),
          x[[3]][[idx[4]]],
          ylab="population",
          xlab="population",
          main="Pairwise Jost's D",
          col = cols,
          las=1)
    D<-as.vector(x[[3]][[idx[4]]])
    D<-as.vector(na.omit(D))
    collab111<-c(round(min(D),3),
                 round(mean(D),3),
                 round(max(D),3))
    
    plotrix::color.legend(nrow(x[[3]][[idx[4]]])/5,
                          nrow(x[[3]][[idx[4]]])/3,
                          nrow(x[[3]][[idx[4]]])/4,
                          nrow(x[[3]][[idx[4]]])/1.2,
                          collab111,
                          cols,
                          gradient="y")
  }
  if(exists("abx", where=".GlobalEnv")==TRUE){
    rm(abx, pos=".GlobalEnv")
  }
  if(exists("collab111", where=".GlobalEnv")==TRUE){
    rm(collab111, pos=".GlobalEnv")
  }
  if(exists("colleer", where=".GlobalEnv")==TRUE){
    rm(colleer, pos=".GlobalEnv")
  }
  if(exists("sp.header", where=".GlobalEnv")==TRUE){
    rm(sp.header, pos=".GlobalEnv")
  }
  par(mfrow = c(1,1))
}