#
# function to make plots off of Splot object
#



makeSplot <- function(Splot,
                      fname.root="Splot",
                      dir="./",
                      overwriteSourcePlot = NA,
                      makeInteractive=TRUE,
                      overrideInteractive=NA,
                      Default=TRUE,
                      header="v3",
                      window.size = "800x1100", # in px
                      returnObj = TRUE
                     
                      ){

  # set up file names
  fname.ps=paste(dir,fname.root,".ps",sep="")
  fname.png=paste(dir,fname.root,".png",sep="")
  fname.jpeg=paste(dir,fname.root,".jpeg",sep="")
  fname.tiff=paste(dir,fname.root,".tif",sep="")
  fname.html=paste(dir,fname.root,".html",sep="")
  fname.pdf=paste(dir,fname.root,".pdf",sep="")
  
  #
  # create static image
  #




  # if single entry check for acceptable type
  if(length(overwriteSourcePlot)==1){

    if(!is.na(overwriteSourcePlot) & ((overwriteSourcePlot != "ps") & (overwriteSourcePlot != "png") & (overwriteSourcePlot != "jpeg") & (overwriteSourcePlot != "tiff") & (overwriteSourcePlot != "pdf"))){
      overwriteSourcePlot = NA
    }
    if(!is.na(overwriteSourcePlot) & ((overwriteSourcePlot != "png") & (overwriteSourcePlot != "jpeg"))){
      overwriteSourcePlot = c("png", overwriteSourcePlot)
    }                                     
    if(is.na(overwriteSourcePlot[1]))  overwriteSourcePlot = Splot$source.plot

    pngF = length(which(overwriteSourcePlot == "png")) > 0
    psF = length(which(overwriteSourcePlot == "ps")) > 0
    tiffF = length(which(overwriteSourcePlot == "tiff")) > 0
    jpegF = length(which(overwriteSourcePlot == "jpeg")) > 0
    pdfF = length(which(overwriteSourcePlot == "pdf")) > 0
       
 
  }else{
    # if vector check acceptable types 
    pngF = length(which(overwriteSourcePlot == "png")) > 0
    psF = length(which(overwriteSourcePlot == "ps")) > 0
    tiffF = length(which(overwriteSourcePlot == "tiff")) > 0
    jpegF = length(which(overwriteSourcePlot == "jpeg")) > 0
    pdfF = length(which(overwriteSourcePlot == "pdf")) > 0
    
    if(!pngF & !jpegF) pngF = TRUE
     
  }

  
  wi = strsplit(Splot$image.size, "x")[[1]][1]
  hi = strsplit(Splot$image.size, "x")[[1]][2]

  startFig = "png"
  if(!pngF) startFig = "jpeg"

  
  if(startFig == "png"){
    png(filename=fname.png, width=as.double(wi), height=as.double(hi), pointsize=Splot$pointsize, res=Splot$res)
    dev.control("enable")
    
    Splot$plot.output = list()
  
    # initiate layout
    if(max(as.vector(Splot$mat),na.rm=TRUE)>1) nf = layout(Splot$mat, respect=TRUE)
    if(Splot$mai.prc) mai.def=par("mai")

    # initialize getLims 
    if(is.null(Splot$plot.lims)){
      plot.lims = list()
      plot.lims$xmins = rep(NA, length(Splot$plot.calls))
      plot.lims$xmaxs = rep(NA, length(Splot$plot.calls))
      plot.lims$ymins = rep(NA, length(Splot$plot.calls))
      plot.lims$ymaxs = rep(NA, length(Splot$plot.calls))
      Splot$plot.lims = plot.lims    
    }
  
    # loop over plot calls to place plots in order or 1:n in layout
    for(i in 1:length(Splot$plot.calls)){
      # set up plot margin
      if(length(Splot$mai.mat)>1){
        if(!Splot$mai.prc) par(mai=Splot$mai.mat[i,])
        if(Splot$mai.prc)  par(mai=Splot$mai.mat[i,]*mai.def)            
      }
      # evaluate plot call
      plt = eval.js(Splot$plot.calls[[i]])

      if(is.null(plt)) Splot$plot.output[i] = NA
      if(!is.null(plt)){
        Splot$plot.output[i] = NA
        class(Splot$plot.output[i]) = "list"
        Splot$plot.output[[i]] = plt
      }
    
      # cycle through all additional plot calls for current plot
      sub.np = length(Splot$plot.extras[[i]])
      for(sp in 1:sub.np){
        if(!is.na(Splot$plot.extras[[i]][[sp]]))
          eval.js(Splot$plot.extras[[i]][[sp]])
      }

      Splot$plot.lims$xmins[i] = par()$usr[1]
      Splot$plot.lims$xmaxs[i] = par()$usr[2]
      Splot$plot.lims$ymins[i] = par()$usr[3]
      Splot$plot.lims$ymaxs[i] = par()$usr[4]
           
    }# end for loop over plot calls


    if(jpegF){
      dev.copy(jpeg, filename=fname.jpeg, width=as.double(wi), height=as.double(hi), pointsize=Splot$pointsize, res=Splot$res)
      dev.off()
    }
    if(tiffF){
      dev.copy(tiff, filename=fname.tiff, width=as.double(wi), height=as.double(hi), pointsize=Splot$pointsize, res=Splot$res)
      dev.off()
    }
    if(psF){
      dev.copy(postscript,file=fname.ps,paper=Splot$ps.paper,width=Splot$ps.width,height=Splot$ps.height,horizontal=FALSE, pointsize=Splot$pointsize)
      dev.off()
    }
    if(pdfF){
      dev.copy2pdf(file=fname.pdf,paper=Splot$ps.paper,width=Splot$ps.width,height=Splot$ps.height,pointsize=Splot$pointsize)
      #dev.off()
    }
    dev.off() 
    

  }else{

    jpeg(filename=fname.jpeg, width=as.double(wi), height=as.double(hi), pointsize=Splot$pointsize, res=Splot$res)
    dev.control("enable")

    
    Splot$plot.output = list()
  
    # initiate layout
    if(max(as.vector(Splot$mat),na.rm=TRUE)>1) nf = layout(Splot$mat, respect=TRUE)
    if(Splot$mai.prc) mai.def=par("mai")

    # initialize getLims 
    if(is.null(Splot$plot.lims)){
      plot.lims = list()
      plot.lims$xmins = rep(NA, length(Splot$plot.calls))
      plot.lims$xmaxs = rep(NA, length(Splot$plot.calls))
      plot.lims$ymins = rep(NA, length(Splot$plot.calls))
      plot.lims$ymaxs = rep(NA, length(Splot$plot.calls))
      Splot$plot.lims = plot.lims    
    }
  
    # loop over plot calls to place plots in order or 1:n in layout
    for(i in 1:length(Splot$plot.calls)){
      # set up plot margin
      if(length(Splot$mai.mat)>1){
        if(!Splot$mai.prc) par(mai=Splot$mai.mat[i,])
        if(Splot$mai.prc)  par(mai=Splot$mai.mat[i,]*mai.def)            
      }
      # evaluate plot call
      plt = eval.js(Splot$plot.calls[[i]])

      if(is.null(plt)) Splot$plot.output[i] = NA
      if(!is.null(plt)){
        Splot$plot.output[i] = NA
        class(Splot$plot.output[i]) = "list"
        Splot$plot.output[[i]] = plt
      }
    
      # cycle through all additional plot calls for current plot
      sub.np = length(Splot$plot.extras[[i]])
      for(sp in 1:sub.np){
        if(!is.na(Splot$plot.extras[[i]][[sp]]))
          eval.js(Splot$plot.extras[[i]][[sp]])
      }

      Splot$plot.lims$xmins[i] = par()$usr[1]
      Splot$plot.lims$xmaxs[i] = par()$usr[2]
      Splot$plot.lims$ymins[i] = par()$usr[3]
      Splot$plot.lims$ymaxs[i] = par()$usr[4]
           
    }# end for loop over plot calls


    if(tiffF){
      dev.copy(tiff, filename=fname.tiff, width=as.double(wi), height=as.double(hi), pointsize=Splot$pointsize, res=Splot$res)
      dev.off()
    }
    if(psF){
      dev.copy(postscript,file=fname.ps,paper=Splot$ps.paper,width=Splot$ps.width,height=Splot$ps.height,horizontal=FALSE, pointsize=Splot$pointsize)
      dev.off()
    }
    if(pdfF){
      dev.copy2pdf(file=fname.pdf,paper=Splot$ps.paper,width=Splot$ps.width,height=Splot$ps.height,pointsize=Splot$pointsize)
      #dev.off()
    }
    dev.off() 
        
  } # end else







  
  


  #
  # interactive webpage
  #
  
  if(makeInteractive){
    
    # check for interactive plots
    if(!is.na(overrideInteractive[1])){
      if(length(overrideInteractive) != Splot$nfig){
        cat(paste("Note: overrideInteractive is not of correct length\n       Length of overrideInteractve:", length(overrideInteractive),"\n       should be equal to the number of figures:", Splot$nfig, "\n       Continuing with originally set interactive plots \n"))
        overrideInteractive = NA            
      }
      if(class(overrideInteractive) != "logical"){
        cat(paste("Note: overrideInteractive is not correct \n       Must be Logical (T/F) vector of length:", Splot$nfig, "\n       Continuing with originally set interactive plots \n"))
        overrideInteractive = NA      
      }      
    }
    # if override is set use override otherwise use stored Iflag object
    if(is.na(overrideInteractive[1])) Ifig = which(Splot$Iflag)
    if(!is.na(overrideInteractive[1])) Ifig = which(overrideInteractive)

    # if there is at least one interactive plot
    if(length(Ifig) !=0){

      # checks to make sure at least one of the labelled interactive plots
      #  does have a mapping 
      Fsum = 0
      for(f in Ifig){
        Fsum = Fsum + length(Splot$iMap[[f]])
      }
      if(Fsum != 0){

     
           ##############################################################
           ##############################################################
           #
           # eventually need to run a check and combined
           #   duplicate coordinates
           #
           # print points before regions?? 
           #   
           ##############################################################
           ##############################################################
      
    
        #
        #  combined image mapping and make .html
        #


        # load header information for html file
        if(header!="v1" &  header!="v2" & header!="v3") header="v3"
        # mapfile header info
        if(header=="v2") sp.header= sp.header2 #data(v2.header)
        if(header=="v1") sp.header= sp.header1 #data(v1.header)
        
        # begin html file
        sink(fname.html)
        # add header information
        if(header=="v3"){
          w.wi = strsplit(window.size, "x")[[1]][1]
          w.hi = strsplit(window.size, "x")[[1]][2]

          cat(paste("<html>\n<head>\n    <title>sendplot</title>\n    <style type=\"text/css\">\n        /* CSS GOES HERE */\n        .plot {border:1px solid #CCCCCC;display:block;overflow:auto;padding:5px;\n               max-width:",w.wi,"px; max-height:",w.hi,"px;}\n    </style>\n</head>\n", sep=""))

          #data(v3.header)
          sp.header=sp.header3
        }

        sp.header=sp.header
        
        for(i in 1:length(sp.header)) cat(sp.header[i],fill=TRUE)
   
        
        # loop over labelled interactive
        # writing to file 
        for(fi in Ifig){
          
          lenList = length(Splot$iMap[[fi]])
          if(lenList != 0){
            
            for(ll in 1:lenList){

              obj = Splot$iMap[[fi]][[ll]]
              DFs = makeCharacter(obj)
              iType = Splot$iType[[fi]][ll]
              if(header=="v1") writeToHTML1(obj, DFs, iType)
              if(header=="v2" | header=="v3") writeToHTML2(obj, DFs, iType)
              
            }
          }
        }

        if(!is.null(Splot$Default)){
          if(Default){
            if(header=="v1") writeDefault1(Splot)
            if(header=="v2" | header=="v3") writeDefault2(Splot)
          }
        }
        
        cat("</map>",fill=TRUE)
        cat("<div class=\"plot\">",fill=TRUE)
        if(!pngF){

          image.name.jpeg=paste(fname.root,".jpeg",sep="")
 
          if(header=="v1")cat("<img border=\"0\" src=\"",image.name.jpeg,"\" usemap=\"img-map\" />",sep="",fill=TRUE)
          if(header=="v2" | header=="v3")cat("<img border=\"0\" src=\"",image.name.jpeg,"\" usemap=\"#img-map\" />",sep="",fill=TRUE)          
        }else{

          image.name.png=paste(fname.root,".png",sep="")
           
          if(header=="v1")cat("<img border=\"0\" src=\"",image.name.png,"\" usemap=\"img-map\" />",sep="",fill=TRUE)
          if(header=="v2" | header=="v3")cat("<img border=\"0\" src=\"",image.name.png,"\" usemap=\"#img-map\" />",sep="",fill=TRUE)
        }
        cat("</div>",fill=TRUE)
        cat("</body>",fill=TRUE)
        cat("</html>",fill=TRUE)
        
        sink()

        if(header=="v1") cat("Note: hyperlinks currenly only work with header=v2 \n")

        
      }else{
        cat("Note:  Plot[s] are labelled as interactive:\n       But none have been mapped\n       Please set using makeImap \n")
      }

    }else{

      # load header information for html file
      if(header!="v1" &  header!="v2" & header!="v3") header="v3"
      # mapfile header info
      if(header=="v2") sp.header= sp.header2 #data(v2.header)
      if(header=="v1") sp.header= sp.header1 #data(v1.header)
        
      # begin html file
      sink(fname.html)
      # add header information
      if(header=="v3"){
        w.wi = strsplit(window.size, "x")[[1]][1]
        w.hi = strsplit(window.size, "x")[[1]][2]

        cat(paste("<html>\n<head>\n    <title>sendplot</title>\n    <style type=\"text/css\">\n        /* CSS GOES HERE */\n        .plot {border:1px solid #CCCCCC;display:block;overflow:auto;padding:5px;\n               max-width:",w.wi,"px; max-height:",w.hi,"px;}\n    </style>\n</head>\n", sep=""))

        #data(v3.header)
        sp.header= sp.header3  
      }

      sp.header=sp.header
        
      for(i in 1:length(sp.header)) cat(sp.header[i],fill=TRUE)

      if(!is.null(Splot$Default)){
        if(Default){
          if(header=="v1") writeDefault1(Splot)
          if(header=="v2" | header=="v3") writeDefault2(Splot)
        }
      }
      
      cat("</map>",fill=TRUE)
      cat("<div class=\"plot\">",fill=TRUE)
      if(!pngF){
        
        image.name.jpeg=paste(fname.root,".jpeg",sep="")
        
        if(header=="v1")cat("<img border=\"0\" src=\"",image.name.jpeg,"\" usemap=\"img-map\" />",sep="",fill=TRUE)
        if(header=="v2" | header=="v3")cat("<img border=\"0\" src=\"",image.name.jpeg,"\" usemap=\"#img-map\" />",sep="",fill=TRUE)          
      }else{
        
        image.name.png=paste(fname.root,".png",sep="")
        
        if(header=="v1")cat("<img border=\"0\" src=\"",image.name.png,"\" usemap=\"img-map\" />",sep="",fill=TRUE)
        if(header=="v2" | header=="v3")cat("<img border=\"0\" src=\"",image.name.png,"\" usemap=\"#img-map\" />",sep="",fill=TRUE)
      }
      cat("</div>",fill=TRUE)
      cat("</body>",fill=TRUE)
      cat("</html>",fill=TRUE)
      
      sink()
      
      if(header=="v1") cat("Note: hyperlinks currenly only work with header=v2 \n")


      
      cat("Note:  No plots are designated interactive \n       Please set using makeImap \n")
    }
    
  }

  
  #if(getLims) returnObj = TRUE
  if(returnObj) return(Splot)
  
}




