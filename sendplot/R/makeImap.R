#
# acts on Splot object
#   gives a single figure interactive mapping
#



makeImap <- function(Splot,
                     figure=1,

                     xy.type=NA, # points, image.midpoints, image.boundaries, image.box, circle, rect, polygon 
                     x.pos,
                     y.pos,

                     x.right.pos=NA,  # if xy.type is rect
                     y.bottom.pos=NA, # if xy.type is rect
                     spot.radius = 5, # can be a vector, each entry in x.pos, # if xy.type=circle/points
         
                     x.labels=NA,
                     y.labels=NA,
                     xy.labels=NA,
                     
                     x.links=NA,
                     y.links=NA,
                     xy.links=NA,

                     asLinks=NA,

                     x.images=NA,
                     y.images=NA,
                     xy.images=NA,

                     sep.chr=":",

                     font.type="Helvetica", # 'Arial, Helvetica, sans-serif'   
                     font.color="black",  # name, or #------
                     font.size="12",      # can specify type px,pt,em, etc.
                     bg.color="#D6E3F6",  # name or #------
                     
                     
                     fname.root="Splot",
                     dir="./",

                     automap.method="mode",
                     
                     bb.clr=NA,
                     bb.cex=2,  # note: if image.size value is high this should be increase inorder for mapping to work
              
                     returnVl=TRUE,
                     saveFlag=FALSE,
                     saveName="Splot.RData",

                     cleanDir = TRUE

                      ){



  # first make normal static plot
  #  if necessary, retrieving xlim/ylims for all plots
  #                stores for future reference b/c limits are static
  if(is.null(Splot$plot.lims)){
    Splot = makeSplot(Splot, fname.root=fname.root, dir=dir, makeInteractive=FALSE, overwriteSourcePlot=c("png", "tiff"), returnObj=TRUE)
  }else{
    Splot = makeSplot(Splot, fname.root=fname.root, dir=dir, makeInteractive=FALSE, overwriteSourcePlot=c("png", "tiff"))
  }

  # make separate file with bound points
  # need to compare to original file 
  boundFileName = paste(fname.root, "Dot", sep="")

  

  ###############################
  ###############################
  #
  #  Start setting up automapping
  #
  ###############################
  ###############################

 # if(automap){
  
    detected = FALSE

    # loops over given colors, or 
    # if NA try: blue, red, black, white, green
    if(is.na(bb.clr[1])) bb.clr =  colors()[c(26,552,24,1,254)]  
    Cl.idx = 1
  
    # while mapping doesn't work 
    while(!detected){
    
      # if there are more colors to try 
      if(Cl.idx <= length(bb.clr)){

        # make file with bound points (always png)
        addBounding(Splot, figure, boundFileName=boundFileName, dir=dir, bb.clr=bb.clr[Cl.idx], bb.cex= bb.cex)
        # try to automatically find bounding points 
        boundingPt = automapPts(Splot, fname.root=fname.root, boundFileName=boundFileName, dir=dir, automap.method=automap.method)

        # if worked end loop
        if(class(boundingPt) != "try-error") detected=TRUE
        # if didn't work print error and try different color 
        if(class(boundingPt) == "try-error"){
          cat("Warning: First attempt at automapping failed \n\n")
          cat(paste(boundingPt))
          cat("\n\n         Attempting with different boundingPt color \n")
          Cl.idx = Cl.idx + 1
        }
      # if there are no more colors to try print Error   
      }else{
        detected = TRUE
        cat("ERROR: Could not map correctly \n")
      }
    }
#  }



  
  ###############################
  ###############################
  #
  #  Make figure mapping 
  #
  ###############################
  ###############################
 

  # if automap performed correctly
  if(class(boundingPt) != "try-error"){


    # get limits for specific figure
    xlim = c(Splot$plot.lims$xmins[figure], Splot$plot.lims$xmaxs[figure])
    ylim = c(Splot$plot.lims$ymins[figure], Splot$plot.lims$ymaxs[figure])

    # check how x and y values should be treated
    #    xy directly
    #    cuts as cuts as in image
    
    if(is.na(xy.type)){
      cat("Note:  xy.type not specified\n       Continuning with x.pos and y.pos as points\n")
      xy.type="points"
    }
    if( (xy.type != "points") & (xy.type != "image.midpoints") & (xy.type != "image.boundaries") & (xy.type != "image.box") &  (xy.type != "circle") &  (xy.type != "rect") &  (xy.type != "polygon") ){
      cat("Note:  xy.type is not acceptable \n       Continuning with x.pos and y.pos as points\n")
      xy.type="points"
    }


    
    # send to makeDF function for creation and alteration of data matrix entries 

    if( (xy.type=="points") | (xy.type=="circle") )  MapObj = makeScatterDF(Splot=Splot, xlim= xlim, ylim=ylim, x.pos=x.pos, y.pos=y.pos,boundingPt=boundingPt, x.labels=x.labels, y.labels=y.labels, xy.labels=xy.labels, x.links=x.links, y.links=y.links, xy.links=xy.links, asLinks=asLinks, x.images=x.images, y.images=y.images, xy.images=xy.images)

    if( (xy.type=="image.midpoints") | (xy.type=="image.boundaries") |  (xy.type=="image.box"))  MapObj = makeImageDF(Splot=Splot, xy.type=xy.type, xlim= xlim, ylim=ylim,x.pos=x.pos, y.pos=y.pos,boundingPt=boundingPt, x.labels=x.labels, y.labels=y.labels, xy.labels=xy.labels, x.links=x.links, y.links=y.links, xy.links=xy.links, asLinks=asLinks, x.images=x.images, y.images=y.images, xy.images=xy.images)

    if( xy.type=="rect")  MapObj = makeRectDF(Splot=Splot, xlim= xlim, ylim=ylim,x.left=x.pos, y.top=y.pos, x.right=x.right.pos, y.bottom=y.bottom.pos,  boundingPt=boundingPt, x.labels=x.labels, y.labels=y.labels, xy.labels=xy.labels, x.links=x.links, y.links=y.links, xy.links=xy.links, asLinks=asLinks, x.images=x.images, y.images=y.images, xy.images=xy.images)

    if( xy.type=="polygon" )  MapObj = makePolyDF(Splot=Splot, xlim= xlim, ylim=ylim, x.pos=x.pos, y.pos=y.pos,boundingPt=boundingPt, x.labels=x.labels, y.labels=y.labels, xy.labels=xy.labels, x.links=x.links, y.links=y.links, xy.links=xy.links, asLinks=asLinks, x.images=x.images, y.images=y.images, xy.images=xy.images)



    # format sep.chr
    if(class(sep.chr)=="character"){
      if(length(sep.chr) == 1) sep.chr = rep(sep.chr,9)
      if(length(sep.chr) < 9) sep.chr = c(sep.chr, rep(sep.chr[length(sep.chr)],9-length(sep.chr)))
      if(length(sep.chr) > 9) sep.chr = sep.chr[1:9]

      temp = list(x.lbl.sep=NA, y.lbl.sep=NA,  xy.lbl.sep=NA, x.lnk.sep=NA, y.lnk.sep=NA, xy.lnk.sep=NA, x.img.sep=NA ,y.img.sep=NA,  xy.img.sep=NA)
      temp.vec = c("x.labels","y.labels","xy.labels","x.links","y.links","xy.links", "x.images","y.images","xy.images")
      for(i in 1:9){
        temp.vl=NA
        eval.js(paste("temp.vl=", temp.vec[i]))
        if(class(temp.vl)!="data.frame"){
          if(length(temp.vl) > 1){
            eval.js(paste("temp$",names(temp)[i]," = '", sep.chr[i], "'", sep=""))
          }else{
            if(class(temp.vl)=="logical") temp[i] = NA
            else eval.js(paste("temp$",names(temp)[i]," = '", sep.chr[i], "'",sep=""))
          }
        }else{
          eval.js(paste("temp$",names(temp)[i]," = rep(sep.chr[i], dim(temp.vl)[2])",sep=""))        
        }
      }
      sep.chr = temp
    }
    if(class(sep.chr) == "list"){
      temp.vec = c("x.lbl.sep","y.lbl.sep","xy.lbl.sep","x.lnk.sep","y.lnk.sep","xy.lnk.sep", "x.img.sep","y.img.sep","xy.img.sep")
      nm = names(sep.chr)
      mis = setdiff(temp.vec, nm)
      if(length(mis) != 0){
        for(m in mis){
          eval.js(paste("sep.chr$", m ,"=NA",sep=""))
        }
      }
      temp.vec2 = c("x.labels","y.labels","xy.labels","x.links","y.links","xy.links", "x.images","y.images","xy.images")
      log.vec = rep(NA,9)
      for(i in 1:9){
        if(eval.js(paste("class(", temp.vec2[i], ")=='data.frame'",sep=""))){
          log.vec[i] = eval.js(paste("class(", temp.vec2[i], ")=='data.frame'",sep=""))
        }else{
          if(eval.js(paste("length(", temp.vec2[i], ") > 1",sep=""))){
            log.vec[i] = TRUE
          }else{
            log.vec[i] = eval.js(paste("!is.na(", temp.vec2[i], ")",sep=""))
          }          
        }
      }
      idx = which(log.vec)
      if(length(idx) != 0){
        for(j in idx){


          if(eval.js(paste("class(", temp.vec2[j],")=='list'", sep=""))) eval.js(paste( temp.vec2[j], "= as.data.frame(", temp.vec2[j], ")", sep=""))
          
          if(eval.js(paste("class(", temp.vec2[j],")=='data.frame'", sep=""))){
            dm.df = eval.js(paste("dim(", temp.vec2[j],")[2]", sep=""))
            ln.sep = eval.js(paste("length(sep.chr$",temp.vec[j],")", sep=""))
            if(dm.df > ln.sep){
              eval.js(paste("sep.chr$",temp.vec[j],"= c(sep.chr$", temp.vec[j], ",", "rep(","sep.chr$",temp.vec[j],"[length(","sep.chr$",temp.vec[j] ,")],",dm.df-ln.sep,"))", sep=""))
            }
            if(dm.df < ln.sep){
              eval.js(paste("sep.chr$",temp.vec[j],"= sep.chr$", temp.vec[j], "[1:",dm.df,"]", sep=""))
              
            }
            cng = eval.js(paste("which(is.na(sep.chr$",temp.vec[j], "))", sep=""))
            eval.js(paste("sep.chr$",temp.vec[j],"[cng] = ':'", sep=""))

          }else{
       
            if(eval.js(paste("length(", temp.vec2[j],") > 1", sep=""))){
              eval.js(paste("sep.chr$",temp.vec[j],"= sep.chr$", temp.vec[j], sep=""))              
            }else{
              if(eval.js(paste("!is.na(", temp.vec2[j], ")",sep=""))) eval.js(paste("sep.chr$",temp.vec[j],"=sep.chr$", temp.vec[j], "[1]",sep=""))
            }
          }          
        }        
      }      
    }
    
    # add to mapobj

     MapObj$sep.chr = sep.chr

    
    
    if(length(MapObj) != 1){

      if((xy.type=="points") | (xy.type=="circle")){
        if(length(spot.radius) == 1) spot.radius = rep(spot.radius, length(x.pos))
        if(length(spot.radius) < length(x.pos))  spot.radius = c(spot.radius, rep(5, (length(x.pos)-length(spot.radius))))
      }
      if(xy.type=="image.midpoints"){
        len = length(x.pos)*length(y.pos)
        if(length(spot.radius) == 1) spot.radius = rep(spot.radius, len)
        if(length(spot.radius) < len)  spot.radius = c(spot.radius, rep(5, (len-length(spot.radius))))
      }
      if(xy.type=="image.boundaries"){
        len = (length(x.pos)-1)*(length(y.pos)-1)
        if(length(spot.radius) == 1) spot.radius = rep(spot.radius, len)
        if(length(spot.radius) < len)  spot.radius = c(spot.radius, rep(5, (len-length(spot.radius))))
      }

      MapObj$spot.radius = spot.radius

      # add display information
      MapObj$font.type  = font.type
      MapObj$font.color = font.color
      MapObj$font.size = font.size
      MapObj$bg.color = bg.color
                   
      
      
      # add MapObj to Splot
      if(length(Splot$iMap[[figure]]) == 0){
        eval.js(paste("Splot$iMaps$Figure",figure,"$MapObj1 = MapObj", sep=""))
      }else{
        len = length(Splot$iMap[[figure]])
        eval.js(paste("Splot$iMaps$Figure",figure,"$MapObj",len+1," = MapObj", sep=""))
      }
      
    
      # add iType to Splot
      if(length(Splot$iTypes[[figure]]) == 0){
        eval.js(paste("Splot$iTypes$Figure",figure," = MapObj$xy.type", sep=""))
      }else{
        len = length(Splot$iTypes[[figure]])
        eval.js(paste("Splot$iTypes$Figure",figure,"[len+1] =MapObj$xy.type", sep=""))
                                        #c(Splot$iTypes$Figure",figure,"[[1]] ,MapObj$xy.type)", sep=""))
      }
  
  
      # update Iflag
      if(!Splot$Iflag[figure]){
        cat("Note: Figure currently not labelled as interactive \n      Changing to Interactive now \n")
        Splot$Iflag[figure]=TRUE
      }
    }


    
    
  }else{# end if automap worked if(class(boundingPt) !="try-error")

    cat("ERROR: Could not automap points \n       Mapping NOT performed \n       Please try again using different bb.cex/bb.clr settings \n       Or using manual mapping \n")
    
  }

  
  if(cleanDir){
    
    if(Splot$platform == "unix"){

      system("mkdir makeImapTempDir", ignore.stderr =TRUE)
      system(paste("mv ",dir,"*Dot* makeImapTempDir/",sep=""), ignore.stderr =TRUE)
      system(paste("mv ",dir,fname.root,".tif makeImapTempDir/", sep=""), ignore.stderr =TRUE)
      system("rm -r makeImapTempDir", ignore.stderr =TRUE)
      
      
    }else{
      
      shell("mkdir makeImapTempDir", mustWork=NA,ignore.stderr =TRUE)
      shell(paste("mv ",dir,"*Dot* makeImapTempDir/",sep=""), mustWork=NA,ignore.stderr =TRUE)
      shell(paste("mv ",dir,fname.root,".tif makeImapTempDir/", sep=""), mustWork=NA,ignore.stderr =TRUE)
      #shell("rm -r makeImapTempDir", mustWork=NA,ignore.stderr =TRUE)
      shell("del /Q makeImapTempDir", mustWork=NA,ignore.stderr =TRUE)
      shell("rd /Q makeImapTempDir", mustWork=NA,ignore.stderr =TRUE)
      
    }

    
  }


  
  #########################################
  #########################################
  #
  #  save  mapping and save/return object
  #
  #########################################
  #########################################

  

  # save and return
  if(saveFlag) save(Splot, file=saveName, compress=TRUE)
  if(returnVl) return(Splot)


}
