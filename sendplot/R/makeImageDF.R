#
# maps x and y values as cuts
#  checks and alters lbls/links matricies
#



makeImageDF <- function(Splot,xy.type,
                        xlim, ylim,
                        x.pos,y.pos,
                        boundingPt,
                        x.labels=NA,
                        y.labels=NA,
                        xy.labels=NA,
                        x.links=NA,
                        y.links=NA,
                        xy.links=NA,
                        asLinks=NA,
                        x.images=NA,
                        y.images=NA,
                        xy.images=NA
                        ){


  up.left = boundingPt$up.left
  low.right = boundingPt$low.right
  
  # create object to return
  MapObj = list()


  if(xy.type == "image.midpoints"){
    if(xlim[1] == xlim[2]){
      xdif = 1
    }else{
      xdif = xlim[2]-xlim[1]
    }
    if(low.right[1] == up.left[1]){
      ptdif = 1
    }else{
      ptdif = low.right[1]-up.left[1]
    }
    if(ylim[1] == ylim[2]){
      ydif = 1
    }else{
      ydif = ylim[2]-ylim[1]
    }
    if(low.right[2] == up.left[2]){
      ptdif2 = 1
    }else{
      ptdif2 = low.right[2]-up.left[2]
    }

    #x.image =  round(up.left[1] + ((x.pos-xlim[1])/(xlim[2]-xlim[1]))*(low.right[1]-up.left[1]))
    #y.image = round(up.left[2] + ((ylim[2]-y.pos)/(ylim[2]-ylim[1]))*(low.right[2]-up.left[2]))
    x.image = round(up.left[1] + ((x.pos-xlim[1])/(xdif))*(ptdif))
    y.image = round(up.left[2] + ((ylim[2]-y.pos)/(ydif))*(ptdif2))



    
    # initiate data frame of info
    dat = data.frame(
      # pixil locations
      pix.x = as.vector(mapply(rep,x=x.image,MoreArgs=list(times=length(y.image)))),
      pix.y = rep(y.image, length(x.image))
      )
    dat2 = data.frame(rep(NA, (length(y.image)*length(x.image))))
    names(dat2) = "tempNA"

   # initiate data frame for images 
    dat3 = data.frame(rep(NA, (length(y.image)*length(x.image))))
    names(dat3) = "tempNA"


    
    MapObj$xy.type = "circle"
  }
  
  if(xy.type == "image.boundaries"){

    # find midpoints
    new.x = cntrs=(x.pos[1:(length(x.pos)-1)]+x.pos[2:length(x.pos)])/2
    new.y = cntrs=(y.pos[1:(length(y.pos)-1)]+y.pos[2:length(y.pos)])/2
    x.image =  round(up.left[1] + ((new.x-xlim[1])/(xlim[2]-xlim[1]))*(low.right[1]-up.left[1]))
    y.image = round(up.left[2] + ((ylim[2]-new.y)/(ylim[2]-ylim[1]))*(low.right[2]-up.left[2]))
    
    # initiate data frame of info
    dat = data.frame(
      # pixil locations
      pix.x = as.vector(mapply(rep,x=x.image,MoreArgs=list(times=length(y.image)))),
      pix.y = rep(y.image, length(x.image))
      )
    dat2 = data.frame(rep(NA, (length(y.image)*length(x.image))))
    names(dat2) = "tempNA"
    dat3 = data.frame(rep(NA, (length(y.image)*length(x.image))))
    names(dat3) = "tempNA"



    
    MapObj$xy.type = "circle"
  }

    
  if(xy.type == "image.box"){

    # assumes x.pos and y.pos are boundaries

    # map boxes low left and upper right
    x.left = rep(c(x.pos[1]:x.pos[length(x.pos)-1]), each=(length(y.pos)-1))
    y.bottom = rep(c(y.pos[1]:y.pos[length(y.pos)-1]), times=(length(x.pos)-1))
    x.right = rep(c(x.pos[2]:x.pos[length(x.pos)]), each=(length(y.pos)-1))
    y.top = rep(c(y.pos[2]:y.pos[length(y.pos)]), times=(length(x.pos)-1))

    # covert to pixels 
    new.x.left =  round(up.left[1] + ((x.left-xlim[1])/(xlim[2]-xlim[1]))*(low.right[1]-up.left[1]))
    new.y.bottom =  round(up.left[2] + ((ylim[2]-y.bottom)/(ylim[2]-ylim[1]))*(low.right[2]-up.left[2]))
    new.x.right = round(up.left[1] + ((x.right-xlim[1])/(xlim[2]-xlim[1]))*(low.right[1]-up.left[1]))
    new.y.top = round(up.left[2] + ((ylim[2]-y.top)/(ylim[2]-ylim[1]))*(low.right[2]-up.left[2]))

    # for matrix dimension checks
    x.image = x.pos[1:(length(x.pos)-1)]
    y.image = y.pos[1:(length(y.pos)-1)]
    
    dat = data.frame(
      # pixil locations
      pix.x.left = new.x.left,
      pix.y.top = new.y.top,
      pix.x.right = new.x.right,
      pix.y.bottom = new.y.bottom      
      )
    dat2 = data.frame(rep(NA, (length(y.image)*length(x.image))))
    names(dat2) = "tempNA"
    dat3 = data.frame(rep(NA, (length(y.image)*length(x.image))))
    names(dat3) = "tempNA"


    MapObj$xy.type = "rect"
  }
  





#######################################################  
#######################################################  
 
  if(FALSE){
  # legacy code
    # this code essential finds boundaries but does nothing
    # because assumes want midpoints (which is given)
    
    wdth=low.right[1]-up.left[1]
    hght=low.right[2]-up.left[2]
    # get min and max x values for image      
    nx=length(x.pos)
    xmin=x.pos[1]-(x.pos[2]-x.pos[1])/2
    xmax=x.pos[nx]+(x.pos[nx]-x.pos[nx-1])/2
    # calculate cuts and center of x values on image
    bndrs=c(xmin,(x.pos[1:(nx-1)]+x.pos[2:nx])/2,xmax)
    cntrs=(bndrs[1:(length(bndrs)-1)]+bndrs[2:length(bndrs)])/2
    # adjust 
    unit.int.wdth=cntrs-xmin
    unit.int.wdth=unit.int.wdth/(xmax-xmin)
    # calculate pixil positions
    x.image=round(unit.int.wdth*wdth+up.left[1])
    # get min and max y values for image
    ny=length(y.pos)
    ymin=y.pos[1]-(y.pos[2]-y.pos[1])/2
    ymax=y.pos[ny]+(y.pos[ny]-y.pos[ny-1])/2
    # calculate cuts and centers of y values on image
    bndrs=c(ymin,(y.pos[1:(ny-1)]+y.pos[2:ny])/2,ymax)
    cntrs=(bndrs[1:(length(bndrs)-1)]+bndrs[2:length(bndrs)])/2
    # adjust
    unit.int.hght=cntrs-ymin
    unit.int.hght=unit.int.hght/(ymax-ymin)
    # calculate pixil positions
    y.image= round((1-unit.int.hght)*hght+up.left[2])

  }

#######################################################  
#######################################################  






  






  
 # xy specific data
 # eval.js(paste("dat$",z.value,"=as.vector(z)",sep=""))

 # x specific data
  cont = TRUE
  x.labels = as.data.frame(x.labels)
  cngName =  grep("if ", names(x.labels))
  names(x.labels)[cngName] = paste("Value", cngName, sep="")
  names(x.labels) = gsub(pattern=" ", replacement=".",names(x.labels))
  if( (dim(x.labels)[1]==1) & (dim(x.labels)[2]==1)){
    if(is.na(x.labels[1,1])) cont = FALSE
  }
  # dimension check
  if(cont){
    if(dim(x.labels)[1] != length(x.image)){
      cont = FALSE
      cat(paste("Warning: x.labels does not have correct dimensions \n   number of rows should equal length(x.image):",length(x.image), "\n   Continuing with x.labels = NA \n", sep=""))
      x.labels = NA      
    }         
  }       
  # if x.labels is not NA continue
  if(cont){
    lev = levels(factor(names(x.labels)))
    num = length(lev)
    for(i in 1:num){
      idx = which(names(x.labels) == lev[i])
      if(length(idx) > 1){
        for(j in 2:length(idx)){
          vch=""
          for(k in 1:(j-1)){
            vch = paste(vch, ".", sep="")
          }
          names(x.labels)[idx[j]] = paste(names(x.labels)[idx[j]], vch, sep="")
        }        
      }
    }    
  
    for(i in 1:dim(x.labels)[2]){
      eval.js(paste("dat$",names(x.labels)[i], "=as.vector(mapply(rep,x=x.labels[,i], MoreArgs=list(times=length(y.image))))", sep=""))
    }
  }
  # y specific data
  cont = TRUE
  y.labels = as.data.frame(y.labels)
  cngName =  grep("if ", names(y.labels))
  names(y.labels)[cngName] = paste("Value", cngName, sep="")
  names(y.labels) = gsub(pattern=" ", replacement=".",names(y.labels))
  if( (dim(y.labels)[1]==1) & (dim(y.labels)[2]==1)){
    if(is.na(y.labels[1,1])) cont = FALSE
  }
  # dimension check
  if(cont){
    if(dim(y.labels)[1] != length(y.image)){
      cont = FALSE
      cat(paste("Warning: y.labels does not have correct dimensions \n   number of rows should equal length(y.image):",length(y.image), "\n   Continuing with y.labels = NA \n", sep=""))
      y.labels = NA      
    }         
  }    
  # if y.labels is not NA continue
  if(cont){
    lev = levels(factor(names(y.labels)))
    num = length(lev)
    for(i in 1:num){
      idx = which(names(y.labels) == lev[i])
      if(length(idx) > 1){
        for(j in 2:length(idx)){
          vch=""
          for(k in 1:(j-1)){
            vch = paste(vch, ".", sep="")
          }
          names(y.labels)[idx[j]] = paste(names(y.labels)[idx[j]], vch, sep="")
        }        
      }
    }        
   
    for(i in 1:dim(y.labels)[2]){
      eval.js(paste("dat$",names(y.labels)[i], "=rep(y.labels[,i],length(x.image))", sep=""))
    }
  }  
  cont = TRUE
  if(is.na(xy.labels[1])) cont = FALSE
  # if xy.labels is not NA continue
  if(cont){
    lev = levels(factor(names(xy.labels)))
    num = length(lev)
    for(i in 1:num){
      idx = which(names(xy.labels) == lev[i])
      if(length(idx) > 1){
        for(j in 2:length(idx)){
          vch=""
          for(k in 1:(j-1)){
            vch = paste(vch, ".", sep="")
          }
          names(xy.labels)[idx[j]] = paste(names(xy.labels)[idx[j]], vch, sep="")
        }        
      }
    }  
    
    for(i in 1:length(xy.labels)){
      # check dimension
      if((dim(xy.labels[[i]])[2] == length(x.image)) & (dim(xy.labels[[i]])[1] == length(y.image))){                  
        eval.js(paste("dat$",names(xy.labels)[i],"=as.vector(xy.labels[[i]])", sep=""))
      }else{
        cat(paste("Warning: at least one of the xy.labels matricies are not of the correct dimension. \n    All should be of the dimension ",length(x.image), " by ", length(y.image), "\n", sep="")) 
      }                  
    }
  }

  #################
  #
  # hyperlinks
  #
  #################

  
  # x specific hyperlinks
  cont = TRUE
  x.links = as.data.frame(x.links)
  cngName =  grep("if ", names(x.links))
  names(x.links)[cngName] = paste("Value", cngName, sep="")
  names(x.links) = gsub(pattern=" ", replacement=".",names(x.links))
  if( (dim(x.links)[1]==1) & (dim(x.links)[2]==1)){
    if(is.na(x.links[1,1])) cont = FALSE
  }    
  # dimension check
  if(cont){
    if(dim(x.links)[1] != length(x.image)){
      cont = FALSE
      cat(paste("Warning: x.links does not have correct dimensions \n   number of rows should equal length(x.image):",length(x.image), "\n   Continuing with x.links = NA \n", sep=""))
      x.links = NA      
    }         
  }
  # if x.links is not NA
  if(cont){
    lev = levels(factor(names(x.links)))
    num = length(lev)
    for(i in 1:num){
      idx = which(names(x.links) == lev[i])
      if(length(idx) > 1){
        for(j in 2:length(idx)){
          vch=""
          for(k in 1:(j-1)){
            vch = paste(vch, ".", sep="")
          }
          names(x.links)[idx[j]] = paste(names(x.links)[idx[j]], vch, sep="")
        }        
      }
    }    
    
    # for each column get links
    for(i in 1:dim(x.links)[2]){           
      eval.js("temp=as.vector(mapply(rep,x=x.links[,i], MoreArgs=list(times=length(y.image))))")
      # for each points link
      for(j in 1:length(temp)){
        tmp = temp[j]
        # if not NA
        if(is.na(tmp)){
          temp[j] = NA
          # split multiple links...assumed seperated by a comma
        }else{
          links = strsplit(tmp, split=",")[[1]]
          new.t = " "
          for(l in 1:length(links)){
            new.t = paste(new.t, paste("<a href=\\'", gsub(links[l], pattern=" ", replacement=""), "\\'> ", paste((names(x.links)[i]),l, sep="."), " </a>", sep=""), sep=",")          
          }
          new.t = gsub(new.t, pattern=" ,", replacement="")
          temp[j] = new.t
        }
      }
      # put correctly syntaxed link in matrix 
      eval.js(paste("dat2$", names(x.links)[i], "=temp", sep=""))    
    }  
  }      
  # y specific hyperlinks
  cont = TRUE
  y.links = as.data.frame(y.links)
  cngName =  grep("if ", names(y.links))
  names(y.links)[cngName] = paste("Value", cngName, sep="")
  names(y.links) = gsub(pattern=" ", replacement=".",names(y.links))
  if( (dim(y.links)[1]==1) & (dim(y.links)[2]==1)){
    if(is.na(y.links[1,1])) cont = FALSE
  }
  # dimension check
  if(cont){
    if(dim(y.links)[1] != length(y.image)){
      cont = FALSE
      cat(paste("Warning: y.links does not have correct dimensions \n   number of rows should equal length(y.image):",length(y.image), "\n   Continuing with y.links = NA \n", sep=""))
      y.links = NA      
    }         
  }
  # if y.links is not NA
  if(cont){

    lev = levels(factor(names(y.links)))
    num = length(lev)
    for(i in 1:num){
      idx = which(names(y.links) == lev[i])
      if(length(idx) > 1){
        for(j in 2:length(idx)){
          vch=""
          for(k in 1:(j-1)){
            vch = paste(vch, ".", sep="")
          }
          names(y.links)[idx[j]] = paste(names(y.links)[idx[j]], vch, sep="")
        }        
      }
    }    
    
    # for each column get links
    for(i in 1:dim(y.links)[2]){
      eval.js("temp=as.vector(rep(y.links[,i],length(x.image)))")
      # for each points link
      for(j in 1:length(temp)){
        tmp = temp[j]
        # if not NA
        if(is.na(tmp)){
          temp[j] = NA
        # split multiple links...assumed seperated by a comma
        }else{
          links = strsplit(tmp, split=",")[[1]]
          new.t = " "
          for(l in 1:length(links)){
            new.t = paste(new.t, paste("<a href=\\'", gsub(links[l], pattern=" ", replacement=""), "\\'> ", paste((names(y.links)[i]),l, sep="."), " </a>", sep=""), sep=",")          
          }
          new.t = gsub(new.t, pattern=" ,", replacement="")
          temp[j] = new.t
        }
      }
      # put correctly syntaxed link in matrix 
      eval.js(paste("dat2$", names(y.links)[i], "=temp", sep=""))    
    }  
  }
  # xy specific hyperlinks
  cont = TRUE
  if(is.na(xy.links[1])) cont = FALSE
  # if xy.links is not NA
  if(cont){
    lev = levels(factor(names(xy.links)))
    num = length(lev)
    for(i in 1:num){
      idx = which(names(xy.links) == lev[i])
      if(length(idx) > 1){
        for(j in 2:length(idx)){
          vch=""
          for(k in 1:(j-1)){
            vch = paste(vch, ".", sep="")
          }
          names(xy.links)[idx[j]] = paste(names(xy.links)[idx[j]], vch, sep="")
        }        
      }
    }    
  
    # for each matrix of links
    for(i in 1:length(xy.links)){
      eval.js("temp=xy.links[[i]]")
      # check dimensions
      if((dim(temp)[2] == length(x.image)) & (dim(temp)[1] == length(y.image))){
        temp = as.vector(temp)
        # for each points link
        for(j in 1:length(temp)){
          tmp = temp[j]
          # if not NA
          if(is.na(tmp)){
            temp[j] = NA
            # split multiple links...assumed seperated by a comma  
          }else{
            links = strsplit(tmp, split=",")[[1]]
            new.t = " "
            for(l in 1:length(links)){
              new.t = paste(new.t, paste("<a href=\\'", gsub(links[l], pattern=" ", replacement=""), "\\'> ", paste((names(xy.links)[i]),l, sep="."), " </a>", sep=""), sep=",")          
            }
            new.t = gsub(new.t, pattern=" ,", replacement="")
            temp[j] = new.t
          }
        }
        # put correctly syntaxed link in matrix 
        eval.js(paste("dat2$", names(xy.links)[i], "=temp", sep=""))    
      }else{
        cat(paste("Warning: at least one of the xy.links matricies are not of the correct dimension. \n    All should be of the dimension ",length(x.image), " by ", length(y.image), "\n", sep="")) 
      }          
    }        
  }      
      
  # get points as Links information 
  contLinks = TRUE
  # if data frame convert to matrix
  if(class(asLinks) == "data.frame") asLinks = as.matrix(asLinks)
  # if matrix convert to vector
  if(class(asLinks) == "matrix") asLinks = as.vector(asLinks)
  # repeat values if necessary
  if(length(asLinks) == length(x.image)) asLinks = rep(asLinks, each=length(y.image))
  if(length(asLinks) == length(y.image)) asLinks = rep(asLinks, length(x.image))
  if((length(asLinks) == 1) & !is.na(asLinks[1])) asLinks = rep(asLinks, (length(x.image)*length(y.image)))
  # convert to character vector for easy access
  asLinks = as.character(asLinks)
  # check dimension
  if((length(asLinks) != (length(x.image)*length(y.image))) & !is.na(asLinks[1])){
    cat("Warning: cannot create points as links \n     length must be equal to x or y or dimensions equal to x*y \n")
    contLinks = FALSE
  }      
  if(length(asLinks) ==1){
    if(is.na(asLinks[1])) contLinks=FALSE
  }

 
  #
  # images
  #


  #
  # x specific data
  #
  contxi = TRUE
  x.images = as.data.frame(x.images)
  cngName =  grep("if ", names(x.images))
  names(x.images)[cngName] = paste("Value", cngName, sep="")
  names(x.images) = gsub(pattern=" ", replacement=".",names(x.images))
  if( (dim(x.images)[1]==1) & (dim(x.images)[2]==1)){
    if(is.na(x.images[1,1])) contxi = FALSE
  }
  # dimension check
  if(contxi){
    if((dim(x.images)[1] != length(x.image))){
      contxi = FALSE
      cat(paste("Warning: x.images does not have correct dimensions \n   number of rows should equal length(x.image):",length(x.image), "\n   Continuing with x.images = NA \n", sep=""))
      x.images = NA
    }
  }
  # if x.images is not NA continue
  if(contxi){
    lev = levels(factor(names(x.images)))
    num = length(lev)
    for(i in 1:num){
      idx = which(names(x.images) == lev[i])
      if(length(idx) > 1){
        for(j in 2:length(idx)){
          vch=""
          for(k in 1:(j-1)){
            vch = paste(vch, ".", sep="")
          }
          names(x.images)[idx[j]] = paste(names(x.images)[idx[j]], vch, sep="")
        }        
      }
    }    
   
    for(i in 1:dim(x.images)[2]){
      eval.js("temp=as.vector(mapply(rep,x=x.images[,i], MoreArgs=list(times=length(y.image))))")
      # for each points link
      for(j in 1:length(temp)){
        tmp = temp[j]
        # if not NA
        if(is.na(tmp)){
          temp[j] = NA
        }else{
          new.ti= paste("<img src=\\'",tmp,"\\'>", sep="")          
          temp[j] = new.ti
        
        }
      }
     # put link in correct syntax into character matrix
      eval.js(paste("dat3$", names(x.images)[i], "=temp", sep=""))    
    }
  }


  

  #
  # y specific data
  #
  contyi = TRUE
  y.images = as.data.frame(y.images)
  cngName =  grep("if ", names(y.images))
  names(y.images)[cngName] = paste("Value", cngName, sep="")
  names(y.images) = gsub(pattern=" ", replacement=".",names(y.images))
  if( (dim(y.images)[1]==1) & (dim(y.images)[2]==1)){
    if(is.na(y.images[1,1])) contyi = FALSE
  }
  # dimension check
  if((dim(y.images)[1] != length(y.image)) & contyi){
    contyi = FALSE
    cat(paste("Warning: y.images does not have correct dimensions \n   number of rows should equal length(y.image):",length(y.image), "\n   Continuing with y.images = NA \n", sep=""))
    y.images = NA
  }      
  # if y.images is not NA continue
  if(contyi){

   lev = levels(factor(names(y.images)))
    num = length(lev)
    for(i in 1:num){
      idx = which(names(y.images) == lev[i])
      if(length(idx) > 1){
        for(j in 2:length(idx)){
          vch=""
          for(k in 1:(j-1)){
            vch = paste(vch, ".", sep="")
          }
          names(y.images)[idx[j]] = paste(names(y.images)[idx[j]], vch, sep="")
        }        
      }
    }    
    
    for(i in 1:dim(y.images)[2]){
      eval.js("temp=as.vector(rep(y.images[,i],length(x.image)))")
      # for each points link
      for(j in 1:length(temp)){
        tmp = temp[j]
        # if not NA
        if(is.na(tmp)){
          temp[j] = NA
        }else{
          new.ti= paste("<img src=\\'",tmp,"\\'>", sep="")          
          temp[j] = new.ti
        
        }
      }
     # put link in correct syntax into character matrix
      eval.js(paste("dat3$", names(y.images)[i], "=temp", sep=""))    
    }
  }


  
  #
  # xy -- assumes in this case that columns are different data vectors of row == nsmpls
  #
  contxyi = TRUE
  if(is.na(xy.images[1])) contxyi = FALSE
  # if xy.images is not NA continue
  if(contxyi){

  lev = levels(factor(names(xy.images)))
    num = length(lev)
    for(i in 1:num){
      idx = which(names(xy.images) == lev[i])
      if(length(idx) > 1){
        for(j in 2:length(idx)){
          vch=""
          for(k in 1:(j-1)){
            vch = paste(vch, ".", sep="")
          }
          names(xy.images)[idx[j]] = paste(names(xy.images)[idx[j]], vch, sep="")
        }        
      }
    }    

    
    for(i in 1:length(xy.images)){
      eval.js("temp=xy.images[[i]]")
      if((dim(temp)[2] == length(x.image)) & (dim(temp)[1] == length(y.image))){
        # for each points link
        temp = as.vector(temp)
        for(j in 1:length(temp)){
          tmp = temp[j]
        # if not NA
          if(is.na(tmp)){
            temp[j] = NA
          }else{
  
            new.ti= paste("<img src=\\'",tmp,"\\'>", sep="")          
            temp[j] = new.ti
            
          }
        }
     # put link in correct syntax into character matrix
        eval.js(paste("dat3$", names(xy.images)[i], "=temp", sep=""))    
      }else{
   cat(paste("Warning: at least one of the xy.images matricies are not of the correct dimension. \n    All should be of the dimension ",length(x.image), " by ", length(y.image), "\n", sep="")) 
    
      }
    }   
  }
       



  

  
  # add to object to return
  if(xy.type != "image.box") MapObj$Pixcoord = paste(as.character(x.image), as.character(y.image), sep=",")
  MapObj$dat = dat
  MapObj$dat2 = dat2
  MapObj$dat3 = dat3
  MapObj$contLinks = contLinks
  MapObj$asLinks = asLinks

  # return object
  return(MapObj)

  
}

