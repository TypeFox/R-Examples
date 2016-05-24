makePolyDF <-function(Splot,
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


  if( (x.pos[1] != x.pos[length(x.pos)]) & (y.pos[1] != y.pos[length(y.pos)]) ) {
    x.pos = c(x.pos, x.pos[1])
    y.pos = c(y.pos, y.pos[1])
  }


  
  # estimate pixil location
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
  #x.new = round(up.left[1] + ((x.pos-xlim[1])/(xlim[2]-xlim[1]))*(low.right[1]-up.left[1]))
  #y.new = round(up.left[2] + ((ylim[2]-y.pos)/(ylim[2]-ylim[1]))*(low.right[2]-up.left[2]))
  x.new = round(up.left[1] + ((x.pos-xlim[1])/(xdif))*(ptdif))
  y.new = round(up.left[2] + ((ylim[2]-y.pos)/(ydif))*(ptdif2))


  vertices = rep(NA, (length(x.new)*2))
  vertices[(1:length(x.new))*2] = y.new
  vertices[(((1:length(x.new))*2)-1)] = x.new

  tempDat = matrix(vertices, ncol=length(vertices))
  
  # initiate data frame for lbls
  dat = as.data.frame(tempDat)
  names(dat) = paste("Coords", c(1:length(names(dat))), sep="")

  # initiate data frame for links 
  dat2 = data.frame(rep(NA, (dim(dat)[1])))
  names(dat2) = "tempNA"

  # initiate data frame for images
  dat3 = data.frame(rep(NA, (dim(dat)[1])))
  names(dat3) = "tempNA"


  #######################
  #######################
  #
  #  check and fill lbls
  #
  #######################
  #######################
 
  #
  # x specific data
  #
  contx = TRUE
  x.labels = as.data.frame(x.labels)
  cngName =  grep("if ", names(x.labels))
  names(x.labels)[cngName] = paste("Value", cngName, sep="")        
  names(x.labels) = gsub(pattern=" ", replacement=".",names(x.labels))
  if( (dim(x.labels)[1]==1) & (dim(x.labels)[2]==1)){
    if(is.na(x.labels[1,1])) contx = FALSE
  }
  # dimension check
  if(contx){
    if((dim(x.labels)[1] != 1) & contx){
      contx = FALSE
      cat(paste("Warning: x.labels does not have correct dimensions \n   number of rows should be 1\n   Different variables should be in columns\n   Continuing with x.labels = NA \n", sep=""))
      x.labels = NA
    }
  }
  # if x.labels is not NA continue
  if(contx){
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
      if(i == 1) z.value = names(x.labels)[i]
      eval.js(paste("dat$",names(x.labels)[i], "=as.vector(x.labels[,i])", sep=""))
    }      
  }
  #
  # y specific data
  #
  conty = TRUE
  y.labels = as.data.frame(y.labels)
  cngName =  grep("if ", names(y.labels))
  names(y.labels)[cngName] = paste("Value", cngName, sep="")        
  names(y.labels) = gsub(pattern=" ", replacement=".",names(y.labels))
  if( (dim(y.labels)[1]==1) & (dim(y.labels)[2]==1)){
    if(is.na(y.labels[1,1])) conty = FALSE
  }
  # dimension check
  if((dim(y.labels)[1] != 1) & conty){
    conty = FALSE
    cat(paste("Warning: y.labels does not have correct dimensions \n   number of rows should equal 1\n   Different variables should be in columns\n   Continuing with y.labels = NA \n", sep=""))
    y.labels = NA
  }      
  # if y.labels is not NA continue
  if(conty){
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
      if((i == 1) & !contx) z.value = names(y.labels)[i]
      eval.js(paste("dat$",names(y.labels)[i], "=as.vector(y.labels[,i])", sep=""))
    }
  }
  #
  # xy -- assumes in this case that columns are different data vectors of row == nsmpls
  #
  contxy = TRUE
  xy.labels = as.data.frame(xy.labels)
  cngName =  grep("if ", names(xy.labels))
  names(xy.labels)[cngName] = paste("Value", cngName, sep="")    
  names(xy.labels) = gsub(pattern=" ", replacement=".",names(xy.labels))
  if( (dim(xy.labels)[1]==1) & (dim(xy.labels)[2]==1)){
    if(is.na(xy.labels[1,1])) contxy = FALSE
  }
  # dimension check
  if(((dim(xy.labels)[1] != 1) | (dim(xy.labels)[1] != 1)) & contxy){
    contxy = FALSE
    cat(paste("Warning: xy.labels does not have correct dimensions \n   number of rows should equal 1\n   Continuing with xy.labels = NA \n", sep=""))
    xy.labels = NA
  }         
  # if xy.labels is not NA continue
  if(contxy){
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
    for(i in 1:dim(xy.labels)[2]){
      if((i == 1) & !contx & !conty) z.value = names(xy.labels)[i]
      eval.js(paste("dat$",names(xy.labels)[i], "=as.vector(xy.labels[,i])", sep=""))
    }
  }
  # if all: x.labels, y.labels, and xy.labels were NA no data to display
  # set up dummy vector with blanks 
  if(!contx & !conty & !contxy){  
    eval.js(paste("dat$value=rep('',dim(dat)[2])",sep=""))
  }  

  #############################
  #############################
  #
  #  check and fill hyperlinks
  #
  #############################
  #############################

  #
  # if x specific hyperlinks
  #
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
    if(dim(x.links)[1] != 1){
      cont = FALSE
      cat(paste("Warning: x.link does not have correct dimensions \n   number of rows should equal 1\n   Different variables should be in columns\n   Continuing with x.links = NA \n", sep=""))
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
      eval.js("temp=as.vector(x.links[,i])")
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
      # put link in correct syntax into character matrix
      eval.js(paste("dat2$", names(x.links)[i], "=temp", sep=""))    
    }  
  }

  #
  # if y specific hyperlinks
  #
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
    if(dim(y.links)[1] != 1){
      cont = FALSE
      cat(paste("Warning: y.link does not have correct dimensions \n   number of rows should equal 1\n   Different variables should be in columns\n   Continuing with y.links = NA", sep=""))
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
      eval.js("temp=as.vector(y.links[,i])")
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
      # put link in correct syntax into character matrix
      eval.js(paste("dat2$", names(y.links)[i], "=temp", sep=""))    
    }  
  }

  #
  # if xy specific hyperlinks
  #
  cont = TRUE
  xy.links = as.data.frame(xy.links)
  cngName =  grep("if ", names(xy.links))
  names(y.links)[cngName] = paste("Value", cngName, sep="") 
  names(xy.links) = gsub(pattern=" ", replacement=".",names(xy.links))
  if( (dim(xy.links)[1]==1) & (dim(xy.links)[2]==1)){
    if(is.na(xy.links[1,1])) cont = FALSE
  }
  # dimension check
  if(((dim(xy.links)[1] != 1) | (dim(xy.links)[1] != 1)) & cont){
    cont = FALSE
    cat(paste("Warning: xy.links does not have correct dimensions \n   number of rows should equal 1\n   Continuing with xy.links = NA", sep=""))
    xy.links = NA
  }
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
   
    # for each column get links
    for(i in 1:length(xy.links)){
      eval.js("temp=as.vector(xy.links[,i])")
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
      # put link in correct syntax into character matrix
      eval.js(paste("dat2$", names(xy.links)[i], "=temp", sep=""))    
    }  
  }

  
  #############################
  #############################
  #
  #  check if points themselves
  #    are hyperlinks
  #
  #############################
  #############################

  
  
  # get points as Links information 
  contLinks = TRUE
  # if data frame convert to matrix
  if(class(asLinks) == "data.frame") asLinks = as.matrix(asLinks)
  # if matrix convert to vector
  if(class(asLinks) == "matrix") asLinks = as.vector(asLinks)
  # if single entry assume same for all points
  if((length(asLinks) == 1) & !is.na(asLinks[1])) asLinks = rep(asLinks, 1)
  # convert to character vector
  asLinks = as.character(asLinks)
  # check dimensions
  if((length(asLinks) != 1) & !is.na(asLinks[1])){
    cat("Warning: cannot create points as links \n     length must be equal to 1 \n")
    contLinks = FALSE
  }
  if(length(asLinks) ==1){
    if(is.na(asLinks[1])) contLinks=FALSE
  }


  ##############################
  #############################
  # check images
  #############################
  #############################


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
    if((dim(x.images)[1] != length(x.pos)) & contxi){
      contxi = FALSE
      cat(paste("Warning: x.images does not have correct dimensions \n   number of rows should be 1\n   Different variables should be in columns\n   Continuing with x.images = NA \n", sep=""))
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
      eval.js("temp=as.vector(x.images[,i])")
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
  if((dim(y.images)[1] != length(y.pos)) & contyi){
    contyi = FALSE
    cat(paste("Warning: y.images does not have correct dimensions \n   number of rows should equal 1\n   Different variables should be in columns\n   Continuing with y.images = NA", sep=""))
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
      eval.js("temp=as.vector(y.images[,i])")
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
  xy.images = as.data.frame(xy.images)
  cngName =  grep("if ", names(xy.images))
  names(xy.images)[cngName] = paste("Value", cngName, sep="")
  names(xy.images) = gsub(pattern=" ", replacement=".",names(xy.images))
  if( (dim(xy.images)[1]==1) & (dim(xy.images)[2]==1)){
    if(is.na(xy.images[1,1])) contxyi = FALSE
  }
  # dimension check
  if(((dim(xy.images)[1] != length(y.pos)) | (dim(xy.images)[1] != length(x.pos))) & contxyi){
    contxyi = FALSE
    cat(paste("Warning: xy.images does not have correct dimensions \n   number of rows should equal 1\n   Continuing with xy.images = NA", sep=""))
    xy.images = NA
  }         
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
   
    for(i in 1:dim(xy.images)[2]){
      eval.js("temp=as.vector(xy.images[,i])")
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
      eval.js(paste("dat3$", names(xy.images)[i], "=temp", sep=""))    
    }
  }
       

























  
  # create object to return
  MapObj = list()
  MapObj$xy.type = "poly"
  MapObj$nCoords = (length(x.pos)*2)
  MapObj$dat = dat
  MapObj$dat2 = dat2
  MapObj$dat3 = dat3
  MapObj$contLinks = contLinks
  MapObj$asLinks = asLinks

  
  # return object
  return(MapObj)
      
}

