
writeCircle.1 <-function(DFs, cdat, ndat, obj){

  sep.chr = DFs$sep.chr

  # if nothing is specified besides coordinates
  # make mock black vector
  if(length(ndat) == 2){
    cdat = cbind(cdat, rep("NA", dim(cdat)[1]))
    ndat = c(ndat, "value")
    sep.chr=c(sep.chr, " ")
  }

  # dat dataframe dimensions (because hyperlinks not active in header = v1
  # loop over each pt.  coordinate  
  for(i in 1:(DFs$orgDatDim[1]) ){
    
    ctmp=paste("<area shape=\"circle\" coords=\"",cdat[i,1],",",cdat[i,2],
      ",",obj$spot.radius[i],"\" onmouseover=\"setData(\'",ndat[3],"&nbsp;&nbsp;",sep.chr[3],"&nbsp;",
      cdat[i,3],sep="")
    
    if(DFs$orgDatDim[2]>3){
      for(j in 4:(DFs$orgDatDim[2])){
        ctmp = paste(ctmp, "<br> ",ndat[j],"&amp;nbsp;&amp;nbsp;",sep.chr[j],"&amp;nbsp;",
          cdat[i,j],sep="")
      }
    }
    
    ctmp = paste(ctmp, "\')\" onMouseOut=\"clearData();\" />",sep="")
    
    # write to file
    cat(ctmp, fill=TRUE)
  } 
  
}


writeRect.1 <- function(DFs, cdat, ndat){
  
  sep.chr = DFs$sep.chr
  
  # if nothing is specified besides coordinates
  # make mock black vector
  if(length(ndat) == 4){
    cdat = cbind(cdat, rep("NA", dim(cdat)[1]))
    ndat = c(ndat, "value")
    sep.chr=c(sep.chr, " ")
  }

  # dat dataframe dimensions (because hyperlinks not active in header = v1
  # loop over each pt.  coordinate  
  for(i in 1:(DFs$orgDatDim[1]) ){
    
    ctmp=paste("<area shape=\"rect\" coords=\"",cdat[i,1],",",cdat[i,2],",",cdat[i,3],",",cdat[i,4],
      " \" onmouseover=\"setData(\'",ndat[5],"&nbsp;&nbsp;",sep.chr[5],"&nbsp;",
      cdat[i,5],sep="")
    
    if(DFs$orgDatDim[2]>5){
      for(j in 6:(DFs$orgDatDim[2])){
        ctmp = paste(ctmp, "<br> ",ndat[j],"&amp;nbsp;&amp;nbsp;",sep.chr[j],"&amp;nbsp;",
          cdat[i,j],sep="")
      }
    }
    
    ctmp = paste(ctmp, "\')\" onMouseOut=\"clearData();\" />",sep="")
    
    # write to file
    cat(ctmp, fill=TRUE)
  } 

}

writePoly.1 <- function(DFs, cdat, ndat, obj){

  nv = obj$nCoords
  sep.chr = DFs$sep.chr
    
  # if nothing is specified besides coordinates
  # make mock black vector
  if(length(ndat) == nv){
    cdat = cbind(cdat, rep("NA", dim(cdat)[1]))
    ndat = c(ndat, "value")
    sep.chr=c(sep.chr, " ")
 }

  # dat dataframe dimensions (because hyperlinks not active in header = v1
  # loop over each pt.  coordinate  
  for(i in 1:(DFs$orgDatDim[1]) ){
    
    ctmp="<area shape=\"poly\" coords=\""

    for(v in 1:(nv-1)){

     ctmp = paste(ctmp, cdat[i,v],",")
    }
    ctmp = paste(ctmp, cdat[i,nv], sep="")
    
    ctmp = paste(ctmp, " \" onmouseover=\"setData(\'",ndat[(nv+1)],"&nbsp;&nbsp;",sep.chr[(nv+1)],"&nbsp;",
      cdat[i,(nv+1)],sep="")
    
    if(DFs$orgDatDim[2]>(nv+1)){
      for(j in (nv+2):(DFs$orgDatDim[2])){
        ctmp = paste(ctmp, "<br> ",ndat[j],"&amp;nbsp;&amp;nbsp;",sep.chr[j],"&amp;nbsp;",
          cdat[i,j],sep="")
      }
    }
    
    ctmp = paste(ctmp, "\')\" onMouseOut=\"clearData();\" />",sep="")
    
    # write to file
    cat(ctmp, fill=TRUE)
  } 
  
}



writeCircle.2 <-function(DFs, cdat, ndat, obj){

  sep.chr = DFs$sep.chr
  
  # if nothing is specified besides coordinates
  # make mock black vector
  if(length(ndat) == 2){
    cdat = cbind(cdat, rep("NA", dim(cdat)[1]))
    ndat = c(ndat, "value")
    sep.chr=c(sep.chr, " ")
  }
  # if nothing static only hyperlinks
  # make mock blank vector
  if(DFs$links.st == 3){
    cdat = cbind(cdat[,1:2], rep("NA", dim(cdat)[1]), cdat[,3:dim(cdat)[2]])
    ndat = c(ndat[1:2], "value", ndat[3:length(ndat)])
    sep.chr=c(sep.chr[1:2]," ", sep.chr[3:length(sep.chr)])
  }
  

  
  # dat dataframe dimensions (because hyperlinks not active in header = v1
  # loop over each pt.  coordinate  
  for(i in 1:(DFs$orgDatDim[1]) ){    
  
    ctmp=paste("<area shape=\"circle\" coords=\"",cdat[i,1],",",cdat[i,2],
      ",",obj$spot.radius[i],"\" onmouseover=\"Tip(\'",ndat[3],"&nbsp;&nbsp;",sep.chr[3],"&nbsp;",
      cdat[i,3],sep="")

   
    if(dim(cdat)[2]>3){
       
      # static values
      if(DFs$orgDatDim[2]>3){
        #for(j in 4:(DFs$links.st-1)){
        for(j in 4:(DFs$image.st-1)){
          ctmp = paste(ctmp, "<br> ",ndat[j],"&amp;nbsp;&amp;nbsp;",sep.chr[j]," &amp;nbsp;",
            cdat[i,j],sep="")
        }
      }

      # hyperlinks 
      linkFlag = FALSE

      if(DFs$orgDat2Dim[2] > 1){
        for(j in DFs$links.st:(dim(cdat)[2])){
          if(!is.na(cdat[i,j])){
            linkFlag = TRUE
            ctmp = paste(ctmp, "<br> ",ndat[j],sep.chr[j], cdat[i,j],sep="")
          }
        }         
      }

      # images
      if(DFs$orgDat3Dim[2] > 1){
        if(DFs$orgDat2Dim[2] > 1){
          end.num = DFs$links.st - 1
        }else{
          end.num = (dim(cdat)[2])
        }
        for(j in DFs$image.st:end.num){
          if(!is.na(cdat[i,j])){
            linkFlag = TRUE
            ctmp = paste(ctmp, "<br> ",ndat[j],sep.chr[j], cdat[i,j],sep="")
          }
        }  
        
      }
      
    }else{
      linkFlag = FALSE
    }

    if(linkFlag) ctmp = paste(ctmp, "\', STICKY,true,CLICKCLOSE,true,CLOSEBTN,false, FONTFACE,\'",obj$font.type,"\', FONTCOLOR, \'",obj$font.color,"\', FONTSIZE, \'",obj$font.size,"\', BGCOLOR, \'",obj$bg.color,"\')\" ",sep="")
    if(!linkFlag) ctmp = paste(ctmp, "\', FONTFACE,\'",obj$font.type,"\', FONTCOLOR, \'",obj$font.color,"\', FONTSIZE, \'",obj$font.size,"\', BGCOLOR, \'",obj$bg.color,"\')\" ",sep="")
    # points as links
    if(obj$contLinks){
      link = obj$asLinks[i]
      if(!is.na(link)){
        ctmp = paste(ctmp, " href=\" ", link, "\" target=\"blank\" ", sep="")
      }
    }
            
    ctmp = paste(ctmp, "  />", sep="")
    # write to file
    cat(ctmp, fill=TRUE)
  }

}



writeRect.2 <-function(DFs, cdat, ndat, obj){
  
  sep.chr = DFs$sep.chr
  
  # if nothing is specified besides coordinates
  # make mock black vector
  if(length(ndat) == 4){
    cdat = cbind(cdat, rep("NA", dim(cdat)[1]))
    ndat = c(ndat, "value")
    sep.chr=c(sep.chr, " ")
  }
  # if nothing static only hyperlinks
  # make mock blank vector
  if(DFs$links.st == 5){
    cdat = cbind(cdat[,1:4], rep("NA", dim(cdat)[1]), cdat[,5:dim(cdat)[2]])
    ndat = c(ndat[1:4], "value", ndat[5:length(ndat)])
    sep.chr=c(sep.chr[1:4]," ", sep.chr[5:length(sep.chr)])
   
  }

  # dat dataframe dimensions (because hyperlinks not active in header = v1
  # loop over each pt.  coordinate  
  for(i in 1:(DFs$orgDatDim[1]) ){    
  
    ctmp=paste("<area shape=\"rect\" coords=\"",cdat[i,1],",",cdat[i,2],",",cdat[i,3],",",cdat[i,4],
      " \" onmouseover=\"Tip(\'",ndat[5],"&nbsp;&nbsp;",sep.chr[5],"&nbsp;",
      cdat[i,5],sep="")
 
    if(dim(cdat)[2]>5){

      # static values
      if(DFs$orgDatDim[2]>5){
        for(j in 6:(DFs$image.st-1)){
          ctmp = paste(ctmp, "<br> ",ndat[j],"&amp;nbsp;&amp;nbsp;",sep.chr[j],"&amp;nbsp;",
            cdat[i,j],sep="")
        }
      }
      # hyperlinks 
      linkFlag = FALSE
      if(DFs$orgDat2Dim[2] > 1){
        for(j in DFs$links.st:(dim(cdat)[2])){
          if(!is.na(cdat[i,j])){
            linkFlag = TRUE
            ctmp = paste(ctmp, "<br> ",ndat[j],sep.chr[j], cdat[i,j],sep="")
          }
        }         
      }


      # images
      if(DFs$orgDat3Dim[2] > 1){
        if(DFs$orgDat2Dim[2] > 1){
          end.num = DFs$links.st - 1
        }else{
          end.num = (dim(cdat)[2])
        }
        for(j in DFs$image.st:end.num){
          if(!is.na(cdat[i,j])){
            linkFlag = TRUE
            ctmp = paste(ctmp, "<br> ",ndat[j],sep.chr[j], cdat[i,j],sep="")
          }
        }  
        
      }
      
    }else{
      linkFlag = FALSE
    }




    
    if(linkFlag) ctmp = paste(ctmp, "\', STICKY,true,CLICKCLOSE,true,CLOSEBTN,false, FONTFACE,\'",obj$font.type,"\', FONTCOLOR, \'",obj$font.color,"\', FONTSIZE, \'",obj$font.size,"\', BGCOLOR, \'",obj$bg.color,"\')\" ",sep="")
    if(!linkFlag) ctmp = paste(ctmp, "\', FONTFACE,\'",obj$font.type,"\', FONTCOLOR, \'",obj$font.color,"\', FONTSIZE, \'",obj$font.size,"\', BGCOLOR, \'",obj$bg.color,"\')\" ",sep="")
    # points as links
    if(obj$contLinks){
      link = obj$asLinks[i]
      if(!is.na(link)){
        ctmp = paste(ctmp, " href=\" ", link, "\" target=\"blank\" ", sep="")
      }
    }
            
    ctmp = paste(ctmp, "  />", sep="")
    # write to file
    cat(ctmp, fill=TRUE)
  }

}

writePoly.2 <- function(DFs, cdat, ndat, obj){

  nv = obj$nCoords
  sep.chr = DFs$sep.chr

  
  # if nothing is specified besides coordinates
  # make mock black vector
  if(length(ndat) == nv){
    cdat = cbind(cdat, rep("NA", dim(cdat)[1]))
    ndat = c(ndat, "value")
    sep.chr=c(sep.chr, " ")    
  }
  # if nothing static only hyperlinks
  # make mock blank vector
  if(DFs$links.st == (nv+1)){
    cdat = cbind(cdat[,1:nv], rep("NA", dim(cdat)[1]), cdat[,(nv+1):dim(cdat)[2]])
    ndat = c(ndat[1:nv], "value", ndat[(nv+1):length(ndat)])
    sep.chr = c(sep.chr[1:nv], " ", sep.chr[(nv+1):length(sep.chr)])
    
  }

  # dat dataframe dimensions (because hyperlinks not active in header = v1
  # loop over each pt.  coordinate  
  for(i in 1:(DFs$orgDatDim[1]) ){    

   
    ctmp="<area shape=\"poly\" coords=\""

    for(v in 1:(nv-1)){

     ctmp = paste(ctmp, cdat[i,v],",")
    }
    ctmp = paste(ctmp, cdat[i,nv], sep="")
    
    ctmp = paste(ctmp, " \" onmouseover=\"Tip(\'",ndat[(nv+1)],"&nbsp;&nbsp;",sep.chr[(nv+1)],"&nbsp;",
      cdat[i,(nv+1)],sep="")

    if(dim(cdat)[2]>(nv+1)){

      # static values
      if(DFs$orgDatDim[2]>(nv+1)){
        for(j in (nv+2):(DFs$image.st-1)){
          ctmp = paste(ctmp, "<br> ",ndat[j],"&amp;nbsp;&amp;nbsp;",sep.chr[j],"&amp;nbsp;",
            cdat[i,j],sep="")
        }
      }
      # hyperlinks 
      linkFlag = FALSE
      if(DFs$orgDat2Dim[2] > 1){
        for(j in DFs$links.st:(dim(cdat)[2])){
          if(!is.na(cdat[i,j])){
            linkFlag = TRUE
            ctmp = paste(ctmp, "<br> ",ndat[j],sep.chr[j], cdat[i,j],sep="")
          }
        }         
      }

      # images
      if(DFs$orgDat3Dim[2] > 1){
        if(DFs$orgDat2Dim[2] > 1){
          end.num = DFs$links.st - 1
        }else{
          end.num = (dim(cdat)[2])
        }
        for(j in DFs$image.st:end.num){
          if(!is.na(cdat[i,j])){
            linkFlag = TRUE
            ctmp = paste(ctmp, "<br> ",ndat[j],sep.chr[j], cdat[i,j],sep="")
          }
        }  
        
      }
  

    }else{
      linkFlag = FALSE
    }
    if(linkFlag) ctmp = paste(ctmp, "\', STICKY,true,CLICKCLOSE,true,CLOSEBTN,false, FONTFACE,\'",obj$font.type,"\', FONTCOLOR, \'",obj$font.color,"\', FONTSIZE, \'",obj$font.size,"\', BGCOLOR, \'",obj$bg.color,"\')\" ",sep="")
    if(!linkFlag) ctmp = paste(ctmp, "\', FONTFACE,\'",obj$font.type,"\', FONTCOLOR, \'",obj$font.color,"\', FONTSIZE, \'",obj$font.size,"\', BGCOLOR, \'",obj$bg.color,"\')\" ",sep="")
    # points as links
    if(obj$contLinks){
      link = obj$asLinks[i]
      if(!is.na(link)){
        ctmp = paste(ctmp, " href=\" ", link, "\" target=\"blank\" ", sep="")
      }
    }
            
    ctmp = paste(ctmp, "  />", sep="")
    # write to file
    cat(ctmp, fill=TRUE)
  }
  
}


