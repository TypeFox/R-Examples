writeDefault1 <- function(Splot){

  data = Splot$Default$data
  data.labels = Splot$Default$data.labels
  links = Splot$Default$links
  links.labels = Splot$Default$links$labels
  asLink = Splot$Default$asLink

  ctmp = "<area shape=\"default\" onmouseover=\"setData(\'"
  cont = TRUE
  if(length(data) == 1){
    if(is.na(data[1])) cont=FALSE
  }
  if(cont){

    ctmp = paste(ctmp, data.labels[1], "&nbsp;&nbsp;:&nbsp;", data[1], sep="")

    if(length(data) > 1){
      for(i in 2:length(data)){
        ctmp = paste(ctmp,"<br>", data.labels[i], "&nbsp;&nbsp;:&nbsp;", data[i], sep="")
      }
    }
  }
  ctmp = paste(ctmp, "\')\" onMouseOut=\"clearData();\" />",sep="")

  cat(ctmp, fill=TRUE)
} 





writeDefault2 <- function(Splot){

  data = Splot$Default$data
  data.labels = Splot$Default$data.labels
  links = Splot$Default$links
  links.labels = Splot$Default$links.labels
  asLink = Splot$Default$asLink

  ctmp = "<area shape=\"default\" onmouseover=\"Tip(\'"
  
  cont = TRUE
  if(length(data) == 1){
    if(is.na(data[1])) cont=FALSE
  }
  if(cont){

    ctmp = paste(ctmp, data.labels[1], "&nbsp;&nbsp;:&nbsp;", data[1], sep="")

    if(length(data) > 1){
      for(i in 2:length(data)){
        ctmp = paste(ctmp,"<br>", data.labels[i], "&nbsp;&nbsp;:&nbsp;", data[i], sep="")
      }
    }
  }

  cont = TRUE
  linkFlag = FALSE
  if(length(links) == 1){
    if(is.na(links[1])) cont=FALSE
  }
  if(cont){
    for(i in 1:length(links)){
      if(!is.na(links[i])){
        linkFlag = TRUE
        ctmp = paste(ctmp, "<br> ",links.labels[i],":", links[i],sep="")
      } 
    }
  }


  if(linkFlag) ctmp = paste(ctmp, "\', STICKY,true,CLICKCLOSE,true,CLOSEBTN,false, FONTFACE,\'",Splot$Default$font.type,"\', FONTCOLOR, \'",Splot$Default$font.color,"\', FONTSIZE, \'",Splot$Default$font.size,"\', BGCOLOR, \'",Splot$Default$bg.color,"\')\" ",sep="")
  if(!linkFlag) ctmp = paste(ctmp, "\', FONTFACE,\'",Splot$Default$font.type,"\', FONTCOLOR, \'",Splot$Default$font.color,"\', FONTSIZE, \'",Splot$Default$font.size,"\', BGCOLOR, \'",Splot$Default$bg.color,"\')\" ",sep="")

  if(!is.na(asLink)) ctmp = paste(ctmp, " href=\" ", asLink, "\" target=\"blank\" ", sep="")

  ctmp = paste(ctmp, "  />", sep="")
  cat(ctmp, fill=TRUE)
  
}
