addDefault <- function(Splot,
                       data=NA,
                       data.labels=NA,
                       links=NA,
                       links.labels=NA,
                       asLink=NA,
                       font.type="Helvetica", # 'Arial, Helvetica, sans-serif'   
                       font.color="black",  # name, or #------
                       font.size="12",      # can specify type px,pt,em, etc.
                       bg.color='#D6E3F6',  # name or #------
                       returnVl=TRUE,
                       saveFlag=FALSE,
                       saveName="Splot.RData"){

  data = as.character(data)
  data.labels = as.character(data.labels)
  links = as.character(links)
  links.labels = as.character(links.labels)
  asLink = as.character(asLink)

  if(is.na(data.labels[1])) data.labels = paste("Value", 1:length(data), sep="")
  if(is.na(links.labels[1])) links.labels = paste("Value", 1:length(links), sep="")

  for(l in 1:length(links)){
    temp = links[l]
    if(is.na(temp)){
      links[l] = NA
    }else{
      
      linkstmp = strsplit(temp, split=",")[[1]]
      new.t = " "
      for(ln in 1:length(linkstmp)){
        new.t = paste(new.t, paste("<a href=\\'", gsub(linkstmp[ln], pattern=" ", replacement=""), "\\'> ", paste(links.labels[l],ln, sep="."), " </a>", sep=""), sep=",")
      }
      new.t = gsub(new.t, pattern=" ,", replacement="")
      links[l] = new.t
      
    }
  }


  
  Default.Obj = list(
    data=data,
    data.labels=data.labels,
    links=links,
    links.labels=links.labels,
    asLink=asLink, 
    font.type=font.type,
    font.color=font.color,
    font.size=font.size,
    bg.color=bg.color
    )
  
  Splot$Default = Default.Obj


  # save and return
  if(saveFlag) save(Splot, file=saveName, compress=TRUE)
  if(returnVl) return(Splot)


  
}
  
                       
                   
