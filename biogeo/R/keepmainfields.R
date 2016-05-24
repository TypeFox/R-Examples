keepmainfields <-
function (dat,ID='',Species='',x='',y='',others='') 
  {
    inputNames <- c(ID,Species,x,y)
    if(all(nchar(inputNames)==0)){
      stop("At least one of 'ID','Species','x' or'y' column names need to be specified")
    }
    reqMinNames <- c("ID","Species","x","y")
    if(any(nchar(inputNames)==0)){
      fm <- which(nchar(inputNames)==0)
    } else {
      fm <- 0
    }
    for(c in 1:length(reqMinNames)){
      if(c%in%fm){
        if(c==1){
          zv <- 1:nrow(dat)
        } else {
          zv <- rep(NA,nrow(dat))
        }
        assign(reqMinNames[c], zv)
      } else {
        assign(reqMinNames[c], dat[,inputNames[c]])
      }
    }
    if(any(nchar(others))==0){
      z1 <- data.frame(ID,Species,x,y)
    } else {
      for(o in 1:length(others)){
        zo <- dat[,others[o]]
        assign(others[o], zo)
      }
      z1 <- data.frame(ID,Species,dat[,others],x,y)
      names(z1) <- c("ID","Species",others,"x","y")
    }
    reqNames <- c("ID","Species","x","y","x_original","y_original","Correction","Modified","Exclude","Reason")
    missingNames <- reqNames[!sapply(reqNames,FUN=function(x) x%in%names(z1))]
    z2 <- data.frame(z1, ID=1:nrow(dat), Species, x=NA, y=NA, x_original = NA, y_original = NA, Correction = "........", 
                    Modified = "01-01-1900 12:01:01", Exclude = 0, Reason = "........", stringsAsFactors = F)
    z2 <- z2[,c(names(z1),missingNames)]
    return(z2)
  }
