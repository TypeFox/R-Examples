

makeCharacter <- function(DF){
  
  dat = DF$dat
  dat2 = DF$dat2
  dat3 = DF$dat3
  

  # update dat into character array to make writing more efficient
  cdat=array(" ",dim=dim(dat))
  ndat=rep(" ",dim(dat)[2])
  for(j in 1:(dim(dat)[2])){
    cdat[,j]=as.character(dat[,j])
    ndat[j]=names(dat)[j]
  }
  ndat = gsub(pattern=".", replacement=" ",ndat, fixed=TRUE)
   
  cdat2=array(" ",dim=dim(dat2))
  ndat2=rep(" ",dim(dat2)[2])
  for(j in 1:(dim(dat2)[2])){
    cdat2[,j]=as.character(dat2[,j])
    ndat2[j]=names(dat2)[j]
  }
  ndat2 = gsub(pattern=".", replacement=" ",ndat2, fixed=TRUE)
   
  cdat3=array(" ",dim=dim(dat3))
  ndat3=rep(" ",dim(dat3)[2])
  for(j in 1:(dim(dat3)[2])){
    cdat3[,j]=as.character(dat3[,j])
    ndat3[j]=names(dat3)[j]
  }
  ndat3 = gsub(pattern=".", replacement=" ",ndat3, fixed=TRUE)

  imageIn=FALSE
  # combined information data frame and image data frame
  if(dim(cdat3)[2] > 1){

    if(dim(cdat)[1] == 1){
      cdat = matrix(c(cdat, cdat3[,2:dim(cdat3)[2]]), nrow=1)
      ndat = c(ndat, ndat3[2:length(ndat3)])
    }else{
      cdat = cbind(cdat, cdat3[,2:dim(cdat3)[2]])
      ndat = c(ndat, ndat3[2:length(ndat3)])
    }
    imageIn = TRUE
  }
 
  
  hypIn = FALSE
  # combined information data frame and hyper link data frame
  if(dim(cdat2)[2] > 1){

    if(dim(cdat)[1] == 1){
      cdat = matrix(c(cdat, cdat2[,2:dim(cdat2)[2]]), nrow=1)
      ndat = c(ndat, ndat2[2:length(ndat2)])
    }else{
      cdat = cbind(cdat, cdat2[,2:dim(cdat2)[2]])
      ndat = c(ndat, ndat2[2:length(ndat2)])
    }
    hypIn = TRUE
  }

  if(dim(cdat3)[2] > 1){
    image.st = (dim(dat)[2]) + 1
  }else{
    image.st = (dim(dat)[2]) + 1
  }

  
  if(dim(cdat3)[2] > 1){
    links.st = (dim(dat)[2] + dim(dat3)[2])
  }else{
    links.st = (dim(dat)[2]) + 1
  }


  # combine sep.chr
  sep.chr = DF$sep.chr
  xy.type = DF$xy.type
  if((xy.type=="points") | (xy.type=="circle")) new.sep = c(NA,NA)
  if( (xy.type=="image.midpoints") | (xy.type=="image.boundaries") |  (xy.type=="image.box"))new.sep = c(NA,NA)
  if( xy.type=="rect")new.sep = c(NA,NA, NA, NA)
  if( xy.type=="poly" ) new.sep = rep(NA,DF$nCoords)
  for(i in 1:9){
    vec = sep.chr[i]
    if(class(vec[[1]]) != "logical") new.sep = c(new.sep, vec[[1]])
    
  }

  
  cDF = list()
  cDF$cdat = cdat
  cDF$ndat = ndat
  cDF$orgDatDim = dim(dat)
  cDF$orgDat2Dim = dim(dat2)
  cDF$orgDat3Dim = dim(dat3)
  cDF$links.st = links.st
  cDF$image.st = image.st
  cDF$imageIn = imageIn
  cDF$hyperIn = hypIn
  cDF$sep.chr = new.sep
  
  return(cDF)
  
}  
