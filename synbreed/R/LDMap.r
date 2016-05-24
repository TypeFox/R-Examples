LDMap <- function(LDmat,gpData,chr=NULL,file=NULL,fileFormat="pdf",onefile=TRUE,...){

    # catch (possible) errors
    if(class(LDmat)!="LDmat") stop("'LDmat' must be of class 'LDmat'")
    if(is.null(chr)) lg <- (1:length(LDmat$LD))[!as.logical(lapply(LDmat$LD,is.null))]
    else lg <- chr
    if (class(gpData) == "gpData"){
      pos <- gpData$map$pos
      names(pos) <- rownames(gpData$map)
    } else { stop("gpData has to be of class gpData!")}
    
     # use LD from input arguement
    ret <- LDmat

    if(!is.null(file) & onefile & fileFormat == "pdf"){
      if(substr(file, nchar(file)-nchar(fileFormat)+1, nchar(file)) != fileFormat | nchar(file) < 5)
        file <- paste(file, ".", fileFormat, sep="")
        pdf(file, onefile=onefile)
    }
 
    for (i in lg){    
      if(!is.null(file)&(fileFormat != "pdf" | fileFormat == "pdf" & !onefile)){
        if(substr(file, nchar(file)-nchar(fileFormat)+1, nchar(file)) != fileFormat | nchar(file) < 5){
          if(length(lg) <2)  
            fileName <- paste(file, ".", fileFormat, sep="") 
          else 
            fileName <- paste(file, "_chr", i, ".", fileFormat, sep="")
        } else {
          if(length(lg)>1)
            fileName <- paste(substr(file, 1, nchar(file)-nchar(fileFormat)-1), "_chr", i, ".", fileFormat, sep="")
          else 
            fileName <- file
        }
        if(fileFormat == "pdf") pdf(fileName)
        else if (fileFormat == "png") png(fileName)
        else stop("not supported file format choosen!")
      }
        color = c("#7F0000","#B30000","#D7301F","#EF6548","#FC8D59","#FDBB84","#FDD49E","#FEE8C8","#FFF7EC")
        #   using function LDheatmap
        MapUnit <- ifelse(gpData$info$map.unit=="cM","genetics","physical")
        LDheatmap(LDmat$LD[[i]], LDmeasure="r", color=color, genetic.distances=pos[rownames(LDmat$LD[[i]])],distances=MapUnit,   
                  geneMapLabelY=0.12, geneMapLabelX=0.35, title=paste("Pairwise LD on chromosome", i),...)
        if(!is.null(file)&(fileFormat != "pdf" | fileFormat == "pdf" & !onefile)){
          dev.off() 
        } else if(is.null(file) & length(lg)>1) readline()
     }
     if(!is.null(file) & onefile & fileFormat == "pdf") dev.off()
  
}
