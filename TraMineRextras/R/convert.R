##prevwd <- setwd("C:/G/Project/Oris/eurasia/figSw")
## function for converting graphic files with the 
## ImageMagick mogrify function

## author: G. Ritschard

## usage examples
#convert.g(path="C:/G/Project/biomining/Traminer/trunk/articles/TypicalLifeCourses/Graphiques")
#convert.g(path="C:/G/Project/biotree/parcoursEric/figSw", from="pdf", to="jpg")

convert.g <- function(path = NULL, 
        fileroot= "*",
        from = "pdf",
        to = "png",
        create.path = TRUE,
        options = NULL)
    {
    
    if ((fileroot %in% c("",".")) | !is.character(fileroot)){
        warning("Invalid 'fileroot' value, no conversion has been done!")
        return(NULL)
    }
    
    v.fileroot <- strsplit(fileroot, "", fixed=TRUE)
    if (v.fileroot[length(v.fileroot)] != ".")
        {
        fileroot <- paste(fileroot,".",sep="")
        }
    
    if (!file.exists(path)){
        warning(paste(path, "not found, no conversion has been done!"))
        return(NULL)
    }

    if (fileroot != "*." & !file.exists(file <- paste(path,"\\\\",fileroot,from,sep=""))){
        warning(paste(file, "not found, no conversion has been done!"))
        return(NULL)
    }


    if(!is.null(path)){
        prevwd <- setwd(path)
    }
    else{
        prevwd <- getwd()
    }
    
    if(create.path){
        to.path <- paste(to," ", sep="")
        dir.create(to, showWarnings = FALSE)
    }
    else{
        to.path <- " "
    }
    
    if(is.null(options)){
        if(to %in% c("png","jpg")){
            options <- "-quality 100 -density 150x150"
        }
        else {
            options <- ""
        }
        
    }
    options <- paste(" ", options, sep="")
    
    fileroot <- paste(" ", fileroot, sep="")
    
    mogr.str <- paste("mogrify -path ", to.path, "-format ", to, options, fileroot, from, sep="")
    shell(mogr.str)
    setwd(prevwd)
    return(mogr.str) 
}

