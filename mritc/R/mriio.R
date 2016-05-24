readMRI <- function(file, dim=NULL, format) {
    if(format == "rawb.gz"){
        if(! is.null(dim)){
            if(!is.vector(dim) || length(dim) != 3)
                stop("The dimensions of the file are wrong.")
        }
        f <- gzfile(file, open="rb")
        on.exit(close(f))
        data <- readBin(f,"integer", prod(dim), size = 1, signed = FALSE)
        array(data, dim)
    }
    else if(format == "analyze"){
        readANALYZE(file)@.Data
    }
    else if(format == "nifti"){
        readNIfTI(file)@.Data
    }
    else
        stop("The 'format' is not supported.")
}

writeMRI <- function(data, file, header=NULL, format) {

    objClass <- class(data)
    if(!((objClass=="array" &&
         (length(dim(data)) ==3 || length(dim(data)) == 4 && dim(data)[4] == 1)) ||
         objClass %in% c("nifti", "anlz")))
        stop("'data' has to be an array of dimension 3 or dimension 4 with the
             forth dimension equal to 1, an object of class 'nifti', or
             an object of class 'anlz'.")
    if(format == "rawb.gz"){
        f <- gzfile(file, open="wb")
        writeBin(as.vector(ifelse(objClass=="array", data, data@.Data)), f, size=1)
        on.exit(close(f))
    }
    else if(format == "analyze"){
        switch(objClass,
               array=writeANALYZE(as(data, "anlz"), file),
               anlz=writeANALYZE(data, file),
               nifti=writeANALYZE(as(data, "anlz"), file))
    }
    else if(format == "nifti"){
        switch(objClass,
               array=writeNIfTI(as(data, "nifti"), file),
               anlz=writeNIfTI(as(data, "nifti"), file),
               nifti=writeNIfTI(data, file))
    }
    else
        stop("The 'format' is not supported.")   
       
}
