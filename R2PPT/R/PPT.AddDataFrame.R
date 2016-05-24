"PPT.AddDataFrame"<-function (ppt, df = stop("Data frame must be specified"), size = c(120,110,480,320), row.names = TRUE, col.names = TRUE){

if(ppt$method=="rcom"){

    if (!comIsValidHandle(ppt$ppt)) 
        stop("Invalid handle for powerpoint application")
    if (!comIsValidHandle(ppt$pres)) 
        stop("Invalid handle for powerpoint presentation")
    if (!comIsValidHandle(ppt$Current.Slide)) 
        stop("Invalid handle for presentation slide. Make sure you add a slide before adding a data frame.")
}


    if (length(size) != 4) 
        stop("size must be a vector of length 4")
    if (!all(is.finite(size))) 
        stop("non-finite values detected in size vector")
    if (!is.data.frame(df)) 
        stop("df must be a data frame.")

    PPTtemp <- paste(tempfile(), "csv", sep = ".")
    PPTtemp <- gsub("/", "\\\\", PPTtemp)
    if (!col.names) {
        write.table(df, file = PPTtemp, row.names = row.names,col.names = col.names, sep = ",")
    }
    else {
        write.csv(df, file = PPTtemp, row.names = row.names)
    }

    #myShapes <- comGetProperty(ppt$Current.Slide, "Shapes")
    myShapes <- ppt$Current.Slide[["Shapes"]]

    #mydf <- comInvoke(myShapes, "AddOLEObject", size[1], size[2], size[3], size[4], "", PPTtemp, 0, "", 0, 0)
    mydf <- myShapes$AddOLEObject(size[1],size[2],size[3],size[4],"", PPTtemp, 0, "", 0, 0) 

    unlink(PPTtemp)
    return(invisible(ppt))
}
