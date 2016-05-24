# $Id: grid.R 45 2008-05-23 15:50:48Z mentus $ 

"HTMLgrid_summary" <- function(
    x,
    file=NULL,
    append=TRUE,
    digits=getOption("R2HTML.format.digits"),
    nsmall = getOption("R2HTML.format.nsmall"), 
    big.mark = getOption("R2HTML.format.big.mark"), 
    big.interval = getOption("R2HTML.format.big.interval"), 
    decimal.mark = getOption("R2HTML.format.decimal.mark"),
    browse=FALSE)
{
    if (!is.data.frame(x)) stop("\n x must be a data.frame")
     # Handle file
     if (is.null(file)) file <- paste(tempfile(),"htm",sep=".")
    # ensure append mode
    cat("\n", file=file,append=append)
    whichnum= sapply(x,is.numeric)
    # numeric variables
    if (any(whichnum)){
          cat("<p>Numeric variables:</p>",file=file,append=TRUE)
          numx <- x[,whichnum,drop=FALSE]
          
          mysummary <- function(vec) {
            out<- c(
            min(vec,na.rm=TRUE),
            quantile(vec,0.25,na.rm=TRUE),
            median(vec,na.rm=TRUE),
            mean(vec,na.rm=TRUE),
            quantile(vec,0.75,na.rm=TRUE),
            max(vec,na.rm=TRUE),
            nobs=sum(!is.na(vec))) 
            names(out) <- c("Min.","First Qu.","Median","Mean","Third Qu.","Max.","N.obs")
            return(out)
          }
          out=do.call("rbind", lapply(numx,FUN=mysummary))
          out=as.data.frame(out)
          out<-cbind("Variable"=rownames(out) ,out)
           # format DF 
          out.formatted <- format(out, digits=digits, nsmall=nsmall, 
            big.mark=big.mark, big.interval=big.interval, 
            decimal.mark=decimal.mark)       
          out.formatted <- as.matrix(out.formatted)
          out.formatted[is.na(out.formatted)] <- "NA"
          out.formatted[is.nan(out.formatted)] <- "NA"
          out.formatted[(out.formatted)==" NA"] <- "NA"
          out.formatted <- as.data.frame(out.formatted)
          HTMLgrid_inline(out.formatted,
            file=file, append=TRUE,digits=digits,
            browse=FALSE,classes=c("character",rep("numeric",7)),
            showimages=FALSE)
      }
    else     cat("<p>No numeric variable.</p>",file=file,append=TRUE)
   # Factor variables
    if (any(!whichnum)){
      cat("<p>Factor variables:</p>",file=file,append=TRUE)
      facx <- x[,!whichnum,drop=FALSE]
      mysummary <- function(vec)
      {
        out <- cbind(nlevels(vec),as.data.frame(table(vec)))
        return(out)
      }
      out=do.call("rbind", lapply(facx,FUN=mysummary))
      out=as.data.frame(out)
      out<-cbind(rownames(out) ,out)
      colnames(out) <- c("Variable.level","N.levels","level","Nobs")
      HTMLgrid_inline(out,file=file,append=TRUE,
          digits=digits,browse=FALSE,includeref=FALSE,,showimages=FALSE)
    } 
    else   cat("<p>No factor variable.</p>",file=file,append=TRUE)
   if (browse) browseURL(file)     
}

### HTMLgrid_summary(x,file=NULL,browse=TRUE)

# HTMLgrid uses the widget from ActiveWidget to display a DataFrame

"HTMLgrid" <- function(
    x, 
    file = HTMLGetFile(), 
    append=TRUE ,
    includeref=FALSE,
    align="center", 
    digits=getOption("R2HTML.format.digits"),
    nsmall = getOption("R2HTML.format.nsmall"),
    big.mark = getOption("R2HTML.format.big.mark"),
    big.interval = getOption("R2HTML.format.big.interval"),
    decimal.mark = getOption("R2HTML.format.decimal.mark"),
    asDF=TRUE,
    browse=FALSE,
    classes=NULL,
    showimages=TRUE
    ){
  if (!is.data.frame(x)) stop("\n x must be a data.frame")

  # Handle file
  if (is.null(file)) file <- paste(tempfile(),"htm",sep=".")

  # Handle append
  cat("\n",file=file,append=TRUE)

  # First grid? Should we include reference to stuff in the file?
   if (includeref) HTMLgrid_references(file)

  # Handle classes
  if (is.null(classes)) classes <- sapply(x,class)
  
  # format DF 
   x.formatted <- format(x, digits=digits, nsmall=nsmall, big.mark=big.mark,
     big.interval=big.interval, decimal.mark=decimal.mark)       
   x.formatted <- as.matrix(x.formatted)
   x.formatted[is.na(x.formatted)] <- "NA"
   x.formatted[is.nan(x.formatted)] <- "NA"
   x.formatted[(x.formatted)==" NA"] <- "NA"
            
   # write it to the same dir as the file, as a text file
   targetdatafile <- paste(dirname(file),
      paste(basename(tempfile()),"txt",sep="."),sep="/")
   write.table(x.formatted,file=targetdatafile,append=FALSE,
        sep="\t",row.names=FALSE,col.names=TRUE)
   
    # write the code that writes the grid (import data)
    txt<-c("")
    if (align=="center") txt <- c(txt,"\n<center>")
    txt<-c(txt,"\n\t<script language='javascript'>")
   
   # function ReadExternalData(URLfic,classes,asDF,ID)
    txt<-c(txt,paste("ReadExternalData('",
      basename(targetdatafile),"', ",
      javascriptArray(classes), ", ",
      as.numeric(asDF),", '",
      basename(tempfile(pattern="grid",tmpdir="")), 
      "',",
      as.numeric(showimages)
       ,");",sep=""))
    
    txt<-c(txt,"\t</script>")
    if (align=="center") txt <- c(txt,"\n</center>")

    cat(paste(txt,collapse="\n"),file=file,append=TRUE)
  if (browse) browseURL(file)
  invisible(return(file))
}



"HTMLgrid_inline" <- function(
   x,
   file = HTMLGetFile(),
   append=TRUE ,
   includeref=FALSE,
   align="center",
   digits=getOption("R2HTML.format.digits"),
   nsmall = getOption("R2HTML.format.nsmall"),
   big.mark = getOption("R2HTML.format.big.mark"),
   big.interval = getOption("R2HTML.format.big.interval"),
   decimal.mark = getOption("R2HTML.format.decimal.mark"),
   asDF=TRUE,
   browse=FALSE,
   classes=sapply(x,class),
   showimages=TRUE){

  # Handle file
  if (is.null(file)) file <- paste(tempfile(),"htm",sep=".")
  # format DF 
   x.formatted <- format(x, digits=digits, nsmall=nsmall, big.mark=big.mark, big.interval=big.interval, decimal.mark=decimal.mark)       
   x.formatted <- as.matrix(x.formatted)
   x.formatted[is.na(x.formatted)] <- "NA"
   x.formatted[is.nan(x.formatted)] <- "NA"
   x.formatted[(x.formatted)==" NA"] <- "NA"

   # First grid? Should we include reference to stuff in the file?
   if (includeref) HTMLgrid_references(file)
  
    # write the code that writes the grid 
    txt<-c("")
    if (align=="center") txt <- c(txt,"\n<center>")
    txt<-c(txt,"\n\t<script language='javascript'>")

    tmpcode=javascriptArrayWrite(x,writeit=FALSE)    
    txt <- c(txt, tmpcode$code)  
    ID <- basename(tempfile(pattern="grid",tmpdir=""))
     txt <- c(txt,paste("var ",ID,"=new Active.Controls.Grid;",sep=""))
     txt <- c(txt,paste(ID,".setId('",ID,"'),",sep=""))
      ### ReadInlineData(data,columns,nrow,ncol,classes,asDF,obj)
     txt<-c(txt,paste("ReadInlineData(",
        tmpcode[1],", ",
        tmpcode[2],",",
        nrow(x),",",
        ncol(x),",",
        javascriptArray(classes),", ",
        as.numeric(asDF),", ",
        ID, ", ",
        as.numeric(showimages),
        ");",sep=""))
    txt<-c(txt,"\t</script>")
    if (align=="center") txt <- c(txt,"\n</center>")
    cat(paste(txt,collapse="\n"),file=file,append=TRUE)
  if (browse) browseURL(file)
  invisible(return(file))
}
###     HTMLgrid_inline(iris,file=NULL,browse=TRUE)




 "HTMLgrid_references" <- function(file=NULL)
 {
      txt <- paste("<link href=\"",getOption("R2HTML.grid.stuffbasepath"),"runtime/styles/xp/grid.css\" rel=\"stylesheet\" type=\"text/css\" ></link>",sep="")
      txt <- c(txt, "<link href=\"gridR2HTML.css\" rel=\"stylesheet\" type=\"text/css\" ></link>")
    txt <- c(txt,paste("\n<script src=\"",getOption("R2HTML.grid.stuffbasepath"), "runtime/lib/grid.js\"></script>",sep=""))
    txt <- c(txt, paste("\n<script src=\"",getOption("R2HTML.grid.stuffbasepath"),"gridR2HTML.js\"></script>",sep=""))
    txt <- paste(txt, collapse="\n")
    if (!is.null(file)) cat(txt,file=file,append=TRUE)
    return(txt)
 }
 


"get.environment" <- function(pos) {
	# Retrieve an environment from a position in search(), or from its name
	# if pos = -1, returns the parent environment of the calling function
	# if the environment is "package:base", We need to use '.BaseNamespaceEnv'
	# which is experimental (!) to retrieve it!
	envir <- if (pos == -1) parent.frame(2) else as.environment(pos)
	#It seems to be normal to have NULL for the base environment!
   if (is.null(envir) && (pos == length(search()) || pos == "package:base")) envir <- .BaseNamespaceEnv
	return(envir)
}

"javascriptArray"=function(x,dimx=length(dim(x)))
  {
  #   Compute the expression
  	xexp <- try(if (inherits(x, "expression")) x else NULL, silent = TRUE)
  	if (inherits(xexp, "try-error") || is.null(xexp)) {
  		xexp <- substitute(x)
  		if (is.character(xexp)) # To make sure that non conventional names will be correctly evaluated, we use backticks!
  			xexp <- parse(text=paste("`", xexp, "`", sep = ""))
  		xexp <- as.expression(xexp)
  	}
    x=eval(xexp,envir=get.environment(-1))
  
  if (dimx==2) {  
    Jarray= paste("['",do.call("paste",c(x,sep="','")),"']",sep="",collapse=",\n")
    Jarray=paste("[",Jarray,",]",sep="\n")
  }
  else   {
       Jarray= paste("[",paste("'",do.call("paste",list(x,sep="','")),"'",sep="",collapse=","),"]",sep="")
  }
  return(Jarray)
}

# write two javascripts arrays definitions
# Return their names
"javascriptArrayWrite" <- function(x,file=NULL,append=TRUE,namex=NULL,insidejavascript=TRUE,writeit=TRUE)
{
    if (!is.data.frame(x)) stop("\n Apply for a data frame...")      
    # Handle file
    if (is.null(file)) file <- paste(tempfile(),"htm",sep=".")
    # Handle append
    cat("\n", file=file,append=append)
    # handle namex
    if (is.null(namex)){
          namex=basename(tempfile("data",""))
    }
    namexcol = paste(namex,"col",sep="")
    txt <- c("\n")
    if (!insidejavascript) txt <- c(txt,"<script language='javascript'>")
    txt <- c(txt, paste("var ", namex, "=", javascriptArray(x),";",sep=""))
    txt <- c(txt, paste("var ", namexcol, "=", javascriptArray(colnames(x)),";",sep=""))
    if (!insidejavascript) txt <- c(txt,"</script>")
    code= paste(txt,collapse="\n")
    if(writeit) cat(code,"\n",file=file,append=TRUE,sep="")
    return(list(namex,namexcol,code=code)) 
}


#
#   x=iris[sample(1:nrow(iris),size=10),]
#   x[2:3,2]<-NA
#   .HTML.file=paste(tempfile(),"htm",sep=".")
#   HTMLgrid(x)
#   browseURL(.HTML.file)



# file.copy(paste("t:/Rpackages/R2HTML/inst/output/",dir("t:/Rpackages/R2HTML/inst/output"),sep=""),tempdir(),overwrite=TRUE)
# source("t:/Rpackages/R2HTML/R/R2HTMLgrid.R")


# Dont forget to gather all things in a ZIP file using DOS zip...
# > zip -r path/name.zip path/*.*
# zip -r R2HTMLstuff.zip *.*
# Then in R
# zip.unpack("t:/rpackages/R2HTML/inst/output/R2HTMLstuff.zip",dest=tempdir())

