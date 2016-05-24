 data(iris)
 data <- iris[sample(1:nrow(iris),size=30),]

 .HTML.file = HTMLInitFile(useGrid=TRUE,useLaTeX=FALSE)
 HTML.title("Iris dataset (again)",1,file=.HTML.file)  
 HTML(as.title("20 random observations displayed with HTMLgrid"),HR=3)
 HTML("Try to click on columns headers to sort them")
 HTMLgrid_inline(data)
 HTML(as.title("A summary of those observations displayed with HTMLgrid"),HR=3)
 
 HTMLgrid_summary(data,file=.HTML.file)
 cat("file:", .HTML.file, "is created")
 browseURL(paste("file://",.HTML.file,sep=""))

