t1<-XLwriteOpen("generic1.xls") 
### Just a meaningless matrix; function converts to data.frame and exports.
XLgeneric(t1,"s1",matrix(1:4,nrow=2))
### Now adding row names, title, etc. Note adding the title shifts the table one row down.
XLgeneric(t1,"s1",matrix(1:4,nrow=2),col1=5,addRownames=TRUE,
          title="Another Meaningless Table",rowTitle="What?",
          rowNames=c("Hey","You!"))

###... and now adding some text
XLaddText(t1,"s1","You can also add text here...",row1=10)
XLaddText(t1,"s1","...or here.",row1=11,col1=8)
XLaddText(t1,"s2",
          "Adding text to a new sheet name will create that sheet!"
          ,row1=2,col1=2)

### A more complicated example, showing how a "flattened" 3-way table might be exported:

carnames=paste(rep(c(4,6,8),each=2),"cylinders",rep(c("automatic","manual"),3))
XLgeneric(t1,'cars',ftable(mtcars$cyl,mtcars$vs,mtcars$am),
          addRownames=TRUE,rowNames=carnames,rowTitle="Engine Type",colNames=c("S","V"))

cat("Look for",paste(getwd(),"generic1.xls",sep='/'),"to see the results!\n")
