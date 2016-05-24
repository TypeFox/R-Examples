book1<-XLwriteOpen("chick1.xls") 
XLoneWay(book1,"Diets",ChickWeight$Diet)
### Now in separate columns, and with a title - note it shifts the table down.
XLoneWay(book1,"Diets",ChickWeight$Diet,combine=FALSE,row1=10,
         rowTitle="Diet",title="Counts by Diet:")
cat("Look for",paste(getwd(),"chick1.xls",sep='/'),"to see the results!\n")
