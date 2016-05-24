
### Contrived example looking at, e.g., the distribution of A-K-Q card counts
### in two partners' Bridge hands


hand1=rhyper(1000,12,40,13)
hand2=rhyper(1000,12-hand1,27+hand1,13)
handNames=c("0-1",2:4,"5 or more")


### The problem is ridiculously symmetric, so I de-symmetrize the presentation slightly:

book3<-XLwriteOpen("hands.xls") 
XLtwoWay(book3,"PartnersAKQcounts",cut(hand1,c(0,2:6,14)-0.5),cut(hand2,c(0,2:5,14)-0.5),
         rowTitle="Hand 1 vs. Hand 2",rowNames=c(handNames[-5],5,"6 or more","Total"),
         colNames=c(handNames,"Total"),header=TRUE)

## Same table, but percents now condition on columns rather than rows, 
## counts/pct header row removed - but a title added:
XLtwoWay(book3,"PartnersAKQcounts",cut(hand1,c(0,2:6,14)-0.5),cut(hand2,c(0,2:5,14)-0.5),
         rowTitle="Hand 1 vs. Hand 2",rowNames=c(handNames[-5],5,"6 or more","Total"),
         colNames=c(handNames,"Total"),header=FALSE,row1=12,sumby=2,
         title="Now Percents are Summed by Column:")

cat("Look for",paste(getwd(),"hands.xls",sep='/'),"to see the results!\n")
