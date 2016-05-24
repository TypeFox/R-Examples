book2<-XLwriteOpen("chick2.xls") 
## Plain-vanilla
XLunivariate(book2,"weightByDiet",ChickWeight$weight,ChickWeight$Diet,
             title="Mean Weights by Diet",rowTitle="Diet")

## Replace mean/SD with median/range, put results beside previous
XLunivariate(book2,"weightByDiet",ChickWeight$weight,ChickWeight$Diet,
             title="Median Weights by Diet",rowTitle="Diet",col1=8,
             fun1=list(fun=roundmedian,name="Median"),fun2=list(fun=rangeString,name="range"))

### You can also do only one statistic... by "killing" one of the functions
XLunivariate(book2,"weightByAge",ChickWeight$weight,ChickWeight$Time,
             title="Mean Weights by Age",rowTitle="Age (Days)",seps=rep("",3),
             fun2=list(fun=emptee,name=""))
cat("Look for",paste(getwd(),"chick2.xls",sep='/'),"to see the results!\n")
