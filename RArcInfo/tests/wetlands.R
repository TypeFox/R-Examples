library(RArcInfo)

datadir<-system.file("exampleData",package="RArcInfo")
infodir<-system.file("exampleData","info",package="RArcInfo")
coveragedir<-system.file("exampleData","wetlands",package="RArcInfo")

#get.bnddata needs the last slash...
infodir<-paste(c(infodir,"/"), collapse="")

#List all the tables

covnames<-get.namesofcoverages(datadir)
tablenames<-get.tablenames(infodir)


#Display the name of the table and its filds
for(i in 1:length(tablenames[[1]]))
{
	print(c("Table: ",tablenames$TableName[i]))

	fields<-get.tablefields(infodir,tablenames$TableName[i])
	print("Fields")
	for(j in 1:length(fields))
		print(fields[[j]][1])

	#Get the data
	if(i==1)
		tabledata<-get.tabledata(infodir,tablenames$TableName[i])
	else
		tabledata<-c(tabledata, get.tabledata(infodir,tablenames$TableName[i]) )
}

#Import data fromsome tables
arc<-get.arcdata(datadir,"wetlands")
pal<-get.paldata(datadir,"wetlands")
lab<-get.labdata(datadir,"wetlands")
cnt<-get.cntdata(datadir,"wetlands")

bnd<-get.bnddata(infodir,"WETLANDS.BND")

print("Plotting all the arcs")
plotarc(arc)

print("Plotting the first ten polygons (in red) on the previous plot")
par(col="red")
plotpal(arc,pal,new=FALSE, index=1:10)

