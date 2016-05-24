require(rtf)
require(lattice)

output<-"test_rtf-package.doc"
png.res<-300

rtf<-RTF(output,width=8.5,height=11,font.size=10,omi=c(1,1,1,1))
addHeader(rtf,title="Test",subtitle="2011-08-15\n")
addPlot(rtf,plot.fun=plot,width=6,height=6,res=300, iris[,1],iris[,2])

# Try trellis plots

# single page trellis objects
addPageBreak(rtf, width=11,height=8.5,omi=c(0.5,0.5,0.5,0.5))

p <- histogram( ~ height | voice.part, data = singer, xlab="Height")
addTrellisObject(rtf,trellis.object=p,width=10,height=7.5,res=png.res)

p <- densityplot( ~ height | voice.part, data = singer, xlab = "Height")
addTrellisObject(rtf,trellis.object=p,width=10,height=7.5,res=png.res)

# multipage trellis object
p2<-xyplot(uptake ~ conc | Plant, CO2, layout = c(2,2))
addTrellisObject(rtf,trellis.object=p2,width=6,height=6,res=png.res)

addPageBreak(rtf, width=6.5,height=10,omi=c(0.5,0.5,0.5,0.5))
addParagraph(rtf,"\n\nHere's a new page with a custom size and margins.\n")

tab<-table(iris$Species,floor(iris$Sepal.Length))
names(dimnames(tab))<-c("Species","Sepal Length")
addParagraph(rtf,"\n\nTable based on a 'table' class object:\n")
addTable(rtf,tab,font.size=10,row.names=TRUE,NA.string="-",col.widths=c(1,rep(0.5,4)) )

addParagraph(rtf,"\n\nTable based on a 'data.frame' class with automatic column widths:\n")
tab<-as.data.frame(head(iris))
colnames(tab)<-gsub("\\."," ",colnames(tab))
addTable(rtf,tab,font.size=9,row.names=FALSE,NA.string="-")


tab<-head(as.data.frame(Seatbelts))
colnames(tab)<-c("car drivers killed","UK Driver Deaths","front-seat passengers killed or seriously injured","rear-seat passengers killed or seriously injured","distance driven (km)","petrol price","number of van drivers","Law in effect that month? (1/0)")
tab<-format(tab,digits=3)
addParagraph(rtf,"\n\nTable based on a 'data.frame' class with automatic column widths (exceeds page content width)\n")
addTable(rtf,tab,font.size=9,row.names=TRUE,NA.string="-")


fake.word<-function(all.upper=F,width.range=c(2,10)) {
	if(all.upper==FALSE) {
		return(paste(sample(c(LETTERS,letters),sample(width.range[1]:width.range[2])[1]),collapse=""))
	} else {
		return(paste(sample(c(LETTERS),sample(width.range[1]:width.range[2])[1]),collapse=""))
	}
}

fake.words<-function(word.max=3,all.upper=F,width.range=c(2,10)) {
	paste(replicate(sample(word.max)[1],fake.word(all.upper=all.upper,width.range=width.range)),collapse=" ")
}


# Single "word" data, first all uppercase, then mixed upper and lower case
ncols<-3
tab<-as.data.frame(matrix(replicate(ncols*5,fake.words(word.max=1,all.upper=T)),ncol=ncols))
colnames(tab)<-paste("C",c(1:ncols),sep="")
addParagraph(rtf,"\n\nAnother table based on a 'data.frame' class with automatic column widths (single-word/all upper case)\n")
addTable(rtf,tab,font.size=10,row.names=TRUE,NA.string="-")

ncols<-3
tab<-as.data.frame(matrix(replicate(ncols*5,fake.words(word.max=1,all.upper=F)),ncol=ncols))
colnames(tab)<-paste("C",c(1:ncols),sep="")
addParagraph(rtf,"\n\nAnother table based on a 'data.frame' class with automatic column widths (single-word/mixed case)\n")
addTable(rtf,tab,font.size=10,row.names=TRUE,NA.string="-")


# Multi-"word" data, first all uppercase, then mixed upper and lower case
ncols<-3
tab<-as.data.frame(matrix(replicate(ncols*5,fake.words(word.max=2,all.upper=T)),ncol=ncols))
colnames(tab)<-paste("C",c(1:ncols),sep="")
addParagraph(rtf,"\n\nAnother table based on a 'data.frame' class with automatic column widths (multi-word/all upper case)\n")
addTable(rtf,tab,font.size=10,row.names=TRUE,NA.string="-")

ncols<-4
tab<-as.data.frame(matrix(replicate(ncols*5,fake.words(word.max=3,all.upper=F)),ncol=ncols))
colnames(tab)<-paste("C",c(1:ncols),sep="")
addParagraph(rtf,"\n\nAnother table based on a 'data.frame' class with automatic column widths (multi-word/mixed case)\n")
addTable(rtf,tab,font.size=10,row.names=TRUE,NA.string="-")

# Try some UTF-8 output
addParagraph(rtf,"\n\nTry some UTF-8 output in a table\n")
addTable(rtf,data.frame(utf="\u2586"),font.size=10,row.names=FALSE,NA.string="-")

# Add a Greek substitutions
addParagraph(rtf,"\n\nTry some Greek letters\n&Alpha; &Beta; &Gamma; &Delta; &Epsilon; &Zeta; &Eta; &Theta; &Iota; &Kappa; &Lambda; &Mu; &Nu; &Xi; &Omicron; &Pi; &Rho; &Sigma; &Tau; &Upsilon; &Phi; &Chi; &Psi; &Omega; &alpha; &beta; &gamma; &delta; &epsilon; &zeta; &eta; &theta; &iota; &kappa; &lambda; &mu; &nu; &xi; &omicron; &pi; &rho; &sigmaf; &sigma; &tau; &upsilon; &phi; &chi; &psi; &omega;\n\n")


addNewLine(rtf)
addNewLine(rtf)
addSessionInfo(rtf)
done(rtf)

# Open in word processor
#system(paste("open",output))


# Test forest plot
tab<-data.frame(
	Label=c("Test1","Test2","Test3"),
	HR=c(1,2,0.45),
	Lower.CI=c(0.5,1.1,0.25),
	Upper.CI=c(2,3.5,0.9), 
	stringsAsFactors=FALSE, 
	check.names=FALSE)

# create forest plots by row
forest.plot.args<-list(xlim=c(0.1,5),width=3.0,height=0.3,cex=1,lwd=0.75,res=300)
tab$"HR Plot (log scale)"<-mapply(rtf.forest.plot,tab$HR,tab$Lower.CI,tab$Upper.CI,
			MoreArgs=forest.plot.args)

# rbind the x-scale to the table in the plot column
xscale<-rtf.forest.plot.xscale(xlim=c(0.1,5),width=3.0,height=0.3,cex=1,
			lwd=0.75,res=300)

tab<-data.frame(lapply(tab, as.character), 
			stringsAsFactors=FALSE, 
			check.names=FALSE)

tab<-rbind(tab,list("","","","",xscale))

# write the RTF output
rtf<-RTF("test.forest.plot.doc",width=8.5,height=11,font.size=10,omi=c(1,1,1,1))
addTable(rtf,tab,col.widths=c(0.75,0.75,0.75,0.75,3))
done(rtf)

# End of forest plot test







