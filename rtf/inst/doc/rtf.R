### R code from vignette source 'rtf.Rnw'

###################################################
### code chunk number 1: rtf.Rnw:36-37
###################################################
library(rtf)


###################################################
### code chunk number 2: rtf.Rnw:43-50
###################################################
output<-"rtf_vignette.doc"   # although this is RTF, we can use the
                             # .doc extension so it opens in MS Word
rtf<-RTF(output,width=8.5,height=11,font.size=10,omi=c(1,1,1,1))

# Other rtf commands here...

done(rtf)                    # writes and closes the file


###################################################
### code chunk number 3: rtf.Rnw:61-63
###################################################
addHeader(rtf,title="Section Header",
	      subtitle="This is the subheading or section text.")


###################################################
### code chunk number 4: rtf.Rnw:68-69
###################################################
addParagraph(rtf,"This is a new self-contained paragraph.\n")


###################################################
### code chunk number 5: rtf.Rnw:75-80
###################################################
startParagraph(rtf)
addText(rtf,"This text was added with the addText command.  ")
addText(rtf,"You can add styled text too.  ",bold=TRUE,italic=TRUE)
addText(rtf,"You must end the paragraph manually.")
endParagraph(rtf)


###################################################
### code chunk number 6: rtf.Rnw:87-88
###################################################
addNewLine(rtf)


###################################################
### code chunk number 7: rtf.Rnw:99-100
###################################################
addParagraph(rtf,"&Alpha; &Beta; &Gamma; &Delta; &Epsilon;\n\n")


###################################################
### code chunk number 8: rtf.Rnw:105-106
###################################################
addParagraph(rtf,"&alpha; &beta; &gamma; &delta; &epsilon;\n\n")


###################################################
### code chunk number 9: rtf.Rnw:112-113
###################################################
addParagraph(rtf,"\\u9829\\3 \\u9829\\3 \\u9829\\3\n\n")


###################################################
### code chunk number 10: rtf.Rnw:133-134
###################################################
addParagraph(rtf,"Normal, \\b this is bold\\b0, normal.\n")


###################################################
### code chunk number 11: rtf.Rnw:138-139
###################################################
addParagraph(rtf,"Normal, {\\b\\i bold-italic}, normal.\n")


###################################################
### code chunk number 12: rtf.Rnw:158-160
###################################################
tab<-as.data.frame(head(iris)) # create a data.frame
colnames(tab)<-gsub("\\."," ",colnames(tab)) # format column names


###################################################
### code chunk number 13: tab1
###################################################
library(xtable)
print(xtable(tab), table.placement = "!htbp")


###################################################
### code chunk number 14: rtf.Rnw:168-169
###################################################
addTable(rtf,tab,font.size=9,row.names=FALSE,NA.string="-")


###################################################
### code chunk number 15: rtf.Rnw:175-177
###################################################
tab<-table(iris$Species,floor(iris$Sepal.Length))
names(dimnames(tab))<-c("Species","Sepal Length")


###################################################
### code chunk number 16: tab2
###################################################
print(xtable(tab), table.placement = "!htbp")


###################################################
### code chunk number 17: rtf.Rnw:185-187
###################################################
addTable(rtf,tab,font.size=10,row.names=TRUE,NA.string="-",
             col.widths=c(1,0.5,0.5,0.5,0.5) )


###################################################
### code chunk number 18: rtf.Rnw:197-198 (eval = FALSE)
###################################################
## addPlot(RTF.object, plot.fun=plot.fun, width=4, height=5, res=300, ...)


###################################################
### code chunk number 19: rtf.Rnw:205-206
###################################################
plot(iris[,1],iris[,2])


###################################################
### code chunk number 20: rtf.Rnw:211-212
###################################################
addPlot(rtf,plot.fun=plot,width=6,height=6,res=300, iris[,1],iris[,2])


###################################################
### code chunk number 21: rtf.Rnw:217-224
###################################################
newPlot<-function() {
	par(pty="s",cex=0.7)      # adjust plot style
	plot(iris[,1],iris[,2])
	abline(h=2.5,v=6.0,lty=2) # add some lines
}

newPlot()


###################################################
### code chunk number 22: rtf.Rnw:229-230
###################################################
addPlot(rtf,plot.fun=newPlot,width=6,height=6,res=300)


###################################################
### code chunk number 23: rtf.Rnw:238-241
###################################################
library(lattice)
p <- histogram( ~ height | voice.part, data = singer, xlab="Height")
print(p)


###################################################
### code chunk number 24: rtf.Rnw:246-247
###################################################
addPlot(rtf,plot.fun=print,width=5,height=5,res=300,p)


###################################################
### code chunk number 25: rtf.Rnw:253-255
###################################################
p2 <- densityplot( ~ height | voice.part, data = singer, xlab = "Height")
print(p2)


###################################################
### code chunk number 26: rtf.Rnw:258-259
###################################################
addTrellisObject(rtf,trellis.object=p2,width=5,height=5,res=300)


###################################################
### code chunk number 27: rtf.Rnw:263-266
###################################################
p3<-xyplot(uptake ~ conc | Plant, CO2, layout = c(2,2))
print(p3) # note this is a multipage lattice plot
          # but Sweave only shows the first plot


###################################################
### code chunk number 28: rtf.Rnw:269-270
###################################################
addTrellisObject(rtf,trellis.object=p3,width=6,height=6,res=300)


###################################################
### code chunk number 29: rtf.Rnw:278-282
###################################################
# plot
library(ggplot2)
mt <- ggplot(mtcars, aes(mpg, wt, colour = factor(cyl))) + geom_point()
print(mt)


###################################################
### code chunk number 30: rtf.Rnw:290-291
###################################################
addPlot(rtf,plot.fun=print,width=5,height=4,res=300, mt)


###################################################
### code chunk number 31: rtf.Rnw:298-299 (eval = FALSE)
###################################################
## addPng(rtf, "foo.png", width=5, height=5)


###################################################
### code chunk number 32: rtf.Rnw:306-312 (eval = FALSE)
###################################################
## addHeader(rtf,"Table of Contents")
## addTOC(rtf)
## 
## addHeader(rtf,"Section 3",TOC.level=1)
## addHeader(rtf,"Section 3A",TOC.level=2)
## addHeader(rtf,"Section 3B",TOC.level=2)


###################################################
### code chunk number 33: rtf.Rnw:319-323
###################################################
addPageBreak(rtf, width=8.5, height=11, omi=c(1,1,1,1))

addSessionInfo(rtf)
done(rtf)


