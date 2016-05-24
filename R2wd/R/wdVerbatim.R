wdVerbatim<-
function (text = "", paragraph = TRUE, fontsize=9,fontname="Courier New",wdapp = .R2wd)
{
   wdsel <- wdapp[["Selection"]]
   wdsel$TypeParagraph()
   savestyle <- wdsel[["Style"]]
   wdsel[["Style"]] <- -67
   wdfmt<- wdsel[["ParagraphFormat"]]
   wdfmt[["LineSpacingRule"]] <- 0
   wdfont<-wdsel[["Font"]]
   wdfont[["Name"]]<-"Courier New"
   wdfont[["Size"]]<-fontsize
   newtext <- paste(text,collapse="\n")
   wdsel$TypeText(newtext)
   wdsel$TypeParagraph()
   wdsel[["Style"]] <- savestyle
}
