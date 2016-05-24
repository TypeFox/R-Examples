wdType <-
function (text = "",
          italic = FALSE,
          alignment = "nothing",
          paragraph = TRUE,
          wdapp = .R2wd)
{
    wdsel <- wdapp[['Selection']]
    tt<-wdsel[['Font']]
    if (italic) tt[['Italic']]<--1 else tt[['Italic']]<-0
    tt<-wdsel[['ParagraphFormat']]
    if (alignment=="left") tt[['Alignment']]<-0
    if (alignment=="center") tt[['Alignment']]<-1
    if (alignment=="right") tt[['Alignment']]<-2
    newtext<-text
    newtext[newtext==""]<-"\n"
    newtext<-paste(newtext,collapse=" ")
    newtext<-gsub("\n ","\n",newtext)
    newtext<-gsub(" \n","\n",newtext)
    wdsel$TypeText(newtext)
    if (paragraph){
        wdsel$TypeParagraph()
    }
}

