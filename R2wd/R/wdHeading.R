wdHeading <-
function (level=1,text = "", paragraph = TRUE, wdapp = .R2wd)
{
    wdsel<-wdapp[['Selection']]
    wdsel[["Style"]]<--1-level
    newtext<-text
    newtext[newtext==""]<-"\n"
    newtext<-paste(newtext,collapse=" ")
    newtext<-gsub("\n ","\n",newtext)
    newtext<-gsub(" \n","\n",newtext)
    wdsel$TypeText(newtext)
    if (paragraph){
        wdsel$TypeParagraph()
        wdsel[["Style"]]<- -1
    }
    invisible()
}

