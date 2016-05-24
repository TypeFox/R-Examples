wdBody <-
function (text = "", paragraph = TRUE, wdapp = .R2wd)
{
    wdsel <- wdapp[['Selection']]
    savestyle<-wdsel[["Style"]]
    wdsel[["Style"]]<--67
    newtext<-text
    newtext[newtext==""]<-"\n"
    newtext<-paste(newtext,collapse=" ")
    newtext<-gsub("\n ","\n",newtext)
    newtext<-gsub(" \n","\n",newtext)
    wdsel$TypeText(newtext)
    if (paragraph){
        wdsel$TypeParagraph()
        wdsel[["Style"]]<-savestyle}
}

