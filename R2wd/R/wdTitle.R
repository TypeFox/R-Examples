wdTitle <-
function(title,
         label=substring(gsub("[, .]","_",paste("text",title,sep="_")),1,16),
         paragraph=TRUE,wdapp=.R2wd){
    wdsel<-wdapp[['Selection']]
    wdsel[["Style"]]<- -63
    newtext<-title
    newtext[newtext==""]<-"\n"
    newtext<-paste(newtext,collapse=" ")
    newtext<-gsub("\n ","\n",newtext)
    newtext<-gsub(" \n","\n",newtext)
    wdsel$TypeText(newtext)
    wdInsertBookmark(label)
    if(paragraph) {
        wdsel$TypeParagraph()
        wdsel[["Style"]]<--1}
    invisible()
}

