wdInsertFootnote<-
function(text="",reference="",wdapp=.R2wd)
{
    wdsel<-wdapp[['Selection']]
    wdsel$TypeText(" ")
    wdsel$MoveLeft()
    tt<-wdsel[['Footnotes']]$Add(wdsel[['Range']], reference,text)
    wdsel$MoveRight()
    wdsel$MoveRight()
    invisible(tt)
}
