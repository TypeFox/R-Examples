ViewHelp <-
function()
{
helpfile<-paste(find.package("dgmb"))
helpfile<-substring(helpfile,1,nchar(helpfile))
helpfile<-paste(helpfile,"/docs/dgmb-manual.pdf")
helpfile<-sub("dgmb ","dgmb", helpfile)
cat ("Viewing dgmb help...otherwise, you should write help(dgmb) in R console\n")
browseURL(helpfile)
}

