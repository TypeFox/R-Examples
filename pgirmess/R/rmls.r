rmls<-function(){
mylist<-select.list(ls(envir = parent.frame()),multiple=TRUE)

 if (length(mylist)!=0) {
    OK<-winDialog("yesno", "Do you want to remove these objects?")
    if (OK=="YES") rm(list=mylist,envir = parent.frame())
    }
}
