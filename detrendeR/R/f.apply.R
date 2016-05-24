f.apply <-
function ( FUN=NULL, simplify = TRUE, ...  ) {
#arrFiles<-as.matrix(sort(choose.files())) #tk_choose.files                     #select n files to be open
arrFiles<-as.matrix(sort(tk_choose.files()))
iLast<-dim(arrFiles)[1]                                                         #to know the nº of files selected == nrow(arrFiles)
if (iLast<1) return(
invisible(tk_messageBox(type="ok", "You should select at least one file!", caption="Problems"))   #to stop if no files were selected
#cat("\nYou should select at least one file!")
)             

sapply(X=arrFiles, FUN=FUN, simplify = simplify, ...)->n.files
if(length(n.files)!=iLast) {
tk_messageBox(type="ok", "There is some kind of problems!", caption="Problems")
#cat("There is some kind of problems!")
return(n.files)}
}

