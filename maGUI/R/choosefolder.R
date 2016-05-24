choose_folder<-function(h,...)
{
	folderchoose=NULL;
	folder<-tclvalue(tkchooseDirectory(title="Select the folder"))
	folderchoose<<-folder
}
