choose_file<-function(h,...)
{
	filechoose=NULL;
	file<-tclvalue(tkgetOpenFile(title="Select the file"))
	filechoose<<-file
}

