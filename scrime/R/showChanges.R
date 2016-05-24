showChanges<-function(pkg="scrime"){
	filename<-file.path(find.package(pkg),"CHANGES")
	file.show(filename)
}

