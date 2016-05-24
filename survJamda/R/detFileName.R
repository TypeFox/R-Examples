detFileName <-
function (file.name)
{	
	sub = substr(file.name,1,3)	
	if(sub == "gse") 
		file.name = toupper(file.name)
	else
		file.name = file.name
	return (file.name) 
}

