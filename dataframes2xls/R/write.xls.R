# Copyright 2008-2009 by Guido van Steen. Consult the LICENSE file of the csv2xls R package for the terms of use. 

# This R script makes use csv2xls version 0.4.2 

handle.pythonscript <- function(string) {
	if (string == "default")
		return (system.file("python", "csv2xls.py", package = "dataframes2xls"))
	else 
		return (string)
}

system.returnvalue = 0 # initialize (nothing is wrong yet)

handle.sheetnames <- function(string , n) {
	tobe.returned = '-s '
	if (string == "default") {
		for (i in 1:n) {
			sheet = paste("Sheet", as.character(i), sep = '')
			if (i != n) {
				tobe.returned = paste(tobe.returned, sheet, ',', sep = '')
			}
			else {
				tobe.returned = paste(tobe.returned, sheet, ' ', sep = '')
			}
		}
	}
	else {
		tobe.returned = paste(tobe.returned, string, ' ', sep = '')
	}
	return (tobe.returned) 
}

handle.standard.argument <- function(arg.name, arg.string) {
	if (arg.string == "default") {
		tobe.returned = ""
	}
	else {
		arg.string = gsub(":::",",", arg.string) # Enable R users to avoid comma's in function arguments 
		tobe.returned = paste(arg.name, arg.string, " ", sep = "")
	}
	return (tobe.returned) 
}

get.header <- function(df) {
	header.elements = c(unlist(names(df)))
	header = ""
	for (element.count in 1:(length(header.elements)-1)) { 
		header = paste(header, '"', header.elements[element.count], '"', ",", sep = "") 
	}
	header = as.character(paste(header, '"', header.elements[length(header.elements)], '"', sep = "") )
} 

write.xls <- function(x, file, sh.names = "default", formats = "default", t.formats = FALSE, fnt.names = "Helvetica", fnt.metr = "default", col.widths = 48, col.names = TRUE, row.names = FALSE, to.floats = "default", python = "python", py.script = "default", sh.return=FALSE) {

	tmp.dir = tempdir()  
	if (t.formats == FALSE) { 
		t.formats = "default"
	}
	if (t.formats == TRUE) { 
		t.formats = "true"
	}
	if (fnt.names == "Helvetica") {
		fnt.names = "default"
	}
	col.widths = as.character(col.widths) 
	if (col.widths == "48") {
		col.widths = "default"
	}
	if (col.names == TRUE) { 
		col.names = "default"
	}
	if (col.names == FALSE) { 
		col.names = "false"
	}
	if (row.names == FALSE) { 
		row.names = "default"
	}
	if (row.names == TRUE) { 
		row.names = "true"
	}

	x.deparsed = (deparse(substitute(x)))
	s = unlist(strsplit(x.deparsed, ",", fixed = TRUE))
	### add by Wei Zhang, 3rd May 2011
	if(any(s==" ")) {
		s <- s[-which(s==" ")]
	}
	#######
	s = gsub("c\\(", "", s)
	s = gsub("\\)", "", s)
	s = gsub(" ", "", s)
	numberof.data.frames = length(strsplit(s,","))
	for (i in 1:(numberof.data.frames)) {
		df.tobewritten = as.data.frame(get(s[i]))
		csv.filename = paste(tmp.dir,"/csvfile",as.character(i),".csv",sep = "")
		if (row.names == "true") { 
			if (col.names == "false") {
				header = ""
			}
			else {
				header = paste ("\"\",", get.header(df.tobewritten), "\n", sep = "")
			}
			cat(header, file = csv.filename, sep = "") 
			write.table(df.tobewritten, csv.filename, append = TRUE, sep = ',', row.names = TRUE, col.names = FALSE)
		}
		if ( (row.names %in% c("false", "default") ) ) {
			if (col.names == "false") {
				header = ""
			}
			else {
				header = paste(get.header(df.tobewritten),"\n",sep = "")
			}
			cat(header, file = csv.filename, sep = "") 
			write.table(df.tobewritten, csv.filename, append = TRUE, sep = ',', row.names = FALSE, col.names = FALSE)
		}
	}

	arg0 = "python " 
	if (python != "python") { 
		arg0 = paste(python, " ", sep = "")
	} 
	test.python = system(paste(arg0, "-c pass", sep = ""))
	if (test.python > 0) {
		print(paste("Does '", python, "' exist, and is it in the path?", sep = "")) 
		system.returnvalue = 1
	}

	arg1 = handle.pythonscript(py.script)
	if (file.exists(arg1) == FALSE) {
		print(paste("Python script ", arg1 ,"does not exist.", sep = ""))
		system.returnvalue = 1
	}
	else {
		arg1 = paste(arg1, " ", sep = "")
	}

	arg2 = '-i '
	for (i in 1:(numberof.data.frames)) {
		csv.filename = paste(tmp.dir, "/csvfile", as.character(i), ".csv", sep = "")
		if (i != numberof.data.frames) {
			arg2 = paste(arg2, csv.filename, ',', sep = '')
		}
		else {
			arg2 = paste(arg2, csv.filename, ' ', sep = '')
		}
	}

	arg3 = paste("-o ", file, " ", sep = "")
	sh.names = gsub(":::", ",", sh.names)  # Enable R users to avoid comma's in function arguments 
	arg4 = handle.sheetnames(sh.names,numberof.data.frames)
	arg5 = handle.standard.argument("-f ", formats)
	arg6 = handle.standard.argument("-t ", t.formats) 
	arg7 = handle.standard.argument("-n ", fnt.names)
	arg8 = handle.standard.argument("-m ", fnt.metr)
	arg9 = handle.standard.argument("-w ", col.widths) 
	arg10 = handle.standard.argument("-r ", row.names)
	arg11 = handle.standard.argument("-x ", col.names)
	arg12 = handle.standard.argument("-c ", to.floats)
	
	csv2xls.cmd = paste(arg0,arg1,arg2,arg3,arg4,arg5,arg6,arg7,arg8,arg9,arg10,arg11,arg12,sep="")
	if (system.returnvalue == 0) { 
		system.returnvalue = system(csv2xls.cmd) 
	}
	if (system.returnvalue != 0) {
		system.returnvalue = 1
	}
	if (sh.return) {
		print(system.returnvalue)
	}
}



