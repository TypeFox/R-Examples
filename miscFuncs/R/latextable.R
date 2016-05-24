##' latextable function
##'
##' A very useful function to create a LaTeX table from a matrix. Rounds numeric entries and also replaces
##' small numbers with standard index form equivalents.
##'
##' To get a backslash to appear, use a double backslash
##'
##' Just copy and paste the results into your LaTeX document.
##'
##' @param x a matrix, or object that can be coerced to a matrix. x can include mixed character and numeric entries. 
##' @param digits see help file for format
##' @param scientific see help file for format
##' @param colnames optional column names set to NULL (default) to automatically use column names of x. NOTE! if rownames is not NULL present, colnames must include an entry for the rownames i.e. it should be a vector of length the number of columns of x plus 1.
##' @param rownames optional row names set to NULL (default) to automatically use row names of x
##' @param caption optional caption, not normally used
##' @param narep string giving replacement for NA entries in the matrix
##' @param laststr string to write at end, eg note the double backslash!!
##' @param ... additional arguments passed to format
##' @return prints the LaTeX table to screen, so it can be copied into reports
##' @export

latextable <- function(x,digits=3,scientific=-3,colnames=NULL,rownames=NULL,caption=NULL,narep=" ",laststr="",...){

    x <- as.matrix(x)

    if(is.null(colnames)){
        colnames <- c("",colnames(x))
    }
    if(is.null(rownames)){
        rownames <- rownames(x)
    }
    


    form <- function(x,...){
        if(is.character(x)){
            x <- as.numeric(x)
        }
        xtxt <- format(x,digits=digits,scientific=scientific,...)
        if(length(grep("e",xtxt))>0){
            spl <- unlist(strsplit(xtxt,"e"))
            xtxt <- paste(spl[1],"$\\times10^{",as.character(as.numeric(spl[2])),"}$",sep="")
        }    
        return(xtxt)
    }

	d <- dim(x)
	write("","")
	write("\\begin{table}[htbp]","")
	write("    \\centering","")
	if(!is.null(caption)){write(paste("    \\caption{",caption,"}",sep=""),"")}
	cs <- "    \\begin{tabular}{"
	cn <- ""
	times <- d[2]
	if(!is.null(rownames)){times <- times + 1}
	for (i in 1:times){
		cs <- paste(cs,"c",sep="")
		if (i<times){cn <- paste(cn,colnames[i]," & ",sep="")}
	}
	cs <- paste(cs,"}",sep="")
	cn <- paste(cn,colnames[times]," \\\\ \\hline",sep="")
	write(cs,"")
	if(!is.null(colnames)){
		if (!any(length(colnames)==c(d[2],d[2]+1))){stop("Incorrect number of column names")}
		write(paste("        ",cn),"")
	}
	for (i in 1:d[1]){
		if (!is.null(rownames)){
			if (is.na(x[i,1])){
				towrite <- paste("        ",rownames[i]," & ",narep," & ",sep="")
			}
			else if (is.numeric(x[i,1])){
				towrite <- paste("        ",rownames[i]," & ",form(x[i,1])," & ",sep="")
			}
			else{
				towrite <- paste("        ",rownames[i]," & ",x[i,1]," & ",sep="")
			}
		}
		else{
			if (is.na(x[i,1])){
				towrite <- paste("        ",narep," & ",sep="")
			}
			else if (is.numeric(x[i,1])){
				towrite <- paste("        ",form(x[i,1])," & ",sep="")
			}
			else{
				towrite <- paste("        ",x[i,1]," & ",sep="")
			}
		}
		if (d[2] > 2){			
			for (j in 2:(d[2]-1)){
				if (is.na(x[i,j])){
					towrite <- paste(towrite,narep," & ",sep="")
				}
				else if (!is.na(as.numeric(x[i,j]))){
					towrite <- paste(towrite,form(x[i,j])," & ",sep="")
				} 
				else{
					towrite <- paste(towrite,x[i,j]," & ",sep="")
				} 
			}
		}	
		if (is.na(x[i,d[2]])){
			towrite <- paste(towrite,narep," \\\\",sep="")
		}
		else if (!is.na(as.numeric(x[i,d[2]]))){
			towrite <- paste(towrite,form(x[i,d[2]])," \\\\",sep="")
		}
		else{
			towrite <- paste(towrite,x[i,d[2]]," \\\\",sep="")
		}
		if (i==d[1] & laststr!=""){
			towrite <- paste(towrite," ", laststr,sep="")
		}
		write(towrite,"")
	}
	write("    \\end{tabular}","")	
	write("\\end{table}","")	
	write("","")
}
