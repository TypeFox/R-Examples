`desctab` <- function(filename=NULL,caption=NULL,label=NULL,csv=FALSE,dcolumn=NULL,booktabs=FALSE,store="default"){
###
# --- variable definition
header <- c("Min.","1st Qu.","Median","Mean","3rd Qu.","Max.","Missing Values")
prev.list <- paste(store,".dcl",sep="")
input_matrix <- eval(lapply(prev.list,as.name)[[1]],estout:::estoutstorage)

# --- dcolumn or center?
if(is.null(dcolumn)){
	dcolumn <- "c"
}

# --- booktabs == TRUE
if(booktabs == TRUE){
	toprule <- "\\toprule\n"
	bottomrule <- "\\bottomrule\n"
	midrule <- "\\midrule\n"
}
else{
	toprule <- "\\hline\\hline\n"
	midrule <- "\\hline\n"
	bottomrule <- "\\hline\\hline\n"
}

# --- csv output
if(csv == TRUE){
	if(! is.null(filename)){
		filename <- paste(filename,".csv",sep="")
		sink(filename)
	}
        if(! is.null(caption)){cat(caption,"\n",sep="")}  
        cat(paste(",",header,collapse=""),"\n",sep="")
        for(i in seq(1:length(input_matrix))){
                write(paste(input_matrix[[i]],collapse=","),file=filename,append=TRUE)
        }
	if(! is.null(filename)){
		sink()
	}
}
# --- Standard Output / TeX
else{
# --- writing body
	if(! is.null(filename)){
		filename <- paste(filename,".tex",sep="")
		sink(filename)
	}
        cat("\\begin{table}[htbp]\n\\centering\n\\begin{tabular}{l*{7}{",dcolumn,"}}\n",toprule,sep="")
        cat((paste("&\t\t",header,collapse="")),"\\\\\n",midrule,sep="")
        for(i in seq(1:length(input_matrix))){
                if(length(input_matrix[[i]]) == 7){
                       input_matrix[[i]] <- append(input_matrix[[i]],"0")
                }
                cat(paste(input_matrix[[i]],collapse="\t\t &"),"\\\\\n",sep="")
} 
        cat(bottomrule,"\\end{tabular}\n")
        if(! is.null(caption)){cat("\\caption{",caption,"}\n",sep="")}  
	if(! is.null(label)){cat("\\label{tab:",label,"}\n",sep="")}
        cat("\\end{table}\n",sep="")
	if(! is.null(filename)){
		sink()
	}
        
        }
# --- function end
}
