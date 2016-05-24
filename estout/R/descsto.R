`descsto` <- function(x=NULL,row=NULL,name=NULL,drop.row=NULL,store="default"){

# --- creating storage lists and vars
# check if descsto was run previously
prev.list <- paste(store,".dcl",sep="")
if(exists(prev.list,envir=estout:::estoutstorage) && length(eval(lapply(prev.list,as.name)[[1]],estout:::estoutstorage)) != 0){
        output <- eval(lapply(prev.list,as.name)[[1]],estout:::estoutstorage)
}
else{
        output <- list()
}

# --- function to drop rows
drop.row.fct <- function(d.r,out){
	if(! is.null(d.r) && length(out) != 0){
		for(i in 1:length(d.r)){
			n.row <- grep(d.r[[i]],out)
			out <- out[-n.row]
		}
	}
	return(output <<- out)
}

# --- grabbing parameters
if(! is.null(x)){
	input <- x
}
else{
	if(! is.null(drop.row) && length(output) != 0){
		drop.row.fct(d.r=drop.row,out=output)
	}
        assign("dcl",output,estout:::estoutstorage)
}

# --- row to be overwritten
if(is.null(row)){   
	row <- length(output)+1
}


#print(length(output)) --------- control!
# --- if is data.frame or column name
if(is(input,"data.frame")){
#print(length(input)) --------- control!
        for(i in seq(1:length(input))){
                column <- c(attributes(input)$names[[i]])
print(column)
                column <- c(column,summary(input[[i]]))
                output[[length(output)+1]] <- column
        }
}
# --- is column name (vector)
else{
        if(is.null(name)){
                return('You want to insert a single column into the output. In this case it is necessary that you provide a name for the column. Please restart the function providing a name using name=!')
                #return()
        }
        else{
                output[[row]] <- c(name,summary(input))
        }
}
# --- dropping names from the new list
if(! is.null(drop.row) && length(output) != 0){
	drop.row.fct(d.r=drop.row,out=output)
}
assign("dcl",output,estout:::estoutstorage)
}
