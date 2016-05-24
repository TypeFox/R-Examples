# an interactive function to perform a function on an object.
# Originally developed to be a clone of unix %rm -i
# July 2013:  refactor to better style, and include evout in the returned list.
#
askrm<-function(items=ls(parent.frame()),fn="rm",ask=TRUE){
       selected<-NA
	   evout<- vector('list')
       j<-1
for (thenam in c(items)){
       if(ask==TRUE){
               prmpt<-paste("Do ",fn," on ",thenam,"? [y/n]")
               readline(prompt=prmpt)->theans
               }  else theans="y"
       if(theans=="y"){
                #as.name() gets rid of quotes...
               evout[[j]] <- eval(call(fn,as.name(thenam)),envir=parent.frame(1))
               cat("performed ", paste(fn, thenam,collapse=" "),'\n')
               selected[j]<-thenam
               j<-j+1
               }
       }
 #keeping track of what happened
 outs<-list(func=fn, selected=selected,evout=evout)
 return(invisible(outs))
 }