internal_dictionary <- function(table,file.res){


  write("\n\n",file.res)
      
  vari.labels<-unlist(sapply(table, function(x) ifelse(is.null(attr(x,"label")),"",attr(x,"label")))) 
  var.names<-names(table)
  missings<-sapply(table, function(x) attr(x,"miss.values"))    
  value.labels.list<-lapply(table,function(x) attr(x,"value.labels"))


  ## variables ##

  vari.labels.sint<-
      c("vari labels",
      paste("\t",var.names,"\t\"",vari.labels,"\"",
      c(rep("",length(vari.labels)-1),"."),sep=""))

  write(paste("\n\n\n********** DICTIONARY: ****************.",sep=""),file.res,append=TRUE)
  write(" ",file.res,append=TRUE)
  write(" ",file.res,append=TRUE)
  write("*** VARI LABELS ***.",file.res,append=TRUE)
  write(vari.labels.sint,file.res,append=TRUE)

  ## values ##

  value.labels.sint=NULL

  if (length(value.labels.list)>0){

    for (j in 1:length(value.labels.list)){
    
      value.labels<-value.labels.list[[j]]

      if (!is.null(value.labels)){
		
		    values<-as.character(value.labels)

		    if (length(grep("\\)",values))==0){     
		
		      if (length(grep("[0-9]-",values))>0) values<-paste("'",values,"'",sep="")
        
          value.labels<-data.frame(labels=names(value.labels),values=values)
 
          aux<-c(paste("value labels",var.names[j]),
              paste("\t",value.labels[,2],"\t\"",value.labels[,1],"\"",
              c(rep("",nrow(value.labels)-1),"."),sep=""))

          value.labels.sint=c(value.labels.sint," ",aux)

          if (length(value.labels.list[[j]])>1){
            table.value.labels<-cbind(c(var.names[j],rep("",nrow(value.labels)-1)),c(vari.labels[j],rep("",nrow(value.labels)-1)),value.labels)
          } else table.value.labels<-c(var.names[j],vari.labels[j],value.labels)

        }

      }
      
    }
    
    write(" ",file.res,append=TRUE)
    write("*** VALUE LABELS ***.",file.res,append=TRUE)
    write(value.labels.sint,file.res,append=TRUE)

  }
	
  ## missings ##

  if (length(missings)>0){  
    missing.sint=NULL
    for (j in 1:length(missings)){
      if (!is.null(missings[[j]])){
        aux.miss<-paste("(",paste(missings[[j]],collapse=" "),")",sep="")
        aux.miss<-gsub(":"," thru ",aux.miss)
        missing.sint<-c(missing.sint,paste("MISSING VALUES ",names(missings)[j]," ",aux.miss,".",sep=""))
      }
    }
    write(" ",file.res,append=TRUE)
    write("*** MISSING VALUES ***.",file.res,append=TRUE)
    write(" ",file.res,append=TRUE)
    write(missing.sint,file.res,append=TRUE)
  }

}
