
## QQQ { not optimized, need to be changed
get.file.basename<- function(filename) {
	tempname<-strsplit( basename(filename),".RData")[[1]][1]
	basename<-paste( strsplit(tempname, "_")[[1]][1],strsplit(tempname, "_")[[1]][2],
                     strsplit(tempname, "_")[[1]][3],strsplit(tempname, "_")[[1]][4],sep="_")
	basename
}

get.project.plate.name<- function(file) {
	tempname<-strsplit( basename(file),".RData")[[1]][1]
	basename<-paste( strsplit(tempname, "_")[[1]][1],strsplit(tempname, "_")[[1]][2],
                     strsplit(tempname, "_")[[1]][3],sep="_")
	basename
}
## QQQ }

.get.all.electrodes<- function(r) {
	plate <- plateinfo(r$layout$array)
	wells <- as.matrix(sort(plate$wells))
	result <- as.vector(apply(wells,c(1,2),function(well) {.get.electrode.layout(r,well)$electrodes}))
	result 
}

.get.electrode.layout<-function(r,well) {
	plateinfo <- plateinfo(r$layout$array)
	d1 <- expand.grid(col=1:plateinfo$n.elec.c,row=1:plateinfo$n.elec.r)
	electrodes <- sort(paste(well,"_", d1[,"row"],d1[,"col"],sep=""))
	layout <- c(plateinfo$n.elec.r, plateinfo$n.elec.c)
	return(list(electrodes  = electrodes, layout = layout))
}




IGM.write.UI.to.log<-function(files=NULL,parameterList, new.file=F ){
  
  if(new.file){
    for( i in 1:length(files) ){
      cur.file=files[i]
      write(file=files[i], " "  , append=F )
    }#end of for
  }
  
  if (!is.null(files)){
    
    for( i in 1:length(files) ){
      cur.file=files[i]

      #write params
      for (j in 1:length(parameterList)){
        write(file=cur.file, " "  , append=T )
        if (length(parameterList[[j]])>1 ){
          write(file=cur.file, names(parameterList)[j]  , append=T )
          for (d in 1:length(parameterList[[j]]) ){
            write(file=cur.file, paste(names(parameterList[[j]])[d], " = ",
                                       parameterList[[j]][d] ) , append=T )
          } # end for 
        } else {
          write(file=cur.file, paste(names(parameterList)[j], " = ",
                                     parameterList[j] ) , append=T )
        }
      }# end of for legnth param list


    }#end of for lenght(files)
    
  }# end of if(!is.null(files))
} #end of .IGM.write.UI.to.log
    
