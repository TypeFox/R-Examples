#' @title Reads zTree output files
#' @description Extracts variables from a zTree output file.
#' @param object a zTree file or a list of files.
#' @param tables the tables of intrest.  
#' @return A list of dataframes, one for each table
#' @examples
#' \dontrun{url <- 
#' zTables <- read.zTree( "131126_0009.xls" , "contracts" )
#' zTables <- read.zTree( c("131126_0009.xls",
#' "131126_0010.xls"), c("globals","subjects", "contracts" ))
#' }
#' 
#' @export
read.zTree <- function(object, tables=c("globals","subjects")) {
  splittable <- function(objname,tables=c("globals","subjects")) {
    getTable <- function(start, stop) {
      if (!is.na(stop) && !is.na(start)) {
        names<-aa2[[start]][-3]
        names[1]<-"Date"
        names[2]<-"Treatment"
        tab<-matrix(nrow=stop-start-1,ncol=length(names))
        colnames(tab)<-names
        for( i in 1:(stop-start-1)) {
          tab[i,] <- aa2[[start+i]][-3]
        }
        for (n in names(tab)) {
          if (is.factor(tab[[n]])) {
            tab[[n]][tab[[n]]=="-"] <- NA
            mm<-mean(as.numeric( levels(tab[[n]])[tab[[n]]]),na.rm=TRUE)
            if (!is.na(mm)) {
              tab[[n]]<-as.numeric( levels( tab[[n]])[tab[[n]]] )
            }
          }
        }
        tab
      }
    }
    
    getTables <- function(name) {
      tab<-NULL
      for (i in which ((splitname==name))) {
        new<-getTable(splitpoints[i],splitpoints[i+1])
        if (length(new)>0) {
          if (is.null(tab)) {
            tab<-new
          } else {
            tab <- merge(tab,new,all=TRUE)
          }
        }
      }
      tab
    }
    cat("reading ",objname,"...\n")
    Tfile<-file(objname,"r")
    aa<-readLines(Tfile)
    close(Tfile)
    aa2<-strsplit(aa,"\t")
    splitpoints<-array()
    splitname<-array()
    splittreat<-array()
    table(splitname)
    splitcols<-array()
    last<-0
    for (i in 1:length(aa2)) {
      if (last==0 || (aa2[[i]][3] != aa2[[i-1]][3])) {
        last<-last+1
        splitpoints[last]<-i
        splitname[last]<-aa2[[i]][3]
        splittreat[last]<-aa2[[i]][2]
        splitcols[last]<-length(aa2[[i]])
      }
      splitpoints[last+1]<-i+1
    }
    
    result<-list()
    do <- intersect(splitname,tables)
    miss <- setdiff(splitname,tables)
                                        #if (length(miss)>0)
    cat ("Skipping:",miss,"\n")
    for (name in do) {
      cat ("Doing:",name,"\n")
      aTable<-getTables(name)
      if (!is.null(aTable)) result[[name]]<-aTable
    }
    result
  }
  
  
  z<-splittable(object[1],tables)
  for (name in object[-1]) {
    a=splittable(name,tables)
    for(t in tables) {
      z[[t]]=merge(z[[t]],a[[t]],all=TRUE)
    }
  }
  z
}

