format_corrector <-
function(table,identif=NULL,force=FALSE,rate.miss.date=0.5){

  options(warn=-1) 

  table<-as.data.frame(table)  

  original.table<-table    

  formats<-rep(NA,ncol(table))
  lost.table<-white.table<-lost.table.numeric<-lost.table.old<-string.table<-list()
                                
  for (j in 1:ncol(table)){ # loop for all variables
    

    if (force || is.null(attr(table[,j],"fixed.formats")) || !attr(table[,j],"fixed.formats")){   

      cat("\n\n-----Fixing variable '",names(table)[j],"'---------\n")
      
      old.format<-attr(table[,j],"format.SPSS")
      
      lost.table[[j]]<-white.table[[j]]<-lost.table.old[[j]]<-string.table[[j]]<-0
      
      var<-as.character(table[,j])    
      var<-trim(var)
      var<-ifelse(var=="",NA,var)  # 
      
      if (sum(var=="NA",na.rm=TRUE)==nrow(table)) var<-rep(NA,nrow(table))
      
      assigned=FALSE

      ### Empty var ###
      
      if (!assigned & sum(is.na(var))==nrow(table)){
        formats[j]<-"F8.2"
        table[,j]<-var
        assigned=TRUE
      }
      
      ### Numeric var ###
      
      if (!assigned){
        aux<-var
        index.blank<-grep(" ",aux)
        aux<-gsub(" ","",aux)
        aux<-gsub(",",".",aux)
        aux<-gsub("[\\.]+",".",aux)  
        aux<-ifelse(aux==".",NA,aux)  
        na.before<-is.na(aux)
        na.after<-is.na(as.double(aux))
        if (all(is.na(aux))){ 
          table[,j]<-aux
          formats[j]<-"F8.2"
          assigned=TRUE              
        }else{
          if (sum(na.before)==sum(na.after)){   
            var<-as.double(aux)
            integer<-trunc(var)
            dec<-var-trunc(var)
  			    num.integer<-max(nchar(as.character(integer[!is.na(integer)])))
    			  num.dec<-max(nchar(format(dec[!is.na(dec)])))-2
    			  num.dec<-ifelse(num.dec<0,0,num.dec)
   			    num.dec<-ifelse(num.dec>16,16,num.dec)  #max 16 dec
    			  formats[j]<-paste("F",num.integer+num.dec+1,".",num.dec,sep="")
            table[,j]<-var
      		  assigned=TRUE
      		  if (length(index.blank)>0){ 
              if (!is.null(identif)){
                white.table[[j]]<-eval(parse(text=paste("data.frame(list(",identif,"=table[index.blank,identif],ara=table[index.blank,j],abans=original.table[index.blank,j]))",sep="")))
              }
              if (is.null(identif)){
                white.table[[j]]<-data.frame(list(row_num=index.blank,ara=table[index.blank,j],abans=original.table[index.blank,j]))
              }  
              cat("\nEmpty fields have appeared which have been removed in the variable for the individuals:\n")
              print(white.table[[j]])
              cat("\n\n")
            }
         } 
        }
      }
   

      ### Date var ###
      
      if (!assigned){

        var.class<-class(table[,j])
        if (is.null(var.class) || !sum(var.class%in%c("dates","times"))==2){ 
        
          na.before<-is.na(var)  
        
          # format  /.
          aux<-var
          aux<-gsub("  "," ",aux)  
          aux<-gsub("\\)","",aux)
          aux<-gsub("\\(","",aux)
          aux<-gsub("-","/",aux)
          aux<-gsub("\\.","/",aux)
                    
          ## day
          first<-unlist(lapply(strsplit(aux,"/"),function(x) x[1]))
          ## month
          second<-unlist(lapply(strsplit(aux,"/"),function(x) x[2]))
            
          second<-tolower(second)
            # English
          second[grep("^jan",second)]<-1
          second[grep("^feb",second)]<-2
          second[grep("^mar",second)]<-3          
          second[grep("^apr",second)]<-4
          second[grep("^may",second)]<-5
          second[grep("^jun",second)]<-6
          second[grep("^jul",second)]<-7
          second[grep("^aug",second)]<-8
          second[grep("^sep",second)]<-9          
          second[grep("^oct",second)]<-10
          second[grep("^nov",second)]<-11
          second[grep("^dec",second)]<-12
             # Spanish
          second[grep("^ene",second)]<-1
          second[grep("^feb",second)]<-2
          second[grep("^mar",second)]<-3          
          second[grep("^abr",second)]<-4
          second[grep("^may",second)]<-5
          second[grep("^jun",second)]<-6
          second[grep("^jul",second)]<-7
          second[grep("^ago",second)]<-8
          second[grep("^sep",second)]<-9          
          second[grep("^oct",second)]<-10
          second[grep("^nov",second)]<-11
          second[grep("^dic",second)]<-12
             # Catalan
          second[grep("^gen",second)]<-1
          second[grep("^feb",second)]<-2
          second[grep("^mar",second)]<-3          
          second[grep("^abr",second)]<-4
          second[grep("^mai",second)]<-5
          second[grep("^jun",second)]<-6
          second[grep("^jul",second)]<-7
          second[grep("^ago",second)]<-8
          second[grep("^set",second)]<-9          
          second[grep("^oct",second)]<-10
          second[grep("^nov",second)]<-11
          second[grep("^dic",second)]<-12          
          
          third.aux<-unlist(lapply(strsplit(aux,"/"),function(x) x[3]))
          third.aux<-trim(third.aux)  
          third<-unlist(lapply(strsplit(third.aux," "),function(x) x[1]))
          if (length(grep(" ",third.aux))>0) fourth<-unlist(lapply(strsplit(third.aux," "),function(x) x[2]))     

          #individuals with ' / /1999', o ' // 'dates remain missing.
          data.missings<-apply(is.na(cbind(as.integer(first),as.integer(second),as.integer(third))),1,any)     
          if (sum(data.missings)>0) first[data.missings]<-second[data.missings]<-third[data.missings]<-rep(NA,sum(data.missings))

          first<-as.integer(first)
          second<-as.integer(second)
          third<-as.integer(third)
          
          if (sum(is.na(first))<nrow(table) && sum(is.na(second))<nrow(table) && sum(is.na(third))<nrow(table)){ 
          
            aux2<-paste(as.character(first),as.character(second),as.character(third),sep="/")
            
            format.data<-matrix(NA,nrow(table),6)

            format.data[,1]<-ifelse(first%in%c(1:31) & second%in%1:12 & (third>31 | third<10),1,0)  
            format.data[,2]<-ifelse(first%in%c(1:12) & second%in%1:31 & (third>31 | third<10),1,0)    
            format.data[,3]<-ifelse((first>31 | first<10) & second%in%1:12 & third%in%1:31,1,0)                
            format.data[,4]<-ifelse((first>31 | first<10) & second%in%1:31 & third%in%1:12,1,0)  
            format.data[,5]<-ifelse(first%in%1:12 & (second>31 | second<10) & third%in%1:31,1,0)   
            format.data[,6]<-ifelse(first%in%1:31 & (second>31 | second<10) & third%in%1:12,1,0)  

            colnames(format.data)<-c("d/m/y","m/d/y","y/m/d","y/d/m","m/y/d","d/y/m")
            possible.data.formats<-colnames(format.data)[apply(format.data==1,2,any)]
          
            format.data<-names(sort(colSums(format.data),decreasing = TRUE))[1] # The most frequent format is selected.
            
            dates<-dates(aux2,format=c(dates=format.data),out.format=c(dates="day-mon-year"))

            
            
            if (inherits(try(get("fourth"),silent=TRUE), "try-error")){ ## format date
              var<-chron(dates.=dates,out.format=c(dates="day-mon-year"))
              table[,j]<-var
              formats[j]<-"DATE11"
              assigned=TRUE        

            } else {     ### format date hour

              fourth<-gsub("\\.",":",fourth)  # separats per ':'.
              time.first<-unlist(lapply(strsplit(fourth,":"),function(x) x[1]))
              time.second<-unlist(lapply(strsplit(fourth,":"),function(x) x[2]))
              time.third<-unlist(lapply(strsplit(fourth,":"),function(x) x[3]))
              
              if (!any(is.na(as.double(time.first)) || is.na(as.double(time.second)) || is.na(as.double(time.third)))){ 
               
                time.first<-as.integer(time.first)
                time.second<-as.integer(time.second)
                time.third<-as.integer(time.third) 
      
                format.time<-NULL
                aux2<-paste(as.character(time.first),as.character(time.second),as.character(time.third),sep=":")
                
                time.first.aux<-time.first[!is.na(time.first)]
                time.second.aux<-time.second[!is.na(time.second)]
                time.third.aux<-time.third[!is.na(time.third)]
                if (all(time.first.aux%in%0:60) & all(time.second.aux%in%0:60) & all(time.third.aux%in%0:24))  {times<-times(aux2,format=c(times = "s:m:h"));format.time=c(format.time,"s:m:h")}
                if (all(time.first.aux%in%0:24) & all(time.second.aux%in%0:60) & all(time.third.aux%in%0:60))  {times<-times(aux2,format=c(times = "h:m:s"));format.time=c(format.time,"h:m:s")}
                if (all(time.first.aux%in%0:24) & all(time.second.aux%in%0:60) & all(time.third.aux%in%0:60))  {times<-times(aux2,format=c(times = "h:s:m"));format.time=c(format.time,"h:s:m")}
                if (all(time.first.aux%in%0:60) & all(time.second.aux%in%0:60) & all(time.third.aux%in%0:24))  {times<-times(aux2,format=c(times = "m:s:h"));format.time=c(format.time,"m:s:h")}  
                if (all(time.first.aux%in%0:60) & all(time.second.aux%in%0:24) & all(time.third.aux%in%0:60))  {times<-times(aux2,format=c(times = "m:h:s"));format.time=c(format.time,"m:h:s")}              
                if (all(time.first.aux%in%0:60) & all(time.second.aux%in%0:24) & all(time.third.aux%in%0:60))  {times<-times(aux2,format=c(times = "s:h:m"));format.time=c(format.time,"s:h:m")}
                
                if (!is.null(format.time)){ 
                  if (length(format.time)>1){
                    format.time<-format.time[1]  # the first one is selected
                    times<-times(aux2,format=c(times = format.time),out.format=c(times="h:m:s"))
                  }
                  var<-chron(dates,times)
                  table[,j]<-var
                  formats[j]<-"DATETIME23"
                  assigned=TRUE
                } else{
                  cat("WARNING: The variable have impossible seconds, minutes or hours so only the days will be stored","\n\n") # remains as character
                  var<-chron(dates)
                  table[,j]<-var
                  formats[j]<-"DATE11"
                  assigned=TRUE
                }
                
              }
            }
          } 

        } else {
          
          var<-table[,j]
            
          if (length(attributes(var)$format)>1){ 
            var<-chron(var,out.format=c(dates="day-mon-year",times="h:m:s"))
            table[,j]<-var
            formats[j]<-"DATETIME23"
            assigned=TRUE
          } else { 
            var<-chron(var,out.format=c(dates="day-mon-year"))
            table[,j]<-var
            formats[j]<-"DATE11"
            assigned=TRUE
          }
        }
        
        if (assigned) na.after<-is.na(var)
        
        if (assigned && sum(na.after)>sum(na.before)){  
          if ((sum(na.after)-sum(na.before))/sum(!na.before)>rate.miss.date){     
            assigned=FALSE
            var<-original.table[,j] 
            cat("Aborted variable conversion due to the loss of more than",rate.miss.date*100,"% of the fields")
          }else{
          
            ######  Observations lost during the date conversion #######
            lost.list<-original.table[na.after!=na.before,j]
            if (!is.null(identif)){
              num_indiv=table[na.after!=na.before,identif]
              lost.table[[j]]<-eval(parse(text=paste("data.frame(list(",identif,"=num_indiv,values=lost.list))",sep="")))
            }
            if (is.null(identif)){              
              num_indiv=which(na.after!=na.before)
              lost.table[[j]]<-data.frame(list(row_num=num_indiv,values=lost.list))
            }
            cat("\n\n",nrow(lost.table),"variable fields were lost during the date format conversion\n\n")
            print(lost.table[[j]])
            cat("\n\n")    
            
          }
        }
        
        if (assigned){
          lost.list.old<-(var<chron(dates.="14/10/1582",format=c(dates="d/m/y")))
          lost.list.old<-ifelse(is.na(lost.list.old),FALSE,lost.list.old)
        }
        if (assigned && sum(lost.list.old)>0){
          if (!is.null(identif)){
            num_indiv=table[lost.list.old,identif]
            lost.table.old[[j]]<-eval(parse(text=paste("data.frame(list(",identif,"=num_indiv,values=original.table[lost.list.old,j]))",sep="")))
          }
          if (is.null(identif)){
            num_indiv=which(lost.list.old)
            lost.table.old[[j]]<-data.frame(list(row_num=num_indiv,values=as.character(original.table[lost.list.old,j])))
          }              
          cat("\n\nFor the following individuals, during the SPSS export the variable fields will be lost due to having older dates than 14/10/1582:\n")
          print(lost.table.old[[j]])
          cat("\n\n")
        }
      }

      
      ### Character var
      
      if (!assigned){ ## remains character
        table[,j]<-as.character(var)
        formats[j]<-paste("A",max(nchar(as.character(table[,j]))),sep="")
        index<-unique(c(grep("\n",table[,j]),grep("\t",table[,j]),grep("\r",table[,j])))
        if (length(index)>0){
          var.old<-table[index,j]  
          table[,j]<-gsub("\n"," ",table[,j]) 
          table[,j]<-gsub("\t"," ",table[,j]) 
          table[,j]<-gsub("\r"," ",table[,j]) 
          if (!is.null(identif)){
            num_indiv=table[index,identif]
            string.table[[j]]<-eval(parse(text=paste("data.frame(list(",identif,"=num_indiv,values=var.old",sep="")))
          }
          if (is.null(identif)){
            num_indiv=index
            string.table[[j]]<-data.frame(list(row_num=num_indiv,values=var.old))
          }
          cat("\n\nInstructions ('\\t'), ('\\n') and ('\\r') have been removed from the following individuals:\n")
          print(string.table[[j]])
          cat("\n\n")          
        }
        assigned=TRUE
      }
      
      ## value labels and vari.label we lost are recovered.
      attr(table[,j],"value.labels")<-attr(original.table[,j],"value.labels")
      attr(table[,j],"label")<-attr(original.table[,j],"label")
      attr(table[,j],"format.SPSS")<-formats[j]
      attr(table[,j],"fixed.formats")<-TRUE     

      cat("\n   The following SPSS format has been assigned:",formats[j],"     \n")
      if (!is.null(old.format)) cat("   The older SPSS format was",old.format,"\n\n")

      if (!is.null(old.format) && (sub("[0-9.]+$","",old.format)=="A" & sub("[0-9.]+$","",formats[j])=="F" & !is.null(attr(table[,j],"value.labels")))){
        cat("\nThe variable had these value labels:\n")
        print(attr(table[,j],"value.labels"))
        attr(table[,j],"value.labels")=structure(as.double(attr(table[,j],"value.labels")),names=names(attr(table[,j],"value.labels")))
        cat("\nand these become:\n")
        print(attr(table[,j],"value.labels"))        
      }  
    }  
    
    if (!inherits(try(get("fourth"),silent=TRUE), "try-error")) rm(fourth) # si exiteix la variable fourth l'elimina.

  } ## END LOOP
    
  options(warn=0)  

  return(table) 	      		   
   
}
