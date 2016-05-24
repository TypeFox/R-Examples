alter_files <-
function(change_params){
#library(operators)
filetype=NULL
for(ft in unique(change_params$filetype)){
        print(ft)
        files=list.files(,paste(ft,"$",sep=""))
        for (file in files) {
                fileorig=paste(file,".unixorig",sep="");
                junk%<%fileorig
                file_change_params=subset(change_params,filetype==ft)
                for( i in 1:length(rownames(file_change_params))){
                        startstr=file_change_params[i,"startstr"]
                        endstr=file_change_params[i,"endstr"]
                        current=file_change_params[i,"current"]
                        param=file_change_params[i,"parameter"]
                        multi=file_change_params[i,"multi"]
                        alter_type=file_change_params[i,"alter_type"]
                        frformat=unlist(strsplit(as.character(file_change_params[i,"frformat"]),","))
                        fwformat=unlist(strsplit(as.character(file_change_params[i,"fwformat"]),","))
                                junkline=grep(param,junk,value=T)
                                fread_tmp=read.fortran(con1<-textConnection(junkline),frformat,as.is=T)
                                close(con1)
                                fwrite_tmp=as.vector(1:length(fread_tmp[1,]))
                                for (j in 1:length(fread_tmp[1,])) {
				   lastj=j
                                   if(is.na(fread_tmp[, j])){lastj=j-1;break}
                                   if(is.numeric(fread_tmp[,j])){
                                      if(alter_type=="percent"){
				          fwrite_tmp[j]=sprintf(fwformat[min(j,length(fwformat))],fread_tmp[,j]*current)
                                      } else {
                                          fwrite_tmp[j]=sprintf(fwformat[min(j,length(fwformat))],as.numeric(current))
                                      } 
                                   } else {
                                      fwrite_tmp[j]=sprintf(fwformat[j],fread_tmp[,j])
                                   }
                                   stringb=paste(fwrite_tmp[1:lastj],sep="",collapse="")
                                   junk=gsub(paste(".*",param,".*",sep=""),stringb,junk)
                                }
                }
                cat(junk,file=file,sep="\n")
        }
}
}

