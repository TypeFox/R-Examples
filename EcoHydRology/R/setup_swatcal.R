setup_swatcal <-
function(change_params){

for(ft in unique(change_params$filetype)){
        print(ft)
        files=list.files(,paste(ft,"$",sep=""))
        for (file in files) {
                junk%<%file
                junk=gsub("\r","",junk)
                file_swatcal=paste(file,".unixorig",sep="");
                cat(junk,file=file_swatcal,sep="\n")
        }
}
}

