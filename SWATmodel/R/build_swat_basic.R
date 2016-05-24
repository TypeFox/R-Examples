# Function for example
# 
#
build_swat_basic<-function(dirname,iyr=2000,nbyr=2,wsarea=45,elev=500,declat=45.7,declon=-76,hist_wx=NULL){
#
# Removing the global variables Note. Delete the first line when past supported versions are past v15.1
#
swat_general<- NULL; rm(swat_general);  # A code dummy
# Correct way to let code know about global variables post 2.15.1
#if(getRversion() >= "2.15.1") {
#   globalVariables(c("swat_general"))
#} 
#
#data(swat_general)
load(paste(path.package("SWATmodel"), "data/swat_general.rda", sep = "/"))
dir.create(dirname)
setwd(dirname)
for (file in names(swat_general)) {
    print(file)
    cat(unlist(swat_general[file]), file = file, sep = "\n")
}

if (length(hist_wx) > 2) {
  tmp_head=paste("Tmp\nLati Not Used\nLong Not Used\nElev        ",sprintf("%5.0f\n",elev),sep="")
  pcp_head=paste("Pcp\nLati Not Used\nLong Not Used\nElev        ",sprintf("%5.0f\n",elev),sep="")
  cat(tmp_head,sprintf("%s%005.1f%005.1f\n",format(hist_wx$DATE,"%Y%j"),hist_wx$TMX,hist_wx$TMN),file="tmp.tmp",sep="")
  cat(pcp_head,sprintf("%s%005.1f\n",format(hist_wx$DATE,"%Y%j"),hist_wx$PRECIP),file="pcp.pcp",sep="")
}

#
# Setting up the filecio file
#
filecio=readLines("file.cio")

if (length(hist_wx) > 2) {
  filecio=sub("pcp1\\.pcp","pcp\\.pcp",filecio)
  filecio=sub("tmp1\\.tmp","tmp\\.tmp",filecio)
}

filecio=sub("\\d\\d\\d\\d(    \\| IYR)",paste(iyr,"\\1",sep=""),filecio)
substr=sprintf("      %10d    | NBYR",nbyr)
filecio=sub("^.*\\| NBYR",substr,filecio,perl=T)
cat(filecio,file="file.cio",sep="\n")
#
# Setting up the sub basin files
#
num_subbasins=length(list.files(pattern="0*.sub"))
for (subfilename in list.files(pattern="0*.sub")){
        subfile=readLines(subfilename)
        substr=sprintf("%12.1f        | SUB_KM",wsarea/num_subbasins)
        subfile=sub("^.*\\| SUB_KM",substr,subfile,perl=T)
        substr=sprintf("%12.1f        | SUB_LAT",declat)
        subfile=sub("^.*\\| SUB_LAT",substr,subfile,perl=T)
        substr=sprintf("%12.1f        | SUB_ELEV",elev)
        subfile=sub("^.*\\| SUB_ELEV",substr,subfile,perl=T)
        cat(subfile,file=subfilename,sep="\n")
}

}
# End build_swat_basic


