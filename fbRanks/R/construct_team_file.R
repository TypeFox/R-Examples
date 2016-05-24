###################################
# Construct a team name resolver from the scores file
#  ignore=c("b00","bu12", "u12", " 00", "'", "-", "[(]2000 boys[)]","2000B","2012","RCL")
###################################
construct.team.file=function(scores.file="scores.csv", add.to.name="", file="team_template.csv", ignore=""){

#read in the scores files (match files)
scores=data.frame()
for(filename in scores.file){
  if(length(scores)==0){ 
    #scores=read.csv(file=filename, colClasses=c("character"),strip.white=TRUE,stringsAsFactors=FALSE)
    scores=read.csv(file=filename,strip.white=TRUE,stringsAsFactors=FALSE)
  }else{
    #next.f=read.csv(file=filename, colClasses=c("character"),strip.white=TRUE)
    next.f=read.csv(file=filename, strip.white=TRUE,stringsAsFactors=FALSE)
    if(any(!(names(next.f) %in% names(scores))))
      scores[names(next.f)[!(names(next.f) %in% names(scores))]]=NA
    if(any(!(names(scores) %in% names(next.f))))
      next.f[names(scores)[!(names(scores) %in% names(next.f))]]=NA
    scores=rbind(scores,next.f)
  } 
}
scores$home.team=str_strip.white(scores$home.team)
scores$away.team=str_strip.white(scores$away.team)

home.names = as.character(scores$home.team)
away.names = as.character(scores$away.team)
orig.team.names = unique(c(home.names, away.names))
orig.team.names = sapply(orig.team.names, str_trim, side="both")


#standardize the names
stan.names=function(x){
 #x=paste(" ",x," ",sep="")
 for(pattern in ignore)
 x=gsub(pattern, " ", x, ignore.case = TRUE)
 x=gsub("  ", " ", x, ignore.case = TRUE)
 x=gsub("  ", " ", x, ignore.case = TRUE)
 x=str_trim(x, side = "both")
 x
}

simpleCap <- function(x) {
    s <- strsplit(x, " ")[[1]]
    paste(toupper(substring(s, 1,1)), tolower(substring(s, 2)),
          sep="", collapse=" ")
}

#get rid of any leading or trailing white space
team.names = sapply(orig.team.names, stan.names)

display.name=c()
alt.names=list()
tmp.team.names = team.names
tmp.orig.team.names = orig.team.names
while(length(tmp.team.names)>0){
 x = tmp.team.names[1]
 dis=c()
 m=agrep(x,tmp.team.names,ignore.case=TRUE,max.distance=0)
 display.name=c(display.name,x)
 alt.names[[x]]=unique(c(tmp.orig.team.names[m]))
 #get rid up used up names
 tmp.team.names=tmp.team.names[-m]
 tmp.orig.team.names=tmp.orig.team.names[-m]
}

#set up the dataframe
max.num.alt.names = max(unlist(lapply(alt.names,length)))
dummy = matrix("",1,2+max.num.alt.names)
colnames(dummy)=c("name","league",paste("alt.name.",1:max.num.alt.names,sep=""))
dummy = as.data.frame(dummy,stringsAsFactors = FALSE)
dummy[1,3:(max.num.alt.names+2)]=""
team.info=dummy
disp.name = stan.names(display.name[1])
if(str_length(disp.name)==dim(str_locate_all(disp.name,"[[:upper:][:digit:][:blank:]]")[[1]])[1])
 disp.name=simpleCap(disp.name)
team.info[1,1]=disp.name
team.info[1,3:(length(alt.names[[1]])+2)]=alt.names[[1]]
for(i in 2:length(display.name)){
team.info = rbind(team.info,dummy)
disp.name = stan.names(display.name[i])
if(str_length(disp.name)==dim(str_locate_all(disp.name,"[[:upper:][:digit:][:blank:]]")[[1]])[1]){
 disp.name=simpleCap(disp.name)
}
team.info$name[i]=disp.name
team.info[i,3:(length(alt.names[[i]])+2)]=alt.names[[i]]
}

ix=sort(team.info$name,index.return=TRUE)$ix
team.info=team.info[ix,,drop=FALSE]
rownames(team.info)=NULL
team.info$name=paste(team.info$name,add.to.name,sep="")
team.info$name=str_proper(team.info$name)
team.info$name=str_replace(team.info$name," Fc"," FC")
team.info$name=str_replace(team.info$name," Sc "," SC ")
team.info$name=str_replace(team.info$name,"Soccer Club","SC")

write.table(team.info, file=file,row.names=FALSE,col.names=FALSE,append=FALSE,sep=",",qmethod="double")
}

