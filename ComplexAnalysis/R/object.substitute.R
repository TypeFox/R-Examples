object.substitute <-
function(func,replaced,replacement){
B.List<-NULL
for(kk in 1:length(func)){
if(length(replaced)!=length(replacement)) {stop("Length of replaced not match that of replacement")}
alph<-c(letters,toupper(letters))
num<-c(0:9,".")
symb<-c("(",")","|","+","-","*","/","~","^",">","<","=",",",":","%")
sqb<-c("[","]")
if(sum(unlist(strsplit(replaced,""))%in%c("[","]"))>0){subset.q<-TRUE}else{subset.q<-FALSE}
if(subset.q==TRUE){sqbX<-sqb}else{sqbX<-NULL}
#split string and indexing
func[kk]<-paste(strsplit(func[kk]," ")[[1]],collapse="",sep="")
td1<-c("|",strsplit(func[kk],"")[[1]],"|")
names(td1)<-seq(1,length(td1),1)
i<-1;j1<-1;j3<-1;Break<-NULL
repeat{
   if(td1[i]%in%alph){
    Break<-c(Break,paste(td1[j3:(i-1)],collapse="",sep=""))
       j1<-i
       repeat{
       if(subset.q==TRUE & td1[j1+1]=="["){ j2<-j1; repeat{if(td1[j2+1]=="]") break; j2<-j2+1};  Break<-c(Break,paste(td1[i:(j2+1)],collapse="",sep="")); j1<-j2+1}
       if(td1[j1+1]%in%c(symb,"[","]")) break
       j1<-j1+1}
       if(td1[j1]!="]"){Break<-c(Break,paste(td1[i:j1],collapse="",sep=""))}
    i<-j1+1;j3<-j1+1
    }else{i<-i+1}
if(i==length(td1)) break
}
Break<-c(Break,paste(td1[(j1+1):i],collapse="",sep=""))
for(i in 1:length(Break)){
for(j in 1:length(replaced)){
if(Break[i]==replaced[j]){Break[i]<-replacement[j]; break}
}}
Break<-paste(Break,collapse="",sep="");Break<-strsplit(Break,"[|]")[[1]];Break<-paste(Break,collapse="",sep="")
B.List<-c(B.List,Break)
}
return(B.List)}
