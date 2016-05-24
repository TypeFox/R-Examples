`survival` <-
function(time, status){

if (is.factor(time)) time<-as.numeric(levels(time))[time]
if (is.factor(status)) status<-as.numeric(levels(status))[status]
if (sum(is.na(time))==length(time)) stop("Time variable is tota de missings\n")
if (sum(is.na(status))==length(status)) stop("Status variable <U+00E9>s tota de missings\n")
if (sum(time<0)>0) stop("Time variable can not take negative values.\n")
if (length(unique(status))>2) stop("Status variable takes more than two different values.\n")
s<-paste(time,status,sep="/")
class(s)<-"character"
s
}

