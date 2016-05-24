functionall<-function(fd,num,setdis,Meth,resultname,original,randgenenum,randtime,randname,presult)
{
distable(fd,num,setdis,Meth,resultname)

randdis(original,randgenenum,randtime,setdis,num,Meth,randname)
MM<-read.table(resultname,header=TRUE)
N<-read.table(randname);
valuep(MM,N,presult)
}
