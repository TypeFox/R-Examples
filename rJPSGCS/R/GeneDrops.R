GeneDrops <-
function(ld.par="ld.par", ped="in.ped", n, complete.data=FALSE, gd.ped="gd.ped", compress.pedfiles=FALSE)
{
  # The function will create file gd.ped in the current working directory

  args<-c(ld.par,ped,n,gd.ped)
  if(complete.data) {
    args<-c(args,"-a")
  }
  if(compress.pedfiles) {
    args<-c(args,"-z")
  }

  .jcall("GeneDrops","V","main", .jarray(args))

}
