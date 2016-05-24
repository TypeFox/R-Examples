annotateplane<-function(Rp=matrix(ncol=4, nrow=3)  , add=TRUE)
{
  if(missing(add)) { add=TRUE }
  if(missing(Rp)) { add=FALSE }

  if(!add)
    {
      S= stressSETUP()
      pstart()

      
      PLOTbox(S$Rax, S$Rbox, axcol= 'green', boxcol= 'purple')
      PLOTplane(S$Rp, planecol="brown")
      
      Rp = S$Rp
      Rbox = S$Rbox
    }

  points(Rp[,1],Rp[,2], col='blue')
  text(Rp[,1],Rp[,2], labels=1:3, pos=1)

}
