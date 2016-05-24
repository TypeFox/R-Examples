annotatebox<-function(Rbox=matrix(ncol=4, nrow=8)  , add=TRUE)
{
if(missing(add)) { add=TRUE }
if(missing(Rbox)) { add=FALSE }



  if(!add)
    {
      S= stressSETUP()
      pstart()

    
     PLOTbox(S$Rax, S$Rbox, axcol= 'green', boxcol= 'purple')
      
      Rbox = S$Rbox
    }

points(Rbox[,1],Rbox[,2], col='blue')
text(Rbox[,1],Rbox[,2], labels=1:8, pos=1)


}
