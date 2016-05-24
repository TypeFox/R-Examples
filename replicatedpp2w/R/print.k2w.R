print.k2w <-
function(x,...){
      if(!is.null(x$nameB)){
         btss<- c(x$btss.i, x$btss.j, x$btss.ij)
         length.p<-length(x$btss.i.res)+1
         r.i<- rank(c(x$btss.i, x$btss.i.res))[1]
         r.j<- rank(c(x$btss.j, x$btss.j.res))[1]
         r.ij<- rank(c(x$btss.ij, x$btss.ij.res))[1]
         p.value<- 1-(c(r.i, r.j, r.ij)/length.p)
         result<- data.frame(btss=btss,p.value=p.value)      
         rownames(result) <- c(x$nameA, x$nameB, paste(x$nameA,x$nameB, sep=" : "))
       }
      if(is.null(x$nameB)){
         btss<- c(x$btss.i)
         length.p<-length(x$btss.i.res)+1
         r.i<- rank(c(x$btss.i, x$btss.i.res))[1]
         p.value<- 1-(r.i/length.p)
         result<- data.frame(btss=btss,p.value=p.value)      
         rownames(result) <- x$nameA
      }

      print(result)
}
