kin.morgan<-function(ped,verbose=FALSE)
{
   v2k <- function (pars) 
   # 5/6/2004
   # this function is in spirit similar to g2a and make.del(mvnmle)
   {
       k <- floor((-1 + sqrt(1 + 8 * length(pars)))/2)
       mymatrix <- diag(1:k)
       if (k > 1) {
           for (i in 1:k) {
               mymatrix[i, 1:i] <- mymatrix[1:i,i] <- pars[1:i]
               pars <- pars[-(1:i)]
           }
       }
       mymatrix
   }
   pedsize<-dim(ped)[1]
   kin<-rep(0,pedsize*(pedsize+1)/2)
   id <- ped[,1]
   father <- ped[,2]
   mother <- ped[,3]
   tid <- c(id,father,mother)
   uid <- unique(tid[tid!=0])
   iid <- match(id,uid)
   fid <- match(father,uid,nomatch=0)
   mid <- match(mother,uid,nomatch=0)
   peddata <- rbind(id,father,mother)
   pedindex <- rbind(iid,fid,mid)
   if (verbose)
   {
      cat("The original pedigree IDs and their indices:\n")
      print(t(rbind(peddata,pedindex)))
   }
   z<-.C("kin_morgan",data=as.integer(peddata),pedsize=as.integer(pedsize),
         pedinex=as.integer(pedindex),kin=as.double(array(kin)),PACKAGE="gap")
   kin.matrix=v2k(z$kin)
   colnames(kin.matrix) <- rownames(kin.matrix) <- id
   list(kin=z$kin,kin.matrix=kin.matrix)
}
