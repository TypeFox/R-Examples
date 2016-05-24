i_info<-function(hvec,itembeta,theta)
#calculates information (Samejima, 1969) for an item i as a function of theta
#
#@input: hvec...number of categories of item
#        itembeta...cumulative item parameters
#        theta ... supporting or sampling points on latent trait
#@output: a list with
#         $c.info...matrix of category information in columns for theta (rows)
#         $i.info...vector of item information at values of theta
#@author: Thomas Rusch
#@date:6.12.2013 Happy Nikolaus! 
#
  {
   if(missing(theta)) theta <-seq(-5,5,0.01)
   p.ih<-function(hvec,itembeta,theta)
   #Calculates p.ih of given item i and the weird expression in its first derivative
   #needs categories given (hvec) and the cumulative item parameters of the item (itembeta)
   #@output: a list with
   #         $p.ih...matrix of probabilities to fall into category h (colums) for given items as a function of theta (rows).
   #         $weird...the weird expression from the derivative
  {
    beta <- c(0,itembeta) #eRm gives itempar with first fixed to zero
    numerator<-exp(outer(hvec,theta)+beta) #Numerator
    tmp<-hvec*numerator
    weird.exp.num <- apply(tmp,2,sum) #numerator of weird expression in the derivative
    denom <- apply(numerator,2,sum) #denominator
    p.ih<-t(numerator)/denom    #categories in column,thetas in rows
    weird.exp<- weird.exp.num/denom #weird expression in derivative
    return(list("p.ih"=p.ih,"weird"=weird.exp))
  }

  ic.derivative<-function(hvec,itembeta,theta)
   {
   #Calculates first derivative of p.ih of given item i, needs number of categories and     cumulative item parameters
   #
   #@output: a list with
   #         $out...first derivative of p.ih with categories h in columns and theta in ro   ws
    f1<-p.ih(hvec,itembeta,theta) #to get p.ih and weird expression
    out <- t(hvec*t(f1$p.ih))-f1$p.ih*f1$weird #first derivative
    return(out)
  }

   tmp <- ic.derivative(hvec,itembeta,theta)#call ic.derivative
   c.info <- tmp^2/p.ih(hvec,itembeta,theta)$p.ih #calculates category info (columns) for all theta(rows)
   i.info <-apply(c.info,1,sum)#calculates item for all theta(rows)
   return(list("c.info"=c.info,"i.info"=i.info))
 }

item_info <- function(ermobject,theta=seq(-5,5,0.01))
##Calculates information (Samejima, 1969) of all items as a function of the latent trait, theta
#        ermobject ... object of class eRm
#        theta ... supporting or sampling points on latent trait
#@output: a list where each element corresponds to an item and contains
#         $c.info...matrix of category information in columns for theta (rows)
#         $i.info...vector of item information at values of theta
#@author: Thomas Rusch
#@date:13.6.2011
#
{
   vec.tmp <- get_item_cats(X=ermobject$X,nitems=dim(ermobject$X)[2],grp_n=dim(ermobject$X)[1])
   betapar <- ermobject$betapar
   veco <- unlist(lapply(vec.tmp,max))
   alloc.list<-vector("list",length(veco))
   hvec.list <- vector("list",length(veco))
   out.list <- vector("list",length(veco))
   for (i in 1:length(veco))
     {
       alloc.list[[i]] <- rep(i,veco[i])
       hvec.list[[i]] <- 0:veco[i]
     }
   uu<-unlist(alloc.list)
   itembeta.list <- split(betapar,uu)
   for (i in 1:length(itembeta.list))
     {
      out.list[[i]] <- i_info(hvec.list[[i]],itembeta.list[[i]],theta) #patch
    }
   return(out.list)
 }

#THANK YOU FOR READING THE SOURCE

