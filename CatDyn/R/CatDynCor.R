CatDynCor <-
function(x, ttl, meths, arr)
  {
   if(length(x) != length(ttl))
     {stop("'x' and 'ttl' should be of the same length")}
   if(class(x) != "list")
     {stop("'x' should be a list with each component an output of function CatDynFit(), class 'catdyn'")}
   if(sum(sapply(x,class) == "catdyn") != length(x))
     {stop("'x' should be a list with each component an output of function CatDynFit(), class 'catdyn'")}
   if(class(ttl) != 'character')
     {stop("'ttl' should be a character vector")} 
   if(length(meths) != length(x))
      {stop("One numerical method for each component of 'x' should be provided")}
   if(length(arr) != 2)
     {stop("'arr' is the number of rows and number of columns to deploy the panels. Passed as is to 'par'")}
   if(arr[1]*arr[2] < length(x))
     {warning("Not all model fits are shown")}
   if(length(x) == 1)
     {
      par(mfrow=arr,mar=c(2,3,1,1),oma=c(2,3,1,1))
      hist(unique(x[[1]]$Model[[meths]]$Cor[-diag(x[[1]]$Model[[meths]]$Cor)]),main=ttl,xlab="",cex.axis=1.5)
     }
   else
     {
      par(mfrow=arr,mar=c(2,3,1,1),oma=c(2,3,1,1))
      for(i in 1:length(x))
        {
         hist(unique(x[[i]]$Model[[meths[i]]]$Cor[-diag(x[[i]]$Model[[meths[i]]]$Cor)]),main=ttl[i],xlab="",cex.axis=1.5)
        }
     }
   mtext(side=1,outer=TRUE,text="Pairwise Correlation Coefficients",cex=1.25,line=1)
   mtext(side=2,outer=TRUE,text="Frequency",cex=1.25,line=1)
  }
