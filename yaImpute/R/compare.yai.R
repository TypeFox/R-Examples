compare.yai = function (...,ancillaryData=NULL,vars=NULL,method="rmsd",scale=TRUE)
{
   if (missing(...)) stop ("... required")
   okClasses <- c("yai","impute.yai")
   for (object in list(...))
      if (length(intersect(class(object),okClasses))==0)
         stop("object classes must be one of ",paste(okClasses,collapse=", "))
   args <- list(...)
   if (length(intersect(method,c("rmsd","cor")))==0) stop("method must be rmsd or cor")
   names(args) <- as.list(substitute(list(...)))[-1]  #who would of guessed that this is the way!
   ans <- list()
   scales <- list()
   tag <- if (method=="rmsd") "rmsdS" else "cor"
   meth <- match(method,c("rmsd","cor"))
   i <- 0
   for (object in list(...))
   {
      i <- i+1
      if (inherits(object,"yai")) object <- impute.yai(object,ancillaryData=ancillaryData,vars=vars,observed=TRUE)
      one <- switch(meth,
                    rmsd.yai(object,vars=vars,scale=scale),
                    cor.yai(object,vars=vars), NULL)
      names(one) <- paste(names(args)[i],tag,sep=".")
      ans[[i]] <- one
      if (meth == 1) scales[[i]] <- attributes(one)$scale
   }
   names(ans) <- names(args)
   ans <- unionDataJoin(ans)
   class(ans) <- c("compare.yai",class(ans))
   if (length(scales) > 0) 
   {
     names(scales) <- names(args)
     ident = TRUE
     if (length(scales) > 1) 
     {
       s1 <- scales[[1]]
       for (i in 2:length(scales))
       {  
         if (!identical (s1,scales[[i]]))
         {
           warning ("not all scale factors are the same.")
           ident = FALSE
           break         
         }
       }
     }
     attr(ans,"scales") <- if (ident) scales[[1]] else scales
   }
   ans
}
