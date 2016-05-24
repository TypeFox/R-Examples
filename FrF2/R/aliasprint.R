aliasprint <- function(design,...){
   ##function to print but not return a design's alias information
   if(!"design" %in% class(design)) stop("design must be of class design")
   lgd <- design.info(design)$aliased[[1]]
   als <- design.info(design)$aliased[-1]
   if (all(names(als)%in% c("main","fi2"))) print(design.info(design)$aliased, quote=FALSE, ...)
   else{
   if (length(als)>1 | length(als[1])>1){
            alprint <- list(legend=lgd, aliases=t(t(sapply(als, "paste", collapse=" = "))))
            colnames(alprint$aliases) <- ""
            rownames(alprint$aliases) <- rep("",nrow(alprint$aliases))
       if (is.null(alprint$legend)) print(alprint$aliases) else{
            if(!is.null(alprint$legend)) names(alprint$legend)<-rep("", length(alprint$legend))
            print(alprint, quote=FALSE, ...)
       }}
   else print(design.info(design)$aliased, quote=FALSE, ...)
   }
}