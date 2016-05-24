combine.select <-
 function(sel1=NULL, sel2=NULL, ..., operator="AND", verbose=TRUE) {
   cl <- match.call()
   sels <- list(sel1, sel2, ...)

   if(any(sapply(sels, function(x) !is.null(x) && !inherits(x, "select"))))
      stop("Invalid atom select(s)")

   rm.inds = sapply(sels, is.null)
   sels = sels[!rm.inds]

   if(length(sels) == 0)
      return(NULL)
   else if(length(sels) == 1)
      return(sels[[1]])

   op.tbl <- c(rep("AND",3), rep("OR",4), rep("NOT",4))
   operator <- op.tbl[match(operator, c("AND","and","&","OR","or","|","+","NOT","not","!","-"))]

   # Print message
   if(verbose) {
      msg <- switch(operator,
        AND = " Intersect of selects",
        OR  = " Union of selects",
        NOT = " Select 2 (, 3, ...) is subtracted from select 1",
        stop("Unknown operation") )
      cat(msg, "\n", sep="")
   }

   sel <- sels[[1]]$atom
   for(i in 2:length(sels)) {
      sel <-
         switch(operator,
         "AND" = intersect(sel, sels[[i]]$atom),
         "OR"  = sort(union(sel, sels[[i]]$atom)),
         "NOT" = setdiff(sel, sels[[i]]$atom),
         stop("Unknown operation") )
   }
   sel <- list(atom = sel, xyz = atom2xyz(sel), call=cl)
   if(verbose) { cat(paste(" *  Selected a total of:", length(sel$atom), "atoms  *\n")) }

   class(sel) = "select"
   return(sel)
}
