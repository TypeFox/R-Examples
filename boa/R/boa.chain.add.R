"boa.chain.add" <-
function(link, lname)
#-------------------------------------------------------------------------------
#-------------------------------------------------------------------------------
{
   added <- FALSE
   if(is.matrix(link) && is.numeric(link)) {
      added <- TRUE
      dn <- dimnames(link)
      if(is.null(dn)) {
         dimnames(link) <- list(1:nrow(link),
                                paste("par", 1:ncol(link), sep = ""))
      } else if(length(dn[[1]]) == 0) {
         dimnames(link)[[1]] <- 1:nrow(link)
      } else if(length(dn[[2]]) == 0) {
         dimnames(link)[[2]] <- paste("par", 1:ncol(link), sep = "")
      }
      master <- boa.chain("master")
      master.support <- boa.chain("master.support")
      # master[[lname]] <- boa.sortparms(link)
      master[[lname]] <- link
      master.support[[lname]] <- matrix(c(-Inf, Inf), nrow = 2,
                                        ncol = ncol(link))
      dimnames(master.support[[lname]]) <- list(c("Min", "Max"),
                                                boa.pnames(master[[lname]]))
      if(boa.chain("work.sync")) {
         boa.chain(master = master, master.support = master.support,
                   work = master, work.support = master.support)
      } else {
         boa.chain(master = master, master.support = master.support)
      }
   } else {
      cat("Warning: object must be a numeric matrix\n")
   }

   return(added)
}
