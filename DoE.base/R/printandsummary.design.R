## controls printing of the design
## especially with structure information
print.design <- function(x,show.order=NULL, group.print=TRUE, std.order=FALSE, ...){
   if (!"design" %in% class(x)) stop("this function works for class design objects only")
      else{ ## do the right thing for class design from package conf.design
         if (is.null(design.info(x)))
             print.data.frame(x, ...)
         else{
             ## now designs generated with suite DoE.base etc.
             ## this else closes at the end of the function
   di <- design.info(x)
   if (std.order) {
        print(cbind(run.order(x)[,c(1,2)],x)[ord(run.order(x)),])
        cat("NOTE: columns run.no.in.std.order and run.no are annotation, not part of the data frame",fill=TRUE)
         if (length(grep("param",di$type))>0 & length(grep("wide",di$type))>0 ){
             cat("Outer array:\n")
             print(di$outer, ...)
      }
    return(invisible())
   }
   if (group.print)
   group.print <- di$type %in% c("full factorial.blocked", "FrF2.blocked", "FrF2.blockedcenter", 
        "FrF2.splitplot", "FrF2.splitplot.folded", "Dopt.blocked", "Dopt.splitplot", "oa.blocked", "pb.blocked")
         # | length(grep("param",di$type)) > 0
   if (is.null(show.order)) 
       show.order <- group.print | di$replications > 1 | di$type=="crossed" | length(grep("param",di$type)) > 0
   if (show.order){
       if (!group.print)
       print(cbind(run.order(x)[,2:3],x), ...)
       else{
          ## provisions for some Dopt types; not yet known whether all of them will exist 
          if (di$type %in% c("full factorial.blocked", "FrF2.blocked", "FrF2.blockedcenter", "Dopt.blocked", "oa.blocked", "pb.blocked"))
             printBy(cbind(run.order(x)[,2:3],x), di$block.name,...)
          if (di$type %in% c("FrF2.splitplot", "FrF2.splitplot.folded","Dopt.splitplot"))
             printBy(cbind(run.order(x)[,2:3],x), names(di$factor.names)[1:di$nfac.WP], ...)
         ## must find something more convenient
         ## for many crossed and parameter designs, group printing is a nuisance only
         ## should not be done for wide designs
         ## may be useful for long designs, but only if there are more than two rows each
         # if (di$type == "crossed" | length(grep("param",di$type)) > 0)
         #    printBy(cbind(run.order(x)[,2:3],x), 
         #                names(di$factor.names)[1:sum(di$cross.nfactors[-length(di$cross.nfactors)])])
       }
       }
   else {
      if (!group.print)
          print(cbind(x), ...)
      else
       {
          if (di$type %in% c("full factorial.blocked", "FrF2.blocked", "FrF2.blockedcenter", "Dopt.blocked"))
             printBy(cbind(x), di$block.name,...)
          if (di$type %in% c("FrF2.splitplot", "FrF2.splitplot.folded","Dopt.splitlot"))
             printBy(cbind(x), names(di$factor.names)[1:di$nfac.WP], ...)
          ## see above (with show.order)
          #if (di$type == "crossed" | length(grep("param",di$type)) > 0)
          #   printBy(cbind(x), 
          #               names(di$factor.names)[1:sum(di$cross.nfactors[-length(di$cross.nfactors)])])
       }}
   cat("class=design, type=", di$type,"\n") 
   if (show.order) 
       cat("NOTE: columns run.no and run.no.std.rp are annotation, not part of the data frame",fill=TRUE)
   if (length(grep("param",di$type))>0 & length(grep("wide",di$type))>0 ){
       cat("Outer array:\n")
       print(di$outer, ...)
      }
    }
    }
}

printBy <- function(data, byvars, ...){
       ### structured printing
       ### currently separates lines by Variable names
       ### would prefer separation by blank line 
       zaehl <- 0
       zeil <- 0
       while (zeil < nrow(data)){
           zaehl <- zaehl + 1
           current <- data[zeil + 1, byvars]
           curz <- zaehl
           aus <- data[-(1:nrow(data)),]
           while (zaehl == curz & zeil < nrow(data)){
             zeil <- zeil + 1
             if (all(data[zeil,byvars] == current))
                 aus <- rbind(aus, data[zeil,])
             else {
                 print(aus, ...)
                 zaehl <- zaehl + 1
                 zeil <- zeil - 1
             }
             if (zeil == nrow(data))
                 print(aus, ...)
           }
       }
}


## brief summary without printout
## long summary with printout
summary.design <- function(object, brief = NULL, quote=FALSE, ...){
##summary.design <- function(object, ...){
   di <- design.info(object)
   
   ## dataframe summary for class design objects from package conf.design
   if (is.null(di)) return(summary.data.frame(object=object, quote=quote, ...))
   
   if (is.null(brief)) 
       if (nrow(object) <= 40 & ncol(object)<=12) brief <- FALSE else brief <- TRUE
   if (is.language(di$creator)){ 
       cat("Call:\n")
       print(di$creator, quote=quote, ...)
       cat("\n")
       }
   else {if (length(class(di$creator))>1)
       cat("design was generated with RcmdrPlugin.DoE\n\n")
       else {
           cat("Multi-step-call:\n")
           print(di$creator, quote=quote, ...)
           cat("\n")}
       }
       cat("Experimental design of type ", di$type,"\n")
       cat(di$nruns, " runs\n")
## handle blocks from ccd differently
## report varying block sizes, if applicable
   blocks <- di$nblocks
   if (is.null(blocks)) blocks <- 1
       if (blocks > 1){
          if (length(grep("ccd",di$type))>0) 
               cat("blocked design with ", blocks, " cube blocks and one star block\n")
          else
              cat("blocked design with ", blocks, " blocks")
          if (!all(di$blocksize==di$blocksize[1])){
              cat("\nVarying block sizes: \n")
              print(di$blocksize)}
              else cat(" of size ", di$blocksize, "\n")
          if (!length(grep("Dopt",di$type))>0){
          if (di$bbreps>1)
             cat("each type of block independently conducted ", di$bbreps, " times\n")
          if (di$wbreps > 1 & !di$repeat.only)
             cat("each run within each block independently conducted ", di$wbreps, " times\n")
          if (di$wbreps > 1 & di$repeat.only)
             cat("each run measured ", di$wbreps, " times (no proper replication)\n")
             }
          if (di$type=="full factorial.blocked"){
             hilf <- factorize(di$nlevels)
             names(hilf) <- Letters[1:di$nfactors]
             hilf <- unlist(hilf)
             if (is.null(colnames(di$block.gen)))
                colnames(di$block.gen) <- names(hilf)
             hilf.primes <- apply(di$block.gen, 1, 
                 function(obj) unique(hilf[!obj==0]))
             for (p in unique(hilf.primes)){
                   chilf <- conf.set(di$block.gen[which(hilf.primes==p),,drop=FALSE], p)
                   cat("\nConfounding of ", p, "-level pseudo-factors with blocks",
                   "\n(each row gives one independent confounded effect):\n")
                   print(chilf)
                   cat("\n")
                 }
          }
       }
    else if (di$replications>1)
      if (di$repeat.only)
         cat(di$replications, " measurements per run (not proper replications)\n")
      else
         cat("each run independently conducted ", di$replications, " times\n")
    
## add white space
    cat("\n")

#   nlevels <- di$nlevels
#   if (is.null(nlevels))
#      nlevels <- sapply(di$factor.names, "length")
#   names(nlevels) <- names(di$factor.names)

#   if (length(unique(nlevels))==1) message(di$nfactors, " factors with ", unique(nlevels), " levels each")
#   else {message(di$nfactors, " factors") 
#         message("numbers of levels:")
#         print(nlevels)
#   }
   pfn <- di$factor.names
   lfn <- max(sapply(pfn, "length"))
   pfn <- lapply(pfn, function(obj) if (length(obj)==lfn) obj else c(obj,rep("",lfn-length(obj))))
   pfn <- as.data.frame(pfn)
   if (all(di$quantitative)){
      if (!"ccd" %in% di$type)
          cat("Factor settings (scale ends):\n")
      else cat("Factor settings (cube):\n")
      }
   else
      cat("Factor settings:\n")
   print(pfn, quote=quote, ...)

   if ("ccd" %in% di$type){
      cat("\nNumbers of cube and star points: \n") 
      print(c(Cube=di$ncube, Star=di$nstar))
      cat("Numbers of center points: \n") 
      print(c(Cube=di$ncenter[1], Star=di$ncenter[2]))
   }

   if (length(grep("Dopt",di$type))>0 | length(grep("lhs",di$type))>0)
      if (!is.null(di$optimality.criteria)){
        cat("\nOptimality criteria:\n ") 
        print(unlist(di$optimality.criteria))
        }

   if (!is.null(response.names(object))){
       cat("\nResponses:\n")
       if (is.null(di$responselist)) print(response.names(object), quote=quote, ...)
       else print(di$responselist, quote=quote, ...)
   } 
   if (length(grep("param",design.info(object)$type))>0 & length(grep("wide",design.info(object)$type))>0 ){
      cat("\nOuter array:\n")
      print(design.info(object)$outer, quote=quote, ...)
      }
   ## alias information for FrF2 designs
   if (substr(di$type,1,4)=="FrF2"){
      cat("\nDesign generating information:\n")
      print(list(legend=di$aliased$legend), quote=quote, ...)
      ### show generator information only if valid, 
      ### i.e. if design was generated with FrF2.version at least 1.1 or 
      ### if it is not a blocked or splitplot design
      ###     other designs should not be problematic
      ###     (blocked designs without blockpick.big should also work, but ...)
      neuver <- FALSE
      if (!is.null(di$FrF2.version) & length(di$FrF2.version)==1)
         if (compareVersion(di$FrF2.version, "1.1") >= 0) neuver <- TRUE
      ## FrF2 version only relevant for single step FrF2 designs
      #if (!is.null(di$FrF2.version) & length(di$FrF2.version) > 1)
      #   if (all(sapply(di$FrF2.version, "compareVersion", "1.1") >= 0)) neuver <- TRUE

      if ((neuver | !(length(grep("blocked",di$type)) > 0 | length(grep("splitplot",di$type)) > 0)) & 
             !(length(grep("param",di$type)) > 0 | length(grep("folded", di$type))>0) )
          print(generators(object), quote=quote, ...)
          

      ## only the legend entry          
      if (length(di$aliased) == 1)
         cat("\nno aliasing of main effects or 2fis among experimental factors\n", fill=TRUE)
      else{
        ## several entries but no aliasing
        if (all(sapply(di$aliased[-1], "length")==0))
           cat("\nno aliasing of main effects or 2fis among experimental factors\n", fill=TRUE)
           else{
             ## relevant alias entries              
              if (all(sapply(di$aliased, "length") >= 1) && length(di$aliased) > 1){
                 ## more than only the legend entry, no NULL
                 cat("\nAlias structure:\n")
                 print(di$aliased[-1], quote=quote, ...)}
                 else {
                     cat("\nAlias structure:\n")
                 if (length(di$aliased$main) > 0) 
                     print(di$aliased[2], quote=quote, ...)
                 if (length(di$aliased$fi2) > 0) 
                     print(di$aliased[3], quote=quote, ...)
                 if (length(di$aliased$fi3) > 0)
                     print(di$aliased[4], quote=quote, ...)
                      }
                 }
             }
      if (di$type=="FrF2.blocked"){
        if (length(di$aliased.with.blocks) > 0){ 
           cat("Aliased with block main effects:\n")
           print(di$aliased.with.blocks, quote=quote, ...)
         }
         else cat("no main effects or 2fis aliased with blocks\n")
       }
   }
   ## what for pb and oa.design?
   if (substr(di$type,1,4)=="oa"){
       cat("Generating Orthogonal Array:\n")
         print(di$generating.oa, quote=quote, ...)
       cat("Selected Columns:\n")
         print(di$selected.columns,...)
       if (di$nfactors <= 15){
          cat("Numbers of generalized words of lengths 3 and 4:\n")
          print(c("3"=length3(object),"4"=length4(object)))}
       else if (di$nfactors <= 30)
          cat("Number of generalized words of length 3: ", length3(object),"\n")
       }
   ## nothing for pb or full factorials

   ## what for rsm designs? 
   nWPs <- di$nWPs
   if (is.null(nWPs) | length(nWPs)==0) nWPs <- 1
        ### nWPs = numeric(0) for folded designs; why?
   if (nWPs > 1){ 
          cat("\nsplit-plot design: ", nWPs, " whole plots\n")
          cat("                 : first ", di$nfac.WP, " factors are whole plot factors\n")
          }
   if (!brief){ 
      cat("\nThe design itself:\n")
      print(object, quote=quote, ...)
   }
}

