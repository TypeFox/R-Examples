print.catlg <- function(x, name="all", nruns="all", nfactors="all", 
                        res.min=3, MaxC2=FALSE, show=10, 
                        gen.letters=FALSE, show.alias=FALSE, ...){
     if (name[1]=="all" & nruns[1]=="all" & nfactors[1]=="all"){
       ## no filtering needed on x, apart from resolution filtering
       if (gen.letters) hilf <- names(Yates) else hilf <- 1:4095
       if (!is.numeric(res.min)) stop("res.min must be a number.")
       if (!(res.min==floor(res.min) & res.min>2)) stop("res.min must be an integer number larger than 2.")
       x <- x[res.catlg(x)>=res.min]
       if (MaxC2) x <- x[sort(nclear.2fis.catlg(x),index.return=TRUE, decreasing=TRUE)$ix]
       alias <- NULL
       for (i in 1:min(length(x),show)) {
           resp <- as.character(as.roman(x[[i]]$res))
           if (is.na(resp)) resp <- "full factorial"
           ## resp allows to extend catlg to include full factorials 
           ## with gen equal to numeric(0)
           ## warning for as.roman would have to be blocked
           ## further warnings in gen.check would have to be blocked
           ##     would there be further issues ???
           ## Sole purpose so far: also show full factorial possibilities with buttons in 
           ##     create design menu of RcmdrPlugin.DoE 
           if (show.alias) alias <- alias3fi(round(log2(x[[i]]$nruns)),x[[i]]$gen,order=2)
           cat("Design: ", names(x)[i], "\n  ",
                x[[i]]$nruns, " runs, ", 
                x[[i]]$nfac, " factors,  \n   Resolution ", resp, "\n") 
           if (!resp=="full factorial")
               cat("   Generating columns: ", hilf[x[[i]]$gen], "\n",
               "  WLP (3plus): ", x[[i]]$WLP[-(1:2)], ", ", x[[i]]$nclear.2fis, " clear 2fis")
               if (length(x[[i]]$all.2fis.clear)>0) if(!x[[i]]$all.2fis.clear[1]=="all"){
                    if (x[[i]]$nfac <= 50) cat("\n Factors with all 2fis clear: ", Letters[x[[i]]$all.2fis.clear])
                    else cat("\n Factors with all 2fis clear: ", paste("F",x[[i]]$all.2fis.clear,sep=""))}
               if (is.list(alias)){ 
                   cat("\n Alias structure: \n")
                   if (x[[i]]$res==3){ 
                         cat("   Main effects aliasing: \n ")
                         print(alias[[1]],quote=FALSE) 
                         }
                    cat("   2fi aliasing: \n ")
                   print(alias[[2]],quote=FALSE) 
               }
               cat("\n")
       }
       }
     else{
       ## filtering x, and recalling print function afterwards
       if (!name[1]=="all") {
          ## names take precedence over other specifications
          if (!is.character(name)) stop("name must be a vector of character strings naming catalogued designs.")
          if (!length(setdiff(name,names(x)))==0) stop("name must contain names of catalogued designs only.")
          x <- x[name]
       }
       else if (!(nruns[1]=="all" & nfactors[1]=="all")){
         if (!nruns[1]=="all") {
            if (!is.numeric(nruns)) stop("nruns must be a numeric vector.")
            if (!all(2^round(log2(nruns))==nruns)) 
                stop("nruns must be a vector of powers of 2.")
            x <- x[nruns.catlg(x) %in% nruns]
          }
         if (!nfactors[1]=="all"){
            if (!is.numeric(nfactors)) stop("nfactors must be a numeric vector.")
            x <- x[nfac.catlg(x) %in% nfactors]
          }
     }
     print(x, res.min=res.min, show=show, MaxC2=MaxC2, gen.letters=gen.letters, show.alias=show.alias, ...)
     }
}
