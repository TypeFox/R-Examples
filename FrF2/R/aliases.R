`aliases` <-
function(fit, code=FALSE, condense=FALSE){
  if (is.null(stats::alias(fit)$Complete)){
  ## there is no complete aliasing whatsoever
  ## in the linear model
      aus <- list(legend=NULL, aliases=NULL)
      fit <- remodel(fit)$model
      ## there is partial aliasing 
      ## (even with -1 and 1 effects coding as done by remodel)
      if (!is.null(stats::alias(fit, partial=TRUE)$Partial))
         stop(paste("This model contains partially aliased effects.", "\n", 
            "You can look at the alias structure with alias(fit, partial=TRUE)."))
      }
  else {
      ## prevents non-admitted designs
      if (!check(fit)) 
        stop(paste("This routine is applicable for 2-level factorial designs", 
             "with not or fully aliased effects only.")) 
      
      ## create list of confounding structures
      ## using built-in function alias
      fit <- remodel(fit)$model  ## recoded factors to -1 and 1
      mm <- model.matrix(fit)
      al <- round(stats::alias(fit)$Complete,5)
          ## columns contain the "master" effects 
          ## that are not confounded with each other
      if (colnames(al)[1]=="(Intercept)") 
           al <- al[,2:ncol(al)]
      confounded <- as.list(colnames(al))
      for (i in 1:length(confounded))
         if (length(which(abs(al[,i])==1)) > 0) {
              ## in case of proper 2-level fractional factorials,
              ##    at most one entry of each row is non-zero;
              ## rows of al can thus be uniquely assigned 
              ##    to one list element of confounded
              posi <- which(abs(al[,i])==1)
              addnam <- rownames(al)[posi]
              addsign <- gsub("1","",as.character(sign(al[posi,i])))
              addnam <- paste(addsign,addnam,sep="")
              confounded[[i]] <- c(confounded[[i]],addnam)
             }
      
      ## replace variable names with codes (I and i not admitted) for brevity, 
      ##        if requested
      ## and add legend to output object
      legend <- NULL
      if (code){
        terms1 <- names(coef(fit))
        if (terms1[1]=="(Intercept)") terms1 <- terms1[-1]
        terms1 <- terms1[attr(terms(fit),"order")==1]
        faclet <- c(LETTERS[-9],letters[-9])
        codes <- faclet[1:length(terms1)]
        legend <- paste(codes,terms1,sep="=")
        for (i in 1:length(terms1)){
          for (j in 1:length(confounded)){
           for (k in 1:length(confounded[[j]])){
              if (confounded[[j]][k] == terms1[i]) 
                confounded[[j]][k] <- codes[i]
              else if (length( grep(paste(terms1[i],":",sep=""), 
                         confounded[[j]][k] )) > 0)
                confounded[[j]][k] <- gsub(paste(terms1[i],":",sep=""), 
                     paste(codes[i],":",sep=""), confounded[[j]][k] )
              else if (length( grep(paste(":",terms1[i],sep=""), 
                         confounded[[j]][k] )) > 0)
                confounded[[j]][k] <- gsub(paste(":",terms1[i],sep=""), 
                         paste(":",codes[i],sep=""), confounded[[j]][k] )
           }
          }
        }
      }  ## end of recoding
      if (condense){
      lang <- sapply(confounded,"length")
      lang <- which(lang>1)
      grad <- lapply(confounded, function(obj) unlist(strsplit(obj[1],":")))
      grad <- sapply(grad, "length")
      aliased <- list(main=NULL, fi2=NULL, fi3=NULL)
      m <- intersect(which(grad==1), lang)
      if (length(m) > 0)
         aliased$main <- sapply(confounded[m],"paste",collapse=" = ")
      m <- intersect(which(grad==2), lang)
      if (length(m) > 0)
         aliased$fi2 <- sapply(confounded[m],"paste",collapse=" = ")
      m <- intersect(which(grad==3), lang)
      if (length(m) > 0)
         aliased$fi3 <- sapply(confounded[m],"paste",collapse=" = ")
      
      aus <- c(legend=list(legend),aliased)
      }
      else aus <- list(legend = legend, aliases = confounded)
      }  ## end of long else
class(aus) <- c("aliases", class(aus))
aus
}

