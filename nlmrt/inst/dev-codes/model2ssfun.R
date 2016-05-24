model2resfun <- function(resformula, pvec, filename) {
   pnames<-names(pvec)
#   cat("pnames:")
#   print(pnames)
   if (is.null(pnames) ) stop("MUST have named parameters in pvec")
   if (is.character(resformula)){
      es<-resformula
   } else {
      tstr<-as.character(resformula) # note ordering of terms!
      es<-paste(tstr[[2]],"~",tstr[[3]],'')
   }
   xx <- all.vars(parse(text=es))
#   cat("xx:")
#   print(xx)
   rp <- match(pnames, xx) # Match names to parameters
# ?? How to ensure there are names?
   xx2 <- c(xx[rp], xx[-rp])
   xxparm<-xx[rp]
#   cat("xx2:")
#   print(xx2)
#   cat("xxparm:")
#   print(xxparm)
   pstr<-"c("
   npar<-length(xxparm)
   if (npar>0) {
      for (i in 1:npar){
         pstr<-paste(pstr,"\"",xxparm[i],"\"", sep='')
         if (i<npar) pstr<-paste(pstr,", ",sep='')
      }
   }
   pstr<-paste(pstr,")",sep='')
#   cat("pstr:")
#   print(pstr)
#   tmp<-readline("...")
   xxvars<-xx[-rp]
   nvar<-length(xxvars)
   vstr<-""
   if (nvar>0) {
      for (i in 1:nvar){
         vstr<-paste(vstr,xxvars[i]," = NULL", sep='')
         if (i<nvar) vstr<-paste(vstr,", ",sep='')
      }
   }
#   cat("vstr:")
#   print(vstr)
#   tmp<-readline("...")
   ff <- vector("list", length(xx2))
   names(ff) <- xx2
   sf<-as.character(resformula)
   parts<-strsplit(as.character(sf), "~")[[1]]
   if (length(parts)!=2) stop("Model expression is incorrect!")
   lhs<-parts[1]
   rhs<-parts[2]
#  And build the residual at the parameters
   resexp<-paste(rhs,"-",lhs, collapse=" ") # build the residuals
   fnexp<-paste("as.numeric(eval(crossprod(",resexp,")))", sep="") ##3
   pparse<-paste("for (i in 1:length(prm) ){\n",
      "joe<-paste(names(prm)[[i]],\"<-\",prm[[i]]);\n",
      " eval(parse(text=joe));\n",
      "}", sep='')
   myfstr<-paste("myfn<-function(prm, ",vstr,"){;\n",
      "if ( is.null(names(prm)) ) { warning('Parameter vector must have names'); return(NA)};\n",
      pparse,";\n ",
      fnexp,"\n }",sep='')
   write(myfstr, file=filename) # write out the file
   myfn<-source(file=filename)$value # This may be inefficient, but ...
   attr(myfn, "source-text")<-myfstr
   return(myfn)      
}


