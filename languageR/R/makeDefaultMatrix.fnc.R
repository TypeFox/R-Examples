`makeDefaultMatrix.fnc` <-
function(model, n = 100,  
  conditioningPred="", 
  conditioningValue=NULL,
  control=NA) {


   coefs = fixef(model)
   ncoefs = length(coefs)
   X = model@pp$X

   nams =  strsplit(names(coefs), ":")

   if (is.character(conditioningValue)) {
     condName = paste(conditioningPred, conditioningValue, sep="")
   } else {
     condName = conditioningPred
   }
   
   m = matrix(0, n, ncoefs)
   #X = unique(model@X)
   rownames(X) = 1:nrow(X)
   for (i in 1:ncoefs) {
     if (names(coefs[i]) == "(Intercept)") {
       m[,i] = rep(1, n)
     } else {  
       v = names(table(X[,names(coefs[i])]))
       if (length(v)==2 & v[1] == "0" & v[2] == "1"){   # dummy vector
         if (condName==names(coefs)[i]) m[,i] = rep(1, length=n)
         else m[,i] = rep(0, length=n)
       } else {
         if (length(nams[[i]])==1) {  # a main effect
           if (condName==names(coefs)[i]) {
             m[,i] = rep(conditioningValue, length=n)
           } else {
             if (regexpr('^poly\\(', names(coefs[i])) > 0) {
               if (regexpr("1$", names(coefs[i])) > 0) {
                 maxval = max(X[X[,i]<median(X[,i]),][,i])
                 maxbelowmedianpos = which(X[,i] == maxval)[1]
               } 
               m[,i] = rep(X[maxbelowmedianpos,names(coefs[i])], length=n)
             } else {
               if (regexpr('^rcs\\(', names(coefs[i])) > 0) {
                 if (regexpr("[^']$", names(coefs[i]))>0) {
                   maxval = max(X[X[,i]<median(X[,i]),][,i])
                   maxbelowmedianpos = which(X[,i] == maxval)[1]
                 }
                 m[,i] = rep(X[maxbelowmedianpos,names(coefs[i])], length=n)
               } else {  
                 m[,i] = rep(median(X[,names(coefs[i])]), length=n)
               }
             }
           }
         } else {  # an interaction, to be updated later by implementInteractions.fnc()
           m[,i] = rep(0, length=n)
         }
       }
     }
   }
   colnames(m) = colnames(X)

   
   if (!is.na(control)[[1]]) {
     controlPredName = as.vector(control[[1]])
     if (!is.element(controlPredName, colnames(m))) {
       stop(paste ("the control predictor name", controlPredName, "is not a valid column name\n", sep=" "))
     } else {
       m[,controlPredName] = rep(as.vector(control[[2]]), nrow(m))
     }
   }

   return(m)
}

