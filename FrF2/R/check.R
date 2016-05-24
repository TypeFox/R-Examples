`check` <-
function(obj){
     check<-TRUE
     ## check whether completely aliased only
     ## works for numeric -1 and 1 coded data only
     mm <- model.matrix(remodel(obj)$model)
     if (any(colnames(mm)=="(Intercept)")) 
             mm <- mm[,-which(colnames(mm)=="(Intercept)")]
     if (!all(round(t(mm)%*%mm, 10) %in% c(0,nrow(mm),-nrow(mm)))) check <- FALSE 
     check
   }

