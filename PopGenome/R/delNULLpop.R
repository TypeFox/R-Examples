delNULLpop <- function(Populations){


 npops        <- length(Populations)
 PPopulations <- list()
 popmissing   <- NULL

 yy <- 1
 for(xx in 1:npops){
   if(length(Populations[[xx]])>0){
   PPopulations[[yy]] <- Populations[[xx]]
   yy <- yy + 1
   }else{
   popmissing <- c(popmissing,xx)
   }
 }

 Populations <- PPopulations
 if(length(popmissing)==0){
   popmissing <- integer(0)
 }
 
return(list(Populations=Populations,popmissing=popmissing))
}