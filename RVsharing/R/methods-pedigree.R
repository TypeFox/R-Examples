setMethod("RVsharing",  signature(data="pedigree", dad.id="missing",mom.id="missing",carriers="character"), function(data,carriers){

  fa.index <- ifelse(data$findex !=0, data$findex, NA)
  fa.id <- data$id[fa.index]
  ma.index <- ifelse(data$mindex !=0, data$mindex, NA)
  ma.id <- data$id[ma.index]
  
    RVsharing.fn( id = data$id, dad.id = fa.id, mom.id = ma.id, carriers=carriers )
})
setMethod("RVsharing",  signature(data="pedigree", dad.id="missing",mom.id="missing",carriers="missing"), function(data){

  fa.index <- ifelse(data$findex !=0, data$findex, NA)
  fa.id <- data$id[fa.index]
  ma.index <- ifelse(data$mindex !=0, data$mindex, NA)
  ma.id <- data$id[ma.index]
  
    RVsharing.fn( id = data$id, dad.id = fa.id, mom.id = ma.id)
})
setMethod("RVsharing",  signature(data="character",dad.id="character",mom.id="character",carriers="character"), function(data, dad.id, mom.id, carriers){
    # Converting 0s to NA  in dad.id and mom.id
    dad.id = ifelse(dad.id=="0",NA,dad.id)
    mom.id = ifelse(mom.id=="0",NA,mom.id)
    RVsharing.fn( id = data, dad.id = dad.id, mom.id = mom.id, carriers=carriers )
})
setMethod("RVsharing",  signature(data="character",dad.id="character",mom.id="character",carriers="missing"), function(data, dad.id, mom.id){
    # Converting 0s to NA  in dad.id and mom.id
    dad.id = ifelse(dad.id=="0",NA,dad.id)
    mom.id = ifelse(mom.id=="0",NA,mom.id)
    RVsharing.fn( id = data, dad.id = dad.id, mom.id = mom.id)
})
setMethod("RVsharing",  signature(data="numeric",dad.id="numeric",mom.id="numeric",carriers="numeric"), function(data, dad.id, mom.id, carriers){
    # Converting 0s to NA  in dad.id and mom.id
    dad.id = ifelse(dad.id==0,NA,dad.id)
    mom.id = ifelse(mom.id==0,NA,mom.id)
    RVsharing.fn( id = data, dad.id = dad.id, mom.id = mom.id, carriers=carriers )
})
setMethod("RVsharing",  signature(data="numeric",dad.id="numeric",mom.id="numeric",carriers="missing"), function(data, dad.id, mom.id){
    # Converting 0s to NA  in dad.id and mom.id
    dad.id = ifelse(dad.id==0,NA,dad.id)
    mom.id = ifelse(mom.id==0,NA,mom.id)
    RVsharing.fn( id = data, dad.id = dad.id, mom.id = mom.id)
})
