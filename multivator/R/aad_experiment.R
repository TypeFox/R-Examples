setClass("experiment",
         representation = representation(
           mm  = "mdm",
           obs = "numeric")
         )
"get_mdm" <- function(x){x@mm}
"get_obs" <- function(x){x@obs}

# there are no occurrences of '@' below this line.

".experiment_valid" <- function(object){
    if(nrow(get_mdm(object)) != length(get_obs(object))){
    return("length mismatch")
  } else {
    return(TRUE)
  }
}

setValidity("experiment" , .experiment_valid)

"experiment" <- function(mm,obs){
  new("experiment",mm=mm,obs=obs)
}
  
".experiment_print" <- function(x){
  as.data.frame(x)
}
    
"print.experiment" <- function(x, ...){
  jj <- .experiment_print(x, ...)
  print(jj)
  return(invisible(jj))
}

setMethod("show", "experiment", function(object){print.experiment(object)})
setMethod("levels","experiment",function(x){levels(get_mdm(x))})

setAs("experiment","mhp",function(from){
  as.mhp(get_mdm(from))
} )

setMethod("as.mhp","experiment",function(x){as(x,"mhp")})

setAs("experiment","data.frame",function(from){
  data.frame(as.data.frame(get_mdm(from)),obs=get_obs(from))
} )

setMethod("as.data.frame",signature="experiment",
          function(x){as(x,"data.frame") })


setMethod("head",signature="experiment",function(x,n=6,...){
  experiment(head(get_mdm(x),n=n,...), head(get_obs(x),n=n,...))
} )

setMethod("head",signature="experiment",function(x,n=6,...){
  experiment(tail(get_mdm(x),n=n,...), tail(get_obs(x),n=n,...))
} )


setMethod("[", signature(x="experiment"),
          function(x,i,j,drop=FALSE){
            experiment(
                       mm  = get_mdm(x)[i,j,drop=drop],
                       obs = get_obs(x)[i]
                       )
          } )
