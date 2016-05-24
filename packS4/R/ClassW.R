#################################################################
###                         ClassW.R                          ###
#################################################################

setClass(
    "ClassW",
    representation(w1="character"),
    contains="ClassU"
)
classW <- function(u1,u2,w1){new("ClassW",u1=u1,u2=u2,w1=w1)}
setMethod("publicA","ClassW",function(object){sqrt(object@u2^5)})
setMethod("plot","ClassW",
    function(x,y){barplot(c(x@u1,x@u2),main=x@w1)}
)
setMethod("[","ClassW",
    function(x,i,j,drop){
        switch(EXP=i,
            "u1"={return(x@u1)},
            "u2"={return(x@u2)},
            "w1"={return(x@w1)}
        )
    }
)

setReplaceMethod("[","ClassW",
    function(x,i,j,value){
        switch(EXP=i,
            "u1"={x@u1<-value},
            "u2"={x@u2<-value},
            "w1"={x@w1<-value}
        )
        validObject(x)
        return(x)
    }
)
