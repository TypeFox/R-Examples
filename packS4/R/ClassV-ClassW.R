##################################################
###              ClassV-ClassW.R
setMethod("publicB",c("ClassV","ClassW"),
    function(objectV,objectW){objectV@v1*objectW@u1})
