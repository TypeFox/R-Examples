setClassUnion("param",c("numeric","NULL"))

setClass("scores", representation(phi = "function", Dphi = "function", param="param"))

### Define methods ###
setGeneric(name="getScores", def=function(object, x) {standardGeneric("getScores")})
setGeneric(name="getScoresDeriv", def=function(object, x) {standardGeneric("getScoresDeriv")})

setMethod("getScores", "scores",
        function(object, x) {

                if(is.null(object@param)) {
                        a<-object@phi(x)
                } else {
                        a<-object@phi(x, object@param)
                }

				ac<-a-mean(a)
				sigma<-sum(ac*ac)/(length(ac)+1)
                ac/sqrt(sigma)

        }
)

setMethod("getScoresDeriv", "scores",
        function(object, x) {
        if(is.null(object@param)) {
                aP<-object@Dphi(x)
				a<-object@phi(x)
        } else {
                aP<-object@Dphi(x, object@param)
				a<-object@phi(x,object@param)
        }

		ac<-a-mean(a)
		sigma<-sum(ac*ac)/(length(ac)+1)
        aP/sqrt(sigma)
		
}
)

