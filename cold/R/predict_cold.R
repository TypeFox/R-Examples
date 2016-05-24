
setClass("cold", representation(coefficients = "matrix", se = "matrix", covariance = "matrix", correlation="matrix", 
		log.likelihood="numeric", message ="integer",n.cases="numeric", ni.cases="numeric", aic="numeric",  
		 Fitted="numeric", Fitted.av="numeric", Time="numeric", model.matrix= "matrix", y.matrix="matrix",
		subset.data="data.frame",final.data="data.frame", y.av="numeric", f.value="factor", data.id="numeric",call="language"))


setMethod(f="predict",signature=c(object="cold"), 

function (object) 
{
    
        mu <- object@Fitted

        res<-log(mu)
res

}
)