
setClass("cold", representation(coefficients = "matrix", se = "matrix", covariance = "matrix", correlation="matrix", 
		log.likelihood="numeric", message ="integer",n.cases="numeric", ni.cases="numeric", aic="numeric",  
		 Fitted="numeric", Fitted.av="numeric", Time="numeric", model.matrix= "matrix", y.matrix="matrix",
		subset.data="data.frame",final.data="data.frame", y.av="numeric", f.value="factor", data.id="numeric",call="language"))

setMethod(f="resid",signature=c(object="cold"), 

function (object, type = c( "pearson","response","null"),...) 
{
    

	type <- match.arg(type)

	
	data<-object@subset.data

	data$y<-as.vector(t(object@y.matrix))  ### IMPORTANTE
	ynew<-data$y

        mu <- object@Fitted


        switch(type, pearson = , response = )


 	if(type=="pearson"||type=="null"){
   
	res = (ynew - mu)/sqrt(mu)}
   
	else
  
        {res = (ynew - mu)}

    res
}
)

