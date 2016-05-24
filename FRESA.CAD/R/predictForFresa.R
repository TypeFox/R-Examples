predictForFresa <-
function (object,testData, predictType = c("prob", "linear")) 
{

	frm <- formula(terms(object));
	classlen=length(class(object))
	
	cobj <- substr(class(object)[classlen], 1, 2);
	switch(cobj,
		co =
		{
			pobj <- object;
			switch(predictType, 
				linear = 
					{		
						out <- predict(pobj,testData,type = 'lp');
					},
				prob = 
					{
						out <- 1.0/(1.0+exp(-predict(pobj,testData,type = 'lp')));
					},
					{
						out <- predict(pobj,testData,type = 'lp');
					}
			)
		},
		fi =
		{
			cf <- object$estimations;
			pred <-.Call("predictForFresaCpp",cf,model.matrix(frm, testData),predictType,object$type);
			out <- pred$prediction;
		},
		{
			predictType <- match.arg(predictType)
			cf <- coef(object)
			s <- is.na(cf);
			if (any(s)) 
			{
				cf[s] <- 0;
			}
			
			switch(predictType, 
				linear = 
					{
					   out <- as.vector(model.matrix(frm, testData)  %*% cf);
					}, 
				prob = 
					{
					  out <- as.vector(model.matrix(frm, testData)  %*% cf);
					  out <- 1.0/(1.0+exp(-out));
					}, 
					{
					  out <- as.vector(model.matrix(frm, testData)  %*% cf);
					}
			)
		}
	)
	s <- is.na(out);
	if (any(s)) 
	{
		cat("Warning NA predictForFresa \n");
		out[s] <- 0;
	}
    return (out)
}
