### Compute predicted values witin a BJ iteration                                    
predval <- function(learner, twin, dat1.glm, b, k, x, s, mselect){
    beta0bj <- betabj <- bdiff <- NULL
    if(learner=="linear.regression"){
        if(twin){
            beta0bj <- attr(coef(dat1.glm), "offset2int")
            betabj <- coef(dat1.glm)     ###for bst
        }
	    else{
                beta0bj <- coef(dat1.glm)[1] + dat1.glm$offset
		betabj <- coef(dat1.glm, which = 1:length(variable.names(dat1.glm)))[-1] 
            }
	Fboost <- predict(dat1.glm)
        b[k,] <- betabj
        if(k == 1) 
            bdiff <- 1000
        else
            bdiff <- sum((b[k,] - b[k-1,])^2)
    }
    else if(learner=="enet"){
        Fboost <- predict(dat1.glm, x, type="fit", s=s, mode="fraction")$fit
    }
    else if(learner %in% c("enet2", "mnet", "snet")){
        Fboost <- predict(dat1.glm, newx=x, type="response", which=mselect)
    } 
    else if(learner %in%c("pspline", "mars")){
        Fboost <- predict(dat1.glm)
    }
    else if(learner=="tree"){
        if(!twin){            
            Fboost <- dat1.glm$fit   #predicted values on training data
        }
	  else{ 
	      Fboost <- predict(dat1.glm)
          }
    }
    RET <- list(beta0bj = beta0bj, betabj = betabj, Fboost = Fboost, b = b, bdiff=bdiff)
    RET	
}
