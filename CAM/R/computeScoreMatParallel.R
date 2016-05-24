computeScoreMatParallel <-
function(rowParents, scoreName, X, selMat, output, node2, i, parsScore, intervMat, intervData)
{
    #the i-th row of rowParents contains possible parents of node2 (we call them "parentsToCheck") 
    parentsToCheck <- rowParents[i,]
    if(output)
    {
        cat("\r compute score entry for regressing",node2,"on",parentsToCheck,"                  \r")
    }
    if(intervData)
    {
        X2 <- X[!intervMat[,node2],]
    } else
    {
        X2 <- X
    }
    
    if(!(node2 %in% parentsToCheck) && (prod(selMat[parentsToCheck,node2]) == TRUE))
    {
        if(scoreName == "SEMSEV")
        {      
            stop("This score does not work. It does not decouple.")
        } else if(scoreName == "SEMIND")
        {
            stop("NOT IMPLEMENTED YET")
        } else if(scoreName == "SEMGAM")
        {
            mod_gam <- train_gam(X2[,parentsToCheck],X2[,node2],pars=parsScore)
            score <- (-log(var(mod_gam$residuals)))
        } else if(scoreName == "SEMLIN")
        {
            mod_gam <- train_linear(X2[,parentsToCheck],X2[,node2])
            score <- (-log(var(mod_gam$residuals)))
        } else if(scoreName == "SEMGP")
        {
            mod_gp <- train_gp(X2[,parentsToCheck],X2[,node2])
            score <- (-log(var(mod_gp$residuals)))
        } else
        {
            stop("I do not know this score function.")
        }
    } else
    {
        score <- (-Inf)
    }
    return(score)
}
