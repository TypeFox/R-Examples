

##########################################################################
# model fit for polytomous item responses
modelfit.cor.poly <- function( data , probs , theta.k , f.qk.yi){
    # create input for tam.modelfit
    tamobj <- list()
    resp.ind <- 1 - is.na(data)
    resp <- data
    resp[ is.na(data) ] <- 0
    tamobj$resp <- as.matrix(resp)
    tamobj$resp.ind <- resp.ind
    tamobj$rprobs <- probs
    tamobj$theta <- vec2mat.sirt(theta.k)
    tamobj$hwt <- f.qk.yi
    # apply tam.modelfit
    tmod <- TAM::tam.modelfit( tamobj)
	# tmod <- tam.modelfit( tamobj)
    return(tmod)
            }
##########################################################################