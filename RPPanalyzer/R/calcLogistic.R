`calcLogistic` <-
function(x,sample.id=c("sample","sample.n"),dilution="dilution"
        ,xVal=NULL, plot=F, detectionLimit=F) {

    ## generate column to identify individual samples
    xi <- create.ID.col(x,sample.id=sample.id)

    ## identify samples in a vector
    id <- (unique(xi[[4]][,"identifier"]))

    ## generate empty matrix to store concentration vals
    vals <- matrix(NA, nrow=length(id), ncol=ncol(xi[[1]])
            ,dimnames=list(id,colnames(xi[[1]])))

    ## for loop over the samples
    for (i in seq(along=id)){

        x.lines <- which(as.character(xi[[4]][,"identifier"])==id[i])

        ## for loop over the analyzed targets
        for (j in 1:ncol(xi[[1]])){

            xvals <- xi[[4]][x.lines,dilution]
            yvals <- xi[[1]][x.lines,j]

            params <- curveFitSigmoid(xvals,yvals)

            vals[i,j] <- params$val
        }
    }

    tempdat <- pick.high.conc(xi,highest=dilution)

    tempid <-unique(as.character(tempdat[[4]][,"identifier"]))

    order.lines <- match(tempid,rownames(vals))
    vals <- vals[order.lines,]

    sampledescription <- tempdat[[4]][match(tempid,tempdat[[4]][,"identifier"]),]

    data.list <- list(expression=vals
            , dummy=vals
            , arraydescription=xi[[3]]
            , sampledescription=sampledescription)
    return(data.list)
}

