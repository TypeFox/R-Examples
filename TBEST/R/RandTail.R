RandTail <-
#function(mydata, mystat, mymethod, mymetric, rand.fun = c("shuffle.column", 
function(mydata, myinput, mystat, mymethod, mymetric, rand.fun = c("shuffle.column",
    "shuffle.block", "define.function"), by.block = NA, metric.args = list(), 
    rand.args = list()) 
{
    #myinput <- mydata$myinput
    ntest <- mydata$nperm
    indextable <- TreeStat(myinput, mystat = mystat, method = mymethod, 
        metric = mymetric, metric.args = metric.args)
    statnames <- mystat   
    nullstat <-vector("list",length(statnames))
    names(nullstat) <- statnames 
    for (i in 1:ntest) {
        if (rand.fun == "shuffle.column") { 
            myrdata <- apply(myinput, 2, sample)
        }else if (rand.fun == "shuffle.block") {
	    if (is.na(by.block[1]))stop("by.block needs to be specified")
            myrdata <- t(myinput)
            myrlist <- by(myrdata, by.block, FUN = byfactor)
            for (j in 1:length(myrlist)) {
                if (j == 1) {
                  myrdata <- myrlist[[j]]
                }
                else {
                  myrdata <- rbind(myrdata, myrlist[[j]])
                }
            }
            myrdata <- t(myrdata)
        }else if (rand.fun == "define.function") {
            define.function <- match.fun(define.function)
            myrand.args <- vector("list", length(rand.args) + 
                1)
            myrand.args[[1]] <- myinput
            if (length(myrand.args) > 1) {
                myrand.args[2:length(myrand.args)] <- rand.args
            }
            myrdata <- do.call(define.function, myrand.args)
        }
        rindextable <- TreeStat(myrdata, mystat = mystat, method = mymethod, 
            metric = mymetric, metric.args = metric.args)
        size <-rindextable[,"clustersize"]
	if(any(statnames!="slb")){	
	for (statname in statnames) {
		rstat <- rindextable[,statname]
		statmax <- max(rstat)
		randomX <- sort(size + 0.5 * rstat/statmax)
		rmatch <- bestmatch(rsize=sort(size),size = indextable[,
            		"clustersize"])
		rmatchl <- bestmatchl(rsize=sort(size),size = indextable[,
                        "clustersize"])
		data <- 2 * statmax * (randomX - sort(size))[rmatch]
        	datal <- 2 * statmax * (randomX - sort(size))[rmatchl]
		mydata <- pmax(data, datal)
		if(i==1){nullstat[[statname]]<-mydata}
		else{nullstat[[statname]]<-c(nullstat[[statname]],mydata)}
    	}
	}else if(statnames=="slb"){
		if(i==1){nullstat[[statnames]]<-rindextable[nrow(rindextable),"slb"]}
		else{nullstat[[statnames]]<-c(nullstat[[statnames]],rindextable[nrow(rindextable),"slb"])}
	}
	}
    return(nullstat)
}
