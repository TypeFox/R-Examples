###################################################
### chunk number 1: 
###################################################

test <- function(data = c("k1", "m5"), range = 157:339){
  res <- list()
  for(i in data){
    disProgObj <- readData(i,week53to52=TRUE)
    disProgObj <- enlargeData(disProgObj)
    survResults <- algo.call(disProgObj, 
                             control = list(
                               list(funcName = "rki1", range = range),
                               list(funcName = "rki2", range = range),
                               list(funcName = "rki3", range = range),
                               list(funcName = "bayes", range = range,alpha=0.05)))
    res[[i]] <- algo.compare(survResults)
    cat("\n\n\n", i, " Res:\n")
    print(compMatrix.writeTable(res[[i]]))
  }
  sum <- algo.summary(res)
  cat("\n\nSummary:\n")
  print(compMatrix.writeTable(sum))
}


###################################################
### chunk number 2: 
###################################################

testSim <- function(p = 0.99, r = 0.01, length = 400, A = 1, alpha = 1, beta = 0,
                                phi = 0, frequency = 1, state = NULL, K, range = 200:400){

        disProgObj <- sim.pointSource(p, r, length, A, alpha, beta,
                                phi, frequency, state, K)
        survResults <- algo.call(disProgObj, control = list(list(funcName = "rki1", range = range)))
        res <- algo.compare(survResults)
        plot(survResults[[1]], "RKI 1", "Simulation")
        print(compMatrix.writeTable(res))

}




###################################################
### chunk number 3: 
###################################################

makePlot <- function(outputpath, data = "k1", method = "rki1", name, disease, range = 157:339){
        disProgObj <- readData(data,week53to52=TRUE)
        disProgObj <- enlargeData(disProgObj)
        res <- algo.call(disProgObj, control = list(list(funcName = method, range = range)))
        pdf(paste(outputpath, data, "_", method, "_plot.pdf", sep=""), width = 10)
                plot(res[[1]],name,disease)
        dev.off()
}



