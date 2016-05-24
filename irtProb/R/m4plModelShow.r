`m4plModelShow`<-
function(x, ...) {
 res <- NULL
 ID  <- 1:dim(x)[1]
 #print("res")
 res <- data.frame(ID=NA,MODEL=NA,LL=NA,AIC=NA,BIC=NA,T=0,SeT=NA,S=0,SeS=NA,C=0,SeC=NA,D=0,SeD=NA)

 #print("T")
 temp1 <- m4plPersonParameters(x=x, model="T", more=TRUE, ...)
 temp2 <- m4plMoreSummary(temp1, out="result")
 temp3 <- data.frame(ID=ID, MODEL="T",
                     temp2[c("LL","AIC","BIC","T","SeT")],
                     S=0,SeS=NA,C=0,SeC=NA,D=0,SeD=NA)
 res   <- rbind(res,temp3)[-1,]

 #print("TS")
 temp1 <- m4plPersonParameters(x=x, model="S", more=TRUE, ...)
 temp2 <- m4plMoreSummary(temp1, out="result")
 temp3 <- data.frame(ID=ID, MODEL="TS",
                     temp2[c("LL","AIC","BIC","T","SeT","S","SeS")],
                     C=0,SeC=NA,D=0,SeD=NA)
 res   <- rbind(res,temp3)
 
 #print("TC")
 temp1 <- m4plPersonParameters(x=x, model="C", more=TRUE, ...)
 temp2 <- m4plMoreSummary(temp1, out="result")
 temp3 <- data.frame(ID=ID, MODEL="TC",
                     temp2[c("LL","AIC","BIC","T","SeT")],
                     S=0,SeS=NA,
                     temp2[c("C","SeC")],
                     D=0,SeD=NA)
 res   <- rbind(res,temp3)

 #print("TD")
 temp1 <- m4plPersonParameters(x=x, model="D", more=TRUE, ...)
 temp2 <- m4plMoreSummary(temp1, out="result")
 temp3 <- data.frame(ID=ID, MODEL="TD",
                     temp2[c("LL","AIC","BIC","T","SeT")],
                     S=0,SeS=NA,C=0,SeC=NA,
                     temp2[c("D","SeD")])
 res   <- rbind(res,temp3)

 #print("TSC")
 temp1 <- m4plPersonParameters(x=x, model="SC", more=TRUE, ...)
 temp2 <- m4plMoreSummary(temp1, out="result")
 temp3 <- data.frame(ID=ID, MODEL="TSC",
                     temp2[c("LL","AIC","BIC","T","SeT","S","SeS","C","SeC")],
                     D=0,SeD=NA)
 res   <- rbind(res,temp3)

 #print("TSD")
 temp1 <- m4plPersonParameters(x=x, model="SD", more=TRUE, ...)
 temp2 <- m4plMoreSummary(temp1, out="result")
 temp3 <- data.frame(ID=ID, MODEL="TSD",
                     temp2[c("LL","AIC","BIC","T","SeT","S","SeS")],
                     C=0,SeC=NA,
                     temp2[c("D","SeD")])
 res   <- rbind(res,temp3)

 #print("TCD")
 temp1 <- m4plPersonParameters(x=x, model="CD", more=TRUE, ...)
 temp2 <- m4plMoreSummary(temp1, out="result")
 temp3 <- data.frame(ID=ID, MODEL="TCD",
                     temp2[c("LL","AIC","BIC","T","SeT")],
                     S=0,SeS=NA,
                     temp2[c("C","SeC","D","SeD")])
 res   <- rbind(res,temp3)

 #print("TSCD")
 temp1      <- m4plPersonParameters(x=x, model="SCD", more=TRUE, ...)
 temp2      <- m4plMoreSummary(temp1, out="result")
 temp3      <- data.frame(ID=ID, MODEL="TSCD",
                          temp2[c("LL","AIC","BIC","T","SeT","S","SeS","C","SeC",
                          "D","SeD")])
 res        <- rbind(res,temp3)
 class(res) <- c("modelShow","data.frame")
 return(res)
 }

