# DIF: COMPARING DIF STATISTICS

genDichoDif <- function (Data, group, focal.names, method, anchor=NULL,match="score",type= "both", criterion="LRT",
    alpha = 0.05, model = "2PL", c = NULL, engine="ltm", discr=1,irtParam = NULL, 
    nrFocal = 2, same.scale = TRUE, purify = FALSE, nrIter = 10, save.output=FALSE, output=c("out","default"))  
{
internalGenDicho<-function(){
    mets <- c("GMH", "genLogistic", "genLord")
    prov.met <- rep(0, length(method))
    for (i in 1:length(method)) {
        if (sum(method[i] == mets) == 1) 
            prov.met[i] <- 1
    }
    if (min(prov.met) == 0) {
        ind <- min((1:length(method))[prov.met == 0])
        RES <- list(NULL, method[ind])
        class(RES) <- "genDichoDif"
        return(RES)
    }
    else {
        if (length(method) == 1) 
            return(selectGenDif(Data = Data, group = group, focal.names = focal.names, match=match,type=type, criterion=criterion,
                method = method, anchor=anchor,alpha = alpha, model = model, c = c, engine=engine,discr=discr,irtParam = irtParam, 
                same.scale = same.scale, purify = purify, nrIter = nrIter,save.output=save.output,output=output))
        else {
            mat <- iters <- conv <- anchor.names<-NULL
            for (met in 1:length(method)) {
                prov <- selectGenDif(Data = Data, group = group, 
                  focal.names = focal.names, match=match,type=type,criterion=criterion, method = method[met], 
                  anchor=anchor,alpha = alpha, model = model, c = c, engine=engine,discr=discr, irtParam = irtParam, 
                  same.scale = same.scale, purify = purify, nrIter = nrIter)
anchor.names<-prov$anchor.names
                mat <- cbind(mat, rep("NoDIF", length(prov[[1]])))
                if (!is.character(prov$DIFitems)) 
                  mat[prov$DIFitems, met] <- "DIF"
                rname <- prov$names
                if (purify) {
                  iters <- c(iters, prov$nrPur)
                  conv <- c(conv, prov$convergence)
                }
            }
            method2 <- method
            method2[method == "GMH"] <- "M.-H."
            method2[method == "genLogistic"] <- "Logistic"
            method2[method == "genLord"] <- "Lord"
            colnames(mat) <- method2
            if (!is.null(rname)) rownames(mat) <- rname
            else {
                rname <- NULL
                for (i in 1:nrow(mat)) rname <- c(rname, paste("Item", 
                  i, sep = ""))
                rownames(mat) <- rname
            }
            RES <- list(DIF = mat, alpha = alpha, method = method,
                type = type, criterion=criterion, model = model, c = c, engine=engine,discr=discr,irtParam = irtParam, 
		    same.scale = same.scale, purification = purify, nrPur = iters, 
                convergence = conv, anchor.names=anchor.names,save.output=save.output,output=output)
            class(RES) <- "genDichoDif"
            return(RES)
        }
    }
}
resToReturn<-internalGenDicho()
if (save.output){
if (output[2]=="default") wd<-paste(getwd(),"/",sep="")
else wd<-output[2]
fileName<-paste(wd,output[1],".txt",sep="")
capture.output(resToReturn,file=fileName)
}
return(resToReturn)

}



# METHODS
print.genDichoDif<-function(x,...){
res <- x
if (is.null(res[[1]])) cat("Error: '",res[[2]],"' is not a correct method!","\n","\n",sep="")
else{
cat("Comparison of DIF detection among multiple groups, using",ncol(res$DIF),"methods","\n","\n")
methods<-colnames(res$DIF)
methods2<-methods
methods2[methods=="M.-H."]<-"Mantel-Haenszel"
methods2[methods=="Logistic"]<-"Logistic regression"
methods2[methods=="Lord"]<-"Lord's chi-square test"
met<-methods2[1]
for (i in 2:min(c(2,length(methods)))) met<-paste(met,methods2[i],sep=", ")
if (length(methods)<=2) cat("Methods used:",met,"\n")
else cat("Generalized methods used: ",met,",","\n",sep="")
if (length(methods)>2){
met2<-methods2[3]
if (length(methods)>3){
for (i in 4:length(methods)) met2<-paste(met2,methods2[i],sep=", ")
}
cat(met2,"\n")
}
cat("\n")
if (is.null(res$anchor.names)) {
itk<-1:nrow(res$DIF)
cat("No set of anchor items was provided", "\n", "\n")
}
else {
if (is.numeric(res$anchor.names)) {
itk<-res$anchor.names
itk.names<-rownames(res$DIF)[itk]
}
else{
itk<-NULL
for (tt in 1:length(res$anchor.names)) itk[tt]<-(1:nrow(res$DIF))[rownames(res$DIF)==res$anchor.names[tt]]
itk.names<-res$anchor.names
}
cat("Anchor items (provided by the user):", "\n")
mm <- cbind(itk.names)
rownames(mm) <- rep("", nrow(mm))
colnames(mm) <- ""
print(mm, quote = FALSE)
cat("\n", "\n")
}
cat("Parameters:","\n")
cat("Significance level: ",res$alpha,"\n",sep="")
if (sum(methods=="Logistic")==1){
if (res$type=="both") cat("DIF effects tested by generalized logistic regression: ",
				   res$type," effects","\n",sep="")
else cat("DIF effects tested by generalized logistic regression: ",res$type,
         " DIF effect","\n",sep="")
cat("DIF flagging criterion:",ifelse(res$criterion=="Wald","Wald test","Likelihood ratio test"),"\n")
}
if (sum(methods=="Lord") ==1) {
cat("Item response model:",res$model,"\n")
if (res$model!="1PL" | res$engine=="ltm") cat("Engine 'ltm' for item parameter estimation","\n")
else cat("Engine 'lme4' for item parameter estimation","\n")
if (res$model=="1PL" & res$engine=="ltm") {
if (is.null(res$discr)) cat("Common discrimination parameter: estimated from 'ltm'","\n")
else cat("Common discrimination parameter: fixed to ",res$discr,"\n",sep="")
}
if (!is.null(res$c)){
if (length(res$c)==1) cat("Common pseudo-guessing value: ",res$c,"\n",sep="")
else {
pg<-cbind(res$c)
rownames(pg)<-res$names
colnames(pg)<-"c"
cat("Common pseudo-guessing values:","\n")
print(pg)
cat("\n")
}
}
}
if (res$purification & is.null(res$anchor.names)) {
cat("Item purification: Yes","\n","\n")
cat("Item purification results:","\n","\n")
co<-rep("Yes",length(res$convergence))
co[!res$convergence]<-"No"
resConv<-data.frame(rbind(co,res$nrPur))
colnames(resConv)<-colnames(res$DIF)
rownames(resConv)<-c("Convergence","Iterations")
print(format(resConv,justify="centre"))
cat("\n")
}
else cat("Item purification: No","\n","\n")
cat("Comparison of DIF detection results:","\n","\n")
nr<-NULL
for (i in 1:nrow(res$DIF)) nr[i]<-paste(length(res$DIF[i,][res$DIF[i,]=="DIF"]),"/",ncol(res$DIF),sep="")
MAT<-cbind(res$DIF,nr)
colnames(MAT)[ncol(MAT)]<-"#DIF"
print(format(MAT,justify="centre"),quote=FALSE)
}
if (!x$save.output) cat("\n","Output was not captured!","\n")
else {
if (x$output[2]=="default") wd<-paste(getwd(),"/",sep="")
else wd<-x$output[2]
fileName<-paste(wd,x$output[1],".txt",sep="")
cat("\n","Output was captured and saved into file","\n"," '",fileName,"'","\n","\n",sep="")
}
}

