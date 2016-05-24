# DIF: COMPARING DIF STATISTICS

dichoDif<-function(Data,group,focal.name,method,anchor=NULL,props=NULL,thrTID=1.5,alpha=0.05,MHstat="MHChisq",correct=TRUE,exact=FALSE,stdWeight="focal",thrSTD=0.1,BDstat="BD",
member.type="group", match="score",type="both",criterion="LRT",model="2PL",c=NULL,engine="ltm",discr=1,irtParam=NULL,same.scale=TRUE,signed=FALSE,purify=FALSE,
nrIter=10,save.output=FALSE, output=c("out","default")) 
{
internalDicho<-function(){
mets<-c("TID","MH","Std","Logistic","BD","Lord","Raju","LRT")
prov.met<-rep(0,length(method))
for (i in 1:length(method)){
if (sum(method[i]==mets)==1) prov.met[i]<-1
}
if (min(prov.met)==0){
ind<-min((1:length(method))[prov.met==0])
RES<-list(NULL,method[ind])
class(RES)<-"dichoDif"
return(RES)
}
else{
if (length(method)==1) return(selectDif(Data=Data,group=group,focal.name=focal.name,method=method,anchor=anchor,props=props,thrTID=thrTID,alpha=alpha,MHstat=MHstat,
correct=correct,exact=exact,stdWeight=stdWeight,thrSTD=thrSTD,BDstat=BDstat,member.type=member.type, match=match,type=type,criterion=criterion,model=model,c=c,engine=engine,discr=discr,irtParam=irtParam,
same.scale=same.scale,signed=signed,purify=purify,nrIter=nrIter, save.output=save.output,output=output))
else{
mat<-iters<-conv<-anchor.names<-NULL
for (met in 1:length(method)){
prov<-selectDif(Data=Data,group=group,focal.name=focal.name,method=method[met],anchor=anchor,props=props,thrTID=thrTID,alpha=alpha,MHstat=MHstat,correct=correct,exact=exact,
stdWeight=stdWeight,thrSTD=thrSTD,BDstat=BDstat,member.type=member.type, match=match,type=type,criterion=criterion,model=model,c=c,engine=engine,discr=discr,irtParam=irtParam,same.scale=same.scale,
signed=signed,purify=purify,nrIter=nrIter)
if (method[met]!="LRT") anchor.names<-prov$anchor.names
if (method[met]=="BD") mat<-cbind(mat,rep("NoDIF",nrow(prov[[1]])))
else mat<-cbind(mat,rep("NoDIF",length(prov[[1]])))
if (!is.character(prov$DIFitems)) mat[prov$DIFitems,met]<-"DIF"
rname<-prov$names
if (purify){
iters<-c(iters,prov$nrPur)
conv<-c(conv,prov$convergence)
} 
}
method2<-method
method2[method=="TID"]<-"T.I.D."
method2[method=="MH"]<-"M-H"
method2[method=="Std"]<-"Stand."
colnames(mat)<-method2
if (!is.null(rname)) rownames(mat)<-rname
else{
rname<-NULL
for (i in 1:nrow(mat)) rname<-c(rname,paste("Item",i,sep=""))
rownames(mat)<-rname
}
RES<-list(DIF=mat,props=props,thrTID=thrTID,correct=correct,exact=exact,alpha=alpha,MHstat=MHstat,stdWeight=stdWeight,thrSTD=thrSTD,
BDstat=BDstat,member.type=member.type, match=match,type=type,criterion=criterion,model=model,c=c,engine=engine,discr=discr,irtParam=irtParam,same.scale=same.scale,
signed=signed,purification=purify,nrPur=iters,convergence=conv, anchor.names=anchor.names,save.output=save.output,output=output)
class(RES)<-"dichoDif"
return(RES)}
}
}
resToReturn<-internalDicho()
if (save.output){
if (output[2]=="default") wd<-paste(getwd(),"/",sep="")
else wd<-output[2]
fileName<-paste(wd,output[1],".txt",sep="")
capture.output(resToReturn,file=fileName)
}
return(resToReturn)
}



# METHODS
print.dichoDif<-function (x, ...) 
{
    res <- x
    if (is.null(res[[1]])) 
        cat("Error: '", res[[2]], "' is not a correct method!", 
            "\n", "\n", sep = "")
    else {
        cat("Comparison of DIF detection results using", ncol(res$DIF), 
            "methods", "\n", "\n")
        methods <- colnames(res$DIF)
        methods2 <- methods
        methods2[methods == "T.I.D."] <- "Transformed item difficulties (TID)"
        methods2[methods == "M-H"] <- "Mantel-Haenszel"
        methods2[methods == "Stand."] <- "Standardization"
        methods2[methods == "Logistic"] <- "Logistic regression"
        methods2[methods == "BD"] <- "Breslow-Day"
        methods2[methods == "Raju"] <- "Raju's area"
        methods2[methods == "Lord"] <- "Lord's chi-square test"
        methods2[methods == "LRT"] <- "Likelihood ratio test"
        cat("Methods used:", "\n")
        for (i in 1:length(methods2)) cat(" ",methods2[i],"\n",sep="")     
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
        cat("Parameters:", "\n")
        cat(" Significance level: ", res$alpha, "\n", sep = "")
        if (sum(methods == "T.I.D.") == 1) 
            cat(" TID threshold:", res$thrTID, "\n")
        if (sum(methods == "Stand.") == 1) 
            cat(" Standardization threshold:", res$thrSTD, "\n")
        if (sum(methods == "M-H") == 1) {
            if (res$MHstat == "MHChisq") 
                MHmet <- "Chi-square statistic"
            else MHmet <- "Log odds-ratio statistic"
            cat(" Mantel-Haenszel DIF statistic:", MHmet, "\n")
            if (res$correct) corr <- "Yes"
            else corr <- "No"
            cat(" Mantel-Haenszel continuity correction:", corr, "\n")
            if (res$exact) cat(" Type of Mantel-Haenszel test: exact test","\n")
            else cat(" Type of Mantel-Haenszel test: asymptotic test","\n")
        }
        if (sum(methods == "Stand.") == 1) {
            stdw <- ifelse(res$stdWeight == "total", "both groups", 
                ifelse(res$stdWeight == "focal", "the focal group", 
                  "the reference group"))
            cat(" Weights for standardized P-DIF statistic: based on", 
                stdw, "\n")
        }
        if (sum(methods == "BD") == 1) {
            if (res$BDstat == "BD") 
                BDmet <- "Breslow-Day statistic"
            else BDmet <- "trend test statistic"
            cat(" Breslow-Day DIF statistic:", BDmet, "\n")
        }
        if (sum(methods == "Logistic") == 1) {
            cat(" Logistic regression DIF statistic:", res$criterion, 
                "statistic", "\n")
            resLog <- ifelse(res$type == "both", "both", ifelse(res$type == 
                "udif", "uniform", "non uniform"))
            resLog2 <- ifelse(res$type == "both", "effects", 
                "effect")
            cat(" DIF effect(s) tested by logistic regression:", 
                resLog, "DIF", resLog2, "\n")
        }
        if (sum(methods == "Lord" | methods == "Raju") >= 1) {
            cat(" Item response model:", res$model, "\n")
            if (res$model != "1PL" | res$engine == "ltm") 
                cat(" Engine 'ltm' for item parameter estimation", 
                  "\n")
            else cat(" Engine 'lme4' for item parameter estimation", 
                "\n")
            if (res$model == "1PL" & res$engine == "ltm") {
                if (is.null(res$discr)) 
                  cat(" Common discrimination parameter: estimated from 'ltm'", 
                    "\n")
                else cat(" Common discrimination parameter: fixed to ", 
                  res$discr, "\n", sep = "")
            }
            if (!is.null(res$c)) {
                if (length(res$c) == 1) 
                  cat(" Common pseudo-guessing value: ", res$c, 
                    "\n", sep = "")
                else {
                  pg <- cbind(res$c)
                  rownames(pg) <- res$names
                  colnames(pg) <- "c"
                  cat(" Common pseudo-guessing values:", "\n")
                  print(pg)
                  cat("\n")
                }
            }
        }
        if (sum(methods == "Raju") == 1) {
            if (res$signed) cat(" Type of Raju's Z statistic: signed area", "\n")
            else cat(" Type of Raju's Z statistic: unsigned area", "\n")
}
        if (res$purification & is.null(res$anchor.names)) {
            cat(" Item purification: Yes", "\n", "\n")
            cat(" Item purification results:", "\n", "\n")
            co <- rep("Yes", length(res$convergence))
            co[!res$convergence] <- "No"
            resConv <- data.frame(rbind(co, res$nrPur))
            colnames(resConv) <- colnames(res$DIF)
            rownames(resConv) <- c("Convergence", "Iterations")
            print(format(resConv, justify = "centre"))
            cat("\n")
        }
        else cat(" Item purification: No", "\n", "\n")
        cat("Comparison of DIF detection results:", "\n", "\n")
        nr <- NULL
        for (i in 1:nrow(res$DIF)) nr[i] <- paste(length(res$DIF[i, 
            ][res$DIF[i, ] == "DIF"]), "/", ncol(res$DIF), sep = "")
        MAT <- cbind(res$DIF, nr)
        colnames(MAT)[ncol(MAT)] <- "#DIF"
        print(format(MAT, justify = "centre"), quote = FALSE)
    }
    if (!x$save.output) 
        cat("\n", "Output was not captured!", "\n")
    else {
        if (x$output[2] == "default") 
            wd <- paste(getwd(), "/", sep = "")
        else wd <- x$output[2]
        fileName <- paste(wd, x$output[1], ".txt", sep = "")
        cat("\n", "Output was captured and saved into file", 
            "\n", " '", fileName, "'", "\n", "\n", sep = "")
    }
}