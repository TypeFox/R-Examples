html.report <-
function(result,dir=getwd(),file="report", CSS="R2HTML", ...){
    target <- HTMLInitFile(outdir=dir,filename=file, Title="Agreement Statistics", ...);
    HTML("<center>Agreement Statistics</center>",file=target);

    result$Conf_Limit$CCC <- round(result$Conf_Limit$CCC,4);
    result$Estimate$CCC <- round(result$Estimate$CCC,4);
    result$Conf_Limit$Precision <- round(result$Conf_Limit$Precision,4);
    result$Estimate$Precision <- round(result$Estimate$Precision,4);
    result$Conf_Limit$Accuracy <- round(result$Conf_Limit$Accuracy,4);
    result$Estimate$Accuracy <- round(result$Estimate$Accuracy,4);
    result$Conf_Limit$CP <- round(result$Conf_Limit$CP,4);
    result$Estimate$CP <- round(result$Estimate$CP,4);
    result$Conf_Limit$RBS <- round(result$Conf_Limit$RBS,2);
    result$Estimate$RBS <- round(result$Estimate$RBS,2);
    result$Conf_Limit$TDI <- round(result$Conf_Limit$TDI,result$Data$dec);
    result$Estimate$TDI <- round(result$Estimate$TDI,result$Data$dec);

    out <- cbind(array(c("Estimate","95% Conf.Limit","Allowance"),c(3,1)),rbind(format(result$Estimate),format(result$Conf_Limit),format(result$Allowance)));
    out <- replace(out, out=="NA", " ");
    HTMLCSS(file=target, CSSfile=CSS);
    HTML(out,file=target);
    plot(result);
    HTMLplot(file=target,GraphDirectory=dir);
    HTMLEndFile();
}

