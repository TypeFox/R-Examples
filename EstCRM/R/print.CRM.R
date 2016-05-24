print.CRM <-
function(x, ...){
cat("***********************************************************************","\n")
cat("EstCRM -- An R Package for Estimating Samejima's Continuous IRT Model Parameters","\n")
cat("         ","Via Marginal Maximum Likelihood Estimation and EM Algorithm","\n")
cat("","\n")
cat("Version 1.4  2015","\n")
cat("","\n")
cat("Cengiz Zopluoglu","\n")
cat("","\n")
cat("University of Miami - Educational and Psychological Studies","\n")
cat("","\n")
cat("c.zopluoglu@miami.edu","\n")
cat("*************************************************************************","\n")
cat("Please use the following citation for any use of the package in any","\n")
cat("publication:","\n")
cat("","\n")
cat("     Zopluoglu, C.(2012). EstCRM: An R package for Samejima's continuous IRT model.","\n") 
cat("         Applied Psychological Measurement,36,",sprintf("%7s","149_150."),"\n") 
cat("","\n")
cat("*************************************************************************","\n")
cat("","\n")
cat("Processing Date: ",date(),"\n")
cat("","\n")
cat("Number of Items: ",nrow(x$param),"\n")
cat("Number of Subjects: ",nrow(x$data),"\n")
cat("","\n")
cat("Item Descriptive Statistics","\n")
cat("","\n")
cat("                  ","Raw Scores","             ",
"Transformed Scores(Z scale)","\n")
cat("            ",sprintf("%6s %6s %6s %6s","Mean","SD","Min","Max"),"        ",
sprintf("%6s %6s","Mean","SD"),"\n")
for(i in 1:nrow(x$param)){
cat("     ",sprintf("%6s %6.2f %6.2f %6.2f %6.2f",colnames(x$data)[i],
x$descriptive[i,1],x$descriptive[i,2],x$descriptive[i,3],x$descriptive[i,4]),
"       ",sprintf("%6.2f %6.2f",x$descriptive[i,5],x$descriptive[i,6]),"\n")
}
cat("","\n")
cat("Iteration was terminated at EM cycle",length(x$iterations),"\n")
cat("","\n")
cat("The difference in loglikelihoods between the last two EM cycles:",round(x$dif,7),"\n")
cat("","\n")
cat("The largest parameter change at the last EM cycle is:",
max(as.data.frame(x$iterations[length(x$iterations)])-as.data.frame(x$iterations[length(x$iterations)-1])),
"\n")
cat("","\n")
cat("Loglikelihood:",x$LL,"\n")
cat("","\n")
cat("Final Item Paramater Estimates:","\n")
cat("","\n")
cat("          ",sprintf("%6s %6s %8s","a","b","alpha"),"\n")
for(i in 1:nrow(x$param)){
cat("     ",sprintf("%6s %6.3f %6.3f %6.3f",colnames(x$data)[i],
x$param[i,1],x$param[i,2],x$param[i,3]),"\n")
}
cat("","\n")
cat("Standard Error of the Final Item Parameter Estimates:","\n")
cat("","\n")
cat("          ",sprintf("%6s %6s %8s","a","b","alpha"),"\n")
for(i in 1:nrow(x$param)){
cat("     ",sprintf("%6s %6.3f %6.3f %6.3f",colnames(x$data)[i],
x$std.err[i,1],x$std.err[i,2],x$std.err[i,3]),"\n")
}
cat("","\n")
}

