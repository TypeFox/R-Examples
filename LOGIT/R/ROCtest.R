# Hilbe, J.M., Practical Guide to Logistic Regression
# Rafael de Souza, Eotvos Lorand Univ. 2015
#' @title  Display ROC curve and related AUC statistic, or sensitivity-specificity plot of glm
#' with binomial family.
#' @description Provides two options following the glm() function with binomial family.
#' 1: Senstivity-specificity plot with optimal cut point statistic
#' 2: ROC plot with Area Under Curve (AUC) statistic
#' @aliases ROCtest
#' @importFrom stats binomial glm predict
#' @import  pROC caret ggplot2 reshape
#' @usage ROCtest(model= model, fold=10, type=c("ROC","Sensitivity"))
#'
#' @format   \describe{
#' \item{x}{
#' The function has three arguments: modelname, folds, type of plot}
#' }
#' @param model  model name
#' @param fold  number of k-folds
#' @param type type of plot
#' @return plot
#' @note ROCtest() must be loaded into memory in order to be effectve. As a function in LOGIT,
#' it is immediately available to a user.
#' @details ROCtest is a post-estimation function for logistic regression, following the use
#' of glm(). Options to display a sensitivity-specificity plot or ROC curve are
#' available.
#' @seealso \code{\link{glm}}
#' @author Rafael de Souza, ELTE, Hungary,
#' Joseph M. Hilbe, Arizona State University.
#'
#' @references Hilbe, Joseph M. (2016), Practical Guide to Logistic Regression, Chapman & Hall/CRC.
#' Hilbe, Joseph M. (2009), Logistic Regression Models, Chapman & Hall/CRC.
#'@examples
#'  library(MASS)
#'  library(LOGIT)
#'  data(R84)
#'  R84$cage <- R84$age - mean(R84$age)
#'  R84$cdoc <- R84$docvis - mean(R84$docvis)
#'  mylogit <- glm(outwork ~ cdoc + female + kids + cage + factor(edlevel),
#'  family=binomial, data=R84)
#'  summary(mylogit)
#'  ROCtest(mylogit, fold=10, type="Sensitivity")
#'  ROCtest(mylogit, fold=10, type="ROC")
#'

#' @keywords models
#' @export
#'


########### Function


ROCtest <- function(model= model, fold=10, type=c("ROC","Sensitivity")) {
# To pass R CMD check
Sensitivity <- NULL
Specificity <- NULL
variable <- NULL
value <- NULL
probcut <- NULL

# Extract information from dataset
  data <- model$data
  response = as.character(model$formula[[2]])
  folds <- createFolds(data[,response], k = fold)
  AUC <- c()
  ROC_all <- c()
  cut <- c()

# Run k-fold and ROC analysis
   for (i in 1:fold){
    training <- data[-folds[[i]], ]
    testing <- data[folds[[i]], ]
    myroc <- glm(model$formula, family = binomial, data = training)


    ROC.a <- data.frame(True = training[,response],predicted = predict(myroc, newdata = training,type = "response"))

    F1 <- roc(ROC.a$True,ROC.a$predicted)

    ROC.b <- data.frame(True = testing[,response],predicted = predict(myroc, newdata = testing,type = "response"))

    ROC.b$class <- ROC.b$predicted
    ROC.b$class[which(ROC.b$class >= coords(F1,x = "best")[1])] <- 1
    ROC.b$class[which(ROC.b$class < coords(F1,x = "best")[1])] <- 0
    F2 <- roc(ROC.b$True,ROC.b$predicted)
    ROC_all <- rbind(ROC_all,ROC.b)



    AUC <- append(AUC,F2$auc)
#    cut<-append(cut,coords(F1,x="best")[1])
  }




# Plot information
GROC <- roc(ROC_all$True,ROC_all$predicted)

cut <- coords(GROC,x = "best",best.method=c("closest.topleft"))[[1]]
#xs<-coords(GROC,x="best")[[2]]
#ys<-coords(GROC,x="best")[[3]]

# Plot ROC curve and Confusion Matrix


if(type == "ROC") {
  g1 <- data.frame(Sensitivity=GROC$sensitivities,Specificity=1-GROC$specificities)


 gg<- ggplot(data=g1,aes(x=Specificity,y=Sensitivity))+geom_line(size=1.5,alpha=0.7,color="blue4")+
    theme(plot.title = element_text(hjust=0.5),axis.title.y=element_text(vjust=0.75),axis.title.x=element_text(vjust=-0.25),
          text = element_text(size=20))+
    annotate("text",size=7, x = 0.75, y = 0.14,color="black", label = paste("AUC: ",round(GROC$auc[1],3),sep="") ) +
    geom_segment(x=0,y=0,xend=1,yend=1,colour="red",size=1,linetype="dashed")+
    xlab("1-Specificity")

return(list(plot=gg,Observed=factor(ROC_all$True),Predicted=factor(ROC_all$class),cut=cut))

}

if(type=="Sensitivity") {

  g2<-data.frame(Sensitivity=GROC$sensitivities,Specificity=GROC$specificities,
                 probcut=GROC$thresholds)
  g2<-melt(g2, id=c("probcut"))

 gg<- ggplot(data=g2,aes(x=probcut,y=value,group=variable,colour=variable))+geom_line(size=1.5,alpha=0.7)+
   geom_point(size=2)+
    theme(plot.title = element_text(hjust=0.5),axis.title.y=element_text(vjust=0.75),axis.title.x=element_text(vjust=-0.25),
          text = element_text(size=20))+
    annotate("text",size=7, x = 0.75, y = 0.15,color="black", label = paste("Cut Point: ",round(mean(cut),3),sep="") ) +
    xlab("Probability cutoff")+ylab("Sensitivity/Specificity")+
   coord_cartesian(xlim=c(0,1),ylim=c(0,1))+scale_y_continuous(breaks=c(0.25,0.5,0.75,1))+
   scale_x_continuous(breaks=c(0,round(cut,2),0.25,0.5,0.75,1))
#   geom_segment(x=cut,y=0,xend=cut,yend=1,colour="red",size=1)

return(list(plot=gg,Observed=factor(ROC_all$True),Predicted=factor(ROC_all$class),cut=cut))

}
gg
}















