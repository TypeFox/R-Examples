# PerfMeas.R 
# November 2011
# Modified july 2012
# February 2014: Corrected the degenerate case in AUC.single and corrected a bug in precision.at.recall.level

library("graph");
library("RBGL");
library("limma");

# dyn.load("PerfMeas.so");

#################################################################
######## 1. AUC measure #########################################
#################################################################

# Function to compute the AUC for a single class
# Input
# pred : numeric vector of  the values  of the predicted labels
# labels : vector of the true labels (0 negative, 1 positive examples)
# Output:
# a numeric value corresponding to the AUC
AUC.single <- function(pred,labels) {
  
  if (length(pred)!=length(labels))
         stop("AUC.single: lengths of true and predicted labels do not match.");
  if (any((labels!=0) & (labels!=1)))
       stop("AUC.single: labels variable must take values 0 or 1");
  if ((all(labels==0)) || (all(labels==1))) # considering the degenerate case where we have labels all equal ...
    return (0);
  if (sum(abs(diff(pred))) == 0)  { # if predictions are all equal then gives the ratio of the smaller class
    n <- length(pred);
	return(min(sum(labels==0)/n, sum(labels==1)/n));
  }
  return (auROC(labels, pred));
}



# Function that computes AUC for each class
# Both the target e predicted matrices have a number of rows equal to the number of examples
# and a number of columns equal to the number of the classes.
# Input:
# target : matrix with the target multilabels
# predicted : a numeric matrix with predicted values
# g : a tree of the multilabels. If g is missing no per.level results are computed
# root : the name of the root node (def. "00")
# Output :
# a list with three elements:
# average :  the average AUC across classes.              
# per.level : a named vector with average  AUC for each level of the hierarchy.
#                       names correspond to levels
# per.class : a named vector with AUC for each class.
#                       Names correspond to classes
AUC.single.over.classes <- function(target, predicted, g, root="00") { 
       n.examples <- nrow(target);
       n.classes <- ncol(target);
       if ((n.examples!=nrow(predicted)) || (n.classes!=ncol(predicted)))
          stop ("AUC.single.over.classes: number of rows or columns do not match between target and predicted classes");
       
       class.names <- colnames(target);
       n.classes <- length(class.names);
       
       # per class AUC
       per.class <- numeric(n.classes);
       names(per.class) <- class.names;
       
       for (i in 1:n.classes) {
         pred.labels <- predicted[,i];
	     true.labels <- target[,i];
         per.class[i] <- AUC.single(pred.labels,true.labels);       
       }
       
       # average over classes precision recall specificity F accuracy
       average <- mean(per.class);
       
	   if (!missing(g)) {
    	 # average per level AUC
    	 levels <- get.all.nodes.by.depth(g,  root=root);
    	 n.levels <- length(levels);
    	 level.names <- 1:n.levels;
    	 per.level <- numeric(n.levels);
    	 names(per.level) <-  level.names;    
    	 for (i in 1:n.levels)  {
           per.level[i] <- mean(per.class[levels[[i]]]);
    	 }
	   } else
	     return(list(average=average, per.level=NULL, per.class=per.class));
       	   
       return(list(average=average, per.level=per.level, per.class=per.class));
}



# Function to compute means across folds of AUC.single.over.classes
# Input:
# y : a list of lists. The compoents of the outer list is a list returned from the method AUC.single.over.classes
# Output:
# a list obtained by averaging the results across folds of the input x. The components are:
# a list with three elements:
# average :  the average AUC across classes.              
# per.level : a named vector with average  AUC for each level of the hierarchy.
#                       names correspond to levels
# per.class : a named vector with AUC for each class.
#                       Names correspond to classes
compute.mean.AUC.single.over.classes <- function (y) {
   n <- length(y);
   average <- y[[1]]$average;
   per.level <- y[[1]]$per.level;
   per.class <- y[[1]]$per.class;
   for (i in 2:n) {
      average <- average + y[[i]]$average;
      per.level <- per.level + y[[i]]$per.level;
      per.class <- per.class + y[[i]]$per.class;     
   }
   average <- average / n;
   per.level <- per.level / n;
   per.class <- per.class / n;
   return(list(average=average, per.level=per.level, per.class=per.class));
} 



############################################################
######## 2. F-score ########################################
############################################################

# Function to compute the F-measure for a single class
# Input
# pred : vector of the predicted labels
# labels : vector of the true labels
# Note that 0  stands for negative and 1 for positive.
# In general the first level is negative and the second positive
# Output:
# a named vector with six elements:
# P: precision
# R : recall (sensitivity)
# S : specificity
# F : F measure
# A : 0/1 loss accuracy
# npos : number of positives

F.measure.single <- function(pred,labels) {
      
   if (length(pred)!=length(labels))
       stop("F.measure: lengths of true and predicted labels do not match.");
   neg.ex <- which(labels <= 0);	
   np <- length(neg.ex);
   pos.ex <- which(labels > 0);
   npos <- length(pos.ex);	 
   TP <- sum(pred[pos.ex] > 0);
   FN <- sum(pred[pos.ex] <= 0);	
   TN <- sum(pred[neg.ex] <= 0);
   FP <- sum(pred[neg.ex] > 0);	           
   acc <- (TP+TN)/length(labels);
   if ((TP+FP) == 0)
     precision <- 0
   else 
     precision <- TP/(TP+FP);
   if ((TP+FN) == 0)
     recall <- 0
   else
     recall <- TP/(TP+FN);
   if ((TN+FP) == 0)
     specificity <- 0
   else
     specificity <- TN/(TN+FP);
   if ((precision+recall) == 0)
      F <- 0
   else
      F = 2 *(precision*recall) / (precision+recall); 

   res <- c(precision,recall,specificity,F,acc, npos);
   names(res) <- c("P", "R", "S", "F", "A", "Pos.");
   return (res);
}


# Function that computes  precision, recall, specificity and F-measure for each class
# Both the target e predicted matrices have a number of rows equal to the number of examples
# and a number of columns equal to the number of the classes.
# Input:
# target : matrix with the target multilabels
# predicted : matrix with the predicted multilabels
# g : a tree of the multilabels. If missing no per level results are computed
# root : the name of the root node (def. "00")
# Output :
# a list with three elements:
# average :  a named vector with the average precision, recall, specificity, F-measure and accuracy across classes.              
# per.level : a named matrix with average  precision, recall, specificity, F-measure and accuracy for each level of the hierarchy.
#                       Named rows correspond to levels,
#                       named columns correspond respectively to precision, recall, specificity, F-measure and accuracy. 
#                       If g is missing non per.level results are computed.
# per.class : a named matrix with  precision, recall, specificity, F-measure and accuracy for each class.
#                       Named rows correspond to classes,
#                       named columns correspond respectively to precision, recall, specificity, F-measure and accuracy
F.measure.single.over.classes <- function(target, predicted, g, root="00") { 
  n.examples <- nrow(target);
  n.classes <- ncol(target);
  if ((n.examples!=nrow(predicted)) || (n.classes!=ncol(predicted)))
     stop ("F.measure.single.over.classes: number of rows or columns do not match between target and predicted classes");

  class.names <- colnames(target);
  n.classes <- length(class.names);
  measures <- c("P", "R", "S", "F", "A", "Pos.");
  n.measures <- length(measures);
  # per class precision recall F accuracy
  per.class <- matrix(numeric(n.classes * n.measures), ncol=n.measures, dimnames=list(class.names, measures));      
  for (i in 1:n.classes) {
    pred.labels <- predicted[,i];
	true.labels <- target[,i];
    per.class[i,] <- F.measure.single(pred.labels,true.labels);       
  }

  # average over classes precision recall specificity F accuracy
  average <- apply(per.class, 2, mean);

  # average per level precision recall F accuracy
  if (!missing(g)) {
	levels <- get.all.nodes.by.depth(g,  root=root);
	n.levels <- length(levels);
	level.names <- 1:n.levels;
	per.level <- matrix(numeric(n.levels * n.measures), ncol=n.measures, dimnames=list(level.names, measures));  	
	for (i in 1:n.levels)  {
      if (length(levels[[i]]) == 1)
		ff <- matrix(per.class[levels[[i]],], nrow=1)
	  else
		ff <- per.class[levels[[i]],];
      per.level[i, ] <- apply(ff, 2, mean);
	}
  } else
    per.level <- NULL;

  return(list(average=average, per.level=per.level, per.class=per.class));
}


# Function to compute means across folds of F.measure.single.over.classes
# Input:
# y : a list of lists. The compoents of the outer list is a list returned from the method F.measure.single.over.classes
# Output:
# a list obtained by averaging the results across folds of the input x. The components are:
# average :  a named vector with the average precision, recall, specificity, F-measure and accuracy across classes across folds            
# per.level : a named matrix with average  precision, recall, specificity, F-measure and accuracy for each level of the hierarchy across folds
#                       Named rows correspond to levels,
#                       named columns correspond respectively to precision, recall, specificity, F-measure and accuracy
# per.class : a named matrix with  precision, recall, specificity, F-measure and accuracy for each class across folds
#                       Named rows correspond to classes,
#                       named columns correspond respectively to precision, recall, specificity, F-measure and accuracy
compute.mean.F.measure.single.over.classes <- function (y) {
   n <- length(y);
   average <- y[[1]]$average;
   per.level <- y[[1]]$per.level;
   per.class <- y[[1]]$per.class;
   for (i in 2:n) {
      average <- average + y[[i]]$average;
      per.level <- per.level + y[[i]]$per.level;
      per.class <- per.class + y[[i]]$per.class;     
   }
   average <- average / n;
   per.level <- per.level / n;
   per.class <- per.class / n;
   return(list(average=average, per.level=per.level, per.class=per.class));
} 



############################################################
######## 3. Precision at a given recall ####################
############################################################

# Function to compute the precision at a given recall level for a single class
# Input:
# scores : vector of the predicted scores in [0,1]
# labels : 0/1 vector of the true labels 
# rec.level: the desired recall level
# Output: 
# res : the precision at the requested recall
# Function to compute the precision at a given recall level for a single class
# Input:
# scores : vector of the predicted scores in [0,1]
# labels : 0/1 vector of the true labels 
# rec.level: the desired recall level
# Output: 
# res : the precision at the requested recall
precision.at.recall.level <- function(scores, labels, rec.level=0.2){
  n<-length(scores); 
  if (n!=length(labels))
	stop("precision.at.recall.level: length of labels and scores does not match");
  if(length(which(labels > 0)) == 0)
	  return(res=0);
  # sort.scores <- sort(scores, decreasing=TRUE);
  scores.ordered <- order(scores, decreasing=TRUE);	
  precision <- recall <- rep(0, n);
  
  res <- .C("prec_recall", as.double(precision), as.double(recall), as.integer(scores.ordered), 
             as.integer(labels), as.integer(n), PACKAGE="PerfMeas");
  
  precision <- res[[1]];
  recall <- res[[2]];
  #cat("Precision: \n"); cat(precision, "\n");
  #cat("Recall: \n"); cat(recall, "\n");
  #cat(res[[4]]);
  
  # choosing the first element in recall > rec.level
  tmp_ind <- which(recall - rec.level >= 0);	  
  # there could be some identical recall values, so chosen.index could be a vector.
  chosen.index <- which(recall == min(recall[tmp_ind]));
  return(res=max(precision[chosen.index]));	
}




# Function that computes precision at a given recall level for each class
# Both the target e predicted matrices have a number of rows equal to the number of examples
# and a number of columns equal to the number of the classes.
# Input:
# target : matrix with the target multilabels
# predicted : a numeric matrix with predicted values
# g : a tree of the multilabels.If g is missing, no hierarchy is considered are per level results are NULL
# rec.level : requested recall level (def. 0.2)
# root : the name of the root node (def. "00")
# Output :
# a list with three elements:
# average :  the average precision at a given recall level across classes.              
# per.level : a named vector with average  precision at a given recall level for each level of the hierarchy.
#                       names correspond to levels
# per.class : a named vector with precision at a given recall level for each class.
#                       Names correspond to classes
precision.at.recall.level.over.classes <- function(target, predicted, g, rec.level=0.2, root="00") { 
  n.examples <- nrow(target);
  n.classes <- ncol(target);
  if ((n.examples!=nrow(predicted)) || (n.classes!=ncol(predicted)))
     stop ("precision.at.recall.level.over.classes: number of rows or columns do not match between target and predicted classes");

  class.names <- colnames(target);
  n.classes <- length(class.names);

  # per class precision.at.recall.level
  per.class <- numeric(n.classes);
  names(per.class) <- class.names;

  for (i in 1:n.classes) {
    pred.labels <- predicted[,i];
	true.labels <- target[,i];
	per.class[i] <- precision.at.recall.level (pred.labels, true.labels, rec.level=rec.level);
  }

  # average over classes precision.at.recall.level
  average <- mean(per.class);

  # average per level precision.at.recall.level
  if (!missing(g)) {
	levels <- get.all.nodes.by.depth(g,  root=root);
	n.levels <- length(levels);
	level.names <- 1:n.levels;
	per.level <- numeric(n.levels);
	names(per.level) <-  level.names;    
	for (i in 1:n.levels)  {
      per.level[i] <- mean(per.class[levels[[i]]]);
	}
  } else
    per.level<-NULL;

  return(list(average=average, per.level=per.level, per.class=per.class));
}


# Function to compute the precision at multiple levels of recall for a single class
# Input:
# scores : vector of the predicted scores in [0,1]
# labels : 0/1 vector of the true labels 
# rec.levels: a vector with the desired recall level (def. 0.1 to 1 by 0.1 step)
# Output: 
# a list with 2 elements:
# - a vector with the precision at different recall levels
# - a vector with the f.score at different recall levels
precision.at.multiple.recall.level <- function(scores, labels, rec.levels=seq(from=0.1, to=1, by=0.1)){
  n<-length(scores); 
  n.levels <- length(rec.levels);
  if (n!=length(labels))
	stop("precision.at.recall.level: length of labels and scores does not match");
  if(length(which(labels > 0)) == 0)
	  return(list(precisions=rep(0,n.levels),f.score=rep(0,n.levels)));
  scores.ordered <- order(scores, decreasing=TRUE);	
  precision <- recall <- rep(0, n);
  res <- .C("prec_recall", as.double(precision), as.double(recall), as.integer(scores.ordered), 
      as.integer(labels), as.integer(n), PACKAGE="PerfMeas");
  
  precision <- res[[1]];
  recall <- res[[2]];
  
  precisions.at.rec.level <- numeric(n.levels);
  names(precisions.at.rec.level) <- rec.levels;
  for (i in 1:n.levels) {
	# choosing the first element in recall > rec.levels
	tmp_ind <- which(recall - rec.levels[i] >= 0);	  
	# there could be some identical recall values, so chosen.index could be a vector.
	chosen.index <- which(recall == min(recall[tmp_ind]));
	precisions.at.rec.level[i] <- max(precision[chosen.index]);
  }
  f.score <- (2 * precisions.at.rec.level * rec.levels)/(precisions.at.rec.level + rec.levels);
  f.score[is.nan(f.score)] <- 0;
  names(f.score) <- rec.levels;
  return(list(precisions=precisions.at.rec.level, f.score=f.score));	
}



# Function to compute the precision at multiple levels of recall for multiple classes
# Input:
# target : 0/1 matrix  of the true labels:  rows are examples, columns classes.
# predicted :  matrix of the predicted scores  with values in [0,1]. Rows are examples, columns classes 
#     S and L must have the same dimension ant the same examples and classes
# rec.levels: a vector with the desired recall level (def. 0.1 to 1 by 0.1 step)
# Output: 
# a list with two elements
# - PXR: a matrix with the precisions at different recall levels: rows are classes, columns 
#       precisions at different recall levels
# - avgPXR: a vector with the the average precisions at different recall levels across classes
precision.at.multiple.recall.level.over.classes <- function(target, predicted, rec.levels=seq(from=0.1, to=1, by=0.1)){
  n.examples <- nrow(predicted);
  n.classes <- ncol(predicted);
  n.rec.level <- length(rec.levels);
  if ((n.examples != nrow(target)) || (n.classes != ncol(target)))
    stop("precision.at.multiple.recall.level.over.classes: dimensions of score and label matrices do not agree");
  classes.names <- colnames(predicted);
  PXR <- matrix(numeric(n.classes * n.rec.level), nrow=n.classes, dimnames=list(classes.names, rec.levels));                
  for (i in 1:n.classes)
    PXR[i,] <- precision.at.multiple.recall.level(predicted[,i], target[,i], rec.levels)$precisions;
  avgPXR <- apply(PXR, 2, mean);
  names(avgPXR) <- rec.levels;
  return(list(avgPXR = avgPXR, PXR = PXR));
}
    
# Function to compute the precision at all recall levels  for a single class
# Input:
# scores : vector of the predicted scores in [0,1]
# labels : 0/1 vector of the true labels 
# Output: 
# a list with 3 elements:
# precision : precision at different thresholds
# recall : recall at different thresholds
# f.score : f.score at different thresholds
precision.at.all.recall.levels <- function(scores, labels){
  n<-length(scores); 
  if (n!=length(labels))
	stop("precision.at.recall.level: length of labels and scores does not match");
  if(length(which(labels > 0)) == 0)
	  return(list(res=0,precision=rep(0,n),recall=rep(0,n)));
  scores.ordered <- order(scores, decreasing=TRUE);	
  precision <- recall <- rep(0, n);
  res <- .C("prec_recall", as.double(precision), as.double(recall), as.integer(scores.ordered), 
            as.integer(labels), as.integer(n), PACKAGE="PerfMeas");
  
  precision <- res[[1]];
  recall <- res[[2]];
  
  f.score <- (2 * precision * recall)/(precision + recall);
  f.score[is.nan(f.score)] <- 0;
  return(list(precision=precision,recall=recall,f.score=f.score));	
}

# Function to compute multiple AUPRC (Area Under Precision and Recall Curves)
# Input:
# z : a list of lists. Each component list is a list returned from precision.at.all.recall.levels
#     that reports precision, recall and f-score results at different levels for different methods or tasks
# comp.precision: boolean. It TRUE (default) the AUPRC is computed otherwise the area under the F-score curve is computed
# Output:
# a named vector with the AUPRC (or the AUFRC) for different methods or tasks
AUPRC <- function(z, comp.precision=TRUE) {
  n <- length(z);
  curve.names <- names(z);
  if (is.null(names(curve.names)))
    curve.names<-as.character(1:n);
  integral <- numeric(n);
  names(integral) <- curve.names;
  for (i in 1:n)  {
    if (comp.precision)
      integral[i] <- trap.rule.integral(z[[i]][[2]], z[[i]][[1]])
    else
      integral[i] <- trap.rule.integral(z[[i]][[2]], z[[i]][[3]])
  }
  return(integral);
}


############################################################
######## 4. Utility functions ##############################
############################################################

# Function to plot multiple precision recall or f-score recall curves
# Input:
# y : a list of lists. Each component list is a list returned from precision.at.all.recall.levels
#     that reports precision, recall and f-score results at different levels for different methods or tasks
# curve.names : names of the compared methods to be reported in the legenda (def: numbers)
# f : file name. If is given, a postscript file is created, otherwise the output is rendered on
#                a window.
# cex.val : magnification value for characters (def. 0.6)
# height : heigth of the graph (def. 9)
# width : width of the graph (def 11)
# col : colors of the lines. 14 colors are given as default, but any vector of color from colors() 
#       can be used. Colors are recycled if length(col) < length(y)
# line.type : type of the line. Any valid vector of integer can be assigned (values between 1 
#             and 6, see lty in par() for details). Values are recycled if length(line.type) < length(y).
#             Def.: 1 (solid lines).
precision.recall.curves.plot <- function(y, range=seq(from=0, to=1, by=0.1), 
       curve.names=1:length(y), cex.val=0.6, f="", height=9, width=11,
	   col=c("black","red1","blue1","green1","darkgrey","brown1","yellow1","orange1",
	   "red4","blue4","green4","lightgrey","brown4","yellow4","orange4"),
	   line.type=1, leg=TRUE, pos=c(range[length(range)-2], range[length(range)]), do.grid=TRUE, plot.precision=TRUE,
	   trap.rule=TRUE)  {
  prec <- rec <- range;  
  n <- length(y);
  x.label<-"Recall";
  if (plot.precision) 
    y.label<-"Precision"
  else
    y.label<-"F-score";
  len.line.type <- length(line.type);
  mm.line <- as.integer(n/len.line.type); 
  m.line <-   n-(mm.line * len.line.type);            
  str.line=integer(0);
  if (mm.line>0)
    for (i in 1:mm.line)
	  str.line <- c(str.line, line.type);
  if (m.line>0)
	str.line <- c(str.line, line.type[1:m.line]);
	
  len.col <- length(col);
  mm.col <- as.integer(n/len.col); 
  m.col <- n-(mm.col * len.col);
  str.col=character(0);
  if (mm.col>0)
    for (i in 1:mm.col)
	  str.col <- c(str.col, col);
  if (m.col>0)
	str.col <- c(str.col, col[1:m.col]);
  
  if (f!="")
    postscript(f, paper="special", height=height, width=width, horizontal=F);
  plot(prec, rec, type="n",  xlab=x.label, ylab=y.label, xaxt="n");
  axis(side=1, at = range, labels = range);
  
  for (i in 1:n) 
    if (plot.precision)
      lines(y[[i]][[2]], y[[i]][[1]], type="l", lty=str.line[i], col=str.col[i])
	else
	  lines(y[[i]][[2]], y[[i]][[3]], type="l", lty=str.line[i], col=str.col[i]) 

  if (leg)  
     legend(x=pos[1], y=pos[2], curve.names, lty=str.line, col=str.col);
  
  if (do.grid)
    grid(lwd=1,col="gray");  
      
  if(f!="")
    dev.off(); 
	
  if (trap.rule) {
    integral <- numeric(n);
	names(integral) <- curve.names;
    for (i in 1:n)  {
	  if (plot.precision)
	    integral[i] <- trap.rule.integral(y[[i]][[2]], y[[i]][[1]])
	  else
	    integral[i] <- trap.rule.integral(y[[i]][[2]], y[[i]][[3]])
    }
	return(integral);
  }
}


# Function to plot multiple performance curves
# It may be used to compute e.g. precision or f-score at given recall levels.
# Input:
# m : a numeric matrix. Rows correspond to different methods and columns to precision or f-score at given recall values
# curve.names : names of the compared methods to be reported in the legenda (def: numbers). 
# f : file name. If is given, a postscript file is created, otherwise the output is rendered on
#                a window.
# cex.val : magnification value for characters (def. 0.6)
# height : heigth of the graph (def. 9)
# width : width of the graph (def 11)
# col : colors of the lines. 14 colors are given as default, but any vector of color from colors() 
#       can be used. Colors are recycled if length(col) < length(y)
# line.type : type of the line. Any valid vector of integer can be assigned (values between 1 
#             and 6, see lty in par() for details). Values are recycled if length(line.type) < length(y).
#             Def.: 1 (solid lines).
# patch.type: symbols to be plotted according to the pch parameter (see par in package graphics).
performance.curves.plot <- function(m, x.range=seq(from=0.1, to=1, by=0.1), y.range=c(0,1),
       curve.names=1:nrow(m), cex.val=0.6, f="", height=9, width=11,
	   col=c("black","red1","blue1","green1","darkgrey","brown1","yellow1","orange1",
	   "red4","blue4","green4","lightgrey","brown4","yellow4","orange4"),
	   line.type=1, patch.type=1:16, leg=TRUE, pos=c(x.range[length(x.range)-2], y.range[2]), do.grid=TRUE,
	   x.label="Recall", y.label="Precision")  {
  n <- nrow(m);
  len.line.type <- length(line.type);
  mm.line <- as.integer(n/len.line.type); 
  m.line <-   n-(mm.line * len.line.type);            
  str.line=integer(0);
  if (mm.line>0)
    for (i in 1:mm.line)
	  str.line <- c(str.line, line.type);
  if (m.line>0)
	str.line <- c(str.line, line.type[1:m.line]);
	
  len.col <- length(col);
  mm.col <- as.integer(n/len.col); 
  m.col <- n-(mm.col * len.col);
  str.col=character(0);
  if (mm.col>0)
    for (i in 1:mm.col)
	  str.col <- c(str.col, col);
  if (m.col>0)
	str.col <- c(str.col, col[1:m.col]);
	
  len.patch <- length(patch.type);
  mm.patch <- as.integer(n/len.patch);
  m.patch <- n - (mm.patch * len.patch);
  str.patch=integer(0);
  if (mm.patch > 0)
    for (i in 1:mm.patch)
	  str.patch <- c(str.patch, patch.type);
  if (m.patch>0)
	str.patch <- c(str.patch, patch.type[1:m.patch]);    
  
  if (f!="")
    postscript(f, paper="special", height=height, width=width, horizontal=F);
  plot(x.range, y=seq(from=y.range[1], to=y.range[2], along.with=x.range), type="n", xlab=x.label, ylab=y.label, xaxt="n", cex=cex.val, cex.axis=cex.val, cex.lab=cex.val);
  axis(side=1, at = x.range, labels = x.range, cex.axis=cex.val);
  
  for (i in 1:n) 
     lines(x.range, m[i,], type="b", lty=str.line[i], col=str.col[i], pch=str.patch[i], cex=cex.val);
  
  if (leg)  
     legend(x=pos[1], y=pos[2], curve.names, lty=str.line, col=str.col, pch=str.patch, cex=cex.val);
  
  if (do.grid)
    grid(lwd=1,col="gray");  
  if(f!="")
    dev.off(); 
}


# It groups set of nodes according to their depth in the graph
# Input:
# g: a graph 
# root : name of the root node (def. 00)
# Output:
# a  list of the nodes, grouped w.r.t. the distance from the root:
# the first element of the list corresponds to the nodes at distance 1, 
# the second to nodes at distance 2 and so on.
get.all.nodes.by.depth <- function(g,  root="00") {
  l<-dijkstra.sp(g, root);
  depth <- max(l$distance);
  levels <- list();
  for (i in 1:depth)
     levels[[i]] <- names(l$distance[l$distance==i])
  return(levels);
}


# Function that implements the trapezoidal rule for integration
# Input:
# x : abscissa values in increasing order
# y : ordinate values
# Output:
# value of the integral
trap.rule.integral <- function (x,y){
    if (length(x) != length(y))
	  stop("trap.rule.integral: length of x and y vectors must match");
	integral_value = 0.0;
    integral_value <- .C("trap_rule", as.double(x), as.double(y), as.integer(length(x)), 
             as.double(integral_value), PACKAGE="PerfMeas")[[4]];
    return (integral_value);
}



.onLoad <- function(libname=.libPaths(), pkgname="PerfMeas")
       library.dynam("PerfMeas", pkgname, libname);
	   
.onAttach <- function(libname=.libPaths(), pkgname="PerfMeas")
               packageStartupMessage("PerfMeas: Performance Measures for ranking and classification tasks.\n")
