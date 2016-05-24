#' decision tree in a recursive way
#'@export
#'@param data the dataset including the response
#'@param Y the name of the response
#'@param modele (optional) vector of names of covariates allowed in the tree
#'@param kill vector of the names to kill (variables won't be used in the tree)
#'@param index to give a number to the plot
#'@param print boolean to print the tree parameters
#'@param plot boolean to plot the tree
#'@param main the main title if plot=TRUE
#'@param sub the subtitle (if NULL it is automatically added)
#'@param lang the language for the automatic subtitle in the plot
#'@return returns the tree as an "rpart" object and the modele as a vector of the names of the covariates the tree could have used (to give as an input of the function).
#'\item{modele}{ vector of the names of the covariates the tree could have used} 
#'\item{tree}{ the regression tree as an "rpart" object} 
#' @examples
#'    \dontrun{

#'    data<-mtcars
#'    require(CorReg)
#' main="Regression tree of cars consumption (in mpg)"
#'   mytree=recursive_tree(data = data,Y ="mpg" ,main=main)
#'    #want to try without cylinder and disp
#'   mytree2=recursive_tree(data = data,Y ="mpg" ,kill=c("cyl","disp"),modele=mytree$modele,main=main)
#'   
#'     }
recursive_tree<-function(data=data,Y="Y",modele=NULL,kill=NULL,index=NULL,print=TRUE,plot=TRUE,main=NULL,sub=NULL,lang=c("en","fr")){
   if (is.null(modele)){
      modele=colnames(data)
      modele=modele[modele!=Y]   
   }
   if(!is.null(kill)){
      modele=modele[!modele%in%kill]
   }
   
   formule=as.formula(paste(Y," ~",paste(modele,collapse="+")))
   arbre=rpart(formule,data)
   if(plot){
      opar<-par()
      par( xpd=NA)
      if(is.null(sub)){
         lang=lang[1]
         if(lang=="fr"){
            sub=paste("oui a gauche, non a droite. Moyenne (global", round(arbre$frame$yval[1],digits=3),") et effectif (global ", arbre$frame$n[1],") ")
            vertical= "la hauteur indique la significativite"        
         }else if (lang=="en"){
            sub=paste("True left, False right. Mean response value (global", round(arbre$frame$yval[1],digits=3),") and effectives (global ", arbre$frame$n[1],") ")
            vertical=" heights indicates significativity"         
         }
      }
      plot(arbre)
      text(arbre, use.n=TRUE)
      title(main = main,ylab=vertical,xlab = sub, col.main = "red", col.lab = gray(.5),cex.main = 1.2, cex.lab = 1.0, font.main = 4, font.lab = 3)
      par( xpd=opar$xpd)
   }
   if(print){
      print(arbre)
   }
   return(list(modele=modele,tree=arbre))
}