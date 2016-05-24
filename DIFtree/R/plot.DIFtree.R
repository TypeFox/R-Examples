#' Plotting of Item focussed Trees
#' 
#' @description
#' Visualization of trees for items with DIF identified by item focussed recursive partitioning 
#' based on the Rasch Model or the Logistic Regression Approach for DIF detection.
#' 
#' @param x Object of class \code{\link[DIFtree]{DIFtree}}
#' @param item Number of the item, for which the tree shall be plotted 
#' @param component Component of the model for which the tree shall be plotted; 
#' can be \code{"intercept"} or \code{"slope"}. For \code{"Rasch"} model only one tree of item difficulties 
#' is available for each DIF item and therefore \code{component} will be ignored. 
#' @param cex.lines Width of branches of the tree
#' @param cex.branches Size of the labels of branches of the tree 
#' @param cex.coefs Size of coefficients given in the terminal nodes of the tree
#' @param cex.main Size of the title of the tree
#' @param title Optional title, which is added to the tree;
#' if \code{title=NULL} the title is the number of the plotted item.
#' @param ... Further arguments passed to or from other methods
#' 
#' @author Moritz Berger <moritz.berger@stat.uni-muenchen.de> \cr \url{http://www.statistik.lmu.de/~mberger/}
#' 
#' @references 
#' Berger, Moritz and Tutz, Gerhard (2015): Detection of Uniform and Non-Uniform Differential Item Functioning 
#' by Item Focussed Trees, Cornell University Library, arXiv:1511.07178
#' 
#' Tutz, Gerhard and Berger, Moritz (2015): Item Focused Trees for the Identification of Items
#' in Differential Item Functioning, Psychometrika, published online, DOI: 10.1007/s11336-015-9488-3 
#' 
#' @seealso \code{\link[DIFtree]{DIFtree}}, \code{\link[DIFtree]{predict.DIFtree}}, \code{\link[DIFtree]{summary.DIFtree}}
#' 
#' @examples 
#' data(data_sim)
#'  
#' Y <- data_sim[,1]
#' X <- data_sim[,-1]
#'  
#' \dontrun{
#'  
#' mod <- DIFtree(Y=Y,X=X,model="Logistic",type="udif",alpha=0.05,nperm=1000,trace=TRUE)
#'  
#' plot(mod,item=1)
#' }
#
#' @method plot DIFtree
#' @export
#' @importFrom plotrix draw.ellipse
#' @importFrom grDevices grey
#' @importFrom graphics lines plot.new plot.window points rect text 

plot.DIFtree <-
function(x, # object of class DIFtree
                         item,
                         component="intercept",
                         cex.lines=2,cex.branches=1,cex.coefs=1,cex.main=1,title=NULL,...){
  
  if(is.null(x$splits)){
    cat("There is no plot available in the case of no DIF item")
  } else{
    
    X <- x$X
    
    model <- which(c("Rasch","Logistic") %in% paste(x$call))
    if(model==1){
      info     <- x$splits[which(x$splits[,"item"]==item),]
      if(nrow(info)==0){
        beta_item <- x$betas_nodif[paste0("beta",item)]
        cat("Item", item, "is no DIF item. There is no tree to plot.\n")
        cat("Estimated item difficulty:", beta_item)
      } else{
        betas_item <- x$coefficients$betas_dif[[which(names(x$coefficients$betas_dif)==item)]]
        ptree(info,item,betas_item,X,cex.lines,cex.main,cex.branches,cex.coefs,title)
      }
    }
    if(model==2){
      type <- which(c("udif","dif","nudif") %in% paste(x$call))
      if(type==1){
        info     <- x$splits[which(x$splits[,"item"]==item),]
        if(nrow(info)==0){
          gamma_item <- x$coefficients$gammas_nodif[paste0("gamma",item)]
          cat("Item", item, "is no DIF item. There is no tree to plot.\n")
          cat("Estimated intercept:", gamma_item)
        } else{
          gammas_item <- x$coefficients$gammas_dif[[which(names(x$coefficients$gammas_dif)==item)]]
          ptree(info,item,gammas_item,X,cex.lines,cex.main,cex.branches,cex.coefs,title)
        }
      }
      if(type==2){
      
        dif_items <- unique(c(x$splits[[1]][,"item"],x$splits[[2]][,"item"]))
        if(is.null(dif_items)){
          cat("There is no plot available in the case of no DIF item")
        } else{
          if(!item %in% dif_items){
            gamma_item <- x$coefficients$gammas_nodif[paste0("gamma",item)]
            alpha_item <- x$coefficients$alphas_nodif[paste0("alpha",item)]
            cat("Item", item, "is no DIF item. There are no trees to plot.\n")
            cat("Estimated intercept:", gamma_item,"\n")
            cat("Estimated slope:", alpha_item)
          } else{
            
            if(!(component %in% c("intercept","slope"))){
              stop(paste("Component",component,"undefined. Must be 'intercept' or 'slope'."))
            }
            
            if(component=="intercept"){
              info     <- x$splits[[1]][which(x$splits[[1]][,"item"]==item),]
              if(is.null(info) || nrow(info)==0){
                gamma_item <- x$coefficients$gammas_dif[[which(names(x$coefficients$gammas_dif)==item)]]
                cat("Item", item, "has no DIF in intercept. There is no tree to plot.\n")
                cat("Estimated intercept:", gamma_item)
              } else{
                gammas_item <- x$coefficients$gammas_dif[[which(names(x$coefficients$gammas_dif)==item)]]
                ptree(info,item,gammas_item,X,cex.lines,cex.main,cex.branches,cex.coefs,title)
              }
            } else{
              info     <- x$splits[[2]][which(x$splits[[2]][,"item"]==item),]
              if(is.null(info) || nrow(info)==0 ){
                alpha_item <- x$coefficients$alphas_dif[[which(names(x$coefficients$alphas_dif)==item)]]
                cat("Item", item, "has no DIF in slope. There is no tree to plot.\n")
                cat("Estimated slope:", alpha_item)
              } else{
                alphas_item <- x$coefficients$alphas_dif[[which(names(x$coefficients$alphas_dif)==item)]]
                ptree(info,item,alphas_item,X,cex.lines,cex.main,cex.branches,cex.coefs,title)
              }
            }
          }
        }
      }
      if(type==3){
        
        info     <- x$splits[which(x$splits[,"item"]==item),]
        if(nrow(info)==0){
          gamma_item <- x$coefficients$gammas_nodif[paste0("gamma",item)]
          alpha_item <- x$coefficients$alphas_nodif[paste0("alpha",item)]
          cat("Item", item, "is no DIF item. There are no trees to plot.\n")
          cat("Estimated intercept:", gamma_item,"\n")
          cat("Estimated slope:", alpha_item)
        } else{
          
          if(!(component %in% c("intercept","slope"))){
            stop(paste("Component",component,"undefined. Must be 'intercept' or 'slope'."))
          }
          
          if(component=="intercept"){
            gammas_item <- x$coefficients$gammas_dif[[which(names(x$coefficients$gammas_dif)==item)]]
            ptree(info,item,gammas_item,X,cex.lines,cex.main,cex.branches,cex.coefs,title)
          } else{
            alphas_item <- x$coefficients$alphas_dif[[which(names(x$coefficients$alphas_dif)==item)]]
            ptree(info,item,alphas_item,X,cex.lines,cex.main,cex.branches,cex.coefs,title)
          }
        }
      }
    }
  }
  invisible(x)
}
