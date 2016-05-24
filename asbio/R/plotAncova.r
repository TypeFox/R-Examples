plotAncova <- function(model, pch = NULL, lty = NULL, col = NULL, leg.loc = "topright", leg.cex = 1, leg.bty = "o", leg.bg = par("bg"), legend.title = NULL, ...) {
  if(model$contrasts!="contr.treatment") stop("Please use treatment contrasts")
  if(ncol(model$model)>3) stop("function can only be used one with one-way ANCOVAs")
    
    x.name <- names(attr(model$terms,"dataClasses"))[attr(model$terms,"dataClasses")=="character"|attr(model$terms,"dataClasses")=="factor"]
    x1 <- model$model[,names(model$model)==x.name]
      cats <- unique(x1); r <- length(cats); n <- dim(model$model)[1]
    
    preds.vec <- names(attr(model$terms,"dataClasses"))[2:3]
    con.name <- preds.vec[preds.vec!=x.name] 
    con <- model$model[,names(model$model)==con.name]
    y <- model$model[,1]
  
  if(is.null(pch))  pch <- as.numeric(as.factor(x1))
  if(is.null(col))  col <- as.numeric(as.factor(x1))
  if(is.null(lty))  lty <- as.numeric(as.factor(x1))
  
  if(length(pch)==1) pch <- rep(pch, n)
  if(length(col)==1) col <- rep(col, n)
  if(length(lty)==1) lty <- rep(lty, n)
   
  mat <- data.frame(y = y, con = con, x = x1, pch, col, lty)
  o <- order(mat$x)
  mat <- mat[o,]
  names <- unique(mat$x)  
  
  leg.pch <- unique(mat$pch)
    if(length(leg.pch)==1) leg.pch <- rep(leg.pch, r) 
  leg.col <- unique(mat$col)
    if(length(leg.col)==1) leg.col <- rep(leg.col, r)
  leg.lty <- unique(mat$lty)
    if(length(leg.lty)==1) leg.lty <- rep(leg.lty, r)
  
      
      if(length(attr(model$terms,"term.labels"))==2){
          cof <- t(coef(model))
          slope <- cof[,colnames(cof)==con.name]
          subcof <- as.matrix(cof[,colnames(cof)!=con.name])
          int <- c(subcof[1],subcof[1]+subcof[2:r])
          
          
          plot(con, y, pch = pch, col = col,...)
          for(i in 1:r){
          abline(int[i],slope,lty=leg.lty[i],col=leg.col[i])
          }
      }

      if(length(attr(model$terms,"term.labels"))==3){
      
          new.model <- lm(mat$y ~ -1 + mat$x + mat$x:mat$con)
          
          cof <- coef(new.model)
          
          int <- cof[1:r]
          slope <- cof[(r + 1) : (2 * r)]
          
          plot(con, y, pch = pch, col = col,...)
          for(i in 1:r){
          abline(int[i],slope[i],lty=leg.lty[i],col=leg.col[i])
          }
      }

  if(leg.loc!="n"){
  legend(leg.loc, legend = names, pch = leg.pch, col = leg.col, lty = leg.lty, cex = leg.cex, bty = leg.bty, bg = leg.bg, title=legend.title) 
  }
  print(model)
}

