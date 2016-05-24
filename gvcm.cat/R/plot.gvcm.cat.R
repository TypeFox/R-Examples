plot.gvcm.cat <-
function(

x,
accuracy=2,
type="path",

individual=FALSE, # for type=="path", "coefs"
xlim,
ylim,
main = NULL,

indent = 0,
color = TRUE,
xscale = "lambda", # or "beta"; for type=="path"
label=TRUE, # for type=="path"; print lambdacv etc in figure?!
intercept=TRUE, # type=="coefs" - shall intercept be added to smooth functions, yes, no.
                # type=="path" - shall intercept be plottet or not
...
)


{

# check 
  # x
  if (!("gvcm.cat" %in% is(x)))
       stop ("x must be a 'gvcm.cat' object. \n")
  # type
  if (!(type %in% c("path", "score", "coefs"))) 
       stop ("type is incorrect. \n")
     
# colors
farben <- c('blue','lightblue','turquoise','lightgreen', 'darkblue',
'darkgreen','green','yellow','gold','orange','darkorange','red','darkred','darkviolet',
'violet','magenta','pink','grey','darkgrey','black')

colorvec <- c(farben, farben, rep('black',84))
    
# defaults wieder herstellen.
# if (exists("individual.paths")) individual <- individual.paths

if (x$method %in% c("AIC", "BIC")){ # forward selection

} else { # penalties...

 if (type=="path"){

    # check 
    if (!is.matrix(x$plot[[1]]))
    {stop ("Input argument 'plot' must be 'TRUE' for plotting coefficient paths. \n")}
    
    if (!is.logical(individual) && !is.character(individual) && !is.vector(individual))
    {stop ("Error in input argument 'individual'. \n")}

    if (length(as.vector(accuracy))!=1 || accuracy < 1 || accuracy > 4)
         stop ("accuracy must be a single integer > 0. \n")
         
    if (missing(xlim)) {xlim <- c(0,1)}
    if (!is.numeric(xlim) || !is.vector(xlim) || length(xlim)!=2 || min(xlim)<0 ) #|| max(xlim)>1 )
         stop ("Error in arguemnt xlim. \n") 
     xlim <- c(min(xlim), max(xlim))     
         
    if (!is.logical(color))     
     {stop ("Error in input argument 'color'. \n")}

    if (indent<0 || !is.numeric(indent)) {indent <- 0}
    if (indent>0.8) {indent <- 0.8}
    
    # definitions
    lambda <- x$tuning[[1]]
    assured.intercept <- x$control$assured.intercept
    f <- x$formula
    number.selectable.parameters <- x$number.selectable.parameters
    index1 <- x$indices[1,]
    index2 <- colSums(x$indices) - colSums(as.matrix(x$indices[c("index1", "index2b"),]))
    control <- x$control
    n.p <- length(index1)

    # prepare path
    path <- x$plot[[1]]
    path[,-1] <- round(path[,-1], digits=accuracy)

    i <- 2  
    while (i<dim(path)[2] && dim(path)[1]-length(reduce(path[,i],x$indices,assured.intercept)$beta) 
            < number.selectable.parameters) {i<-i+1}
    path <- path[,1:i]
    loml <- sum(abs(path[which(rep(index2, index1)!=0),1]))
    path <- path[,order(as.numeric(colnames(path)),decreasing=TRUE)]
    
    # other definitions
    lambdas <- as.numeric(colnames(path))
    lambda.upper <- lambdas[1]
        # lambda
        x. <- 1 - lambdas/lambda.upper
        s <- 1- (lambda/lambda.upper)
    if (xscale %in% c("beta", "beta.abs")) {
        x. <- colSums(abs(path)[which(rep(index2, index1)!=0),])/loml 
        s <- sum(abs(x$coefficients[which(rep(index2, index1)!=0)]))/loml
    } 
    if (color==TRUE){
        line.type <- - abs(rep(index2, times=index1)) + 2
        if (intercept==FALSE) colorvec[1] <- gray(1)
        colour <- rep(colorvec[1:n.p], times=index1) # rep(1:n.p, times=index1)
    } else {
        line.type <- rep(rep(c(2,5,1,4,6,1), length.out=length(index1)),times=index1) 
        colour  <- if (intercept) {rep(gray( (0:(n.p-1)) / n.p),times=index1)} else {rep(gray(c(1, (0:(n.p-2))/n.p)),times=index1)}
        # (0=blank, 1=solid (default), 2=dashed, 3=dotted, 4=dotdash, 5=longdash, 6=twodash) 
    }
    
    # zoooom
    path <- path[,min(which(x.>=min(xlim))):max(which(x.<=max(xlim)))]
    x. <- x.[min(which(x.>=min(xlim))):max(which(x.<=max(xlim)))]

    # plot options
#    lambda.max.print = TRUE
#    lambda.opt.plot = TRUE
#    lambda.opt.print = TRUE
    
    anteil.b <- 1/(3*22)
    range.x <- max(xlim)-min(xlim)
    coefs.left <- 4*anteil.b*range.x
    
    # remove NA's (e.g. if ML-estimate does not exisit)
    if (sum(is.na(path))>0) 
       warning("Some estimates contain NAs; set to zero for plotting. \n")
    for (k in 1:ncol(path)) path[which(is.na(path[,k])),k] <- 0 
     
    # individual=TRUE    
    if (!is.logical(individual) || individual==TRUE) {
        lambda.max.print <- FALSE
        lambda.opt.print <- FALSE
        e <- c(1,cumsum(index1)[1:(n.p-1)]+1)
        index <- matrix(c(e, e+1, cumsum(index1)),nrow=n.p,ncol=3,byrow=FALSE)
        if (!is.logical(individual)) {
            tf <- terms(f,specials="v", "p")
            arg <- strsplit(gsub(" ", "", x$formula[length(x$formula)]), "\\+")  # intercept
            w <- if (grepl("v\\(1",arg)==TRUE) 1 else 0                          # intercept
            coefs.names <- if (w==0) c("1")  else c()                            # 
            coefs.names <- c(coefs.names, gsub(" ", "", as.character(attr(tf,"term.labels"))))
            individual <- gsub(" ", "", individual)
            which.coefs <- c()
            for (i in 1:length(coefs.names)) {
               if (coefs.names[i] %in% individual) {
                  which.coefs <- c(which.coefs, i)
               }
            }
        } else { 
            which.coefs <- 1:nrow(index)
        }
        if(length(which.coefs)<4){
           mfrows<-c(1,length(which.coefs))}else{mfrows<-c(2,ceiling(length(which.coefs)/2))
        }
        par(mfrow=mfrows)
        main <- ""# unlist(strsplit(as.character(f[3]),"\\+"))
    } else {
        index <- matrix(c(1, 2, sum(index1)),ncol=3, nrow=1)
        which.coefs <- 1:nrow(index)
    }
        
       
    for (j in which.coefs){
    
        if (!missing(ylim)) {
             if (!is.numeric(ylim) || !is.vector(ylim) || length(xlim)!=2)
                 stop ("Error in arguemnt ylim. \n")  
             ylimj <- c(min(ylim), max(ylim))
           } else {
             ylimj <- c(floor(min(path[index[j,1]:index[j,3],])), 
             ceiling(max(path[index[j,1]:index[j,3],])))
           } 
        heigth <- ylimj[2] - ylimj[1]
    
        # legend
          legend.y <- path[index[j,1],dim(path)[2]]+0.05*heigth
          if (index[j,1] < index[j,3]) {
            for (i in index[j,2]:index[j,3]) {
                    legend.y <- c(legend.y, path[i,dim(path)[2]]+0.05*heigth)
                }
          }
          # legend.x
            legend.x <- rep(1-coefs.left, times=length(legend.y))
            if (length(legend.y)>1) {
                legend.y.ordered <- legend.y[order(legend.y,decreasing=TRUE)]
                diff.legend <- (legend.y.ordered[1:(length(legend.y)-1)] - legend.y.ordered[2:(length(legend.y))])/heigth
                for (i in 1:length(diff.legend)){
                    if (diff.legend[i] <0.02180){
                        if(legend.x[order(legend.y,decreasing=TRUE)[i]] < 1-coefs.left+2*indent){
                            legend.x[order(legend.y,decreasing=TRUE)[i+1]] <- legend.x[order(legend.y,decreasing=TRUE)[i]] + indent
                        } else {
                            legend.x[order(legend.y,decreasing=TRUE)[i+1]] <- legend.x[order(legend.y,decreasing=TRUE)[i]] - 2*indent            
                        }
                    }
                }
            }
          # correct margin
          anteil <- anteil.b * max(nchar(names(x$coefficients)[[1]])) +  (indent!=0)*.9*(max(legend.x) - 1)*(1-indent) + (1-indent!=0)*.05
          right.margin <- anteil*range.x/(1-anteil)
#          print( right.margin)
          const. <- 0
          
        # plot
        xlab. <- expression(paste(1-lambda/lambda[max]))
        if (xscale=="lambda.abs") xlab. <- expression(paste(lambda))
        if (xscale=="beta") xlab. <- expression(paste(abs(beta)/abs(beta[max])))
        if (xscale=="beta.abs") xlab. <- expression(paste(abs(beta)))
        
        matplot(x.,path[index[j,1],], type="l", lty=line.type[index[j,1]],    # lty=1 => solid
                col = colour[index[j,1]], lwd=2,                              # lty=2 => dashed
                main = main[j], axes = FALSE,                                 # lty=3 => punkte
                xlim = c(min(xlim),max(xlim) + right.margin),
                ylim = ylimj,
                xlab= xlab., 
                ylab= "estimated coefficient"
                )
        legend(x=legend.x[1]+const., y=legend.y[1],
               legend=rownames(path)[index[j,1]],box.lty=0, text.col=colour[index[j,1]])
        if (index[j,1] < index[j,3]) {
            for (i in index[j,2]:index[j,3]) {
                matlines(x.,path[i,], type="l", lty= line.type[i] ,col= colour[i],  lwd=2)
                legend(x=legend.x[i-index[j,2]+2]+const., y=legend.y[i-index[j,2]+2],
                       legend=rownames(path)[i],box.lty=0, text.col=colour[i]) 
            } 
        }

        
        # x$lambda => dotted line
        if (label==TRUE){
              if(s>0){
                matlines(x=c(s, s), y=ylimj, type="l", lty=3, lwd=1)
              } else {s<-1} 
            }
        if (label==TRUE){
            legend("bottomleft", legend=
#                    c(expression(paste(lambda[CV], " = ")), (round(lambda, digits=2))) , 
                    bquote(lambda[CV] == .(eval(round(lambda, digits=2)))),
                    box.lty=0, horiz=TRUE)
            }
        
        # max.lambda
        if (label==TRUE){
            legend("topleft", legend=
#                   c(expression(paste(lambda[max], " = ")),(round(lambda.upper, digits=2))), 
                   bquote(lambda[max] == .(eval(round(lambda.upper, digits=2)))),
                   box.lty=0, horiz=TRUE)
            }
            
        # the good look
        xlim. <- xlim
        if (xscale=="beta.abs") xlim. <- xlim*loml
        if (xscale=="lambda.abs") xlim. <- (1-xlim)*lambda.upper
        
        box()
        axis(1, at=seq(from=min(xlim),to=max(xlim),length.out=6), 
                labels=as.character(seq(from=xlim.[1], to=xlim.[2],length.out=6)), 
                col = "black", lty = 1, lwd = 1 )          
        axis(2, #at=seq(from = ylimj[1], to = ylimj[2], length.out=heigth),
                #labels=as.character(seq(from = ylimj[1], to = ylimj[2], length.out=heigth)),
                col = "black", lty = 1, lwd = 1 )    
        
        }

}
# type score ###################################################################
if (type=="score"){
    # check x$plot
    if (is.numeric(x$plot[[2]])==FALSE){stop ("type='score' requires cross-validation of lambda. \n")}
    
    # defintions
    score <- x$plot[[2]]
#    ph <- as.numeric(rownames(score))
#    phi <- x$tuning[[2]]
    lambda <- x$tuning[[1]]
    
    lambdas <- if (is.matrix(score)) {
                as.numeric(colnames(score))
               } else {
                as.numeric(names(score))
               }
    l <- lambdas[order(lambdas)]
    z <- as.vector(score)[order(lambdas)]
       
    if(missing(xlim)){xlim <- c(floor(min(l)),ceiling(max(l)))}
    if (!is.numeric(xlim) || !is.vector(xlim) || length(xlim)!=2)
         stop ("Error in arguemnt xlim. \n")  
    if(missing(ylim)){ylim <- c(floor(min(z)),ceiling(max(z)))}
    if (!is.numeric(ylim) || !is.vector(ylim) || length(ylim)!=2)
         stop ("Error in arguemnt ylim. \n")  
    
    matplot(l, z, type="l", lty=1, 
            col = 1,                         
            main = main, axes = FALSE, 
            xlim = xlim, 
            ylim = ylim,
            xlab= expression(paste(lambda)), ylab= "score"
            )
    
    matlines(x=c(lambda, lambda), y=ylim, type="l", lty=3, lwd=1)
    
    # phi flexible => add rest    
#    if (length(ph)>1) {
#        j <- 2
#        for (i in rest) {
#            z <- score[i,]
#            matlines(l, z, type="l", lty= 1 ,col=j)
#            j <- j+1
#           }
#        label <- factor(1:length(ph), levels= 1:length(ph), labels = as.character(ph[c(best,rest)]) )
#        legend("bottomright", levels(label), fill=1:length(ph))    
#       }
       
    # good look
    box()
    axis(1, at=seq(from = floor(min(l)), to = ceiling(max(l)), by = 1),
            labels=as.character(seq(from = min(xlim), to = max(xlim), by = 1)),
            col = "black", lty = 1, lwd = 1 )    
    axis(2, #at=seq(from = floor(min(z)), to = ceiling(max(z)), by = 1),
            #labels=as.character(seq(from = min(ylim), to = max(ylim), by = 1)),
            col = "black", lty = 1, lwd = 1 )    
    
} 

# type coefs ###################################################################
if (type=="coefs"){

plotstair <- function(y, x, coef.matrix, coef.matrix.oml, xlim, xlab, color, main) {
  nlevel <- length(levels(x)) # incl. reference category
  plot(0:nlevel, c(0, coef.matrix.oml, rev(coef.matrix.oml)[1]),
       type="s", col=gray(.3), lwd=2, lty=ifelse(color,1,2),
       axes = FALSE, xlim=c(0,nlevel),  
       ylim=c(min(0, coef.matrix, coef.matrix.oml), max(0, coef.matrix, coef.matrix.oml)),
       xlab=paste("Levels of ", xlab, sep=""), ylab="estimated coefficients", main=main)
  matlines(0:nlevel, c(0, coef.matrix, rev(coef.matrix)[1]),
       type="s", col=ifelse(color, colorvec[2], 1), lty=1, lwd = 2)
  box()
  axis(1, at=seq(from = .5, to = nlevel-.5, by = 1),
          labels=levels(x),
          col = "black", lty = 1, lwd = 1 )
  axis(2, col = "black", lty = 1, lwd = 1 )
  legend("bottomleft", col=c(gray(.3), ifelse(color, colorvec[2], 1)), legend=c("ML", "pen."), bty="n", lwd=2, lty=c(ifelse(color,1,2),1))
}

plotstair. <- function(y, x, coef.matrix, coef.matrix.oml, xlim, xlab, color, main) {
  nlevel <- length(levels(x)) # incl. reference category
  plot(0:(nlevel), c(coef.matrix.oml, rev(coef.matrix.oml)[1]),
       type="s", col=gray(.3), lwd=2, lty=ifelse(color,1,2),
       axes = FALSE, xlim=c(0,nlevel),
       ylim=c(min(0, coef.matrix, coef.matrix.oml), max(0, coef.matrix, coef.matrix.oml)),
       xlab=paste("Levels of ", xlab, sep=""), ylab="estimated coefficients", main=main)
  matlines(0:nlevel, c(coef.matrix, rev(coef.matrix)[1]),
       type="s", col=ifelse(color, colorvec[2], 1), lty=1, lwd = 2)
  box()
  axis(1, at=seq(from = .5, to = nlevel-.5, by = 1),
          labels=levels(x),
          col = "black", lty = 1, lwd = 1 )
  axis(2, col = "black", lty = 1, lwd = 1 )
  legend("bottomleft", col=c(gray(.3), ifelse(color, colorvec[2], 1)), legend=c("ML", "pen."), bty="n", lwd=2, lty=c(ifelse(color,1,2),1))
}

plotlines <- function(y, x, u, coef.matrix, int.matrix, xlim, ylim, xlab, ylab, color, main) {
#  plot(x, y, col=gray(.4), xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, main=main)
  plot(0, 0, col=gray(1), xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, main=main)
  farbe <- if(color) colorvec[1:length(coef.matrix)] else gray( (0:length(coef.matrix)) / length(coef.matrix))
  linie <- if(color) rep(1, length(coef.matrix)) else 1:length(coef.matrix)
  for (i in 1:length(coef.matrix)) {
       abline(int.matrix[i], coef.matrix[i], col=farbe[i], lwd=2, lty=linie[i])
  }
  legend("bottomleft", legend=levels(u), col=farbe, bty="n", lwd=2, lty=linie)
}

plotsp <- function(y, x, knots, coef.matrix, int.matrix, xlim, ylim, xlab, ylab, color, main, intercept) {

  farbe <- if(color) colorvec[1:length(int.matrix)] else gray( (0:length(int.matrix)) / length(int.matrix))
  linie <- if(color) rep(1, length(int.matrix)) else 1:length(int.matrix)

  # ylim
  if (any(is.na(ylim))) {
      ymax <- c()
      ymin <- c()
      for (i in 1:length(int.matrix)) {
           ymax <- c(ymax, max(sp(x, knots=knots) %*% coef.matrix + ifelse(intercept, int.matrix[i], 0)))
           ymin <- c(ymin, min(sp(x, knots=knots) %*% coef.matrix + ifelse(intercept, int.matrix[i], 0)))
           }
      ylim <- c(min(ymin), max(ymax))
      }
  # data points
  plot(jitter(x), rep(min(ylim), length(x)), col=gray(.4), xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, pch=3, main=main)
  
  # abline knots
  mini <- min(x)
  maxi <- max(x)
  dis <- (maxi - mini)/(knots - 7) # 5, by = ((to - from)/(length.out - 1))
  knoten <- seq(from=mini-3*dis, to=maxi+3*dis, length.out=knots)
  for(i in 1:length(knoten)) abline(v=knoten[i], lty=3, col=gray(.4))
  
  # plot function
  for (i in 1:length(int.matrix)) {
       lines(x[order(x)], (sp(x, knots=knots) %*% coef.matrix +
              ifelse(intercept, int.matrix[i], 0))[order(x)], lwd=2, col=farbe[i], lty=linie[i])
  }
  if (length(int.matrix) > 1 && intercept) {
      spaltennamen <- if (is.null(colnames(int.matrix))) 1:length(int.matrix) else colnames(int.matrix)
      legend("bottomleft", legend=spaltennamen,
             lwd=2, col = farbe, bty="n", lty=linie)
  }
}

plotvspline <- function(y, x, u, knots, coef.matrix, int.matrix, xlim, ylim, xlab, ylab, color, main, intercept) {

  farbe <- if(color) colorvec[1:length(levels(u))] else gray( (0:length(levels(u))) / length(levels(u)))
  linie <- if(color) rep(1, length(levels(u))) else 1:length(levels(u))

  # ylim
  if (any(is.na(ylim))) {
      ymax <- c()
      ymin <- c()
      for (i in 1:length(int.matrix)) {
           ymax <- c(ymax, max(sp(x, knots=knots) %*% coef.matrix + ifelse(intercept, int.matrix[i], 0)))
           ymin <- c(ymin, min(sp(x, knots=knots) %*% coef.matrix + ifelse(intercept, int.matrix[i], 0)))
           }
      ylim <- c(min(ymin), max(ymax))
      }
  # data points
  plot(jitter(x), rep(min(ylim), length(x)), col=gray(.4), xlab=xlab, ylab=ylab, xlim=xlim, ylim=ylim, pch=3, main=main)

  # abline knots
  mini <- min(x)
  maxi <- max(x)
  dis <- (maxi - mini)/(knots - 7) # 5, by = ((to - from)/(length.out - 1))
  knoten <- seq(from=mini-3*dis, to=maxi+3*dis, length.out=knots)
  for(i in 1:length(knoten)) abline(v=knoten[i], lty=3, col=gray(.4))

  # plot functions
  des <- vspline(x, u, knots=knots)
  for (i in 1:length(levels(u))) {
       level <- levels(u)[i]
       xu <- x[which(u==level)]
       von <- nrow(coef.matrix)*(i - 1) + 1
       bis <- nrow(coef.matrix)* i
       lines(xu[order(xu)], (des[which(u==level),  von:bis] %*% coef.matrix[,i] +
            ifelse(length(int.matrix)>1, ifelse(intercept, int.matrix[i], 0), ifelse(intercept, int.matrix[1], 0)))[order(xu)],
            col=farbe[i], lwd=2, lty=linie[i])

  }
  legend("topright", legend=levels(u), lwd=2, col = farbe, bty="n", lty=linie)
}
 
     # definitions
     formula <- x$formula
     family <- x$family
     link <- family$linkfun
     indices <- x$indices
     indexs <- as.matrix(indices[c("index1", "index2", "index3", "index4", "index5", "index6", "index7", "index8", "index9"), ])
     data <- x$data
     Terms <- x$terms
     label <- attr(Terms,"term.labels")
     variables <- as.character(attr(Terms,"variables"))[-c(1)]
     ncoefs <- ncol(indices)
     special <- c("v", "p", "grouped", "grouped.fused", "sp", "SCAD", "elastic", "vspline")
     m <- model.frame(formula=terms(formula, specials=special, data=data), data)

     r <- ifelse(indices[1,1]!=1, 0, -1)  # korrektur beziehung indices[,i], label[i]
           # falls v(1,u): r = 0 <=> indices[,i]   = label[i]
           # falls 1     : r = -1 <=> indices[,i] = label[i+r]


     # plotables
     plotables <- c()
     for (i in 1:ncol(indices)){
          variable <- which(variables==label[i+r])
#          if (indices["index2",i]!=0 && i>1) plotables <- c(plotables, i) # v
          if (indices["index2",i]!=0) plotables <- c(plotables, i) # v
          if (indices["index5",i]!=0) plotables <- c(plotables, i) # grouped.fused
          if (indices["index6",i]!=0) plotables <- c(plotables, i) # sp
          if (indices["index9",i]!=0) plotables <- c(plotables, i) # vspline
          if (indices["index3",i]!=0){                             # p
              u <- model.frame(~eval(parse(text=sub(")", "", gsub(" ", "",
                   strsplit(strsplit(label[i+r], "p\\(")[[1]][2], ",")[[1]][1])))), data)
              erms <- attr(u, "terms")
              if (attr(erms, "dataClasses")[[1]]!="numeric"){plotables <- c(plotables, i)}
          }
          if (indices["index4",i]!=0){                             # grouped
              u <- model.frame(~eval(parse(text=sub(")", "", gsub(" ", "",
                   strsplit(label[i+r], "grouped\\(")[[1]][2])))), data)
              erms <- attr(u, "terms")
              if (attr(erms, "dataClasses")[[1]]!="numeric"){plotables <- c(plotables, i)}
          }
          if (indices["index7",i]!=0){                             # SCAD
              u <- model.frame(~eval(parse(text=sub(")", "", gsub(" ", "",
                   strsplit(strsplit(label[i+r], "SCAD\\(")[[1]][2], ",")[[1]][1])))), data)
              erms <- attr(u, "terms")
              if (attr(erms, "dataClasses")[[1]]!="numeric"){plotables <- c(plotables, i)}
          }
          if (indices["index8",i]!=0){                             # elastic
              u <- model.frame(~eval(parse(text=sub(")", "", gsub(" ", "",
                   strsplit(strsplit(label[i+r], "elastic\\(")[[1]][2], ",")[[1]][1])))), data)
              erms <- attr(u, "terms")
              if (attr(erms, "dataClasses")[[1]]!="numeric"){plotables <- c(plotables, i)}
          }
     }

    # individual    
    if (!is.logical(individual)) {

        #tf <- terms(formula,specials="v", "p")
        arg <- strsplit(gsub(" ", "", formula[length(formula)]), "\\+")  # intercept
        coefs.names <- if (r==-1) c("1")  else c()                            # 
        coefs.names <- c(coefs.names, gsub(" ", "", as.character(attr(Terms,"term.labels"))))
        individual <- gsub(" ", "", individual)
        which.coefs <- c()
        for (i in 1:length(coefs.names)) {
           if (coefs.names[i] %in% individual) {
              which.coefs <- c(which.coefs, i)
           }
        }
        plotables <- which.coefs[which(which.coefs %in% plotables)]        
    } 

     # par
     main <- if (length(plotables)>1) "" else main
     mfrows <- if (length(plotables)<4) c(1, length(plotables)) else c(2,ceiling(length(plotables)/2))
     par(mfrow=mfrows)

     # loop
     if (length(plotables)!=0) {
     
     y <- model.extract(m, "response")
     if (is.factor(y)==TRUE){y <- as.numeric(y)-1}
     if (!is.null(dim(y)[2]) && family$family=="binomial") {
          y <- y[,1]/(y[,1]+y[,2])
     }
     if (family$family=="poisson") y[y==0] <- y[y==0] + .01
     if (family$family=="binomial") {
         y[y==0] <- .001
         y[y==1] <- .999
         }
     y <- link(y)

     for (i in plotables) { # plotables bzgl. indices
     
          typespecial <- which(indexs[,i]!=0)
          typespecial <- typespecial[-1]
          vonbis <- (sum(indexs["index1",1:i])-indexs["index1",i]+1):sum(indexs["index1",1:i])
          ki <- indexs["index1",i]

          if (typespecial %in% c(3, 4, 7, 8)) { # treppe
              # p (3), grouped (4), SCAD (7), elastic (8) (intercept mit v wird nicht geplottet!!)
              if (typespecial==3) X <- data[, which(colnames(data) == sub(")", "", gsub(" ", "", strsplit(strsplit(label[i+r], "p\\(")[[1]][2], ",")[[1]][1])))] # p, group, SCAD, elastic => ref kategorie
              if (typespecial==4) X <- data[, which(colnames(data) == sub(")", "", gsub(" ", "", strsplit(strsplit(label[i+r], "grouped\\(")[[1]][2], ",")[[1]][1])))] # p, group, SCAD, elastic => ref kategorie
              if (typespecial==7) X <- data[, which(colnames(data) == sub(")", "", gsub(" ", "", strsplit(strsplit(label[i+r], "SCAD\\(")[[1]][2], ",")[[1]][1])))] # p, group, SCAD, elastic => ref kategorie
              if (typespecial==8) X <- data[, which(colnames(data) == sub(")", "", gsub(" ", "", strsplit(strsplit(label[i+r], "elastic\\(")[[1]][2], ",")[[1]][1])))] # p, group, SCAD, elastic => ref kategorie
              nlevel <- length(levels(u)) # incl reference category

              coef.matrix <- x$coefficients[vonbis]
              coef.matrix.oml <- x$coefficients.oml[vonbis]
              
              xlimj <- c(0,nlevel)

              plotstair(y, X, coef.matrix, coef.matrix.oml, xlimj, label[i+r], color, main)

          }

          if (typespecial %in% c(2, 5) && i > 1) { # mehrere geraden
              # v ohne intercept (2), grouped.fused (5)
              u <- data[, which(colnames(data) == sub(")", "", gsub(" ", "", strsplit(label[i+r], "\\,")[[1]][2])) )]
              X <- data[, which(colnames(data) == sub("grouped.fused(", "", sub("v(", "", sub("x=", "", gsub(" ", "", strsplit(label[i+r], "\\,")[[1]][1])), fixed=TRUE), fixed=TRUE))]
              coef.matrix <- x$coefficients[vonbis]
              index.u <- which(colnames(data) == sub(")", "", gsub(" ", "", strsplit(label[i+r], "\\,")[[1]][2])))    ###
              int.matrix <- if (indices[1,1]>1) {
                   index.u.int <- which(colnames(data) == gsub(")", "", gsub(" ", "", strsplit(names(m)[2], "\\,")[[1]][2])))
                   if (index.u == index.u.int) {
                       matrix(x$coefficients[1:indices[1,1]], nrow=1)
                       } else {
                       matrix(x$coefficients[1], nrow=1)
                       warning("Plot assumes x$coefficients[1] as intercept.", call. = FALSE)
                       }
                   } else {matrix(x$coefficients[1], nrow=1)}
              if (!intercept) {int.matrix <- matrix(0, ncol=levels(u), nrow=1)}

              # X numeric
              if (!is.factor(X)){
              if (!missing(xlim)) {
                     if (!is.numeric(xlim) || !is.vector(xlim) || length(xlim)!=2)
                         stop ("Error in arguemnt xlim. \n")
                     xlimj <- c(min(xlim), max(xlim))
              } else {
                     xlimj <- c(floor(min(X)), ceiling(max(X)))
              }

              if (!missing(ylim)) {
                     if (!is.numeric(ylim) || !is.vector(ylim) || length(xlim)!=2)
                         stop ("Error in arguemnt ylim. \n")
                     ylimj <- c(min(ylim), max(ylim))
              } else {
                     can <- c(min(coef.matrix)*min(X), min(coef.matrix)*max(X), max(coef.matrix)*min(X), max(coef.matrix)*max(X))
                     can <- c(can + min(int.matrix), can + max(int.matrix))
                     ylimj <- c(floor(min(can)), ceiling(max(can)))
              }
              xlab <- gsub("v(", "", gsub("grouped.fused(", "", gsub(")", "", label[i+r], fixed=TRUE), fixed=TRUE), fixed=TRUE)
              plotlines (y, X, u, coef.matrix, int.matrix, xlimj, ylimj, xlab, "predictor", color, main)
              } else {
              # X binary
              nlevel <- length(levels(u)) # incl reference category
              coef.matrix.oml <- x$coefficients.oml[vonbis]
              xlimj <- c(-1,1)
              plotstair.(y, u, coef.matrix, coef.matrix.oml, xlimj, label[i+r], color, main)
              }
          }
         if (typespecial %in% c(2, 5) && i==1) { # var intercept
              u <- data[, which(colnames(data) == sub(")", "", gsub(" ", "", strsplit(label[i+r], "\\,")[[1]][2])) )]
              coef.matrix <- x$coefficients[vonbis]
              nlevel <- length(levels(u)) # incl reference category
              coef.matrix.oml <- x$coefficients.oml[vonbis]
              xlimj <- c(-1,1)
              plotstair.(y, u, coef.matrix, coef.matrix.oml, xlimj, label[i+r], color, main)
          }

          if (typespecial == 6) { # sp (6)

              # variable <- which(variables==label[i+r])
              X <- data[, which(colnames(data) == sub(")", "", sub("sp(", "", sub("x=", "", gsub(" ", "", strsplit(label[i+r], "\\,")[[1]][1])), fixed=TRUE), fixed=TRUE))]
              knots <- ki + 4 + 1
              coef.matrix <- matrix(x$coefficients[vonbis], ncol=1, byrow=FALSE)
              int.matrix <- matrix(x$coefficients[1:indices[1,1]], nrow=1)

              if (!missing(xlim)) {
                     if (!is.numeric(xlim) || !is.vector(xlim) || length(xlim)!=2)
                         stop ("Error in arguemnt ylim. \n")
                     xlimj <- c(min(xlim), max(xlim))
              } else {
                     xlimj <- c(floor(min(X)), ceiling(max(X)))
              }

              if (!missing(ylim)) {
                     if (!is.numeric(ylim) || !is.vector(ylim) || length(ylim)!=2)
                         stop ("Error in arguemnt ylim. \n")
                     ylimj <- c(min(ylim), max(ylim))
              } else {
                     ylimj <- NA
              }
              xlab <- gsub("sp(", "", gsub(")", "", label[i+r], fixed=TRUE), fixed=TRUE)
              plotsp(y, X, knots, coef.matrix, int.matrix, xlimj, ylimj, xlab, label[i+r], color, main, intercept)

          }

          if (typespecial == 9) {  # vspline (9)

              # variable <- which(variables==label[i+r])
              levels.u <- indices["index2b",i] # index2b = number levels of u;
              k.spline <- ki/levels.u # k.spline = number of coefficients per spline and level
              X <- data[, which(colnames(data) == sub("vspline(", "", sub("x=", "", gsub(" ", "", strsplit(label[i+r], "\\,")[[1]][1])), fixed=TRUE))]
              index.u <- which(colnames(data) == sub(")", "", gsub(" ", "", strsplit(label[i+r], "\\,")[[1]][2])))
              u <- data[, index.u]
              knots <- k.spline + 4 + 1
              coef.matrix <- matrix(x$coefficients[vonbis], ncol=levels.u, nrow=k.spline, byrow=FALSE)
              int.matrix <- if (indices[1,1]>1) {
                   index.u.int <- which(colnames(data) == gsub(")", "", gsub(" ", "", strsplit(names(m)[2], "\\,")[[1]][2])))
                   if (index.u == index.u.int) {
                       matrix(x$coefficients[1:indices[1,1]], nrow=1)
                       } else {
                       matrix(x$coefficients[1], nrow=1)
                       warning("Plot assumes x$coefficients[1] as intercept.", call. = FALSE)
                       }
                   } else {matrix(x$coefficients[1:indices[1,1]], nrow=1)}

              if (!missing(xlim)) {
                     if (!is.numeric(xlim) || !is.vector(xlim) || length(xlim)!=2)
                         stop ("Error in arguemnt xlim. \n")
                     xlimj <- c(min(xlim), max(xlim))
              } else {
                     xlimj <- c(floor(min(X)), ceiling(max(X)))
              }

              if (!missing(ylim)) {
                     if (!is.numeric(ylim) || !is.vector(ylim) || length(ylim)!=2)
                         stop ("Error in arguemnt ylim. \n")
                     ylimj <- c(min(ylim), max(ylim))
              } else {
                     ylimj <- NA
              }
              xlab <- gsub("vspline(", "", gsub(")", "", label[i+r], fixed=TRUE), fixed=TRUE)
              plotvspline(y, X, u, knots, coef.matrix, int.matrix, xlimj, ylimj, xlab, label[i+r], color, main, intercept)
          }

     }
     }

  } # type="coefs"

} # method lqa

}


