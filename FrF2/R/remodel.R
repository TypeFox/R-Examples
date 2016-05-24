`remodel` <-
function(obj){   
     ## obj is a linear model object
     ## check whether as many coefficients as terms plus intercept
     ## i.e. two-level factors only
     if (length(coef(obj)) > length(attr(terms(obj),"term.labels")) + 
         sign(attr(terms(obj),"intercept"))) 
         stop("Only 2-level-designs are covered by this routine.")
     mod <- obj$model
     term.ord <- attr(terms(obj),"order")
     nmain <- length(which(term.ord==1))
     intcol <- attr(attr(mod,"terms"),"intercept")
     respcol <- attr(attr(mod,"terms"),"response")
     respnam <- colnames(mod)[respcol]
     labs <- lapply(vector("list",nmain),function(sp){c("-","+")})
     xmod <- mod[,-c(intcol,respcol)] #[,rowSums(attr(attr(mod,"terms"),"factors"))>0]
       # numeric -1 1 columns instead of factors 
       # because otherwise prediction difficult for intermediate (0) level
  #   if (any(sapply(xmod,"is.factor")) || "aov" %in% class(obj)) {
       hilf <- 0
       for (i in 1:length(xmod)) {
           ## character variables are already factors because of lm object
           if (is.factor(xmod[[i]])) 
              labs[[i]] <- levels(xmod[[i]])
              else labs[[i]] <- range(xmod[[i]])
              xmod[[i]] <- as.numeric(xmod[[i]])
              mu <- mean(xmod[[i]])
              r <- range(xmod[[i]])
              r <- r[2]-r[1]
              xmod[[i]] <- 2*(xmod[[i]]-mu)/r
  #         }
           }
       mod[,-c(intcol,respcol)] <- xmod
       obj <- lm(terms(obj),data=data.frame(mod))
  #     }
       ## return list with linear model for further calculations
       ## and labels for annotating plots
       list(model=obj,labs=labs)
}

