predict.repolr <-
function(object, newdata = NULL, se.fit = FALSE, robust.var = TRUE,
        type = c("link", "response", "terms"), ...){
  type <- match.arg(type)
  if (missing(newdata) || is.null(newdata)) {
    mm <- model.matrix(object)
    pred.names <- paste(object$id, 1:object$max.id, sep=".")
  } else {
   nsubject.comp <- sum(table(newdata[,which(names(newdata) == object$subjects)])/length(object$times))
   if(nsubject.comp != round(nsubject.comp)){
    stop("newdata: need data.frame with complete times for all (unique) subjects")
   } else {
    pred.names <-  paste(rep(newdata[,which(names(newdata) == object$subjects)], each = length(object$times)),
                     1:object$max.id, sep=".")
    newdata[,which(names(newdata) == object$subjects)] <- rep(1:nsubject.comp, each = length(object$times))
   }
   resp.i <- match(as.character(terms(object$formula)[[2]]), names(newdata))
   if(is.na(resp.i) == TRUE){
    newdata <- data.frame(newdata, resp.var = rep(1, dim(newdata)[1]))
    names(newdata)[dim(newdata)[2]] <- as.character(terms(object$formula)[[2]])
   }
   exdata <- ord.expand(space = object$poly.mod$space, formula = object$orig.formula, 
                times = object$times, poly = object$poly.mod$poly, data = newdata, 
                subjects = object$subjects, categories = object$categories)
   mm <- model.matrix(exdata$formula, data = exdata$data)
  }
  if(type == "link"){
    pred <- as.numeric(coef(object) %*% t(mm))
    names(pred) <- pred.names
    if(se.fit == TRUE){
     se.funct <- function(y, object, Xmat, sel.vcov) {
            if(sel.vcov == TRUE){
              sqrt(Xmat[y,] %*% object$robust.var %*% Xmat[y,])
            } else {
              sqrt(Xmat[y,] %*% object$naive.var %*% Xmat[y,])
            }
            }
     se <- sapply(1:dim(mm)[1], se.funct, Xmat = mm, sel.vcov = robust.var,
               object = object, simplify = TRUE)
     names(se) <- pred.names
    }
  } else if (type == "response"){
    se.fit <- FALSE
    ilink.funct <- binomial(link="logit")$linkinv
    pred <- ilink.funct(as.numeric(coef(object) %*% t(mm)))
    names(pred) <- pred.names
   } else if (type == "terms"){  
    se.fit <- FALSE
    pred <- mm * matrix(rep(coef(object), times = dim(mm)[1]),
             nrow = dim(mm)[1], byrow = TRUE)
    rownames(pred) <- pred.names
  }
  if (se.fit == FALSE){
   pred
  } else {
   list(fit = pred, se.fit = se)
  }
}
