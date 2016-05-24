make.mon.ui <- function(x, type = "coeff", contr = NULL){
    if (!type %in% c("coeff","mean")) stop("type must be either coeff or mean")
    if (type=="coeff"){
      if (!is.factor(x)) stop("For type coeff, x must be a factor.")
      if (is.null(contr)){ 
        if (!is.character(attr(x, "contrasts"))) contr <- "contr.treatment"
        else contr <- attr(x, "contrasts")
        }
      if (!contr %in% c("contr.treatment", "contr.SAS", "contr.diff", "contr.sum"))
        stop("invalid value for contr")
      levels <- levels(x)
      lenglev <- length(levels)

      if (!all(contrasts(x)==eval(parse(text=paste(contr, "(levels)", sep="")))))
                 stop("x is not coded with the chosen contrast type.")
      if (contr=="contr.treatment"){
          ui <- diag(1, lenglev-1)
          for (i in 2:(lenglev-1))
             ui[i,i-1] <- -1
      }
      if (contr=="contr.SAS"){
          ui <- diag(-1, lenglev-1)
          for (i in 2:(lenglev-1))
             ui[i-1,i] <- 1
      }
      if (contr=="contr.diff"){
          ui <- diag(1, lenglev-1)
      }
      if (contr=="contr.sum"){
          ui <- diag(-1, lenglev-1)
          for (i in 2:(lenglev-1))
             ui[i-1,i] <- 1
          ui[lenglev-1,] <- rep(-1,lenglev-1)
          ui[lenglev-1,lenglev-1] <- -2
      }
    }
    else{
        if (!(is.numeric(x))) stop("For type mean, x must be the dimension of the mean vector.")
        if (!length(x)==1) stop("For type mean, x must be the dimension of the mean vector.")
        if (!x%%1==0) stop(stop("x must be the dimension of the mean vector."))
        if (x<=1) stop("x must be larger than 1.")
        ui <- diag(1,x)
          for (i in 2:x)
             ui[i,i-1] <- -1
        ui <- ui[-1,]
    }
    ui
}
