NCovered.io <- function(par,model,data,group=TRUE,...){
  # see docs in NCovered
  #
  # computes abundance in covered region for io model object
  #
  # par   - parameter values (used when computing derivatives wrt 
  #         parameter uncertainty)
  # model - ddf model object
  # group - if TRUE computes group abundance and if FALSE individual abundance
  #
  # result - abundance estimate
  #
  # Functions Used: predict (predict.io), compute.Nht 

  # set pars if we are differentiating
  # then extract fitted values
  if(!is.null(par)){
    model$mr$mr$coefficients <- par[1:length(model$mr$mr$coefficients)]
    model$ds$par <- par[(length(model$mr$mr$coefficients)+1):length(par)]
    model$ds$ds$aux$ddfobj <- assign.par(model$ds$ds$aux$ddfobj,model$ds$par)
    fitted <- predict(model,compute=TRUE)$fitted
  }else{
    fitted <- model$fitted
  }

  if(!group){
    size <- model$data$size[model$data$observer==1&model$data$object %in%
                            as.numeric(names(model$fitted))]
    Nhat <- sum(compute.Nht(fitted,FALSE,size))
  }else{
    Nhat <- sum(compute.Nht(fitted,TRUE,size=NULL))
  }

  return(Nhat)
}
