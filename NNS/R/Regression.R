#' VN Regression
#'
#' Generates a nonlinear regression based on partial moment quadrant means.
#'
#' @param x Independent Variable
#' @param y Dependent Variable
#' @param order Controls the number of partial moment quadrant means.  Defaults to smaller order to avoid overfitting
#' @param point.est Returns the fitted value for any value of the independent variable
#' @param location Sets the legend location within the plot
#' @param print.values Defaults to FALSE, set to TRUE in order to return all fitted values for independent variable
#' @param print.equation Defaults to FALSE, set to TRUE in order to return the local coefficients
#' @keywords nonlinear regression
#' @author Fred Viole, OVVO Financial Systems
#' @references Viole, F. and Nawrocki, D. (2013) "Nonlinear Nonparametric Statistics: Using Partial Moments"
#' \url{http://amzn.com/1490523995}
#' @examples
#' set.seed(123)
#' x<-rnorm(100); y<-rnorm(100)
#' \dontrun{VN.reg(x,y)}
#' @export


VN.reg = function (x, y,
            order=max(2,ceiling(VN.dep.reg(x,y,2)[2]*ceiling(log10(length(x))))),
            point.est = NULL,
            location = 'top',
            print.values = FALSE,
            print.equation = FALSE){

  temp_df = data.frame(x=x, y=y)
  temp_df[,'temp_part'] = 'p'
  temp_df[,'master_part'] = 'p'

  regression.points = data.frame(matrix(ncol = 2))
  Regression.Coefficients = data.frame(matrix(ncol=3))

  names(Regression.Coefficients) = c('Coefficient','X Lower Range','X Upper Range')



  if(order<1){return("Please Increase the Order Specification")}
  order=order-1

    for(i in 0:(order)){

      for(item in unique(temp_df$master_part)){
        tmp_xbar = mean(temp_df[temp_df$master_part == item,'x'])
        tmp_ybar = mean(temp_df[temp_df$master_part == item, 'y'])



        temp_df[temp_df$x >= tmp_xbar & temp_df$y >= tmp_ybar & temp_df$master_part == item,'temp_part'] = paste(temp_df[temp_df$x >= tmp_xbar & temp_df$y >= tmp_ybar & temp_df$master_part == item,'master_part'], 1, sep = '')
        temp_df[temp_df$x < tmp_xbar & temp_df$y >= tmp_ybar & temp_df$master_part == item,'temp_part'] = paste(temp_df[temp_df$x < tmp_xbar & temp_df$y >= tmp_ybar & temp_df$master_part == item,'master_part'], 2, sep = '')
        temp_df[temp_df$x >= tmp_xbar & temp_df$y < tmp_ybar & temp_df$master_part == item,'temp_part'] = paste(temp_df[temp_df$x >= tmp_xbar & temp_df$y < tmp_ybar & temp_df$master_part == item,'master_part'], 3, sep = '')
        temp_df[temp_df$x < tmp_xbar & temp_df$y < tmp_ybar & temp_df$master_part == item,'temp_part'] = paste(temp_df[temp_df$x < tmp_xbar & temp_df$y < tmp_ybar & temp_df$master_part == item,'master_part'], 4, sep = '')


        ### order + 1 to account for 'p'
        if(nchar(item)==order+1){

        regression.points[item,] = cbind(tmp_xbar,tmp_ybar)

        }


      }

      temp_df[,'master_part'] = temp_df[, 'temp_part']


    #}


  }

  min.range = min(na.omit(regression.points[,1]))
  max.range = max(na.omit(regression.points[,1]))


  Dynamic.average.min = mean(y[x<min.range])
  Dynamic.average.max = mean(y[x>max.range])

  ###Endpoints
  if(length(x[x<min.range])>0){
    if(VN.dep.reg(x,y,2)[2]<.5){
      x0 = Dynamic.average.min} else {
        x0 = y[x==min(x)]} }  else {x0 = y[x==min(x)]}

  if(length(x[x>max.range])>0){
    if(VN.dep.reg(x,y,2)[2]<.5){x.max = Dynamic.average.max} else {x.max = y[x==max(x)]}}  else { x.max = y[x==max(x)]}


  regression.points[1,2] = x0
  regression.points[1,1] = min(x)

  regression.points[length(regression.points[,2])+1,2] = x.max
  regression.points[length(regression.points[,1]),1] = max(x)



  ###Regression Equation

  regression.points = na.omit(regression.points[order(regression.points),])


  q=length(regression.points[,1])



  for(i in 1:q){

      rise = regression.points[i+1,2] - regression.points[i,2]
      run = regression.points[i+1,1] - regression.points[i,1]

      Regression.Coefficients[i,] = cbind((rise/run),regression.points[i,1],regression.points[i+1,1])
      Regression.Coefficients[q,] = cbind(1,regression.points[i,1],regression.points[i,1]+1e-10)
  }

  Regression.Coefficients= na.omit(Regression.Coefficients)

 ### Fitted Values
  p = length((Regression.Coefficients)[,1])


  fitted = numeric()
  fitted.new = numeric()

  for (i in 1:p){

      z=(which(x>=Regression.Coefficients[i,2] & x<Regression.Coefficients[(i),3]))

      z.diff = ((x[z]- Regression.Coefficients[i,2])*Regression.Coefficients[i,1])+regression.points[i,2]


      if(is.null(point.est)){point.est.y = NULL} else{

          if(!is.null(point.est) && point.est>=Regression.Coefficients[i,2] && point.est<Regression.Coefficients[i,3]){ point.est.y = (point.est - Regression.Coefficients[i,2])*(Regression.Coefficients[i,1])+regression.points[i+ceiling(abs(p-q)/2),2]}

            else{if(!is.null(point.est) && point.est<Regression.Coefficients[1,2]){
     point.est.y = ((point.est - Regression.Coefficients[1,2])*(Regression.Coefficients[1,1]))+(regression.points[1,2])
      }

                else{if(!is.null(point.est) && point.est>Regression.Coefficients[p,2]){point.est.y = ((point.est - Regression.Coefficients[(p-0),2])*(Regression.Coefficients[(p-1),1]))+(regression.points[(p+ceiling(abs(p-q)/2)),2])
      }
     }
    }
  }


   fitted.new =  cbind(z,z.diff)


   fitted = rbind(fitted,fitted.new)
   fitted = fitted[order(fitted[,1]),]

 }

  if(print.equation==TRUE){
    print(regression.points)
    print(Regression.Coefficients)
  }

  Values = (cbind(x,Fitted=fitted[,2],Actual=y,Difference=fitted[,2]-(y)))

  MSE = mean((fitted[,2]-y)^2)

  R=cor(fitted[,2],y)
  R2=R^2

  R2.adj = 1 - (((1-R2)*length(fitted))/(length(fitted)-p-1))

  ###Plotting and regression equation
  xmin= min(c(point.est,x))
  xmax= max(c(point.est,x))
  ymin= min(c(point.est.y,y))
  ymax= max(c(point.est.y,y))
  plot(x,y,xlim=c(xmin,xmax),ylim=c(ymin,ymax),col='steelblue', xlab = "X",ylab="Y",main=paste0("Order = ",order+1))


  ### Plot Regression points and fitted values and legend
  points(na.omit(regression.points[order(regression.points),]),col='red',pch=19)
  lines(na.omit(regression.points[order(regression.points),]),col='red',lwd=2,lty = 2)


  if(!is.null(point.est)){ points(point.est,point.est.y, col='green',pch=18)
  legend(location, bty="n", y.intersp = 0.75,legend=c(paste("R2",format(R2,digits=4)),
    paste("R2 Adjusted",format(R2.adj,digits=4)),paste("Predictors",(p-1)),
    paste("Point Estimate",point.est),
    paste("Fitted Value",format(point.est.y,digits = 6))
  ))}

  if(is.null(point.est)){
    legend(location, bty="n", y.intersp = 0.75,legend=c(paste("R2",format(R2,digits=4)),paste("R2 Adjusted",format(R2.adj,digits=4)),paste("Predictors",(p-1))))
  }

  if(!is.null(point.est)){
  if(point.est>max(x)) segments(point.est,point.est.y,regression.points[p+ceiling(abs(p-q)/2),1],regression.points[p+ceiling(abs(p-q)/2),2],col="green",lty=2)
  if(point.est<min(x)) segments(point.est,point.est.y,regression.points[1,1],regression.points[1,2],col="green",lty=2)
  }


### Print Values

  if(print.values ==FALSE){
    if(is.null(point.est)){
    return(c("Predictors" = p-1,"R2"=R2,"R2 Adjusted"=R2.adj))
    }}

  if(print.values ==FALSE){
    if(!is.null(point.est)) {
    print(c("Predictors" = p-1,"R2"=R2,"R2 Adjusted"=R2.adj))
    return(c(Point=point.est, Fitted.value=point.est.y))
    }}

  if(print.values ==TRUE){
    if(is.null(point.est)){
    print(Values)
    return(c("Predictors" = p-1,"R2"=R2,"R2 Adjusted"=R2.adj))
    }}


  if(print.values ==TRUE){
    if(!is.null(point.est)) {
      print(Values)
      print(c("Predictors" = p-1,"R2"=R2,"R2 Adjusted"=R2.adj))
      return(c(Point=point.est, Fitted.value=point.est.y))
    }}



}
