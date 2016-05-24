modelAssumptions <-
function(Model, Type='NORMAL', ID=row.names(Model$model), one.page=TRUE)
#Provides diagnositic graphs and score tests to evaluate assumptions of normality, constant variance and linearity
#model: an model object from lm
#type: NORMAL, CONSTANT, or LINEAR
#ID:  Use to identify points.  Default = row.names(thedata).  NULL = no identification
#one.page:  put all graphs on one page

#Revision history
#2011-02-28:  added print to display ncvTest(), JJC


{
  switch(toupper(Type),
  
     NORMAL =
     {
     
        if (one.page) {
            dev.new(width=14,height=7)
            par(mfrow = c(1,2))} 
        else  {
            dev.new(width=7,height=7,record=TRUE)
        }
        par(cex.lab=1.5, cex.axis=1.2, lwd=2)   #make pretty graphs
          
        qqPlot(Model, labels=FALSE, sim=TRUE, main='Quantile-Comparison Plot to Assess Normality', xlab='t Quantiles', ylab ='Studentized Residuals')  ##Quantile-comparison plot of studentize residuals vs. t-distribution
        plot(density(rstudent(Model)),main = 'Density Plot to Assess Normality of Residuals', xlab='Studentized Residual')       #Non-parametric density plot of studentized residuals
        zx <- seq(-4, 4, length.out=100)
        lines(zx,dnorm(zx, mean=0, sd=sd(rstudent(Model))),lty = 2, col="blue")
        cat('Descriptive Statistics for Studentized Residuals\n')
        describe(rstudent(Model))  
                                   
     },
  
     CONSTANT =
     {
        if (one.page) {
          dev.new(width=14,height=7)
          par(mfrow = c(1,2))} 
        else  {
          dev.new(width=7,height=7,record=TRUE)
        }
        par(cex.lab=1.5, cex.axis=1.2, lwd=2)   #make pretty graphs
        
        plot(rstudent(Model) ~ fitted.values(Model), main='Studentized Residuals vs. Fitted Values', xlab='Fitted Values', ylab='Studentized Residuals')   #Studentized residuals vs. Fitted values
        abline(h=0,lty=2, col="blue")
        print(spreadLevelPlot(Model)) #log(abs(studentized residuals) vs. log(fitted values).
        cat('\n\n')
        print(ncvTest(Model))
     },
     
     LINEAR =
     {
        dev.new(width=7,height=7,record=TRUE)
        par(cex.lab=1.5, cex.axis=1.2, lwd=2)   #make pretty graphs
        crPlots(Model, ask=TRUE)     
     },
     {print('Valid options for type: normal, constant, linear')}     #OTHERWISE
  )#end switch
  
  
  #Provide Global test of linear Model assumptions
  print(gvlma(Model))
}

