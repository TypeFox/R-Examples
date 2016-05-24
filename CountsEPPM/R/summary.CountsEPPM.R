summary.CountsEPPM <-
function(object, ...) { cat('\n')
          cat('Model type:',object$model.type,'\n')
          cat('Model     :',object$model,'\n')
          if ((object$model=='Faddy distribution fixed b') | 
              (object$model=='negative binomial fixed b') | 
              (object$model=='general fixed b')) { 
             cat('b fixed at:',object$fixed.b,'\n') } # end if output$model
          cat('Link for mean         :','log','\n')   
# Introduction of link information for scale factor and variance in Version 1.01
          if (object$model.type=='mean and variance') { 
             if (object$scale.factor.model=='yes') { 
                          cat('Link for scale factor :','log','\n') 
                 } else { cat('Link for variance     :','log','\n') }}   
          offsetid_mean     <- sum(object$offset.mean)
          offsetid_variance <- sum(object$offset.variance)
          if (offsetid_mean!=0) { 
             cat('non zero offsets in linear predictor for mean','\n') }   
          if (offsetid_variance!=0) { 
              cat('non zero offsets in linear predictor for variance','\n') }   
          if (object$model.type=='mean and variance') { 
             if (object$scale.factor.model=='no') { 
                 cat('variance model fitted','\n') 
                    } else {
                 cat('scale factor model model fitted','\n') }}   
          if (is.na(object$estses[1,3])==TRUE) { 
              cat('Determinant of hessian matrix is zero','\n') 
              cat('or the hessian matrix is ill conditioned.','\n') }

          if (is.na(object$loglikelihood)==FALSE) {
             if ((object$model=="Faddy distribution") | 
                 (object$model=="Faddy distribution fixed b")) { 
                  npar   <- nrow(object$estses)
                  nparm1 <- npar - 1
               if (object$model=="Faddy distribution") {
                  wk.c <- round(object$estses[nparm1,2],digits=7)
                  } else {
                  wk.c <- round(object$estses[npar,2],digits=7) }
               if (wk.c==1) {
                 cat('Boundary for c of 1 has been reached','\n') 
                 cat('hence its se is set to NA.','\n') }}} # end is.na 
               
# Introduced in Version 1.01 to deal with the missing name se for 
# the single parameter situation.
          names(object$estses) <- c('name','Estimates','se')

# standard normal z values, p values and asterisks 
          zvalue <- object$estses$'Estimates' / object$estses$'se'
          pvalue <- 2*(1-pnorm(abs(zvalue)))
          asterisks <- pvalue
          asterisks <- sapply(1:length(pvalue), function(i) 
                    if (is.finite(asterisks[i])==FALSE) { asterisks[i] <- '   '
                    } else { if (asterisks[i]<0.001) { asterisks[i] <- '***'
                           } else { if (asterisks[i]<0.01) { asterisks[i] <- '** '
                                  } else { if (asterisks[i]<0.05) { asterisks[i] <- '*  '
                                     } else { if (asterisks[i]<0.10) { asterisks[i] <- '.  '
                                        } else { asterisks[i] <- '   ' }}}}} )

          coeff.table <- data.frame(object$estses, zvalue, pvalue, asterisks)
          cat('Parameter estimates and se\'s','\n')
          names(coeff.table) <- c('name','Estimate','Std. Error','z value','Pr(>|z|)',' ')
          wks <- length(object$estses[1])
          if (wks) {
              print.data.frame(coeff.table, row.names=FALSE)
              cat("\n","Signif. codes: '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ","\n")
                   } else { cat("No coefficients in model)\n\n") }

          if (is.finite(object$loglikelihood)==FALSE) {
                  cat('loglikelihood is NA or -Inf suggesting an inappropriate model','\n')
                                                      } else { 
                  cat('\n','Log-likelihood:',object$loglikelihood,'on', nrow(object$estses),
                       'Df', '\n', sep=" ")
                 } # end of is.finite(object$loglikelihood
                                   }
