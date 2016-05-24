simexaft <-
function(formula=formula(data),data=parent.frame(),SIMEXvariable,repeated =FALSE,repind=list(),err.mat=err.mat,B=50,lambda=seq(0,2,0.1),extrapolation="quadratic",dist="weibull")
{


       ############################ check the input of SIMEXvariable #######################################
       colname=colnames(data)
       SIMEXvariable=unique(SIMEXvariable)
       nSIMEXvariable=length(SIMEXvariable)

       if(!is.character(SIMEXvariable) | nSIMEXvariable>length(colname)){

               stop("Invalid SIMEXvariable object")

               }
       if(!all(SIMEXvariable %in% colname)){

               stop("SIMEXvariable must selected from the data")

               }


       #################### check the input of repeated, err.mat, repind ##################################
       if (!(repeated == FALSE) & !(repeated == TRUE) ) {

               stop("Repeated indicator should only be 'TRUE' or 'FALSE'. ")

               }

       if(repeated==FALSE){

                      err.mat=as.matrix(err.mat)

                      if(!is.numeric(err.mat) | any(err.mat < 0)){

                             stop("Invalide err.mat object, err.mat must be a square symmetric numeric matrix")

                             }

                      if (nrow(err.mat) != ncol(err.mat)) {

                             stop("err.mat must be a square matrix")

                             }

                      if (length(SIMEXvariable) != nrow(err.mat)) {

                            stop("SIMEXvariable and err.mat have non-conforming size")

                            }


                      SSigma <- err.mat

                      dimnames(SSigma) <- NULL

                      if (!isTRUE(all.equal(SSigma, t(SSigma)))) {

                             warning("err.mat is numerically not symmetric")

                             }
                }


       else if(repeated==TRUE){

                    if(length(SIMEXvariable) != length(repind)){

                             stop("SIMEXvariable and repind have non-conforming size")

                             }

                   }




       ##################### specifies the parameter B #########################
       if(length(B)!=1){

                stop("B must be positive integer")

                }

       if(!is.numeric(B) | B<=0 ){

                stop("B must be positive integer")

                }

        else{

                B=ceiling(B)

                }




         ############################  lambda ################################
         if(!is.vector(lambda) |!is.numeric(lambda)){

                   stop(":Invalide lambda object")

                   }

        if (any(lambda < 0)) {

               warning("Lambda should be positive values. Negative values will be ignored",call. = FALSE)

               lambda <- lambda[lambda >= 0]

               }




        ##### specifies the regression form used in the extrpolation step #####
        extrapolation = substr(extrapolation, 1, 4)

        if(!is.character(extrapolation) | length(extrapolation)!=1){

                     warning("Invalid extrapolation object. Using: quadatic\n\n",call.=FALSE)
                     
                     extrapolation="quad"
                     }



        ##############################  prepare reading data from data.frame############


        ndata=nrow(data)


        nformula=length(attr(terms(formula),"term.labels"))+1
        nlambda=length(lambda)

        #the matrixes to save the estimates of the coefficents, variance #
        A1=matrix(data=NA,nlambda,nformula)

        A2=matrix(data=NA,nlambda,nformula)

        A3=matrix(data=NA,nlambda,nformula)

        #the matrix of estimate#
        theta=matrix(data=NA,B,nformula)

        colnames(theta)=c("Intercept",attr(terms(formula),"term.labels"))
        p.names=colnames(theta)

        #all estimates corresonding to each lambda#
        theta.all=vector(mode="list",nlambda)


        #####################use SIMEX method to analysis the data#########
        k=1
        while(k <= length(lambda))
        {
                 ## the coefficients estimates of the kth sample##
                 w=numeric()

                 ##the variance of the coefficent estimate of the kth sample##
                 v=numeric()

                 ##the variance estimate of the kth sample##
                 omega=numeric()


                 temp=data

                 #the matrix of variance estimate#
                 estivarB=matrix(data=NA,B,nformula)


                 #the vector of the estimate scale#
                 estiscaleB=matrix(data=NA,B,ncol=1)

                 ### simulation step of the SIMEX algorithm ###
                 for(r in 1:B)
                 {

                      if(repeated==FALSE){

                               temp[SIMEXvariable]=data[SIMEXvariable]+sqrt(lambda[k])*rmvnorm(ndata,rep(0,length(SIMEXvariable)),err.mat)

                           }


                       else{
                                ## use repeated data set to estimate the varaince of the measurement error ###
                                ## generate the contrasts and get the corresponding W ##
                                constrast=list()

                                for(i in 1:nSIMEXvariable){

                                    n.i=length(repind[[i]])

                                    z.i=rnorm(n.i, 0, 1)

                                    constrast[[i]]=(z.i-mean(z.i))/sqrt(sum((z.i-mean(z.i))^2))

                                    mean.i=apply(temp[repind[[i]]],1,mean)

                                    temp[SIMEXvariable[i]]= mean.i + sqrt(lambda[k]/n.i)*as.matrix(temp[repind[[i]]])%*%as.vector(constrast[[i]])


                                    }


                           }

                           ### use survreg to fit the artifical simulated data to the the AFT model ###
                           re = survreg(formula=formula,data=temp,dist=dist)

                           ##obtain the coefficients estimates w, associated variance estimate omega and the scale estimate scale##
                           scale=re$scale
                           w=re$coefficients
                           omega=diag(re$var)[1:nformula]
                           theta[r,]=w
                           estivarB[r,]=omega
                           estiscaleB[r,]=scale
                 }

                 ###for each lambda value calculate the mean of the each estimates of B samples ###
                 ##the coefficients estimates and the variance of the estiamtes##
                 w=apply(theta,2,FUN=mean)
                 v=apply(theta,2,FUN=var)


                 ##the variance estimates##
                 omega =apply(estivarB,2,FUN=mean)


                 ##the difference of the variance estimates and the varaince of the coefficient estimates, use this to estimate the variance##
                 tau=omega-v


                 ### remove the sample with negative tau###
                 if(all(tau>0)){
                       A1[k,] = w
                       A2[k,] = tau
                       A3[k,] = apply(estiscaleB,2,FUN=mean)
                       theta.all[[k]]=theta
                       k = k + 1
                       }


        }

       ### save all the estimates in the theta.all###
       theta=matrix(unlist(theta.all),nrow=B)
       theta.all=list()

       for (i in 1:nformula){

        theta.all[[p.names[i]]] <- data.frame(theta[,seq(i, nformula * nlambda, by = nformula)])
       }



      #####################extrapolation step of the SIMEX algorithm, fit the regression model#####################
       if(extrapolation=="line"){

                    result1=linearextrapolation(A1,A2,A3,lambda)


                     }


        else if(extrapolation=="quad"){

                    result1=quadraticextrapolation(A1,A2,A3,lambda)


                     }


        else stop("extrapolation method must be linear or quadratic")




        ####### outputs###########################
        estimate=result1$reg1
        names(estimate)=p.names
        se=sqrt(result1$reg2)
        names(se)=p.names
        scalereg=result1$scalereg
        pvalue=2*(1-pnorm(abs(estimate/se)))



        if(repeated==FALSE){
        erg=list(coefficients=estimate,se=se,scalereg=scalereg,pvalue=pvalue,lambda=lambda, B=B, formula=formula, err.mat=err.mat, extrapolation=extrapolation,SIMEXvariable=SIMEXvariable,theta=theta.all)

        }


        else{
        erg=list(coefficients=estimate,se=se,scalereg=scalereg,pvalue=pvalue,lambda=lambda, B=B, formula=formula, extrapolation=extrapolation,SIMEXvariable=SIMEXvariable,repind=repind,theta=theta.all)

        }


        class(erg)<-("simexaft")


        return(erg)

}
