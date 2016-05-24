

dynamicCorrelation <- function(dataFrame,
                               depVar = c('resp1', 'resp2', 'resp3', 'resp4', 'resp5'),
                               indepVar = 'time',
                               subjectVar = 'subject',
                               function.choice = c(1),
					 width.range = c(((range(dataFrame[[indepVar]]))[2]-(range(dataFrame[[indepVar]]))[1])/4, 
							((range(dataFrame[[indepVar]]))[2]-(range(dataFrame[[indepVar]]))[1])/4),
					 width.place = c(NA, NA), 
					 boundary.trunc = c(0, 0),
                               lag.input = c(),
                               byOrder = c(),
                               by.deriv.only = TRUE,
                               full.lag.output = FALSE)

	 ## arguments are:
       ##               dataFrame       ( The dataFrame that contains the dataset)
	 ##			dep.var         ( dependents or y-axis variables)
	 ##			indep.var       ( independent or x-axis variable, usually the recording time for the data)
	 ##			function.choice ( a vector of any length to indicate which derivatives
	 ##				            are desired for each curve smoothing; 1st is 0th derivative, 2nd is 1st
	 ##				            derivative, 3rd is 2nd derivative, etc; 1=yes, 0=no
       ##                                 e.g. c(1,0,1) means we would like to look at derivative 0 and derivative 2 )		
	 ##			width.range     ( bandwidth (size) for local regression in curve smoothing, it can be a list that specifies different bandwidth 
       ##                                 vector for each response, or a single width vector that is used among all the responses )	 
	 ##			width.place     ( endpoints for width change, it can be a list that specifies difference endpoints for each response, or a single
       ##                                 vector of endpoints that is used among all the responses )
	 ##		      boundary.trunc  ( indicate the boundary of indep.var that should be truncated after smoothing)
       ##               lag.input       ( values of lag to be considered, can be a vector of requested lags)
       ##               byOrder         ( a vector that specifies the order of the variables and derivatives to be the leading variable in the calculation 
       ##                                 of the dynamic correlations)
       ##               by.deriv.only   ( If true, the inter-dynamic correlations between different derivatives are not computed, which 
       ##                                 can same computation time)
       ##               full.lag.output ( If true, the dynamic correlation values for each response and derivative will be stored in vectors corresponding
       ##                                 to different lag values, which enables plotting the correlations vs lag values. All the vectors are stored in the 
       ##                                 returned attribute resultMatrix )
{  
      ## set up the variables
      attach(dataFrame)
      vec <-  unique(dataFrame[[subjectVar]])
      dep.var= dataFrame[depVar]
      subject.var = dataFrame[[subjectVar]]
      indep.var=dataFrame[[indepVar]]
       
      #maybe used to find the number of days to be used later
      min<-c()
      for(i in 1:length(vec))
      { min <- c(min, max(indep.var[subject.var==vec[i]]))
      }
      limit<- min(min)
      
      ## To check which derivatives to be considered
      funcVar<-c()
      for(i in 1:length(function.choice))
      { if(function.choice[i]==1)
        { funcVar<-c(funcVar, i)
        }
      }

      ## check parameters
      if(length(byOrder)!=0)
      {  largest<-length(depVar)*length(funcVar)
         if(max(byOrder)>largest)
         { stop("The specified index order of leading variables is over range")
         }
         if(length(byOrder)!=largest && length(byOrder)!=0)
         { stop("The index order of leading variables is not completely specified")
         }
      }
    
      ## Make the width.range and the width.place parameters to be lists of the same length as the number of responses
      if(!is.list(width.range) && !is.list(width.place))
      {  widthRange <- width.range
         widthPlace <- width.place  
         width.range <- width.place <- vector(mode='list', length(depVar))
         for(i in 1:length(depVar))
         { width.range[[i]] <- widthRange
           width.place[[i]] <- widthPlace
         }
      }else if(!is.list(width.range) && is.list(width.place))
      {  widthRange <- width.range
         width.range <- vector(mode='list', length(depVar))
         for(i in 1:length(depVar))
         {  width.range[[i]] <- widthRange  
         } 
         if(length(width.place)!= length(depVar))
         {  stop("If width.place is a list, it needs to have the same number of components as the number of responses")
         }
      }else if(is.list(width.range) && !is.list(width.place))
      {  widthPlace <- width.place
         width.place <- vector(mode='list', length(depVar))
         for( i in 1:length(depVar))
         {  width.place[[i]] <- widthPlace
         }
         if(length(width.range) != length(depVar))
         {  stop("If width.range is a list, it needs to have the same number of components as the number of responses")
         }
      }else
      {  if(length(width.place)!= length(depVar))
         {  stop("If width.place is a list, it needs to have the same number of components as the number of responses")
         }
         if(length(width.range) != length(depVar))
         {  stop("If width.range is a list, it needs to have the same number of components as the number of responses")
         }
      }

      smoothedCurves<-vector(mode='list', length(vec)) 
      ## curve smoothing 
      for (k in 1:length(vec))
      {  #cat('subject=', vec[k], ' ')
         indep<- indep.var[subject.var==vec[k]]
         smoothedCurves[[k]]<-vector(mode='list', length(depVar))
         for( i in 1:length(depVar))
         {  #cat('Dependent variable = ', depVar[i], '\n') 
            if( is.na(width.place[[i]][1]))
            {  width.place[[i]][1]<-min(indep[!is.na(indep)])
            }
            if( is.na(width.place[[i]][2]))
            {  width.place[[i]][2]<- max(indep[!is.na(indep)])
            }
            ## compute the bandwidth for each point to be outputted for all subjects
            points.use = seq(boundary.trunc[1],ceiling(max(indep[!is.na(indep)])-boundary.trunc[2]), by=1)
            size <- ceiling(max(indep[!is.na(indep)])-boundary.trunc[2]) - boundary.trunc[1]+1
            band <-c()
            for( count in 1:size)
            {  if(abs(width.range[[i]][1] - width.range[[i]][2]) < 1e-5)  
													### for single-bandwidth smoothing
	             width <- width.range[[i]][1]
	         else										### for >1 bandwidth smoothing
	         {
	             widthrange <- width.range[[i]][2] - width.range[[i]][1]
	             if(points.use[count] < width.place[[i]][1])
                       width <- width.range[[i]][1]
	             else if(points.use[count] >= width.place[[i]][1] & points.use[count] <= width.place[[i]][2])
                       width <- width.range[[i]][1] + 
		            ((points.use[count]-width.place[[i]][1])/
		            (width.place[[i]][2]-width.place[[i]][1])) * widthrange
	             else width <- width.range[[i]][2]
	         }
               band<-c(band, width)
            }
            smoothedCurves[[k]][[i]]<-vector(mode='list', length(funcVar))
            for( j in 1: length(funcVar))
            { #cat('deriative function index= ', funcVar[j], '\n')
              temp<- lpepa(x=indep.var[subject.var==vec[k]][!is.na(dep.var[[i]][subject.var==vec[k]])], y=dep.var[[i]][subject.var==vec[k]][!is.na(dep.var[[i]][subject.var==vec[k]])], 
              bandwidth=band,  deriv =(funcVar[j]-1), n.out=size,  x.out= seq(boundary.trunc[1], ceiling(max(indep[!is.na(indep)])-boundary.trunc[2])),
              order=funcVar[j], var= FALSE)$est
              smoothedCurves[[k]][[i]][[j]]<- temp
            }
         }
      } 
       
       ###############################
       ####    Obtain column means for responses 
 
       max.len<- (ceiling(limit)-boundary.trunc[2])-boundary.trunc[1]+1
       meanMatrix<-vector(mode='list', length(depVar))
       for( i in 1:length(depVar))
       {  meanMatrix[[i]]<-vector(mode='list', length(funcVar))
          #cat('Dependent Variable = ', depVar[i], '\n')
          for( j in 1:length(funcVar))
          {  #cat('deriative function index= ', funcVar[j], '\n')
             forEachDerivMean <- rep(NA, max.len)
             for( m in 1:max.len)
             { for (n in 1:length(vec))
               {  if(n==1)
                  { total <- 1
                    forEachDerivMean[m] <- smoothedCurves[[n]][[i]][[j]][m]
                  }else if( n > 1 )
                  { total <- total+1
                    forEachDerivMean[m] <- (1/total)*((total-1)*forEachDerivMean[m] + smoothedCurves[[n]][[i]][[j]][m])      
                  }
               }
             }            
             meanMatrix[[i]][[j]] <- forEachDerivMean
          }
       }


       ####### now, standardize data within each curve, and generate
       #######  cross-correlation between estimated curves and derviatives

       ####### standardize the response estimates and derivatives and 
       #######  place all value (including grid values) into one large list 
        
       correst.all.standz <- vector(mode='list', length(vec))
       for( k in 1:length(vec))
       {  #cat('subject= ', k)
          correst.all.standz[[k]]<-list()
          correst.all.standz[[k]][[1]] <- list(boundary.trunc[1]:(boundary.trunc[1]+max.len-1))
          ## grid values above
          for( i in 2:(length(depVar)+1))
          {  correst.all.standz[[k]][[i]]<- vector(mode='list', length(funcVar))
             for( j in 1:length(funcVar))
             {  temp <- smoothedCurves[[k]][[i-1]][[j]][1:max.len]
                correct <- temp - meanMatrix[[i-1]][[j]]
                mean.correct <- mean(correct)
                sd.correct <- sqrt(var(correct))
                correst.all.standz[[k]][[i]][[j]]<- list((correct-mean.correct)/(sqrt((max.len-1)/max.len)*sd.correct))
             }
          }
       }
       
       
       ## calculate weighted correlations from correst.all.standz

       weights.vec <- rep(NA, length(vec))
       time.extend <- (boundary.trunc[1] + max.len - 1) + boundary.trunc[1] 
       
       if(by.deriv.only==TRUE)##Then we do not include inter-dynamic-correlation compuatations
       {  cov.mtx.listz <- vector(mode='list', length(funcVar))
          cov.wgt.mtxz <- vector(mode='list', length(funcVar))
          dim <- length(depVar)
          for( deriv in 1:length(funcVar))
          {  cov.mtx.listz[[deriv]] <- vector(mode='list', length(vec))
             for( i in 1:length(vec))
             {  # use days with actual observations for weights
                weights.vec[i] <- length(subject.var[(subject.var==vec[i]) & 
							(indep.var >= 0) & (indep.var <= time.extend)])
                 cov.mtx.listz[[deriv]][[i]] <- matrix(nrow=dim, ncol=dim)
                 m <- max.len
                 diag(cov.mtx.listz[[deriv]][[i]]) <- 1
                 for(j in 1: dim)
                 {  for (k in (j+1):dim)
                    { if(k <= dim)
                      {   cov.mtx.listz[[deriv]][[i]][j,k] <- cov.mtx.listz[[deriv]][[i]][k,j] <- 
                                 ((1/m)*sum(correst.all.standz[[i]][[(j+1)]][[deriv]][[1]]
                            *correst.all.standz[[i]][[(k+1)]][[deriv]][[1]]))
                      }
                    }
                 }
             }
             cov.wgt.mtxz[[deriv]] <-matrix(0,dim,dim)
             for( i in 1:length(vec))
             { cov.wgt.mtxz[[deriv]] <-cov.wgt.mtxz[[deriv]] + (weights.vec[i] * cov.mtx.listz[[deriv]][[i]])
             }
             cov.wgt.mtxz[[deriv]] <- (1/sum(weights.vec)) * cov.wgt.mtxz[[deriv]]
             names<-c()
             for(i in 1:length(dep.var))
             {  name<-depVar[i]
                name2<-paste(name, (deriv-1))
                names<-c(names, name2)   
             }
             dimnames(cov.wgt.mtxz[[deriv]]) <- list(names,names)
          }
       }else
       {  cov.mtx.listz <- list()
          dim <- length(dep.var)*length(funcVar)
          for( i in 1:length(vec))
          {  # use days with actual observations for weights
             weights.vec[i] <- length(subject.var[(subject.var==vec[i]) & (indep.var >= 0) & (indep.var <= time.extend)])
             cov.mtx.listz[[i]] <- matrix(nrow=dim, ncol=dim)
             m <- max.len
             diag(cov.mtx.listz[[i]]) <- 1
             for(j in 1: dim)
             {  for (k in (j+1):dim)
                { if( k <= dim)
                  {   m1<- (j-1)%/%length(funcVar)+2
                      m2<-(j-1)%%length(funcVar)+1
                      s1<-(k-1)%/%length(funcVar)+2
                      s2<-(k-1)%%length(funcVar)+1
                      cov.mtx.listz[[i]][j,k] <- cov.mtx.listz[[i]][k,j] <- 
                             ((1/m)*sum(correst.all.standz[[i]][[m1]][[m2]][[1]]
                            *correst.all.standz[[i]][[s1]][[s2]][[1]]))
                  }
                }
             }
          }
          
          cov.wgt.mtxz <-matrix(0,dim,dim)
          for( i in 1:length(vec))
          { cov.wgt.mtxz <-cov.wgt.mtxz + (weights.vec[i] * cov.mtx.listz[[i]])
          }
          cov.wgt.mtxz <- (1/sum(weights.vec)) * cov.wgt.mtxz
          names<-c()
          for( i in 1:length(depVar))
          {   name<-depVar[i]
              for( j in 1:length(funcVar))
              {  name2<-paste(name, (funcVar[j]-1))
                 names<-c(names, name2)
              }
          }
          dimnames(cov.wgt.mtxz) <- list(names,names)
       }
       cov.wgt.mtxz

     ##############################################################################
     # now incorporate testing for lags which maximize correlations, 
     #   using new lag shift, with common weight function 
     ##############################################################################
     if(length(lag.input)!=0)
     {   if(by.deriv.only == FALSE)
         {dim<-length(depVar)*length(funcVar)
          cov.lag.mtx.listz <-list()
          for(i in 1:length(vec))
          {  #cat('subject = ', vec[i])
             cov.lag.mtx.listz[[i]] <-list()
           
             ##corrected curves for subject i
             dep.correct<-list()
             for(m in 1:length(depVar))
             {  dep.correct.function<-list()
                for(n in 1:length(funcVar))
                { dep.correct.function[[n]]<- smoothedCurves[[i]][[m]][[n]][1:max.len]-meanMatrix[[m]][[n]]
                }
                dep.correct[[m]]<-dep.correct.function
             }
           
             for( j in 1:length(lag.input))
             {  #cat(' grid pos. =', lag.input[j])
                cov.lag.mtx.listz[[i]][[j]]<- matrix(nrow=dim, ncol=dim)
                diag(cov.lag.mtx.listz[[i]][[j]]) <- 1
                if( lag.input[j] >= 0)
                { lag.end <- max.len - lag.input[j]
                  lag.beg <- 1 + lag.input[j]
                  m.support <- lag.end
                  ##standardize each curve according to the lag value and the order of the leading variables, then put into the matrix
                  if(length(byOrder)==0)
                  {  correst.standz<-list()
                     index<-1
                     mean<- mean(dep.correct[[1]][[1]][1:lag.end])
                     sd <- sqrt(var(dep.correct[[1]][[1]][1:lag.end]))
                     correst.standz[[index]] <- (dep.correct[[1]][[1]][1:lag.end]-mean)/
                                                (sqrt((m.support-1)/m.support)*sd)
                     for(index in 2:dim )
                     {  if(dim >=2)
                        {  dep<-(index-1)%/%length(funcVar)+1
                           deriv<-(index-1)%%length(funcVar)+1
                           mean<- mean(dep.correct[[dep]][[deriv]][lag.beg:max.len])
                           sd<- sqrt(var(dep.correct[[dep]][[deriv]][lag.beg:max.len]))
                           correst.standz[[index]] <- (dep.correct[[dep]][[deriv]][lag.beg:max.len]-mean)/
                                                  (sqrt((m.support-1)/m.support)*sd)
                           cov.lag.mtx.listz[[i]][[j]][1,index]<-cov.lag.mtx.listz[[i]][[j]][index,1] <-
                                        (1/m.support)*sum(correst.standz[[1]]*correst.standz[[index]])
                        }
                     }
                     for( d in 2:(dim-1)  )
                     {  if( dim >=3)
                        {  dep<- (d-1)%/%length(funcVar)+1
                           deriv<-(d-1)%%length(funcVar)+1
                           mean<- mean(dep.correct[[dep]][[deriv]][1:lag.end])
                           sd <- sqrt(var(dep.correct[[dep]][[deriv]][1:lag.end]))
                           correst.standz[[d]] <- (dep.correct[[dep]][[deriv]][1:lag.end]-mean)/
                                                (sqrt((m.support-1)/m.support)*sd)
                           for( e in (d+1):dim)
                          {  cov.lag.mtx.listz[[i]][[j]][d,e]<-cov.lag.mtx.listz[[i]][[j]][e,d] <-
                                       (1/m.support)*sum(correst.standz[[d]]*correst.standz[[e]])
                          }
                        }
                     }
                  ##end length of byorder=0
                  }else
                  {  correst.standz<-list()
                     index<-byOrder[1]
                     dep<-(index-1)%/%length(funcVar)+1   ## %/% get the integral party if divided
                     deriv<-(index-1)%%length(funcVar)+1  ## %% get the module
                     mean<- mean(dep.correct[[dep]][[deriv]][1:lag.end])
                     sd <- sqrt(var(dep.correct[[dep]][[deriv]][1:lag.end]))
                     correst.standz[[index]] <- (dep.correct[[dep]][[deriv]][1:lag.end]-mean)/
                                                (sqrt((m.support-1)/m.support)*sd)
                     for( ord in 2:dim )
                     {  if(dim >=2)
                        {  index<-byOrder[ord]
                           dep<-(index-1)%/%length(funcVar)+1
                           deriv<-(index-1)%%length(funcVar)+1
                           mean<- mean(dep.correct[[dep]][[deriv]][lag.beg:max.len])
                           sd<- sqrt(var(dep.correct[[dep]][[deriv]][lag.beg:max.len]))
                           correst.standz[[index]] <- (dep.correct[[dep]][[deriv]][lag.beg:max.len]-mean)/
                                                  (sqrt((m.support-1)/m.support)*sd)
                           cov.lag.mtx.listz[[i]][[j]][byOrder[1],index]<-cov.lag.mtx.listz[[i]][[j]][index,byOrder[1]] <-
                                            (1/m.support)*sum(correst.standz[[byOrder[1]]]*correst.standz[[index]])
                        }
                     }
                     for( ord in 2:(dim-1)  )
                     {  if( dim >=3)
                        {  d<-byOrder[ord]
                           dep<- (d-1)%/%length(funcVar)+1
                           deriv<-(d-1)%%length(funcVar)+1
                           mean<- mean(dep.correct[[dep]][[deriv]][1:lag.end])
                           sd <- sqrt(var(dep.correct[[dep]][[deriv]][1:lag.end]))
                           correst.standz[[d]] <- (dep.correct[[dep]][[deriv]][1:lag.end]-mean)/
                                                (sqrt((m.support-1)/m.support)*sd)
                           for( ord2 in (ord+1):dim)
                           {  e<-byOrder[ord2]
                              cov.lag.mtx.listz[[i]][[j]][d,e]<-cov.lag.mtx.listz[[i]][[j]][e,d] <-
                                          (1/m.support)*sum(correst.standz[[d]]*correst.standz[[e]])
                           }
                        }
                     }
                  }
                 
               }else if(lag.input[j]<0)
               {  lag.end <- max.len+lag.input[j]
                  lag.beg <- 1 - lag.input[j]
                  m.support<- lag.end
                  ##standardize each curve according to the lag value and the order of leading variables specified first, then put into the matrix
                  if(length(byOrder)==0)
                  {  correst.standz<-list()
                     index<-1
                     mean<- mean(dep.correct[[1]][[1]][lag.beg:max.len])
                     sd <- sqrt(var(dep.correct[[1]][[1]][lag.beg:max.len]))
                     correst.standz[[index]] <- (dep.correct[[1]][[1]][lag.beg:max.len]-mean)/
                                                (sqrt((m.support-1)/m.support)*sd)
                     for(index in 2:dim)
                     {  if(dim>=2)
                        {  dep<-(index-1)%/%length(funcVar)+1
                           deriv<-(index-1)%%length(funcVar)+1
                           mean<- mean(dep.correct[[dep]][[deriv]][1:lag.end])
                           sd<- sqrt(var(dep.correct[[dep]][[deriv]][1:lag.end]))
                           correst.standz[[index]] <- (dep.correct[[dep]][[deriv]][1:lag.end]-mean)/
                                                  (sqrt((m.support-1)/m.support)*sd)
                           cov.lag.mtx.listz[[i]][[j]][1,index]<-cov.lag.mtx.listz[[i]][[j]][index,1] <-
                           (1/m.support)*sum(correst.standz[[1]]*correst.standz[[index]])

                        }
                     }
                     for( d in 2:(dim-1))
                     {  if( dim >=3)
                        {  dep<- (d-1)%/%length(funcVar)+1
                           deriv<-(d-1)%%length(funcVar)+1
                           mean<- mean(dep.correct[[dep]][[deriv]][lag.beg:max.len])
                           sd <- sqrt(var(dep.correct[[dep]][[deriv]][lag.beg:max.len]))
                           correst.standz[[d]] <- (dep.correct[[dep]][[deriv]][lag.beg:max.len]-mean)/
                                                (sqrt((m.support-1)/m.support)*sd)
                     
                           for( e in (d+1):dim)
                           {  cov.lag.mtx.listz[[i]][[j]][d,e]<-cov.lag.mtx.listz[[i]][[j]][e,d] <-
                                            (1/m.support)*sum(correst.standz[[d]]*correst.standz[[e]])
                           }
                      
                        } 
                     }
                  }else ##byorder!=0
                  {  correst.standz<-list()
                     index<-byOrder[1]
                     dep<-(index-1)%/%length(funcVar)+1
                     deriv<-(index-1)%%length(funcVar)+1
                     mean<- mean(dep.correct[[dep]][[deriv]][lag.beg:max.len])
                     sd <- sqrt(var(dep.correct[[dep]][[deriv]][lag.beg:max.len]))
                     correst.standz[[index]] <- (dep.correct[[dep]][[deriv]][lag.beg:max.len]-mean)/
                                                (sqrt((m.support-1)/m.support)*sd)
                     for(ord in 2:dim)
                     {  index<-byOrder[ord]
                        if(dim>=2)
                        {  dep<-(index-1)%/%length(funcVar)+1
                           deriv<-(index-1)%%length(funcVar)+1
                           mean<- mean(dep.correct[[dep]][[deriv]][1:lag.end])
                           sd<- sqrt(var(dep.correct[[dep]][[deriv]][1:lag.end]))
                           correst.standz[[index]] <- (dep.correct[[dep]][[deriv]][1:lag.end]-mean)/
                                                  (sqrt((m.support-1)/m.support)*sd)
                           cov.lag.mtx.listz[[i]][[j]][1,index]<-cov.lag.mtx.listz[[i]][[j]][index,1] <-
                                        (1/m.support)*sum(correst.standz[[1]]*correst.standz[[index]])
                        }
                     }
                     for( ord in 2:(dim-1))
                     {  if( dim >=3)
                        {  d<-byOrder[ord]
                           dep<- (d-1)%/%length(funcVar)+1
                           deriv<-(d-1)%%length(funcVar)+1
                           mean<- mean(dep.correct[[dep]][[deriv]][lag.beg:max.len])
                           sd <- sqrt(var(dep.correct[[dep]][[deriv]][lag.beg:max.len]))
                           correst.standz[[d]] <- (dep.correct[[dep]][[deriv]][lag.beg:max.len]-mean)/
                                                (sqrt((m.support-1)/m.support)*sd)
                      
                           for( ord2 in (ord+1):dim)
                           {  e<-byOrder[ord2]
                              cov.lag.mtx.listz[[i]][[j]][d,e]<-cov.lag.mtx.listz[[i]][[j]][e,d] <-
                                           (1/m.support)*sum(correst.standz[[d]]*correst.standz[[e]])
                           }
                      
                        }
                     }  

                  }##byorder!=0

                }#end of negative if
             }#end of j loop
          }#end of i loop
          
           
          ## calculate overall weighted covariance matrix for each lag;
          ## then determine which lag produces highest magnitude correlation for each unique matrix entry
          ## -> even with lags, we will leave some wieghts as for non-lagged approach
        
          cov.lag.wgt.mtx.listz <- list()
          for( j in 1:length(lag.input))
          {  cov.lag.wgt.mtx.listz[[j]]<-matrix(0, nrow=dim, ncol=dim)
             for( i in 1:length(vec))
             {  cov.lag.wgt.mtx.listz[[j]] <- cov.lag.wgt.mtx.listz[[j]]+(weights.vec[i]*cov.lag.mtx.listz[[i]][[j]])
             }
             cov.lag.wgt.mtx.listz[[j]]<-(1/sum(weights.vec))*cov.lag.wgt.mtx.listz[[j]]
          }

          ## find highest correlation amongst all lags for each individual 
          cov.lag.wgt.mtxz <-matrix(0, nrow=dim, ncol=dim)
          lag.max.cov.mtxz<-matrix(NA, nrow=dim, ncol=dim)
          for( j in 1:length(lag.input))
          {  for( k in 1:(dim-1))
              {  for( l in (k+1):dim)
                 {  if( abs(cov.lag.wgt.mtx.listz[[j]][k,l]) > abs(cov.lag.wgt.mtxz[k,l]))
                    {  cov.lag.wgt.mtxz[k,l]<- cov.lag.wgt.mtxz[l,k] <- cov.lag.wgt.mtx.listz[[j]][k,l]
                       lag.max.cov.mtxz[k,l]<- lag.max.cov.mtxz[l,k] <- lag.input[j]
                    }
                 }
              }
          }
          diag(cov.lag.wgt.mtxz) <-1 
          names<-c()
          for(i in 1:length(dep.var))
          {  name<-depVar[i]
             for(j in 1:length(funcVar))
             {  name2<-paste(name, (funcVar[j]-1))
                names<-c(names, name2)
             }
          }
          dimnames(cov.lag.wgt.mtxz) <-dimnames(lag.max.cov.mtxz)<-
          list(names, names)
          round(cov.lag.wgt.mtxz, 3)
          lag.max.cov.mtxz

       }else 
       {  ## by.deriv.only==TRUE
          cov.lag.mtx.listz <- vector(mode='list', length(vec))
          dim<-length(depVar)
          for(i in 1:length(vec))
          { #cat('subject = ', vec[i])
            cov.lag.mtx.listz[[i]] <-vector(mode='list', length(funcVar))
            ##corrected curves for subject i
            dep.correct<-list()
            for(m in 1:length(depVar))
            {  dep.correct.function<-list()
               for(n in 1:length(funcVar))
               {  dep.correct.function[[n]]<- smoothedCurves[[i]][[m]][[n]][1:max.len]-meanMatrix[[m]][[n]]
               }
               dep.correct[[m]]<-dep.correct.function
            }
            for(deriv in 1:length(funcVar))
            {  cov.lag.mtx.listz[[i]][[deriv]] <- vector(mode='list', length(lag.input))
               for( j in 1:length(lag.input))
               {  #cat(' grid pos. =', j)
                  cov.lag.mtx.listz[[i]][[deriv]][[j]]<- matrix(nrow=dim, ncol=dim)
                  diag(cov.lag.mtx.listz[[i]][[deriv]][[j]])<-1  
                  if( lag.input[j]>=0)
                  {  lag.end <- max.len - lag.input[j]
                     lag.beg <- 1 + lag.input[j]
                     m.support <- lag.end
                     ##standardize each curve according to the lag value and the order of the leading variables, then put into the matrix
                     if(length(byOrder) == 0)
                     {  correst.standz<-list()
                        index<-1
                        mean<- mean(dep.correct[[index]][[deriv]][1:lag.end])
                        sd <- sqrt(var(dep.correct[[index]][[deriv]][1:lag.end]))
                        correst.standz[[index]] <- (dep.correct[[index]][[deriv]][1:lag.end]-mean)/
                                                 (sqrt((m.support-1)/m.support)*sd)
                        for(index in 2:dim )
                        {  if(dim >=2)
                           {  mean<- mean(dep.correct[[index]][[deriv]][lag.beg:max.len])
                              sd<- sqrt(var(dep.correct[[index]][[deriv]][lag.beg:max.len]))
                              correst.standz[[index]] <- (dep.correct[[index]][[deriv]][lag.beg:max.len]-mean)/
                                                  (sqrt((m.support-1)/m.support)*sd)
                              cov.lag.mtx.listz[[i]][[deriv]][[j]][1,index]<-cov.lag.mtx.listz[[i]][[deriv]][[j]][index,1] <-
                                              (1/m.support)*sum(correst.standz[[1]]*correst.standz[[index]])
                           }
                        }
                        for( d in 2:(dim-1)  )
                        {  if( dim >=3)
                           {  mean<- mean(dep.correct[[d]][[deriv]][1:lag.end])
                              sd <- sqrt(var(dep.correct[[d]][[deriv]][1:lag.end]))
                              correst.standz[[d]] <- (dep.correct[[d]][[deriv]][1:lag.end]-mean)/
                                                (sqrt((m.support-1)/m.support)*sd)
                      
                              for( e in (d+1):dim)
                              {  cov.lag.mtx.listz[[i]][[deriv]][[j]][d,e]<-cov.lag.mtx.listz[[i]][[deriv]][[j]][e,d] <-
                                               (1/m.support)*sum(correst.standz[[d]]*correst.standz[[e]])
                              }
                           }
                        }
                        ##end length of byorder=0
                     }else
                     {  ##sort the byorder 
                        byOrderPerD<-c()
                        for(ite in 1: length(byOrder))
                        {  dep<-(byOrder[ite]-1)%/%length(funcVar)+1
                           deriv2<-(byOrder[ite]-1)%%length(funcVar)+1
                           if( deriv2 == deriv)
                           {  byOrderPerD<-c(byOrderPerD, dep)
                           }
                        }
                        ##check
                        if(length(byOrderPerD)!=length(depVar))
                        { stop("**********Error for the length of the order of variables for each derivatives")
                        }
                        correst.standz<-list()
                        index<-byOrderPerD[1]
                        mean<- mean(dep.correct[[index]][[deriv]][1:lag.end])
                        sd <- sqrt(var(dep.correct[[index]][[deriv]][1:lag.end]))
                        correst.standz[[index]] <- (dep.correct[[index]][[deriv]][1:lag.end]-mean)/
                                                (sqrt((m.support-1)/m.support)*sd)
                        for( ord in 2:dim )
                        {  if(dim >=2)
                           {  index<-byOrderPerD[ord]
                              mean<- mean(dep.correct[[index]][[deriv]][lag.beg:max.len])
                              sd<- sqrt(var(dep.correct[[index]][[deriv]][lag.beg:max.len]))
                              correst.standz[[index]] <- (dep.correct[[index]][[deriv]][lag.beg:max.len]-mean)/
                                                  (sqrt((m.support-1)/m.support)*sd)
                              cov.lag.mtx.listz[[i]][[deriv]][[j]][byOrderPerD[1],index]<-cov.lag.mtx.listz[[i]][[deriv]][[j]][index,byOrderPerD[1]] <-
                                          (1/m.support)*sum(correst.standz[[byOrderPerD[1]]]*correst.standz[[index]])
                           }
                        }
                        for( ord in 2:(dim-1)  )
                        {  if( dim >=3)
                           {  d<-byOrderPerD[ord]
                              mean<- mean(dep.correct[[d]][[deriv]][1:lag.end])
                              sd <- sqrt(var(dep.correct[[d]][[deriv]][1:lag.end]))
                              correst.standz[[d]] <- (dep.correct[[d]][[deriv]][1:lag.end]-mean)/
                                                     (sqrt((m.support-1)/m.support)*sd)                     
                              for( ord2 in (ord+1):dim)
                              {  e<-byOrderPerD[ord2]
                                 cov.lag.mtx.listz[[i]][[deriv]][[j]][d,e]<-cov.lag.mtx.listz[[i]][[deriv]][[j]][e,d] <-
                                                                            (1/m.support)*sum(correst.standz[[d]]*correst.standz[[e]])
                              }
                           }
                        }
                     }
                  }else if(lag.input[j]<0)
                  {  lag.end <- max.len+lag.input[j]
                     lag.beg <- 1 - lag.input[j]
                     m.support<- lag.end
                     ##standardize each curve according to the lag value and the order of leading variables specified first, then put into the matrix
                     if(length(byOrder)==0)
                     {  correst.standz<-list()
                        index<-1
                        mean<- mean(dep.correct[[index]][[deriv]][lag.beg:max.len])
                        sd <- sqrt(var(dep.correct[[index]][[deriv]][lag.beg:max.len]))
                        correst.standz[[index]] <- (dep.correct[[index]][[deriv]][lag.beg:max.len]-mean)/
                                                   (sqrt((m.support-1)/m.support)*sd)
                        for(index in 2:dim)
                        {  if(dim>=2)
                           {  mean<- mean(dep.correct[[index]][[deriv]][1:lag.end])
                              sd<- sqrt(var(dep.correct[[index]][[deriv]][1:lag.end]))
                              correst.standz[[index]] <- (dep.correct[[index]][[deriv]][1:lag.end]-mean)/
                                                         (sqrt((m.support-1)/m.support)*sd)
                              cov.lag.mtx.listz[[i]][[deriv]][[j]][1,index]<-cov.lag.mtx.listz[[i]][[deriv]][[j]][index,1] <-
                                                         (1/m.support)*sum(correst.standz[[1]]*correst.standz[[index]])
                           }
                        }
                        for( d in 2:(dim-1))
                        {  if( dim >=3)
                           {  mean<- mean(dep.correct[[d]][[deriv]][lag.beg:max.len])
                              sd <- sqrt(var(dep.correct[[d]][[deriv]][lag.beg:max.len]))
                              correst.standz[[d]] <- (dep.correct[[d]][[deriv]][lag.beg:max.len]-mean)/
                                                     (sqrt((m.support-1)/m.support)*sd)
                              for( e in (d+1):dim)
                              {  cov.lag.mtx.listz[[i]][[deriv]][[j]][d,e]<-cov.lag.mtx.listz[[i]][[deriv]][[j]][e,d] <-
                                                                            (1/m.support)*sum(correst.standz[[d]]*correst.standz[[e]])
                              }
                           }
                        }
                     }else##byOrder not default
                     {  ##sort the byorder 
                        byOrderPerD<-c()
                        for(ite in 1: length(byOrder))
                        {  dep<-(byOrder[ite]-1)%/%length(funcVar)+1
                           deriv2<-(byOrder[ite]-1)%%length(funcVar)+1
                           if(deriv2==deriv)
                           {  byOrderPerD<-c(byOrderPerD, dep)
                           }
                        }
                        ##check
                        if(length(byOrderPerD)!=length(depVar))
                        {  stop("**********Error for the length of the order of variables for each derivatives")
                        }
                        correst.standz<-list()
                        index<-byOrderPerD[1]
                        mean<- mean(dep.correct[[index]][[deriv]][lag.beg:max.len])
                        sd <- sqrt(var(dep.correct[[index]][[deriv]][lag.beg:max.len]))
                        correst.standz[[index]] <- (dep.correct[[index]][[deriv]][lag.beg:max.len]-mean)/
                                                   (sqrt((m.support-1)/m.support)*sd)
                        for(ord in 2:dim)
                        {  index<-byOrderPerD[ord]
                           if(dim>=2)
                           {  mean<- mean(dep.correct[[index]][[deriv]][1:lag.end])
                              sd<- sqrt(var(dep.correct[[index]][[deriv]][1:lag.end]))
                              correst.standz[[index]] <- (dep.correct[[index]][[deriv]][1:lag.end]-mean)/
                                                         (sqrt((m.support-1)/m.support)*sd)
                              cov.lag.mtx.listz[[i]][[deriv]][[j]][byOrderPerD[1],index]<-cov.lag.mtx.listz[[i]][[deriv]][[j]][index,byOrderPerD[1]] <-
                                                          (1/m.support)*sum(correst.standz[[1]]*correst.standz[[index]])
                           }
                        }
                        for( ord in 2:(dim-1))
                        {  if( dim >=3)
                           {  d<-byOrderPerD[ord]
                              mean<- mean(dep.correct[[d]][[deriv]][lag.beg:max.len])
                              sd <- sqrt(var(dep.correct[[d]][[deriv]][lag.beg:max.len]))
                              correst.standz[[d]] <- (dep.correct[[d]][[deriv]][lag.beg:max.len]-mean)/
                                                     (sqrt((m.support-1)/m.support)*sd)
                              for( ord2 in (ord+1):dim)
                              {  e<-byOrderPerD[ord2]
                                 cov.lag.mtx.listz[[i]][[deriv]][[j]][d,e]<-cov.lag.mtx.listz[[i]][[deriv]][[j]][e,d] <-
                                                                            (1/m.support)*sum(correst.standz[[d]]*correst.standz[[e]])
                              }                     
                           }
                        }
                     }##byOrder not default
                  }#end of negative lag if
               }#end of j loop
            }#end of deriv function.choice
         }#end of i loop

         ## calculate overall weighted covariance matrix for each lag;
         ## then determine which lag produces highest magnitude correlation for each unique matrix entry
         ## -> even with lags, we will leave some wieghts as for non-lagged approach
        
         cov.lag.wgt.mtx.listz <- vector(mode='list', length(funcVar))
         for( deriv in 1:length(funcVar))
         {  cov.lag.wgt.mtx.listz[[deriv]]<-vector(mode='list', length(lag.input))
            for( j in 1:length(lag.input))
            {  cov.lag.wgt.mtx.listz[[deriv]][[j]]<-matrix(0, nrow=dim, ncol=dim)
               for( i in 1:length(vec))
               {  cov.lag.wgt.mtx.listz[[deriv]][[j]] <- cov.lag.wgt.mtx.listz[[deriv]][[j]]+(weights.vec[i]*cov.lag.mtx.listz[[i]][[deriv]][[j]])               
               }
               cov.lag.wgt.mtx.listz[[deriv]][[j]]<-(1/sum(weights.vec))*cov.lag.wgt.mtx.listz[[deriv]][[j]]
            }
         }

         ## find highest correlation amongst all lags for each individual 
         cov.lag.wgt.mtxz <-vector(mode='list', length(funcVar))
         lag.max.cov.mtxz<-vector(mode='list', length(funcVar))
         for(deriv in 1:length(funcVar))
         {  cov.lag.wgt.mtxz[[deriv]] <-matrix(0, nrow=dim, ncol=dim)
            lag.max.cov.mtxz[[deriv]]<-matrix(NA, nrow=dim, ncol=dim)
            for( j in 1:length(lag.input))
            {  for( k in 1:(dim-1))
               {  for( l in (k+1):dim)
                  {  if( abs(cov.lag.wgt.mtx.listz[[deriv]][[j]][k,l]) > abs(cov.lag.wgt.mtxz[[deriv]][k,l]))
                     {  cov.lag.wgt.mtxz[[deriv]][k,l]<- cov.lag.wgt.mtxz[[deriv]][l,k] <- cov.lag.wgt.mtx.listz[[deriv]][[j]][k,l]
                        lag.max.cov.mtxz[[deriv]][k,l]<- lag.max.cov.mtxz[[deriv]][l,k] <- lag.input[j]
                     }
                  }
               }
            }
            diag(cov.lag.wgt.mtxz[[deriv]]) <-1
            names<-c()
            for(i in 1:length(dep.var))
            {  name<-depVar[i] 
               name2<-paste(name, (funcVar[deriv]-1))
               names<-c(names, name2)   
            }
            dimnames(cov.lag.wgt.mtxz[[deriv]]) <-dimnames(lag.max.cov.mtxz[[deriv]])<-
            list(names, names)
         }
         cov.lag.wgt.mtxz
         lag.max.cov.mtxz    
      }##end of if statement for by.deriv.only==TRUE      
   }## end of length(lag.input!=0)

      ##############################################################################################################
      # Set up the list to be returned
      ##############################################################################################################
      detach(dataFrame)
      if( by.deriv.only ==FALSE)
      {  if( full.lag.output ==FALSE)
         {  if(length(lag.input)==0)
            {  return(list(dynCorrMatrix=cov.wgt.mtxz))
            }else
            {  return(list(dynCorrMatrix=cov.wgt.mtxz, lag.input=lag.input, max.dynCorr=cov.lag.wgt.mtxz, max.dynCorrLag=lag.max.cov.mtxz))
            }
         }else ## full.lag.output==TRUE
         {  if(length(lag.input)==0)
            {  return(list(dynCorrMatrix=cov.wgt.mtxz))
            }else
            {  dim <- length(depVar)*length(funcVar)
               numRow <- dim *(dim-1)/2
               result <- matrix(0, numRow, length(lag.input))
               namesRow <- c()
               namesCol <- c()
               for( m in 1:length(lag.input))
               { lagName<-paste(lag.input[m])
                 namesCol <- c(namesCol, lagName)
               }
               row <- 1
               for( i in 1:dim)
               {  for( j in (i+1):dim)
                  {  if( dim >= (i+1))
                     {  dep1 <- (i-1)%/%length(funcVar)+1
                        deriv1 <- (i-1)%%length(funcVar)+1
                        dep2 <- (j-1)%/%length(funcVar)+1
                        deriv2 <- (j-1)%%length(funcVar)+1
                        name <- depVar[dep1]
                        name2 <- depVar[dep2]
                        nameT <- paste(name, (funcVar[deriv1]-1), ' x ', name2, (funcVar[deriv2]-1))
                        namesRow <- c(namesRow, nameT)
                        for( m in 1:length(lag.input))
                        {  result[row, m] <- cov.lag.wgt.mtx.listz[[m]][i, j]
                        }
                        row <- row + 1
                     }
                  }
               }
               dimnames(result)<-list(namesRow, namesCol)
               return(list(dynCorrMatrix=cov.wgt.mtxz, lag.input=lag.input, lagResultMatrix=result, max.dynCorr=cov.lag.wgt.mtxz, max.dynCorrLag=lag.max.cov.mtxz))
            }  
         }   
      }else
      {  if( full.lag.output==FALSE)
         {  if(length(lag.input)==0)
            {  return(list(dynCorrMatrix=cov.wgt.mtxz))
            }else
            {  return(list(dynCorrMatrix=cov.wgt.mtxz, lag.input=lag.input, max.dynCorr=cov.lag.wgt.mtxz, max.dynCorrLag=lag.max.cov.mtxz))
            }

         }else
         {  if(length(lag.input)==0)
            {  return(list(dynCorrMatrix=cov.wgt.mtxz))
            }else
            {  dim <- length(depVar)
               numRows <- dim*(dim-1)/2
               namesCol <- c()
               for( m in 1:length(lag.input))
               { lagName<-paste(lag.input[m])
                 namesCol <- c(namesCol, lagName)
               }
               result <- vector(mode='list', length(funcVar))
               for ( k in 1:length(funcVar))
               {  result[[k]]<- matrix(0, numRows, length(lag.input))
                  namesRow <- c()
                  row <- 1
                  for( i in 1:dim)
                  {  for( j in (i+1):dim)
                     {  if( dim >= i+1)
                        {  name <- paste(depVar[i], funcVar[k]-1, 'x', depVar[j], funcVar[k]-1)
                           namesRow <- c(namesRow, name)
                           for ( m in 1:length(lag.input))
                           {  result[[k]][row, m] <- cov.lag.wgt.mtx.listz[[k]][[m]][i, j]
                           } 
                           row <- row + 1
                        }
                     } 
                  }
                  dimnames(result[[k]])<- list(namesRow, namesCol)
               }
               return(list(dynCorrMatrix=cov.wgt.mtxz, lag.input=lag.input, lagResultMatrix=result, max.dynCorr=cov.lag.wgt.mtxz, max.dynCorrLag=lag.max.cov.mtxz))
            }  
         }
      }

}   
            



             
       


       

       