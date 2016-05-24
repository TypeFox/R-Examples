

bootstrapCI<- function(dataFrame,
				   depVar= c('resp1', 'resp2', 'resp3', 'resp4', 'resp5'),
                           indepVar='time',
                           subjectVar= 'subject',                           
				   function.choice =c(1), 
                           width.range=c(((range(dataFrame[[indepVar]]))[2]-(range(dataFrame[[indepVar]]))[1])/4, ((range(dataFrame[[indepVar]]))[2]-(range(dataFrame[[indepVar]]))[1])/4),
			         width.place=c(NA, NA), 
				   boundary.trunc= c(0, 0),
                           byOrder=c(),
                           max.dynCorrLag=0,  
                           B=100, 
                           percentile=c(0.025, 0.975),
                           by.deriv.only = TRUE,
                           seed = 299)

	 ## arguments are:
       ##               dataFrame       ( the dataFrame that contains the dataset )
	 ##			depVar          ( dependents or y-axis variables )
	 ##			indepVar        ( independent or x-axis variable, usually the recording time for the data )
       ##               subjectVar      ( column name of the individuals )
	 ##			function.choice ( a vector of length 3 to indicate which derivatives
	 ##				            are desired for each curve smoothing; 1st is 0th derivative, 2nd is 1st
	 ##				            derivative, 3rd is 2nd derivative, etc; 1=yes, 0=no.
       ##                                 e.g. c(1,0,1) means we would like to look at derivative 0 and derivative 2 )	
	 ##			width.range     ( bandwidth (size) for local regression in curve smoothing, it can be a list that specifies different bandwidth 
       ##                                 vector for each response, or a single width vector that is used among all the responses )	 
	 ##			width.place     ( endpoints for width change, it can be a list that specifies difference endpoints for each response, or a single
       ##                                 vector of endpoints that is used among all the responses )
	 ##		      boundary.trunc  ( indicate the boundary of indep.var that should be truncated after smoothing )
       ##               byOrder         ( a vector that specifies the order of the variables and derivatives to be the leading variable in the calculation )
       ##               max.dynCorrLag       ( a specified lag value that we would like to analyze at )
       ##               B               ( The number of samples used for the bootstrap )
       ##               percentile      ( The percentile used to construct confidence intervals for dynamic correlations )
       ##               by.deriv.only   ( If true, the inter-dynamic correlations between different derivatives are not computed, which 
       ##                                 can save computation time )
       ##               seed            ( The seed used in generating the samples )
{
     ## set up the variables
     attach(dataFrame)
     dep.var= dataFrame[depVar]   
     indep.var=dataFrame[[indepVar]]
     subject.var=dataFrame[[subjectVar]]
     vec <-  unique(subject.var)
  
     #maybe used to find the number of days to be used later
     min<-c()
     for(i in 1:length(vec))
     { min <- c(min, max(indep.var[subject.var==vec[i]]))
     }
     limit<- min(min)
     
     ## To check which derivatives to be considered
     funcVar<-c()
     for(i in 1:length(function.choice))
     {  if(function.choice[i]==1)
        { funcVar<-c(funcVar, i)
        }
     }

      ## check parameters
      if(length(byOrder)!=0)
      {  largest<-length(depVar)*length(funcVar)
         if(max(byOrder)>largest)
         { stop("The specified index order of leading variables is over range")
         }
         if(length(byOrder)!=largest)
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
  
      ## Do smoothing on the curves
      smoothedCurves <-vector(mode='list', length(vec))
      for(k in 1:length(vec))
      {  # cat('subject=',vec[k], ' ')
         indep<- indep.var[subject.var==vec[k]]
         smoothedCurves[[k]]<-vector(mode='list', length(depVar))
         for(i in 1:length(depVar))
         {  if( is.na(width.place[[i]][1]))
            {  width.place[[i]][1]<-min(indep[!is.na(indep)])
            }
            if ( is.na(width.place[[i]][2]))
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

           temp<- lpepa(x=indep.var[subject.var==vec[k]][!is.na(dep.var[[i]][subject.var==vec[k]])], y=dep.var[[i]][subject.var==vec[k]][!is.na(dep.var[[i]][subject.var==vec[k]])],
                     bandwidth=band,  deriv = 0, n.out=size,  x.out= seq(boundary.trunc[1], ceiling(max(indep[!is.na(indep)])-boundary.trunc[2])),
                     order=1, var= FALSE)$est
           smoothedCurves[[k]][[i]] <- vector(mode='list', length=2)
           smoothedCurves[[k]][[i]][[1]]<-seq(boundary.trunc[1], ceiling(max(indep[!is.na(indep)])-boundary.trunc[2]))       
           smoothedCurves[[k]][[i]][[2]]<-temp                    
        }
    }  
    
     ## Compute the residuals
     days.out<-boundary.trunc[1]:(ceiling(limit)-boundary.trunc[2])
     resid<-vector(mode='list', length(vec))
     smooth<-vector(mode='list', length(vec))

     for( k in 1:length(vec))
     {  resid[[k]]<-vector(mode='list', length(depVar))
        smooth[[k]]<-vector(mode='list', length(depVar))
        
        for(i in 1:length(depVar))
        {  indep<- indep.var[subject.var==vec[k]][!is.na(dep.var[[i]][subject.var==vec[k]])]   
           temp.obs<- indep[indep>=boundary.trunc[1]& indep<=ceiling(max(indep))-boundary.trunc[2]]
           dep.v<- dep.var[[i]][subject.var==vec[k]][!is.na(dep.var[[i]][subject.var==vec[k]])]
           depValue<-dep.v[indep>=boundary.trunc[1]& indep<=ceiling(max(indep[!is.na(indep)]))-boundary.trunc[2]]
           temp.obs<-round(temp.obs, 0)
           ##to form the smoothed curves according to temp.obs
           smooth[[k]][[i]] <- vector(mode='numeric', length(temp.obs))
           for(j in 1:length(temp.obs))
           {  smooth[[k]][[i]][[j]] <-  smoothedCurves[[k]][[i]][[2]][is.element(smoothedCurves[[k]][[i]][[1]], temp.obs[j])]
           }   
           resid[[k]][[i]]<- depValue - smooth[[k]][[i]]              
        }
     }

     rm(temp.obs)
     
     ### step 2: get random sample (w/ replacement) of subjects, then, within each sample, 
     ###				get random sample (w/ replacement) of residuals, then 
     ###				re-create "original" data using original smoothed curves and these residuals
     ### step 3: repeat step 2 B times (hopefully, we can have B = 1000)
   
     ## a list of all the samples for each of B bootstraps
     bs.seq.samp<-vector(mode='list', length=B)
     set.seed(seed)
     indexVec<-1:length(vec)
     for(b in 1:B)
     { bs.seq.samp[[b]]<-sort(sample(indexVec, replace=TRUE))
     }
     ##recreate list after sampling with replacement on the time stamps
     recreate.B.list<-vector(mode='list', length(B))
     for(b in 1:B)
     {  # cat(' B is ', b)
        vecSample<-bs.seq.samp[[b]] 
        recreate.B.list[[b]]<-vector(mode='list', length(vecSample))
        for(j in 1:length(vecSample))
        {  recreate.B.list[[b]][[j]]<-vector(mode='list', length(depVar))
           for(i in 1:length(depVar))
           {  recreate.B.list[[b]][[j]][[i]]<- sample(resid[[vecSample[j]]][[i]], replace=TRUE)+ smooth[[vecSample[j]]][[i]]
           }  
        }
     }

     ### step 4, smoothing of recreated data; note that I will append this recreated data
     # 		with original data from days 0 to 25 and days 206 to 230 simply to borrow that 
     #		information for smoothing in the interior  

     days.out<-boundary.trunc[1]:(ceiling(limit)-boundary.trunc[2])
     temp.days<-tempData<-vector(mode='list', length(B))
      
     ## recreate data
     for(b in 1:B)
     {  # cat(' bs = ', b, '\n')
        vecSample<-bs.seq.samp[[b]]
        tempData[[b]]<-vector(mode='list', length(vecSample))
        temp.days[[b]]<-vector(mode='list', length(vecSample))

        for(j in 1:length(vecSample))
        {  tempData[[b]][[j]]<-vector(mode='list', length(depVar))
           temp.days[[b]][[j]] <- vector(mode='list', length(depVar))
           for(i in 1:length(depVar))
           {  # cat(' depVar = ', depVar[i], '\n') 
              indep<- indep.var[subject.var==vec[vecSample[j]]][!is.na(dep.var[[i]][subject.var==vec[vecSample[j]]])]
              temp.1<- indep[indep>=0& indep< boundary.trunc[1]]
              temp.2<- indep[indep>ceiling(max(indep[!is.na(indep)]))-boundary.trunc[2]& indep<=max(indep[!is.na(indep)])]
              temp.days[[b]][[j]][[i]]<- indep[indep>=0& indep<=max(indep[!is.na(indep)])]
            
              depValue<-dep.var[[i]][subject.var==vec[vecSample[j]]]
              depV1<-depValue[is.element(indep.var[subject.var==vec[vecSample[j]]], temp.1)]
              depV2<-depValue[is.element(indep.var[subject.var==vec[vecSample[j]]], temp.2)]
              if(length(temp.2)>0)
              {  tempData[[b]][[j]][[i]]<- c(depV1, recreate.B.list[[b]][[j]][[i]], depV2)
              }else if(length(temp.1)>0)
              {  tempData[[b]][[j]][[i]]<-c(depV1, recreate.B.list[[b]][[j]][[i]])
              }else
              {  tempData[[b]][[j]][[i]]<-recreate.B.list[[b]][[j]][[i]]
              }
           }##depVar
        }##vecSample
     }##B
    
     # perform smoothing and get correlation matrix for each B and subject
     len<-length(vec)
     cov.wgt.mtx.B <- vector(mode='list', length=B)
     
     ###############################################################################
     # Compute Dynamic Correlation matrix for each B
     ###############################################################################
     for (b in 1:B)
     {  cat(' B =', b, '\n')
        vecSample<-bs.seq.samp[[b]]
        smoothedCurves2 <-vector(mode='list', length(vecSample)) 
        for(k in 1:length(vecSample))
        {  # cat('subject=', vec[vecSample[k]], ' ')
           
           smoothedCurves2[[k]] <-vector(mode='list', length(depVar))     
           for(i in 1:length(depVar))
           {  # cat('Dependent variable = ', depVar[i], '\n')
             indep<- temp.days[[b]][[k]][[i]]
             if( is.na(width.place[[i]][1]))
             {  width.place[[i]][1]<-min(indep[!is.na(indep)])
             }
             if ( is.na(width.place[[i]][2]))
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
             smoothedCurves2[[k]][[i]] <-vector(mode='list', length(funcVar))
             for( j in 1: length(funcVar))
              {  # cat('deriative function index= ', funcVar[j], '\n')
                 tempFunc<- lpepa(x=indep, y=tempData[[b]][[k]][[i]], 
                           bandwidth=band,  deriv =(funcVar[j]-1), n.out=size,  x.out= seq(boundary.trunc[1], (ceiling(max(indep[!is.na(indep)]))-boundary.trunc[2])),
                           order=funcVar[j], var= FALSE)$est
                 smoothedCurves2[[k]][[i]][[j]]<-tempFunc
              }
           }   
        }
    
        ## smoothedCurves
        ################################
        ######  Obtain column means for responses
     
        max.len<- (ceiling(limit)-boundary.trunc[2])-boundary.trunc[1]+1
        meanMatrix<- list() 
        for(i in 1:length(depVar))
        {  forEachDepMean<-list()
           for( j in 1:length(funcVar))
           {  ##cat('deriative function index= ', funcVar[j], '\n')
              forEachDerivMean<-rep(NA, max.len)
              for( m in 1:max.len)
              {  for( n in 1:length(vecSample))
                 {  if( n==1)
                    {  total <-1 
                       forEachDerivMean[m] <- smoothedCurves2[[n]][[i]][[j]][m]
                    }else if(n>1)
                    {  total<- total+1
                       forEachDerivMean[m] <- (1/total)*((total-1)* forEachDerivMean[m] +  smoothedCurves2[[n]][[i]][[j]][m])
                    }
                 }
              }  
              forEachDepMean[[j]]<-forEachDerivMean
           }  
           meanMatrix[[i]]<-forEachDepMean
        }     

        if(by.deriv.only == FALSE)
        {  dim<-length(depVar)*length(funcVar)
           cov.lag.mtx.listz <-list()
           weights.vec <- rep(NA, length(vecSample))
           time.extend <- (boundary.trunc[1] + max.len - 1) + boundary.trunc[1] 
           for(i in 1:length(vecSample))
           {  # cat('subject = ', vec[vecSample[i]])
              cov.lag.mtx.listz[[i]] <-list()
              # use days with actual observations for weights
              weights.vec[i] <- length(subject.var[(subject.var==vec[vecSample[i]]) & 
							(indep.var >= 0) & (indep.var <= time.extend)])
              ##corrected curves for subject i
              dep.correct<-list()
              for(m in 1:length(depVar))
              {  dep.correct.function<-list()
                 for(n in 1:length(funcVar))
                 { dep.correct.function[[n]]<- smoothedCurves2[[i]][[m]][[n]][1:max.len]-meanMatrix[[m]][[n]]
                 }
                 dep.correct[[m]]<-dep.correct.function
              }
           
              cov.lag.mtx.listz[[i]]<- matrix(nrow=dim, ncol=dim)
              diag(cov.lag.mtx.listz[[i]])<-1
              if( max.dynCorrLag >= 0)
              {   lag.end <- max.len -  max.dynCorrLag
                  lag.beg <- 1 +  max.dynCorrLag
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
                           cov.lag.mtx.listz[[i]][1,index]<-cov.lag.mtx.listz[[i]][index,1] <-
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
                          {  cov.lag.mtx.listz[[i]][d,e]<-cov.lag.mtx.listz[[i]][e,d] <-
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
                           cov.lag.mtx.listz[[i]][byOrder[1],index]<-cov.lag.mtx.listz[[i]][index,byOrder[1]] <-
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
                              cov.lag.mtx.listz[[i]][d,e]<-cov.lag.mtx.listz[[i]][e,d] <-
                                          (1/m.support)*sum(correst.standz[[d]]*correst.standz[[e]])
                           }
                        }
                     }
                  }#end of byorder is specified
                 
               }else if( max.dynCorrLag<0)
               {  lag.end <- max.len+ max.dynCorrLag
                  lag.beg <- 1 -  max.dynCorrLag
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
                           cov.lag.mtx.listz[[i]][1,index]<-cov.lag.mtx.listz[[i]][index,1] <-
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
                           {  cov.lag.mtx.listz[[i]][d,e]<-cov.lag.mtx.listz[[i]][e,d] <-
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
                           cov.lag.mtx.listz[[i]][1,index]<-cov.lag.mtx.listz[[i]][index,1] <-
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
                              cov.lag.mtx.listz[[i]][d,e]<-cov.lag.mtx.listz[[i]][e,d] <-
                                           (1/m.support)*sum(correst.standz[[d]]*correst.standz[[e]])
                           }
                      
                        }
                     }  

                  }##byorder!=0
                }#end of negative if
          }#end of i loop
          
           
          ## calculate overall weighted covariance matrix for each lag;
          ## then determine which lag produces highest magnitude correlation for each unique matrix entry
          ## -> even with lags, we will leave some wieghts as for non-lagged approach
        
          cov.lag.wgt.mtx<- matrix(0, nrow=dim, ncol=dim)
          for( i in 1:length(vecSample))
          {  cov.lag.wgt.mtx <- cov.lag.wgt.mtx+(weights.vec[i]*cov.lag.mtx.listz[[i]])
          }
          cov.lag.wgt.mtx<-(1/sum(weights.vec))*cov.lag.wgt.mtx

          diag(cov.lag.wgt.mtx) <-1 
          names<-c()
          for(i in 1:length(dep.var))
          {  name<-depVar[i]
             for(j in 1:length(funcVar))
             {  name2<-paste(name, (funcVar[j]-1))
                names<-c(names, name2)
             }
          }
          dimnames(cov.lag.wgt.mtx)<-list(names, names)
          cov.wgt.mtx.B[[b]]<-cov.lag.wgt.mtx
        

       }else  
       {  ## by.deriv.only==TRUE
          weights.vec <- rep(NA, length(vecSample))
          time.extend <- (boundary.trunc[1] + max.len - 1) + boundary.trunc[1] 
          cov.lag.mtx.listz <- vector(mode='list', length(vec))
          dim<-length(depVar)
          for(i in 1:length(vecSample))
          { # cat('subject = ', vec[vecSample[i]])
            # use days with actual observations for weights
            weights.vec[i] <- length(subject.var[(subject.var==vec[vecSample[i]]) & 
							(indep.var >= 0) & (indep.var <= time.extend)]) 
            cov.lag.mtx.listz[[i]] <-vector(mode='list', length(funcVar))
            ##corrected curves for subject i
            dep.correct<-list()
            for(m in 1:length(depVar))
            {  dep.correct.function<-list()
               for(n in 1:length(funcVar))
               {  dep.correct.function[[n]]<- smoothedCurves2[[i]][[m]][[n]][1:max.len]-meanMatrix[[m]][[n]]
               }
               dep.correct[[m]]<-dep.correct.function
            }
            for(deriv in 1:length(funcVar))
            {     cov.lag.mtx.listz[[i]][[deriv]]<- matrix(nrow=dim, ncol=dim)
                  diag(cov.lag.mtx.listz[[i]][[deriv]])<-1  
                  if(  max.dynCorrLag>=0)
                  {  lag.end <- max.len -  max.dynCorrLag
                     lag.beg <- 1 +  max.dynCorrLag
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
                              cov.lag.mtx.listz[[i]][[deriv]][1,index]<-cov.lag.mtx.listz[[i]][[deriv]][index,1] <-
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
                              {  cov.lag.mtx.listz[[i]][[deriv]][d,e]<-cov.lag.mtx.listz[[i]][[deriv]][e,d] <-
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
                              cov.lag.mtx.listz[[i]][[deriv]][byOrderPerD[1],index]<-cov.lag.mtx.listz[[i]][[deriv]][index,byOrderPerD[1]] <-
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
                                 cov.lag.mtx.listz[[i]][[deriv]][d,e]<-cov.lag.mtx.listz[[i]][[deriv]][e,d] <-
                                                                            (1/m.support)*sum(correst.standz[[d]]*correst.standz[[e]])
                              }
                           }
                        }
                     }
                  }else if( max.dynCorrLag<0)
                  {  lag.end <- max.len+ max.dynCorrLag
                     lag.beg <- 1 -  max.dynCorrLag
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
                              cov.lag.mtx.listz[[i]][[deriv]][1,index]<-cov.lag.mtx.listz[[i]][[deriv]][index,1] <-
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
                              {  cov.lag.mtx.listz[[i]][[deriv]][d,e]<-cov.lag.mtx.listz[[i]][[deriv]][e,d] <-
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
                              cov.lag.mtx.listz[[i]][[deriv]][byOrderPerD[1],index]<-cov.lag.mtx.listz[[i]][[deriv]][index,byOrderPerD[1]] <-
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
                                 cov.lag.mtx.listz[[i]][[deriv]][d,e]<-cov.lag.mtx.listz[[i]][[deriv]][e,d] <-
                                                                            (1/m.support)*sum(correst.standz[[d]]*correst.standz[[e]])
                              }                     
                           }
                        }
                     }##byOrder not default
                  }#end of negative lag if
            }#end of deriv function.choice
         }#end of i loop

         ## calculate overall weighted covariance matrix for each lag;
         ## then determine which lag produces highest magnitude correlation for each unique matrix entry
         ## -> even with lags, we will leave some wieghts as for non-lagged approach
        
         cov.lag.wgt.mtx<- vector(mode='list', length(funcVar))
         for( deriv in 1:length(funcVar))
         {  cov.lag.wgt.mtx[[deriv]]<-matrix(0, nrow=dim, ncol=dim)
            for( i in 1:length(vecSample))
            {  cov.lag.wgt.mtx[[deriv]] <- cov.lag.wgt.mtx[[deriv]]+(weights.vec[i]*cov.lag.mtx.listz[[i]][[deriv]])               
            }
            cov.lag.wgt.mtx[[deriv]]<-(1/sum(weights.vec))*cov.lag.wgt.mtx[[deriv]]
     
         }

         for(deriv in 1:length(funcVar))
         {  names<-c()
            for(i in 1:length(dep.var))
            {  name<-depVar[i] 
               name2<-paste(name, (funcVar[deriv]-1))
               names<-c(names, name2)   
            }
            dimnames(cov.lag.wgt.mtx[[deriv]]) <-list(names, names)
         }
         cov.wgt.mtx.B[[b]]<-cov.lag.wgt.mtx
      }##end by.deriv.only=TRUE 
   }##end b


   if(by.deriv.only==FALSE)
   {  # first, calculate mean correlation across BS samples
      dim<-length(depVar)*length(funcVar)
      mean.cov.wgt.mtx.B <- matrix(0, nrow=dim, ncol=dim) 
      for(b in 1:B)
      { mean.cov.wgt.mtx.B <- mean.cov.wgt.mtx.B + cov.wgt.mtx.B[[b]]
      }
      mean.cov.wgt.mtx.B <- mean.cov.wgt.mtx.B / B
   }else
   {  ##by.deriv.only==TRUE
      mean.cov.wgt.mtx.B <- vector(mode='list', length(funcVar))
      dim<-length(depVar)
      for(deriv in 1:length(funcVar))
      {  mean.cov.wgt.mtx.B[[deriv]]<-matrix(0, nrow=dim, ncol=dim)
         for(b in 1:B)
         {  mean.cov.wgt.mtx.B[[deriv]] <- mean.cov.wgt.mtx.B[[deriv]] + cov.wgt.mtx.B[[b]][[deriv]]
         }
         mean.cov.wgt.mtx.B[[deriv]] <- mean.cov.wgt.mtx.B[[deriv]]/B
      } 
   }  
   mean.cov.wgt.mtx.B  

   ### note, as would seem expected, BS mean correlations are downward biased (towards 0)
   ###		from the true observed correlations (except for cer-trf pair), 
   ###		almost surely due to the intra-residual  
   ###		BS adjustment that would cause the correlations to not be quite as strong
   ###		in absolute value 
   # second, generate ordered (empirical) distribution of correlations,  
   #	and obtain BS interval estimates based on percentile method
   detach(dataFrame)
   if(by.deriv.only==FALSE)
   {  dim <- length(depVar)*length(funcVar)
      corr<-vector(mode='list', dim)
      for( i in 1:dim)
      {  corr[[i]]<-vector(mode='list', length(dim))
         for( j in 1:dim)
         {  corr[[i]][[j]]<-c(cov.wgt.mtx.B[[1]][i,j])
            for( b in 2:B)
            {  corr[[i]][[j]]<- c(corr[[i]][[j]], cov.wgt.mtx.B[[b]][i,j])
            }
         }
      }
      CI<-matrix(0, dim*(dim-1)/2, 5) 
      count<-1
      namesRow<-c()
      for( i in 1:dim)
      {  for ( j in (i+1):dim)
         {  if( dim >= (i+1) )
            {  rowVar <- (i-1)%/%length(funcVar)+1
               rowDer <- (i-1)%%length(funcVar)+1
               colVar <- (j-1)%/%length(funcVar)+1
               colDer <- (j-1)%%length(funcVar)+1
               temp<- quantile(corr[[i]][[j]], percentile)
               CI[count,3] <- round(temp[[1]], 4)
               CI[count,4] <- round(temp[[2]], 4)
               CI[count,2] <- max.dynCorrLag
               CI[count,1] <- " "
               if( mean.cov.wgt.mtx.B[i,j] > 0 )
               {  CI[count,5] <- round(2*sum(corr[[i]][[j]]<=0)/B,5)
               }else
               {  CI[count,5] <- round(2*sum(corr[[i]][[j]]>=0)/B,5)
               }
               count <- count + 1
               name<-paste(depVar[rowVar], funcVar[rowDer]-1, 'x', depVar[colVar], funcVar[colDer]-1)
               namesRow <-c(namesRow, name)
            }
         }  
      } 
      perc1 <- paste(percentile[1]*100, '%')
      perc2 <- paste(percentile[2]*100, '%')
      namesCol <- c(' ', 'lag', perc1, perc2, 'BS p-value')
      dimnames(CI) <- list(namesRow, namesCol)
      CI<-as.data.frame(CI)
      return(list(quantilesMatrix = CI))
   }else
   {  ## by.deriv.only=TRUE
      dim <- length(depVar)
      corr <- vector(mode='list', length(funcVar))
      for( deriv in 1:length(funcVar))
      {  corr[[deriv]]<-vector(mode='list', dim)
         for( i in 1:dim)
         {  corr[[deriv]][[i]]<-vector(mode='list', length(dim))
            for( j in 1:dim)
            {  corr[[deriv]][[i]][[j]] <-c(cov.wgt.mtx.B[[1]][[deriv]][i,j])
               for(b in 2:B)
               {  corr[[deriv]][[i]][[j]] <- c(corr[[deriv]][[i]][[j]], cov.wgt.mtx.B[[b]][[deriv]][i,j])
               }
            }
         }
      }
      CI<-vector(mode='list', length(funcVar))
      perc1 <- paste(percentile[1]*100, '%')
      perc2 <- paste(percentile[2]*100, '%')
      namesCol <- c(' ', 'lag', perc1, perc2, 'BS p-value')
      for( deriv in 1:length(funcVar))
      {  CI[[deriv]]<-matrix(0, dim*(dim-1)/2, 5)
         count<-1
         namesRow<-c()
         for( i in 1:dim)
         {  for( j in (i+1):dim)
            {  if( dim >= (i+1))
               {  temp <- quantile(corr[[deriv]][[i]][[j]], percentile)
                  CI[[deriv]][count, 3] <- round(temp[[1]], 4)
                  CI[[deriv]][count, 4] <- round(temp[[2]], 4)
                  CI[[deriv]][count, 2] <- max.dynCorrLag
                  CI[[deriv]][count, 1] <- " "
                  if( mean.cov.wgt.mtx.B[[deriv]][i,j] > 0 )
                  {  CI[[deriv]][count, 5] <- round(2*sum(corr[[deriv]][[i]][[j]]<=0)/B,5)
                  }else
                  {  CI[[deriv]][count, 5] <- round(2*sum(corr[[deriv]][[i]][[j]]>=0)/B,5)
                  }
                  count <- count + 1
                  name <- paste(depVar[i], funcVar[deriv]-1, 'x', depVar[j], funcVar[deriv]-1)
                  namesRow <- c(namesRow, name)
               }
            }
         }
         dimnames(CI[[deriv]]) <- list(namesRow, namesCol)
         CI[[deriv]]<-as.data.frame(CI[[deriv]])
      }
      return(list(quantilesMatrix = CI))
   }  
}  