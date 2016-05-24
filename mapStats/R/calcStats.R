calcStats <- function (d, var, d.geo.var, stat=c("mean","quantile"), quantiles=c(.5,.75),  
                       by.var=NULL, wt.var = NULL, cell.min = 0) 
    {

    
   stat <- match.arg(tolower(stat), c("mean", "quantile", "total","var", "sd"), several.ok = TRUE)

   #check if all variables defined
   
   for (g in c(var, d.geo.var, by.var, wt.var)) { if (! g %in% colnames(d)) { stop (paste("Variable",g,"not defined.")) } }
   
   #list to contain all variable statistics lists
   all_stats <- list()

    #create dummy weight vector if no weight vector specified
    if (is.null(wt.var)) {
          d$weight.vector.default <- 1
          wt.var <- "weight.vector.default"
     }

    #make the by variables factors just in case
    for (k in c(d.geo.var, by.var)) {

           d[, k] <- as.factor(d[, k])      
 
    }

    #check if there is only one level of the by-variable:
    
    if (! is.null(by.var)) {
	#number of factor levels for each by variable
        len <- unlist(with(d, lapply(by.var, function(x) { length(levels(d[,x])) }))) 
         
        if (all(len==1)) { 
            by.var <- NULL 
            warning("None of the by-variables have more than one level.\nAnalysis will continue without by-variables.\n")
            }
                   
        else if (any( len ==1)) {

            #if some (but not all) of the by-variables have length one

            by.var <- by.var[ len > 1 ]
            warning(paste("The following by-variables are omitted from analyis\nbecause they have only one level:",
                      paste(by.var[ len == 1 ], collapse=", "), sep="\n"))
              
            }

     }#finish checking by-var


                


   
    
    by.form <- as.formula( paste("~", paste( c(d.geo.var, by.var), collapse=" + ")))

    #create formulas for weights and analysis variable
    wt.form <- as.formula(paste("~", wt.var))



    #option for multiple variables:

    for (x in var) {

      

       #list to contain matrices of statistics
       v.stats <- list()


       #eliminate observations with missing values for the variables
       #the fact that the by-variables were made factors before keeps 
       #the levels even though we may have missing combinations

       d_tmp <- na.omit(d[, c(wt.var, x, d.geo.var, by.var)])

       var.form <- as.formula(paste("~", x))

       #now recast the resulting matrix  so
       #have one row per d.geo.var level for merging  
                  
       if (! is.null(by.var)) { cast.form <- as.formula(paste(d.geo.var, paste(by.var, collapse=" + "), sep = "~")) }  

       #now define structure for the data design
       data.des <- svydesign(ids = ~1, data = d_tmp, weights = wt.form)

       #here we create names for the columns
       #this is necessary to avoid common characters like . or _

  
       if (is.null(by.var)) { cnames <- c(d.geo.var, x) }
       else {
          levs <- list()
          for (v in rev(by.var)) { levs[[v]] <- levels(d_tmp[,v]) }
          cnames <- as.data.frame(expand.grid(levs)[,length(levs):1])
          cnames <- apply(cnames, 1, paste, collapse="__.__")
          cnames <- c(d.geo.var, cnames)

         }  


      #now loop through each statistic and calculate

       for (k in stat) {
           svyfunc <- match.fun(paste("svy", ifelse(k %in% c("sd","var"), "var", k), sep = ""))
        
           #create matrix of statistics

       
           if (k == "quantile") {
            
           #extra processing for quantiles
            stat.mat <- svyby(formula = var.form, by = by.form, design = data.des, 
                              FUN = svyfunc, keep.var = FALSE, na.rm = TRUE, drop.empty.groups = FALSE,
                              quantiles=quantiles, method="constant")
            
            #create statistic names
            
            qnames <- paste("Q",100*quantiles, sep="")
            qnames[ qnames=="Q0" ] <- "Minimum"
            qnames[ qnames=="Q50" ] <- "Median"
            qnames[ qnames=="Q100" ] <- "Maximum"
            
            #break out each quantile from matrix
              for (q in 1:length(quantiles)) {
                  qn <- ifelse(length(quantiles)>1, paste("statistic", q, sep=""), "statistic")
              
                 mat <- stat.mat[, c(d.geo.var, by.var, qn)]
                 if ( ! is.null(by.var)) { mat <-  reshape2::dcast(data=mat, formula=cast.form, value.var=qn) }                     
                
                  colnames(mat) <- cnames
                  v.stats[[ qnames[q] ]] <- mat

               }

            }
          
           #mean/total/var/sd

           else {

              mat <- svyby(formula = var.form, by = by.form, design = data.des, 
                                FUN = svyfunc, keep.var = FALSE, na.rm = TRUE, drop.empty.groups = FALSE)
            
                 
             if (! is.null(by.var)) { mat <- reshape2::dcast(data=mat, formula=cast.form, value.var="statistic") }
             colnames(mat) <- cnames

              #standard deviation

              if (k == "sd" ) { mat[ , 2:ncol(mat) ] <- sqrt(mat[ , 2:ncol(mat) ] )  }
               
              v.stats[[Hmisc::capitalize( ifelse(k=="var", "Variance", ifelse(k=="sd", "Standard deviation", k))  )]] <- mat

              }

          }#finish looping over statistics.  Now some extra processing


           #table of frequency counts

           t.c <- as.data.frame(table(d_tmp[, c(d.geo.var, rev(by.var)) ])) 
           if (! is.null(by.var)) {
                      
              t.c <- reshape2::dcast(data=t.c, formula=cast.form, value.var="Freq")
            }
         

           #now loop over all statistics and make NA all the ones with NAs
	     if (cell.min >0 ) {
              colrange <- 2:ncol(t.c)
              NAfreqs <- (t.c[, colrange ] < cell.min)
           
              for (k in names(v.stats)) {
                  
                  mat <- v.stats[[ k ]]
                  mat[ , colrange ][ NAfreqs ] <- NA
                  v.stats[[ k ]][ , colrange ] <- mat[ , colrange ]
                }
              }
               
           
           v.stats[["Freqs"]] <- t.c  
           attributes(v.stats)$variable <- x
           all_stats[[ x ]] <- v.stats

        }#finish looping over variables


    all_stats  
    }




     
