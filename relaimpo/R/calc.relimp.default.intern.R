"calc.relimp.default.intern" <-
function (object, x = NULL, ..., type = "lmg", diff = FALSE, rank = TRUE, rela = FALSE, always = NULL, 
        groups = NULL, groupnames=NULL, weights=NULL, design=NULL, WW=NULL, test=FALSE, ynam=NULL, ngroups=NULL) 
{
## ngroups is used for obtaining adequate ordering of effects in case of interactions including
## factors (or groups)
## uses the original numbering of effects

    # Author and copyright holder: Ulrike Groemping

    #This routine is distributed under GPL version 2 or newer.
    #The text of this license can be found at http://www.gnu.org/copyleft/gpl.html.

    # object can be one of the following:
    #  - the variance-covariance matrix of y with all regressor variables (y in first position),
    #    which needs to be square, at least 3x3 (= 2 regressors), and positive definite
    #   - a response variable, 
    #   - a data frame with the response in the first column (no factors allowed)
    #   - a matrix with the response in the first column 
    # x is
    #   - a matrix or data frame of x-variables if y is a response variable 
    #   - NULL otherwise
    # type is a character listing (vector incl. single character string, list, matrix) 
    #    that chooses one or more of the available relative importance metrics (in any order)
    #    The selected metrics will always appear in the result in the order given below:
    #    "lmg","pmvd","last","first","betasq","pratt" (without "pmvd" for the global version)


    # diff is a logical requesting differences between metrics if TRUE
    # rank is a logical requesting ranks of metrics if TRUE (largest=1, smallest=p)
    # rela is a logical requesting normalization to a percentage scale (sum 100%) if TRUE
        # if FALSE, all metrics are relative to var(Y), lmg, pmvd and pratt sum to R^2, 
        # first, last and betasq do not sum to R^2 but neither to 100%
    # always is a vector of integer numbers or a character vector given column positions 
    #    of regressors that are always kept in the model, i.e. are adjusted for in all calculations
    #    numbers refer to the position of the regressors in the variable list (1 refers to response),
    #    character strings to variable names
    # groups gives variable combinations that are to be treated as a group
    # groupnames gives a name for each variable combination
    # weights gives a vector of observation weights (all equally weighted, if NULL)
    # design gives a design according to package survey (either weights OR design only)
    # WW is determined from calc.relimp.lm in case of interactions (German: WechselWirkungen)
    # test is set to TRUE by routines that use calc.relimp.default.intern for error checking and preparation
    # ynam is used by calc.relimp.default in case data are given as y and x



    y<-object
    ogroups <- groups
    ogroupnames <- groupnames
    xcall <- x
    
    hasfacs <- FALSE
    if (is.null(x) && is.data.frame(y)) hasfacs <- any(sapply(y,is.factor))
    if (!is.null(x) && is.data.frame(x)) hasfacs <- any(sapply(x,is.factor))
    
    #error control and initialization

    nobs <- NULL
    
    if (!is.logical(rank)) 
        stop("rank must be a logical")
    if (!is.logical(diff)) 
        stop("diff must be a logical")
    if (!is.logical(rela)) 
        stop("rela must be a logical")

    # alltype is set in zzz.R
    alltype <- alltype()
    if (!all(type %in% alltype) && (length(alltype) == 6 || !("pmvd" %in% 
        type))) 
        stop("invalid type requested")
    if (!all(type %in% alltype) && length(alltype) == 5 && "pmvd" %in% 
        type) 
        stop("pmvd is not a valid type in this version of relaimpo, obtain the non-US version, if you want to use pmvd")

    ## interim value of p, since numbers in always and groups would refer to these columns
    ## in case of factors, this is not identical to the later p
    if (is.null(x) && (is.data.frame(y) || is.matrix(y))) 
    p <- ncol(y)-1 
    if (!is.null(x) && (is.data.frame(x) || is.matrix(x)))
    p <- ncol(x) 
                   

    if (!is.null(always) && !length(always) == length(unique(always)))
        stop("Duplicate elements in always")
    if (!is.null(always) && length(always) > p - 2)
        stop("always has too many elements, less than two regressors left")
    if (!is.null(always) && any(c("car","genizi") %in% type)) 
         stop("Metrics genizi and car do not work with always.") 

    if (is.numeric(always) && any(!always %in% 2:(p+1)))
        stop("Numbers in always must be between 2 and p+1, corresponding to the position of regressors in cov(Y,X1,...,Xp).")

    if (is.null(x)) {
       if (!is.data.frame(y) && !is.matrix(y)) 
           stop("If x is NULL, then object must be a data frame or a matrix.")
       if (is.matrix(y) && !is.numeric(y))
           stop("Matrix object must be numeric.")
       names <- colnames(y)
       # not a covariance matrix
       if (!(ncol(y) == nrow(y))) {
             n <- nrow(y)     
             nomiss <- complete.cases(y)
             y <- y[nomiss,] 
             nobs <- nrow(y)
             if (!nobs==n) 
                 warning(paste((n-nobs), "observations deleted due to missing"))
             if (!(nobs > ncol(y))) 
                stop("Too few complete observations for estimating this model")
       }
    }
    
    if (!is.null(x)) {
       if ( (is.matrix(y) && ncol(y)>1) || !is.numeric(y) )
           stop(paste("object must be a numeric vector or one-column matrix,",
              "\n", "if x is also given.",sep=""))
       if (!(length(y) == nrow(x))) 
            stop("number of rows in object and x MUST be identical")
       if (is.matrix(x) && !is.numeric(x))
           stop("Matrix x must be numeric.")
        n <- length(y)
        nomiss <- complete.cases(y,x)
        nobs <- sum(nomiss)
        if (!nobs==n) 
              warning(paste((n-nobs), "observations deleted due to missing"))
        if (!(nobs > ncol(x)+1)) 
              stop("Too few complete observations for estimating this model")
         ## ynam <- deparse(substitute(object)) happened in calling program
         if (is.matrix(y)  & !is.null(colnames(y))) ynam <- colnames(y)
         if (is.null(colnames(x))) 
            colnames(x) <- paste("X", 1:ncol(x), sep = "")
         names <- c(ynam, colnames(x))
         y <- y[nomiss]
         x <- x[nomiss,]
        }
            if (is.character(always) && !all(always %in% names[2:(p+1)]))
                stop("Names in always must come from the names of the regressors.")
            if (is.character(always)) always <- which(names %in% always)

       if ((is.null(x) && !(ncol(y) == nrow(y))) || !is.null(x)) {
             ## not a covariance matrix, so that n is available and weights make sense
             wt <- rep(1/n, n)
             if (!is.null(weights)) {
                if (!(length(weights)==n && is.numeric(weights)) ) 
                stop("weights must be a numeric vector with an entry for each observation.")
                if (min(weights)<0) 
                stop("All weights must be non-negative.")
                if (!is.null(design)) stop("You can EITHER specify weights OR a design.")
                wt <- weights
                }
             if (!is.null(design)) {
                if (!"survey.design" %in% class(design))
                stop("design must be an object of class survey.design.")
                if (!nrow(design$variables)==n) 
                    stop (paste("The number of rows of the dataframe in the design must be identical", "\n",
                     "to the number of data rows given to the function directly (after omitting missing values and applying subset requests).",
                     sep=""))
                wt <- weights(design)
                if (!length(wt) == n) stop("Length of object and length of design do not match.")
                }
             wt <- wt[nomiss]  }
    
    #### error checks for groups-related things
    #### groups is transferred into a list of vectors
    if (!is.null(groups)) 
    {
      if (any(c("betasq","pratt","car","genizi") %in% type)) 
         stop("Metrics betasq, pratt, genizi and car do not work with groups.") 
      if (!is.list(groups)) {
        if ((is.numeric(groups) || is.integer(groups)) && !all(groups %in% 2:(p+1)))
            stop(paste("Numbers in groups must refer to columns 2 to ", p+1, " in Y,X1,...,Xp.", sep=""))
        if (length(groups) <= 1) 
            stop("groups must list groups of more than one regressor.")
        if (!(is.character(groups) || is.numeric(groups))) stop("groups must be numeric or character")
        groups<-list(groups)
        }
        ## make list into numeric list
        groups <- lapply(groups, function(obj){        
                  if (is.numeric(obj) || is.integer(obj)) {
                     if (!all(obj %in% 2:(p+1)))
                     stop(paste("Numbers in elements of groups must refer to columns 2 to ", p+1, " in Y,X1,...,Xp.", sep=""))
                     obj <- names[obj]
                     }
                     
                     else if (!is.character(obj)) stop("groups must be numeric or character") 
                     if (is.character(obj)) { 
                          if (!all(obj %in% names[2:(p+1)])) {
                          if (hasfacs) {
                            if (is.null(x)) {
                              faclevs <- lapply(y,function(obj) levels(obj) )
                              names(faclevs)<-colnames(y)
                            }
                            if (!is.null(x)) {
                              faclevs <- lapply(x,function(obj) levels(obj))
                              names(faclevs)<-colnames(x)
                            }
                            faclevs <- faclevs[sapply(faclevs,function(obj) length(obj)>1)]
                            for (i in 1:length(faclevs)) {
                                 if (!is.null(faclevs[[i]]))
                                  faclevs[[i]]<- paste(names(faclevs)[i],faclevs[[i]],sep="")
                                 }
                             kommtvor <- sapply(faclevs,function(fobj) {
                                   obj[!obj %in% names[2:(p+1)]] %in% fobj
                                   })
                             if (all(rowSums(kommtvor))) {
                                  ## map factorlevels to factor names again
                                  for (i in 1:length(faclevs)) {
                                     obj[which(obj %in% faclevs[[i]])] <- names(faclevs)[i]
                                     obj <- unique(obj)
                                     }
                          }
                              else stop("Names in groups are not all from regressor columns." )
                          }
                          else stop("names in groups are not all from regressor columns." )
                          }}
                     obj
                     })
        if (!is.null(groupnames)) 
             groupnames <- groupnames[which(sapply(groups,function(obj) length(obj)>1))]         
        groups <- groups[which(sapply(groups,function(obj) length(obj)>1))]         
        if (length(groups)==0) {groups <- NULL; groupnames <- NULL}
                     }
        ogroups <- groups
        ogroupnames <- groupnames
      ## now groups is a character list, 
      ## and mismatch between factor levels in groups and data still containing unsplit factors 
      ## can be handled 
    if (!is.null(ogroups)) 
    {
        ogroups <- lapply(ogroups, function(obj){        
              if(is.character(obj)) obj <- which(names %in% obj)
              obj } )
        # now ogroups is a numeric list
             if (!length(list2vec(ogroups))==length(unique(list2vec(ogroups))))
              stop("Overlapping groups are not permitted!")
             if (any(always %in% list2vec(ogroups))) 
             stop("groups must not refer to regressors that also occur in always.")

    
       if (!is.null(ogroupnames)) {
          if (!length(ogroups)==length(ogroupnames)) 
           stop(paste("groupnames must have one entry for each group.", "\n", 
               "There are", length(ogroups), "groups and", length(ogroupnames), "group names.")) }
    }

    if ("subset" %in% names(list(...)))
        warning("subset selection has been ignored.")    

    ## following block originally by ML
    # convert factors to binary coded data
    g_facnames <- character(0) 
    if (is.null(x) && is.data.frame(y) && any(as.logical(lapply(y,is.factor)))) {
         ## groups is list or NULL, make it character list
        if (!is.null(groups)) 
           groups <- lapply(groups, function(obj){
               if (is.numeric(obj)) colnames(y)[obj] else obj
               })
        ## make ogroups contain the correct column numbers
        if (!is.null(groups)) 
           ogroups <- lapply(groups, function(obj){
               which(colnames(y) %in% obj)
               })
        if (is.numeric(always)) always <- colnames(y[always])
        # Modelmatrix erzeugen
        y.modelm <- model.matrix(lm(y))
        y <- cbind(y[,1], as.data.frame(y.modelm)[,-1])
        colnames(y)[1] <- colnames(object)[1]
        # Gruppenvariablen hinzufugen, wenn notig
        g_names <- attr(y.modelm, "dimnames")[[2]]
        g_assign <- attr(y.modelm, "assign")
        g_facpos <- as.data.frame(table(g_assign))
        for (k in 2:length(g_facpos[,1])) {
          ## k=1 is Intercept
          fac.name <- colnames(object)[g_facpos[k,1]]
          fac.varn <- g_names[g_assign == g_facpos[k,1]]
          if (g_facpos[k,2]==1) {
              colnames(y)[g_assign == g_facpos[k,1]]<-fac.name
              ## make groups findable in case of two levels only
              if (!is.null(groups)) groups <- lapply(groups, function(obj){
                    if (fac.varn %in% obj) obj[obj==fac.varn]<-fac.name
                    obj
                 })
               next
              }
          if (any(always %in% fac.name)) {
               always <- c(always[always != fac.name], fac.varn)
                 ## darf nachher nicht mehr als Faktor auftauchen
               next
               }
          ## take care of the case with groups including a factor with more than two levels
           if (!is.null(groups)) {
                   schonda <- FALSE
                   for (ll in 1:length(groups)) {
                        if (fac.name %in% groups[[ll]]) {
                           schonda <- TRUE
                           groups[[ll]] <- c(groups[[ll]][groups[[ll]] != fac.name], fac.varn)
                           ogroups[[ll]] <- unique(c(setdiff(ogroups[[ll]],which(fac.name %in% names)),
                                     which(g_names %in% groups[[ll]])))
                        }
                      }
                 if (schonda) next
              }
             g_facnames <- c(g_facnames,fac.name)
             if (is.null(groups)) groups <- list(fac.varn)
             else groups[[length(groups)+1]] <- fac.varn
              }
            if (!is.null(groups)) groups <- lapply(groups, function(obj){which(colnames(y) %in% obj)})
            if (!is.null(always)) always <- which(colnames(y) %in% always)
            }

    if (!is.null(x) && is.data.frame(x) && any(as.logical(lapply(x,is.factor)))) {
        if (!is.null(groups)) 
           groups <- lapply(groups, function(obj){
               if (is.numeric(obj)) colnames(x)[obj-1] else obj
               })
        if (is.numeric(always)) always <- colnames(x[always-1])
        # Modelmatrix erzeugen
        x.model <- model.matrix(lm(y ~ ., data = x))
        x <- as.data.frame(x.model[,-1])
        ## make ogroups contain the correct column numbers
        if (!is.null(groups)) 
           ogroups <- lapply(groups, function(obj){
               which(colnames(xcall) %in% obj)+1
               })
        # weitere Variablen
        g_names <- attr(x.model, "dimnames")[[2]]
        g_assign <- attr(x.model, "assign")
        g_facpos <- as.data.frame(table(g_assign))
        for (k in 2:length(g_facpos[,1])) {
            fac.name <- colnames(xcall)[as.numeric(as.character(g_facpos[k,1]))]
            fac.varn <- g_names[g_assign == g_facpos[k,1]]
            if (g_facpos[k,2]==1) {
                colnames(x)[(g_assign == g_facpos[k,1])[-1]]<-fac.name
              ## make groups findable in case of two levels only
              if (!is.null(groups)) groups <- lapply(groups, function(obj){
                    if (fac.varn %in% obj) obj[obj==fac.varn]<-fac.name
                    obj
                 })
                next
                }
            if (any(always %in% fac.name)) {
               always <- c(always[always != fac.name], fac.varn)
                 ## darf nachher nicht mehr als Faktor auftauchen
               next
               }
          ## take care of the case with groups including a factor with more than two levels
           if (!is.null(ogroups)) {
                   schonda <- FALSE
                   for (ll in 1:length(groups)) {
                        if (fac.name %in% groups[[ll]]) {
                           schonda <- TRUE
                           groups[[ll]] <- c(groups[[ll]][groups[[ll]] != fac.name], fac.varn)
                           ogroups[[ll]] <- unique(c(setdiff(ogroups[[ll]],which(fac.name %in% names)),
                                     which(g_names %in% groups[[ll]])))
                        }
                      }
              if (schonda) next
              }
            g_facnames <- c(g_facnames,fac.name)
            if (is.null(groups)) groups <- list(fac.varn)
            else groups[[length(groups)+1]] <- fac.varn
                }
            if (!is.null(groups)) groups <- lapply(groups, function(obj){which(colnames(x) %in% obj)+1})
            if (!is.null(always)) always <- which(colnames(x) %in% always)+1
        # now groups is a list of numerics
        # groups could be NULL if only two-level-factors
            }

    
    ## previous block by ML    

    #error control and initialization continued

    if (is.null(x)) {
       # not a covariance matrix
       if (!(ncol(y) == nrow(y))) {
             if (!(nobs > ncol(y))) 
                stop("Too few complete observations for estimating this model")
             ## change UG 1.3: die cov-Funktion gewichtet
             covg <- cov.wt(y, wt=wt)$cov
       }
       else covg <- y
    }

    if (!is.null(x)) {
        if (!(nobs > ncol(x)+1)) 
              stop("Too few complete observations for estimating this model")
        ## change UG 1.3: die cov-Funktion gewichtet
        covg <- cov.wt(cbind(y,x), wt=wt)$cov
        if (is.null(colnames(covg))) colnames(covg) <- c("y", paste("X", 1:ncol(x), sep = ""))
        if (colnames(covg)[1]=="y") colnames(covg)[1] <- ynam       
        if (is.null(colnames(x))) colnames(covg)[2:ncol(covg)] <- paste("X", 1:ncol(x), sep = "")
        }

    # now, after treatment of the data frame, all covg need to be square matrices
    if (!is.matrix(covg)) stop(paste("If object is square, ", "\n", 
          "it must be a covariance matrix (of type matrix)."))

    # check covariance matrix properties
    hilf <- eigen(covg, only.values = T)$values
    if (is.complex(hilf) || min(hilf) <= 0)
        stop(paste("covg must be", "\n", 
          "a positive definite covariance matrix", "\n",
          "or a data matrix / data frame with linearly independent columns."))

    # no of regressors
    p <- ncol(covg) - 1
    g <- p
    names <- colnames(covg)
    if (is.null(names)) 
        names <- c("y", paste("X", 1:p, sep = ""))

    ## make sure that all group lists are numeric
    if (!is.null(groups)) groups <- lapply(groups, function(obj){
        if (is.character(obj)) which(names %in% obj) else obj})

    #initialise output matrices
    wahr <- matrix(0, 1, p)
    #vector for setdiff
    alle <- 1:(p + 1)
    var.y <- covg[1,1]
    covall <- covg
    andere <- alle

    # adjust out variables that are supposed to always stay in the model
    if (!is.null(always))
         {if (is.character(always)) always <- which(names %in% always)
             covall <- covg
             andere <- setdiff(alle, always)
             if (length(always)==1)
             covg <- covg[andere, andere] - covg[andere,always]%o%covg[always,andere]/covg[always,always]
             if (length(always)>1)
             covg <- covg[andere, andere] - covg[andere,always]%*%solve(covg[always,always],covg[always,andere])
                 var.y.distrib <- covg[1,1]  
             alwaysnam <- names[always]
             R2.always <- (var.y - covg[1,1])/var.y
         }
    g <- p - length(always)
    
    ## redo here because of factor-based groups
    if (!is.null(groups)) {
       if (any(c("betasq","pratt","genizi","car") %in% type)) 
          stop("Metrics betasq, pratt, genizi and car do not work with groups or factors.")} 

    if (!is.null(ogroups)) {
    if (!is.null(ogroupnames)) 
       groupnames <- c(ogroupnames,g_facnames) 
       else groupnames <- c(paste("G", 1:length(ogroups), sep=""), g_facnames)

       ## groupdocu will support printing the meaning of groups
       groupdocu <- list(groupnames, lapply(groups, function(obj){names[obj]}))
                 ### use correct columns of x-matrix
       groupdocu[[1]] <- append(groupdocu[[1]],  as.list(names[setdiff(2:(p+1), c(list2vec(groups),always))]))
       groupdocu[[2]] <- append(groupdocu[[2]],  as.list(names[setdiff(2:(p+1), c(list2vec(groups),always))]))

       ## groups will be used for picking appropriate elements from covariance matrix
       ## g ist number of groups
       groupnames <- c(groupnames, names[setdiff(alle, c(1,list2vec(groups),always))]) 
       groups <- append(groups, as.list(setdiff(alle, c(1,list2vec(groups),always))))

       g <- length(groups)
       }       
    else {
        if (length(g_facnames)>0) {
            groupnames <- g_facnames
            groupdocu <- list(groupnames, lapply(groups, function(obj){names[obj]}))
                 ### use correct columns of x-matrix
            groupdocu[[1]] <- append(groupdocu[[1]],  as.list(names[setdiff(2:(p+1), c(list2vec(groups),always))]))
            groupdocu[[2]] <- append(groupdocu[[2]],  as.list(names[setdiff(2:(p+1), c(list2vec(groups),always))]))

             ## groups will be used for picking appropriate elements from covariance matrix
             ## g ist number of groups
             groupnames <- c(groupnames, names[setdiff(alle, c(1,list2vec(groups),always))]) 
             groups <- append(groups, as.list(setdiff(alle, c(1,list2vec(groups),always))))

            g <- length(groups)   
        }
    }


##    ## always-Manipulationen, die nach groups-Abschnitt erfolgen sollen
    if (!is.null(always)) {
        names <- names[andere]
        p <- p - length(always)
        alle <- 1:(p+1)
    if (!is.null(groups)) 
       groups <- lapply(groups, function(obj){
           obj - rowSums(matrix(obj,length(obj),length(always),byrow=F) >  
                  matrix(always,length(obj),length(always),byrow=T))
         } )
##    ## always==obj impossible because of error checking
    }



    # start by calculating all conditional variances of y
    # construct all subsets with function nchoosek (included in relaimpo, taken from package vsn),
    # calculate conditional variances as top left corner of appropriate conditional covariance matrix
    # and write them to vectors in same position as with index matrices

    #initialise lists of index matrices and variance rows
    
    #first element of list unconditional, list initialized in full length
    indices <- rep(list(0), g + 1)
    variances <- rep(list(covg[1, 1]), g + 1)
    betas <- matrix(NA,p,g)
    betas[,g] <- solve(covg[-1,-1],covg[-1,1])

    #conditioning on all variables, i.e. var=s^2
    indices[[g + 1]] <- matrix(1:g, g, 1)
    variances[[g + 1]] <- covg[1:1] - covg[1, 2:(p + 1)] %*% 
        solve(covg[2:(p + 1), 2:(p + 1)], covg[2:(p + 1), 1])

    hilf <- varicalc(type, alle, covg, p, indices, variances, betas, g, groups, ngroups, WW)

    indices <- hilf$indices
    variances <- hilf$variances
    betas <- hilf$betas

    #output R-squared in order to show the total that is subdivided
    if (!is.null(always)) ausgabe <- new("relimplm", var.y = var.y, 
        R2 = as.numeric(1 - variances[[g + 1]]/var.y), 
        R2.decomp = as.numeric(variances[[1]] - variances[[g + 1]])/var.y)
    if (is.null(always)) ausgabe <- new("relimplm", var.y = as.numeric(variances[[1]]), 
        R2 = as.numeric(1 - variances[[g + 1]]/variances[[1]]), 
        R2.decomp = as.numeric(1 - variances[[g + 1]]/variances[[1]]))

    if (!is.null(groups)) names <- c(names[1],as.character(groupnames))
    if ("lmg" %in% type) 
        {
        if (is.null(WW))
        ausgabe <- lmgcalc(ausgabe, g, indices, variances, rank, 
            diff, rela, var.y)
        else
        ausgabe <- lmgcalcWW(ausgabe, g, indices, variances, rank, 
            diff, rela, var.y, WW, ngroups=ngroups)
        names(ausgabe@lmg)<-names[2:(g+1)]
        if (rank) names(ausgabe@lmg.rank)<-names[2:(g+1)]
        }
    if ("pmvd" %in% type) 
        {
        ausgabe <- pmvdcalc(ausgabe, g, indices, variances, rank, 
            diff, rela)
        names(ausgabe@pmvd)<-names[2:(g+1)]
        if (rank) names(ausgabe@pmvd.rank)<-names[2:(g+1)]
        }
   if ("last" %in% type) 
        {
        ausgabe <- lastcalc(ausgabe, g, variances, rank, diff, 
            rela, var.y)
        names(ausgabe@last)<-names[2:(g+1)] 
        if (rank) names(ausgabe@last.rank)<-names[2:(g+1)]
        }
   if ("first" %in% type) 
        {
        ausgabe <- firstcalc(ausgabe, g, variances, rank, diff, 
            rela, var.y)
        names(ausgabe@first)<-names[2:(g+1)]
        if (rank) names(ausgabe@first.rank)<-names[2:(g+1)]
        }
    if ("betasq" %in% type) 
        {
        ausgabe <- betasqcalc(ausgabe, covg, g, variances, rank, 
            diff, rela, var.y)
        names(ausgabe@betasq)<-names[2:(g+1)]
        if (rank) names(ausgabe@betasq.rank)<-names[2:(g+1)]
        }
    if ("pratt" %in% type) 
        {
        ausgabe <- prattcalc(ausgabe, covg, g, rank, diff, rela, var.y)
        names(ausgabe@pratt)<-names[2:(g+1)]
        if (rank) names(ausgabe@pratt.rank)<-names[2:(g+1)]
        }

    if ("genizi" %in% type) 
        {
        ausgabe <- genizicalc(ausgabe, covg, g, rank, diff, rela, var.y)
        names(ausgabe@genizi)<-names[2:(g+1)]
        if (rank) names(ausgabe@genizi.rank)<-names[2:(g+1)]
        }

    if ("car" %in% type) 
        {
        ausgabe <- carcalc(ausgabe, covg, g, rank, diff, rela, var.y)
        names(ausgabe@car)<-names[2:(g+1)]
        if (rank) names(ausgabe@car.rank)<-names[2:(g+1)]
        }

    #ausgabe contains (in this order) var.y, R2, lmg, rank.lmg, diff.lmg, 
    #                        pmvd, rank.pmvd, diff.pmvd,  (non-US version only)
    #                        last, rank.last, diff.last, first, rank.first, diff.first,
    #                                 betasq, rank.betasq, diff.betasq, pratt, rank.pratt, diff.pratt,
    #                                 genizi, genizi.rank, genizi.diff, car, rank.car, diff.car    
    # as far as requested by the call
    # default: R2, lmg, rank.lmg
    # in addition, some logicals and names are included

    slot(ausgabe, "rela") <- rela
    slot(ausgabe, "namen") <- names
    if (!is.null(nobs)) slot(ausgabe, "nobs") <- nobs
     slot(ausgabe, "always") <- always
    if (!is.null(always)) slot(ausgabe, "alwaysnam") <- alwaysnam
    slot(ausgabe, "type") <- alltype[which(alltype %in% type)]
           # this cryptic approach ensures the correct order of types
           # and makes it possible for type to be a list without generating an error
    ## change UG 1.3: incorporate call
    slot(ausgabe, "call") <- sys.call(which = 1)
    if (is.null(groups) & !is.null(betas)) {
       rownames(betas) <- names[-1]
       colnames(betas) <- paste(1:g,c("X",rep("Xs",g-1)),sep="")
    }
    else {
       if (!is.null(betas)) {
         rownames(betas) <- colnames(covall)[-c(1,always)]
         colnames(betas) <- paste(1:g,c("group",rep("groups",g-1)),sep="")
       }
    }
    slot(ausgabe, "ave.coeffs") <- betas
    if (!is.null(groups)) slot(ausgabe, "groupdocu") <- groupdocu
    if (test)
       {
       ausgabe <- new("relimplmtest", ausgabe)
       if (is.null(x)) {
           if (!nrow(y)==ncol(y)) {
           ausgabe@daten <- as.matrix(y)
           ausgabe@wt <- wt
           }}
       else {
           ausgabe@daten <- as.matrix(cbind(y,x))
           colnames(ausgabe@daten)[1]<-ynam
           ausgabe@wt <- wt
           }
       }
    return(ausgabe)
}


