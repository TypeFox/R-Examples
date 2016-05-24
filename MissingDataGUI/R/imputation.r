##' Impute the missing data with the method selected under the
##' condition.
##'
##' This function provides eight methods for imputation with
##' categorical varaibles as conditions.
##'
##' The imputation methods: This list displays all the imputation
##' methods. Users can make one selection. (1) 'Below 10%' means
##' NA's of one variable will be replaced by the value which
##' equals to the minimum of the variable minus 10% of the range.
##' For categorical variables, NA's are treated as a new category.
##' Under this status the selected conditioning variables are
##' ignored. If the data are already imputed, then this item
##' will show the imputed result. (2) 'Simple' will create three
##' tabs: Median, Mean, and Random Value. 'Median' means NA's will be
##' replaced by the median of this variable (omit NA's).
##' 'Mean' means NA's will be replaced by the mean of the
##' variable (omit NA's). The median does not apply to the nominal
##' variable, neither does the mean to the categorical variable.
##' In these cases the mode (omit NA's) is provided.
##' 'Random Value' means NA's will be replaced by any values
##' of this variable (omit NA's) which are randomly selected.
##' (3) 'Neighbor' contains two methods: 'Average Neighbor' and
##' 'Random Neighbor'. 'Average Neighbor' will replace the
##' NA's by the mean of the nearest neighbors. 'Random Neighbor'
##' substitutes the missing for a random sample of the k nearest
##' neighbors. The number of neighbors is default to 5, and
##' can be changed by argument knn. The Neighbor methods
##' require at lease one case to be complete, at least two
##' variables to be selected, and no factor/character variables.
##' The ordered factors are treated as integers. The method will
##' return the overall mean or a global random sample value
##' if the observation only contains NA's.
##' (4) 'MI:areg' uses function \code{\link[Hmisc]{aregImpute}}
##' from package \pkg{Hmisc}. It requires at lease one case to
##' be complete, and at least two variables to be selected.
##' (5) 'MI:norm' uses function \code{\link[norm]{imp.norm}}
##' from package \pkg{norm}. It requires all selected variables
##' to be numeric(at least integer), and at least two variables
##' to be selected. Sometimes it cannot converge, then the
##' programme will leave NA's without imputation.
##' (6) 'MI:mice' uses the \pkg{mice} package. The methods of
##' the variables containing NA's must be attached with argument
##' method. If not, then default methods are used.
##' (7) 'MI:mi' employes the \pkg{mi} package.
##' @param origdata A data frame whose missing values need to be
##' imputed. This data frame should be selected from the missing
##' data GUI.
##' @param method The imputation method selected from the missing
##' data GUI. Must be one of 'Below 10%','Simple','Neighbor',
##' 'MI:areg','MI:norm','MI:mice','MI:mi'. If method='MI:mice',
##' then the methods of the variables containing NA's must be attached
##' with argument method. If not, then default methods are used.
##' @param vartype A vector of the classes of origdata. The length is
##' the same as the number of columns of origdata. The value should be
##' from "integer", "numeric", "logical", "character", "factor", and "ordered".
##' @param missingpct A vector of the percentage of missings of the
##' variables in origdata. The length is the same as the number of
##' columns of origdata. The values should be between 0 and 1.
##' @param condition A vector of categorical variables. The dataset
##' will be partitioned based on those variables, and then the
##' imputation is implemented in each group. There are no missing
##' values in those variables. If it is null, then there is no
##' division. The imputation is based on the whole dataset.
##' @param knn number of the neighbors.
##' @param mi.n number of the imputation sets for multiple imputation
##' @param mi.seed random number seed for multiple imputation
##' @param row_var A column name (character) that defines the ID of rows.
##' @return The imputed data frame with the last column being the row
##' number from the original dataset. During the procedure of the
##' function, rows may be exchanged, thus a column of row number could
##' keep track of the original row number and then help to find the
##' shadow matrix.
##' @export
##' @author Xiaoyue Cheng <\email{xycheng@@unomaha.edu}>
imputation = function(origdata, method, vartype=NULL, missingpct=NULL, condition=NULL, knn=5, mi.n=3, mi.seed=1234567, row_var=NULL){
    if (is.null(origdata)) return(NULL)
    if (is.null(vartype)) vartype=unname(sapply(origdata,function(x)class(x)[1]))
    if (is.null(missingpct)) missingpct=unname(sapply(origdata,function(x)mean(is.na(x))))
    row_NO = if (is.null(row_var)) {1:nrow(origdata)} else {origdata[,row_var]}
    origdata = origdata[,setdiff(colnames(origdata),row_var),drop=FALSE]
    n=ncol(origdata)
    if (n==1 && substr(method,1,3) == 'MI:'){
        gmessage('You only selected one variable. Cannot apply the multiple imputation.', icon = "error")
        return(NULL)
    }
    if (!is.null(condition)) {
      if (sum(is.na(condition))) {
        condition = sapply(condition,as.character)
        condition[is.na(condition)] = "NaN"
        condition = data.frame(condition)
      }
        dat = by(data.frame(origdata,row_id=1:nrow(origdata)),condition,imputation,method=method,vartype=vartype,missingpct=NULL,condition=NULL,knn=knn,mi.n=mi.n,mi.seed=mi.seed,row_var='row_id',simplify=FALSE)
        dat = dat[!sapply(dat,is.null)]
        k1=length(dat[[1]])
        k2=length(dat)
        res=list()
        for (j in 1:k1) res[[j]]=lapply(dat,function(x) x[[j]])
        names(res)=names(dat[[1]])
        res=lapply(res, function(x){
            tmp=NULL
            for (i in 1:k2) {
                tmp=rbind(tmp,x[[i]])
            }
            return(tmp)
        })
        res = lapply(res,function(x)x[order(x$row_number),c(colnames(origdata),'row_number')])
        return(res)
    }

    dat = list(d1=origdata)
    origshadow = is.na(origdata[,1:n,drop=FALSE])
    if (method == 'Below 10%') {
        for (i in 1:n) {
            if (vartype[i] %in% c('integer','numeric','logical') &
                    as.numeric(as.character(missingpct))[i]>0) {
                dat$d1[origshadow[,i],i] = min(origdata[,i], na.rm=TRUE)*1.1-
                    max(origdata[,i],na.rm=TRUE)/10
            }
            if (vartype[i] == "character") {
                dat$d1[origshadow[,i],i] = 'NAN'
            }
            if (vartype[i] %in% c("factor","ordered") &&
                    as.numeric(as.character(missingpct))[i]>0) {
                dat$d1[,i]=factor(dat$d1[,i],levels = c('NAN',levels(factor(dat$d1[,i]))))
                dat$d1[origshadow[,i],i] = 'NAN'
            }
        }
        names(dat)="Below 10%"
    }
    else if (method == 'Simple') {
      set.seed(mi.seed)
        dat$d3=dat$d2=dat$d1
        for (i in 1:n) {
          if (sum(is.na(dat$d1[,i]))){
            if (vartype[i] %in% c('integer','numeric')) {
              dat$d1[is.na(dat$d1[,i]),i] = median(dat$d1[,i], na.rm=TRUE)
              dat$d2[is.na(dat$d2[,i]),i] = mean(dat$d2[,i], na.rm=TRUE)
            } else if (vartype[i] %in% c('logical','factor','character')){
              biggroup = names(sort(table(na.omit(dat$d2[,i])),decreasing=TRUE))[1]
              dat$d1[origshadow[,i],i] = biggroup
              dat$d2[origshadow[,i],i] = biggroup
            } else {
              grouplevel = levels(dat$d1[,i])
              dat$d1[origshadow[,i],i] = grouplevel[round(median(as.integer(dat$d1[,i]), na.rm=TRUE))]
              dat$d2[origshadow[,i],i] = names(sort(table(na.omit(dat$d2[,i])),decreasing=TRUE))[1]
            }
            fill=sample(dat$d3[!origshadow[,i],i], sum(origshadow[,i]), replace = TRUE)
            dat$d3[origshadow[,i],i] = fill
          }
        }
        names(dat)=c('Median','Mean','Random Value')
    }
    else if (method == 'Neighbor') {
        set.seed(mi.seed)
        if (n==1 || sum(complete.cases(origdata))==0 || any(c("factor","character") %in% vartype)) {
            if (n==1) warning_message='You only selected one variable. Cannot apply the neighbor imputation. Only the random value imputation is given.'
            if (sum(complete.cases(origdata))==0) warning_message="All the observations have missing values. Cannot find the neighbors."
            if (any(c("factor","character") %in% vartype)) warning_message='Cannot impute by the neighbor with one or more factor or character variables.'
            gmessage(warning_message, icon='warning')
            return(NULL)
        } else {
            dat$d2=dat$d1
            CmpltDat = dat$d1[complete.cases(dat$d1),]
            orderedfactor = which(vartype=="ordered")
            if (length(orderedfactor)) {
              for (j in orderedfactor) {
                CmpltDat[,j]=as.integer(CmpltDat[,j])
                dat$d1[,j]=as.integer(dat$d1[,j])
                dat$d2[,j]=as.integer(dat$d2[,j])
              }
            }
            for (i in which(!complete.cases(dat$d1))){
                usecol = which(!is.na(dat$d1[i,]))
                if (length(usecol)>0){
                    NNdat = CmpltDat
                    a = rbind(dat$d1[i,], NNdat)[,usecol,drop=FALSE]
                    b = apply(a,2,scale)
                    NNdat$distance = dist(b)[1:nrow(NNdat)]
                    kNNdat = NNdat[order(NNdat$distance,decreasing=FALSE),][1:min(knn,nrow(NNdat)),]
                    dat$d1[i,-usecol] = colMeans(kNNdat[,1:n][,-usecol,drop=FALSE])
                    dat$d2[i,-usecol] = kNNdat[,1:n][sample(1:nrow(kNNdat),1),-usecol]
                } else {
                    dat$d1[i,] = sapply(origdata, mean, na.rm=TRUE)
                    dat$d2[i,] = sapply(origdata, function(x)sample(na.omit(x),1))
                }
            }
            if (length(orderedfactor)) {
              for (j in orderedfactor){
                dat$d1[,j] = round(dat$d1[,j])
                l = levels(origdata[,j])
                dat$d1[,j] = ordered(l[dat$d1[,j]],levels=l)
                dat$d2[,j] = ordered(l[dat$d2[,j]],levels=l)
              }
            }
            names(dat)=c('Average Neighbor','Random Neighbor')
        }
    }
    else if (method == 'MI:norm') {
        if (any(c('factor','character','ordered') %in% vartype)) {
            gmessage("Not every variable is numeric. Cannot impute missing values under the multivariate normal model.", icon = "error")
            return(NULL)
        } else {
          library_gui('norm')
            s = norm::prelim.norm(as.matrix(origdata))
            thetahat = norm::em.norm(s)
            norm::rngseed(mi.seed)
            theta = norm::da.norm(s,thetahat,steps=10,showits=TRUE)
            for (i in 1:mi.n){
              dat[[i]] = norm::imp.norm(s,theta,as.matrix(origdata))
            }
            names(dat)=paste('norm',1:mi.n)
            if (any(sapply(dat,function(x){any(c(Inf,NaN) %in% x)}))) {
                gmessage("The algorithm doesn't converge. Return the original data with missing values.", icon = "warning")
                dat = list(original=origdata)
            }
        }
    }
    else if (method == 'MI:areg') {
        if (sum(complete.cases(origdata))==0) {
            gmessage('All the observations have missing values. Cannot impute by Hmisc::aregImpute.', icon = "warning")
            return(NULL)
        } else {
          library_gui('Hmisc')
            formula0 = as.formula(paste('~ ',paste(names(origdata),collapse=' + ')))
            set.seed(mi.seed)
            f = Hmisc::aregImpute(formula0, data=origdata, n.impute=mi.n)
            tmpres = f$imputed
            tmpres = tmpres[!sapply(tmpres,is.null)]
            for (i in 1:mi.n) {
              dat[[i]]=origdata
              for (j in 1:length(tmpres)){
                dat[[i]][rownames(tmpres[[j]]),names(tmpres)[j]]=tmpres[[j]][,i]
              }
            }
            names(dat)=paste('areg',1:mi.n)
        }
    }
    else if (method == 'MI:mice') {
      library_gui('mice')
      if (is.null(attr(method,'method'))) attr(method,'method')=vector("character", length = ncol(data))
        f = mice::mice(origdata, method=attr(method,'method'), m=mi.n, printFlag=FALSE, seed=mi.seed)
        tmpres = f$imp
        for (j in 1:mi.n){
          dat[[j]] = origdata
          for (i in which(!sapply(tmpres,is.null))) {
            dat[[j]][rownames(tmpres[[i]]),names(tmpres)[i]]=tmpres[[i]][,j]
          }
        }
        names(dat)=paste('mice',1:mi.n)
    }
    else if (method == 'MI:mi') {
      library_gui('mi')
        f = mi::mi(origdata, n.chains=mi.n, seed=mi.seed)
        dat = mi::complete(f, include_missing = FALSE)
        names(dat)=paste('mi',1:mi.n)
    }
    res = lapply(dat, function(x) data.frame(x, row_number=row_NO))
    return(res)
}
