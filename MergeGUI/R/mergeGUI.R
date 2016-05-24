##' Obtain the intersection of a list of vectors.
##' Function "intersect" in the base package can only intersect two
##' vectors. The function "intersect2" is designed to obtain the
##' intersection and the difference for more than two vectors. The
##' input should be a list whose elements are the vectors, and the
##' outputs include the intersection of all vectors and a list whose
##' elements are the input vectors substracting the intersection.
##' Besides, intersect2 allows the labels of the vectors. If a list of
##' labels is given in the input, then the outputs will also include a
##' matrix of labels which match the intersection for the vectors, and
##' a list of labels which match the left part of the vectors.
##'
##' @param vname A list of labels.
##' @param simplifiedname A list of vectors to make the intersection.
##' Each element in the list has the same length as the corresponding
##' element in vname. Default to be vname. If simplifiedname is not
##' vname, then it works as the real vectors to match, and vname is
##' like the labels of simplifiedname. If simplifiedname is the same
##' as vname, then the returned value simpleuniq=uniq.
##' @return The outputs are 'public', 'individual', 'uniq', and
##' 'simpleuniq'.  'public' is a vector of the intersection of
##' 'simplifiedname'.  'individual' is a matrix with the original
##' colnames matched to 'public' in all files.  'simpleuniq' is a list
##' of the left part of 'simplifiedname' if we pick 'public' out.
##' 'uniq' is a list of the left part of 'vname' if we pick
##' 'individual' out.
##' @author Xiaoyue Cheng <\email{xycheng@@iastate.edu}>
##' @exportPattern "^[^\\.]"
##' @export
##' @examples
##' a = list(x1=c("label11","label12"),
##'    x2=c("label21","label22","label23"),
##'	x3=c("label31","label32"))
##' b = list(x1=c(1,2),x2=c(3,1,2),x3=c(2,1))
##' intersect2(a,b)
##'
intersect2 = function(vname, simplifiedname=vname) {
    
    s = as.vector(simplifiedname[[1]])
    for (i in 2:length(simplifiedname)) {
        s = intersect(as.vector(simplifiedname[[i]]), s)
    }
    v1 = matrix(nrow = length(s), ncol = length(vname))
    v2 = vname
    v3 = simplifiedname
    if (length(s) > 0) {
        for (i in 1:length(vname)) {
            for (j in 1:length(s)) {
                if (s[j] %in% simplifiedname[[i]]) {
                    v1[j, i] = vname[[i]][which(simplifiedname[[i]] ==
                        s[j])[1]]
                    tmp = vname[[i]][-which(simplifiedname[[i]] ==
                        s[j])[1]]
                    v2[[i]] = intersect(v2[[i]], tmp)
                    v3[[i]] = v3[[i]][-which(v3[[i]] == s[j])[1]]
                }
            }
        }
    }
    return(list(public = s, individual = v1, uniq = v2, simpleuniq = v3))
}

##' Short the names from a template.
##' The merging GUI is designed to merge data from different
##' files. But sometimes the file names are too long to be displayed
##' in the GUI. Hence this function is used to short the basenames by
##' removing the same beginning letters of each name. Hence the output
##' is a character vector whose elements will not start with the same
##' letter.
##'
##' @param namevector A character vector.
##' @return A character vector which cuts the first several same
##' letters from the input.
##' @author Xiaoyue Cheng <\email{xycheng@@iastate.edu}>
##' @export
##' @examples
##' simplifynames(c("abc234efg.csv","abc234hfg.csv"))
##' simplifynames(c("12345","54321"))
##' simplifynames(c("aeiou","aerial"))
##'
simplifynames=function(namevector) {
    n=max(nchar(namevector))
    for (i in 1:n){
        if (!all(substr(namevector,i,i)==substr(namevector,i,i)[1])){
            newnamevec=substring(namevector,i)
            return(newnamevec)
        }
    }
    return(namevector)
}

##' Detect the classes of the variables.
##'
##' This function gives an initial guess of the classes of each variable in the merged data.
##'
##' @param nametable.class A matrix of the matched variable names. The
##' number of columns is equal to the number of files. Each row
##' represents a variable that is going to be merged. Any elements
##' except NA in nametable.class must be the variable names in
##' dataset.class.
##' @param dataset.class The dataset list. The length of the
##' list is equal to the number of files, and the order of the
##' list is the same as the order of columns in nametable.class.
##' @return A vector matching the rows of 'nametable.class'. The value
##' includes NA if any variable are only NA's.
##' @author Xiaoyue Cheng <\email{xycheng@@iastate.edu}>
##' @export
##' @examples
##' a=data.frame(aa=1:5, ab=LETTERS[6:2], ac=as.logical(c(0,1,0,NA,0)))
##' b=data.frame(b1=letters[12:14],b2=3:1)
##' dat=list(a,b)
##' name=matrix(c("ab","aa","ac","b1","b2",NA),ncol=2)
##' var.class(name,dat)
##'
var.class = function(nametable.class, dataset.class) {
    varclass = rep("NA", nrow(nametable.class))
    for (i in 1:nrow(nametable.class)) {
        notNAcolset = which(!is.na(nametable.class[i, ]))
        notNAcolclass = c()
        for (k in notNAcolset) {
            if (sum(!is.na(dataset.class[[k]][, nametable.class[i, k]])) > 0)
                notNAcolclass = c(notNAcolclass,class(dataset.class[[k]][,nametable.class[i, k]]))
        }
        if (length(notNAcolclass) > 0) {
            varclass[i] = if ('factor' %in% notNAcolclass) {"factor"} else {
                if ('character' %in% notNAcolclass) {'character'} else {
                    if ('numeric' %in% notNAcolclass) {'numeric'} else {'integer'}
                }
            }
        }
    }
    return(varclass)
}

##' Compute the misclassification rate for each variable.
##' When merging data from several datasets, it is meaningful to
##' detect whether the matched variables from different files have
##' different centers. The function computes the misclassification
##' rate variable by variable using classification tree (the rpart
##' package). It will firstly merge the dataset by the given
##' nametable.class, then use rpart for each variable to seperate the
##' data without any covariates and compute the misclassification
##' rate.
##'
##' @param nametable.class A matrix of the matched variable names. The
##' number of columns is equal to the number of files.
##' The column names are required.  Each row
##' represents a variable that is going to be merged. Any elements
##' except NA in nametable.class must be the variable names in
##' dataset.class.
##' @param dataset.class The dataset list. The length of the
##' list is equal to the number of files, and the order of the
##' list is the same as the order of columns in nametable.class.
##' @param name.class A character vector of variable names. The length
##' of the vector must be equal to the number of rows in
##' nametable.class. Since the variable names in nametable.class may
##' not be consistent, name.class is needed to name the variables.
##' @param varclass A character vector of variable classes. The length
##' of the vector must be equal to the number of rows in
##' nametable.class. All the classes should be in "numeric",
##' "integer", "factor", and "character". Default to be null, then
##' it will be determined by \code{\link{var.class}}.
##' @return A vector of the misclassification rate. The rate is
##' between 0 and 1, or equal to 9 if one of more groups only have
##' NA's.
##' @author Xiaoyue Cheng <\email{xycheng@@iastate.edu}>
##' @export
##' @examples
##' a=data.frame(aa=1:5, ab=LETTERS[6:2], ac=as.logical(c(0,1,0,NA,0)))
##' b=data.frame(b1=letters[12:14],b2=3:1)
##' dat=list(a,b)
##' name=matrix(c("ab","aa","ac","b1","b2",NA),ncol=2)
##' colnames(name)=c("a","b")
##' newname=c("letter","int","logic")
##' scale_rpart(name,dat,newname)
##'
scale_rpart = function(nametable.class, dataset.class, name.class,varclass=NULL) {
    if (is.null(varclass)) {
        varclass = var.class(nametable.class,dataset.class)
    }
    rows = unlist(lapply(dataset.class, nrow))
    selectedvariables = which(varclass %in% c('numeric','integer','logical','factor'))
    mergedata = matrix(nrow = sum(rows), ncol = length(selectedvariables) + 1)
    colnames(mergedata) = c("source", name.class[selectedvariables])
    for (i in 1:length(dataset.class)) {
        tmp = matrix(c(rep(colnames(nametable.class)[i], rows[i]),
                       rep(NA, rows[i] * length(selectedvariables))), nrow = rows[i])
        colnames(tmp) = c("source", nametable.class[selectedvariables, i])
        tmp[, na.omit(nametable.class[selectedvariables, i])] = as.matrix(dataset.class[[i]])[,
                                                                                              na.omit(nametable.class[selectedvariables, i])]
        mergedata[(cumsum(rows) - rows + 1)[i]:cumsum(rows)[i],
                  ] = tmp
    }
    mergedata=as.data.frame(mergedata)
    mergedata$source=factor(mergedata$source)
    num = which(varclass[selectedvariables] %in% c('integer','numeric'))
    fac = which(varclass[selectedvariables] %in% c('logical','factor'))
    if (length(num)!=0) {
        if (length(num)>1) {
            mergedata[,num+1] = sapply(mergedata[,num+1],function(avec){as.numeric(as.character(avec))},simplify = FALSE)
        } else {
            mergedata[,num+1] = as.numeric(as.character(mergedata[,num+1]))
        }
    }
    if (length(fac)!=0) {
        if (length(fac)>1) {
            mergedata[,fac+1]=sapply(mergedata[,fac+1],as.factor,simplify = FALSE)
        } else {
            mergedata[,fac+1] = as.factor(mergedata[,fac+1])
        }
    }
    
    res = rep(NA, nrow(nametable.class))
    group = mergedata$source
    if (length(levels(group))==2) {
        for (i in 2:ncol(mergedata)) {
            if (!(varclass[i-1] %in% c("factor","character") & length(unique(mergedata[,i]))>10) & sum(tapply(mergedata[,i],mergedata[,1],function(x){!all(is.na(x))}))>=2) {	
                fit_rpart = rpart::rpart(mergedata$source~mergedata[,i], control=c(maxdepth=1))
                tmperror = weighted.mean(residuals(fit_rpart), 1/table(group)[group[!is.na(mergedata[,i])]])
                res[name.class==colnames(mergedata)[i]] = round(tmperror,3)
            }
        }
    } else {
        for (i in 2:ncol(mergedata)) {
            if (!(varclass[i-1] %in% c("factor","character") & length(unique(mergedata[,i]))>10) & sum(tapply(mergedata[,i],mergedata[,1],function(x){!all(is.na(x))}))>=2) {
                tmperror2 = c()
                for (j in 1:length(levels(group))) {
                    tmpdat = data.frame(file=factor(group==levels(group)[j]),mergedata[,i])
                    if (any(tapply(tmpdat[,2],tmpdat[,1],function(x){all(is.na(x))}))){
                        tmperror2[j] = NA
                    } else {
                        fit_rpart = rpart::rpart(file~., data=tmpdat, control=c(maxdepth=1))
                        tmperror2[j] = weighted.mean(residuals(fit_rpart), 1/table(tmpdat[,1])[tmpdat[!is.na(mergedata[,i]),1]])
                    }                    
                }
                tmperror = min(tmperror2,na.rm=TRUE)
                res[name.class==colnames(mergedata)[i]] = round(tmperror,3)
            }
        }
    }
    return(as.character(res))
}

##' Compute the p-values of the Kolmogorov-Smirnov tests between
##' different sources for each variable.
##' This function is used to detect whether the matched variables from
##' different files have different distributions. For each variable,
##' it will compute the pairwise KS-test p-values among the sources,
##' then report the lowest p-value as the indice for this variable.
##'
##' @param nametable.class A matrix of the matched variable names. The
##' number of columns is equal to the number of files. Each row
##' represents a variable that is going to be merged. Any elements
##' except NA in nametable.class must be the variable names in
##' dataset.class.
##' @param dataset.class The dataset list. The length of the
##' list is equal to the number of files, and the order of the
##' list is the same as the order of columns in nametable.class.
##' @param name.class A character vector of variable names. The length
##' of the vector must be equal to the number of rows in
##' nametable.class. Since the variable names in nametable.class may
##' not be consistent, name.class is needed to name the variables.
##' @param varclass A character vector of variable classes. The length
##' of the vector must be equal to the number of rows in
##' nametable.class. All the classes should be in "numeric",
##' "integer", "factor", and "character". Default to be null, then
##' it will be determined by \code{\link{var.class}}.
##' @return A vector of p-values from the KS-test for each
##' variable.The p-values are between 0 and 1, or equal to 9 if one of
##' more groups only have NA's.
##' @author Xiaoyue Cheng <\email{xycheng@@iastate.edu}>
##' @export
##' @examples
##' a=data.frame(aa=1:5, ab=LETTERS[6:2], ac=as.logical(c(0,1,0,NA,0)))
##' b=data.frame(b1=letters[12:14],b2=3:1)
##' dat=list(a,b)
##' name=matrix(c("ab","aa","ac","b1","b2",NA),ncol=2)
##' colnames(name)=c("a","b")
##' newname=c("letter","int","logic")
##' scale_kstest(name,dat,newname)
##'
scale_kstest = function(nametable.class, dataset.class, name.class,varclass=NULL) {
    if (is.null(varclass)) {
        varclass = var.class(nametable.class,dataset.class)
    }
    rows = unlist(lapply(dataset.class, nrow))
    selectedvariables = which(varclass %in% c('numeric','integer'))
    mergedata = matrix(nrow = sum(rows), ncol = length(selectedvariables) + 1)
    colnames(mergedata) = c("source", name.class[selectedvariables])
    for (i in 1:length(dataset.class)) {
        tmp = matrix(c(rep(colnames(nametable.class)[i], rows[i]),
                       rep(NA, rows[i] * length(selectedvariables))), nrow = rows[i])
        colnames(tmp) = c("source", nametable.class[selectedvariables, i])
        tmp[, na.omit(nametable.class[selectedvariables, i])] = as.matrix(dataset.class[[i]])[,
                                                                                              na.omit(nametable.class[selectedvariables, i])]
        mergedata[(cumsum(rows) - rows + 1)[i]:cumsum(rows)[i],
                  ] = tmp
    }
    mergedata=as.data.frame(mergedata)
    mergedata$source=factor(mergedata$source)
    if (ncol(mergedata)>2) {
        mergedata[,2:ncol(mergedata)]=sapply(mergedata[,2:ncol(mergedata)], function(avec){as.numeric(as.character(avec))},simplify = FALSE)
    } else {
        mergedata[,2]=as.numeric(as.character(mergedata[,2]))
    }
    
    scaleclass = rep(NA, nrow(nametable.class))
    for (i in 2:ncol(mergedata)) {
        tmpdat = mergedata[,c(1,i)]
        sig = c()
        for (j in 1:length(levels(tmpdat[,1]))) {
            a = scale(tmpdat[tmpdat[,1]==levels(tmpdat[,1])[j],2], scale=FALSE)
            b = scale(tmpdat[tmpdat[,1]!=levels(tmpdat[,1])[j],2], scale=FALSE)
            if (all(is.na(a)) | all(is.na(b))) {
                sig[j] = 1
            } else {
                sig[j] = ks.test(a, b)$p.value
            }
        }
        scaleclass[selectedvariables][i-1] = min(sig, na.rm=TRUE)
    }
    return(as.character(round(scaleclass,3)))
}

##' Chi-square tests for the counts of missing and non-missing.
##' This function is used to detect whether the matched variables from
##' different files have different missing patterns. For each
##' variable, it will firstly count the missing and non-missing values
##' among the sources, and then form a contingency table. The p-value
##' of Chi-square test is computed from the contingency table and
##' finally reported for the variable.
##'
##' @param nametable.class A matrix of the matched variable names. The
##' number of columns is equal to the number of files. Each row
##' represents a variable that is going to be merged. Any elements
##' except NA in nametable.class must be the variable names in
##' dataset.class.
##' @param dataset.class The dataset list. The length of the
##' list is equal to the number of files, and the order of the
##' list is the same as the order of columns in nametable.class.
##' @param name.class A character vector of variable names. The length
##' of the vector must be equal to the number of rows in
##' nametable.class. Since the variable names in nametable.class may
##' not be consistent, name.class is needed to name the variables.
##' @return A vector of p-values from the Chisquare-test for the
##' missings of each variable. The p-values are between 0 and 1.
##' @author Xiaoyue Cheng <\email{xycheng@@iastate.edu}>
##' @export
##' @examples
##' a=data.frame(aa=1:5, ab=LETTERS[6:2], ac=as.logical(c(0,1,0,NA,0)))
##' b=data.frame(b1=letters[12:14],b2=3:1)
##' dat=list(a,b)
##' name=matrix(c("ab","aa","ac","b1","b2",NA),ncol=2)
##' colnames(name)=c("a","b")
##' newname=c("letter","int","logic")
##' scale_missing(name,dat,newname)
##'
scale_missing = function(nametable.class, dataset.class, name.class) {
    rows = unlist(lapply(dataset.class, nrow))
    mergedata = matrix(nrow = sum(rows), ncol = length(name.class) + 1)
    colnames(mergedata) = c("source", name.class)
    for (i in 1:length(dataset.class)) {
        tmp = matrix(c(rep(colnames(nametable.class)[i], rows[i]),rep(NA, rows[i] * length(name.class))), nrow = rows[i])
        colnames(tmp) = c("source", nametable.class[, i])
        tmp[, na.omit(nametable.class[, i])] = as.matrix(dataset.class[[i]])[,na.omit(nametable.class[, i])]
        mergedata[(cumsum(rows) - rows + 1)[i]:cumsum(rows)[i],] = tmp
    }
    mergedata=as.data.frame(mergedata)
    mergedata$source=factor(mergedata$source)
    
    chi2testout = c()
    missingclass = rep(NA, nrow(nametable.class))
    for (i in 2:ncol(mergedata)) {
        tmpdat = mergedata[,c(1,i)]
        missingcount = matrix(0, nrow=length(levels(tmpdat[,1])), ncol=2)
        for (j in 1:length(levels(tmpdat[,1]))) {
            missingcount[j,1] = sum(is.na(tmpdat[tmpdat[,1]==colnames(nametable.class)[j],2]))
            missingcount[j,2] = rows[j] - missingcount[j,1]
        }
        if (sum(missingcount[,1])==0 | sum(missingcount[,2])==0) {
            missingclass[i-1] = 1
        } else {
            missingclass[i-1] = chisq.test(missingcount)$p.value
            if (missingclass[i-1]=='NaN') missingclass[i-1]=NA
        }
    }
    return(as.character(round(missingclass,3)))
}

##' The Merging GUI.
##' This function will start with an starting interface, allowing 1)
##' selecting several data files; 2) doing the next command with more
##' than one files. There are two commands which could be selected:
##' match the variables, match the cases by the key variable. In the
##' matching-variable interface the user can 1) check the matching of
##' the variables among files and switch the variable names if they
##' are wrongly matched; 2) look at the numerical and graphical
##' summaries for the selected variables, or the dictionary for
##' selected factor varibles; 3) observe the misclassification rate,
##' KS-test p-values and Chi-square test p-values for each variable,
##' which helps to determine whether any transformation is needed for
##' the variable; (For each variable, the user may want to know
##' whether it could distinguish the sources correctly. So the
##' misclassification rate is calculated through the tree model.
##' KS-test is used to check whether any variable has different
##' distributions for different sources. And the Chi-square test is
##' useful when the user is interested in the pattern of missing
##' values among the sources.) 4) change the name or class for any
##' variable; 5) export the merged dataset and the summary for it. In
##' the matching-case interface the user can determine a primary key
##' for each data file and then merge the cases by the key.
##'
##' The merging GUI consists of four tabs. In the preferences tab,
##' user can choose whether the numerical p-values or the flag symbols
##' are displayed in the summary tab; whether the y-scales are free
##' for different data files when drawing the plots faceted by the
##' sources. In the checking tab, each data file has a list of
##' variable names, and the GUI will automatically arrange the order
##' of variable names to align the same names in one row. The user can
##' switch the order of the variables in one file's list. It is
##' possible to undo, redo, or reset the matching. In the summary tab,
##' there is a list of variable names on the left which corresponds to
##' the checking tab. The misclassification rate, KS-test p-values and
##' Chi-square test p-values for each variable may also be presented
##' with the variable names. On the top right there are three buttons:
##' Numeric summary, Graphical summary, and Dictionary. And the
##' results could be shown below the buttons. For the graphical
##' summary, histogram or barchart will be shown if a single variable
##' is selected. A scatterplot will be drawn if two numeric or two
##' factor varaibles are chosen. Side-by-side boxplots will be
##' presented when one numeric and one factor varaibles are
##' selected. A parallel coordinate plot is shown when all the
##' variables selected are numeric and there are more than two
##' variables. If more than two variables are chosen but the classes
##' of the variables are mixed, i.e. some are numeric, some are factor
##' or character, then histograms and barcharts will be drawn
##' individually. All the plots are facetted by the source. In the
##' export tab the user could select all or none variables by click
##' the buttons or choose several varaibles by Ctrl+Click. Then the
##' export button will export the merged data and the numeric
##' summaries of the selected variables into two csv files.
##'
##' @param ... names of the data frames to read
##' @param filenames A vector of csv file names with their full paths.
##' @param unit whether the test of the difference among the group centers is on or off
##' @param distn whether the test of the difference among the group distributions is on or off
##' @param miss whether the test of the difference among the group missing patterns is on or off
##' @return NULL
##' @author Xiaoyue Cheng <\email{xycheng@@iastate.edu}>
##' @import ggplot2 gWidgetsRGtk2 cairoDevice
##' @export
##' @examples
##' if (interactive()) {
##' MergeGUI()
##' 
##' csvnames=list.files(system.file("doc",package="MergeGUI"),pattern = "\\.csv$")
##' files=system.file("doc",csvnames,package="MergeGUI")
##' MergeGUI(filenames=files)
##' 
##' data(iris)
##' setosa=iris[iris$Species=="setosa",1:4]
##' versicolor=iris[iris$Species=="versicolor",1:4]
##' virginica=iris[iris$Species=="virginica",1:4]
##' MergeGUI(setosa,versicolor,virginica)
##' }
##'
MergeGUI = function(..., filenames=NULL, unit=TRUE, distn=TRUE, miss=TRUE) {
    mergegui_env = new.env()
    
    mergefunc = function(h, ...) {
        
        undo = function(h, ...) {
            #####-----------------------------------------------------#####
            ##  The following buttons are used for switching variables.  ##
            ##  undo button.                                             ##
            #####-----------------------------------------------------#####
            mergegui_env$idx <- mergegui_env$idx - 1
            if (mergegui_env$idx == 0) {
                gmessage("You can not undo anymore!")
                mergegui_env$idx <- 1
            }
            for (i in 1:n) {
                gt2[[i]][,] = data.frame(namecode=rownames(mergegui_env$hstry1[[mergegui_env$idx]]),mergegui_env$hstry1[[mergegui_env$idx]][, i, drop=FALSE],stringsAsFactors = FALSE)
            }
            mergegui_env$redo.indicate <- 1
            mergegui_env$gt4[,] = mergegui_env$hstry2[[mergegui_env$idx]]
            mergegui_env$gt5[,] = mergegui_env$hstry3[[mergegui_env$idx]]
        }
        
        redo = function(h, ...) {
            #####-----------------------------------------------------#####
            ##  The following buttons are used for switching variables.  ##
            ##  redo button.                                             ##
            #####-----------------------------------------------------#####
            if (mergegui_env$redo.indicate == 0) {
                gmessage("There is nothing to redo.")
                return()
            }
            mergegui_env$idx <- mergegui_env$idx + 1
            if (mergegui_env$idx > length(mergegui_env$hstry1)) {
                gmessage("You can not redo anymore!")
                mergegui_env$idx <- length(mergegui_env$hstry1)
            }
            for (i in 1:n) {
                gt2[[i]][,] = data.frame(namecode=rownames(mergegui_env$hstry1[[mergegui_env$idx]]),mergegui_env$hstry1[[mergegui_env$idx]][, i, drop=FALSE],stringsAsFactors = FALSE)
            }
            mergegui_env$gt4[,] = mergegui_env$hstry2[[mergegui_env$idx]]
            mergegui_env$gt5[,] = mergegui_env$hstry3[[mergegui_env$idx]]
        }
        
        reset = function(h, ...) {
            #####-----------------------------------------------------#####
            ##  The following buttons are used for switching variables.  ##
            ##  reset button.                                            ##
            #####-----------------------------------------------------#####
            for (i in 1:n) {
                gt2[[i]][,] = data.frame(namecode=rownames(nametable),nametable[, i, drop = F],stringsAsFactors = FALSE)
            }
            mergegui_env$redo.indicate = 1
            mergegui_env$gt4[, ] = mergegui_env$hstry2[[1]]
            mergegui_env$gt5[, ] = mergegui_env$hstry3[[1]]
            if (mergegui_env$hstry4[[length(mergegui_env$hstry4)]]==length(mergegui_env$hstry4)) {
                indicator1 = all(mergegui_env$hstry1[[length(mergegui_env$hstry1)]]==mergegui_env$hstry1[[1]],na.rm=TRUE)
                indicator2 = all(mergegui_env$hstry2[[length(mergegui_env$hstry2)]][,1:3]==mergegui_env$hstry2[[1]][,1:3],na.rm=TRUE)
                indicator3 = all(mergegui_env$hstry3[[length(mergegui_env$hstry3)]]==mergegui_env$hstry3[[1]],na.rm=TRUE)
                mergegui_env$idx = ifelse(all(indicator1,indicator2,indicator3),length(mergegui_env$hstry1),length(mergegui_env$hstry1)+1)
            } else {
                mergegui_env$idx = length(mergegui_env$hstry1)+1
            }            
            mergegui_env$hstry1[[mergegui_env$idx]] = mergegui_env$hstry1[[1]]
            mergegui_env$hstry2[[mergegui_env$idx]] = mergegui_env$hstry2[[1]]
            mergegui_env$hstry3[[mergegui_env$idx]] = mergegui_env$hstry3[[1]]
            mergegui_env$hstry4[[mergegui_env$idx]] = mergegui_env$idx
            svalue(check141,index=TRUE) = 1:3
        }
        
        VariableOptions = function(h, ...) {
            #####------------------------------------------------------#####
            ##  VariableOptions is the handler when double clicking gt4.  ##
            ##  It gives a new window for                                 ##
            ##          editing the attributes of variables.              ##
            #####------------------------------------------------------#####
            gt4input0 = gwindow("Attributes", visible = T, width = 300,
                                height = 200)
            gt4input = ggroup(horizontal = FALSE, container = gt4input0,
                              expand = TRUE)
            gt4input1 = ggroup(container = gt4input, expand = TRUE)
            gt4input2 = ggroup(container = gt4input, expand = TRUE)
            gt4input4 = ggroup(container = gt4input, horizontal = FALSE, expand = TRUE)
            gt4input3 = ggroup(container = gt4input, expand = TRUE)
            
            gt4input11 = glabel("Name:", container = gt4input1)
            gt4input12 = gedit(text = svalue(mergegui_env$gt4), container = gt4input1,
                               expand = TRUE)
            gt4input21 = glabel("Class:", container = gt4input2)
            gt4input22 = gcombobox(union(mergegui_env$gt4[svalue(mergegui_env$gt4, index = TRUE), 3], c("integer", "numeric", "character", "factor")),
                                   container = gt4input2, expand = TRUE)
            
            gt4input31 = gbutton("Ok", container = gt4input3, expand = TRUE,
                                 handler = function(h, ...) {
                                     if (svalue(gt4input12) != "") {
                                         mergegui_env$gt4[svalue(mergegui_env$gt4, index = TRUE), 2] = svalue(gt4input12)
                                         mergegui_env$gt4[svalue(mergegui_env$gt4, index = TRUE), 3] = svalue(gt4input22)
                                         mergegui_env$gt5[mergegui_env$gt4[svalue(mergegui_env$gt4, index = TRUE), 1], 2] = mergegui_env$gt4[svalue(mergegui_env$gt4, index = TRUE), 2]
                                         mergegui_env$gt5[mergegui_env$gt4[svalue(mergegui_env$gt4, index = TRUE), 1], 3] = mergegui_env$gt4[svalue(mergegui_env$gt4, index = TRUE), 3]
                                         mergegui_env$hstry2[[mergegui_env$idx]] <- mergegui_env$gt4[,]
                                         mergegui_env$hstry3[[mergegui_env$idx]] <- mergegui_env$gt5[,]
                                         dispose(gt4input0)
                                     }
                                     else {
                                         gmessage("Variable name could not be empty!")
                                     }
                                 })
            gt4input32 = gbutton("Cancel", container = gt4input3,
                                 expand = TRUE, handler = function(h, ...) {
                                     dispose(gt4input0)
                                 })
        }
        
        smmry = function(h, ...) {
            #####---------------------------------#####
            ##  smmry is the handler of gbcombo431.  ##
            ##  (gbutton: Numeric Summary)           ##
            #####---------------------------------#####
            graphics.off()
            name.select = svalue(mergegui_env$gt4, index = TRUE)
            
            if (length(name.select) == 0) {
                gmessage("Please select the variables!")
                return()
            }
            name.table = matrix(nrow = length(name.select), ncol = n)
            for (i in 1:n) {
                name.table[, i] = gt2[[i]][mergegui_env$gt4[name.select,1],2]
            }
            name.intersect = as.vector(svalue(mergegui_env$gt4))
            name.class = mergegui_env$gt4[name.select, 3]
            summarytable = list()
            for (i in 1:length(name.select)) {
                if (name.class[i] != "NA") {
                    if (name.class[i] == "numeric" | name.class[i] ==
                        "integer") {
                        summarytable[[i]] = matrix(NA, ncol = n, nrow = 7,
                                            dimnames = list(c("size", "NA#s", "mean",
                                            "std", "min", "median", "max"),
                                            simplifynames(gsub('.csv','',basename(gtfile)))))
                        names(summarytable)[i] = name.intersect[i]
                        for (j in 1:n) {
                            if (!is.na(name.table[i, j])) {
                                tmpdata = dataset[[j]][, name.table[i, j]]
                                summarytable[[i]][1, j] = length(tmpdata)
                                summarytable[[i]][2, j] = sum(is.na(tmpdata))
                                summarytable[[i]][3, j] = mean(tmpdata, na.rm = TRUE)
                                summarytable[[i]][4, j] = sd(tmpdata, na.rm = TRUE)
                                summarytable[[i]][5, j] = min(tmpdata, na.rm = TRUE)
                                summarytable[[i]][6, j] = median(tmpdata, na.rm = TRUE)
                                summarytable[[i]][7, j] = max(tmpdata, na.rm = TRUE)
                            }
                        }
                        summarytable[[i]] = data.frame(t(summarytable[[i]]))
                        summarytable[[i]][,1] = as.integer(as.character(summarytable[[i]][,1]))
                        summarytable[[i]][,2] = as.integer(as.character(summarytable[[i]][,2]))
                        summarytable[[i]][,3] = as.character(round(summarytable[[i]][,3]),3)
                        summarytable[[i]][,4] = as.character(round(summarytable[[i]][,4]),3)
                        summarytable[[i]][,5] = as.character(round(summarytable[[i]][,5]),3)
                        summarytable[[i]][,6] = as.character(round(summarytable[[i]][,6]),3)
                        summarytable[[i]][,7] = as.character(round(summarytable[[i]][,7]),3)
                        summarytable[[i]] = cbind(File=rownames(summarytable[[i]]),summarytable[[i]])
                    }
                    else {
                        summarytable[[i]] = matrix(NA, ncol = n, nrow = 10,
                                            dimnames = list(c("size", "NA#s", "levels",
                                            "matched levels", "top 1 level", "amount 1",
                                            "top 2 level", "amount 2",
                                            "top 3 level", "amount 3"),
                                            simplifynames(gsub('.csv','',basename(gtfile)))))
                        names(summarytable)[i] = name.intersect[i]
                        matchedlevels = list()
                        for (j in 1:n) {
                            if (!is.na(name.table[i, j])) {
                                if (sum(!is.na(dataset[[j]][, name.table[i, j]])) > 0) {
                                    matchedlevels[[j]] = names(table(dataset[[j]][,
                                                         name.table[i, j]], useNA = "no"))
                                }
                                else {
                                    matchedlevels[[j]] = NA
                                }
                            }
                            else {
                                matchedlevels[[j]] = NA
                            }
                        }
                        mtch = intersect2(matchedlevels, matchedlevels)
                        for (j in 1:n) {
                            if (!is.na(name.table[i, j])) {
                                tmpdata = dataset[[j]][, name.table[i, j]]
                                tmptable = sort(table(tmpdata, useNA = "no"),
                                                decreasing = TRUE)
                                summarytable[[i]][1, j] = length(tmpdata)
                                summarytable[[i]][2, j] = sum(is.na(tmpdata))
                                summarytable[[i]][3, j] = length(tmptable)
                                summarytable[[i]][4, j] = length(mtch$public)
                                if (length(tmptable) > 0) {
                                    summarytable[[i]][5, j] = names(tmptable)[1]
                                    summarytable[[i]][6, j] = tmptable[1]
                                }
                                if (length(tmptable) > 1) {
                                    summarytable[[i]][7, j] = names(tmptable)[2]
                                    summarytable[[i]][8, j] = tmptable[2]
                                }
                                if (length(tmptable) > 2) {
                                    summarytable[[i]][9, j] = names(tmptable)[3]
                                    summarytable[[i]][10, j] = tmptable[3]
                                }
                            }
                        }
                        summarytable[[i]] = data.frame(t(summarytable[[i]]))
                        summarytable[[i]][,1] = as.integer(as.character(summarytable[[i]][,1]))
                        summarytable[[i]][,2] = as.integer(as.character(summarytable[[i]][,2]))
                        summarytable[[i]][,3] = as.integer(as.character(summarytable[[i]][,3]))
                        summarytable[[i]][,4] = as.integer(as.character(summarytable[[i]][,4]))
                        summarytable[[i]][,6] = as.integer(as.character(summarytable[[i]][,6]))
                        summarytable[[i]][,8] = as.integer(as.character(summarytable[[i]][,8]))
                        summarytable[[i]][,10] = as.integer(as.character(summarytable[[i]][,10]))
                        summarytable[[i]] = cbind(File=rownames(summarytable[[i]]),summarytable[[i]])
                    }
                }
                else {
                    summarytable[[i]] = matrix(NA, nrow = n, ncol = 1,
                                               dimnames = list(simplifynames(gsub('.csv','',basename(gtfile))), NULL))
                    names(summarytable)[i] = name.intersect[i]
                }
            }
            
            gbcombo441 = list()
            
            delete(mergegui_env$group43, mergegui_env$group45)
            mergegui_env$group45 <- ggroup(container=mergegui_env$group43,expand = TRUE, use.scrollwindow = TRUE)
            gbcombo44 <- glayout(container = mergegui_env$group45,expand = TRUE, use.scrollwindow = TRUE)
            
            for (i in 1:length(name.select)){
                gbcombo44[i*2-1, 1] = glabel(names(summarytable)[i],container=gbcombo44)
                gbcombo44[i*2, 1, expand = TRUE] = gbcombo441[[i]] = gtable(summarytable[[i]], container = gbcombo44)
            }
            
        }
        
        graph = function(h, ...) {
            #####---------------------------------#####
            ##  graph is the handler of gbcombo432.  ##
            ##  (gbutton: Graphic Summary)           ##
            #####---------------------------------#####
            graphics.off()
            delete(mergegui_env$group43, mergegui_env$group45)
            mergegui_env$group45 <- ggroup(container=mergegui_env$group43,expand = TRUE, use.scrollwindow = TRUE)
            gbcombo44 <- glayout(container = mergegui_env$group45,expand = TRUE, use.scrollwindow = TRUE)
            
            yscale = svalue(radio121)
            name.select = svalue(mergegui_env$gt4, index = TRUE)
            if (length(name.select)==0) {
                gmessage("Please select one variables!")
                return()
            }
            if (length(name.select)==1) {
                name.table = rep(NA, n)
                for (i in 1:n) {
                    name.table[i] = gt2[[i]][mergegui_env$gt4[name.select,1], 2]
                }
                name.intersect = as.character(mergegui_env$gt4[name.select,2])
                name.class = as.character(mergegui_env$gt4[name.select, 3])
                mergedata = data.frame(source = rep(simplifynames(gsub('.csv','',basename(gtfile))),
                                                    rows))
                
                is.num = FALSE
                if (name.class != "NA") {
                    if (name.class %in% c("numeric","integer")) {
                        is.num = TRUE
                        tmp.num = c()
                        for (i in 1:n) {
                            if (!is.na(name.table[i])) {
                                tmp.num = c(tmp.num, dataset[[i]][, name.table[i]])
                            }
                            else {
                                tmp.num = c(tmp.num, rep(NA, rows[i]))
                            }
                        }
                        mergedata = data.frame(mergedata, as.numeric(tmp.num))
                        mergedata[,1] = reorder(mergedata[,1], mergedata[,2], median, na.rm=TRUE)
                    }
                    else {
                        tmp.chr = rep("na", sum(rows))
                        for (i in 1:n) {
                            if (!is.na(name.table[i])) {
                                tmp.chr[(cumsum(rows) - rows + 1)[i]:cumsum(rows)[i]] = as.character(dataset[[i]][,
                                                                                                                  name.table[i]])
                            }
                            else {
                                tmp.chr[(cumsum(rows) - rows + 1)[i]:cumsum(rows)[i]] = rep(NA,
                                                                                            rows[i])
                            }
                        }
                        levelorder = names(sort(table(tmp.chr), decreasing = FALSE))
                        mergedata = cbind(mergedata, factor(tmp.chr,
                                                            levels = levelorder))
                    }
                }
                else {
                    mergedata = cbind(mergedata, rep(NA, sum(rows)))
                }
                colnames(mergedata) = c("source", name.intersect)
                
                gbcombo44[1, 1, expand = TRUE] = gbcombo442 = ggroup(container = gbcombo44, use.scrollwindow = TRUE)
                gbcombo4421 = ggraphics(container = gbcombo442,
                                        height = ifelse(is.num, 75 * 3 * n, 75 * 6),  expand = TRUE)
                
                if (yscale=="regular y scale") {
                    eval(parse(text = paste("print(qplot(", name.intersect,
                                            ",data=mergedata,facets=", ifelse(is.num,
                                                                              "source~.)", "~source)+coord_flip()"), ")", collapse = "")))
                } else {
                    eval(parse(text = paste("print(qplot(", name.intersect,
                                            ",data=mergedata,geom='histogram')+facet_wrap(~source, scales = 'free_y', ncol = 1)",
                                            ifelse(is.num, "", "+coord_flip()"), ")", collapse = "")))
                }
            }
            if (length(name.select)==2) {
                name.table = matrix(NA, ncol=n, nrow=2)
                for (i in 1:n) {
                    name.table[,i] = gt2[[i]][mergegui_env$gt4[name.select,1], 2]
                }
                name.intersect = as.character(mergegui_env$gt4[name.select,2])
                name.class = as.character(mergegui_env$gt4[name.select, 3])
                mergedata = data.frame(source = rep(simplifynames(gsub('.csv','',basename(gtfile))),rows))
                for (j in 1:2) {
                    tmp.num = c()
                    for (i in 1:n) {
                        if (!is.na(name.table[j,i])) {
                            tmp.num = c(tmp.num, as.character(dataset[[i]][, name.table[j,i]]))
                        } else {
                            tmp.num = c(tmp.num, rep(NA, rows[i]))
                        }
                    }
                    mergedata = cbind(mergedata,tmp.num)
                    colnames(mergedata)[j+1] = name.intersect[j]
                    eval(parse(text=paste("mergedata[,j+1] = as.",name.class[j],"(as.character(mergedata[,j+1]))",sep="")))
                }
                mergedata = data.frame(mergedata)
                colnames(mergedata)[1]="source"
                mergedata$source = factor(mergedata$source)
                
                gbcombo44[1, 1, expand = TRUE] = gbcombo442 = ggroup(container = gbcombo44, use.scrollwindow = TRUE)
                gbcombo4421 = ggraphics(container = gbcombo442, height = 75 * 3 * n,  expand = TRUE)
                
                if (all(name.class %in% c("integer","numeric"))){
                    if (yscale=="regular y scale") {
                        eval(parse(text = paste("print(qplot(", name.intersect[1],",", name.intersect[2],",data=mergedata,geom='point',facets=source~., alpha=I(0.6)))", sep = "")))
                    } else {
                        eval(parse(text = paste("print(qplot(", name.intersect[1],",", name.intersect[2],",data=mergedata,geom='point', alpha=I(0.6))+facet_wrap(~source, scales = 'free', ncol = 1))", sep = "")))
                    }
                } else {
                    if (all(name.class %in% c("factor","character"))){
                        if (yscale=="regular y scale") {
                            eval(parse(text = paste("print(qplot(", name.intersect[1],",", name.intersect[2],",data=mergedata,geom='point', position=position_jitter(w=0.2,h=0.2), facets=source~., alpha=I(0.6)))", sep = "")))
                        } else {
                            eval(parse(text = paste("print(qplot(", name.intersect[1],",", name.intersect[2],",data=mergedata,geom='point', position=position_jitter(w=0.2,h=0.2), alpha=I(0.6)) + facet_wrap(~source,scales='free',ncol=1))", sep = "")))
                        }
                    } else {
                        if (yscale=="regular y scale") {
                            eval(parse(text = paste("print(qplot(", name.intersect[which(name.class %in% c("factor","character"))],",", name.intersect[which(name.class %in% c("integer","numeric"))],",data=mergedata,geom=c('boxplot','point'),facets=source~.)+coord_flip())", sep = "")))
                        } else {
                            eval(parse(text = paste("print(qplot(", name.intersect[which(name.class %in% c("factor","character"))],",", name.intersect[which(name.class %in% c("integer","numeric"))],",data=mergedata,geom=c('boxplot','point'))+facet_wrap(~source,scales='free',ncol=1)+ coord_flip())", sep = "")))
                        }
                    }
                }
            }
            if (length(name.select)>2) {
                z = length(name.select)
                name.table = matrix(NA, ncol=n, nrow=z)
                for (i in 1:n) {
                    name.table[,i] = gt2[[i]][mergegui_env$gt4[name.select,1], 2]
                }
                name.intersect = as.character(mergegui_env$gt4[name.select,2])
                name.class = as.character(mergegui_env$gt4[name.select, 3])
                mergedata = data.frame(source = rep(simplifynames(gsub('.csv','',basename(gtfile))),rows))
                
                for (j in 1:z) {
                    tmp.num = c()
                    for (i in 1:n) {
                        if (!is.na(name.table[j,i])) {
                            tmp.num = c(tmp.num, as.character(dataset[[i]][, name.table[j,i]]))
                        } else {
                            tmp.num = c(tmp.num, rep(NA, rows[i]))
                        }
                    }
                    mergedata = cbind(mergedata,tmp.num)
                    colnames(mergedata)[j+1] = name.intersect[j]
                    eval(parse(text=paste("mergedata[,j+1] = as.",name.class[j],"(as.character(mergedata[,j+1]))",sep="")))
                }
                mergedata = data.frame(mergedata)
                colnames(mergedata)[1]="source"
                mergedata$source = factor(mergedata$source)
                
                if (sum(name.class %in% c("integer","numeric"))<z) {
                    for (i in 1:z) {
                        is.num = name.class[i] %in% c("integer","numeric")
                        
                        gbcombo44[i, 1, expand = TRUE] = ggraphics(container = gbcombo44, height = ifelse(is.num, 75 * 3 * n, 75 * 6),  expand = TRUE)
                        
                        if (yscale=="regular y scale") {
                            eval(parse(text = paste("print(qplot(", name.intersect[i], ",data=mergedata,facets=", ifelse(is.num, "source~.)", "~source)+coord_flip()"), ")", sep="")))
                        } else {
                            eval(parse(text = paste("print(qplot(", name.intersect[i], ",data=mergedata, geom='histogram')+ facet_wrap(~source, scales='free_y', ncol=1)", ifelse(is.num, "", "+coord_flip()"), ")", sep="")))
                        }
                    }
                } else {
                    gbcombo44[1, 1, expand = TRUE] = ggraphics(container = gbcombo44, expand = TRUE)
                    print(ggpcp(mergedata,vars=names(mergedata)[2:(z+1)]) + geom_line() + facet_wrap(~source, ncol=1))
                }
            }
        }
        
        dict = function(h, ...) {
            #####--------------------------------#####
            ##  dict is the handler of gbcombo432.  ##
            ##  (gbutton: Dictionary)               ##
            #####--------------------------------#####
            graphics.off()
            delete(mergegui_env$group43, mergegui_env$group45)
            mergegui_env$group45 <- ggroup(container=mergegui_env$group43,expand = TRUE, use.scrollwindow = TRUE)
            gbcombo44 <- glayout(container = mergegui_env$group45,expand = TRUE, use.scrollwindow = TRUE)
            gbcombo44[1, 1, expand = TRUE] = gbcombo443 = gtext(container = gbcombo44, expand = TRUE,
                                                                use.scrollwindow = TRUE)
            
            name.select = svalue(mergegui_env$gt4, index = TRUE)
            if (length(name.select) == 0) {
                gmessage("Please select the variables!")
                gbcombo443[, ] = data.frame(VarName = character(0),
                                            Level = integer(0), Label = character(0), stringsAsFactors = FALSE)
                return()
            }
            name.table = matrix(nrow = length(name.select), ncol = n)
            for (i in 1:n) {
                name.table[, i] = gt2[[i]][mergegui_env$gt4[name.select,1], 2]
            }
            name.intersect = as.vector(svalue(mergegui_env$gt4))
            name.class = mergegui_env$gt4[name.select, 3]
            
            dictionary = list()
            dictlength = matrix(0, nrow = length(name.intersect),
                                ncol = n)
            for (i in 1:length(name.intersect)) {
                dictionary[[i]] = list()
                names(dictionary)[i] = name.intersect[i]
                if (name.class[i] == "factor") {
                    for (j in 1:n) {
                        dictionary[[i]][[j]] = levels(factor(dataset[[j]][,
                                                                          name.table[i, j]]))
                        dictlength[i, j] = length(dictionary[[i]][[j]])
                    }
                }
            }
            
            dictlist = list()
            for (i in 1:length(name.intersect)) {
                if (max(dictlength[i, ]) != 0) {
                    dictlist[[i]] = matrix(NA, nrow = max(dictlength[i,
                                                                     ]), ncol = n)
                    names(dictlist)[i] = name.intersect[i]
                    rownames(dictlist[[i]]) = 1:nrow(dictlist[[i]])
                    levelintersect = intersect2(dictionary[[i]],dictionary[[i]])
                    for (j in 1:n) {
                        if (dictlength[i, j]>0) {
                            dictlist[[i]][1:dictlength[i, j], j] = c(levelintersect$individual[,j],levelintersect$uniq[[j]])
                        }
                    }
                    colnames(dictlist[[i]]) = simplifynames(gsub('.csv','',basename(gtfile)))
                }
                else {
                    dictlist[[i]] = "Not a factor"
                    names(dictlist)[i] = name.intersect[i]
                }
            }
            
            if (sum(dictlength) == 0) {
                gmessage("All the variables selected are not factor variables.")
                svalue(gbcombo443) = ""
                return()
            }
            else {
                svalue(gbcombo443) = capture.output(noquote(dictlist))
            }
        }
        
        changetest = function(h,...) {
            #####-----------------------------------#####
            ##  changetest is the handler of radio131  ##
            ##  (gradio: Flag for variables)           ##
            #####-----------------------------------#####
            flagsym = svalue(radio131)
            
            if (flagsym=="Do not show p-values or flags") {
                newgt4 = mergegui_env$gt4[,1:3]
                delete(mergegui_env$group42, mergegui_env$gt4)
                mergegui_env$gt4 <- gtable(newgt4, multiple = T, container = mergegui_env$group42, expand = TRUE, chosencol = 2)
                addhandlerdoubleclick(mergegui_env$gt4, handler = VariableOptions)
                return()
            }
            
            gt4col1 = rownames(mergegui_env$gt4)
            if (!exists("namepanel",where=mergegui_env)) {
                mergegui_env$namepanel = nametable
                mergegui_env$name_intersection_panel[,2]=mergegui_env$gt4[gt4col1,2]
            } else {
                checknamepanel=c()
                for (i in 1:n) {
                    checknamepanel[i]=all(mergegui_env$namepanel[,i]==gt2[[i]][,2], na.rm = TRUE)
                    if (!checknamepanel[i]) mergegui_env$namepanel[,i]<-gt2[[i]][,2]
                }
                checknamepanel[n+1]=all(mergegui_env$name_intersection_panel[,3]==mergegui_env$gt4[order(gt4col1),3], na.rm = TRUE)
                if (!all(checknamepanel)){
                    mergegui_env$nameintersection <- mergegui_env$gt4[order(gt4col1),2]
                    mergegui_env$name_intersection_panel <- data.frame(mergegui_env$gt4[order(gt4col1),1:3],stringsAsFactors = FALSE)
                    colnames(mergegui_env$name_intersection_panel) <- c("Namecode", "Variables", "Class")
                    if (unit) mergegui_env$name_intersection_panel$Unit <- scale_rpart(mergegui_env$namepanel, dataset, mergegui_env$nameintersection, mergegui_env$gt4[order(gt4col1),3])
                    
                    if (distn) mergegui_env$name_intersection_panel$Dist <- scale_kstest(mergegui_env$namepanel, dataset, mergegui_env$nameintersection, mergegui_env$gt4[order(gt4col1),3])
                    if (miss) mergegui_env$name_intersection_panel$Miss <- scale_missing(mergegui_env$namepanel, dataset, mergegui_env$nameintersection)
                }
                mergegui_env$name_intersection_panel[,2]=mergegui_env$gt4[order(gt4col1),2]
            }
            delete(mergegui_env$group42, mergegui_env$gt4)
            
            if (flagsym=="Show the flag symbol") {
                alphalevel = as.numeric(svalue(text133))
                if (is.na(alphalevel) || alphalevel<=0 || alphalevel>=1){
                    gmessage("Invalid alpha-level. Coerce to 0.05.")
                    svalue(text133) = 0.05
                    alphalevel = 0.05
                }
                flag1 = !is.na(mergegui_env$name_intersection_panel[,-(1:3)])
                flag2 = sapply(mergegui_env$name_intersection_panel[,-(1:3)],
                               function(avec){
                                   as.numeric(as.character(avec)) <= alphalevel
                               })
                flag = flag1 & flag2
                newgt4 = mergegui_env$name_intersection_panel
                for (j in 4:ncol(newgt4)) {
                    newgt4[,j] = as.character(newgt4[,j])
                    newgt4[flag[,j-3],j] <- "X"
                    newgt4[!flag[,j-3],j] <- ""
                }
                mergegui_env$gt4 <- gtable(newgt4[gt4col1,], multiple = T, container = mergegui_env$group42, expand = TRUE, chosencol = 2)
                addhandlerdoubleclick(mergegui_env$gt4, handler = VariableOptions)
            } else {
                mergegui_env$gt4 <- gtable(mergegui_env$name_intersection_panel[gt4col1,], multiple = T,container = mergegui_env$group42, expand = TRUE, chosencol = 2)
                addhandlerdoubleclick(mergegui_env$gt4, handler = VariableOptions)
            }
            
        }
        
        changematching = function(h,...) {
            #####---------------------------------------#####
            ##  changematching is the handler of check141  ##
            ##  (gcheckboxgroup: View mode)                ##
            #####---------------------------------------#####
            viewmode = svalue(check141)
            if (length(viewmode)==0) {
                gmessage("You have to see something there. Please check at least one box.")
                return()
            }
            vistable = mergegui_env$hstry1[[mergegui_env$hstry4[[mergegui_env$idx]]]]
            visname = rownames(vistable)
            vispart = c()
            if ("Matched variables" %in% viewmode){
                tmppart = visname[grep("Part1-1-",visname)]
                if (length(tmppart)>0) vispart = c(vispart, tmppart)
            }
            if ("Partial-matched variables" %in% viewmode){
                tmppartidx = c(grep("Part1-1-",visname),grep(paste("Part",n,"-",sep=""),visname))
                if (length(tmppartidx)) {
                    tmppart = visname[-tmppartidx]
                    vispart = c(vispart, tmppart)
                }
            }
            if ("Unmatched variables" %in% viewmode){
                tmppart = visname[grep(paste("Part",n,"-",sep=""),visname)]
                if (length(tmppart)>0) vispart = c(vispart, tmppart)
            }
            if (length(vispart)==0) {
                gmessage("No rows are selected. Please check one more box.")
                return()
            }
            for (i in 1:n) {
                gt2[[i]][,] = gt2[[i]][1:length(vispart),]
                gt2[[i]][,1] = vispart
                gt2[[i]][,2] = vistable[vispart, i]
            }
            mergegui_env$idx = mergegui_env$idx + 1
            mergegui_env$hstry1[[mergegui_env$idx]] = vistable[vispart,]
            mergegui_env$hstry2[[mergegui_env$idx]] = mergegui_env$gt4[,]
            mergegui_env$hstry3[[mergegui_env$idx]] = mergegui_env$gt5[,]
            mergegui_env$hstry4[[mergegui_env$idx]] = ifelse(length(viewmode)==3,mergegui_env$idx,mergegui_env$hstry4[[mergegui_env$idx-1]])
        }
        
        watchdatafunc = function(h, ...) {
            #####-------------------------------------------------------#####
            ##  watchdatafunc is a function to export the merged dataset.  ##
            ##  For the selected checkboxs, we export the corresponding    ##
            ##          variables from all files.                          ##
            ##  The public name for the selected variable is the shortest  ##
            ##          name of that variable among different files.       ##
            ##  mergedata is a matrix to save the merged dataset.          ##
            ##  We should write 'xxx.csv'                                  ##
            ##          when we export mergedata and save the file.        ##
            #####-------------------------------------------------------#####
            name.select = svalue(mergegui_env$gt5, index = TRUE)
            if (length(name.select) == 0) {
                gmessage("Please select the variables!")
                return()
            }
            txtpb = txtProgressBar(min=0,max=1,width = 40,style=3)
            name.table = matrix(nrow = length(name.select), ncol = n)
            for (i in 1:n) {
                name.table[, i] = gt2[[i]][name.select, 2]
            }
            colnames(name.table)=gsub("\\.csv$","",basename(gtfile))
            name.intersect = as.vector(svalue(mergegui_env$gt5))
            name.class = mergegui_env$gt5[name.select, 3]
            mergedata = matrix(nrow = sum(rows), ncol = nrow(name.table) + 1)
            colnames(mergedata) = c("source", name.intersect)
            mergedatadictionary = data.frame(namecode=gt2[[1]][name.select, 1],
                                             newname=name.intersect,
                                             class=name.class,
                                             name.table, stringsAsFactors=FALSE)
            rownames(mergedatadictionary) = name.intersect
            setTxtProgressBar(txtpb, 0.05)
            for (i in 1:n) {
                tmp = matrix(c(rep(gsub("\\.csv$","",basename(gtfile[i])), rows[i]),
                               rep(NA, rows[i] * nrow(name.table))), nrow = rows[i])
                colnames(tmp) = c("source", name.table[, i])
                tmp[, na.omit(name.table[, i])] = as.matrix(dataset[[i]])[,
                                                  na.omit(name.table[, i])]
                mergedata[(cumsum(rows) - rows + 1)[i]:cumsum(rows)[i],] = tmp
                mergedatadictionary[,paste(colnames(name.table)[i],"index",sep="_")]=NA
                mergedatadictionary[,3+n+i]=sapply(name.table[,i],function(x){
                    ifelse(is.na(x),NA,which(colnames(dataset[[i]])==x))
                })
                setTxtProgressBar(txtpb, (0.05+0.4*i/n))
            }
            
            mergedatasummary = matrix(c(colnames(name.table), rows), nrow = n,
                                      dimnames = list(colnames(name.table), c("source","size")))
            for (i in 1:length(name.select)) {
                if (name.class[i] != "NA") {
                    if (name.class[i] == "numeric" | name.class[i] ==
                        "integer") {
                        if (name.class[i] == "numeric") {
                            mergedata[, i + 1] = as.numeric(mergedata[, i + 1])
                        }
                        if (name.class[i] == "integer") {
                            mergedata[, i + 1] = as.integer(mergedata[, i + 1])
                        }
                        datasummary = matrix(NA, nrow = n, ncol = 6,
                                             dimnames = list(colnames(name.table), 
                                             paste(name.intersect[i], 
                                             c("NA#s", "mean", "std", "min", "median", "max"),
                                             sep = ".")))
                        for (j in 1:n) {
                            if (!is.na(name.table[i, j])) {
                                tmpdata = dataset[[j]][, name.table[i,
                                                                    j]]
                                datasummary[j, 1] = sum(is.na(tmpdata))
                                datasummary[j, 2] = mean(tmpdata, na.rm = TRUE)
                                datasummary[j, 3] = sd(tmpdata, na.rm = TRUE)
                                datasummary[j, 4] = min(tmpdata, na.rm = TRUE)
                                datasummary[j, 5] = median(tmpdata, na.rm = TRUE)
                                datasummary[j, 6] = max(tmpdata, na.rm = TRUE)
                            }
                            #setTxtProgressBar(txtpb, 0.45+(i-1)/n*0.5+j/n*0.5/n)
                        }
                        mergedatasummary = cbind(mergedatasummary,
                                                 datasummary)
                    }
                    else {
                        datasummary = matrix(NA, nrow = n, ncol = 9,
                                             dimnames = list(colnames(name.table),
                                             paste(name.intersect[i], 
                                             c("NA#s", "levels", "matched_levels",
                                             "top1_level", "amount_1",
                                             "top2_level", "amount_2",
                                             "top3_level", "amount_3"), sep = ".")))
                        matchedlevels = list()
                        for (j in 1:n) {
                            if (!is.na(name.table[i, j])) {
                                if (sum(!is.na(dataset[[j]][, name.table[i, j]])) > 0) {
                                    matchedlevels[[j]] = names(table(dataset[[j]][,
                                                         name.table[i, j]], useNA = "no"))
                                }
                                else {
                                    matchedlevels[[j]] = NA
                                }
                            }
                            else {
                                matchedlevels[[j]] = NA
                            }
                            #setTxtProgressBar(txtpb, 0.45+(i-1)/n*0.5+j/n*0.25/n)
                        }
                        mtch = intersect2(matchedlevels, matchedlevels)
                        for (j in 1:n) {
                            if (!is.na(name.table[i, j])) {
                                tmpdata = dataset[[j]][, name.table[i, j]]
                                tmptable = sort(table(tmpdata, useNA = "no"),
                                                decreasing = TRUE)
                                datasummary[j, 1] = sum(is.na(tmpdata))
                                datasummary[j, 2] = length(tmptable)
                                datasummary[j, 3] = length(mtch$public)
                                if (length(tmptable) > 0){
                                    datasummary[j, 4] = names(tmptable)[1]
                                    datasummary[j, 5] = tmptable[1]
                                }
                                if (length(tmptable) > 1) {
                                    datasummary[j, 6] = names(tmptable)[2]
                                    datasummary[j, 7] = tmptable[2]
                                }
                                if (length(tmptable) > 2) {
                                    datasummary[j, 8] = names(tmptable)[3]
                                    datasummary[j, 9] = tmptable[3]
                                }
                            }
                            #setTxtProgressBar(txtpb, 0.45+(i-0.5)/n*0.5+j/n*0.25/n)
                        }
                        mergedatasummary = cbind(mergedatasummary, datasummary)
                    }
                }
                else {
                    mergedatasummary = cbind(mergedatasummary, 
                                             matrix(NA, ncol = 1, nrow = n, 
                                             dimnames = list(NULL, name.intersect[i])))
                }
                setTxtProgressBar(txtpb, 0.45+i/n*0.5)
            }
            
            setTxtProgressBar(txtpb, 1)
            if (!is.na(gf <- gfile(type = "save"))) {
                if (regexpr("\\.csv$",gf) %in% c(-1,1)) {
                    gf = paste(gf,".csv",sep="")
                }
                write.csv(mergedata, file = gf, row.names = FALSE)
                summarylocation = sub("\\.csv$", "_summary.csv", gf)
                write.table(t(mergedatasummary), file = summarylocation,
                            sep=",", col.names = FALSE)
                dictionarylocation = sub("\\.csv$", "_dictionary.csv", gf)
                write.csv(mergedatadictionary, file = dictionarylocation,
                          row.names = FALSE)
                gmessage("The files are merged!")
            }
        }
        
        if (exists("combo2",where=mergegui_env,inherits=FALSE)) {
            #if (!isExtant(mergegui_env$combo2)) {
            dispose(mergegui_env$combo2)
        }
        #####-------------------------------------------------------------------#####
        ##  Import the selected files.                                             ##
        ##  'dataset' is a list to save all the data from different files.         ##
        ##  'rows' is  a vector to save the number of observations for each file.  ##
        ##  'vname' is  a list to save the original colnames of the dataset.       ##
        ##  'simplifiedname' is  a list to save the simplified name.               ##
        ##          (delete the filenames in the colnames if they have)            ##
        ##  'vname' & 'simplifiedname' are 1-1 projections,                        ##
        ##          although 'simplifiedname' haves repeated names.                ##
        #####-------------------------------------------------------------------#####
        dataset <- list()
        simplifiedname <- list()
        vname <- list()
        
        if (length(svalue(gt)) == 0) {
            n <- length(gt[])
            gtfile <- gt[]
        } else {
            n <- length(svalue(gt))
            gtfile <- svalue(gt)
        }
        
        rows <- rep(0, n)
        if (n<2) {
            warning('The input data set is not enough. More files are needed.')
            return()
        }
        if (n>9) {
            warning('Too many data sets! Please limit the number to 9.')
            return()
        }
        
        for (i in 1:n) {
            dataset[[i]] <- if (length(grep("\\.csv$",gtfile[i]))) {
                read.csv(file = gtfile[i], header = T)
            } else { gtdata[[i]] }
            rows[i] <- nrow(dataset[[i]])
            vname[[i]] <- colnames(dataset[[i]])
            J = gregexpr(".csv.", vname[[i]])
            loc = c()
            for (j in 1:length(vname[[i]])) {
                loc[j] = max(J[[j]]) + 5
            }
            loc[which(loc == 4)] = 1
            simplifiedname[[i]] = substring(vname[[i]], loc)
        }
        
        #####----------------------------------------------------#####
        ##  Now we are going to generate two stuffs:                ##
        ##          'nameintersect' and 'nametable'.                ##
        ##  'nametable' is a matrix with all matched and            ##
        ##          unmatched variable names ('vname').             ##
        ##  'nameintersect' is a vector of mutually exclusive       ##
        ##          and collectively exhaustive names,              ##
        ##          with 'simplifiedname' at the beginning,         ##
        ##          and the special 'vname' at the end.             ##
        ##  length(nameintersect) == nrow(nametable)                ##
        ##  'Part1-1-' is the intersection for all the n files.     ##
        ##  'Part(n-i+1)' is the intersection for all combination   ##
        ##	        of the n-i+1 files.                             ##
        ##  'Partn-i-' is the left part of the i-th files,          ##
        ##          it cannot be intersected with any other files.  ##
        #####----------------------------------------------------#####
        a = intersect2(vname, simplifiedname)
        nameintersect = a$public
        tmpuniq = a$uniq
        tmpsimpleuniq = a$simpleuniq
        nametable = a$individual
        if (nrow(nametable) != 0) {
            tmpname = paste("Part1-1-", sprintf("%03d", 1:nrow(nametable)), sep = "")
            rownames(nametable) = tmpname
        }
        
        for (i in max((n - 1), 2):2) {
            combnmatrix = combn(1:n, i)
            for (j in 1:ncol(combnmatrix)) {
                tmpintersect = intersect2(tmpuniq[combnmatrix[, j]],
                                          tmpsimpleuniq[combnmatrix[, j]])
                tmptable = matrix(NA, ncol = n, nrow = length(tmpintersect$public))
                tmptable[, combnmatrix[, j]] = tmpintersect$individual
                if (nrow(tmptable) != 0) {
                    tmpname = paste("Part", n - i + 1, "-", j, "-", 
                                    sprintf("%03d",1:length(tmpintersect$public)), sep = "")
                    rownames(tmptable) = tmpname
                }
                nametable <- rbind(nametable, tmptable)
                nameintersect <- c(nameintersect, tmpintersect$public)
                tmpuniq[combnmatrix[, j]] = tmpintersect$uniq
                tmpsimpleuniq[combnmatrix[, j]] = tmpintersect$simpleuniq
            }
        }
        nameintersect <- c(nameintersect, unlist(tmpuniq))
        
        for (i in 1:n) {
            tmptable = matrix(NA, ncol = n, nrow = length(tmpuniq[[i]]))
            tmptable[, i] = tmpuniq[[i]]
            if (nrow(tmptable) != 0) {
                tmpname = paste("Part", n, "-", i, "-",
                                sprintf("%03d",1:nrow(tmptable)), sep = "")
                rownames(tmptable) = tmpname
            }
            nametable <- rbind(nametable, tmptable)
        }
        colnames(nametable) <- simplifynames(gsub('.csv','',basename(gtfile)))
        
        #####-------------------------------#####
        ##  New window for matching variables  ##
        #####-------------------------------#####
        mergegui_env$combo2 = gwindow("Matched Variables", visible = T, width = 900, height = 600)
        tab = gnotebook(container = mergegui_env$combo2)
        if (exists("name_intersection_panel",where=mergegui_env)){
            rm("name_intersection_panel",envir=mergegui_env)
        }
        
        #####-----------------------------------------#####
        ##  In the first tab we can:                     ##
        ##  (1) Determine the scaling way.               ##
        ##  (2) Determine whether show p-values or not.  ##
        #####-----------------------------------------#####
        group11 = ggroup(horizontal = FALSE, container = tab, label = "Preferences", expand = T)
        frame12 = gframe("Scaling of histograms",container = group11, horizontal = FALSE)
        radio121 = gradio(c("regular y scale","relative y scale"),container = frame12)
        frame13 = gframe("Flag for variables",container = group11, horizontal = TRUE)
        radio131 = gradio(c("Show p-values","Show the flag symbol","Do not show p-values or flags"), container = frame13)
        label132 = glabel('"alpha-level" = ',container=frame13)
        text133 = gedit("0.05",container=frame13, width=10)
        if (!unit & !distn & !miss) {
            svalue(radio131) = "Do not show p-values or flags"
            enabled(frame13) = FALSE
            visible(label132) = FALSE
            visible(text133) = FALSE
        } else {
            addHandlerChanged(radio131, handler=changetest)
        }
        frame14 = gframe("View mode of the matching tab",container = group11, horizontal = FALSE)
        check141 = gcheckboxgroup(c("Matched variables","Partial-matched variables","Unmatched variables"), checked = TRUE, container = frame14, handler = changematching)
        
        #####----------------------------------------------#####
        ##  In the second tab we can:                         ##
        ##  (1) Switch the variable names in the same gtable. ##
        ##  (2) Go back or go forth or reset the matching.    ##
        #####----------------------------------------------#####
        group21 = ggroup(horizontal = FALSE, container = tab, label = "Matching", expand = T)
        group22 = ggroup(container = group21, use.scrollwindow = TRUE, expand = T)
        group2 = list()
        gt2 <- list()
        
        mergegui_env$hstry1 <- list()
        mergegui_env$hstry1[[1]] <- nametable
        
        mergegui_env$name_intersection_panel <- data.frame(
            Namecode=rownames(nametable), 
            Variables=nameintersect,
            Class=var.class(nametable,dataset),
            stringsAsFactors = FALSE)
        if (unit) mergegui_env$name_intersection_panel$Unit = scale_rpart(nametable,dataset,nameintersect)
        if (distn) mergegui_env$name_intersection_panel$Dist = scale_kstest(nametable,dataset,nameintersect)
        if (miss) mergegui_env$name_intersection_panel$Miss = scale_missing(nametable,dataset,nameintersect)
        mergegui_env$hstry2 <- list()
        mergegui_env$hstry2[[1]] <- mergegui_env$name_intersection_panel
        
        Matched = substr(rownames(nametable),5,regexpr('-',rownames(nametable))-1)
        FileMatched = as.character((n+1)-as.integer(Matched))
        mergegui_env$hstry3 <- list()
        mergegui_env$hstry3[[1]] <- data.frame(mergegui_env$name_intersection_panel[,1:3],FileMatched)
        
        mergegui_env$hstry4 <- list()
        mergegui_env$hstry4[[1]] <- 1
        
        mergegui_env$idx <- 1
        mergegui_env$redo.indicate <- 0
        
        for (i in 1:n) {
            group2[[i]] = ggroup(horizontal = FALSE, container = group22,
                                 expand = T)
            gt2[[i]] <- gtable(data.frame(namecode=rownames(nametable),nametable[, i, drop = F],stringsAsFactors = FALSE), chosencol = 2, container = group2[[i]], expand = TRUE)
            addHandlerKeystroke(gt2[[i]], handler = function(h,...){})
            tag(gt2[[i]], "prev.idx") <- svalue(gt2[[i]], index = TRUE)
            tag(gt2[[i]], "toggle") <- FALSE
            tag(gt2[[i]], "idx") <- i
            addhandlerclicked(gt2[[i]], handler = function(h, ...) {
                gt.tmp = h$obj
                prev.idx = tag(gt.tmp, "prev.idx")
                gt.tmp.svalue = paste(svalue(gt.tmp))
                if (length(prev.idx) == 1 && tag(gt.tmp, "toggle") && 
                length(gt.tmp.svalue)>0 && 
                gt.tmp.svalue!=paste(gt.tmp[prev.idx, 2])) {
                    tmp = gt.tmp[prev.idx, 2]
                    gt.tmp[prev.idx, 2] = svalue(gt.tmp)
                    gt.tmp[svalue(gt.tmp, index = TRUE), 2] = tmp
                    mergegui_env$idx <- mergegui_env$idx + 1
                    mergegui_env$hstry1[[mergegui_env$idx]] <- mergegui_env$hstry1[[mergegui_env$idx - 1]]
                    mergegui_env$hstry1[[mergegui_env$idx]][, tag(gt.tmp, "idx")] <- gt.tmp[,2]
                    if (tag(gt.tmp, "idx") == 1) {
                        tmpgt4 = mergegui_env$gt4[mergegui_env$gt4[,1]==gt2[[i]][prev.idx,1], 2:3]
                        mergegui_env$gt4[mergegui_env$gt4[,1]==gt2[[i]][prev.idx,1], 2:3] = mergegui_env$gt4[mergegui_env$gt4[,1]==gt2[[i]][svalue(gt.tmp, index = TRUE),1], 2:3]
                        mergegui_env$gt4[mergegui_env$gt4[,1]==gt2[[i]][svalue(gt.tmp, index = TRUE),1], 2:3] = tmpgt4
                    }
                    mergegui_env$gt5[, 2] <- mergegui_env$gt4[order(mergegui_env$gt4[,1]), 2]
                    mergegui_env$gt5[, 3] <- mergegui_env$gt4[order(mergegui_env$gt4[,1]), 3]
                    mergegui_env$hstry2[[mergegui_env$idx]] <- mergegui_env$gt4[,]
                    mergegui_env$hstry3[[mergegui_env$idx]] <- mergegui_env$gt5[,]
                    if (length(svalue(check141))==3) {mergegui_env$hstry4[[mergegui_env$idx]] <- mergegui_env$idx} else {mergegui_env$hstry4[[mergegui_env$idx]] <- mergegui_env$hstry4[[mergegui_env$idx-1]]}
                    mergegui_env$redo.indicate <- 0
                }
                tag(gt.tmp, "toggle") = !tag(gt.tmp, "toggle")
                tag(gt.tmp, "prev.idx") = svalue(gt.tmp, index = TRUE)
            })
        }
        group23 <- ggroup(container = group21)
        gbcombo21 <- gbutton("Undo", container = group23, handler = undo,
                             expand = TRUE)
        gbcombo22 <- gbutton("Redo", container = group23, handler = redo,
                             expand = TRUE)
        gbcombo23 <- gbutton("Reset", container = group23, handler = reset,
                             expand = TRUE)
        
        #####------------------------------------------------#####
        ##  In the third tab we can:                            ##
        ##  (1) Watch and change the name or type of variables. ##
        ##  (2) Numeric or graphic summary.                     ##
        ##  (3) Dictionary for factor variables.                ##
        #####------------------------------------------------#####
        group41 = ggroup(container = tab, label = "Summary", expand = T)
        mergegui_env$group42 <- ggroup(container = group41, use.scrollwindow = TRUE, expand = T)
        
        mergegui_env$gt4 <- gtable(mergegui_env$name_intersection_panel, multiple = T, container = mergegui_env$group42, expand = TRUE, chosencol = 2)
        addhandlerdoubleclick(mergegui_env$gt4, handler = VariableOptions)
        mergegui_env$group43 <- ggroup(horizontal = FALSE, container = group41,
                                       expand = TRUE)
        group44 = ggroup(horizontal = TRUE, container = mergegui_env$group43)
        gbcombo431 <- gbutton("Numeric summary", container = group44,
                              handler = smmry, expand = TRUE)
        gbcombo432 <- gbutton("Graphical summary", container = group44,
                              handler = graph, expand = TRUE)
        gbcombo433 <- gbutton("Dictionary", container = group44, handler = dict,
                              expand = TRUE)
        mergegui_env$group45 <- ggroup(container = mergegui_env$group43, expand = TRUE, use.scrollwindow = TRUE)
        
        #####------------------------------------------------#####
        ##  In the fourth tab we can:                           ##
        ##  (1) Select all or none variables.                   ##
        ##  (2) Export the data.                                ##
        #####------------------------------------------------#####
        group51 = ggroup(container = tab, label = "Export", expand = T)
        group52 = ggroup(container = group51, use.scrollwindow = TRUE,
                         expand = T)
        mergegui_env$gt5 <- gtable(data.frame(mergegui_env$name_intersection_panel[,1:3],FileMatched), multiple = T, container = group52,
                                   expand = TRUE, chosencol = 2)
        addhandlerclicked(mergegui_env$gt5,handler=function(h,...){
            svalue(gbcombo57) = paste("Currently you select",length(svalue(mergegui_env$gt5)),"variables.",sep=" ")
        })
        group53 = ggroup(horizontal = FALSE, container = group51,
                         expand = TRUE)
        gbcombo51 <- gbutton("Select All", container = group53, handler = function(h,
                                                                                   ...) {
            svalue(mergegui_env$gt5, index = TRUE) = 1:length(nameintersect)
            svalue(gbcombo57) = paste("Currently you select all",length(nameintersect),"variables.",sep=" ")
            focus(mergegui_env$gt5)
        })
        gbcombo52 <- gbutton("Clear All", container = group53, handler = function(h,
                                                                                  ...) {
            svalue(mergegui_env$gt5) = NULL
            svalue(gbcombo57) = "Currently you select 0 variable."
        })
        gbcombo55 <- gbutton("Export the matched data", container = group53,
                             handler = watchdatafunc)
        gbcombo56 <- glabel(paste("The complete merged data have ",sum(rows)," rows and ",
                                  length(nameintersect)," columns."),container=group53)
        gbcombo57 <- glabel(paste("Currently you select 0 variable."),container=group53)
        
        svalue(tab)=1
    }
    
    mergeID = function(h, ...) {
        ####################################################
        # mergeID is a function to merge the observations. #
        ####################################################
        
        watchIDfunc = function(h, ...) {
            #####------------------------------------------------------------------------------------#####
            ##  watchIDfunc is a function to export the merged dataset.                                 ##
            ##  key is a vector of the selected primary keys. We checked the validity of the key first. ##
            ##  keyID is the merged ID for the observations.                                            ##
            ##  mergeIDdata is the matrix that merged all the files by the keyID.                       ##
            ##  We should write 'xxx.csv' when we export mergeIDdata and save the file.                 ##
            #####------------------------------------------------------------------------------------#####
            keyID = c()
            vcolumn = rep(0, n)
            key = c()
            for (i in 1:n) {
                key[i] = svalue(gt3[[i]])
                if (sum(duplicated(as.character(dataset[[i]][, key[i]]))) >
                    0) {
                    gmessage(paste(key[i], "could not be the primary key for",
                                   basename(gtfile[i]), "because it has repeated items. Please choose another key."))
                    return()
                }
                keyID = union(keyID, dataset[[i]][, key[i]])
                vcolumn[i] = length(vname[[i]]) - 1
            }
            mergeIDdata = matrix(NA, nrow = length(keyID), ncol = sum(vcolumn) +
                1, dimnames = list(keyID))
            mergeIDdata[, 1] = keyID
            mergeIDcolnames = c(key[1], 1:sum(vcolumn) + 1)
            for (i in 1:n) {
                mergeIDdata[as.character(dataset[[i]][, key[i]]),
                            cumsum(c(2, vcolumn))[i]:cumsum(c(1, vcolumn))[i +
                                1]] = as.matrix(dataset[[i]][, setdiff(vname[[i]],
                                                                       key[i])])
                mergeIDcolnames[cumsum(c(2, vcolumn))[i]:cumsum(c(1,
                                                                  vcolumn))[i + 1]] = paste(basename(gtfile[i]),
                                                                                            ".", colnames(dataset[[i]][, setdiff(vname[[i]],
                                                                                                                                 key[i])]), sep = "")
            }
            colnames(mergeIDdata) = mergeIDcolnames
            
            if (!is.na(gf <- gfile(type = "save"))) {
                if (regexpr("\\.csv$",gf) %in% c(-1,1)) {
                    write.csv(mergeIDdata, file = paste(gf,".csv",sep=""), row.names = FALSE)
                } else {
                    write.csv(mergeIDdata, file = gf, row.names = FALSE)
                }
                gmessage("The files are merged!")
            }
        }
        
        dataset <- list()
        vname <- list()
        simplifiedname <- list()
        if (length(svalue(gt)) == 0) {
            n <- length(gt[])
            gtfile <- gt[]
        }
        else {
            n <- length(svalue(gt))
            gtfile <- svalue(gt)
        }
        rows <- rep(0, n)
        for (i in 1:n) {
            dataset[[i]] <- if (length(grep("\\.csv$",gtfile[i]))) {
                read.csv(file = gtfile[i], header = T)
            } else { gtdata[[i]] }
            rows[i] <- nrow(dataset[[i]])
            vname[[i]] <- colnames(dataset[[i]])
            J = gregexpr(".csv.", vname[[i]])
            loc = c()
            for (j in 1:length(vname[[i]])) {
                loc[j] = max(J[[j]]) + 5
            }
            loc[which(loc == 4)] = 1
            simplifiedname[[i]] = substring(vname[[i]], loc)
        }
        
        a = intersect2(vname, simplifiedname)
        tmpuniq = a$uniq
        tmpnametable = a$individual
        tmpvname = list()
        for (i in 1:n) {
            tmpvname[[i]] = c(tmpnametable[, i], tmpuniq[[i]])
        }
        
        #####-----------------------------------------#####
        ##  In this GUI we can:                          ##
        ##  Select the primary keys for different files. ##
        #####-----------------------------------------#####
        combo3 <- gwindow("Matched Primary Key", visible = TRUE)
        group3 <- ggroup(horizontal = FALSE, container = combo3)
        gt3 <- list()
        
        for (i in 1:n) {
            gl3 = glabel(paste("Select the Primary Key from", basename(gtfile[i])),
                         container = group3)
            gt3[[i]] <- gcombobox(tmpvname[[i]], container = group3,
                                  expand = T)
        }
        gbcombo3 = gbutton("Match by the Key", container = group3,
                           handler = watchIDfunc)
        
    }
    
    #####---------------------#####
    ##  First GUI:  Open files.  ##
    #####---------------------#####
    gtdata=list(...)
    mycall <- as.list(match.call()[-1])
    if (!missing(filenames)) mycall <- mycall[names(mycall)!='filenames']
    if (!missing(unit)) mycall <- mycall[names(mycall)!='unit']
    if (!missing(distn)) mycall <- mycall[names(mycall)!='distn']
    if (!missing(miss)) mycall <- mycall[names(mycall)!='miss']
    datasets = as.character(unlist(mycall))
    
    combo <- gwindow("Combination", visible = TRUE)
    group <- ggroup(horizontal = FALSE, container = combo)
    if (is.null(filenames) & is.null(datasets)) {
        f.list <- matrix(nrow = 0, ncol = 1, dimnames = list(NULL, "File"))
    } else {
        f.list <- matrix(c(datasets,filenames), ncol = 1, dimnames = list(NULL, "File"))
    }
    gt <- gtable(f.list, multiple = T, container = group, expand = TRUE)
    gb1 <- gbutton("Open", container = group, handler = function(h, ...) gt[,] = union(gt[,],na.omit(gfile(multiple=TRUE))))
    gb2 <- gbutton("Match the Variables", container = group,
                   handler = mergefunc)
    gb3 <- gbutton("Match by the Key", container = group,
                   handler = mergeID)
}
