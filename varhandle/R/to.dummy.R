## Function Description:
##     This function turns a categorical character column into multiple binary
##     columns. eg: change the gender (male/female) into male column and female
##     column.

to.dummy <- function(v=NULL, prefix=NULL){
    #----[ checking the input ]----#
    if(is.null(prefix)){
        stop("The input \"prefix\" is missing. This will be added to the begining of each column name to avoid conflicts with other column names.")
    }else if(length(prefix)!=1 | nchar(prefix)==0 | class(prefix)!="character"){
        stop("The input \"prefix\" should be a character vector with length of 1 and character number more than 0.")
    }
    
    if(is.null(v)){
        stop("The input \"v\" is missing. It should be a vector with categories that you are willing to create dummy variables from. It can be a factor, character or numeric vector.")
    }
    
    
    #----[ pre-processing ]----#
    ## convert to character vector if the input is factor
    if(class(v)=="factor"){
        v_levels <- levels(v)
        # convert the factor to character vector
        v <- as.character(v)
    }else{
        # find tha NAs and turn them into character to have a separate column for NAs
        v[which(is.na(v))] <- "NA"
        # get the categories
        v_levels <- names(table(v, useNA = "ifany"))
    }
    
    
    #----[ processing ]----#
    # go through categories one by one
    for(i in 1:length(v_levels)){
        # create a logical vector which has 1 for places that has the category
        assign(x=paste("v", i, sep=""), value=as.numeric(v==v_levels[i]))
    }


    # create a cbind command and run it. It attaches the variables generated in the for loop above.
    df <- eval(parse(text=paste("cbind(", paste('v', 1:i, sep='', collapse = ", "), ")", collapse="", sep="")))
    # strip the white space from begining and end of the name and the middle white space with "_"
    factor_levels <- gsub("\\s+", "_", gsub("^\\s+|\\s+$", "", v_levels))
    # if one of the levels are "", we should create a name for it, so we use "BLANK"
    factor_levels[which(factor_levels=="")] <- "BLANK"
    # set the colnames
    colnames(df) <- paste(prefix, factor_levels, sep=".")
    return(df)
}