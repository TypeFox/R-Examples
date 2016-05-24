#---------------------------------------------------------------------------
#
#   Methods for generic summary() for InclusionZoneGrid class...
#     (1)  inclusionZoneGrid class
#
#Author...									Date: 23-sept-2010
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#



#================================================================================
#  method for class InclusionZoneGrid...
#
setMethod('summary',
          signature(object = 'InclusionZoneGrid'),
function(object,
         ...
        )
{
#------------------------------------------------------------------------------
    cat('\nObject of class:', class(object))
    .StemEnv$underLine(60)
    if(!is.na(object@description))
      cat(object@description, fill=60)
    .StemEnv$underLine(60, prologue='')

    cat('\nInclusionZone class:', class(object@iz))
    cat('\n  units of measurement: ', object@iz@units)    

    grid = object@grid
    cat('\n\nGrid class:', class(grid))
    cat('\nNumber of grid cells =', ncell(grid))
    cat('\nCell dimensions: (nrows=',nrow(grid), ', ncol=',ncol(grid),')', sep='')
    cat('\nGrid cell values**...\n')
    gridValues = getValues(grid)
    gvt = data.frame(table(gridValues, useNA='ifany'))
    print(gvt)
    cat('**Note: data slot values get swapped with zero-valued grid cells as necessary.\n')

    
    cat('\nPer unit area estimates in the data slot (for cells inside IZ only)...\n')
    idx = ifelse(is.na(gridValues), FALSE, TRUE)    #mask out background cells
    df = object@data[idx,]
    #NAs in biomass and carbon will create extra entries in summary, so delete them as necessary...
    jj = apply(df, 2, function(x) all(is.na(x))) #also remove any other NA column with no estimates
    naNames = names(jj)[jj]
    if(length(naNames) > 0)                            #note, we could end up with just one column
      df = df[,-match(naNames, names(jj)), drop=FALSE] #which will convert to vector by default!!

    print(apply(df, 2, summary), digits=4)
    if(!any(is.na(match(c('biomass','carbon'), naNames)))) {
      cat('--Note: either biomass or carbon (or both) had all NAs because no conversion')
      cat('\n        factor was supplied, these columns have been deleted above.')
    }
    if(length(setdiff(naNames,c('biomass','carbon'))) > 0)
      cat('\n        Other columns with all NAs have also been removed.')

    cat('\n\n  Encapulating bounding box...\n')
    print(object@bbox)
    
    cat('\n')
    
    return(invisible())
}   #summary for 'InclusionZoneGrid'
) #setMethod


