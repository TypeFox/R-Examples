#---------------------------------------------------------------------------
#
#   Methods for coercion from some "Stem" class object to other forms.
#
#   1. from=downLogs, to=data.frame
#   2. from=standingTrees, to=data.frame (26-Oct-2011)  
#
#
#Author...									Date: 21-Oct-2010
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#
#   this method just converts from a 'downLogs' object to a data frame in
#   the form of the result returned from sampleLogs()...
#
setAs('downLogs', 'data.frame',
      function(from) {
        dlogs = from
        numLogs = length(dlogs@logs)
        nCols = length(.StemEnv$sampleLogsNames)
        df = data.frame(matrix(NA, nrow=numLogs, ncol=nCols))
        colnames(df) = .StemEnv$sampleLogsNames

        for(i in seq_len(numLogs)) {
          df[i,'species'] = dlogs@logs[[i]]@species
          df[i,'logLen'] = dlogs@logs[[i]]@logLen
          if(dlogs@units == .StemEnv$msrUnits$metric)
            cf = .StemEnv$m2cm
          else
            cf = .StemEnv$ft2in
          df[i,'buttDiam'] = dlogs@logs[[i]]@buttDiam * cf
          df[i,'topDiam'] = dlogs@logs[[i]]@topDiam * cf
          df[i,'solidType'] = dlogs@logs[[i]]@solidType
          df[i,'x'] = dlogs@logs[[i]]@location@coords[,'x']
          df[i,'y'] = dlogs@logs[[i]]@location@coords[,'y']
          df[i,'logAngle'] = dlogs@logs[[i]]@logAngle
          df[i,'logAngle.D'] = dlogs@logs[[i]]@logAngle * .StemEnv$rad2deg
        }

        return(df)
      } #function
)   #setAs




#
#   this method just converts from a 'standingTrees' object to a data frame in
#   the form of the result returned from sampleTrees()...
#
setAs('standingTrees', 'data.frame',
      function(from) {
        strees = from
        numTrees = length(strees@trees)
        nCols = length(.StemEnv$sampleTreesNames)
        df = data.frame(matrix(NA, nrow=numTrees, ncol=nCols))
        colnames(df) = .StemEnv$sampleTreesNames

        for(i in seq_len(numTrees)) {
          df[i,'species'] = strees@trees[[i]]@species
          df[i,'height'] = strees@trees[[i]]@height
          if(strees@units == .StemEnv$msrUnits$metric)
            cf = .StemEnv$m2cm
          else
            cf = .StemEnv$ft2in
          df[i,'dbh'] = strees@trees[[i]]@dbh * cf
          df[i,'topDiam'] = strees@trees[[i]]@topDiam * cf
          df[i,'solidType'] = strees@trees[[i]]@solidType
          df[i,'x'] = strees@trees[[i]]@location@coords[,'x']
          df[i,'y'] = strees@trees[[i]]@location@coords[,'y']
        }

        return(df)
      } #function
)   #setAs 'standingTrees'

