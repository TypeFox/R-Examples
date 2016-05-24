#---------------------------------------------------------------------------
#
#   Methods for coercion from some "izContainer" class object to other forms.
#
#   1. from=downLogIZs, to=downLogs
#   2. from=standingTreeIZs, to=standingTrees  
#
#
#Author...									Date: 12-Dec-2011
#	Jeffrey H. Gove
#	USDA Forest Service
#	Northern Research Station
#	271 Mast Road
#	Durham, NH 03824
#	jhgove@unh.edu
#	phone: 603-868-7667	fax: 603-868-7604
#---------------------------------------------------------------------------
#
#   this method just converts from a 'downLogIZs' izContainer object to a
#   'downLogs' StemContainer object...
#
setAs('downLogIZs', 'downLogs',
      function(from) {
        izs = from@iZones
        numLogs = length(izs)
        logs = vector('list', numLogs)  
        for(i in seq_len(numLogs))
          logs[[i]] = izs[[i]]@downLog
        Logs = downLogs(logs)

        return(Logs)

      } #function
)   #setAs


#
#   this method just converts from a 'standingTreeIZs' izContainer object
#   to a 'standingTrees' StemContainer object...
#
setAs('standingTreeIZs', 'standingTrees',
      function(from) {
        izs = from@iZones
        numTrees = length(izs)
        trees = vector('list', numTrees)  
        for(i in seq_len(numTrees))
          trees[[i]] = izs[[i]]@standingTree
        Trees = standingTrees(trees)

        return(Trees)

      } #function
)   #setAs
