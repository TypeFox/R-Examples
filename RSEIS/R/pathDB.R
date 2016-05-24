pathDB<-function(DB, path1="", path2="" )
  {
    
##############   fix  the data base by change the path locations

    
    ##  files move to a new location (different computer or server) ,
    ##  but retain the same basic file structure,
    ###  need only change the path names

    ###  for example:
 ###  "/data/wadati/soju/SEISMIC_DATA/Reventador2005/rev05/SEGY/R256.01/05.256.12.24.06.9024.6"
    ###  change to:
 ###  "/mnt/SEISMIC_DATA/Reventador2005/rev05/SEGY/R256.01/05.256.12.24.06.9024.6"
 ### i.e.  want to change:   path1="/data/wadati/soju" with path2="/mnt"

    
    fns = DB$fn

   newfns=  sub(path1, path2, fns, perl=TRUE)

    DB$fn = newfns
    invisible(DB)


  }
