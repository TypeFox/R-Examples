## ----load_package--------------------------------------------------------

## Install package from source, then load package into workspace
install.packages( type = "source", pkgs = "spThin_0.1.0.tar.gz", repos = NULL )
library( spThin )


## ------------------------------------------------------------------------
data( Heteromys_anomalus_South_America )
head( Heteromys_anomalus_South_America )

## ------------------------------------------------------------------------
table( Heteromys_anomalus_South_America$REGION )

## ------------------------------------------------------------------------
thinned_dataset_full <-
  thin( loc.data = Heteromys_anomalus_South_America, 
        lat.col = "LAT", long.col = "LONG", 
        spec.col = "SPEC", 
        thin.par = 10, reps = 100, 
        locs.thinned.list.return = TRUE, 
        write.files = TRUE, 
        max.files = 5, 
        out.dir = "hanomalus_thinned_full/", out.base = "hanomalus_thinned", 
        write.log.file = TRUE,
        log.file = "hanomalus_thinned_full_log_file.txt" )

## ------------------------------------------------------------------------
plotThin( thinned_dataset_full )

## ------------------------------------------------------------------------
thinned_dataset_mainland <-
  thin( loc.data = Heteromys_anomalus_South_America[ which( Heteromys_anomalus_South_America$REGION == "mainland" ) , ], 
        lat.col = "LAT", long.col = "LONG", 
        spec.col = "SPEC", 
        thin.par = 10, reps = 100, 
        locs.thinned.list.return = TRUE, 
        write.files = TRUE, 
        max.files = 5, 
        out.dir = "hanomalus_thinned_mainland/", out.base = "hanomalus_thinned", 
        write.log.file = TRUE,
        log.file = "hanomalus_thinned_mainland_log_file.txt" )

## ------------------------------------------------------------------------
thinned_dataset_trin <-
  thin( loc.data = Heteromys_anomalus_South_America[ which( Heteromys_anomalus_South_America$REGION == "trin" ) , ], 
        lat.col = "LAT", long.col = "LONG", 
        spec.col = "SPEC", 
        thin.par = 10, reps = 10, 
        locs.thinned.list.return = TRUE, 
        write.files = TRUE, 
        max.files = 5, 
        out.dir = "hanomalus_thinned_trin/", out.base = "hanomalus_thinned", 
        write.log.file = TRUE,
        log.file = "hanomalus_thinned_trin_log_file.txt" )

## ------------------------------------------------------------------------
thinned_dataset_mar <-
  thin( loc.data = Heteromys_anomalus_South_America[ which( Heteromys_anomalus_South_America$REGION == "mar" ) , ], 
        lat.col = "LAT", long.col = "LONG", 
        spec.col = "SPEC", 
        thin.par = 10, reps = 10, 
        locs.thinned.list.return = TRUE, 
        write.files = TRUE, 
        max.files = 5, 
        out.dir = "hanomalus_thinned_mar/", out.base = "hanomalus_thinned", 
        write.log.file = TRUE,
        log.file = "hanomalus_thinned_mar_log_file.txt" )

## ------------------------------------------------------------------------
thinned_dataset_tobago <-
  thin( loc.data = Heteromys_anomalus_South_America[ which( Heteromys_anomalus_South_America$REGION == "tobago" ) , ], 
        lat.col = "LAT", long.col = "LONG", 
        spec.col = "SPEC", 
        thin.par = 10, reps = 10, 
        locs.thinned.list.return = TRUE, 
        write.files = TRUE, 
        max.files = 5, 
        out.dir = "hanomalus_thinned_tobago/", out.base = "hanomalus_thinned", 
        write.log.file = TRUE,
        log.file = "hanomalus_thinned_tobago_log_file.txt" )

