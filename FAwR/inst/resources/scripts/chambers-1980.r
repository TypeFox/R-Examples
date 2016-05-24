#	$Id: chambers-1980.r 3684 2010-03-17 19:38:33Z hamannj $	

## this file contains R interface functions to the chambers-1980.so shared library
## the library is a single set of implimentation (*.c) and header (*.h) file
## and is created by issuing the following commands into the bash prompt:
## $ R CMD SHLIB chambers_1980.c

## /****************************************************************************/
## /* Charles J. Chambers 1980. Emperical growth and yield tables for the      */
## /* Douglas-fir zone.  DNR Report No. 41.  Washington Department of Natural  */
## /* Resources. Olympia WA. 98504                                             */
## /****************************************************************************/

## /*
##   Emperical yield functions for Douglas-fir from:

##          Min  Max     Mean      SD
##   Age     9    99    38.05   17.35
##   SI     38   165   107.39   21.42 (Kings 1966)
##   TPA               178.06  111.48
##   BA                135.93   78.37
##   CFV              4463.1  3353.87
##   QMD                12.35    3.36
##   normBA              0.82    0.38

## */

## chambers.1980.nba <- function(site.index,total.age)
## chambers.1980.ntpa.dnr <- function(site.index,total.age,pnba)
## chambers.1980.ntpa <- function(site.index,total.age,pnba)
## chambers.1980.adbh.dnr <- function(site.index,total.age,pnba)
## chambers.1980.adbh <- function(site.index,total.age,pnba)
## chambers.1980.atarif <- function(site.index,total.age,pnba)
## chambers.1980.cubic.foot.volume <- function(site.index,total.age,pnba)
## chambers.1980.scribner.16.volume <- function(site.index,total.age,pnba)
## chambers.1980.scribner.32.volume <- function(site.index,total.age,pnba)
## chambers.1980.atdba <- function(dbh,ba)
## chambers.1980.bfcfr <- function(dbh,tarif)
## chambers.1980.smh <- function(site.index,total.age)
## chambers.1980.nbag <- function(total.age,basal.area,site.index)


chambers.1980.nba <- function(site.index,total.age) {
  ret.val <- .C("chambers_1980_normal_basal_area",
                as.double(site.index),
                as.double(total.age),
                pred.nba = as.double(0) )$pred.nba
  ret.val
}


chambers.1980.ntpa.dnr <- function(site.index,total.age,pnba) {
  ret.val <- .C("chambers_1980_normal_tpa_dnr",
                as.double(site.index),
                as.double(total.age),
                as.double(pnba),
                ntpa.dnr = as.double(0) )$ntpa.dnr
  ret.val
}

## void chambers_1980_normal_tpa(
##     double   *site_index,
##     double   *total_age,
##     double   *pnba,
##     double   *ntpa );

chambers.1980.ntpa <- function(site.index,total.age,pnba) {
  ret.val <- .C("chambers_1980_normal_tpa",
                as.double(site.index),
                as.double(total.age),
                as.double(pnba),
                ntpa = as.double(0) )$ntpa
  ret.val
}

## void chambers_1980_avg_dbh_dnr(
##     double *site_index,
##     double *total_age,
##     double *pnba,
##     double *avg_dbh_dnr );

chambers.1980.adbh.dnr <- function(site.index,total.age,pnba) {
  ret.val <- .C("chambers_1980_avg_dbh_dnr",
                as.double(site.index),
                as.double(total.age),
                as.double(pnba),
                adbh.dnr = as.double(0) )$adbh.dnr
  ret.val
}


## /* double chambers_1980_avg_dbh( */
## /*     double *site_index, */
## /*     double *total_age, */
## /*     double *pnba                  ); */

## void chambers_1980_avg_dbh(
##     double *site_index,
##     double *total_age,
##     double *pnba,
##     double *avg_dbh );

chambers.1980.adbh <- function(site.index,total.age,pnba) {
  ret.val <- .C("chambers_1980_avg_dbh",
                as.double(site.index),
                as.double(total.age),
                as.double(pnba),
                adbh = as.double(0) )$adbh
  ret.val
}

## /* double chambers_1980_average_stand_tarif( */
## /*     double *site_index,  */
## /*     double *total_age,  */
## /*     double *pnba                  ); */

## void chambers_1980_average_stand_tarif(
##     double *site_index,
##     double *total_age,
##     double *pnba,
##     double *avg_stand_tarif );

chambers.1980.atarif <- function(site.index,total.age,pnba) {
  ret.val <- .C("chambers_1980_average_stand_tarif",
                as.double(site.index),
                as.double(total.age),
                as.double(pnba),
                atarif = as.double(0) )$atarif
  ret.val
}


## /* double chambers_1980_cubic_foot_volume( */
## /*     double *site_index,  */
## /*     double *total_age,  */
## /*     double *pnba                  ); */

## void chambers_1980_cubic_foot_volume(
##     double *site_index,
##     double *total_age,
##     double *pnba,
##     double *cubic_foot_volume );


chambers.1980.cubic.foot.volume <- function(site.index,total.age,pnba) {
  ret.val <- .C("chambers_1980_cubic_foot_volume",
                as.double(site.index),
                as.double(total.age),
                as.double(pnba),
                cfv = as.double(0) )$cfv
##  print( ret.val )
  ret.val
}

##chambers.1980.cubic.foot.volume( 120.0, 50.0, 1.0 )




## /* these next two functions still don't check out */
## /* double chambers_1980_scribner_16_volume( */
## /*     double *site_index,  */
## /*     double *total_age,  */
## /*     double *pnba                  );  */

## void chambers_1980_scribner_16_volume(
##     double *site_index,
##     double *total_age,
##     double *pnba,
##     double *scibner_16_volume );

chambers.1980.scribner.16.volume <- function(site.index,total.age,pnba) {
  ret.val <- .C("chambers_1980_scribner_16_volume",
                as.double(site.index),
                as.double(total.age),
                as.double(pnba),
                s16vol = as.double(0) )$s16vol
  ret.val
}


## /* double chambers_1980_scribner_32_volume( */
## /*     double *site_index,  */
## /*     double *total_age,  */
## /*     double *pnba                  ); */

## void chambers_1980_scribner_32_volume(
##    double *site_index,
##    double *total_age,
##    double *pnba,
##    double *scribner_32_volume );

chambers.1980.scribner.32.volume <- function(site.index,total.age,pnba) {
  ret.val <- .C("chambers_1980_scribner_32_volume",
                as.double(site.index),
                as.double(total.age),
                as.double(pnba),
                s32vol = as.double(0) )$s32vol
  ret.val
}

## /* double chambers_1980_cubic_volume_dbh_basal_area( */
## /*     double   *dbh, */
## /*     double   *basal_area          ); */

## void chambers_1980_cubic_volume_dbh_basal_area(
##     double   *dbh,
##     double   *basal_area,
##     double *cubic_volume_dbh_basal_area);

chambers.1980.scribner.32.volume <- function(site.index,total.age,pnba) {
  ret.val <- .C("chambers_1980_scribner_32_volume",
                as.double(site.index),
                as.double(total.age),
                as.double(pnba),
                s32vol = as.double(0) )$s32vol
  ret.val
}

    

## /* double chambers_1980_average_tarif_dbh_basal_area( */
## /*     double   *dbh, */
## /*     double   *basal_area          ); */

## void chambers_1980_average_tarif_dbh_basal_area(
##     double   *dbh,
##     double   *basal_area,
##     double *average_tarif_dbh_basal_area );

chambers.1980.atdba <- function(dbh,ba) {
  ret.val <- .C("chambers_1980_average_tarif_dbh_basal_area",
                as.double(dbh),
                as.double(ba),
                atdba = as.double(0) )$atdba
  ret.val
}





## /* double chambers_1980_board_foot_cubic_foot_ratio( */
## /*     double   *dbh, */
## /*     double   *tarif               ); */

## void chambers_1980_board_foot_cubic_foot_ratio(
##     double   *dbh,
##     double   *tarif,
##    double *board_foot_cubic_foot_ratio);

chambers.1980.bfcfr <- function(dbh,tarif) {
  ret.val <- .C("chambers_1980_board_foot_cubic_foot_ratio",
                as.double(dbh),
                as.double(tarif),
                bfcfr = as.double(0) )$bfcfr
  ret.val
}



## /* double chambers_1980_stand_mean_height( */
## /*     double   *site_index, */
## /*     double   *total_age          ); */

## void chambers_1980_stand_mean_height(
##    double   *site_index,
##    double   *total_age,
##    double *stand_mean_height );

chambers.1980.smh <- function(site.index,total.age) {
  ret.val <- .C("chambers_1980_stand_mean_height",
                as.double(site.index),
                as.double(total.age),
                smh = as.double(0) )$smh
  ret.val
}

## /* double chambers_1980_net_basal_area_growth( */
## /*    double   *total_age, */
## /*    double   *basal_area, */
## /*    double   *site_index ); */

## void chambers_1980_net_basal_area_growth(
##    double   *total_age,
##    double   *basal_area,
##    double   *site_index,
##    double *net_basal_area_growth );

chambers.1980.nbag <- function(total.age,basal.area,site.index) {
  ret.val <- .C("chambers_1980_stand_mean_height",
                as.double(total.age),
                as.double(basal.area),
                as.double(site.index),
                nbag = as.double(0) )$nbag
  ret.val
}


## this function generates a data.frame that contain a few metrics of interest.
chambers.1980 <- function( ages=1:100, site=125.0, pnba=1.0 ) {
  ret.val <- matrix( 0, length(ages), 5 )
  for( i in ages ) {    
    res <- c( i, 
             chambers.1980.adbh(site,i,pnba), ## table 5, page 6
             chambers.1980.smh( site,i ), ## table 16, page 10
             chambers.1980.nba(site,i ), ## table 1, page 4
             chambers.1980.ntpa( site,i,pnba ) ## table 3, page 5
             )
    ret.val[i,] <- res
  }
  ret.val <- as.data.frame( ret.val )
  names( ret.val ) <- c("age","qmd","tht","ba","expf" )
  ret.val
}
