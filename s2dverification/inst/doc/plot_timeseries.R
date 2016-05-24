#!/usr/bin/env Rscript

library(s2dverification)
library(ncdf4)
args <- commandArgs(TRUE)

comptrend <- TRUE
compcor <- TRUE
comprms <- TRUE
compspread <- TRUE
plotano <- TRUE

var <- args[1]    # sie/sia/siv/tos/tas/prlr/ohc/lohc/mohc/uohc/amoc
pole <- args[2]   # N/S only for sia/sie
nltimemax <- 124  # number of leadtimes max in the experiments (in months)
nltimeout <- 60   # number of leadtimes to postprocess(in months)
lstexpid <- c('i00k', 'b02p') # list of ids
grid <- '320x160' # atmospheric grid for tos/prlr (ocean/land only)
vertnem <- 'L42'  # Number of vertical levels in nemo
mon0 <- 11        # initial month
year0 <- 1960     # first start date
yearf <- 2005     # last start date
intsdate <- 5     # interval between start dates
runmeanlen <- 12  # length of the window for running mean (in months)

obs <- switch(var, 'sia' = c('HadISST'), 'sie' = c('HadISST'), 
              'tas' = c('NCEP', 'ERA40'), 'tos' = c('ERSST', 'HadISST'),
              'prlr' = c('CRU', 'GPCC'), 'ohc' = c('NEMOVAR_S4'),
              'mohc' = c('NEMOVAR_S4'), 'uohc' = c('NEMOVAR_S4'),
              'lohc' = c('NEMOVAR_S4'), 'amoc' = c('NEMOVAR_S4'),
              'siv' = 'PIOMAS')
toptitle2 <- switch(var, 'sia' = "sea ice area", 'sie' = "sea ice extent",
                    'siv' = "sea ice volume", 'tas' = "global T2m", 
                    'tos' = "global SST (60S-65N)", 
                    'prlr' = 'land precipitation (60S-65N)', 
                    'ohc' = "global ocean heat content", 
                    'lohc' = 'global 800m-bottom ocean heat content',
                    'mohc' = 'global 350m-800m ocean heat content',
                    'uohc' = 'global 0-350m ocean heat content',
                    'amoc' = 'Atlantic Overturning Streamfunction (40-55N, 1-2km)'
                    )
ytitle1 <- switch(var, 'sia' = "Millions km2", 'sie' = "Millions km2", 
                  'siv' = 'Thousands km3', 'tas' = 'K', 'tos' = 'K', 
                  'prlr' = 'mm/day', 'ohc' = '10e22 Joules', 
                  'lohc' = '10e22 Joules', 'mohc' = '10e22 Joules',
                  'uohc' = '10e22 Joules', 'amoc' = 'Sv')

syears <- seq(year0, yearf, intsdate)
imon2 <- paste("0", as.character(mon0), sep = "")
sdates <- paste(as.character(syears), substr(imon2, nchar(imon2) - 1, 
                nchar(imon2)), '01', sep = "")
toptitle1 <- paste(switch(pole, 'N' = "Arctic", 'S' = "Antarctic", ""),
                   toptitle2)

savename <- paste(var, switch(pole, 'N' = paste('_', pole, sep = ''),
                  'S' = paste('_', pole, sep = ''), ''), sep = '')
for (expid in lstexpid ) {
  savename <- paste(savename, '_', expid, sep = '')
}
if (file.exists(paste(savename, '.sav', sep = ''))) {
  load(paste(savename, '.sav', sep = ''))
} else {
  if (var == 'prlr' | var == 'tos' ) {
    fnc <- nc_open(paste('/esnas/exp/ecearth/land_sea_mask_', grid, '.nc',
                     sep = ''))
    mask <- ncvar_get(fnc, 'LSM')
    nc_close(fnc)
    if (var == 'prlr') {
      fnc <- nc_open('/esnas/obs/dwd/gpcc_combined1x1_v6/constant_fields/land_sea_mask.nc'
                       )
      mask_gpcc <- ncvar_get(fnc, 'lsm')
      nc_close(fnc)
      fnc <- nc_open('/esnas/obs/cru/mask_cru_land.nc')
      mask_cru <- ncvar_get(fnc, 'pre')
      nc_close(fnc)
      #fnc <- nc_open('/esnas/obs/noaa/gpcp_v2.2/constant_fields/land_sea_mask.nc'
      #                 )
      #mask_gpcp <- ncvar_get(fnc, 'LSM')
      #nc_close(fnc)
      lstmaskobs <- list(mask_cru, mask_gpcc)
    } else {
      mask <- 1 - mask
      lstmaskobs <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 
                         NULL)
    }
    lstmask <- list()
    for (iexp in 1:length(lstexpid)) {
       lstmask[[iexp]] <- mask
    }
  } else {
    lstmask <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, 
                    NULL, NULL, NULL, NULL, NULL) 
    lstmaskobs <- list(NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL, NULL,
                       NULL)
  }
  latbnd <- switch(var, 'tos' = c(-60, 65), 'prlr' = c(-60, 65), c(-90, 90))
  varname <- switch(var, 'sia' = paste(var, pole, sep = ''), 'sie' = paste(var, 
             pole, sep = ''), 'siv' = paste(var, pole, sep = ''), 
             'ohc' = 'heatc', 'uohc' = switch(vertnem, 'L42' = '0-315_heatc', 
                                                       'L46' = '0-322_heatc',
                                                       'L75' = '0-271_heatc'),
                              'mohc' = switch(vertnem, 'L42' = '373-657_heatc',
                                                       'L46' = '382-735_heatc',
                                                       'L75' = '301-773_heatc'),
                              'lohc' = switch(vertnem, 'L42' = '800-5350_heatc',
                                                       'L46' = '855-5875_heatc',
                                                       'L75' = '857-5902_heatc'),
              'amoc' = 'moc_40N55N_1-2km', var)
  if (is.na(match('b02p', lstexpid)) == TRUE) { 
    lstload <- lstexpid
  } else {
    lstload <- lstexpid[-match('b02p', lstexpid)]
  }       
  toto1 <- Load(varname, lstload, obs, sdates, latmin = latbnd[1],
                latmax = latbnd[2], nleadtime = nltimemax, 
                leadtimemax = nltimeout, maskmod = lstmask, 
                maskobs = lstmaskobs)
  if (is.na(match('b02p', lstexpid)) == FALSE) {
    toto1bis <- Load(varname, 'b02p', obs = NULL, '19501101', 
                     latmin = latbnd[1], latmax = latbnd[2], maskmod = lstmask,
                     maskobs = lstmaskobs)
    toto1ter <- Histo2Hindcast(toto1bis$mod, '19501101', sdates, 
                               nleadtimesout = nltimeout)
    toto1beta <- array(dim = c(dim(toto1$mod)[1] + dim(toto1ter)[1], 
                       max(dim(toto1$mod)[2], dim(toto1ter)[2]),
                       dim(toto1$mod)[3:4]))
    toto1beta[1:dim(toto1$mod)[1], 1:dim(toto1$mod)[2], , ] <- toto1$mod
    toto1beta[(dim(toto1$mod)[1] + 1):(dim(toto1$mod)[1] + dim(toto1ter)[1]), 
              1:dim(toto1ter)[2], , ] <- toto1ter
    toto1$mod <- toto1beta
    lstexpid <- c(lstload, 'b02p')
  }
  if (var == 'prlr') {
    toto1$mod <- toto1$mod * 1000 * 3600 * 24
    toto1$obs <- toto1$obs * 1000 * 3600 * 24
  }
  if (var == 'ohc' | var == 'lohc' | var == 'mohc' | var == 'uohc') {
    toto1$mod <- toto1$mod / 1e22
    toto1$obs <- toto1$obs / 1e22
  }
  if (var == 'sia' | var=='sie' | var=='siv') {
    toto1$mod <- toto1$mod/1000
    if (var == 'siv') {
      toto1$obs <- toto1$obs/1000
    }
  }
  save(toto1, file = paste(savename, '.sav', sep = ''))
}

toto2a <- Clim(toto1$mod, toto1$obs, memb = TRUE)
toto2b_ano_exp <- Ano(toto1$mod, InsertDim(toto2a$clim_exp, 
                                           3, dim(toto1$mod)[3]) )
toto2b_ano_obs <- Ano(toto1$obs, InsertDim(toto2a$clim_obs,
                                           3, dim(toto1$obs)[3]) )
toto3 <- Smoothing(toto2b_ano_exp, runmeanlen, 4)
toto4 <- Smoothing(toto2b_ano_obs, runmeanlen, 4)
suf <- switch(pole, 'N' = paste('_', pole, sep = ''), 'S' = paste('_', pole,
              sep = ''), '')
PlotAno(toto1$mod, toto1$obs, sdates, toptitle = paste(lstexpid, toptitle1),
        ytitle = c(ytitle1, ytitle1, ytitle1), legends = obs, biglab = F, 
        fileout = paste(var, '_', lstexpid, suf, '.eps', sep = ''))
PlotAno(Smoothing(toto1$mod, runmeanlen, 4), 
        Smoothing(toto1$obs, runmeanlen, 4), sdates, 
        toptitle = paste("smoothed", lstexpid, toptitle1),
        ytitle = c(ytitle1, ytitle1, ytitle1), legends = obs, biglab = F, 
        fileout = paste(var, '_', lstexpid, suf, '_smoothed.eps', sep = ''))

if (plotano) {
  PlotAno(toto3, toto4, sdates, toptitle = paste("smoothed", lstexpid,
          toptitle1, "anomalies"), ytitle = c(ytitle1, ytitle1, ytitle1), 
          legends = obs, biglab = F, fileout = paste(var, '_', lstexpid,suf, 
          '_ano.eps', sep = ''))
  PlotClim(toto2a$clim_exp, toto2a$clim_obs, toptitle = paste(switch(pole,
           'N' = "Arctic", 'S' = "Antarctic", ""), toptitle2, "climatologies"),
           ytitle = ytitle1, monini = mon0, listexp = lstexpid, listobs = obs,
           biglab = F, fileout = paste(savename, '_clim.eps', sep = ''))
} 

if (compspread) {
  toto5 <- toto3 - InsertDim(Mean1Dim(toto3, 2, narm = T), 2, dim(toto3)[2])
  toto6 <- Spread(toto5, c(2, 3))
  PlotVsLTime(toto6$iqr, toptitle = paste("InterQuartile Range", toptitle1),
              ytitle = ytitle1, monini = mon0, listexp = lstexpid, biglab = F,
              fileout = paste("IQR_", savename, ".eps", sep = ''))
  PlotVsLTime(toto6$maxmin, toptitle = paste("Maximum-Minimum for", toptitle1),
              ytitle = ytitle1, monini = mon0, listexp = lstexpid, biglab = F, 
              fileout = paste("MaxMin_", savename, ".eps", sep = ''))
  PlotVsLTime(toto6$sd, toptitle = paste("Standard Deviation for", toptitle1),
              ytitle = ytitle1, monini = mon0, listexp = lstexpid, biglab = F,
              fileout = paste("SD_", savename, ".eps", sep = ''))
  PlotVsLTime(toto6$mad, toptitle = paste("Median Absolute Deviation for",
              toptitle1), ytitle = ytitle1, monini = mon0, listexp = lstexpid,
              biglab = F, fileout = paste("Mad_", savename, ".eps", sep = ''))
}

if (compcor) {
  cor <- Corr(Mean1Dim(toto3, 2), Mean1Dim(toto4, 2), 1, 2, compROW = 3, 
              limits = c(ceiling((runmeanlen + 1) / 2), 
              nltimeout - floor(runmeanlen / 2)))
  PlotVsLTime(cor, toptitle = paste("Correlations for", toptitle1), 
              ytitle = "correlation", monini = mon0, limits = c(-1, 2),
              listexp = lstexpid, listobs = obs, biglab = F, 
              hlines = c(-1, 0, 1), fileout = paste("cor_", savename, ".eps",
              sep = ''))
}

if (comprms) {
  rms <- RMS(Mean1Dim(toto3, 2), Mean1Dim(toto4, 2), 1, 2, compROW = 3, 
             limits = c(ceiling((runmeanlen + 1) / 2), 
             nltimeout - floor(runmeanlen / 2)))
  PlotVsLTime(rms, toptitle = paste("RMSE for", toptitle1), ytitle = ytitle1,
              monini = mon0, listexp = lstexpid, listobs = obs, biglab = F,
              fileout = paste("rms_", savename, ".eps", sep = ""))
}

if (comptrend) {
  trends <- Consist_Trend(Mean1Dim(toto3, 2), Mean1Dim(toto4, 2), intsdate / 12)
  PlotVsLTime(trends$trend, toptitle = paste("Trend for", toptitle1), 
              ytitle = paste(ytitle1, "/ year"), monini = mon0,
              listexp = c(lstexpid, obs), biglab = F, fileout = paste("trend_",
              savename, ".eps", sep = ""))
}

rm(list = ls())
quit()