#!/usr/bin/env Rscript

library(s2dverification)
args <- commandArgs(TRUE)

comptrend <- T    # Trend as a function of the start date for each leadtime
compcor <- T      # Correlation Coefficient
comprms <- T      # Root Mean Square Error
comprmsss <- T    # Root Mean Square Skill Score
compratrms <- T   # Ratio RMSE expid1 / expid2

var <- args[1]    # tos/tas/prlr
season <- args[2] # Year/DJF/MAM/JJA/SON
ly <- args[3]     # ly1/ly2-5/ly6-9 for Forecast year 1 / years 2 to 5 / years 
                  # 6 to 9
nltimemax <- 124  # number of leadtimes max in the experiments (in months)
lstexpid <- c('i00k','b02p') # list of ids
mon0 <- 11        # initial month
year0 <- 1960     # first start date
yearf <- 2005     # last start date
intsdate <- 5     # interval between start dates

obs <- switch(var, 'tas' = 'GHCNERSSTGISS', 'tos' = 'ERSST', 'prlr' = 'GPCC')
syears <- seq(year0, yearf, intsdate)
imon2 <- paste("0", as.character(mon0), sep = "")
sdates <- paste(as.character(syears), substr(imon2, nchar(imon2) - 1, 
                nchar(imon2)), '01', sep = "")

savename <- paste(var, '_', season, '_', ly, sep = '')
for (expid in lstexpid ) {
  savename <- paste(savename, '_', expid, sep = '')
}
savename <- paste(savename, '.sav', sep = '')

if (file.exists(savename)) {
  load(savename)
} else {
  if (is.na(match('b02p', lstexpid)) == TRUE) {
    lstload <- lstexpid
  } else {
    lstload <- lstexpid[-match('b02p', lstexpid)]
  }
  toto <- Load(var, lstload, obs,sdates, nleadtime = nltimemax,
               leadtimemin = switch(ly, 'ly1' = 1, 'ly2-5' = 13, 'ly6-9' = 61),
               leadtimemax = switch(ly, 'ly1' = switch(season, 'SON' = 13, 12),
               'ly2-5' = switch(season, 'SON' = 61, 60), 
               'ly6-9' = switch(season, 'SON' = 109, 108)), output = 'lonlat')
  if (is.na(match('b02p', lstexpid)) == FALSE) {
    toto1bis <- Load(var, 'b02p', obs = NULL, '19501101', output = 'lonlat')
    toto1ter <- Histo2Hindcast(toto1bis$mod, '19501101', paste(as.character(
                syears + switch(ly, 'ly1' = 0, 'ly2-5' = 1, 'ly6-9' = 5)),
                substr(imon2, nchar(imon2) - 1, nchar(imon2)), '01', sep = ""),
                nleadtimesout = switch(ly, 'ly1' = switch(season, 'SON' = 13,
                12), switch(season, 'SON' = 49, 48)))
    toto1beta <- array(dim = c(dim(toto$mod)[1] + dim(toto1ter)[1], 
                 max(dim(toto$mod)[2], dim(toto1ter)[2]), dim(toto$mod)[3:6]))
    toto1beta[1:dim(toto$mod)[1], 1:dim(toto$mod)[2], , , , ] <- toto$mod
    toto1beta[(dim(toto$mod)[1] + 1):(dim(toto$mod)[1] + dim(toto1ter)[1]),
              1:dim(toto1ter)[2], , , , ] <- toto1ter
    toto$mod <- toto1beta
    lstexpid <- c(lstload, 'b02p')
  }
  toto_exp <- InsertDim(Mean1Dim(Season(toto$mod, 4, mon0, switch(season,
                        'Year' = mon0, 'DJF' = 12, 'MAM' = 3, 'JJA' = 6,
                        'SON' = 9), switch(season, 
                        'Year' = (mon0 + 12 - 2) %% 12 + 1, 'DJF' = 2, 
                        'MAM' = 5, 'JJA' = 8, 'SON' = 11)), 4), 4, 1)
  toto_obs <- InsertDim(Mean1Dim(Season(toto$obs, 4, mon0, switch(season,
                        'Year' = mon0, 'DJF' = 12, 'MAM' = 3, 'JJA' = 6,
                        'SON' = 9), switch(season, 
                        'Year' = (mon0 + 12 - 2) %% 12 + 1, 'DJF' = 2, 
                        'MAM' = 5, 'JJA' = 8, 'SON' = 11)), 4), 4, 1)
  if (var == 'prlr') {
    toto$mod <- toto$mod * 1000 * 3600 * 24
    toto$obs <- toto$obs * 1000 * 3600 * 24
  }
  toto=list(mod=toto_exp,obs=toto_obs,lat=toto$lat,lon=toto$lon)
  save(toto,file=savename)
}

clims <- Clim(toto$mod, toto$obs)
ano_exp <- Ano(toto$mod, clims$clim_exp)
ano_obs <- Ano(toto$obs, clims$clim_obs)

if (compcor) {
  cor <- Corr(Mean1Dim(ano_exp, 2), Mean1Dim(ano_obs, 2), 1, 2)
  cols <- c("dodgerblue4", "dodgerblue1", "forestgreen", "yellowgreen",
            "white", "white", "yellow", "orange", "red", "saddlebrown")
  lims <- seq(-1, 1, 0.2)
  for (jexp in 1:length(lstexpid)) {
    flag <- array(F, dim = dim(cor[jexp, 1, 2, 1, , ]))
    flag[which(cor[jexp, 1, 2, 1, , ] > cor[jexp, 1, 4, 1, , ])] <- T
    postscript(paste('CorCoef2d_', var, '_', lstexpid[jexp], '_', season, '_',
               ly, '.eps', sep = ''))
    PlotEquiMap(cor[jexp, 1, 2, 1, , ], toto$lon, toto$lat, 
                toptitle = paste('Correlation Coefficient', lstexpid[jexp],
                switch(season, 'Year' = 'annual', season), switch(var, 
                'tas' = 'near surface temperature', 
                'tos' = 'sea surface temperature', 'prlr' = 'precipitation'),
                switch(var, 'tas' = 'GHCNv2+ERSSTv3b+GISTEMP', 
                'tas' = 'ERSSTv3b', 'prlr' = 'GPCC'), switch(ly, 
                'ly1' = 'Year1', 'ly2-5' = 'Year2-5', 'ly6-9' = 'Year6-9')),
                sizetit = 0.8, brks = lims, cols = cols, colNA = switch(var,
                'prlr' = 'white', grey(0.4)), filled.continents = switch(var,
                'tas' = F, 'tos' = T, 'prlr' = F), dots = t(flag), intylat = 45)
    dev.off()
  }
}

if (comprms) {
  rmse <- RMS(Mean1Dim(ano_exp, 2), Mean1Dim(ano_obs, 2), 1, 2)
  cols <- rev(c("dodgerblue4", "dodgerblue1", "forestgreen", "yellowgreen",
                "white", "white", "yellow", "orange", "red", "saddlebrown"))
  lims <- c(0, 0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.7, 1, 1.5, 2)
  lims <- switch(var, 'tas' = lims * 2, 'tos' = lims * 2, lims)
  rmse[which(rmse > max(lims))] <- max(lims)
  for (jexp in 1:length(lstexpid)) {
    postscript(paste('RMSE2d_', var, '_', lstexpid[jexp], '_', season, '_', ly,
               '.eps', sep = ''))
    PlotEquiMap(rmse[jexp, 1, 2, 1, , ], toto$lon, toto$lat, 
                toptitle = paste('RMSE', lstexpid[jexp], switch(season, 
                'Year' = 'annual', season), switch(var, 
                'tas' = 'near surface temperature', 
                'tos' = 'sea surface temperature', 'prlr' = 'precipitation'),
                switch(var, 'tas' = 'GHCNv2+ERSSTv3b+GISTEMP', 
                'tas' = 'ERSSTv3b', 'prlr' = 'GPCC'), switch(ly, 
                'ly1' = 'Year1', 'ly2-5' = 'Year2-5', 'ly6-9' = 'Year6-9')),
                sizetit = 0.8, brks = lims, cols = cols, colNA = switch(var, 
                'prlr' = 'white', grey(0.4)), filled.continents = switch(var,
                'tas' = F, 'tos' = T, 'prlr' = F), intylat = 45)
    dev.off()
  }
}

if (comprmsss) {
  rmsss <- RMSSS(Mean1Dim(ano_exp, 2), Mean1Dim(ano_obs, 2), 1, 2)
  cols <- c("dodgerblue4", "dodgerblue1", "forestgreen", "yellowgreen", 
            "white", "white", "yellow", "orange", "red", "saddlebrown")
  lims <- seq(-1, 1, 0.2)
  for (jexp in 1:length(lstexpid)) {
    flag <- array(F, dim = dim(rmsss[jexp, 1, 2, 1, , ]))
    flag[which(rmsss[jexp, 1, 2, 1, , ] < 0.05)] <- T
    rmsss[which(-1 > rmsss)] = -1
    postscript(paste('RMSSS2d_', var, '_', lstexpid[jexp], '_', season, '_', ly,
               '.eps', sep = ''))
    PlotEquiMap(rmsss[jexp, 1, 1, 1, , ], toto$lon, toto$lat, 
                toptitle = paste('RMSSS', lstexpid[jexp], switch(season,
                'Year' = 'annual', season), switch(var, 
                'tas' = 'near surface temperature', 
                'tos' = 'sea surface temperature', 'prlr' = 'precipitation'), 
                switch(var, 'tas' = 'GHCNv2+ERSSTv3b+GISTEMP', 
                'tas' = 'ERSSTv3b', 'prlr' = 'GPCC'), switch(ly, 
                'ly1' = 'Year1', 'ly2-5' = 'Year2-5', 'ly6-9' = 'Year6-9')),
                sizetit = 0.8, brks = lims, cols = cols, colNA = switch(var,
                'prlr' = 'white', grey(0.4)), filled.continents = switch(var,
                'tas' = F, 'tos' = T, 'prlr' = F), dots = t(flag), intylat = 45)
    dev.off()
  }
}

if (compratrms) { 
  ratrms <- RatioRMS(Mean1Dim(ano_exp, 2)[1, , 1, , ], 
                     Mean1Dim(ano_exp, 2)[2, , 1, , ], 
                     Mean1Dim(ano_obs, 2)[1, , 1, , ], 1)
  cols <- rev(c("dodgerblue4", "dodgerblue1", "forestgreen", "yellowgreen",
                "white", "white", "yellow", "orange", "red", "saddlebrown"))
  lims <- c(0, 0.5, 0.8, 0.9, 0.95, 1, 1.05, 1.1, 1.2, 2, 6)
  flag <- array(F, dim = dim(ratrms[1, , ]))
  flag[which(ratrms[2, , ] < 0.05)] <- T
  postscript(paste('Rati_RMSE2d_', var, '_', lstexpid[1], '_', lstexpid[2], 
                   '_', season, '_', ly, '.eps', sep = ''))
  PlotEquiMap(ratrms[1, , ], toto$lon, toto$lat, toptitle = paste('RMSE',
              lstexpid[1], '/ RMSE', lstexpid[2], switch(season, 
              'Year' = 'annual', season), switch(var, 
              'tas' = 'near surface temperature', 
              'tos' = 'sea surface temperature', 'prlr' = 'precipitation'),
              switch(var, 'tas' = 'GHCNv2+ERSSTv3b+GISTEMP', 'tas' = 'ERSSTv3b',
              'prlr' = 'GPCC'), switch(ly, 'ly1' = 'Year1', 'ly2-5' = 'Year2-5',
              'ly6-9' = 'Year6-9')), sizetit = 0.8, brks = lims, cols = cols,
              colNA = switch(var, 'prlr' = 'white', grey(0.4)), 
              filled.continents = switch(var, 'tas' = F, 'tos' = T, 'prlr' = F),
              dots = t(flag), intylat = 45)
  dev.off()
}
