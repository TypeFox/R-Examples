## ----install rEDM from package file, eval = FALSE------------------------
#  install.packages(file.choose(), type = "source", repos = NULL)

## ----install rEDM from GitHub, eval = FALSE------------------------------
#  library(devtools)
#  install_github("ha0ye/rEDM")

## ----load tentmap data---------------------------------------------------
library(rEDM)
data(tentmap_del)
head(tentmap_del)

## ----lib and pred for tentmap--------------------------------------------
lib <- c(1, 100)
pred <- c(201, 500)

## ----simplex on tentmap--------------------------------------------------
ts <- tentmap_del
simplex_output <- simplex(ts, lib, pred)

## ----rho vs. E for tentmap, tidy = TRUE, fig.width = 5, fig.height = 3.5----
par(mar = c(4,4,1,1), mgp = c(2.5, 1, 0))
plot(simplex_output$E, simplex_output$rho, type = "l", xlab = "Embedding Dimension (E)", ylab = "Forecast Skill (rho)")

## ----simplex varying tp for tentmap--------------------------------------
simplex_output <- simplex(ts, lib, pred, E = 2, tp = 1:10)

## ----rho vs. tp for tentmap, tidy = TRUE, fig.width = 5, fig.height = 3.5----
par(mar = c(4,4,1,1))
plot(simplex_output$tp, simplex_output$rho, type = "l", xlab = "Time to Prediction (tp)", ylab = "Forecast Skill (rho)")

## ----smap for tentmap----------------------------------------------------
smap_output <- s_map(ts, lib, pred, E = 2)

## ----rho vs. theta for tentmap, tidy = TRUE, fig.width = 5, fig.height = 3.5----
par(mar = c(4,4,1,1), mgp = c(2.5, 1, 0))
plot(smap_output$theta, smap_output$rho, type = "l", xlab = "Nonlinearity (theta)", ylab = "Forecast Skill (rho)")

## ----load block_3sp data-------------------------------------------------
data(block_3sp)
head(block_3sp)

## ----block_lnlp for block_3sp, tidy = TRUE, warning = FALSE--------------
lib <- c(1, NROW(block_3sp))
pred <- c(1, NROW(block_3sp))

block_lnlp_output <- block_lnlp(block_3sp, lib = lib, pred = pred, columns = c(1,2,4), target_column = 1, stats_only = FALSE, first_column_time = TRUE)

## ----block_lnlp with named columns, tidy = TRUE, eval = FALSE------------
#  block_lnlp_output <- block_lnlp(block_3sp, lib = lib, pred = pred, columns = c("x_t", "x_t-1", "y_t"), target_column = "x_t", stats_only = FALSE, first_column_time = TRUE)

## ----observed vs predicted for block_lnlp, tidy = TRUE, fig.width = 4, fig.height = 4----
observed <- block_lnlp_output[[1]]$model_output$obs
predicted <- block_lnlp_output[[1]]$model_output$pred

par(mar = c(4,4,1,1), pty = "s")
plot_range <- range(c(observed, predicted), na.rm = TRUE)
plot(observed, predicted, xlim = plot_range, ylim = plot_range, xlab = "Observed", ylab = "Predicted")
abline(a = 0, b = 1, lty = 2, col = "blue")

## ----sardine anchovy ccm, tidy = TRUE, warning = FALSE, cache = TRUE-----
data(sardine_anchovy_sst)
anchovy_xmap_sst <- ccm(sardine_anchovy_sst, E = 3, lib_column = "anchovy", target_column = "np_sst", lib_sizes = seq(10, 80, by = 10), random_libs = FALSE)
sst_xmap_anchovy <- ccm(sardine_anchovy_sst, E = 3, lib_column = "np_sst", target_column = "anchovy", lib_sizes = seq(10, 80, by = 10), random_libs = FALSE)

## ----sardine anchovy ccm plot, tidy = TRUE, fig.width = 5, fig.height = 3.5----
a_xmap_t_means <- ccm_means(anchovy_xmap_sst)
t_xmap_a_means <- ccm_means(sst_xmap_anchovy)

par(mar = c(4,4,1,1), mgp = c(2.5, 1, 0))
plot(a_xmap_t_means$lib_size, pmax(0, a_xmap_t_means$rho), type = "l", col = "red", xlab = "Library Size", ylab = "Cross Map Skill (rho)", ylim = c(0, 0.4))
lines(t_xmap_a_means$lib_size, pmax(0, t_xmap_a_means$rho), col = "blue")
legend(x = "topleft", legend = c("anchovy xmap SST", "SST xmap anchovy"), col = c("red", "blue"), lwd = 1, inset = 0.02, cex = 0.8)

## ----load e054 data, tidy = TRUE-----------------------------------------
data(e054_succession)
head(e054_succession)

# separate time column from data
composite_ts <- e054_succession[,-c(1:5)]

# normalize each time series
n <- NCOL(composite_ts)
blocks <- paste(e054_succession$Exp, e054_succession$OldField, e054_succession$Transect)
blocks_index <- sort(unique(blocks))
for(j in 1:n) {
    for(i in 1:length(blocks_index)) {
        subs <- which(blocks == blocks_index[i])
        composite_ts[subs,j] <- (composite_ts[subs,j] - mean(composite_ts[subs,j])) / sd(composite_ts[subs,j])
        }
    }
composite_ts <- cbind(year = e054_succession$Year, composite_ts)

## ----library for e054----------------------------------------------------
# make composite library
segments <- NULL
startpos <- 1
for(i in 2:NROW(composite_ts)) {
    if(composite_ts$year[i] < composite_ts$year[i-1]) {
        segments <- rbind(segments, c(startpos, i))
        startpos <- i+1
        }
    }
segments <- rbind(segments, c(max(segments)+1, NROW(composite_ts)))

composite_lib <- segments[1:floor(NROW(segments)/2),]
composite_pred <- segments[(floor(NROW(segments)/2)+1):NROW(segments),]

## ----precipitation variable for e054-------------------------------------
precip_ts <- unique(e054_succession[,c("Year", "precipmm")])
precip_ts <- precip_ts[order(precip_ts$Year),]

## ----simplex for e054, tidy = TRUE, warning = FALSE, fig.width = 6, results = "hold"----
par(mar = c(4,4,1,1), mfrow = c(2,2), mgp = c(2,1,0))
varlst <- colnames(composite_ts)[2:4]
simplex_output_list <- NULL

for(i in 1:length(varlst)) {
    simplex_output_list[[i]] <- simplex(composite_ts[,varlst[i]], lib = composite_lib, pred = composite_pred, E = 2:6)
    plot(simplex_output_list[[i]]$E, simplex_output_list[[i]]$rho, type = "l", xlab = "Embedding Dimension (E)", ylab = "Forecast Skill (rho)", main = varlst[i])
    }

simplex_output_list[[4]] <- simplex(precip_ts, E = 2:6, silent = TRUE)
names(simplex_output_list) <- c(varlst, "precipmm")
plot(simplex_output_list[[4]]$E, simplex_output_list[[4]]$rho, type = "l", xlab = "Embedding Dimension (E)", ylab = "Forecast Skill (rho)", main = "Precip")

## ----get best E----------------------------------------------------------
bestE <- sapply(simplex_output_list, function(simplex_output) {
    simplex_output$E[which.max(simplex_output$rho)]
    })
bestE

## ----smap for e054, cache = TRUE, tidy = TRUE, warning = FALSE, fig.width = 6----
par(mar = c(4,4,1,1), mfrow=c(2,2), mgp = c(2,1,0))
smap_output_list <- NULL

for(i in 1:length(varlst)) {
    smap_output_list[[i]] <- s_map(composite_ts[,c("year", varlst[i])], lib = composite_lib, pred = composite_pred, E = bestE[i])
    plot(smap_output_list[[i]]$theta, smap_output_list[[i]]$rho, type = "l", xlab = "Nonlinearity (theta)", ylab = "Forecast Skill (rho)", main=varlst[i])
    }

smap_output_list[[4]] <- s_map(precip_ts, lib = c(1,24), pred = c(1,24), E = bestE[4])
plot(smap_output_list[[4]]$theta, smap_output_list[[4]]$rho, type = "l", xlab = "Nonlinearity (theta)", ylab = "Forecast Skill (rho)", main="Precip")

## ----make blocks for e054------------------------------------------------
n <- NROW(composite_ts)

#Make lags
block_e54 <- data.frame(year=composite_ts$year)
block_e54$AR_tm <- composite_ts$Agropyron.repens
block_e54$AR_tm1 <- c(NA, block_e54$AR_tm[-n])
block_e54$AR_tm2 <- c(NA, block_e54$AR_tm1[-n])
block_e54$AR_tm3 <- c(NA, block_e54$AR_tm2[-n])

block_e54$SS_tm <- composite_ts$Schizachyrium.scoparium
block_e54$SS_tm1 <- c(NA, block_e54$SS_tm[-n])
block_e54$SS_tm2 <- c(NA, block_e54$SS_tm1[-n])
block_e54$SS_tm3 <- c(NA, block_e54$SS_tm2[-n])

block_e54$ML_tm <- composite_ts$Miscellaneous.litter
block_e54$ML_tm1 <- c(NA, block_e54$ML_tm[-n])
block_e54$ML_tm2 <- c(NA, block_e54$ML_tm1[-n])
block_e54$ML_tm3 <- c(NA, block_e54$ML_tm2[-n])

block_e54$PR_tm <- composite_ts$precipmm
block_e54$PR_tm1 <- c(NA, block_e54$PR_tm[-n])
block_e54$PR_tm2 <- c(NA, block_e54$PR_tm1[-n])
block_e54$PR_tm3 <- c(NA, block_e54$PR_tm2[-n])

#Remove overlaps from other plots
startyear <- 2001
for(i in 2:nrow(block_e54)) {
    if(block_e54$year[i]<block_e54$year[i-1]) {
        startyear <- block_e54$year[i]
        }
    if(block_e54$year[i]==startyear) {
        block_e54[i,c("AR_tm1", "SS_tm1", "ML_tm1", "PR_tm1")] <- NA
        block_e54[i,c("AR_tm2", "SS_tm2", "ML_tm2", "PR_tm2")] <- NA
        block_e54[i,c("AR_tm3", "SS_tm3", "ML_tm3", "PR_tm3")] <- NA
        }
    if(block_e54$year[i]==(startyear+1)) {
        block_e54[i,c("AR_tm2", "SS_tm2", "ML_tm2", "PR_tm2")] <- NA
        block_e54[i,c("AR_tm3", "SS_tm3", "ML_tm3", "PR_tm3")] <- NA
        }
    if(block_e54$year[i]==(startyear+2)) {
        block_e54[i,c("AR_tm3", "SS_tm3", "ML_tm3", "PR_tm3")] <- NA
        }
    }
head(block_e54)

## ----block_lnlp for e054, warning = FALSE, tidy = TRUE-------------------
block_lnlp_output_AR <- block_lnlp(block_e54, lib = composite_lib, pred = composite_pred, columns = c("AR_tm", "AR_tm1", "AR_tm2", "AR_tm3"), target_column = "AR_tm", stats_only = FALSE)

block_lnlp_output_ML <- block_lnlp(block_e54, lib = composite_lib, pred = composite_pred, columns = c("AR_tm", "AR_tm1", "AR_tm2", "AR_tm3", "ML_tm1", "PR_tm"), target_column = "AR_tm", stats_only = FALSE)

block_lnlp_output_SS <- block_lnlp(block_e54, lib = composite_lib, pred = composite_pred, columns = c("AR_tm", "AR_tm1", "AR_tm2", "AR_tm3", "ML_tm1", "PR_tm", "SS_tm3"), target_column = "AR_tm", stats_only = FALSE)

## ----block_lnlp plot for e054, fig.width = 6, tidy = TRUE, warning = FALSE----
observed_AR <- block_lnlp_output_AR[[1]]$model_output$obs
predicted_AR <- block_lnlp_output_AR[[1]]$model_output$pred

observed_ML <- block_lnlp_output_ML[[1]]$model_output$obs
predicted_ML <- block_lnlp_output_ML[[1]]$model_output$pred

observed_SS <- block_lnlp_output_SS[[1]]$model_output$obs
predicted_SS <- block_lnlp_output_SS[[1]]$model_output$pred

par(mar = c(4,4,1,1), pty = "s")
plot_range <- range(c(observed_AR, predicted_AR), na.rm = TRUE)
plot(observed_AR, predicted_AR, xlim = plot_range, ylim = plot_range, xlab = "Observed", ylab = "Predicted")
abline(a = 0, b = 1, lty = 2, col = "darkgrey")

abline(lm(predicted_AR~observed_AR), col="black", lty=3)
points(observed_ML, predicted_ML, pch=2, col="red")
abline(lm(predicted_ML~observed_ML), col="red", lty=3)
points(observed_SS, predicted_SS, pch=3, col="blue")
abline(lm(predicted_SS~observed_SS), col="blue", lty=3)

## ---- ccm on e054, cache = TRUE, warning = FALSE, tidy = TRUE, fig.width = 6----
#A. repens:
ar_xmap_prec <- ccm(composite_ts, lib = segments, pred = segments, E = bestE[4], lib_column = "Agropyron.repens", target_column = "precipmm", lib_sizes = seq(1, 1000, by = 100), num_samples = 10)
prec_xmap_ar <- ccm(composite_ts, lib = composite_lib, pred = composite_pred, E = bestE[1], lib_column = "precipmm", target_column = "Agropyron.repens", lib_sizes = seq(1, 1000, by = 100), num_samples = 10)

a_xmap_p_means <- ccm_means(ar_xmap_prec)
p_xmap_a_means <- ccm_means(prec_xmap_ar)

#S. scoparium:
ss_xmap_prec <- ccm(composite_ts, lib = segments, pred = segments, E = bestE[4], lib_column = "Schizachyrium.scoparium", target_column = "precipmm", lib_sizes = seq(1, 1000, by = 100), num_samples = 10)
prec_xmap_ss <- ccm(composite_ts, lib = segments, pred = segments, E = bestE[2], lib_column = "precipmm", target_column = "Schizachyrium.scoparium", lib_sizes = seq(1, 1000, by = 100), num_samples = 10)

s_xmap_p_means <- ccm_means(ss_xmap_prec)
p_xmap_s_means <- ccm_means(prec_xmap_ss)

#Plot output
par(mar = c(4,4,1,1), mgp = c(2.5, 1, 0))
plot(s_xmap_p_means$lib_size, pmax(0, s_xmap_p_means$rho), type = "l", col = "red", xlab = "Library Size", ylab = "Cross Map Skill (rho)", ylim = c(0, 0.2))
lines(p_xmap_s_means$lib_size, pmax(0, p_xmap_s_means$rho), col = "blue")
legend(x = "topleft", legend = c("A. repens xmap Precip.", "Precip. xmap A. repens"), col = c("red", "blue"), lwd = 1, inset = 0.02, bty = "n")

plot(a_xmap_p_means$lib_size, pmax(0, a_xmap_p_means$rho), type = "l", col = "orange", xlab = "Library Size", ylab = "Cross Map Skill (rho)", ylim = c(0, 0.2))
lines(p_xmap_a_means$lib_size, pmax(0, p_xmap_a_means$rho), col = "blue")
legend(x = "topleft", legend = c("S. scoparium xmap Precip.", "Precip. xmap S. scoparium"), col = c("orange", "blue"), lwd = 1, inset = 0.02, bty = "n")

## ----load e120 data, tidy = TRUE-----------------------------------------
data(e120_biodiversity)
head(e120_biodiversity)

# separate time column from data
composite_ts <- e120_biodiversity[,c(7:9,12)]

# normalize each time series
n <- NCOL(composite_ts)
blocks <- e120_biodiversity$Plot
blocks_index <- sort(unique(blocks))
for(j in 1:n) {
    for(i in 1:length(blocks_index)) {
        subs <- which(blocks == blocks_index[i])
        composite_ts[subs,j] <- (composite_ts[subs,j] - mean(composite_ts[subs,j])) / sd(composite_ts[subs,j])
        }
    }

composite_ts <- cbind(year = e120_biodiversity$Year, composite_ts)

## ----make composite ts for e120------------------------------------------
#make composite library
segments <- NULL
startpos <- 1
for(i in 2:nrow(composite_ts)) {
    if(composite_ts$year[i] < composite_ts$year[i-1]) {
        segments <- rbind(segments, c(startpos, i))
        startpos <- i+1
        }
    }
segments <- rbind(segments, c(max(segments)+1, nrow(composite_ts)))

#Choose random segments for prediction
set.seed(2312)
rndlib <- sort(sample(1:nrow(segments), round(nrow(segments)/2,0), rep=FALSE))
composite_lib <- segments[rndlib,]
composite_pred <- segments[-rndlib,]

## ----precip for e120-----------------------------------------------------
precip_ts <- unique(e120_biodiversity[,c("Year", "SummerPrecip.mm.")])
precip_ts <- precip_ts[order(precip_ts$Year),]

## ----simplex on e120, tidy = TRUE, cache = TRUE, warning = FALSE, fig.width = 6----
par(mar = c(4,4,1,1), mfrow=c(2,2), mgp = c(2.5, 1, 0))
varlst <- colnames(composite_ts)[2:4]
simplex_output_list <- NULL

for(i in 1:length(varlst)) {
    simplex_output_list[[i]] <- simplex(composite_ts[,c("year", varlst[i])], lib=composite_lib, pred=composite_pred, E = c(2:6))
    plot(simplex_output_list[[i]]$E, simplex_output_list[[i]]$rho, type = "l", xlab = "Embedding Dimension (E)", ylab = "Forecast Skill (rho)", main=varlst[i])
    }

simplex_output_list[[4]] <- simplex(precip_ts, lib = c(1,7), pred = c(1,7), E = c(2:5), silent = TRUE)
names(simplex_output_list) <- c(varlst, "precipmm")
plot(simplex_output_list[[4]]$E, simplex_output_list[[4]]$rho, type = "l", xlab = "Embedding Dimension (E)", ylab = "Forecast Skill (rho)", main="Precip")

## ----best E for e120-----------------------------------------------------
bestE <- sapply(simplex_output_list, function(simplex_output) {
    simplex_output$E[which.max(simplex_output$rho)]
    })
bestE

## ----smap on e120, cache = TRUE, warning = FALSE, tidy = TRUE, fig.width = 6----
par(mar = c(4,4,1,1), mfrow=c(2,2), mgp = c(2.5, 1, 0))
smap_output_list <- NULL

for(i in 1:length(varlst)) {
    smap_output_list[[i]] <- s_map(composite_ts[,c("year", varlst[i])], lib = composite_lib, pred = composite_pred, E = bestE[i], silent = TRUE)
    plot(smap_output_list[[i]]$theta, smap_output_list[[i]]$rho, type = "l", xlab = "Nonlinearity (theta)", ylab = "Forecast Skill (rho)", main = varlst[i])
    }

smap_output_list[[4]] <- s_map(precip_ts, E = bestE[4], silent = TRUE)
plot(smap_output_list[[4]]$theta, smap_output_list[[4]]$rho, type = "l", xlab = "Nonlinearity (theta)", ylab = "Forecast Skill (rho)", main = "Precip")

## ----make block for e120, tidy = TRUE------------------------------------
n <- NROW(composite_ts)

#Make lags
block_data <- data.frame(year=composite_ts$year)
block_data$AB_tm <- composite_ts$AbvBioAnnProd
block_data$AB_tm1 <- c(NA, block_data$AB_tm[-n])
block_data$AB_tm2 <- c(NA, block_data$AB_tm1[-n])
block_data$AB_tm3 <- c(NA, block_data$AB_tm2[-n])

block_data$NO_tm <- composite_ts$noh020tot
block_data$NO_tm1 <- c(NA, block_data$NO_tm[-n])
block_data$NO_tm2 <- c(NA, block_data$NO_tm1[-n])
block_data$NO_tm3 <- c(NA, block_data$NO_tm2[-n])

block_data$IV_tm <- composite_ts$invrichness
block_data$IV_tm1 <- c(NA, block_data$IV_tm[-n])
block_data$IV_tm2 <- c(NA, block_data$IV_tm1[-n])
block_data$IV_tm3 <- c(NA, block_data$IV_tm2[-n])

block_data$PR_tm <- composite_ts$SummerPrecip.mm
block_data$PR_tm1 <- c(NA, block_data$PR_tm[-n])
block_data$PR_tm2 <- c(NA, block_data$PR_tm1[-n])
block_data$PR_tm3 <- c(NA, block_data$PR_tm2[-n])

#Remove overlaps from other plots
startyear <- 1996
for(i in 2:nrow(block_data)) {
    if(block_data$year[i]<block_data$year[i-1]) {
        startyear <- block_data$year[i]
        }
    if(block_data$year[i]==startyear) {
        block_data[i,c("AB_tm1", "NO_tm1", "IV_tm1", "PR_tm1")] <- NA
        block_data[i,c("AB_tm2", "NO_tm2", "IV_tm2", "PR_tm2")] <- NA
        block_data[i,c("AB_tm3", "NO_tm3", "IV_tm3", "PR_tm3")] <- NA
        }
    if(block_data$year[i]==(startyear+1)) {
        block_data[i,c("AB_tm2", "NO_tm2", "IV_tm2", "PR_tm2")] <- NA
        block_data[i,c("AB_tm3", "NO_tm3", "IV_tm3", "PR_tm3")] <- NA
        }
    if(block_data$year[i]==(startyear+2)) {
        block_data[i,c("AB_tm3", "NO_tm3", "IV_tm3", "PR_tm3")] <- NA
        }
    }
head(block_data[,1:5],20)

## ----block_lnlp for e120, tidy=TRUE, warning = FALSE, cache = TRUE-------
block_lnlp_output_AB <- block_lnlp(block_data, lib = composite_lib, pred = composite_pred, columns = c("AB_tm", "AB_tm1", "AB_tm2"), target_column = 1, stats_only = FALSE, first_column_time = TRUE)

block_lnlp_output_ABPR <- block_lnlp(block_data, lib = composite_lib, pred = composite_pred, columns = c("AB_tm", "AB_tm1", "AB_tm2", "PR_tm", "PR_tm1", "PR_tm2"), target_column = 1, stats_only = FALSE, first_column_time = TRUE)

block_lnlp_output_ABPR <- block_lnlp(block_data, lib = composite_lib, pred = composite_pred, columns = c("AB_tm", "AB_tm1", "AB_tm2", "PR_tm", "PR_tm1", "PR_tm2"), target_column = 1, stats_only = FALSE, first_column_time = TRUE)

## ----block_lnlp on e120, tidy = TRUE, warning = FALSE, fig.width = 4, fig.height = 4----
observed_AB <- block_lnlp_output_AB[[1]]$model_output$obs
predicted_AB <- block_lnlp_output_AB[[1]]$model_output$pred

observed_ABPR <- block_lnlp_output_ABPR[[1]]$model_output$obs
predicted_ABPR <- block_lnlp_output_ABPR[[1]]$model_output$pred

par(mar = c(4,4,1,1), pty = "s", mgp = c(2.5, 1, 0))
plot_range <- range(c(observed_AB, predicted_AB), na.rm = TRUE)
plot(observed_AB, predicted_AB, xlim = plot_range, ylim = plot_range, xlab = "Observed", ylab = "Predicted")
abline(a = 0, b = 1, lty = 2, col = "darkgrey", lwd=2)
abline(lm(predicted_AB~observed_AB), col="black", lty=3, lwd=2)

points(observed_ABPR, predicted_ABPR, pch=2, col="red")
abline(lm(predicted_ABPR~observed_ABPR), col="red", lty=3, lwd=2)

## ----ccm on e120, cache = TRUE, warning = FALSE, tidy = TRUE, fig.width = 5, fig.height = 3.5----
# A. repens:
no_xmap_inv <- ccm(composite_ts, lib=segments, pred=segments, E = bestE[4], lib_column = "noh020tot", target_column = "invrichness", lib_sizes = c(seq(5, 55, by=2), seq(55, 400, by=50)), num_samples = 100, silent = TRUE)
inv_xmap_no <- ccm(composite_ts, lib=composite_lib, pred=composite_pred, E = bestE[1], lib_column = "invrichness", target_column = "noh020tot", lib_sizes = c(seq(5, 55, by=2), seq(55, 400, by=50)), num_samples = 100, silent = TRUE)

n_xmap_i_means <- data.frame(ccm_means(no_xmap_inv), sd.rho = with(no_xmap_inv, tapply(rho, lib_size, sd)))
i_xmap_n_means <- data.frame(ccm_means(inv_xmap_no), sd.rho = with(inv_xmap_no, tapply(rho, lib_size, sd)))

# S. scoparium:
ab_xmap_inv <- ccm(composite_ts, lib=segments, pred=segments, E = bestE[4], lib_column = "AbvBioAnnProd", target_column = "invrichness", lib_sizes = c(seq(5, 55, by=2), seq(55, 400, by=50)), num_samples = 100, silent = TRUE)
inv_xmap_ab <- ccm(composite_ts, lib=segments, pred=segments, E = bestE[2], lib_column = "invrichness", target_column = "AbvBioAnnProd", lib_sizes = c(seq(5, 55, by=2), seq(55, 400, by=50)), num_samples = 100, silent = TRUE)

a_xmap_i_means <- data.frame(ccm_means(ab_xmap_inv), sd.rho=with(ab_xmap_inv, tapply(rho, lib_size, sd)))
i_xmap_a_means <- data.frame(ccm_means(inv_xmap_ab), sd.rho=with(inv_xmap_ab, tapply(rho, lib_size, sd)))

# Plot output
par(mar = c(4,4,1,1))
plot(n_xmap_i_means$lib_size, pmax(0, n_xmap_i_means$rho), type = "l", col = "red", xlab = "Library Size", ylab = "Cross Map Skill (rho)", ylim = c(0, 0.6), lwd=2)
lines(i_xmap_n_means$lib_size, pmax(0, i_xmap_n_means$rho), col = "blue", lwd=2)
legend(x = "topleft", legend = c("Nitrate xmap Inv. Richness", "Inv. Richness xmap Nitrate"), col = c("red", "blue"), lwd = 2, inset = 0.02, bty="n", cex = 0.8)
abline(h=0, lty=3, col="darkgrey", lwd=2)

#Add CI's
lines(n_xmap_i_means$lib_size, n_xmap_i_means$rho+n_xmap_i_means$sd.rho, col = "red", lty=2, lwd=2)
lines(n_xmap_i_means$lib_size, n_xmap_i_means$rho-n_xmap_i_means$sd.rho, col = "red", lty=2, lwd=2)
lines(i_xmap_n_means$lib_size, i_xmap_n_means$rho+i_xmap_n_means$sd.rho, col = "blue", lty=2, lwd=2)
lines(i_xmap_n_means$lib_size, i_xmap_n_means$rho-i_xmap_n_means$sd.rho, col = "blue", lty=2, lwd=2)

plot(a_xmap_i_means$lib_size, pmax(0, a_xmap_i_means$rho), type = "l", col = "orange", xlab = "Library Size", ylab = "Cross Map Skill (rho)", ylim = c(0, 0.6), lwd=2)
lines(i_xmap_a_means$lib_size, pmax(0, i_xmap_a_means$rho), col = "blue", lwd=2)
legend(x = "topleft", legend = c("Abv. Biomass xmap Inv. Richness", "Inv. Richness xmap Abv. Biomass"), col = c("orange", "blue"), lwd = 2, inset = 0.02, bty="n", cex = 0.8)
abline(h=0, lty=3, col="darkgrey", lwd=2)

#Add CI's
lines(a_xmap_i_means$lib_size, a_xmap_i_means$rho+a_xmap_i_means$sd.rho, col = "orange", lty=2, lwd=2)
lines(a_xmap_i_means$lib_size, a_xmap_i_means$rho-a_xmap_i_means$sd.rho, col = "orange", lty=2, lwd=2)
lines(i_xmap_a_means$lib_size, i_xmap_a_means$rho+i_xmap_a_means$sd.rho, col = "blue", lty=2, lwd=2)
lines(i_xmap_a_means$lib_size, i_xmap_a_means$rho-i_xmap_a_means$sd.rho, col = "blue", lty=2, lwd=2)

## ----thrips data---------------------------------------------------------
data(thrips_block)
colnames(thrips_block)

## ----thrips plot, echo = FALSE, fig.width = 6, fig.height = 7------------
par(mar = c(4,4,1,1), mfrow = c(4,1), mgp = c(2.5,1,0))
time_dec <- thrips_block$Year + (thrips_block$Month)/12
plot(time_dec, thrips_block$Thrips_imaginis, type='l', col = 'green', ylab = 'Thrips')
plot(time_dec, thrips_block$maxT_degC, type='l', col = 'red', ylab = 'maxT (oC)')
plot(time_dec, thrips_block$Rain_mm, type='l', col = 'blue', ylab = 'Rain (mm)')
plot(time_dec, thrips_block$Season, type='l', col = 'magenta', ylab = 'Season')

## ----univariate thrips, warning = FALSE----------------------------------
ts <- thrips_block$Thrips_imaginis
lib <- c(1, length(ts))
pred <- c(1, length(ts))
simplex_output <- simplex(ts, lib, pred, tau = 1)

## ----rho vs. e for thrips, echo=FALSE, fig.width = 5, fig.height = 3.5, tidy = TRUE----
par(mar = c(4,4,1,1))
plot(simplex_output$E, simplex_output$rho, type = "l", xlab = "Embedding Dimension (E)", ylab = "Forecast Skill (rho)")

## ----smap for thrips, warning = FALSE------------------------------------
smap_output <- list()
smap_output[[1]] <- s_map(ts, lib, pred, E = 4)
smap_output[[2]] <- s_map(ts, lib, pred, E = 8)

## ----rho vs. theta for thrips, echo=FALSE, tidy = TRUE, fig.width = 6, fig.height = 3.5----
par(mar = c(4,4,1,1), mfrow = c(1,2))
plot(smap_output[[1]]$theta, smap_output[[1]]$rho, type = "l", xlim=c(0,4), xlab = "Nonlinearity (theta)", ylab = "Forecast Skill (rho)")
plot(smap_output[[2]]$theta, smap_output[[2]]$rho, type = "l", xlim=c(0,4), xlab = "Nonlinearity (theta)", ylab = "Forecast Skill (rho)")

## ----compute ccm matrix for thrips, results='hold', tidy=TRUE, cache = TRUE, warning = FALSE----
ncol <- dim(thrips_block)[2]-2
M_rho <- array(NA,dim=c(ncol,ncol),dimnames=list(colnames(thrips_block[3:6]),colnames(thrips_block[3:6])))

for (i in 1:ncol){
    for (j in 1:ncol){
        if (i!=j){
            out_temp <- ccm(thrips_block,E=8,lib_column=2+i,target_column=2+j,
                            lib_sizes = dim(thrips_block)[1],replace=FALSE, silent = TRUE)
            M_rho[i,j] <- out_temp$rho
            } 
        }
    }

## ----compute corr matrix for thrips, tidy=TRUE---------------------------
M_corr <- array(NA,dim=c(ncol,ncol),dimnames=list(colnames(thrips_block[3:6]),colnames(thrips_block[3:6])))

for (i in 1:ncol){
    for (j in 1:ncol){
        if (i!=j){
            cf_temp <- ccf(x=thrips_block[,2+i], y=thrips_block[,2+j], type = "correlation", lag.max = 6, plot = FALSE)$acf
            M_corr[i,j] <- max(abs(cf_temp))
            }
        }
}

## ----xmap matrix for thrips, echo=FALSE----------------------------------
head(M_rho)

## ----corr matrix for thrips, echo=FALSE----------------------------------
head(M_corr)

## ----ccm on thrips, results='hide', tidy=TRUE, cache = TRUE, warning=FALSE----
thrips_xmap_maxT <- ccm(thrips_block, E = 8, random_libs = TRUE, lib_column = "Thrips_imaginis", target_column = "maxT_degC", lib_sizes = seq(10, 75, by = 5), num_samples = 300)
maxT_xmap_thrips <- ccm(thrips_block, E = 8, random_libs = TRUE, lib_column = "maxT_degC", target_column = "Thrips_imaginis", lib_sizes = seq(10, 75, by = 5), num_samples = 300)

xmap_means <- list(ccm_means(thrips_xmap_maxT),ccm_means(maxT_xmap_thrips))

## ----ccm plot, echo=FALSE, tidy = TRUE, fig.width = 5, fig.height = 3.5----
par(mar = c(4,4,1,1), mgp = c(2.5, 1, 0))
plot(xmap_means[[1]]$lib_size, pmax(0, xmap_means[[1]]$rho), type = "l", col = "red",  xlab = "Library Size", ylab = "Cross Map Skill (rho)", ylim = c(0, 1))
lines(xmap_means[[2]]$lib_size, pmax(0, xmap_means[[2]]$rho), col = "blue")
abline(h=M_corr['Thrips_imaginis','maxT_degC'], col = "black", lty = 2)
legend(x = "bottomright", legend = c("Thrips xmap maxT", "maxT xmap Thrips"), col = c("red", "blue"), lwd = 1, inset = 0.02)

## ----ccm on thrips and rainfall, results='hide', tidy=TRUE, cache = TRUE, warning = FALSE----
thrips_xmap_Rain <- ccm(thrips_block, E = 8, random_libs = TRUE, lib_column = "Thrips_imaginis", target_column = "Rain_mm", lib_sizes = seq(10, 75, by = 5), num_samples = 300)
Rain_xmap_thrips <- ccm(thrips_block, E = 8, random_libs = TRUE, lib_column = "Rain_mm", target_column = "Thrips_imaginis", lib_sizes = seq(10, 75, by = 5), num_samples = 300, silent = TRUE)

xmap_means <- list(ccm_means(thrips_xmap_Rain),ccm_means(Rain_xmap_thrips))

## ----rainfall and thrips ccm plot, echo=FALSE, tidy = TRUE, fig.width = 5, fig.height = 3.5----
par(mar = c(4,4,1,1), mgp = c(2.5, 1, 0))
plot(xmap_means[[1]]$lib_size, pmax(0, xmap_means[[1]]$rho), type = "l", col = "red", xlab = "Library Size", ylab = "Cross Map Skill (rho)", ylim = c(0, 1))
lines(xmap_means[[2]]$lib_size, pmax(0, xmap_means[[2]]$rho), col = "blue")
abline(h=M_corr['Thrips_imaginis','Rain_mm'], col = 'black', lty = 2)
legend(x = "topleft", legend = c("Thrips xmap Rain", "Rain xmap Thrips"), 
 col = c("red", "blue"), lwd = 1, inset = 0.02)

## ----ccm on thrips and season, results='hide', tidy=TRUE, cache = TRUE, warning = FALSE----
thrips_xmap_Season <- ccm(thrips_block, E = 8, random_libs = TRUE, lib_column = "Thrips_imaginis", target_column = "Season", lib_sizes = seq(10, 75, by = 5), num_samples = 300)
Season_xmap_thrips <- ccm(thrips_block, E = 8, random_libs = TRUE, lib_column = "Season", target_column = "Thrips_imaginis", lib_sizes = seq(10, 75, by = 5), num_samples = 300)

xmap_means <- list(ccm_means(thrips_xmap_Season), ccm_means(Season_xmap_thrips))

## ----season and thrips ccm plot, echo=FALSE, fig.width = 5, fig.height = 3.5, tidy = TRUE----
par(mar = c(4,4,1,1), mgp = c(2.5, 1, 0))
plot(xmap_means[[1]]$lib_size, pmax(0, xmap_means[[1]]$rho), type = "l", col = "red", xlab = "Library Size", ylab = "Cross Map Skill (rho)", ylim = c(0, 1))
lines(xmap_means[[2]]$lib_size, pmax(0, xmap_means[[2]]$rho), col = "blue")
abline(h=M_corr['Thrips_imaginis','Season'], col = 'black', lty = 2)
legend(x = "bottomright", legend = c("Thrips xmap Season", "Season xmap Thrips"), col = c("red", "blue"), lwd = 1, inset = 0.02)

## ----seasonal surrogates for thrips, cache = TRUE, warning = FALSE, tidy = TRUE----
num_surr <- 1000
surr_maxT <- make_surrogate_data(thrips_block$maxT_degC, method = "seasonal", T_period = 12, num_surr = num_surr)
surr_Rain <- make_surrogate_data(thrips_block$Rain_mm, method = "seasonal", T_period = 12, num_surr = num_surr)

rho_surr <- data.frame(maxT = numeric(num_surr), Rain = numeric(num_surr))

for (i in 1:num_surr) {
    rho_surr$maxT[i] <- ccm(cbind(thrips_block$Thrips_imaginis, surr_maxT[,i]), E = 8, lib_column = 1, target_column = 2, lib_sizes = NROW(thrips_block), replace = FALSE)$rho
    
    rho_surr$Rain[i] <- ccm(cbind(thrips_block$Thrips_imaginis, surr_Rain[,i]), E = 8, lib_column = 1, target_column = 2, lib_sizes = NROW(thrips_block), replace = FALSE)$rho
    }

## ----significance of randomization test----------------------------------
(sum(M_rho['Thrips_imaginis','Rain_mm'] < rho_surr$Rain) + 1) / (length(rho_surr$Rain) + 1)
(sum(M_rho['Thrips_imaginis','maxT_degC'] < rho_surr$maxT) + 1) / (length(rho_surr$maxT) + 1)

## ----load and composite sockeye data-------------------------------------
data(sockeye_returns)

# separate time column from data
time <- sockeye_returns$year
sockeye_returns <- sockeye_returns[,-1]

# normalize each time series
n <- NCOL(sockeye_returns)
for(j in 1:n)
    {
    sockeye_returns[,j] <- (sockeye_returns[,j] - mean(sockeye_returns[,j])) / sd(sockeye_returns[,j])
    }

# make composite time series
composite_ts <- data.frame(year = time, 
                           returns = stack(sockeye_returns)$value)

## ----composite lib and pred for sockeye, tidy = TRUE---------------------
# make composite library
segments <- cbind(seq(from = 1, by = NROW(sockeye_returns), length.out = n),  seq(from = NROW(sockeye_returns), by = NROW(sockeye_returns), length.out = n))
composite_lib <- segments[1:5,]
composite_pred <- segments[6:9,]

## ----simplex for composite sockeye, tidy = TRUE, fig.width = 6-----------
simplex_output <- simplex(composite_ts, composite_lib, composite_pred)
par(mar = c(4,4,1,1), mgp = c(2.5, 1, 0))
plot(simplex_output$E, simplex_output$rho, type = "l",  xlab = "Embedding Dimension (E)", ylab = "Forecast Skill (rho)")

## ----smap for composite sockeye, tidy = TRUE, fig.width = 6--------------
smap_output <- s_map(composite_ts, composite_lib, composite_pred, E = 8)
par(mar = c(4,4,1,1), mgp = c(2.5, 1, 0))
plot(smap_output$theta, smap_output$rho, type = "l", xlab = "Nonlinearity (theta)", ylab = "Forecast Skill (rho)")

