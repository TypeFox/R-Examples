## ----captures,echo=FALSE-------------------------------------------------
library(knitr)
opts_chunk$set(fig.lp="figure:",tidy=TRUE, tidy.opts=list(width.cutoff=60))

fig1_cap <- "The amplification curves were generated with the 
\\textsl{AmpSim} function. All Cqs are unique since random 
values were added to the starting Cq of 25. The parameter $noise = 
0.03$ adds some scatter to the data used for plotting the amplification curves."

fig1_scap <- "Simulation of a qPCR experiment using \\textsl{AmpSim} 
function."


fig2_cap <- "Working principle of the \\textsl{th}.\\textsl{cyc} function.
The function provides two modes (\\textbf{A)} is the linear 
regression. \\textbf{B)} Quadratic regression for the calculation of the 
Cq. In both cases the highest R squared value determines how
many left and right neighbors are used above and the below the defined threshold level ."

fig2_scap <- "Working principle of the \\textsl{th}.\\textsl{cyc} function."

fig3_cap <- "Application of the \\textsl{th}.\\textsl{cyc} function for the analysis of 
ccPCR data. The \\textsl{CPP} function was used to pre-process the data. 
Subsequently, the data were analyzed using the \\textsl{th}.\\textsl{cyc} function 
in the linear regression mode. The threshold level ($r = 50$) was 
identical for all data. The Cq is shown in minutes. The range used 
for the calculation of the Cq is indicated in red. Negative curves are 
automatically excluded from the analysis if the 90\\% percentile is lower 
or equal to the threshold level ($r$)."
fig3_scap <- "Application of the \\textsl{th}.\\textsl{cyc} function for the analysis 
of ccPCR data."

fig4_cap <- "Application of the \\textsl{CPP} and \\textsl{th}.\\textsl{cyc} 
functions. \\textbf{A)} The raw data of the VIMCFX96\\_60 dataset were 
\\textbf{B)} pre-processed with the \\textsl{CPP} function and finally plotted. 
The parameter $trans$ was set to $TRUE$, which leads to a linear trend correction 
and base-lining. By default a Savitzky-Golay filter was used to smooth the data. 
The data were normalized between 0 and 1 ($method.norm = 'minm'$). \\textbf{C)} 
All Cqs were calculated with \\textsl{th}.\\textsl{cyc} function. The Cq for the 
raw data was $17.25 \\pm 0.5$ (at $r = 2575$) and $17.1 \\pm 0.1$ (at $r = 0.1$) 
for the pre-processed data. Our results indicate that the dispersion of the Cq 
values was slightly lower for the pre-processed data."
fig4_scap <- "Application of the \\textsl{CPP} and 
\\textsl{th}.\\textsl{cyc} functions."


fig5_cap <- "Comparison of the normalization methods with the \\textsl{CPP} 
function. The VIMCFX96\\_60 dataset (96-well plate cycler, Bio-Rad CFX96, 
EvaGreen detection) was used. \\emph{(A)} Plot of raw data for all amplification 
curves. The signals are superimposed to circa 2200 RFU and the inter-sample 
baseline and plateau shift is high. Note the positive trend 
(\\textcolor{red}{\\textendash}, fitted with an ordinary least squares method) 
in the background range of cycles 1 to 15. All subsequent plots were processed 
with the \\textsl{CPP} function. By default, the curves are base-lined, smoothed 
(Savitzky-Golay smoother) and the slope corrected by a linear regression ($trans 
= TRUE$). \\emph{(B)} base-lined raw data, \\emph{(C)} \\emph{Min-Max 
normalization}, \\emph{(D)} \\emph{Max normalization}, \\emph{(E)} 
\\emph{lugn-normalization} with a cut off 3\\% 
and \\emph{(F)} \\emph{zscore-normalization}."

fig5_scap <- "Comparison of the normalization methods with the \\textsl{CPP} function."

fig6_cap <- "Amplification standard curve simulation and regression analysis. 
\\emph{(A)} \\textsl{AmpSim} was used to synthesize a qPCR experiment of six 
dilutions (three replicates per dilution) standard samples. The Cqs were 
determined by the $SDM$ method (solid black vertical lines). \\emph{(B)} \\textsl{effcalc} 
was used to automatically perform a linear regression (decadic logarithm 
of input concentration versus the Cq). The 95\\% confidence interval is shown be 
the light-blue solid lines."
fig6_scap <- "Amplification standard curve simulation and regression 
analysis."

fig7_cap <- "Calculation of the amplification efficiency.
Data of a VideoScan HCU dilution experiment (C54 dataset) 
were analyzed. \\emph{(A)} Visualization of the raw data. One of the three 
dilutions contains a missing value due to a sensor error. \\emph{(B, top 
panel)} The \\textsl{CPP} function was used to baseline, remove the missing value 
(\\textcolor{red}{\\textendash}) and smooth 
(\\textcolor{black}{\\textendash}, \\textcolor{red}{\\textendash}, 
\\textcolor{green}{\\textendash}) the raw data. \\emph{(B, bottom panel)}. 
The Cqs ($SDM$ by the \\textsl{diffQ} function) were then determined. An amplification 
efficiency of 87.3~\\% was calculated with the \\textsl{effcalc} function."

fig7_scap <- "Calculation of the amplification efficiency."

fig8_cap <- "Imputation of missing values in the data for generating amplification curves. 
\\emph{(A)} Raw data were generated using the \\textsl{AmpSim} 
simulation function. \\emph{(B)} A missing value was introduced in the 
transition phase. The missing value was imputed either by \\emph{(C)} 
linear approximation or \\emph{(D)} a cubic spline approximation. The 
spline approximation nearly reconstituted the original curve."

fig8_scap <- "Imputation of missing values in the data for generating amplification curves."

fig9_cap <- "Quantification cycle (Cq) by the second derivative maximum 
method. Raw data (\\textbullet) were generated by the \\textsl{AmpSim} 
function. The inflection point is the point where 
the slope is maximum and the curvature is zero. The first derivative of the 
amplification curve has a first derivative maximum ($FDM$) at the 
inflection point. The second derivative maximum method ($SDM$) needs to 
differentiate a curve to the second order prior to quantification. The second 
derivative exhibits a zero-crossing at the $FDM$. The function $y = f(x)$ is 
numerically derived by the five-point stencil. This method is
assumption free regarding the function $f$. \\textsl{inder} calculates the 
approximate $SDM$. The $SDM$ is calculated from a derived cubic spline. 
Similarly, the first approximate derivative maximum ($FDM$), second derivative 
minimum ($SDm$), and approximate second derivative center ($SDC$, geometric mean 
of $SDM$ and $SDm$) are available. $FDM$, $SDm$ and $SDC$ values can be used to 
further characterize the amplification process."

fig9_scap <- "Quantification cycle (Cq) by the second derivative maximum 
method."

fig10_cap <- "Plot of all data from C127EGHP and calculate the $SDM$ (Second 
Derivative Maximum) values with the \\textsl{diffQ2} function.
\\emph{(A)} Plot the samples detected with EvaGreen and 
\\emph{(B)} shows the same samples detected with the Hydrolysis probe for 
\\textsl{MLC-2v}. \\emph{(C)} Stripchart of the Cq values (\\textbullet) with the 
median (\\textcolor{red}{\\textendash}) and the median absolute deviation 
(\\textcolor{blue}{\\textendash~\\textendash}). 
This result indicates, that the variance derived from detection with hydrolysis
probes is higher than for the samples detected with EvaGreen. Note: the $inder$ parameter is set as TRUE."
  
fig10_scap <- "Plot of all data from C127EGHP and calculate the $SDM$ (Second 
Derivative Maximum) values with the \\textsl{diffQ2} function."

fig11_cap <- "Amplification curve profiles from the Bio-Rad iQ5 thermo 
cycler 
for the human gene \\textit{HPRT1}.
\\emph{(A)} The \\textsl{C126EG595} dataset was used with 
96 replicates of equal starting numbers of template molecules. Vertical 
lines represent the Cq ($SDM$ method) determined by the $inder$ method on 
amplification curves fitted with a 5-parameter curve function. Curves with 
Cqs less than 14.5 are indicated in red (\\textcolor{red}{\\textendash}). 
\\emph{(B)} Second derivatives of the amplification curves. Note that 
after differentiation all inter sample baseline and plateau shifts are 
similar. \\emph{(C)} Histogram (class width = 0.05 Cq) of the Cq values 
($SDM$). Cqs were mainly at circa 15.7 (N = 80) while some amplification 
curves had a Cq less than 15.5 (N = 16)."
fig11_scap <- "Amplification curve profiles from the Bio-Rad iQ5 thermo 
cycler 
for the human gene \\textit{HPRT1}."

fig12_cap <- "Signal analysis using the VIMCFX96\\_60 dataset (96-well 
plate cycler (Bio-Rad CFX96)). All cycles (ROI: 1 -- 40) were analyzed by 
the \\textsl{MFIaggr} function. The density plot (right upper panel) and 
quantile-quantile analysis (right lower panel) show no normal 
distribution. Owing to the sigmoidal structure of the curve, the density function can be considered bimodal."

fig12_scap <- "Signal analysis using the VIMCFX96\\_60 dataset (96-well 
plate cycler (Bio-Rad CFX96))."

fig13_cap <- "Helicase-dependent amplification (HDA) of \\textsl{vimentin} (\textit{vim}) in the 
VideoScan Platform. The HDA was performed at 65 degrees Celsius. Three 
concentrations of input DNA (D1, D2, D3) were used. The amplification curves 
were smoothed by a moving average (windowsize 3) and base-lined by a robust 
linear regression by computing MM-type regression estimator. The 
\\textsl{th}.\\textsl{cyc} function was used to determine the time required to 
reach the threshold level of 0.05 (--) arbitrary fluorescence units." 

fig13_scap <- "Helicase-dependent amplification (HDA) of \\textsl{vimentin} (\textit{vim})."

fig14_cap <- "The \\textsl{plotCurves} function. Notice: The function plots many 
curves on one plot in separate cells allowing for quick assessment. Missing 
values were artificially introduced at random positions to selected curves of the 
VIMCFX96\\_60 data set (solid black line). A colored box (topleft of each plot) 
indicates the sample name and if the data contain missing values. The red rug 
indicates the position of the missing values. The red lined shows the 
amplification curve after unsupervised pre-processing (using an instance of 
\\textsl{CPP})."
fig14_scap <- "The \\textsl{plotCurves} function."

fig15_cap <- "Use of \\textsl{MFIaggr} to test for heteroskedacity. The data were 
aggregated with the \\textsl{MFIaggr} function and assigned to the object $res$. 
The standard deviation was transformed to the variance. The plot shows the cycle 
dependent variance measured at 60 degrees Celsius (annealing phase; A, B) and 69 
degrees Celsius (elongation phase, C, D) of 96 qPCR replicate amplification 
curves. The first cycles 1 to 10 and next the cycles 1 to 40 were analyzed. The 
Breusch-Pagan test of the \\textsl{MFIaggr} confirmed the heteroskedasticity in the 
amplification curve data. The VIMCFX96\\_60 and VIMCFX96\\_69 datasets were 
used."

fig15_scap <- "Use of \\textsl{MFIaggr} to test for heteroskedasticity."

fig16_cap <- "The \\textsl{inder} function calculates numeric derivatives on smoothed 
data, which results in data points not observable in reality. The rounder 
function averages such result to the real values of cycle number. An 
amplification curve was simulated with the \\textsl{AmpSim} function."

fig16_scap <- "Application of rounder to average numeric derivatives to the real 
values of cycle number."

fig17_cap <- "\\textsl{lm.coefs} a function to compute coefficients for linear 
models. The function is a wrapper for functions to perform normal (least 
squares) and robust linear regression. If the computation of the robust linear regression fails, 
then \\textsl{lm.coefs} performs a linear regression using the least squares 
method. lmrob, MM-type estimators for linear regression; rq, quantile regression 
fit; least, least squares linear regression; rfit, Rank-based estimates of 
regression coefficients. m, slope; n, asymmetry."

fig17_scap <- "\\textsl{lm.coefs} a function to compute coefficients for linear 
models."

fig18_cap <- "\\textsl{bg}.\\textsl{max} tries to estimate the range between the 
background and the plateau phase of an amplification reaction. \\emph{(A)} in 
absence and \\emph{(B)} presence of noise. The data were simulated with the 
\\textsl{AmpSim} function."

fig18_scap <- "\\textsl{bg}.\\textsl{max} tries to estimate the range between the 
background and the plateau phase of an amplification reaction"

fig19_cap <- "Application of the \\textsl{bg}.\\textsl{max} function. 
Amplification curve data from a capillary convective PCR were used \\emph{(A)} 
as raw data and \\emph{(B)} pre-processed (smoothed (moving average, window size 
3), base-lined and trend corrected (robust MM-estimator)) with the \\textsl{CPP} 
function. The output of the \\textsl{CPP} function was used by 
\\textsl{bg}.\\textsl{max} to detected the start and the end of the 
amplification reaction. The start and end were reliably estimated (range between 
``bg.stop'' and ``amp.stop''). There was no significant difference between raw 
and pre-processed data."

fig19_scap <- "Application of the \\textsl{bg}.\\textsl{max} function to detect 
the start and end of an amplification reaction in a capillary convective PCR."

fig20_cap <- "Inspection of the reps384 dataset. The reps384 dataset was used 
for the analysis of the impact of imputed missing values. Three areas of the 
curve data were defined as ``Linear phase'' (red, cycle 1 -- 10), ``Exponential 
phase'' (blue, cycle 11 -- 33), ``Plateau phase'' (green, cycle 34 -- 40)."

fig20_scap <- "Inspection of the reps384 dataset."

fig21_cap <- "Analysis and interpretation of real-time amplification curves. 
\\emph{(A)} Before procession: The fluorescence values are plotted against the cycle, which results 
in sigmoidal shaped amplification curves 
(\\textcolor{black}{\\textendash},~\\textcolor{red}{\\textendash}). Measurements may 
occasionally contain missing values (``NA'', \\textcolor{red}{\\textendash}) and 
outliers (orange circle, \\textcolor{black}{\\textendash}), due to noise 
introduced by the detection system or sensor errors. Outliers are often 
present in the first cycle due to sensor adjustments. The signal difference 
between the background phase (first cycles) and the plateau phase (last cycles) 
is quantifiable as signal-to-noise ratio (SNR). The SNR between different 
samples (e.g., \\textcolor{black}{\\textendash} and 
\\textcolor{black}{\\textendash}) can vary. For interpretation it is recommended 
to normalize the data. Negative samples (\\textcolor{blue}{\\textendash}) need 
to be identified automatically. \\emph{(B)} Pre-processed raw data, where NAs 
were imputed and the noise slightly removed. The curves were adjusted to have 
the same baseline and plateau level. The quantification cycles (Cq) of the 
positive reactions are determined in the exponential phase (``threshold method'' 
is used in this example). Negative samples are automatically set to zero." 

fig21_scap <- "Analysis and interpretation of real-time amplification curves."

fig22_cap <- "Smoother and filter methods of the \\emph{chipPCR} package. 
\\emph{(A)} Raw data were generated using the \\textsl{AmpSim} simulation 
function. \\emph{(B)} The difference of the raw data to the smoothed data was 
plotted. ``savgol'' (Savitzky-Golay smoothing), ``lowess'' (locally-weighted 
polynomial regression), ``mova3'' (moving average with window size of 3), 
``smooth'' (cubic smoothing spline), ``spline'' (Interpolating cubic spline), 
``supsmu'' (Friedman's SuperSmoother), ``whit1'' (weighted Whittaker smoothing 
with a finite difference penalty of order 1), ``whit2'' (weighted Whittaker 
smoothing with a finite difference penalty of order 2). The ``savgol'', 
``smooth'', ``spline'' ``whit1'' , and ``whit2'' nearly preserved the original 
curve. The other functions resulted in alterations in the transition phases of 
the amplification curve. Optimized time series smoother, like the Kalman filter 
\\citep{Tusell_2010}, are not yet integrated in this package."

fig22_scap <- "Smoother and filter methods of the \\emph{chipPCR} package."

fig23_cap <- "Sample for analysis using the \\textsl{MFIaggr} function. 
The VIMCFX96\\_60 dataset (96-well plate cycler (Bio--Rad CFX96)) was used. 
Either all or a subset of the cycles (ROI: 1 -- 10) or all cycles (ROI: 1 -- 40) 
(Figure~\\ref{figure:MFIaggr_all}) were analyzed. The density plot (right upper 
panel) and quantile-quantile analysis (right lower panel) show no normal 
distribution. Due to the sigmoidal structure of the curve, the density function 
can be considered bimodal."

fig23_scap <- "Signal analysis by the \\textsl{MFIaggr} function."

fig24_cap <- "Analysis of two PCR conditions with the \\textsl{MFIaggr} function. 
The VIMCFX96\\_60 and VIMCFX96\\_69 datasets (96-well plate cycler (Bio--Rad 
CFX96)) were used. A subset of cycles (ROI: 2 -- 10) was analyzed for a qPCR 
measured during the annealing phase and elongation phase. The density 
plot (right upper panel) and quantile-quantile analysis (right lower panel) show 
no normal distribution. The measurements during the annealing phase and 
elongation phase show a similar distribution but differ in their mean fluorescence intensity (MFI) levels."

fig24_scap <- "Analysis of two PCR condition with the \\textsl{MFIaggr} function."

fig25_cap <- "Calculation of the amplification efficiency (AE) from a qPCR 
experiment. DNA of \\textsl{vimentin} (\\emph{A}) and \\textsl{MLC-2v} (\\emph{B}) was diluted and 
amplified in a Roche Light Cycler 1.5 (C60.amp dataset). Cqs ($SDM$) were 
calculated by the \\textsl{inder} function and analyzed with \\textsl{effcalc}. 
The AE was (\\emph{C}) 95.1~\\% for \\textsl{vimentin} and  
(\\emph{D}) 94.1~\\%  for \\textsl{MLC-2v}."

fig25_scap <- "Calculation of the amplification efficiency for the C60.amp dataset."


## ----load_data,message=FALSE,results='asis'------------------------------
# Load chipPCR 
require(chipPCR)
# Load package for table formatting
require(xtable)
# Print table
print(xtable(head(C60.amp[, 1L:5]), caption = "First five cycles of imported data."))

## ----problems,warning=FALSE,message=FALSE,fig.show='hold',fig.cap=fig21_cap,fig.scap=fig21_scap,warning=FALSE,fig.width = 11, fig.height = 8----
# Use AmpSim to generate amplification curves with 40 cycles
# and different Cq's.
res.pos <- AmpSim(cyc = 1:40, noise = TRUE, b.eff = -12, nnl = 0.02)
res.pos[5, 2] <- res.pos[5, 2] * 6

res.low <- AmpSim(cyc = 1:40, noise = TRUE, b.eff = -20, bl = 0.5, 
		  ampl = 0.58, Cq = 33)
# Add missing value to res.low at cycle 31
res.low[31, 2] <- NA

res.neg <- AmpSim(cyc = 1:40, b.eff = -0.1, bl = 0.05, ampl = 0.4, Cq = 1, 
		  noise = FALSE, nnl = 0.5)
		      
res.pos.CPP <- cbind(1:40, CPP(res.pos[, 1], res.pos[, 2], 
		     bg.outliers = TRUE, smoother = TRUE, method = "smooth", 
		      method.norm = "minm", method.reg = "lmrob")$y)
		      
res.low.NA <- cbind(1:40, CPP(res.low[, 1], res.low[, 2], smoother = TRUE, 
		    method = "smooth", bg.outliers = TRUE, method.norm = "minm", 
		    method.reg = "lmrob")$y)
		      
res.neg.exc <- cbind(1:40, amptester(res.neg[, 2]))

par(mfrow = c(1,2), las = 0, bty = "n", cex.axis = 1.5, cex.lab = 1.5, 
    font = 2, cex.main = 1.8, oma = c(1,1,1,1))
plot(NA, NA, xlim = c(1,40), ylim = c(0, max(res.pos[, 2])), xlab = "Cycle", 
     ylab = "Raw fluorescence")
mtext("A", cex = 2, side = 3, adj = 0, font = 2)

lines(res.pos, lwd = 2)
lines(res.low, col = 2, lwd = 2)
arrows(38, min(res.low[, 2], na.rm = TRUE), 38, max(res.low[, 2], 
       na.rm = TRUE), code=3, lwd=3, angle=90, col="grey")
text(38, max(res.low[, 2], na.rm = TRUE) * 0.7,"SNR", cex=1.2)

arrows(29,0.42,31,0.51, lwd=2)
text(29,0.38, "NA", cex=1.2)

points(res.pos[5, 1], res.pos[5, 2], pch=21, cex=4, lwd=5, col="orange")
text(res.pos[5, 1], res.pos[5, 2] * 1.2, "Outlier", cex=1.2)

lines(res.neg, col = 4, lwd = 2)
text(20, mean(res.neg[, 2]) *0.9, "No amplification", cex=1.2, col = 
"blue")


plot(NA, NA, xlim = c(1,40), ylim = c(0, max(res.pos[, 2])), xlab = "Cycle", 
     ylab = "Pre-processed fluorescence")
abline(h = 0.03, lty = 2, lwd = 2)
mtext("B", cex = 2, side = 3, adj = 0, font = 2)

lines(res.pos.CPP, lwd = 2)
lines(res.low.NA, col = 2, lwd = 2)
lines(res.neg.exc, col = 4, lwd = 2)

legend(1, 1, c("Positive (outlier removed)", "Positive (scaled)", 
       "Negative", "Threshold line nof Cq"), 
       col = c("black", "red", "blue", "black"), lty = c(1,1,1,2), 
       lwd = 2, bty = "n")
lines(c(15.1,15.1), c(-1,0.03), lwd = 2, col = "black")
text(14, 0.06, "Cq")
lines(c(28.5,28.5), c(-1,0.03), lwd = 2, col = "red")
text(27, 0.06, "Cq", col = "red")

## ----AmpSim_GUI,eval = FALSE---------------------------------------------
#  # Load the shiny package (chipPCR should already be loaded).
#  # Run from a R console following commands.
#  require(shiny)
#  
#  # Invoke the shiny AmpSim app in the default web browser.
#  runApp(paste(find.package("chipPCR")[1],"/AmpSim.gui", sep = ""))
#  
#  # Alternatively call shiny app AmpSim from gist
#  runGist('https://gist.github.com/michbur/e1def41598f1d0c1e2e6')

## ----AmpSim_random,warning=FALSE,message=FALSE,fig.show='hold',fig.cap=fig1_cap,fig.scap=fig1_scap,fig.width = 11, fig.height = 8----
# Draw an empty plot for 40 cycles with user defined parameters.

par(las = 0, bty = "n", oma = c(.5,.5,.5,.5))
plot(NA, NA, xlim = c(1,40), ylim = c(0,1.1), xlab = "Cycle", ylab = "RFU")
colors <- rainbow(8)

# Create eight amplification curves. The approximate Cqs are synthesized 
# as temporary Cqs by adding a random value to a starting Cq of 25. Note: 
# ``noise'' is set TRUE with a level of nnl = 0.03. This adds some scatter 
# to the amplification curves.

sim <- sapply(1L:8, function(i) {
  Cq.tmp <- 25 + rnorm(1) * 5
  
  tmp <- AmpSim(1:40, Cq = Cq.tmp, noise = TRUE, nnl = 0.03)
  lines(tmp, col = colors[i], lwd = 2)
  
  # Add the approximate Cq values to the plot
  text(3, 1 - i / 10, paste("Cq ", round(Cq.tmp, 2)), col = colors[i])
})

## ----humanrater_gui,eval = FALSE-----------------------------------------
#  # Create a set of data to be analyzed by humanrater.
#  # The function AmpSim creates amplification curves which follow a nearly
#  # optimal sigmoidal curve shape or just noise.
#  
#  testdata <- data.frame(1:35,
#                         AmpSim(Cq = 15, noise = TRUE)[, 2],
#                         AmpSim(Cq = 25, noise = TRUE)[, 2],
#                         rnorm(35),
#                         AmpSim(Cq = 35, noise = TRUE)[, 2],
#                         rnorm(35),
#                         AmpSim(Cq = 45, noise = TRUE)[, 2])
#  
#  # Use testdata as input for humanrater and assign the results to the
#  # object human.test.
#  # check testdata for significance of amplification in two repeats.
#  
#  human.test1 <- humanrater(testdata, repeats = 2)

## ----MFIaggr_intro,warning=FALSE,message=FALSE,fig.show='hold', fig.cap=fig23_cap,fig.scap=fig23_scap,fig.width = 11, fig.height = 8----
par(las = 0, bty = "n", cex.axis = 1.2, cex.lab = 1.2, 
    font = 2, cex.main = 1.2, oma = c(1,1,1,1))

plot(MFIaggr(VIMCFX96_60[, 1], VIMCFX96_60[, 2:ncol(VIMCFX96_60)], 
     llul = c(1,10)), CV = FALSE)

# plot(MFIaggr(VIMCFX96_60[, 1], VIMCFX96_60[, 2:ncol(VIMCFX96_60)], 
#      llul = c(1,40)), CV = FALSE)

## ----MFIaggr_conditions,warning=FALSE,message=FALSE,fig.show='hold', fig.cap=fig24_cap,fig.scap=fig24_scap,fig.width = 11, fig.height = 8----
par(las = 0, bty = "n", cex.axis = 1.2, cex.lab = 1.2, 
    font = 2, cex.main = 1.2, oma = c(1,1,1,1))

plot(x = MFIaggr(VIMCFX96_60, fluo = c(2L:10)), y = MFIaggr(VIMCFX96_69, fluo = c(2L:10)))

## ----MFIaggr_all,fig.show='hold',fig.cap=fig12_cap,fig.scap=fig12_scap----
plot(MFIaggr(VIMCFX96_60[, 1], VIMCFX96_60[, 2:ncol(VIMCFX96_60)], 
             llul = c(1,40)), CV = FALSE)

## ----MFIaggr_heteroskedasticity,fig.show='hold',fig.cap=fig15_cap,fig.scap=fig15_scap,fig.width = 11, fig.height = 8----
par(mfrow = c(2,2), bty = "n")
# Create a helper function "hsk.test" to analyze the heteroskedasticity
# and the variance.
hsk.test <- function(x, y, llul = c(1,15), main = "") {
  res <- MFIaggr(x, y, llul = llul)
  head(res)
  plot(res[, 1], res[, 3]^2, xlab = "Cycle", 
       ylab = "Variance of refMFI", 
       xlim = llul, ylim = c(min(res[llul[1]:llul[2], 3]^2), 
                             max(res[llul[1]:llul[2], 3]^2)), main = main, 
       pch = 19, type = "b")
  abline(v = llul, col = "grey", lty = 2, lwd = 2)
  legend("top", paste0("Breusch-Pagan test p-value: \n", 
                       format(summary(res, print = FALSE)[14], digits = 2)), bty = "n")
}

hsk.test(VIMCFX96_60[, 1], VIMCFX96_60[, 2:ncol(VIMCFX96_60)], 
         llul = c(1,15),
         main = "ROI Cycle 1 to 15\nAnnealing phase") 
mtext("A", cex = 2, side = 3, adj = 0)

hsk.test(VIMCFX96_60[, 1], VIMCFX96_60[, 2:ncol(VIMCFX96_60)], llul = 
           c(1,40), main = "ROI Cycle 1 to 40\nAnnealing phase") 
mtext("B", cex = 2, side = 3, adj = 0)

hsk.test(VIMCFX96_69[, 1], VIMCFX96_69[, 2:ncol(VIMCFX96_69)], llul = 
           c(1,15), main = "ROI Cycle 1 to 15\nElongation phase") 
mtext("C", cex = 2, side = 3, adj = 0)

hsk.test(VIMCFX96_69[, 1], VIMCFX96_69[, 2:ncol(VIMCFX96_69)], llul = 
           c(1,40), main = "ROI Cycle 1 to 40\nElongation phase") 
mtext("D", cex = 2, side = 3, adj = 0)

## ----plotCurves,fig.show='hold',fig.cap=fig14_cap,fig.scap=fig14_scap,warning=FALSE,fig.width = 11, fig.height = 8----
y <- VIMCFX96_60[, 2L:9]
# Introduce some missing values.
y[c(10, 22, 3, 25, 26, 15, 27, 23, 4), c(5, 7, 4, 2, 1)] <- NA

# Show plot with raw data and missing values (black line) and show 
# plots with pre-processed data and imputed missing values (red line).
plotCurves(VIMCFX96_60[, 1], y, nrow = 2, type = "l", CPP = TRUE)

## ----workflow,fig.show='hold',fig.cap=fig4_cap,fig.scap=fig4_scap,warning=FALSE,fig.width = 11, fig.height = 8----
layout(matrix(c(1,2,3,3), 2, 2, byrow = TRUE), respect = TRUE)

par(las = 0, bty = "n", oma = c(.5,.5,.5,.5))

th.cyc.raw <- apply(VIMCFX96_60[, -1], 2, function(i) {
  th.cyc(VIMCFX96_60[, 1], i, r = 2575)[1,1]})

res.CPP <- apply(VIMCFX96_60[, -1], 2, function(i) {
  CPP(VIMCFX96_60[, 1], i, trans = TRUE, 
      method.norm = "minm")[["y.norm"]]})

th.cyc.CPP <- apply(res.CPP, 2, function(i) {
  th.cyc(VIMCFX96_60[, 1], i, r = 0.1)[1,1]})

matplot(VIMCFX96_60[, -1], type = "l", pch = 19, col = 1, lty = 1, 
        xlab = "Cycle", ylab = "Raw fluorescence", main = "Raw")
abline(h = 2575, lty = 2)
mtext("A", cex = 1.2, side = 3, adj = 0, font = 2)

matplot(res.CPP, type = "l", pch = 19, col = 1, lty = 1, xlab = "Cycle", 
        ylab = "Fluorescence", main = "CPP")
abline(h = 0.1, lty = 2)
mtext("B", cex = 1.2, side = 3, adj = 0, font = 2)

boxplot(data.frame(Raw = th.cyc.raw, CPP = th.cyc.CPP), ylab = "Cq (Ct)", 
        notch = TRUE)
mtext("C", cex = 1.2, side = 3, adj = 0, font = 2)

## ----fixNA_data,fig.show='hold',fig.cap=fig20_cap,fig.scap=fig20_scap,message=FALSE,fig.width = 11, fig.height = 8----
library(qpcR)
library(chipPCR)
cols <- adjustcolor(2:4, 0.6)
plot(NA, NA, xlim = c(1,45), ylim = c(min(reps384[, -1]), max(reps384[, -1])), 
     col = 1, pch = 19, type = "b", xlab = "Cycle", ylab = "Fluorescence")
rect(0.8, min(reps384[, -1]), 10.2, max(reps384[, -1]), border = NA, col = cols[1])
rect(10.8, min(reps384[, -1]), 33.2, max(reps384[, -1]), border = NA, col = cols[2])
rect(33.8, min(reps384[, -1]), 45, max(reps384[, -1]), border = NA, col = cols[3])

apply(reps384[, -1], 2, function(i) lines(reps384[, 1], i))

## ----fixNA,fig.show='hold',fig.cap=fig8_cap,fig.scap=fig8_scap-----------
# Simulation of an ideal amplification curve with 40 cycles
# The other parameter of the AmpSim function are identical to
# the default.

res <- AmpSim(cyc = 1:40)

# Introduce a missing value (cycle 18) in the transition between 
# the background and the exponential phase.

res.NA <- res
res.NA[18, 2] <- NA

# Helper function to highlight the position of the missing value.
abliner <- function(x1 = 17.5, x2 = 18.5, y1 = 0.09, y2 = 0.14) {
  abline(v = c(x1, x2), col = "red")
  abline(h = c(y1, y2), col = "red")
}

par(las = 0, mfrow = c(2,2), bty = "n")
plot(res, xlab = "Cycles", ylab = "refMFI", type = "b", pch = 20, 
     main = "Without NA")
abliner()
mtext("A", cex = 1.2, side = 3, adj = 0, font = 2)
res.NA.linear <- fixNA(res.NA[, 1], res.NA[, 2], spline = FALSE, 
                       verbose = FALSE)

plot(res.NA, xlab = "Cycles", ylab = "refMFI", type = "b", pch = 20, 
     main = "With NA during transition")
abliner()
mtext("B", cex = 1.2, side = 3, adj = 0, font = 2)

res.NA.spline <- fixNA(res.NA[, 1], res.NA[, 2], spline = TRUE, 
                       verbose = FALSE)

plot(res.NA.linear, xlab = "Cycles", ylab = "refMFI", type = "b", 
     pch = 20, main = "Linear imputed\n NA")
abliner()
mtext("C", cex = 1.2, side = 3, adj = 0, font = 2)

plot(res.NA.spline, xlab = "Cycles", ylab = "refMFI", type = "b", 
     pch = 20, main = "Spline imputed\n NA")
abliner()
mtext("D", cex = 1.2, side = 3, adj = 0, font = 2)
par(mfrow = c(1,1))

## ----fixNA_simulation,echo=TRUE,eval=FALSE-------------------------------
#  # Evaluation of fixNA by introducing missing values and comparing efficiency measures
#  # for imputed and raw data
#  
#  library(qpcR)
#  library(chipPCR)
#  
#  linear_range <- 1L:10
#  exponential_range <- 11L:33
#  plateau_range <- 34L:45
#  
#  raw.eff <- t(sapply(2L:ncol(reps384), function(i) {
#    fit.raw <- pcrfit(reps384, cyc = 1, fluo = i)
#    c(cpD2.raw = efficiency(fit.raw, type = "cpD2", plot = FALSE)[["cpD2"]],
#      Cy0.raw = efficiency(fit.raw, type = "Cy0", plot = FALSE)[["Cy0"]],
#      background = mean(reps384[1L:5, i]))
#  }))
#  
#  #calculate raw efficiency measures
#  raw.df <- as.vector(apply(raw.eff, 2, function(i)
#    c(mean(i), sd(i))))
#  
#  #reformat raw data
#  raw.df <- cbind(data.frame("none", 0, "-"), t(raw.df))
#  colnames(raw.df) <- c("Imputation", "NA.number", "Phase",
#                        "cpD2.mean", "cpD2.sd",
#                        "Cy0.mean", "Cy0.sd",
#                        "background.mean", "background.sd")
#  
#  #helper function introducing NA into data, fixing them and calculating efficiency
#  gen.fix <- function(number_points, phase, spline){
#    #copy raw data
#    temp_dat <- reps384[, -1]
#  
#    #introduce NA(s)
#    for(i in 1L:ncol(temp_dat))
#      temp_dat[sample(phase, 1), i] <- NA
#  
#    #impute missing values using fixNA
#    fixed_dat <- sapply(1L:ncol(temp_dat), function(i)
#      fixNA(reps384[, 1], temp_dat[, i], spline = spline))
#  
#    #compute efficiency measures
#    t(sapply(2L:ncol(fixed_dat), function(i) {
#      fit.fix <- pcrfit(cbind(reps384[, 1], fixed_dat), cyc = 1, fluo = i)
#      c(cpD2.fix = efficiency(fit.fix, type = "cpD2", plot = FALSE)[["cpD2"]],
#        Cy0.fix = efficiency(fit.fix, type = "Cy0", plot = FALSE)[["Cy0"]],
#        background = mean(fixed_dat[1L:5, i]))
#    }))
#  }
#  
#  #calculate results for linear imputation
#  res.linear <- lapply(c(1, 3), function(number_of_points)
#    lapply(list(linear_range, exponential_range, plateau_range), function(phase)
#      apply(gen.fix(number_of_points, phase, FALSE), 2, function(i)
#        c(mean(i), sd(i)))))
#  linear.df <- data.frame(num = unlist(lapply(c(1, 3), rep, 3)),
#                          region = rep(c("linear", "exponential", "plateau"), 2),
#                          do.call(rbind, lapply(unlist(res.linear, recursive = FALSE), as.vector)))
#  
#  
#  #calculate results for spline imputation
#  res.spline <- lapply(c(1, 3), function(number_of_points)
#    lapply(list(linear_range, exponential_range, plateau_range), function(phase)
#      apply(gen.fix(number_of_points, phase, TRUE), 2, function(i)
#        c(mean(i), sd(i)))))
#  spline.df <- data.frame(num = unlist(lapply(c(1, 3), rep, 3)),
#                          region = rep(c("linear", "exponential", "plateau"), 2),
#                          do.call(rbind, lapply(unlist(res.linear, recursive = FALSE), as.vector)))
#  
#  #join and format results
#  sim.res <- cbind(unlist(lapply(c("linear", "spline"), rep, 6)),
#                   rbind(linear.df, spline.df))
#  colnames(sim.res) <- c("Imputation", "NA.number", "Phase",
#                         "cpD2.mean", "cpD2.sd",
#                         "Cy0.mean", "Cy0.sd",
#                         "background.mean", "background.sd")
#  
#  #join simulation results with results for raw data
#  fixNA.evaluation <- rbind(raw.df, sim.res)
#  xtable(fixNA.evaluation[, c(1L:7)], digit = 4)
#  
#  xtable(fixNA.evaluation[, c(1L:3, 8L:9)], digit = 4)

## ----smoothing_intro,fig.show='hold',fig.cap=fig22_cap,fig.scap=fig22_scap,warning=FALSE,fig.width = 11, fig.height = 8----
# Simulate and amplification curve with the AmpSim function
tmp <- AmpSim(cyc = 1:35, bl = 0)

par(las = 0, bty = "n", cex.axis = 1.5, cex.lab = 1.5, 
    font = 2, cex.main = 1.8, oma = c(1,1,1,1), fig = c(0,1,0.55,1))
plot(tmp, type = "b", col = 1, pch = 20, xlab = "", ylab = "RFU", 
      main = "Raw data")
mtext("A", cex = 2, side = 3, adj = 0, font = 2)

# Apply all (parameter method = "all") smoothers/filter with the default 
# setting to the amplification curve of the object tmp. Smoothers / Filters:
# Savitzky-Golay smoothing filter
# locally-weighted polynomial regression
# moving average, windowsize 3
# cubic spline smooth
# standard cubic spline smooth
# Friedman's SuperSmoother
# weighted Whittaker smoothing with first order finite difference penalty
# weighted Whittaker smoothing with a second order finite difference penalty
res <- smoother(tmp[, 1], tmp[, 2], method = "all", CPP = FALSE)

# Calculate the difference between the ideal curve (tmp) and the smoothed curves
# (res) and assign the results to the object res.out
res.out <- cbind(cycle = tmp[, 1], tmp[, 2] - res)

# Plot the smoothed curves
par(fig = c(0,1,0,0.65), new = TRUE)
plot(NA, NA, type = "b", col = 2, pch = 15, xlim = c(1,35), ylim = c(-0.1,0.1),
     xlab = "Cycle", ylab = "delta refMFI (raw - smoothed)",
     main = "Smoothed / Filtered data")
     
mtext("B", cex = 2, side = 3, adj = 0, font = 2) 
legend(1.5, 0.1, ncol = 2, colnames(res.out[, 2:9]), pch = 15:22, 
       lwd = 2, col = c(2:9))

# Plot the results.
tmp <- sapply(2:9, function(i) {
      lines(res.out[, 1], res.out[, i], type = "b", col = i, pch = i + 13)
     }
)

## ----bgmax,fig.cap=fig18_cap,fig.scap=fig18_scap,fig.show='hold',message=FALSE,results='hide',fig.width = 11, fig.height = 8----
par(las = 0, mfrow = c(2,1), bty = "n", oma = c(.5,.5,.5,.5))

res <- AmpSim(cyc = 1:40, Cq = 25)
plot(res, xlim = c(1,40), ylim = c(-0.1,1), xlab = "Cycles", 
     ylab = "refMFI", 
     main = "Background Range Estimation\n in the Absence of Noise", 
     type = "b", pch = 20)
background <- bg.max(res[, 1], res[, 2])
mtext("A", cex = 2, side = 3, adj = 0, font = 2)

points(background[, 3], col = "red", type = "b", pch = 20)
points(background[, 4], col = "blue", type = "b", pch = 20)
abline(v = background@bg.start)
text(background@bg.start, 0.2, "Background start", pos = 4)
abline(v = background@bg.stop, col = "blue")
text(background@bg.stop, 0.25, "Background stop", pos = 4, 
     col = "blue")
abline(v = background@amp.stop, col = "green")
text(background@amp.stop, 0.3, "Plateau transition", pos = 4, 
     col = "green")
legend(4, 1, c("Raw data", "First derivative", "Second derivative"), 
       pch = rep(20,3), col = c(1,2,4), bty = "n")

res <- AmpSim(cyc = 1:40, Cq = 25, noise = TRUE)
plot(res, xlim = c(1,40), ylim = c(-0.1,1), xlab = "Cycles", 
     ylab = "refMFI", 
     main = "Background Range Estimation\n in the Presence of Noise", 
     type = "b", pch = 20)
mtext("B", cex = 2, side = 3, adj = 0, font = 2)
background <- bg.max(res[, 1], res[, 2])

points(background[, 3], col = "red", type = "b", pch = 20)
points(background[, 4], col = "blue", type = "b", pch = 20)
abline(v = background@bg.start)
text(background@bg.start, 0.2, "Background start", pos = 4)
abline(v = background@bg.stop, col = "blue")
text(background@bg.stop, 0.25, "Background stop", pos = 4, col = "blue")
abline(v = background@amp.stop, col = "green")
text(background@amp.stop, 0.3, "Plateau transition", pos = 4, col = 
       "green")
legend(4, 1, c("Raw data", "First derivative", "Second derivative"), 
       pch = rep(20,3), col = c(1,2,4), bty = "n")
par(mfrow = c(1,1))

## ----bgmax_ccPCR,fig.cap=fig19_cap,fig.scap=fig19_scap,fig.show='hold',message=FALSE,results='hide',warning=FALSE,fig.width = 11, fig.height = 8----
# Set parameter for the plot.
par(mfrow = c(2,1), las = 0, bty = "n")

# Use of bg.max for time-dependent measurements. Amplification curves 
# from the capillaryPCR dataset were processed in a loop. The results of 
#  bg.max are added to the plot. 

colors <- rainbow(8)

plot(NA, NA, xlim = c(0,75), ylim = c(-200,1300), xlab = "Time (min)", 
     ylab = "Voltage (micro V)", main = "ccPCR - Raw data")
mtext("A", cex = 1.5, side = 3, adj = 0)
for (i in c(1,3,5,7)) {
  x <- capillaryPCR[1L:750, i]
  y <- capillaryPCR[1:750, i + 1]
  res.bg <- summary(bg.max(x, y))
  lines(x, y, type = "b", pch = 20, col = colors[i], cex = 0.5)
  lines(c(res.bg[2], res.bg[2], res.bg[4], res.bg[4]), 
        c(-150, -50, -150, -50), col = colors[i], lwd = 1.5)
  text(10, 1200 - i * 50, 
       paste("bg.start: ", res.bg[1], ", bg.stop: ", res.bg[2], 
             ", amp.stop: ", res.bg[4]), col = colors[i], cex = 0.6)
}

plot(NA, NA, xlim = c(0,75), ylim = c(-200,1300), xlab = "Time (min)", 
     ylab = "Voltage (micro V)", main = "ccPCR - Pre-processed")
mtext("B", cex = 1.5, side = 3, adj = 0)
for (i in c(1,3,5,7)) {
  x <- capillaryPCR[1L:750, i]
  y <- CPP(capillaryPCR[1L:750, i], capillaryPCR[1:750, i + 1], 
           method = "mova", trans = TRUE, bg.range = c(1,105), 
           bg.outliers = TRUE)[["y.norm"]]
  res.bg <- summary(bg.max(x, y))
  lines(x, y, type = "b", pch = 20, col = colors[i], cex = 0.5)
  lines(c(res.bg[2], res.bg[2], res.bg[4], res.bg[4]), 
        c(-150, -50, -150, -50), col = colors[i], lwd = 1.5)
  text(10, 1200 - i * 50, 
       paste("bg.start: ", res.bg[1], ", bg.stop: ", res.bg[2], 
             ", amp.stop: ", res.bg[4]), col = colors[i], cex = 0.6)
}

## ----normalization,fig.show='hold',fig.cap=fig5_cap,fig.scap=fig5_scap,fig.width = 11, fig.height = 8----
par(mfrow = c(2,3), las = 0, bty = "n", oma = c(.5,.5,.5,.5))
tmp <- VIMCFX96_60

plot(NA, NA, xlim = c(1,40), ylim = c(0, 6000), xlab = "Cycle", 
     ylab = "RFU", main = "Raw data")
mtext("A", cex = 1.2, side = 3, adj = 0, font = 2) 
lin <- apply(tmp[, -1], 2, function(x) lines(tmp[, 1], x))
abline(lm(rowMeans(tmp[2:10, 2L:ncol(tmp)]) ~ tmp[2:10, 1]), col = 2)

plot(NA, NA, xlim = c(1,40), ylim = c(0, 3300), xlab = "Cycle", 
     ylab = "RFU", main = "Baselined data")
mtext("B", cex = 1.2, side = 3, adj = 0, font = 2) 
lin <- apply(tmp[, -1], 2, function(x) lines(tmp[, 1], CPP(tmp[, 1], x, 
                                                           method.norm = "none")$y))


plot(NA, NA, xlim = c(1,40), ylim = c(0, 1.15), xlab = "Cycle", 
     ylab = "RFU", main = "MinMax-Normalization")
mtext("C", cex = 1.2, side = 3, adj = 0, font = 2) 
lin <- apply(tmp[, -1], 2, function(x) lines(tmp[, 1], CPP(tmp[, 1], x, 
                                                           method.norm = "minm")$y))

plot(NA, NA, xlim = c(1,40), ylim = c(0, 1.15), xlab = "Cycle", 
     ylab = "RFU", main = "Max-Normalization")
mtext("D", cex = 1.2, side = 3, adj = 0, font = 2) 
lin <- apply(tmp[, -1], 2, function(x) lines(tmp[, 1], CPP(tmp[, 1], x,, 
                                                           method.norm = "max")$y))

plot(NA, NA, xlim = c(1,40), ylim = c(0, 1.15), xlab = "Cycle", 
     ylab = "RFU", main = "luqn-Normalization")
mtext("E", cex = 1.2, side = 3, adj = 0, font = 2) 
lin <- apply(tmp[, -1], 2, function(x) lines(tmp[, 1], CPP(tmp[, 1], x, 
                                                           method.norm = "luqn", qnL = 0.03)$y))

plot(NA, NA, xlim = c(1,40), ylim = c(-1.5, 1.5), xlab = "Cycle", 
     ylab = "RFU", main = "zscore-Normalization")
mtext("F", cex = 1.2, side = 3, adj = 0, font = 2) 
lin <- apply(tmp[, -1], 2, function(x) lines(tmp[, 1], CPP(tmp[, 1], x, 
                                                           method.norm = "zscore")$y))

## ----lmcoef,fig.show='hold',fig.cap=fig17_cap,fig.scap=fig17_scap,fig.width = 11, fig.height = 8----
par(bty = "n")
plot(VIMCFX96_69[, 1], VIMCFX96_69[, 2], type = "l", xlab = "Cycle", 
     ylab = "Fluorescence")
rect(1,0,10,5000)
method <- c("lmrob", "rq", "least", "rfit")
for (i in 1:4) {
  tmp <- lm.coefs(VIMCFX96_69[1:10, 1], VIMCFX96_69[1:10, 2], 
                  method.reg = method[i])
  text(9, 3200 - i * 100, paste(method[i], ":", "m: ", 
                                round(tmp[1,1], 4), "n: ", round(tmp[2,1], 3)))
  abline(a = tmp[1, 1], b = tmp[2, 1], col = i + 1, lwd = 1.5)
}
legend("right", c("Data", "lmrob", "rq", "least", "rfit"), lty = 1, 
       col = 1:5, cex = 0.95)

## ----rounder_show,warning=FALSE,message=FALSE,fig.show='hold',fig.cap=fig16_cap,fig.scap=fig16_scap,fig.width = 11, fig.height = 8----
# Simulate an amplification curve with 40 cycles using the AmpSim 
# function.
isPCR <- AmpSim(cyc = 1:40)

# Use inder to calculate the derivatives and assign the results to the 
# object res
res <- inder(isPCR)

# Process res by rounder and assign the results to the object rd
rd <- rounder(res)

# Print details of res and rd. Due to the internal use of interpolating 
# splines in inder are the number of elements in the object res the n-th 
# time of the raw data. In this case 200 virtual instead of 40 real cycles.
head(res)
summary(res)

head(rd)
# summary(rd)

## ----SDM,fig.show='hold',fig.cap=fig9_cap,fig.scap=fig9_scap,fig.width = 11,fig.height = 8----
# Use AmpSim to generate an amplification curve with 40 cycles
# and an approximate Cq of 20 and assign it to the object isPCR.
# isPCR is an object of the class "data.frame".
isPCR <- AmpSim(cyc = 1:40, Cq = 20)

# Invoke the inder function for the object isPCR to interpolate 
# the derivatives of the simulated data as object res. The Nip 
# parameter was set to 5. This leads to smoother curves. res is
# an object of the class "der".
res <- inder(isPCR, Nip = 5)

# Plot the object res and add descriptions to the elements.

par(las = 0, bty = "n", oma = c(.5,.5,.5,.5))

plot(isPCR, xlab = "Cycle", ylab = "RFU", ylim = c(-0.15,1),
     main = "", type = "b", pch = 20, lwd = 2)
colors <- rainbow(4)
# Add graphical elements for the derivatives and the calculated
# Cq values FDM, SDM, SDm and SDC.

lines(res[, "x"], res[, "d1y"], col = "blue", lwd = 2)
lines(res[, "x"], res[, "d2y"], col = "red", lwd = 2)

# Fetch the Cq values from res with the summary function
summ <- summary(res, print = FALSE)

abline(v = summ, col = colors, lwd = 2)
text(15, 0.3, paste("FDM ~ ", round(summ["FDM"], 2)), 
     cex = 1.1, col = colors[1])
text(15, 0.2, paste("SDM ~ ", round(summ["SDM"], 2)), 
     cex = 1.1, col = colors[2])
text(15, - 0.1, paste("SDm ~ ", round(summ["SDm"], 2)), 
     cex = 1.1, col = colors[3])
text(15, 0.7, paste("SDC ~ ", round(summ["SDC"], 2)), 
     cex = 1.1, col = colors[4])

legend(1.1, 0.9, c("Raw data", "First derivative", "Second derivative"), 
       col = c(1,4,2), lty = c(2,1,1), bty = "n")

# Summary of the object res.
summ

## ----inder,fig.show='hold',fig.cap=fig10_cap,fig.scap=fig10_scap,message=FALSE,results='hide'----
# Plot all data from C127EGHP and calculate the SDM (Second Derivative 
# Maximum) values with the diffQ2() function (Note: the inder parameter
# is set as TRUE)
# first plot the samples detected with EvaGreen and next the samples 
# detected with the Hydrolysis probe
require(MBmca)

pointer <- function (x, pos = 1, w = 5, stat = TRUE){
  xx <- pos + rep(seq(-0.1, 0.1, length.out = w), ceiling(length(x)/w))
  yy <- sort(x)
  points(xx[1:length(yy)], yy, pch = 19)
  
  if (stat == TRUE)
    x.median <- median(x, na.rm = T)
  x.mad <- mad(x, na.rm = T) * 2
  param <- c(length= 0, code = 3, pch = 15, cex = 1.2)
  arrows(xx[1] * 0.98, x.median, tail(xx, 1) * 1.02, 
         x.median, param, lwd = 3, col = 2)
  arrows(xx[1] * 1.01, x.median + x.mad, tail(xx, 1) * 0.99, 
         x.median + x.mad, param, lwd = 2, lty = 2, col = 4)
  arrows(xx[1] * 1.01, x.median - x.mad, tail(xx, 1) * 0.99, 
         x.median - x.mad, param, lwd = 2, lty = 2, col = 4)
}

amp.liner <- function(range, input, colors = "black") {
  sapply(range, function(i) {
    lines(input[, 2], input[, i], col = colors, pch = 19)
    tmpP <- mcaSmoother(input[, 2], input[, i])
    SDM <- diffQ2(tmpP, inder = TRUE)[["xTm1.2.D2"]][1]
    abline(v = SDM)
    SDM
  }
  )
}

layout(matrix(c(1,3,2,3), 2, 2, byrow = TRUE), respect = TRUE)
par(las = 0, bty = "n")
plot(NA, NA, xlim = c(1,40), ylim = c(0,10), xlab = "Cycle", 
     ylab = "Fluorescence", main = "EvaGreen")
mtext("A", cex = 1.1, side = 3, adj = 0, font = 2)

EG <- amp.liner(range = 3L:34, input = C127EGHP)

plot(NA, NA, xlim = c(1,40), ylim = c(0,10), xlab = "Cycle", 
     ylab = "Fluorescence", main = "Hydrolysis probe")
mtext("B", cex = 1.1, side = 3, adj = 0, font = 2)

HP <- amp.liner(range = 35L:66, input = C127EGHP)

plot(NA, NA, xlim = c(0.8,2.2), ylim = c(13,14), xaxt = "n", 
     xlab = "", ylab = "Cq (SDM, diffQ2)")
text(c(1.05,2), c(13.05,13.05), c("EG", "HP"), cex = 1.2)
mtext("C", cex = 1.1, side = 3, adj = 0, font = 2)
pointer(EG, pos = 1, w = 8)
pointer(HP, pos = 2, w = 8)

## ----inder_fit,fig.cap=fig11_cap,fig.scap=fig11_scap,fig.show='hold',message=FALSE,fig.width = 11,fig.height = 8,results='hide'----
fit.amp <- function(cyc, fluo, plot = FALSE) {
  
  ampl <- quantile(fluo, 0.999)
  bl <- quantile(fluo, 0.001)
  Cq <- round(mean(cyc))
  b.eff <- 1
  
  fit <- nls(fluo ~ bl + ampl / (1 + exp(- (cyc - Cq) / b.eff)), 
             start = list(Cq = Cq, b.eff = b.eff, ampl = ampl, 
                          bl = bl)
  )
  
  res.pred <- data.frame(cyc, predict(fit))
  res <- inder(res.pred[, 1], res.pred[, 2])
  if (plot) {
    lines(res[, 1], res[, 4])
  }
  # SDM
  summary(res)[2]
}

tmp <- C126EG595

out <- apply(tmp[, -1], 2, function(x) fit.amp(tmp[, 1], x))

layout(matrix(c(1,2,1,3), 2, 2, byrow = TRUE))

plot(NA, NA, xlim = c(1,40), ylim = c(min(tmp[, 2L:97]), 
                                      max(tmp[, 2L:97])), 
                                      xlab = "Cycle", 
                                      ylab = "Raw fluorescence")
mtext("A", cex = 1.2, side = 3, adj = 0, font = 2)
for (i in 2L:97) {
  lines(tmp[, 1], tmp[, i], col = ifelse(out[i - 1] < 15.5, "red", 
                                         "black"), lwd = 2)
}
abline(v = out)

plot(NA, NA, xlab = "Cycle", ylab = "RFU''(Cycle)", main = "", 
     xlim = c(0,40), ylim = c(-850, 850))
abline(v = 15.5, lty = 2)
invisible(apply(tmp[, -1], 2, function(x) {
  fit.amp(tmp[, 1], x, plot = TRUE)
}
))
mtext("B", cex = 1.2, side = 3, adj = 0, font = 2)

hist(out, xlab = "Cq (SDM)", main = "", 
     breaks = seq(14.8, 15.8, 0.05), col = rainbow(96))
abline(v = 15.5, lty = 2)
mtext("C", cex = 1.2, side = 3, adj = 0, font = 2)

## ----thcyc,warning=FALSE,message=FALSE,fig.show='hold',fig.cap=fig2_cap,fig.scap=fig2_scap,fig.width = 11, fig.height = 8----
# Raw data from the VIMCFX96_69 dataset.
# Cycles x and Fluoresce values y
x <- VIMCFX96_69[, 1]
y <- VIMCFX96_69[, 2]

par(mfrow = c(2,1), las = 0, bty = "n")

# Plot the raw data
plot(x, y, xlab = "Cycle", ylab = "Fluo", main = "Linear regression", 
     pch = 19)
mtext("A", cex = 1.3, side = 3, adj = 0) 
# Calculate the Cq (Ct) value
res <- th.cyc(x, y, r = 2400, linear = TRUE)
lines(res@input, col = 2, lwd = 2)

# Threshold fluorescence value
abline(h = res[2], col = 3)

# Calculated Ct value
abline(v = res[1], col = 4)
legend("topleft", paste("Cq (Ct) = ", round(res[1], 3)))

plot(x, y, xlab = "Cycle", ylab = "Fluo", main = "Quadratic regression", 
     pch = 19)
mtext("B", cex = 1.3, side = 3, adj = 0) 

# Calculate the Ct value
res <- th.cyc(x, y, r = 2400, linear = FALSE)
lines(res@input, col = 2, lwd = 2)

# Threshold fluorescence value
abline(h = res[2], col = 3)

# Calculated Ct value
abline(v = res[1], col = 4)
legend("topleft", paste("Cq (Ct) = ", round(res[1], 3)))

## ----thcyc_ccPCR,fig.show='hold',fig.cap=fig3_cap,fig.scap=fig3_scap,fig.width = 11, fig.height = 8----
# Application of the th.cyc method to determine the Cq from a continuous
# amplification reaction.
par(las = 0, bty = "n", oma = c(.5,.5,.5,.5))
       
plot(NA, NA, xlim = c(0,80), ylim = c(0,1200), xlab = "Time (min)", 
     ylab = "Voltage [micro V]", main = "ccPCR - Pre-processed Data")
mtext("B", cex = 2, side = 3, adj = 0)
# Threshold level "r" (50 micro Volts)
for (i in c(1,3,5,7)) {
  y.tmp <- CPP(capillaryPCR[, i], capillaryPCR[, i + 1], trans = TRUE, bg.range = c(1,150))$y.norm
  Ct.tmp <- th.cyc(capillaryPCR[, i], y.tmp, r = 80, linear = FALSE)
  abline(v = Ct.tmp[1])
  text(Ct.tmp[1] * 1.125, 1200, paste(round(Ct.tmp[1], 1), "\nmin"), cex = 0.8)
  lines(capillaryPCR[, i], y.tmp, type = "b", pch = 20 - i) 
  points(Ct.tmp@input, col = "red", pch = 19)
}
abline(h = 80)
legend("topleft", c("Run 1", "Run 2", "Run 3", "Control"), 
       pch = c(19, 17, 15, 13), lwd = 1.3, bty = "n")

## ----HDA,fig.show='hold',fig.cap=fig13_cap,fig.scap=fig13_scap,message=FALSE,warning=FALSE,fig.width = 11, fig.height = 8----
par(mfrow = c(2,1), bty = "n")
plot(NA, NA, xlim = c(0,5000), ylim = c(0.45,0.8), xlab = "Time (sec)", 
     ylab = "Fluorescence", main = "HDA - Raw data")
mtext("A", cex = 2, side = 3, adj = 0)
lines(C85[, 2], C85[, 3], type = "b", col = 2, pch = 20)
lines(C85[, 4], C85[, 5], type = "b", col = 4, pch = 20)
lines(C85[, 6], C85[, 7], type = "b", col = 6, pch = 20)
legend("bottomright", c("D1, 1x", "D2, 1:10", "D3, 1:100"), col = c(2,4,6), 
       pch = rep(20,3), bty = "n")

plot(NA, NA, xlim = c(0,2000), ylim = c(0,0.4), xlab = "Time (sec)", 
     ylab = "Fluorescence", main = "HDA - Pre-processed data")
mtext("B", cex = 2, side = 3, adj = 0)
legend("topleft", c("D1, 1x", "D2, 1:10", "D3, 1:100"), col = c(2,4,6), 
       pch = rep(20,3), bty = "n")

# Define the parameters for the pre-processing by CPP and the th.cyc 
# function.
# smoothing method
sm <- "mova"

# manual range for background
br <- c(2,10)

# time range for analysis
xr <- 3L:200

# method for baseline normalization
lrg <- "least"

# threshold level for the th.cyc function
r <- 0.025
# Calculate in a loop the Cq values (Cycle threshold method) and add the
# calculated time (in minutes) to the plot.
for (i in c(2,4,6)) {
  y.tmp <- CPP(C85[xr, i], C85[xr, i + 1], method = sm, bg.range = br, 
               trans = TRUE)$y.norm
  Ct.tmp <- th.cyc(C85[xr, i], y.tmp, r = r, linear = FALSE)
  abline(v = Ct.tmp[1], col = "grey")
  lines(C85[xr, i], y.tmp, col = i, lwd = 2)
  points(Ct.tmp@input, col = "red", pch = 19)
  text(Ct.tmp[1] * 1.1, 0.36, paste(round(Ct.tmp[1]/60, 1), "\nmin"))
}

# Show the fluorescence value, which defines the threshold.
abline(h = r, lty = 2)

## ----AmpSim_effcalc,fig.show='hold',fig.cap=fig6_cap,fig.scap=fig6_scap,fig.width = 11,fig.height = 8,message=FALSE----

# Load MBmca package (v. 0.0.3-3 or later)
require(MBmca)

# Create an graphic device for two empty plots.
par(mfrow = c(1,2))
plot(NA, NA, xlim = c(1,45), ylim = c(0.01,1.1), xlab = "Cycles", 
     ylab = "Fluorescence", main = "")
mtext("A", cex = 1.1, side = 3, adj = 0, font = 2)

# Create a sequence of "targeted" Cq values (Cq.t) between 15 and 34 
# cycles.

Cq.t <- rep(seq(15, 34, 3.5), 3)

# In-silico experiment set up: Define the levels for the decadic dilutions
# with concentrations from 100 to 0.001 (six steps) as three replicates.

dilution <- rep(10^(2:-4), 3)

# Create an empty matrix for the results of the concentration
# dependent Cq values.

ma.out <- matrix(data = NA, nrow = 45, ncol = length(Cq.t))

# Use AmpSim to simulate amplification curves at different concentrations. 
# The simulation is performed with the addition of some noise. This 
# generates unique (non-reproducible) amplification curves, even under 
# identical parameter settings.

Cq.out <- vector()

# Simulate a qPCR reaction with AmpSim for 45 cycles and some noise.

for (i in 1L:18) {
  ma.out[1:45, i] <- AmpSim(cyc = c(1:45), b.eff = -50, bl = 0.001, 
                            ampl = 1, Cq = Cq.t[i], noise = TRUE, 
                            nnl = 0.02)[, 2]
  lines(1:45, ma.out[, i])
  tmpP <- mcaSmoother(1:45, ma.out[, i])
  # Calculate the pseudo Second Derivative Maximum (SDM) (Cq) using 
  # the diffQ2 function from the MBmca package.
  Cq.tmp <- diffQ2(tmpP, inder = TRUE)$xTm1.2.D2[1]
  abline(v = Cq.tmp)
  Cq.out <- c(Cq.out, Cq.tmp)
}

# Assign the calculated Cqs to the corresponding concentrations.
tmp <- data.frame(dilution[1:6], Cq.out[1:6], Cq.out[7:12], Cq.out[13:18])

# Determine the amplification efficiency by using the effcalc function.
plot(effcalc(tmp[, 1], tmp[, 2:4]), CI = TRUE)
mtext("B", cex = 1.1, side = 3, adj = 0, font = 2) 

## ----CPP_C54,fig.show='hold',fig.cap=fig7_cap,fig.scap=fig7_scap,fig.width = 11,fig.height = 8,message=FALSE----
require(MBmca)
par(las = 0, bty = "n", oma = c(.5,.5,.5,.5))
par(fig = c(0,0.5,0,1), new = TRUE)
plot(NA, NA, xlim = c(1,55), ylim = c(0, 0.7), xlab = "Cycle", 
     ylab = "refMFI", main = "Raw data")
just_line <- apply(C54[, c(2:4)], 2, function(y) lines(C54[, 1], y))
mtext("A", cex = 1.2, side = 3, adj = 0, font = 2) 

par(fig = c(0.5,1,0.5,1), new = TRUE)
plot(NA, NA, xlim = c(1,55), ylim = c(0, 0.55), xlab = "Cycle", 
     ylab = "refMFI", main = "pre-processed data")
mtext("B", cex = 1.2, side = 3, adj = 0, font = 2) 

D1 <- cbind(C54[1:35, 1], CPP(C54[1:35, 1], C54[1:35, 2], trans = TRUE, 
                              bg.range = c(1,8))[["y.norm"]])
D2 <- cbind(C54[1:45, 1], CPP(C54[1:45, 1], C54[1:45, 3], trans = 
                                TRUE)[["y.norm"]])
D3 <- cbind(C54[1:55, 1], CPP(C54[1:55, 1], C54[1:55, 4], trans = 
                                TRUE)[["y.norm"]])

lines(D1, col = 1)
lines(D2, col = 2)
lines(D3, col = 3)

dilution <- c(1E0, 1E-3, 1E-6)
Cq.D1 <- diffQ2(D1, inder = TRUE)[["xTm1.2.D2"]][1]
Cq.D2 <- diffQ2(D2, inder = TRUE)[["xTm1.2.D2"]][1]
Cq.D3 <- diffQ2(D3, inder = TRUE)[["xTm1.2.D2"]][1]

res.dil <- data.frame(dilution, rbind(Cq.D1, Cq.D2, Cq.D3))
par(fig = c(0.5,1,0,0.5), new = TRUE)
plot(effcalc(res.dil[, 1], res.dil[, 2]), res.fit = NULL)

## ----effcalc_output,echo=FALSE,results='asis',message=FALSE--------------
print(xtable(effcalc(res.dil[, 1], res.dil[, 2]), 
      caption = "Output of the effcalc function.", label = "table:effcalc_output"))

## ----effcalc_VIM_MLC,fig.show='hold',fig.cap=fig25_cap,fig.scap=fig25_scap,fig.width = 11, fig.height = 8----
colors <- rep(rainbow(7), each = 2)
par(mfrow = c(2,2))

plot(NA, NA, xlim = c(0,44), ylim = c(0, 6), xlab = "Cycles", ylab = "RFU")
legend(0, 6, colnames(C60.amp[, 4L:17]), ncol = 2, col = colors[1:14], 
       pch = 19, bty = "n")
mtext("A", cex = 1.2, side = 3, adj = 0, font = 2)
SDM.vim <- sapply(4L:17, function(i) {
  lines(C60.amp[, 1], C60.amp[, i], col = colors[i - 3])
  SDM <- summary(inder(C60.amp[, 1], C60.amp[, i]), print = FALSE)[2]
}
)

plot(NA, NA, xlim = c(0,44), ylim = c(0, 4), xlab = "Cycles", ylab = "RFU")
legend(0, 4, colnames(C60.amp[, 18L:31]), ncol = 2, col = colors[1:14], 
       pch = 19, bty = "n")
mtext("B", cex = 1.2, side = 3, adj = 0, font = 2)
SDM.mlc2v <- sapply(18L:31, function(i) {
  lines(C60.amp[, 1], C60.amp[, i], col = colors[i - 17])
  SDM <- summary(inder(C60.amp[, 1], C60.amp[, i]), print = FALSE)[2]
}
)

#create vector of dilutions
dil <- sort(rep(10^(0L:-6), 2), TRUE)

res <- cbind(dil, SDM.vim, SDM.mlc2v)

plot(effcalc(res[, 1], res[, 2]))
mtext("C", cex = 1.2, side = 3, adj = 0, font = 2)

plot(effcalc(res[, 1], res[, 3]))
mtext("D", cex = 1.2, side = 3, adj = 0, font = 2)

## ----summ.datasets,echo=FALSE,results='asis'-----------------------------
load("datdf.RData")

datdf <- cbind(datdf, t(sapply(datdf[["dat.names"]], function(i) 
  dim(get(substr(i, 1, nchar(i) - 1))))))
tmp <- apply(datdf, 1, function(i) {
  safe.name <- i[1]
  if(grepl("_", safe.name, fixed = TRUE)) {
    safe.name <- sub("_", ".", i, fixed = TRUE)
  }
  cat("\\item Dataset: ", safe.name, "\n\\begin{itemize}\n\\item Dataset type: ", 
      i[3], "\n\\item Description: ", i[2], "\n\\item Number of variables: ", i[4], 
      "\n\\item Number of measurements: ", i[5], "\n\\end{itemize}\n")
})

