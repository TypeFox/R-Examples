### R code from vignette source 'speaq.Rnw'
### Encoding: UTF-8

###################################################
### code chunk number 1: Read_data_input
###################################################
library(speaq)

#Generate a simulated NMR data set for this experiment
res=makeSimulatedData();
X=res$data;
groupLabel=res$label;


###################################################
### code chunk number 2: Unaligned_spectral_plots
###################################################
drawSpec(X);


###################################################
### code chunk number 3: Peak_detection
###################################################
cat("\n detect peaks....");
startTime <- proc.time();
peakList <- detectSpecPeaks(X,
    nDivRange = c(128),                
    scales = seq(1, 16, 2),
    baselineThresh = 50000,
    SNR.Th = -1,
    verbose=FALSE
);

endTime <- proc.time();
cat("Peak detection time:",(endTime[3]-startTime[3])/60," minutes");


###################################################
### code chunk number 4: Reference_finding
###################################################

cat("\n Find the spectrum reference...")
resFindRef<- findRef(peakList);
refInd <- resFindRef$refInd;

cat("\n Order of spectrum for reference  \n");
for (i in 1:length(resFindRef$orderSpec))
{
    cat(paste(i, ":",resFindRef$orderSpec[i],sep=""), " ");
    if (i %% 10 == 0) cat("\n")
}
    
cat("\n The reference is: ", refInd);



###################################################
### code chunk number 5: Spectral_alignment
###################################################
# Set maxShift
maxShift = 50;

Y <- dohCluster(X,
                peakList = peakList,
                refInd = refInd,
                maxShift  = maxShift,
                acceptLostPeak = TRUE, verbose=FALSE);



###################################################
### code chunk number 6: Spectral_segment_alignment
###################################################
segmentInfoMat=matrix(data=c(100,200,0,0,0,
                      450,680,1,0,50),nrow=2,ncol=5,byrow=TRUE
                      )
colnames(segmentInfoMat)=c("begin","end","forAlign","ref","maxShift")
segmentInfoMat

Yc <- dohClusterCustommedSegments(X,
                                 peakList = peakList,
                                 refInd = refInd,
                                 maxShift  = maxShift,
                                 acceptLostPeak = TRUE,
                                 segmentInfoMat = segmentInfoMat,
                                 minSegSize = 128,
                                 verbose=FALSE)
                                 


###################################################
### code chunk number 7: Aligned_spectral_plots
###################################################
drawSpec(Y);


###################################################
### code chunk number 8: Aligned_spectral_plots_limited_height
###################################################
drawSpec(Y,
        startP=450,
        endP=680,
        highBound = 5e+5,
        lowBound = -100);


###################################################
### code chunk number 9: Aligned_spectral_plots_customized
###################################################
drawSpec(Yc);


###################################################
### code chunk number 10: Quantitative_analysis
###################################################
N = 100;
alpha = 0.05;

# find the BW-statistic
BW = BWR(Y, groupLabel);

# create sampled H0 and export to file
H0 = createNullSampling(Y, groupLabel, N = N,verbose=FALSE)

#compute percentile of alpha
perc = double(ncol(Y));
alpha_corr = alpha/sum(returnLocalMaxima(Y[2,])$pkMax>50000);
for (i in 1 : length(perc)){    
    perc[i] = quantile(H0[,i],1-alpha_corr, type = 3);
}


###################################################
### code chunk number 11: drawBW_1
###################################################
drawBW(BW, perc,Y, groupLabel = groupLabel)


###################################################
### code chunk number 12: drawBW_2
###################################################
drawBW(BW, perc, Y ,startP=450, endP=680, groupLabel = groupLabel)


