# Plotting and data analysis functions

FTGramImage <- function(sig, dt, ft, time.span = NULL, freq.span = NULL, amp.span = NULL, blur = NULL, taper = 0.05, scaling = "none", grid=TRUE, colorbar=TRUE, backcol=c(0, 0, 0), colormap=NULL, pretty=FALSE, ...)
{
	#Plots a Fourier spectrogram
	#INPUTS
        #    SIG is the signal to analyze
        #    DT is the sample rate (must be constant)
	#    FT is the Fourier transform input parameters, adopted from Jonathan Lees' code in RSEIS
	#        FT$NFFT is the fft length
	#        FT$NS is the number of samples in a window
        #        FT$NOV is the number of samples to overlap
        #    TIME.SPAN is the time span to plot, NULL plots everything
        #    FREQ.SPAN is the frequency span to plot (<=max frequency in spectrogram), NULL plots everything up to the Nyquist frequency
	#    AMP.SPAN is the amplitude range to plot.  NULL plots everything.
        #    BLUR is a list of parameters for a Gaussian image smoothing kernel, if desired.  If not null then
        #        BLUR$SIGMA - Standard deviation of Gaussian kernel.  If a 2 element vector, then the kernel has independent coordinate
        #        BLUR$BLEED - Whether to allow blur to bleed out of the domain of the image 
        #    TAPER is the cosine taper factor (amount of the signal to apply the taper to, must be < 0.5)
        #    SCALING determines whether to apply a logarithmic (log), or square root (sqrt) scaling to the amplitude data
        #    GRID is a boolean asking whether to display grid lines
        #    COLORBAR is a boolean asking whether to plot an amplitude colorbar
        #    BACKCOL is a 3 element vector of RGB values for the background of the spectrogram, based on a 0 to 255 scale: [red, green, blue]
        #    COLORMAP is an R palette object determining how the spectrogram colors should look
        #    PRETTY is a boolean asking whether to adjust axis labels so that they're pretty (TRUE) or give the exactly specified time and frequency intervals (FALSE)
        #    OPTIONAL PARAMETERS
        #       TRACE.FORMAT is the format of the trace minima and maxima in sprintf format
        #       IMG.X.FORMAT is the format of the X axis labels of the image in sprintf format
        #       IMG.Y.FORMAT is the format of the Y axis labels of the image in sprintf format
        #       COLORBAR.FORMAT is the format of the colorbar labels in sprintf format   
        #       CEX.LAB is the font size of the image axis labels
        #       CEX.COLORBAR is the font size of the colorbar
        #       CEX.TRACE is the font size of the trace axis labels
        #       IMG.X.LAB is the X - axis label of the image, it defaults to "time"
        #       IMG.Y.LAB is the Y - axis label of the image, it defaults to "frequency" 
        #       MAIN gives the figure a title.
        #OUTPUTS
        #    IMG is the spectrogram	
	opts = list(...)

        if(!"img.x.lab" %in% names(opts))
        {
            opts$img.x.lab = "time"
        }

        if(!"img.y.lab" %in% names(opts))
        {
            opts$img.y.lab = "frequency"
        }

        if(is.null(time.span))
        {
                time.span=c(dt, length(sig) * dt)
        }

        if(time.span[2] > length(sig) * dt)
        {
                time.span[2]= length(sig) * dt
                warning("The requested spectrogram is longer than the actual signal.")
        }

        if(is.null(freq.span))
        {
                freq.span=c(0, 1/(dt * 2))
        }
        if(freq.span[2] > 1 / (dt * 2))
        {
                freq.span[2] = 1 / (dt * 2)
                warning("Requested maximum frequency is higher than the Nyquist frequency.")
        }

        sig = sig[(time.span[1]/dt):(time.span[2]/dt)]
        tt = (seq_len(length(sig)) * dt) + time.span[1]

	ev=EvolutiveFFT(sig, dt, ft, freq.span, taper) #Calculate the Fourier spectrogram
        ev$tt = tt

        if(is.null(amp.span))
        {    
             amp.span = c(min(ev$z[ev$z>-Inf]), max(ev$z[ev$z<Inf]))
        }
        img.xvec = ev$x + time.span[1]
        img.yvec = seq(freq.span[1], freq.span[2], by = ev$y[2] - ev$y[1])
        img = list(z = array(0, dim = c(length(img.xvec), length(img.yvec))),
              x = img.xvec, y = img.yvec)
        img$z[,img.yvec >= min(ev$y) & img.yvec <= max(ev$y)] = ev$z[,ev$y >= freq.span[1] & ev$y <= freq.span[2]]

        if(scaling == "ln") #Scale by natural log
        {
            img$z[img$z == 0] = NA
            img$z = log(img$z)
            amp.span <- log(amp.span)
        }

        if(scaling == "log") #Log 10 scale
        {
            img$z[img$z == 0] = NA
            img$z = log10(img$z)
            amp.span <- log10(amp.span)
        }

        if(scaling == "sqrt") #Take the square root
        {
            img$z = sqrt(img$z)
            amp.span <- sqrt(amp.span)
        }

        trace = list()
        trace$sig = ev$original.signal[ev$tt >= time.span[1] & ev$tt <= time.span[2]]
        trace$tt = ev$tt[ev$tt >= time.span[1] & ev$tt <= time.span[2]]
       
        window = ft$ns / (length(tt[tt >= min(img$x) & tt <= max(img$x)]))
        HHTPackagePlotter(img, trace, amp.span, blur = blur, opts$img.x.lab, opts$img.y.lab, window = window, colormap = colormap, backcol = backcol, pretty = pretty, grid = grid, colorbar = colorbar, opts = opts)
 
        invisible(img)

}		

HHRender <- function(hres, dt, dfreq, time.span = NULL, freq.span = NULL, scaling = "none", combine.imfs = TRUE, verbose = TRUE)
{
	#Renders a spectrogram of EMD or Ensemble EMD (EEMD) results.
	#INPUTS
	#	HRES is a matrix of data generated by EEMD.COMPILE or the output of HHTRANSFORM
	#	it represents a set on all time/frequency/amplitude points from the given EEMD run
        #       DT is the time resolution of the spectrogram.  Currently, if there is a hres$dt field, DT must be greater than or equal to hres$dt.
        #       this prevents subsample resolution.
        #       DFREQ is the frequency resolution of the spectrogram
        #       TIME.SPAN is the portion of the signal to include.  NULL means the whole signal.
        #       FREQ.SPAN is the frequency range to calculate the spectrum over c(MIN, MAX).  NULL means capture the full frequency spectrum of the signal.
        #       SCALING determines whether to plot frequency as log 10 ("log") or linear ("none")
        #       COMBINE.IMFS will combine all the IMFs into one image, saving space and time for HHGramImage if TRUE.  If FALSE, keep them separate for individual plotting options for HHGramImage.
        #       VERBOSE prints out status messages (i.e. IMF 1 COMPLETE!)
	#OUTPUTS
	#	HGRAM is a spectrogram matrix ready to be plotted by HHGRAM.IMAGE
	#Danny Bowman
	#UNC Chapel Hill

        hgram = hres

        if(scaling == "log")
        {
            hres$hinstfreq = log10(hres$hinstfreq)
        }
        else if (scaling != "none")
        {
             warning("Did not recognize scaling request \"", scaling, ".\" Reverting to linear frequency (scaling = \"none\").")
        }
       
        #Deal with logarithms of 0
        hres$hamp[hres$hinstfreq == -Inf] = 0
        hres$hinstfreq[hres$hinstfreq == -Inf] = 0

        if(is.null(freq.span))
        {
             freq.span = c(min(hres$hinstfreq), max(hres$hinstfreq))
        }
	
        if(!"trials" %in% names(hres))
 	{
		hres$trials=1
                hres$hinstfreq = array(hres$hinstfreq, dim = c(dim(hres$hinstfreq), 1))
                hres$hamp = array(hres$hamp, dim = c(dim(hres$hamp), 1))
	}
        
        if("dt" %in% names(hres))
        {
             if(hres$dt > dt) #We don't want to have to interpolate between samples
             {
                 warning(paste("The time resolution", sprintf("%.2e", dt), "is lower than the sample rate", sprintf("%.2e", hres$dt), "of the time series.  This may introduce time gaps in the spectrogram."))
             }
             if("tt" %in% names(hres))
             {
                 warning("Input data has both DT (sample rate) and TT (sample times) components.  Component TT will be used to calculate the spectrogram")
                 hgram$tt = hres$tt
             }
             else
             {
                 hgram$tt = seq_len(length(hres$original.signal)) * hres$dt
             }
        }

        if(is.null(time.span))
        {
           time.span = c(min(hgram$tt), max(hgram$tt))
        }

        if(!(("tt" %in% names(hres)) | ("dt" %in% names(hres))))
        {
            warning("Neither DT (sample rate) nor TT (sample times) were specified in the input data.  Assuming DT is 1...")
            hgram$tt = seq_len(length(hres$original.signal))
        } 

        if(time.span[2]>max(hgram$tt))
        {
                time.span[2]=max(hgram$tt)
                warning("Requested time window is longer than the actual signal.")
        }
        t.ind = which(hgram$tt >= time.span[1] & hgram$tt <= time.span[2])
        hgram$tt = hgram$tt[t.ind] 
        hres$hinstfreq = array(hres$hinstfreq[t.ind,,], dim = c(length(hgram$tt), hres$nimf, hres$trials))
        hres$hamp = array(hres$hamp[t.ind,,], dim = c(length(hgram$tt), hres$nimf, hres$trials))
        hres$original.signal = hres$original.signal[t.ind]
       
        grid = list()
        grid$x = hgram$tt
        grid$y = seq(from = freq.span[1], to = freq.span[2] + dfreq, by = dfreq)
        if(combine.imfs)
        {
            imf.dim = 1
        }
        else{
            imf.dim = hres$nimf
        }
        hgram$z=array(0,dim=c(length(grid$x),length(grid$y), imf.dim))
        hgram$cluster=hgram$z #Shows how many times a given grid node has data.
        for(i in seq(hres$nimf))
        {
            x = array(c(rep(hgram$tt,hres$trials), hres$hinstfreq[,i,]), dim = c(length(hgram$tt)*hres$trials, 2))
            imf.img = fields::as.image(hres$hamp[,i,], grid = grid, x = x)
            imf.img$z[is.na(imf.img$z)] = 0
            imf.img$weights[is.na(imf.img$weights)] = 0
            if(combine.imfs)
            {
                hgram$z[,,1] = hgram$z[,,1] + imf.img$z
                hgram$cluster[,,1] = hgram$cluster[,,1] + imf.img$weights
            }
            else{
                hgram$z[,,i] = imf.img$z
                hgram$cluster[,,i] = imf.img$weights
            }
            if(verbose)
            {
                print(paste("IMF", i, "COMPLETE!"))
            }
        }
        hgram$combine.imfs = combine.imfs 
        hgram$hinstfreq = hres$hinstfreq
        hgram$hamp = hres$hamp 
        hgram$original.signal = hres$original.signal
        hgram$x = imf.img$x 
        hgram$y = imf.img$y
	hgram$dfreq=dfreq
	hgram$dt=hres$dt
        hgram$scaling = scaling
	invisible(hgram) #Return the spectrogram structure.
}

HHSpectrum <- function(hres, dfreq, freq.span = NULL, time.span = NULL, scaling = "none", verbose = TRUE)
{
    #Calculate the Hilbert spectrogram of a signal contained in HRES (returned by HHTRANSFORM or EEMD.COMPILE)
    #INPUTS
    #       HRES is a matrix of data generated by EEMD.COMPILE or the output of HHTRANSFORM
    #       it represents a set on all time/frequency/amplitude points from the given EEMD run
    #       DFREQ is the frequency resolution of the spectrogram
    #       FREQ.SPAN is the frequency range to calculate the spectrum over c(MIN, MAX).  NULL means capture the full frequency spectrum of the signal.
    #       TIME.SPAN is the time span to calculate the spectrum over c(MIN, MAX).  NULL means use the entire signal
    #       SCALING determines whether to calculate frequency as log 10 ("log") or linear ("none")
    #       VERBOSE prints out status messages (i.e. IMF 1 COMPLETE!)
    #OUTPUTS
    #    HSPEC is the Hilbert spectrum of the signal, separated by IMF.

   if(is.null(time.span))
   {
       dt = max(hres$tt) - min(hres$tt)
   }
   else {
       dt = time.span[2] - time.span[1]
   }

   hgram = HHRender(hres, dt, dfreq, freq.span = freq.span, time.span = time.span, scaling = scaling, combine.imfs = FALSE, verbose = TRUE)

   amps = array(0, dim = dim(hgram$z)[2:3])

   for(i in seq(hres$nimf))
   {
       amps[, i] = apply(hgram$z[, , i], 2, sum) 
   }

  hspec = list(amplitude = amps, frequency = hgram$y, original.signal = hgram$original.signal, dt = dt, tt=hres$tt, dfreq = dfreq)
  invisible(hspec)
} 

HHGramImage <- function(hgram,time.span = NULL,freq.span = NULL, amp.span = NULL, blur = NULL, clustergram = FALSE, cluster.span=NULL, imf.list = NULL, fit.line = FALSE, scaling = "none", grid=TRUE, colorbar=TRUE, backcol=c(0, 0, 0), colormap=NULL, pretty=FALSE, ...)
{
	#Plots a spectrogram of the EEMD processed signal as an image.	
	#INPUTS
	#	HGRAM is the subsetted spectrogram  from HH.RENDER.
	#		HGRAM$X is time
	#		HGRAM$Y is frequency
	#		HGRAM$Z is amplitude normalized to trials
	#		HGRAM$CLUSTER is a matrix containing integer values corresponding to the number of times a signal was recorded in a given spectrogram cell during EEMD
	#		The more often the signal is recorded, the more likely it is that the signal is real and not noise
	#		HGRAM$TRIALS is the number of times EEMD was run to generate signal
	#		HGRAM$ORIGINAL.SIGNAL is the original seismogram (without added noise)
        #               HGRAM$TT is the sample times
	#	TIME.SPAN is the time span to plot, NULL plots everything
	#	FREQ.SPAN is the frequency span to plot (<=max frequency in spectrogram), NULL plots everything
	#	AMP.SPAN is the amplitude span to plot, everything below is set to black, everything above is set to max color, NULL scales to range in signal
        #       BLUR is a list of parameters for a Gaussian image smoothing kernel, if desired.  If not null then
        #           BLUR$SIGMA - Standard deviation of Gaussian kernel.  If a 2 element vector, then the kernel has independent coordinate
        #           BLUR$BLEED - Whether to allow blur to bleed out of the domain of the image 

        #	CLUSTERGRAM tells the code to plot the signal amplitude (FALSE) or the number of times data occupies a given pixel (TRUE).
	#	CLUSTER.SPAN plots only the parts of the signal that have a certain number of data points per pixel [AT LEAST, AT MOST] this only applies to EEMD with multiple trials.
        #       IMF.LIST is a list of IMFs to plot on the spectrogram.  If NULL, plot all IMFs.
        #       IMF.SUM can be set to show the sum of IMFs shown in the spectrogram plotted as a red line against the original trace
        #       SCALING determines whether to apply a logarithmic (log), or square root (sqrt) scaling to the amplitude data, default is "none"
	#	GRID is a boolean asking whether to display grid lines
	#	COLORBAR is a boolean asking whether to plot an amplitude colorbar
        #       BACKCOL is a 3 element vector of RGB values for the background of the spectrogram, based on a 0 to 255 scale: [red, green, blue]
        #       COLORMAP is an R palette object determining how the spectrogram colors should look
        #       PRETTY is a boolean asking whether to adjust axis labels so that they're pretty (TRUE) or give the exactly specified time and frequency intervals (FALSE)
	#OPTIONAL PARAMETERS
        #       TRACE.FORMAT is the format of the trace minima and maxima in sprintf format
        #       IMG.X.FORMAT is the format of the X axis labels of the image in sprintf format
        #       IMG.Y.FORMAT is the format of the Y axis labels of the image in sprintf format
        #       COLORBAR.FORMAT is the format of the colorbar labels in sprintf format   
        #       CEX.LAB is the font size of the image axis labels
        #       CEX.COLORBAR is the font size of the colorbar
        #       CEX.TRACE is the font size of the trace axis labels
        #       IMG.X.LAB is the X - axis label of the image, it defaults to "time"
        #       IMG.Y.LAB is the Y - axis label of the image, it defaults to "frequency"
        #OUTPUTS
        #     IMG is the spectrogram returned as an image

        opts = list(...)
 
        if(!"img.x.lab" %in% names(opts))
        {
            opts$img.x.lab = "time"
        }
        
        if(!"img.y.lab" %in% names(opts))
        {
            opts$img.y.lab = "frequency"
        }
  
        #Subset by IMFs
        if(is.null(imf.list))
        {
            if(hgram$combine.imfs)
            {
                imf.list = seq(1)
            }
            else{

                imf.list = seq(hgram$nimf)
            }
        }
        else
        {
            if(hgram$combine.imfs)
            {
                warning("The IMFs were combined when HHRender was run on this data (combine.imfs = TRUE). Individual IMF spectrograms cannot be plotted - the image you see is the combined IMFs.  Rerun HHRender with combined.imfs = FALSE if you want the ability to plot single IMFs using HHGramImage.")
                imf.list = seq(1)
            }
            if(max(imf.list) > hgram$nimf)
            {
                warning("Requested more IMFs than are present in the actual EMD results!")
                imf.list = imf.list[imf.list < hgram$nimf]
            }
        }   
   
	if(is.null(time.span))
	{
		time.span=c(min(hgram$tt), max(hgram$tt))
	}
	
	if(time.span[2]>max(hgram$tt))
	{
		time.span[2]=max(hgram$tt)
		warning("Requested time window is longer than the actual signal.")
	}
	
	if(is.null(freq.span))
	{
		freq.span=c(min(hgram$y), max(hgram$y))
	}
	if(freq.span[2]>max(hgram$hinstfreq))
	{
		freq.span[2]=max(hgram$y)
		warning("Requested frequency window is higher than maximum frequency in the spectrogram.")
	}

        if(fit.line)
        {
             if(hgram$combine.imfs)
             {
                 warning("User requested the IMF.SUM option but the spectrogram data indicates that the IMFs were combined when HHRender was run (combine.imfs = TRUE).  The IMF sum will still be plotted but the spectrogram will display all the IMFs in the signal.")
             }
             fit.line = rowSums(hgram$averaged.imfs[hgram$x >= time.span[1] & hgram$x <= time.span[2], imf.list])
        }
        else
        {
            fit.line = NULL
        }

        img = list()
        img$x = hgram$x[hgram$x >= time.span[1] & hgram$x <= time.span[2]]
        img$y = hgram$y[hgram$y >= freq.span[1] & hgram$y <= freq.span[2]]
        if(hgram$combine.imfs)
        {
            cluster = hgram$cluster[hgram$x >= time.span[1] & hgram$x <= time.span[2], hgram$y >= freq.span[1] & hgram$y <= freq.span[2],imf.list]
        }
        else{
            cluster = apply(hgram$cluster[hgram$x >= time.span[1] & hgram$x <= time.span[2], hgram$y >= freq.span[1] & hgram$y <= freq.span[2],imf.list], c(1, 2), sum)
        }

        #Determine if we are plotting clustering or amplitudes

        if(clustergram)
        {
            img$z = cluster
        }
        else
        {
            if(hgram$combine.imfs)
            {
                img$z = hgram$z[hgram$x >= time.span[1] & hgram$x <= time.span[2], hgram$y >= freq.span[1] & hgram$y <= freq.span[2],imf.list]
            }
            else 
            {
                img$z = apply(hgram$z[hgram$x >= time.span[1] & hgram$x <= time.span[2], hgram$y >= freq.span[1] & hgram$y <= freq.span[2],imf.list], c(1, 2), sum)
            }
        }

        if(!is.null(cluster.span))
        {
            img$z[cluster <= cluster.span[1] |  cluster >= cluster.span[2]] = 0
        }
       

        if(is.null(amp.span))
        {
             if(scaling == "log")
             {
                 amp.span = c(min(img$z[img$z>0]), max(img$z))
             }
             else
             {
                 amp.span = c(min(img$z), max(img$z))
             }
        }

        if(scaling == "log") #Log 10 scale
        {
            img$z = log10(img$z)
            amp.span = log10(amp.span)
        }

        if(scaling == "sqrt") #Take the square root
        {
            img$z = sqrt(img$z)
            amp.span = sqrt(amp.span)
        }

        trace = list()
        trace$sig = hgram$original.signal[hgram$tt >= time.span[1] & hgram$tt <= time.span[2]]
        trace$tt = hgram$tt[hgram$tt >= time.span[1] & hgram$tt <= time.span[2]]

        HHTPackagePlotter(img, trace, amp.span, opts$img.x.lab, opts$img.y.lab, blur = blur, fit.line = fit.line, colormap = colormap, backcol = backcol, pretty = pretty, grid = grid, colorbar = colorbar, opts = opts)
    
        invisible(img)
}

HHSpecPlot <- function(hspec, freq.span = NULL, scaling = "none", imf.list = NULL, show.total = TRUE, show.fourier = FALSE, scale.fourier = FALSE, show.imfs = FALSE, legend = TRUE, ...)
{
    #Plot the Hilbert spectrum, optionally as individual IMFs, optionally with the scaled Fourier spectrum for comparison
    #INPUTS
    #    HSPEC is the Hilbert spectrogram returned by HHSPECTRUM
    #    FREQ.SPAN is the frequencies to plot, NULL means plot everything
    #    SCALING whether to take the base 10 logarithm of amplitude ("log") or square root of amplitude ("sqrt")  or do nothing ("none")
    #    IMF.LIST means only include these IMFS, NULL includes all of them
    #    SHOW.TOTAL means show the sum of the IMF Hilbert spectra
    #    SHOW.IMFS means plot individual IMFs
    #    SHOW.FOURIER determines whether you want a Fourier spectrum for comparison (TRUE) or not (FALSE)
    #    SCALE.FOURIER scales the Fourier spectrum line to the Hilbert spectrum line if TRUE.  Defaults to FALSE.
    #    LEGEND asks whether to plot a legend.  Additional options will place the legend where you want it.
    #ADDITIONAL OPTIONS
    #    XLAB is the X axis label
    #    YLAB is the Y axis label
    #    LEGEND.LOCATION determines where to put the legend.
    #    TOTAL.COL is the color of the ensemble Hilbert spectrum
    #    TOTAL.LWD is the line weight of the ensemble Hilbert spectrogram
    #    LOTAL.LTY is the line type of the ensemble Hilbert spectrogram
    #    IMF.COLS sets the color of each IMF - a vector with length IMF.LIST    
    #    IMF.LWD is the line weight for the IMFs as a vector with length IMF.LIST
    #    IMF.LTY is the line type for the IMFs as a vector with length IMF.LIST
    #    FOURIER.COL is the color of the Fourier spectrum line
    #    FOURIER.LTY is the line type of the Fourier spectrum line
    #    FOURIER.LWD is the line weight of the Fourier spectrum line

    if(!(show.total | show.imfs | show.fourier))
    {
        stop("Nothing to plot!  Set at least one of SHOW.TOTAL, SHOW.IMFS, or SHOW.FOURIER to TRUE.")
    }

    opts = list(...)

    if(!(scaling == "log" | scaling == "sqrt" | scaling == "none"))
    {
        warning(paste("Did not recognize requested scaling: \"", scaling, "\".  Options are \"log\" (base 10 logarithm), \"sqrt\" (square root), or \"none\""))
        scaling = "none"
    }
    
    if(is.null(freq.span))
    {
        freq.span = c(0, max(hspec$frequency))
    }
   
    hspec$amplitude = hspec$amplitude[hspec$frequency >= freq.span[1] & hspec$frequency<= freq.span[2],]
    hspec$frequency = hspec$frequency[hspec$frequency >= freq.span[1] & hspec$frequency<= freq.span[2]]

    if(!"legend.location" %in% names(opts) & legend)
    {
        opts$legend.location = "topright"
    }


    if(!"total.col" %in% names(opts))
    {
        opts$total.col = "red"
    }

    if(!"total.lwd" %in% names(opts))
    {
        opts$total.lwd = 1
    }
    
    if(!"total.lty" %in% names(opts))
    {
        opts$total.lty = 1
    }

    if(!"xlab" %in% names(opts))
    {
        opts$xlab = "frequency"
    }

    if(!"ylab" %in% names(opts))
    {
        if(scaling != "none")
        {
            opts$ylab = paste(scaling, "amplitude")
        }
        else
        {
             opts$ylab = "amplitude"
        }
    }
    
    if(is.null(imf.list))
    {
        imf.list = seq(dim(hspec$amplitude)[2])
    }

    if(!"imf.cols" %in% names(opts))
    {
        if(show.total)
        {
            opts$imf.cols = rainbow(length(imf.list), start = 1/6, end = 5/6)
        }
        else
        {
            opts$imf.cols = rainbow(length(imf.list), start = 0, end = 5/6)
        }
    }
   
    if(!"imf.lwd" %in% names(opts))
    {
        opts$imf.lwd = rep(1, length(imf.list))
    }

    if(!"imf.lty" %in% names(opts))
    {
        opts$imf.lty = rep(1, length(imf.list))
    }

    if(!"fourier.col" %in% names(opts))
    {
        opts$fourier.col = "black"
    }

    if(!"fourier.lty" %in% names(opts))
    {
        opts$fourier.lty = 1
    }
   
    if(!"fourier.lwd" %in% names(opts))
    {
        opts$fourier.lwd = 1
    }

    if(!"main" %in% names(opts))
    {
        opts$main = ""
    } 

    pmin = Inf
    pmax = -Inf

    if(show.imfs)
    {
        imf.amp = hspec$amplitude[, imf.list]
        pmin = min(imf.amp[imf.amp>0])
        pmax = max(imf.amp)
    }

    if(show.total)
    {
        if(length(imf.list)>1)
        {
            total.amp = apply(hspec$amplitude[,imf.list], 1, sum)
        }
        else
        {
            total.amp = hspec$amplitude[,imf.list]
        }
        if(max(total.amp) > pmax)
        {
            pmax = max(total.amp[total.amp > 0])
        }
        if(min(total.amp) < pmin)
        {
            pmin = min(total.amp[total.amp > 0])
        }
    }

     if(show.fourier)
     {
        fourier.freqs = seq(0, 1/(mean(diff(hspec$tt)) * 2), length.out = length(hspec$original.signal)-1)
        fspec = Mod(fft(hspec$original.signal - mean(hspec$original.signal)))[1:length(hspec$original.signal)/2][fourier.freqs >= freq.span[1] & fourier.freqs <= freq.span[2]]
        if(scale.fourier)
        {
            fspec = fspec * pmax/max(fspec)
        }
        if(max(fspec) > pmax)
        {
            pmax = max(fspec)
        }
        if(min(fspec[fspec > 0]) < pmin)
        {
            pmin = min(fspec[fspec > 0])
        }
    }
    
    if(scaling == "log")
    {
        pmax = log10(pmax)
        pmin = log10(pmin)
    }

    if(scaling == "sqrt")
    {
        pmax = sqrt(pmax)
        pmin = sqrt(pmin)
    }
    
    plot(c(min(hspec$frequency), max(hspec$frequency)), c(pmin, pmax), type = "n", xlab = opts$xlab, ylab = opts$ylab, main = opts$main)

    if(show.imfs)
    {
       for(k in seq_len(length(imf.list)))
       {

           amp = imf.amp[,k]
           if(scaling == "log")
           {
              amp = log10(amp)
           }

           if(scaling == "sqrt")
           {
               amp = sqrt(amp)
           } 
           lines(hspec$frequency[amp > -Inf], amp[amp > -Inf], col = opts$imf.cols[k], lwd = opts$imf.lwd[k], lty = opts$imf.lty[k])
       }
    }
   
    if(show.total) 
    {
        if(scaling == "log")
        {
            total.amp = log10(total.amp)
        }

        if(scaling == "sqrt")
        {
            total.amp = sqrt(total.amp)
        }

        lines(hspec$frequency, total.amp, lwd = opts$total.lwd, lty = opts$total.lty, col = opts$total.col)
    }

    if(show.fourier)
    {

        if(scaling == "log")
        {
            fspec = log10(fspec)
        }

        if(scaling == "sqrt")
        {
            fspec = sqrt(fspec)
        }

        lines(fourier.freqs[fourier.freqs >= freq.span[1] & fourier.freqs <= freq.span[2]], fspec, 
            lty = opts$fourier.lty, lwd = opts$fourier.lwd, col = opts$fourier.col)
    }

    if(legend)
    {
        legend.labs = c()
        legend.cols = c()
        legend.lty = c()
        legend.lwd = c()
        if(show.total)
        {
            legend.labs = c(legend.labs, "Total Hilbert")
            legend.cols = c(legend.cols, opts$total.col)
            legend.lty = c(legend.lty, opts$total.lty)
            legend.lwd = c(legend.lwd, opts$total.lwd) 
        }
        if(show.imfs) 
        {
            legend.labs = c(legend.labs, paste(rep("IMF", length(imf.list)), imf.list))
            legend.cols = c(legend.cols, opts$imf.cols)
            legend.lty = c(legend.lty, opts$imf.lty)
            legend.lwd = c(legend.lwd, opts$imf.lwd)
        }

        if(show.fourier)
        {
            legend.labs = c(legend.labs, "Fourier")
            legend.cols = c(legend.cols, opts$fourier.col)
            legend.lty = c(legend.lty, opts$fourier.lty[1])
            legend.lwd = c(legend.lwd, opts$fourier.lwd[1])
        }
        legend(opts$legend.location, legend = legend.labs, lty = legend.lty, lwd = legend.lwd, col = legend.cols)
     }
}



HHTPackagePlotter <- function(img, trace, amp.span, img.x.lab, img.y.lab, blur = NULL, fit.line = NULL, window = NULL, colormap = NULL, backcol = c(0, 0, 0), pretty = FALSE, grid = TRUE, colorbar = TRUE, opts = list())
{
    #Plots images and time series for Hilbert spectra, Fourier spectra, and cluster analysis.
    #This function is internal to the package and users should not be calling it.
    #
    #INPUTS
    #    IMG is the image portion of the figure
    #        IMG$X is the columns
    #        IMG$Y is the rows
    #        IMG$Z is the image data
    #    TRACE is the time series to plot at the top of the figure
    #        TRACE$SIG is the time series
    #        TRACE$TT is the time of each sample
    #    AMP.SPAN are the maximum and minimum values of the image.
    #    IMG.X.LAB is the label of the X axis of the image
    #    IMG.Y.LAB is the label of the Y axis of the image
    #    BLUR is a list of parameters for a Gaussian image smoothing kernel, if desired.  If not null then
    #        BLUR$SIGMA - Standard deviation of Gaussian kernel.  If a 2 element vector, then the kernel has independent coordinate
    #        BLUR$BLEED - Whether to allow blur to bleed out of the domain of the image 
    #    IMF.SUM is a red line on the time series plot showing the sum of the plotted IMFs, if available
    #        IMF.SUM$SIG is the summed IMFS
    #        IMF.SUM$TT is the time of each sample.  We assume all IMFS have equivalent timing.
    #    WINDOW is the length of the Fourier window, if applicable
    #    COLORMAP is the colormap to use for the image
    #    BACKCOL is the background color of the image
    #    PRETTY allows for nice axis labels
    #    GRID draws a grid on the image
    #    COLORBAR puts a colorbar corresponding to the range of values on the image
    #
    #    OPTS    OTHER POSSIBLE OPTIONS
    #        OPTS$TRACE.FORMAT is the format of the trace minima and maxima in sprintf format
    #        OPTS$IMG.X.FORMAT is the format of the X axis labels of the image in sprintf format
    #        OPTS$IMG.Y.FORMAT is the format of the Y axis labels of the image in sprintf format
    #        OPTS$COLORBAR.FORMAT is the format of the colorbar labels in sprintf format   
    #        OPTS$CEX.LAB is the font size of the image axis labels
    #        OPTS$CEX.COLORBAR is the font size of the colorbar
    #        OPTS$CEX.TRACE is the font size of the trace axis labels
    #        OPTS$TRACE.COL is the color of the trace
    #        OPTS$IMF.SUM.COL is the color of the IMF sums (if shown)
    #        OPTS$PRETTY.X.N is the number of pretty divisions on the X axis
    #        OPTS$PRETTY.Y.N is the number of pretty divisions on the Y axis

    #Configure parameters
    
    if(!"trace.format" %in% names(opts))
    {
        opts$trace.format = "%.1e"
    }
 
    if(!"img.x.format" %in% names(opts))
    {
        opts$img.x.format = "%.2f"
    }
   
    if(!"img.y.format" %in% names(opts))
    {
        opts$img.y.format = "%.2f"  
    }
  
    if(!"colorbar.format" %in% names(opts))
    {
         opts$colorbar.format = "%.1e"
    }
 
    if(!"cex.main" %in% names(opts))
    {
        opts$cex.main = 1
    }
    
    if(!"cex.trace" %in% names(opts))
    {
        opts$cex.trace = opts$cex.main * 0.75
    }

    if(!"cex.colorbar" %in% names(opts))
    {
        opts$cex.colorbar = opts$cex.main * 0.75
    }

    if(!"cex.lab" %in% names(opts))
    {
        opts$cex.lab = opts$cex.main
    }

    if(!"fit.line.col" %in% names(opts))
    {
        opts$fit.line.col = "red"
    }
   
    if(!"trace.col" %in% names(opts))
    {
        opts$trace.col = "black"
    }

    if(!"pretty.x.n" %in% names(opts))
    {
        opts$pretty.x.n = 10
    }

    if(!"pretty.y.n" %in% names(opts))
    {
        opts$pretty.y.n = 5
    }


    if(pretty)
    {   #Get nice divisions
        pretty.x = pretty(img$x, n=opts$pretty.x.n)
        pretty.y = pretty(img$y, n=opts$pretty.y.n) 
        #pretty.x = pretty.x[pretty.x <= max(img$x) & pretty.x >= min(img$x)]
        #pretty.y = pretty.y[pretty.y <= max(img$y) & pretty.y >= min(img$y)]
        if(!is.null(window))
        {
             window = window * ((max(img$x) - min(img$x))/(max(pretty.x) - min(pretty.x)))
        }
        img$z = img$z[img$x <= max(pretty.x) & img$x >= min(pretty.x), img$y <= max(pretty.y) & img$y >= min(pretty.y)]
        img$x = img$x[img$x <= max(pretty.x) & img$x >= min(pretty.x)]
        img$y = img$y[img$y <= max(pretty.y) & img$y >= min(pretty.y)]
        img.x.labels=sprintf(opts$img.x.format, pretty.x)
        img.y.labels=sprintf(opts$img.y.format, pretty.y)
        trace$sig = trace$sig[trace$tt >= min(pretty.x) & trace$tt<= max(pretty.x)]
        trace$tt = trace$tt[trace$tt >= min(pretty.x) & trace$tt<= max(pretty.x)]
        cat("Adjusting Time and Frequency limits to nice looking numbers (the \"pretty\" option is currently set to TRUE)\n")
    }    
    else 
    {    
       img.x.labels=sprintf(opts$img.x.format, seq(min(img$x), max(img$x), length.out = 10))
       img.y.labels=sprintf(opts$img.y.format, seq(min(img$y), max(img$y), length.out=5))
    }    

    if(is.null(colormap))
    {
        colormap=rainbow(500,start=0,end=5/6)
    }

    colorbins = length(colormap)
    
    plot(c(-0.15,1),c(-0.15,1),type="n",axes=FALSE,xlab="", ylab="") # Set up main plot window
 
    #Plot TRACE

    sig = trace$sig - mean(trace$sig)
    trace.y=0.75
    trace.x=0
    trace.yspan=0.10
    trace.xspan=0.9
    trace.at=seq(trace.y,trace.y+trace.yspan,length.out=2)
    trace.labels=c(min(trace$sig), max(trace$sig))
    trace.scale=trace.yspan/(max(sig)-min(sig))
    tt.scale=trace.xspan/(max(trace$tt) - min(trace$tt))
    axis(4,pos=trace.x+trace.xspan,at=trace.at, labels=c("",""), cex.axis=opts$cex.trace)
    lines((trace$tt - min(trace$tt)) * tt.scale + trace.x, trace.y + (sig - min(sig)) * trace.scale, col = opts$trace.col)
    if(!is.null(fit.line))
    {
         lines(((trace$tt - min(trace$tt))*tt.scale+trace.x), (trace.y + (fit.line - min(sig)) * trace.scale), col = opts$fit.line.col)
    }
    rect(trace.x, trace.y, trace.x+trace.xspan, trace.y+trace.yspan)

    #Plot IMAGE

    image.y=0
    image.x=0
    image.yspan=0.75
    image.xspan=0.9
    pixel.width = image.xspan/(length(img$x) * 2)
    pixel.height = image.yspan/(length(img$y) * 2)
    image.xvec=seq(image.x + pixel.width, image.x+image.xspan - pixel.width, length.out=length(img$x))
    image.yvec=seq(image.y + pixel.height, image.y+image.yspan - pixel.height, length.out=length(img$y))
    img.x.at=seq(image.x,image.x+image.xspan,length.out=length(img.x.labels))
    img.y.at=seq(image.y,image.y+image.yspan, length.out=length(img.y.labels))
    rect(image.x,image.y,image.x+image.xspan,image.y+image.yspan,col=rgb(red=backcol[1], green=backcol[2], blue=backcol[3], maxColorValue=255))

    #Add blur, if requested   
    if(!is.null(blur)) {
        if(!"sigma" %in% names(blur)) {
            stop("Please provide a standard deviation value when using the \"blur\" option.")
        } else {
            if("bleed" %in% names(blur)) {
                bleed <- blur$bleed
            } else {
                bleed <- TRUE
            }
        }
        tmp.im <- spatstat::as.im(list(x = image.xvec, y = image.yvec, z = as.matrix(img$z)))
        z <- t(spatstat::blur(tmp.im, sigma = blur$sigma, bleed = bleed)$v)
    } else {
        z <- img$z
    }

    z[z<amp.span[1]] = NA
    z[z>amp.span[2]] = amp.span[2]
    z[z == 0]        = NA


    image(image.xvec,image.yvec, z, zlim = amp.span, col=colormap,add=TRUE)
    axis(2, pos=image.x, at=img.y.at,labels=img.y.labels, cex.axis=opts$cex.lab)
    axis(1,pos=image.y, at=img.x.at,labels=img.x.labels, cex.axis=opts$cex.lab)

    #Plot Fourier window, if applicable
    
    if(!is.null(window))
    {
        rwidth = trace.xspan * window 
        rect(trace.x + trace.xspan - rwidth, trace.y + trace.yspan, trace.x + trace.xspan, trace.y + trace.yspan + 0.01, col = "black")
    }


    #Plot GRID
    if(grid)
    {
        line.color=rgb(red=100, green=100, blue=100, maxColorValue=255)
        line.type=3
        for(k in 2:(length(img.x.at)-1))
        {
            lines(c(img.x.at[k], img.x.at[k]), c(image.y, trace.y+trace.yspan), col=line.color, lty=line.type)
        }

        for(k in 2:(length(img.y.at)-1))
        {
            lines(c(image.x, image.x+image.xspan), c(img.y.at[k], img.y.at[k]), col=line.color, lty=line.type)
        }
    }

    #Plot COLORBAR

    if(colorbar)
    {
        color.x=image.x+image.xspan+0.015
        color.xspan=0.025
        color.y=image.y+image.yspan-0.20
        color.yspan=0.10
        color.xvec=c(color.x,color.x+color.xspan)
        color.yvec=seq(color.y, color.y+color.yspan, length.out=colorbins)
        color.at=seq(color.y,color.y+color.yspan,length.out=2)
        colorbar.matrix=array(seq_len(colorbins),dim=c(1, colorbins))
        image(color.xvec, color.yvec, colorbar.matrix, col=colormap, axes=FALSE, add=TRUE)
    }


    #Plot TEXT


    text(trace.x + trace.xspan + 0.03, trace.y, srt=90, sprintf(opts$trace.format,trace.labels[1]), cex=opts$cex.trace)
    text(trace.x + trace.xspan + 0.03, trace.y+trace.yspan, srt=90, sprintf(opts$trace.format, trace.labels[2]), cex=opts$cex.trace)
    text(image.x-0.095, image.y+image.yspan/2, srt=90, img.y.lab, cex=opts$cex.lab)
    text(image.x+image.xspan/2, image.y-0.1, img.x.lab, cex=opts$cex.lab)
    if("main" %in% names(opts))
    {
        text(trace.x+trace.xspan/2, trace.y+trace.yspan+0.05,opts$main, cex=opts$cex.main)
    }
    if(colorbar)
    {
        text(color.x+0.015, color.y-0.0125, sprintf(opts$colorbar.format, amp.span[1]), cex=opts$cex.colorbar)
        text(color.x+0.015, color.y+color.yspan+0.0125, sprintf(opts$colorbar.format, amp.span[2]), cex=opts$cex.colorbar)
    }
 
}

PlotIMFs <-function(sig, time.span = NULL, imf.list = NULL, original.signal = TRUE, residue = TRUE, fit.line=FALSE, lwd=1, cex=1, ...)
{
    #Better IMF plotter
    #This function plots IMFs on the same plot so they can be checked for mode mixing or other problems.
    #It plots shifted traces in a single window
    #INPUTS
    #    SIG is the signal data structure returned by EEMD or EMD analysis
    #    Note that SIG$AVERAGED.IMFS will be plotted instead of SIG$IMF, likewise SIG$AVERAGED.RESIDUE takes precedence
    #    over SIG$RESIDUE, if both exist.
    #        SIG$IMF is a N by M array where N is the length of the signal and M is the number of IMFs
    #        SIG$ORIGINAL.SIGNAL is the original signal before EEMD
    #        SIG$RESIDUE is the residual after EMD
    #        SIG$DT is the sample rate
    #    TIME.SPAN is a 2 element vector giving the time range to plot
    #    IMF.LIST is the IMFs to plot
    #    ORIGINAL.SIGNAL is a boolean asking if you are going to plot the original signal also (defaults to be on top)
    #    RESIDUE is a boolean asking if you are going to plot the residual (defaults to be on bottom)
    #	 FIT.LINE is a boolean asking if you want to plot a line showing the sum of IMFs on top of the original signal (to check how the selected IMFs reconstruct the original signal)
    #	 LWT is the line weight (for plotting figures)
    #    CEX is the size of text (for plotting figures)
    #    ... other parameters to pass to main plotting function
   
    opts <- list(...)

    if(!"xlab" %in% names(opts)) {
        opts$xlab <- "Time (s)"
    }

    if(!"ylab" %in% names(opts)) {
        opts$ylab <- ""
    }
 
    if(is.null(time.span))
    {
        time.span = c(min(sig$tt), max(sig$tt))
    }
   
    if(is.null(imf.list))
    {
        imf.list = 1:sig$nimf
    }
 
    if("averaged.imfs" %in% names(sig))
    {
        sig$imf=sig$averaged.imfs
    }

    if("averaged.residue" %in% names(sig))
    {
        sig$residue=sig$averaged.residue
    }


    time.ind = which(sig$tt >= time.span[1] & sig$tt <= time.span[2])
    tt = sig$tt[time.ind]
    
    plot(c(0, 1), c(0, 1), type="n", axes=FALSE, xlab=opts$xlab, ylab=opts$ylab, cex.lab=cex)
    
    yax.labs=c()
    snum=length(imf.list)+residue+original.signal
    sp=1/snum # Spacing of subplots

    if(original.signal)
    {
         snum=snum+1
         os=sig$original.signal[time.ind]-mean(sig$original.signal[time.ind])
         scale=max(abs(os)) 
    }
    else
    {
        scale=max(abs(sig$imf[time.ind,imf.list]))
    }
    
    if(residue)
    {
        snum=snum+1
        res=sig$residue[time.ind]-mean(sig$residue[time.ind])
	res=res*(sp/(2*scale))
        yax.labs=append(yax.labs,"Residue")
    }

    
    trace.pos=sp/2 #Where the trace starts on the plot
    imfs=sig$imf*(sp/(scale*2))
    ts=(tt-min(tt))*(1/(time.span[2]-time.span[1]))

    if(residue)
    {
        lines(ts, res+trace.pos, lwd=lwd)
        trace.pos=trace.pos+sp
    } 
    
    for(k in rev(imf.list))
    {
       lines(ts, imfs[time.ind,k]+trace.pos, lwd=lwd)
       trace.pos=trace.pos+sp
       yax.labs=append(yax.labs, paste("IMF",k))
    }
    if(original.signal)
    {
        lines(ts, os*(sp/(2*scale))+trace.pos, lwd=lwd)
        yax.labs=append(yax.labs,"Signal")
        if(fit.line)
        {
            if(length(imf.list)>1)
            {
                fline=rowSums(imfs[time.ind,imf.list])
            }
            else
            {
                fline=imfs[time.ind,imf.list]
            }
            if(residue)
            {
                fline=fline+res
            }
            lines(ts, fline+trace.pos, lwd=lwd, col="red")
        }
    }
    xax.labs=pretty(seq(min(tt), max(tt), length.out=11))
    axis(1, pos=0, at=seq(0,1, length.out=length(xax.labs)), labels=xax.labs, cex.axis=cex)
    axis(2, pos=0, at=seq(sp/2, 1, by=sp), labels=yax.labs, cex.axis=cex)
    segments(c(0,0,1, 0), c(0, 1, 1, 0), c(0,1, 1, 1), c(1,1, 0, 0), lwd=lwd) 
}
