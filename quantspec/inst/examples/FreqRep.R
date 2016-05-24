Y       <- rnorm(32)
freq    <- 2*pi*c(0:31)/32 
levels  <- c(0.25,0.5,0.75)
cFT     <- clippedFT(Y, freq, levels)

plot(cFT)

# Get values for all Fourier frequencies and all levels available.
V.all    <- getValues(cFT)

# Get values for every second frequency available
V.coarse <- getValues(cFT, frequencies = 2*pi*c(0:15)/16, levels = levels)

# Trying to get values on a finer grid of frequencies than available will
# yield a warning and then all values with frequencies closest to that finer
# grid.
V.fine   <- getValues(cFT, frequencies = 2*pi*c(0:63)/64, levels = levels)

# Finally, get values for the available Fourier frequencies from [0,pi] and
# only for tau=0.25
V.part   <- getValues(cFT, frequencies = 2*pi*c(0:16)/32, levels = c(0.25))

# Alternatively this can be phrased like this:
V.part.alt <- getValues(cFT, frequencies = freq[freq <= pi], levels = c(0.25))
