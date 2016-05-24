preWaveform <- function(freq, duration, from, xunit, samp.rate){
    if (!is.numeric(duration) || duration <= 0 || length(duration) != 1)
        stop("'duration' must be a positive numeric of length 1")
    if (!is.numeric(from) || from < 0 || length(from) != 1)
        stop("'from' must be a positive numeric of length 1")
    if (!is.numeric(samp.rate) || samp.rate < 0 || length(samp.rate) != 1)
        stop("'samp.rate' must be a positive numeric of length 1")
    if(!is.numeric(freq) || freq <= 0 || length(freq) != 1)
        stop("'freq' must be a positive numeric of length 1")
    if(xunit == "time"){
        duration <- duration * samp.rate
        from <- from * samp.rate
    }
    return(c(duration = round(duration), from = round(from)))
}

postWaveform <- function(channel, samp.rate, bit, stereo, pcm = FALSE, ...){
    if(!is.numeric(bit) || length(bit)!=1 || (!bit %in% c(0,1,8,16,24,32,64)))
        stop("'bit' must be an integer of length 1 in {0,1,8,16,24,32,64}")
    if(bit == 8)
        channel <- channel + 127
    if(stereo && !is.matrix(channel))
        channel <- matrix(channel, ncol = 2, nrow = length(channel))
    Wobj <- Wave(channel, samp.rate = samp.rate, 
                 bit = if(bit %in% 0:1) 32 else bit, pcm = pcm, ...)
    normalize(Wobj, unit = as.character(bit), center = FALSE)
}

silence <- function(duration = samp.rate, from = 0, samp.rate = 44100, bit = 1, 
                    stereo = FALSE, xunit = c("samples", "time"), ...){
    xunit <- match.arg(xunit)
    durFrom <- preWaveform(freq = 1, duration = duration, from = from, 
        xunit = xunit, samp.rate = samp.rate)
    channel <- rep(0, durFrom["duration"])
    postWaveform(channel = channel, samp.rate = samp.rate, 
        bit = bit, stereo = stereo, ...)
}

sine <- function(freq, duration = samp.rate, from = 0, samp.rate = 44100, bit = 1, 
                 stereo = FALSE, xunit = c("samples", "time"), ...){
    xunit <- match.arg(xunit)
    durFrom <- preWaveform(freq = freq, duration = duration, from = from, 
        xunit = xunit, samp.rate = samp.rate)
    channel <- sin(2 * pi * freq * (durFrom["from"]:(sum(durFrom)-1)) / samp.rate)
    postWaveform(channel = channel, samp.rate = samp.rate, 
        bit = bit, stereo = stereo, ...)
}

sawtooth <- function(freq, duration = samp.rate, from = 0, samp.rate = 44100, bit = 1, 
                     stereo = FALSE, xunit = c("samples", "time"), reverse = FALSE, ...){
    xunit <- match.arg(xunit)
    durFrom <- preWaveform(freq = freq, duration = duration, from = from, 
        xunit = xunit, samp.rate = samp.rate)
    channel <- seq(durFrom["from"], 2*freq*sum(durFrom), 
        length = durFrom["duration"]) %% 2 - 1
    if(!is.logical(reverse) || length(reverse) != 1)
        stop("'reverse' must be a logical value of length 1")
    if(reverse) channel <- rev(channel)
    postWaveform(channel = channel, samp.rate = samp.rate, 
        bit = bit, stereo = stereo, ...)
}

square <- function(freq, duration = samp.rate, from = 0, samp.rate = 44100, bit = 1, 
                   stereo = FALSE, xunit = c("samples", "time"), up = 0.5, ...){
    xunit <- match.arg(xunit)
    durFrom <- preWaveform(freq = freq, duration = duration, from = from, 
        xunit = xunit, samp.rate = samp.rate)
    if(!is.numeric(up) || length(up) != 1 || max(abs(up)) > .5)
        stop("'up' must be a numeric in [-0.5, 0.5] of length 1")
    channel <- sign(seq(durFrom["from"], freq*sum(durFrom), 
                        length = durFrom["duration"])
                    %% 1 - 1 + up)
    postWaveform(channel = channel, samp.rate = samp.rate, 
        bit = bit, stereo = stereo, ...)
}

noise <- function(kind = c("white", "pink", "power", "red"), duration = samp.rate,
                  samp.rate = 44100, bit = 1, stereo = FALSE, 
                  xunit = c("samples", "time"), alpha = 1, ...){
    xunit <- match.arg(xunit)
    kind <- match.arg(kind)
    if(kind != "power" && !missing(alpha))
        warning("alpha ignored if noise kind is not 'power'")
    durFrom <- preWaveform(freq = 1, duration = duration, from = 0, 
        xunit = xunit, samp.rate = samp.rate)
    N <- durFrom["duration"] * (stereo + 1)
    channel <- 
        switch(kind,
            white = rnorm(N),
            pink = TK95(N, alpha = 1),
            power = TK95(N, alpha = alpha),
            red = TK95(N, alpha = 1.5)
        )
    channel <- matrix(channel, ncol = (stereo + 1))
    postWaveform(channel = channel, samp.rate = samp.rate, 
        bit = bit, stereo = stereo, ...)
}




### Power law noise generator by Timmer & Koenig (1995)
### contributed by Anita Thieler
### alpha: decay of Power law (pink noise: alpha=1, red noise: alpha=1.5)

TK95 <- function(N, alpha = 1){ 
    f <- seq(from=0, to=pi, length.out=(N/2+1))[-c(1,(N/2+1))] # Fourier frequencies
    f_ <- 1 / f^alpha # Power law
    RW <- sqrt(0.5*f_) * rnorm(N/2-1) # for the real part
    IW <- sqrt(0.5*f_) * rnorm(N/2-1) # for the imaginary part
    fR <- complex(real = c(rnorm(1), RW, rnorm(1), RW[(N/2-1):1]), 
                  imaginary = c(0, IW, 0, -IW[(N/2-1):1]), length.out=N)
     # Those complex numbers that are to be back transformed for Fourier Frequencies 0, 2pi/N, 2*2pi/N, ..., pi, ..., 2pi-1/N 
     # Choose in a way that frequencies are complex-conjugated and symmetric around pi 
     # 0 and pi do not need an imaginary part
    reihe <- fft(fR, inverse=TRUE) # go back into time domain
    return(Re(reihe)) # imaginary part is 0
}

pulse <- function(freq, duration = samp.rate, from = 0, samp.rate = 44100, bit = 1, 
                  stereo = FALSE, xunit = c("samples", "time"), width = 0.1,
                  plateau = 0.2, interval = 0.5, ...){
    xunit <- match.arg(xunit)
    if((width < 0) || (width > 1)) stop("Parameter 'width' must be between 0 and 1.")
    if((interval < 0) || (interval > 1)) stop("Parameter 'interval' must be between 0 and 1.")
    if((plateau  < 0) || (plateau > 1)) stop("Parameter 'interval' must be between 0 and 1.")
    durFrom <- preWaveform(freq = freq, duration = duration, from = from,
                           xunit = xunit, samp.rate = samp.rate)
    x <- freq * (durFrom["from"]:(sum(durFrom)-1)) / samp.rate
    channel <- .C("pulsewav", as.integer(length(x)),
                  as.double(width), as.double(interval), as.double(plateau),
                  as.double(x), y = double(length(x)))$y
    postWaveform(channel = channel, samp.rate = samp.rate,
                 bit = bit, stereo = stereo, ...)
}