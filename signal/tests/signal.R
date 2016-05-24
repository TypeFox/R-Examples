library("signal")
load("savedTests.Rdata")

x1 <- 1000/(10000/2)
x2 <- 1200/(10000/2)
x3 <- seq(0, 2, by=0.01)
x4 <- seq(0, 1, len=100)
x5 <- c(0, 0.3, 0.3, 0.6, 0.6, 1)
x6 <- c(0, 0, 1, 1/2, 0, 0)
x7 <- seq(0, 2/3-0.0001, length=200) 
x9 <- sin(2*pi*(0:10)/5)
x10 <- seq(0,11,length=500)
x12 <- c(0,4,5,6,8,10) 
x13 <- sin(2*pi*seq(0, 10, length=500)/5) 
x14 <- c(0.09026579, 0.00000000, 0.27079736, 0.00000000, 0.27079736, 0.00000000, 0.09026579)
x15 <- c(1.000000e+00, 1.110223e-16, -6.905559e-01, 4.440892e-16, 8.018905e-01, -1.665335e-16, -3.892083e-01)
x16 <- c(0:4, 6:10)
x17 <- c(0:1, 3:10)
x18 <- sin(2*pi*(0:5)/5)
x19 <- seq(0, pi, length=5)
x20 <- sin(2*pi*seq(0, 10, by=0.05)/5)
x21 <- seq(1, 4, by=2)
x22 <- sin(2*pi*x4*2.3)
x23 <- 2*pi*50   
x27 <- seq(0, 2, by=0.001)

# an
Test2 <- 120*an(30) + 125*an(-160)
# Arma
Test3 <- Arma(b = c(1,2,1)/3, a = c(1,1))
# bartlett
Test4 <- bartlett(1)
Test5 <- bartlett(2)
Test6 <- bartlett(51)
# bilinear
Test7 <- bilinear(Sz=signal:::ncauer(3, 40, 5), 2)                
Test8 <- bilinear(ellip(ellipord(x1, x2, 0.5, 29)), 2) 
Test9 <- bilinear(butter(4, 0.1, type="pass", plane="z"), 2)
# blackman
Test10 <- blackman(1)
Test11 <- blackman(2)
Test12 <- blackman(51)
# boxcar
Test13 <- boxcar(2)
Test14 <- boxcar(51)
# butter
Test15 <- butter(4, 0.1, type="pass", plane="z") 
Test16 <- butter(buttord(x1, x2, 0.5, 29))    # buttord.Rd
Test17 <- butter(5, 0.1)  # cheby1.Rd
# buttord
Test18 <- buttord(x1, x2, 0.5, 29)
# cheb
Test19 <- signal:::cheb(49, cosh(1/49 * acosh(1/10^(-5)))*cos((pi*0:49)/50))  
# cheb1ord
Test20 <- cheb1ord(x1, x2, 0.5, 29)
# chebwin
Test21 <- chebwin(50, 100)
# cheby1
Test22 <- cheby1(cheb1ord(x1, x2, 0.5, 29))   # cheb1ord.Rd
Test23 <- cheby1(5, 3, 0.1)
Test24 <- cheby1(5, 0.5, 0.5)
Test25 <- cheby1(3, 3, 2*c(1000,3000)/8000, 'stop')  # grpdelay.R
# cheby2
Test26 <- cheby2(5, 20, 0.5)  # cheby1.Rd
# chirp
Test27 <- chirp(seq(0, 0.6, len=5000))
Test28 <- chirp(seq(0, 5, by=0.001))
Test29 <- chirp(seq(-2, 15, by=0.001), 400, 10, 100, 'quadratic')
Test30 <- chirp(seq(0, 5, by=1/8000), 200, 2, 500, "logarithmic")
Test31 <- chirp(x27, 0, 2, 500)   # specgram.R
Test32 <- chirp(x3, 2, 0.5, 10, 'quadratic') + sin(2*pi*x3*0.4)  # decimate.Rd
# conv
Test33 <- conv(c(1,2,3), c(1,2))
Test34 <- conv(c(1,2), c(1,2,3))
Test35 <- conv(c(1,-2), c(1,2))
# decimate
Test36 <- decimate(chirp(x3, 2, 0.5, 10, 'quadratic') + sin(2*pi*x3*0.4), 4)
# ellip
Test37 <- ellip(5, 3, 40, 0.1)
Test38 <- ellip(ellipord(x1, x2, 0.5, 29))    # ellipord.Rd
# ellipke
Test39 <- signal:::ellipke(c(0.0, 0.01, 0.1, 0.5, 0.9, 0.99, 1.0)) # test
# ellipord
Test40 <- ellipord(c(0.1, 0.2), 0.4, 1, 90)
Test41 <- ellipord(x1, x2, 0.5, 29)
# FilterOfOrder
# fftfilt
Test42 <- set.seed(1)
Test43 <- fftfilt(rep(1, 10)/10, x22 + 0.25*rnorm(length(x4)))   # example with random numbers, set.seed
# filter
Test44 <- set.seed(1)
Test45 <- filter(butter(3, 0.1), x22 + 0.25*rnorm(length(x4)))   # example with random numbers, set.seed
Test46 <- filter(MedianFilter(7), x22 + 0.25*rlnorm(length(x4), 0.5))      # medfilt1.Rd
# filtfilt
Test47 <- set.seed(1)
Test48 <- filtfilt(butter(3, 0.1), x22 + 0.25*rnorm(length(x4))) # example with random numbers, set.seed
# fir1
Test49 <- fir1(40, 0.3)
Test50 <- fir1(10, c(0.3, 0.5), "stop")
Test51 <- fir1(10, c(0.3, 0.5), "pass")
Test52 <- fir1(15, c(0.2, 0.5), "stop")
Test53 <- fir1(15, c(0.2, 0.5), "stop", scale = 'noscale')
Test54 <- fir1(2, 0.5, 'low',  hanning, scale = TRUE)
Test55 <- fir1(2, 0.5, 'low', "hanning", 'scale')
Test56 <- fir1(2, 0.5, 'low', hamming(3), 'scale')
Test57 <- fir1(10, 0.5, scale='noscale')
Test58 <- fir1(10, 0.5, 'low', 'hamming', 'noscale')
Test59 <- fir1(10, 0.5, 'high')
Test60 <- fir1(10, 0.5, 'high', 'hamming', 'scale')
Test61 <- fir1(10, 0.5, window = 'boxcar')
Test62 <- fir1(10, 0.5, 'low', 'boxcar', 'scale')
Test63 <- fir1(10, 0.5, window='hanning', scale='scale')
Test64 <- fir1(10, 0.5, scale='scale', window='hanning', type='low')
Test65 <- fir1(10, 0.5, window='hanning', scale='noscale')
Test66 <- fir1(10, 0.5, scale='noscale', window='hanning', 'low')
Test67 <- fir1(10, 0.5, window=boxcar(11), c())
Test68 <- fir1(40, 2*(0.3*8000/2)/8000)    # grpdelay.R
# fir2
Test69 <- fir2(100, x5, x6)
Test70 <- fir2(50, x5, x6, 512, 0)
Test71 <- fir2(50, x5, x6, 512, 25.6)
Test72 <- fir2(50, x5, x6, 512, 51.2)
Test73 <- fir2(20, c(x7, 2/3, 1), c((1 - (x7/2/3)^2)^(-1/4), 0, 0))
Test74 <- fir2(50, c(x7, 2/3, 1), c((1 - (x7/2/3)^2)^(-1/4), 0, 0))     
Test75 <- fir2(200, c(x7, 2/3, 1), c((1 - (x7/2/3)^2)^(-1/4), 0, 0)) 
# flattopwin
Test76 <- flattopwin(1, sym = 'periodic')
Test77 <- flattopwin(2, sym = 'symmetric')
Test78 <- flattopwin(2, sym = 'periodic')
Test79 <- flattopwin(51, sym = 'symmetric')
Test80 <- flattopwin(51, sym = 'periodic')       
# fractdiff
Test81 <- signal:::fractdiff(c(1,2,3), 1)
Test82 <- signal:::fractdiff(c(1,2,3), 0)
Test83 <- signal:::fractdiff(c(1,2,3), -0.5)  
Test84 <- try(signal:::fractdiff(c(1,2,3), -2))
Test85 <- try(signal:::fractdiff(1,1))
# freqs
Test86 <- unclass(freqs(c(1,2), c(1,1), seq(0, 4, length=128)))               
# freqz
Test87 <- unclass(freqz(c(0.292893218813452, 0.585786437626905, 0.292893218813452), c(1, 0, 0.171572875253810), 32)) # test
Test88 <- unclass(freqz(c(1,1,1)/3, 1, 32, 'whole', plot=FALSE)) # test
Test89 <- unclass(freqz(c(1,1,1)/3, 1, 16, 'half')) # test
Test90 <- unclass(freqz(c(1,1,1)/3, 1, 16, Fs = 320)) # test
Test91 <- unclass(freqz(c(1,1,1)/3, 1, (0:15)*10, Fs = 320)) # test
Test92 <- unclass(freqz(c(1,1,1)/3, 1, 32, 'whole', 320)) # test
Test93 <- unclass(freqz(c(1, 0, -1), c(1, 0, 0, 0, 0.25)))
Test94 <- unclass(freqz(butter(5, 0.1)))  # cheby1.Rd
Test95 <- unclass(freqz(cheby1(5, 3, 0.1)))   # cheby1.Rd
Test96 <- unclass(freqz(cheby1(5, 0.5, 0.5))) # cheby1.Rd
Test97 <- unclass(freqz(cheby2(5, 20, 0.5)))  # cheby1.Rd
Test98 <- unclass(freqz(butter(5, 0.1)))      # ellip.Rd
Test99 <- unclass(freqz(ellip(5, 3, 40, 0.1)))    # ellip.Rd
Test100 <- unclass(freqz(fir1(40, 0.3)))                   # fir1.Rd
Test101 <- unclass(freqz(fir1(10, c(0.3, 0.5), "stop")))   # fir1.Rd
Test102 <- unclass(freqz(fir1(10, c(0.3, 0.5), "pass")))   # fir1.Rd
Test103 <- unclass(freqz(fir2(100, x5, x6)))       # fir2.Rd
Test104 <- unclass(freqz(fir2(50, x5, x6, 512, 0)))    # fir2.R
Test105 <- unclass(freqz(fir2(50, x5, x6, 512, 25.6))) # fir2.R
Test106 <- unclass(freqz(fir2(50, x5, x6, 512, 51.2))) # fir2.R
Test107 <- unclass(freqz(fir2(20, c(x7, 2/3, 1), c((1 - (x7/2/3)^2)^(-1/4), 0, 0))))   # fir2.R
Test108 <- unclass(freqz(fir2(50, c(x7, 2/3, 1), c((1 - (x7/2/3)^2)^(-1/4), 0, 0))))   # fir2.R
Test109 <- unclass(freqz(fir2(200, c(x7, 2/3, 1), c((1 - (x7/2/3)^2)^(-1/4), 0, 0))))  # fir2.R
# gausswin
Test110 <- gausswin(51, 5)
Test111 <- gausswin(2) 
Test112 <- gausswin(2, 5)
# grpdelay
Test113 <- unclass(grpdelay(c(1, 0.9), 1, 512, 'whole', 1))
Test114 <- unclass(grpdelay(poly(c(1/0.9*exp(1i*pi*0.2), 0.9*exp(1i*pi*0.6))), poly(c(0.9*exp(-1i*pi*0.6), 1/0.9*exp(-1i*pi*0.2))), 512, 'whole', 1))
Test115 <- unclass(grpdelay(c(0,1)))
Test116 <- unclass(grpdelay(c(0,1), 1)) 
Test117 <- unclass(grpdelay(c(0,1), 1, 4)) # test
Test118 <- unclass(grpdelay(c(0,1), 1, 4, 'whole')) # test
Test119 <- unclass(grpdelay(c(0,1), 1, 4, Fs = 0.5)) # test
Test120 <- unclass(grpdelay(c(0,1), 1, 4, 'whole', 1)) # test
Test121 <- unclass(grpdelay(c(1, -0.9i), 1, 4, TRUE, 1)) # test
Test122 <- unclass(grpdelay(1, c(1, 0.9), 4)) # test
Test123 <- unclass(grpdelay(c(1,2), c(1, 0.5, 0.9), 4)) # test
Test124 <- unclass(grpdelay(c(1,2), c(1, 0.5, 0.25), 4)) # test
Test125 <- unclass(grpdelay(conv(c(1,2), c(0.25, 0.5, 1)), 1, 4)) # test
# hamming
Test128 <- hamming(1)
Test129 <- hamming(2)
Test130 <- hamming(5)
Test131 <- hamming(51)
# hanning
Test132 <- hanning(1)
Test133 <- hanning(2)
Test134 <- hanning(5)
Test135 <- hanning(51)
# ifft
Test136 <- ifft(fft(1:4))
Test137 <- ifft(fft(signal:::postpad(c(1,2,3), 4)) * fft(signal:::postpad(c(1,2), 4)))      
Test138 <- ifft(fft(signal:::postpad(c(1,-2), 3)) * fft(signal:::postpad(c(1,2), 3)))      
# impz
Test139 <- unclass(impz(butter(5, 0.3)))                                                  
Test140 <- unclass(impz(ellip(5, 0.5, 30, 0.3)))                                          
# interp1
Test141 <- interp1(0:10, x9, x10, 'linear', extrap = TRUE)
Test142 <- interp1(0:10, x9, x10, 'spline', extrap = TRUE)
Test143 <- interp1(0:10, x9, x10, 'pchip', extrap = TRUE)
Test144 <- interp1(0:10, x9, x10, 'cubic', extrap = TRUE)
Test145 <- interp1(0:10, x9, x10, 'nearest', extrap = TRUE)
Test146 <- interp1(x12, sin(2*pi*x12/5), x13, 'linear')
Test147 <- interp1(x12, sin(2*pi*x12/5), x13, 'spline')
Test148 <- interp1(x12, sin(2*pi*x12/5), x13, 'cubic')
Test149 <- interp1(x12, sin(2*pi*x12/5), x13, 'nearest')
Test150 <- interp1(x16, sin(2*pi*x16/5), x13, 'linear')
Test151 <- interp1(x16, sin(2*pi*x16/5), x13, 'spline')
Test152 <- interp1(x16, sin(2*pi*x16/5), x13, 'cubic')
Test153 <- interp1(x16, sin(2*pi*x16/5), x13, 'nearest')
Test154 <- interp1(0:10, x9, x13, 'linear')
Test155 <- interp1(0:10, x9, x13, 'spline')
Test156 <- interp1(0:10, x9, x13, 'cubic')
Test157 <- interp1(0:10, x9, x13, 'nearest')
Test158 <- interp1(x17, sin(2*pi*x17/5), x13, 'linear')
Test159 <- interp1(x17, sin(2*pi*x17/5), x13, 'spline')
Test160 <- interp1(x17, sin(2*pi*x17/5), x13, 'cubic')
Test161 <- interp1(x17, sin(2*pi*x17/5), x13, 'nearest')
Test162 <- interp1(0:10, x9, x20, 'linear')
Test163 <- interp1(0:10, x9, x20, 'spline')
Test164 <- interp1(0:10, x9, x20, 'cubic')
Test165 <- interp1(0:10, x9, x20, 'nearest')
Test166 <- interp1(0:5, x18, c(min(0:5)-1, max(0:5)+1))    # test
Test167 <- interp1(0:5, x18, 0:5, 'nearest') # test
Test168 <- interp1(0:5, x18, c(-1, max(0:5)+1)) # test
Test169 <- interp1(0:5, x18, 0:5, 'linear')   # test
Test170 <- interp1(0:5, x18, 0:5, 'cubic') # test
Test171 <- interp1(0:5, x18, 0:5, 'spline')   # test
Test172 <- interp1(1:5, seq(3, 11, by=2), c(0,6), 'linear', 'extrap') # test
Test173 <- interp1(0:5, x18, c(-1, max(0:5)+1), 'linear', 5)   # test
Test174 <- interp1(1:2, 1:2, 1.4, 'nearest')  # test
Test175 <- interp1(1:4, 1:4, 1.4, 'cubic')    # test
Test176 <- interp1(1:3, 1:3, 1.4, 'spline')   # test
Test177 <- interp1(x21, x21, 1.4, 'nearest')   # test
Test178 <- interp1(seq(1, 8, by=2), seq(1, 8, by=2), 1.4, 'cubic') # test
Test179 <- interp1(seq(1, 6, by=2), seq(1, 6, by=2), 1.4, 'spline')    # test
Test180 <- interp1(x21, x21, c(0, 1, 1.4, 3, 4), 'linear') # test
Test181 <- interp1(1:2, 1:2, 1.4, 'linear') # test
Test182 <- interp1(t(0:5), t(x18), c(), 'nearest') # test isempty
Test183 <- interp1(0:5, x18, c(), 'nearest')    # test isempty
Test184 <- interp1(t(0:5), t(x18), c(), 'linear') # test isempty
Test185 <- interp1(0:5, x18, c(), 'linear')    # test isempty
Test186 <- interp1(t(0:5), t(x18), c(), 'cubic') # test isempty
Test187 <- interp1(0:5, x18, c(), 'cubic')    # test isempty
Test188 <- interp1(t(0:5), t(x18), c(), 'spline') # test isempty
Test189 <- interp1(0:5, x18, c(), 'spline')    # test isempty
# interp
x26 <- chirp(x3, 2, 0.5, 10, 'quadratic') + sin(2*pi*x3*0.4)
Test190 <- interp(x26[seq(1, length(x26), by=4)], 4, 4, 1)
Test191 <- interp(1, 4, 4, 1)
# kaiser
Test192 <- kaiser(2, 5)
Test193 <- kaiser(2, 10)
Test194 <- kaiser(101, 2)
Test195 <- kaiser(101, 10)
Test196 <- kaiser(101, 50)
# kaiserord
Test197 <- kaiserord(c(1200,1500), c(1,0), c(0.1, 0.1), 11025)
Test198 <- kaiserord(c(1000,1500), c(0,1), c(0.1, 0.1), 11025)
Test199 <- kaiserord(c(1000,1200,3000,3500), c(0,1,0), 0.1, 11025)
Test200 <- kaiserord(100 * c(10,13,15,20,30,33,35,40), c(1,0,1,0,1), 0.05, 11025)
# logseq
Test201 <- signal:::logseq(1, 100, n=500)
# Ma
Test202 <- Ma(c(1,2,1)/3)
# MedianFilter
Test203 <- MedianFilter(7)  
# mkpp see pchip
# ncauer
Test204 <- signal:::ncauer(3, 40, 5)
# pchip
Test205 <- pchip(0:10, x9, x10)
Test206 <- pchip(x17, sin(2*pi*x17/5), seq(0, 10, length=500))
Test207 <- pchip(0:10, x9, seq(0, 10, length=500))
m <- diff(cbind(sin(x19), cos(x19))) / (x19[2]-x19[1])
b <- cbind(sin(x19), cos(x19))[1:4,]
pp <- signal:::mkpp(x19, cbind(as.vector(m), as.vector(b)))
Test208 <- signal:::ppval(pp, x19)   
# poly
Test209 <- poly(c(1,-1))
Test210 <- poly(roots(1:3))
Test211 <- poly(matrix(1:9, 3, 3))    
Test212 <- poly(c(1/0.9*exp(1i*pi*0.2), 0.9*exp(1i*pi*0.6)))   # grpdelay.Rd
Test213 <- poly(c(0.9*exp(-1i*pi*0.6), 1/0.9*exp(-1i*pi*0.2))) # grpdelay.Rd
# polyval
Test214 <- polyval(c(1,0,-2), 1:3)
# postpad
Test215 <- signal:::postpad(c(1,2,3), 4)
Test216 <- signal:::postpad(c(1,2), 4)
Test217 <- signal:::postpad(c(1,-2), 3)
Test218 <- signal:::postpad(c(1,2), 3)
# ppval see pchip
# remez
Test219 <- remez(15, c(0, 0.3, 0.4, 1), c(1,1,0,0))        
# resample
Test220 <- resample(sin(2*pi*(0:10)/5), 1, 0.05)
# roots
Test221 <- roots(1:3)
Test222 <- poly(roots(1:3))   
Test223 <- roots(1:3, method="eigen") 
# sftrans
Test224 <- sftrans(signal:::ncauer(3, 40, 5), 0.1, FALSE)       
Test225 <- sftrans(signal:::ncauer(3, 40, 5), 0.1, TRUE)         
Test226 <- sftrans(bilinear(Sz=signal:::ncauer(3, 40, 5), 2)$zero, bilinear(Sz=signal:::ncauer(3, 40, 5), 2)$pole, 
                   bilinear(Sz=signal:::ncauer(3, 40, 5), 2)$gain, 2, stop = FALSE)
Test228 <- sftrans(ellip(ellipord(x1, x2, 0.5, 29)), 0.1)
Test229 <- sftrans(butter(4, 0.1, type="pass", plane="z"), 0.1)

# sgolay
x24 <- t(0:(2^12-1))/(2^12)
x25 <- x24[2]-x24[1]
d <- 1+exp(-3*(x24-0.5))
dd <- -3*exp(-3*(x24-0.5))
d2d <- 9*exp(-3*(x24-0.5))
d3d <- -27*exp(-3*(x24-0.5)) 
x <- d*sin(x23*x24)
dx <- dd*sin(x23*x24) + x23*d*cos(x23*x24)
d2x <- (d2d-x23^2*d)*sin(x23*x24) + 2*x23*dd*cos(x23*x24)
d3x <- (d3d-3*x23^2*dd)*sin(x23*x24) + (3*x23*d2d-x23^3*d)*cos(x23*x24)

Test230 <- sgolayfilt(x, sgolay(8, 41, 0, x25))
Test231 <- sgolayfilt(x, sgolay(8, 41, 1, x25))
Test232 <- sgolayfilt(x, sgolay(8, 41, 2, x25))
Test233 <- sgolayfilt(x, sgolay(8, 41, 3, x25))
# sgolayfilt
Test234 <- sgolayfilt(c(rep(0, 15), rep(10, 10), rep(0, 15)))
Test235 <- sgolayfilt(cos(2*pi*seq(0, 1, by=0.01)*3), 3, 5)  # demo
# sinc
Test236 <- signal:::sinc(c(1,2,3))
# specgram
Test237 <- unclass(specgram(chirp(x27, 0, 2, 500)))  # test
Test238 <- unclass(specgram(chirp(seq(-2, 15, by=0.001), 400, 10, 100, 'quadratic')))
Test239 <- unclass(specgram(chirp(seq(0, 5, by=1/8000), 200, 2, 500, "logarithmic"), Fs = 8000))
Test240 <- unclass(specgram(chirp(x27, 0, 2, 500), 2^ceiling(log2(abs(ceiling(100)))), 1000, ceiling(100), ceiling(100)-ceiling(20)))
# spencer
Test241 <- set.seed(1)
Test242 <- spencer(x22 + 0.25*rnorm(length(x4)))
# triang
Test243 <- triang(1)   # test
Test244 <- triang(2)   # test
Test245 <- triang(3)   # test
Test246 <- triang(4)   # test
Test247 <- triang(51)
# unwrap
Test248 <- unwrap(c(seq(0, 2*pi, length=500), seq(0, 2*pi, length=500)))   
# wav
# Zpg
# zplane
# no values returned


all.equal(Test2, savedTest2)
all.equal(Test3, savedTest3)
all.equal(Test4, savedTest4)
all.equal(Test5, savedTest5)
all.equal(Test6, savedTest6)
all.equal(Test7, savedTest7)
all.equal(Test8, savedTest8)
all.equal(Test9, savedTest9)
all.equal(Test10, savedTest10)
all.equal(Test11, savedTest11)
all.equal(Test12, savedTest12)
all.equal(Test13, savedTest13)
all.equal(Test14, savedTest14)
all.equal(Test15, savedTest15)
all.equal(Test16, savedTest16)
all.equal(Test17, savedTest17)
all.equal(Test18, savedTest18)
all.equal(Test19, savedTest19)
all.equal(Test20, savedTest20)
all.equal(Test21, savedTest21)
all.equal(Test22, savedTest22)
all.equal(Test23, savedTest23)
all.equal(Test24, savedTest24)
all.equal(Test25, savedTest25)
all.equal(Test26, savedTest26)
all.equal(Test27, savedTest27)
all.equal(Test28, savedTest28)
all.equal(Test29, savedTest29)
all.equal(Test30, savedTest30)
all.equal(Test31, savedTest31)
all.equal(Test32, savedTest32)
all.equal(Test33, savedTest33)
all.equal(Test34, savedTest34)
all.equal(Test35, savedTest35)
all.equal(Test36, savedTest36)
all.equal(Test37, savedTest37)
all.equal(Test38, savedTest38)
all.equal(Test39, savedTest39)
all.equal(Test40, savedTest40)
all.equal(Test41, savedTest41)
all.equal(Test42, savedTest42)
all.equal(Test43, savedTest43)
all.equal(Test44, savedTest44)
all.equal(Test45, savedTest45)
all.equal(Test46, savedTest46)
all.equal(Test47, savedTest47)
all.equal(Test48, savedTest48)
all.equal(Test49, savedTest49)
all.equal(Test50, savedTest50)
all.equal(Test51, savedTest51)
all.equal(Test52, savedTest52)
all.equal(Test53, savedTest53)
all.equal(Test54, savedTest54)
all.equal(Test55, savedTest55)
all.equal(Test56, savedTest56)
all.equal(Test57, savedTest57)
all.equal(Test58, savedTest58)
all.equal(Test59, savedTest59)
all.equal(Test60, savedTest60)
all.equal(Test61, savedTest61)
all.equal(Test62, savedTest62)
all.equal(Test63, savedTest63)
all.equal(Test64, savedTest64)
all.equal(Test65, savedTest65)
all.equal(Test66, savedTest66)
all.equal(Test67, savedTest67)
all.equal(Test68, savedTest68)
all.equal(Test69, savedTest69)
all.equal(Test70, savedTest70)
all.equal(Test71, savedTest71)
all.equal(Test72, savedTest72)
all.equal(Test73, savedTest73)
all.equal(Test74, savedTest74)
all.equal(Test75, savedTest75)
all.equal(Test76, savedTest76)
all.equal(Test77, savedTest77)
all.equal(Test78, savedTest78)
all.equal(Test79, savedTest79)
all.equal(Test80, savedTest80)
all.equal(Test81, savedTest81)
all.equal(Test82, savedTest82)
all.equal(Test83, savedTest83)
all.equal(Test84, savedTest84)
all.equal(Test85, savedTest85)
all.equal(Test86, savedTest86)
all.equal(Test87, savedTest87)
all.equal(Test88, savedTest88)
all.equal(Test89, savedTest89)
all.equal(Test90, savedTest90)
all.equal(Test91, savedTest91)
all.equal(Test92, savedTest92)
all.equal(Test93, savedTest93)
all.equal(Test94, savedTest94)
all.equal(Test95, savedTest95)
all.equal(Test96, savedTest96)
all.equal(Test97, savedTest97)
all.equal(Test98, savedTest98)
all.equal(Test99, savedTest99)
all.equal(Test100, savedTest100)
all.equal(Test101, savedTest101)
all.equal(Test102, savedTest102)
all.equal(Test103, savedTest103)
all.equal(Test104, savedTest104)
all.equal(Test105, savedTest105)
all.equal(Test106, savedTest106)
all.equal(Test107, savedTest107)
all.equal(Test108, savedTest108)
all.equal(Test109, savedTest109)
all.equal(Test110, savedTest110)
all.equal(Test111, savedTest111)
all.equal(Test112, savedTest112)
all.equal(Test113, savedTest113)
all.equal(Test114, savedTest114)
all.equal(Test115, savedTest115)
all.equal(Test116, savedTest116)
all.equal(Test117, savedTest117)
all.equal(Test118, savedTest118)
all.equal(Test119, savedTest119)
all.equal(Test120, savedTest120)
all.equal(Test121, savedTest121)
all.equal(Test122, savedTest122)
all.equal(Test123, savedTest123)
all.equal(Test124, savedTest124)
all.equal(Test125, savedTest125)
all.equal(Test128, savedTest128)
all.equal(Test129, savedTest129)
all.equal(Test130, savedTest130)
all.equal(Test131, savedTest131)
all.equal(Test132, savedTest132)
all.equal(Test133, savedTest133)
all.equal(Test134, savedTest134)
all.equal(Test135, savedTest135)
all.equal(Test136, savedTest136)
all.equal(Test137, savedTest137)
all.equal(Test138, savedTest138)
all.equal(Test139, savedTest139)
all.equal(Test140, savedTest140)
all.equal(Test141, savedTest141)
all.equal(Test142, savedTest142)
all.equal(Test143, savedTest143)
all.equal(Test144, savedTest144)
all.equal(Test145, savedTest145)
all.equal(Test146, savedTest146)
all.equal(Test147, savedTest147)
all.equal(Test148, savedTest148)
all.equal(Test149, savedTest149)
all.equal(Test150, savedTest150)
all.equal(Test151, savedTest151)
all.equal(Test152, savedTest152)
all.equal(Test153, savedTest153)
all.equal(Test154, savedTest154)
all.equal(Test155, savedTest155)
all.equal(Test156, savedTest156)
all.equal(Test157, savedTest157)
all.equal(Test158, savedTest158)
all.equal(Test159, savedTest159)
all.equal(Test160, savedTest160)
all.equal(Test161, savedTest161)
all.equal(Test162, savedTest162)
all.equal(Test163, savedTest163)
all.equal(Test164, savedTest164)
all.equal(Test165, savedTest165)
all.equal(Test166, savedTest166)
all.equal(Test167, savedTest167)
all.equal(Test168, savedTest168)
all.equal(Test169, savedTest169)
all.equal(Test170, savedTest170)
all.equal(Test171, savedTest171)
all.equal(Test172, savedTest172)
all.equal(Test173, savedTest173)
all.equal(Test174, savedTest174)
all.equal(Test175, savedTest175)
all.equal(Test176, savedTest176)
all.equal(Test177, savedTest177)
all.equal(Test178, savedTest178)
all.equal(Test179, savedTest179)
all.equal(Test180, savedTest180)
all.equal(Test181, savedTest181)
all.equal(Test182, savedTest182)
all.equal(Test183, savedTest183)
all.equal(Test184, savedTest184)
all.equal(Test185, savedTest185)
all.equal(Test186, savedTest186)
all.equal(Test187, savedTest187)
all.equal(Test188, savedTest188)
all.equal(Test189, savedTest189)
all.equal(Test190, savedTest190)
all.equal(Test191, savedTest191)
all.equal(Test192, savedTest192)
all.equal(Test193, savedTest193)
all.equal(Test194, savedTest194)
all.equal(Test195, savedTest195)
all.equal(Test196, savedTest196)
all.equal(Test197, savedTest197)
all.equal(Test198, savedTest198)
all.equal(Test199, savedTest199)
all.equal(Test200, savedTest200)
all.equal(Test201, savedTest201)
all.equal(Test202, savedTest202)
all.equal(Test203, savedTest203)
all.equal(Test204, savedTest204)
all.equal(Test205, savedTest205)
all.equal(Test206, savedTest206)
all.equal(Test207, savedTest207)
all.equal(Test208, savedTest208)
all.equal(Test209, savedTest209)
all.equal(Test210, savedTest210)
all.equal(Test211, savedTest211)
all.equal(Test212, savedTest212)
all.equal(Test213, savedTest213)
all.equal(Test214, savedTest214)
all.equal(Test215, savedTest215)
all.equal(Test216, savedTest216)
all.equal(Test217, savedTest217)
all.equal(Test218, savedTest218)
all.equal(Test219, savedTest219)
all.equal(Test220, savedTest220)
all.equal(Test221, savedTest221)
all.equal(Test222, savedTest222)
all.equal(Test223, savedTest223)
all.equal(Test224, savedTest224)
all.equal(Test225, savedTest225)
all.equal(Test226, savedTest226)
all.equal(Test228, savedTest228)
all.equal(Test229, savedTest229)
all.equal(Test230, savedTest230, tolerance = .Machine$double.eps ^ 0.33)
all.equal(Test231, savedTest231, tolerance = .Machine$double.eps ^ 0.33)
all.equal(Test232, savedTest232, tolerance = .Machine$double.eps ^ 0.33)
all.equal(Test233, savedTest233, tolerance = .Machine$double.eps ^ 0.33)
all.equal(Test234, savedTest234)
all.equal(Test235, savedTest235)
all.equal(Test236, savedTest236)
all.equal(Test237, savedTest237)
all.equal(Test238, savedTest238)
all.equal(Test239, savedTest239)
all.equal(Test240, savedTest240)
all.equal(Test241, savedTest241)
all.equal(Test242, savedTest242)
all.equal(Test243, savedTest243)
all.equal(Test244, savedTest244)
all.equal(Test245, savedTest245)
all.equal(Test246, savedTest246)
all.equal(Test247, savedTest247)
all.equal(Test248, savedTest248)
