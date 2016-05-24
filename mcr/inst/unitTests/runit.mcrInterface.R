# Test cases for mcrInterface.R
# 
# Author(s): Fabian Model & André Schützenmeister
###############################################################################


cat("\n\nmcrInterface.R method comparison test cases\n\n")

test.mcreg.call <- function() {
	
	## Some test data
	data.x <- 1:10                                  # regular data x with n=10 elements
	data.y1 <- c(1,3,2,4,6,5,7,9,8,10)              # regular data y with n=10 elements
	data.y2 <- data.y1[-1]                          # data y with only 9 elements
	data.y3 <- data.y1
	data.y3[1] <- data.y3[1] * (-1)                 # data y with negative value
	data.y4 <- as.list(data.y1)                     # data y is list of n=10 elements
    data.y5 <- data.y1
    data.y5[5] <- NA                                # data y with an NA element
    data.y6 <- data.y1
    data.y6[c(3,5)] <- 0                            # data with zero-valued elements
    data.y7 <- c(1, NA, 2, NA, NA, NA, NA, 8, 9, 10)
    data.x7 <- c(NA, 2, 3, 4, NA, NA, NA, NA, NA, NA)
    
    data(creatinine)
    set.seed(19061978)                               
    index <- sample(1:nrow(creatinine), size=30)
    crea.data <- creatinine[index,]                  # no NAs for this choice of seed
    crea.data <- as.matrix(crea.data)
    
    coef.mat <- function(val)
    {
        return(matrix(val, ncol=4, dimnames=list(c("Intercept", "Slope"), c("EST", "SE", "LCI", "UCI"))))
    }
    
    ## Calls with various combinations of parameters checking whether a particular parameter combination is functioning 
    
    checkTrue( class(mcreg(data.x, data.y1, method.reg="LinReg", method.ci="analytical")) == "MCResultAnalytical" )
    checkEquals(getCoefficients(mcreg(crea.data, method.reg="LinReg", method.ci="analytical")), 
                coef.mat(c(-0.11419408, 1.07967176, 0.06192228, 0.04633909, -0.24103611, 0.98475043, 0.01264796, 1.17459309)))    
    checkTrue( class(mcreg(data.x, data.y1, method.reg="LinReg", method.ci="jackknife")) == "MCResultJackknife" )
    checkEquals(getCoefficients(mcreg(crea.data, method.reg="LinReg", method.ci="jackknife")),
                coef.mat(c(-0.11419408, 1.07967176, 0.04268357, 0.03558434, -0.20162740, 1.00678054, -0.02676075, 1.15256298)))
    checkTrue( class(mcreg(data.x, data.y1, method.reg="LinReg", method.ci="bootstrap", method.bootstrap="Student" )) == "MCResultResampling" )
    checkEquals(getCoefficients(mcreg(crea.data, method.reg="LinReg", method.ci="bootstrap", method.bootstrap="Student", rng.seed=42)),
            coef.mat(c( -0.114194075, 1.079671759, 0.056547580, 0.051085212, -0.230026543, 0.975028446, 0.001638393, 1.184315071)))
    checkTrue( class(mcreg(data.x, data.y1, method.reg="LinReg", method.ci="bootstrap", method.bootstrap="quantile" )) == "MCResultResampling" )
    checkEquals(getCoefficients(mcreg(crea.data, method.reg="LinReg", method.ci="bootstrap", method.bootstrap="quantile", rng.seed=42)),
            coef.mat(c(-0.11419408, 1.07967176, NA, NA, -0.23507611, 0.98209947, -0.00924089, 1.19567913)))
    checkTrue( class(mcreg(data.x, data.y1, method.reg="LinReg", method.ci="bootstrap", method.bootstrap="tBoot" )) == "MCResultResampling" )
    checkEquals(getCoefficients(mcreg(crea.data, method.reg="LinReg", method.ci="bootstrap", method.bootstrap="tBoot", rng.seed=42)),
            coef.mat(c(-0.114194075, 1.079671759, 0.061922279, 0.046339095, -0.201760275, 0.984794026, 0.001749079, 1.150777929)))
    checkTrue( class(mcreg(data.x, data.y1, method.reg="LinReg", method.ci="bootstrap", method.bootstrap="BCa" )) == "MCResultBCa" )
    checkEquals(getCoefficients(mcreg(crea.data, method.reg="LinReg", method.ci="bootstrap", method.bootstrap="BCa", rng.seed=42)),
            coef.mat(c(-0.11419408, 1.07967176, NA, NA, -0.21305026, 0.95790412, 0.03054473, 1.18398325)))
    checkTrue( class(mcreg(data.x, data.y1, method.reg="LinReg", method.ci="nestedbootstrap", nsamples=100, nnested=10)) == "MCResultResampling" )
    checkEquals(getCoefficients(mcreg(crea.data, method.reg="LinReg", method.ci="nestedbootstrap", method.bootstrap="tBoot", nsamples=100, nnested=10, rng.seed=42)),
            coef.mat(c(-0.11419408, 1.07967176, 0.06192228, 0.04633909, -0.22908877, 0.95105827, 0.06943129, 1.16451986)))
    
    checkTrue( class(mcreg(data.x, data.y1, method.reg="WLinReg", method.ci="analytical")) == "MCResultAnalytical" )
    checkEquals(getCoefficients(mcreg(crea.data, method.reg="WLinReg", method.ci="analytical")),
                coef.mat(c(-0.13798869, 1.10026004, 0.07234441, 0.06857747, -0.28617950, 0.95978546, 0.01020212, 1.24073463)))
    checkTrue( class(mcreg(data.x, data.y1, method.reg="WLinReg", method.ci="jackknife")) == "MCResultJackknife" )
    checkEquals(getCoefficients(mcreg(crea.data, method.reg="WLinReg", method.ci="jackknife")),
                coef.mat(c(-0.13798869, 1.10026004, 0.07381674, 0.06670479, -0.28919543, 0.96362148, 0.01321805, 1.23689860)))
    checkTrue( class(mcreg(data.x, data.y1, method.reg="WLinReg", method.ci="bootstrap", method.bootstrap="Student" )) == "MCResultResampling" )
    checkEquals(getCoefficients(mcreg(crea.data, method.reg="WLinReg", method.ci="bootstrap", method.bootstrap="Student", rng.seed=42)),
            coef.mat(c(-0.13798869, 1.10026004, 0.07378377, 0.06704120, -0.28912789, 0.96293237, 0.01315051, 1.23758772)))
    checkTrue( class(mcreg(data.x, data.y1, method.reg="WLinReg", method.ci="bootstrap", method.bootstrap="quantile" )) == "MCResultResampling" )
    checkEquals(getCoefficients(mcreg(crea.data, method.reg="WLinReg", method.ci="bootstrap", method.bootstrap="quantile", rng.seed=42)),
            coef.mat(c(-0.1379887, 1.1002600, NA, NA, -0.2655639, 0.9637718, 0.0225326, 1.2244854)), tolerance=10e-7)
    checkTrue( class(mcreg(data.x, data.y1, method.reg="WLinReg", method.ci="bootstrap", method.bootstrap="tBoot" )) == "MCResultResampling" )
    checkEquals(getCoefficients(mcreg(crea.data, method.reg="WLinReg", method.ci="bootstrap", method.bootstrap="tBoot", rng.seed=42)),
            coef.mat(c(-0.13798869, 1.10026004, 0.07234441, 0.06857747, -0.27572839, 0.95242262, 0.02770668, 1.22187735)))
    checkTrue( class(mcreg(data.x, data.y1, method.reg="WLinReg", method.ci="bootstrap", method.bootstrap="BCa" )) == "MCResultBCa" )
    checkEquals(getCoefficients(mcreg(crea.data, method.reg="WLinReg", method.ci="bootstrap", method.bootstrap="BCa", rng.seed=42)),
            coef.mat(c(-0.13798869, 1.10026004, NA, NA, -0.25305387, 0.94567707, 0.04962065, 1.21199073)))
    checkTrue( class(mcreg(data.x, data.y1, method.reg="WLinReg", method.ci="nestedbootstrap", nsamples=100, nnested=10)) == "MCResultResampling" )
    checkEquals(getCoefficients(mcreg(crea.data, method.reg="WLinReg", method.ci="nestedbootstrap", method.bootstrap="tBoot", nsamples=100, nnested=10, rng.seed=42)),
            coef.mat(c(-0.13798869, 1.10026004, 0.07234441, 0.06857747, -0.29955409, 0.88248461, 0.13016791, 1.25830595)))
    
    checkTrue( class(mcreg(data.x, data.y1, method.reg="PaBa", method.ci="analytical")) == "MCResultAnalytical" )
    checkEquals(getCoefficients(mcreg(crea.data, method.reg="PaBa", method.ci="analytical")),
                coef.mat(c(-0.18462385, 1.15775171, NA, NA, -0.34714286, 1.03611901, -0.05895187, 1.28571429)))
    checkTrue( class(mcreg(data.x, data.y1, method.reg="PaBa", method.ci="jackknife")) == "MCResultJackknife" )
    checkEquals(getCoefficients(mcreg(crea.data, method.reg="PaBa", method.ci="jackknife")),
                coef.mat(c(-0.18462385, 1.15775171, 0.05744392, 0.08229692, -0.30229238, 0.98917412, -0.06695532, 1.32632930)))
    checkTrue( class(mcreg(data.x, data.y1, method.reg="PaBa", method.ci="bootstrap", method.bootstrap="Student" )) == "MCResultResampling" )
    checkEquals(getCoefficients(mcreg(crea.data, method.reg="PaBa", method.ci="bootstrap", method.bootstrap="Student", rng.seed=42)),
            coef.mat(c(-0.18462385, 1.15775171, 0.06010908, 0.06114394, -0.30775171, 1.03250403, -0.06149599, 1.28299938)))
    checkTrue( class(mcreg(data.x, data.y1, method.reg="PaBa", method.ci="bootstrap", method.bootstrap="quantile" )) == "MCResultResampling" )
    checkEquals(getCoefficients(mcreg(crea.data, method.reg="PaBa", method.ci="bootstrap", method.bootstrap="quantile", rng.seed=42)),
            coef.mat(c(-0.18462385, 1.15775171, NA, NA, -0.33680348, 1.03703704, -0.07961955, 1.28662420)))
    checkTrue( class(mcreg(data.x, data.y1, method.reg="PaBa", method.ci="bootstrap", method.bootstrap="BCa" )) == "MCResultBCa" )
    checkEquals(getCoefficients(mcreg(crea.data, method.reg="PaBa", method.ci="bootstrap", method.bootstrap="BCa", rng.seed=42)),
            coef.mat(c(-0.18462385, 1.15775171, NA, NA, -0.28131318, 1.03703704, -0.05041045, 1.28611704)))
    checkTrue( class(mcreg(data.x, data.y1, method.reg="PaBa", method.ci="nestedbootstrap", nsamples=100, nnested=10)) == "MCResultResampling" )
    checkEquals(getCoefficients(mcreg(crea.data, method.reg="PaBa", method.ci="nestedbootstrap", method.bootstrap="tBoot", nsamples=100, nnested=10, rng.seed=42)),
            coef.mat(c(-0.18462385, 1.15775171, 0.05937457, 0.05985386, -0.26279985, 1.06673390, -0.08560833, 1.26997589)))
    
    checkTrue( class(mcreg(data.x, data.y1, method.reg="Deming", method.ci="analytical")) == "MCResultAnalytical" )
    checkEquals(getCoefficients(mcreg(crea.data, method.reg="Deming", method.ci="analytical")),
                coef.mat(c(-0.15173632, 1.11002937, 0.06366338, 0.04764203, -0.28214483, 1.01243909, -0.02132781, 1.20761965)))
    checkTrue( class(mcreg(data.x, data.y1, method.reg="Deming", method.ci="jackknife")) == "MCResultJackknife" )
    checkEquals(getCoefficients(mcreg(crea.data, method.reg="Deming", method.ci="jackknife")),
                coef.mat(c(-0.15173632, 1.11002937, 0.06473892, 0.05667105, -0.28434799, 0.99394399, -0.01912465, 1.22611475)))
    checkTrue( class(mcreg(data.x, data.y1, method.reg="Deming", method.ci="bootstrap", method.bootstrap="Student" )) == "MCResultResampling")
    checkEquals(getCoefficients(mcreg(crea.data, method.reg="Deming", method.ci="bootstrap", method.bootstrap="Student", rng.seed=42)),
            coef.mat(c(-0.15173632, 1.11002937, 0.06108003, 0.05516052, -0.27685309, 0.99703817, -0.02661955, 1.22302057)))
    checkTrue( class(mcreg(data.x, data.y1, method.reg="Deming", method.ci="bootstrap", method.bootstrap="quantile" )) == "MCResultResampling")
    checkEquals(getCoefficients(mcreg(crea.data, method.reg="Deming", method.ci="bootstrap", method.bootstrap="quantile", rng.seed=42)),
            coef.mat(c(-0.15173632, 1.11002937, NA, NA, -0.31043006, 1.04167991, -0.06488571, 1.26377263)))
    checkTrue( class(mcreg(data.x, data.y1, method.reg="Deming", method.ci="bootstrap", method.bootstrap="tBoot" )) == "MCResultResampling")
    checkEquals(getCoefficients(mcreg(crea.data, method.reg="Deming", method.ci="bootstrap", method.bootstrap="tBoot", rng.seed=42)),
            coef.mat(c(-0.15173632, 1.11002937, 0.06366338, 0.04764203, -0.24164877, 1.00554677, -0.03062687, 1.17645376)))
    checkTrue( class(mcreg(data.x, data.y1, method.reg="Deming", method.ci="bootstrap", method.bootstrap="BCa" )) == "MCResultBCa")
    checkEquals(getCoefficients(mcreg(crea.data, method.reg="Deming", method.ci="bootstrap", method.bootstrap="BCa", rng.seed=42)),
            coef.mat(c(-0.15173632, 1.11002937, NA, NA, -0.26557764, 1.01120835, -0.02684544, 1.22902264)))
    checkTrue( class(mcreg(data.x, data.y1, method.reg="Deming", method.ci="nestedbootstrap", nsamples=100, nnested=10)) == "MCResultResampling" )
    checkEquals(getCoefficients(mcreg(crea.data, method.reg="Deming",method.ci="nestedbootstrap", method.bootstrap="tBoot", nsamples=100, nnested=10, rng.seed=42)),
            coef.mat(c(-0.15173632, 1.11002937, 0.06366338, 0.04764203, -0.24414475, 0.97408516, 0.03280226, 1.18498067)))
    
    checkTrue( class(mcreg(data.x, data.y1, method.reg="WDeming", method.ci="jackknife")) == "MCResultJackknife" )
    checkEquals(getCoefficients(mcreg(crea.data, method.reg="WDeming", method.ci="jackknife")),
                coef.mat(c(-0.20833788, 1.15884356, 0.06247617, 0.06288783, -0.33631451, 1.03002368, -0.08036125, 1.28766344)))
    checkTrue( class(mcreg(data.x, data.y1, method.reg="WDeming", method.ci="bootstrap", method.bootstrap="Student" )) == "MCResultResampling")
    checkEquals(getCoefficients(mcreg(crea.data, method.reg="WDeming", method.ci="bootstrap", method.bootstrap="Student", rng.seed=42)),
            coef.mat(c(-0.20833788, 1.15884356, 0.06278168, 0.06259293, -0.33694033, 1.03062775, -0.07973543, 1.28705936)))
    checkTrue( class(mcreg(data.x, data.y1, method.reg="WDeming", method.ci="bootstrap", method.bootstrap="quantile" )) == "MCResultResampling")
    checkEquals(getCoefficients(mcreg(crea.data, method.reg="WDeming", method.ci="bootstrap", method.bootstrap="quantile", rng.seed=42)),
            coef.mat(c(-0.20833788, 1.15884356, NA, NA, -0.32584428, 1.03839501, -0.07738795, 1.28487957)))
    checkTrue( class(mcreg(data.x, data.y1, method.reg="WDeming", method.ci="bootstrap", method.bootstrap="BCa" )) == "MCResultBCa")
    checkEquals(getCoefficients(mcreg(crea.data, method.reg="WDeming", method.ci="bootstrap", method.bootstrap="BCa", rng.seed=42)),
            coef.mat(c(-0.20833788, 1.15884356, NA, NA, -0.31981990, 1.03325692, -0.07200329, 1.27792318)))
    checkTrue( class(mcreg(data.x, data.y1, method.reg="WDeming", method.ci="nestedbootstrap", nsamples=100, nnested=10)) == "MCResultResampling" )
    checkEquals(getCoefficients(mcreg(crea.data, method.reg="WDeming", method.ci="nestedbootstrap", method.bootstrap="tBoot", nsamples=100, nnested=10, rng.seed=42)),
            coef.mat(c(-0.20833788, 1.15884356, 0.06874857, 0.06809192, -0.36464491, 0.97591236, -0.01304833, 1.32763234)))
	
    ## Calls with identical x and y coordinates 
    
    for(mreg in c("LinReg","WLinReg","Deming","PaBa"))
    {
        checkTrue( class(mcreg(data.x, data.x, method.reg=mreg, method.ci="analytical")) == "MCResultAnalytical" )
        checkEquals(getCoefficients(mcreg(data.x,data.x, method.reg=mreg, method.ci="analytical"))[,-2], 
                    coef.mat(c(0,1,0,0,0,1,0,1))[,-2])       
    }

    for(mreg in c("LinReg","WLinReg","Deming","WDeming","PaBa"))
    {
        checkTrue( class(mcreg(data.x, data.x, method.reg=mreg, method.ci="jackknife")) == "MCResultJackknife" )
        checkEquals(getCoefficients(mcreg(data.x,data.x, method.reg=mreg, method.ci="jackknife"))[,-2], 
                coef.mat(c(0,1,0,0,0,1,0,1))[,-2])
    }
    
    for(mreg in c("LinReg","WLinReg","Deming","WDeming","PaBa"))
    {
        checkTrue( class(mcreg(data.x, data.x, method.reg=mreg, method.ci="bootstrap",method.bootstrap.ci="quantile")) == "MCResultResampling" )
        checkEquals(getCoefficients(mcreg(data.x,data.x, method.reg=mreg, method.ci="bootstrap",method.bootstrap.ci="quantile"))[,-2], 
                coef.mat(c(0,1,0,0,0,1,0,1))[,-2])
    }
    
	## Calls with incorrect method options
    
	checkException(mcreg(data.x, data.y1, method.reg="something"))
	checkException(mcreg(data.x, data.y1, method.ci="something"))
    checkException(mcreg(data.x, data.y1, method.ci="bootstrap", method.bootstrap.ci="something"))
	
	## Calls with incorrect numeric arguments
    
    checkException(mcreg(data.x, data.y1, method.reg="Deming", error.ratio=-1))
    checkException(mcreg(data.x, data.y1, method.reg="Deming", error.ratio="L"))
	checkException(mcreg(data.x, data.y1, method.reg="Deming", error.ratio=as.numeric(NA)))
    checkException(mcreg(data.x, data.y1, method.reg="WDeming", iter.max=0))
    checkException(mcreg(data.x, data.y1, method.reg="WDeming", iter.max=numeric(0)))
    checkException(mcreg(data.x, data.y1, method.reg="WDeming", iter.max=23.23))
    checkException(mcreg(data.x, data.y1, method.reg="WDeming", iter.max=as.numeric(NA)))
    checkException(mcreg(data.x, data.y1, method.reg="WDeming", iter.max="L"))
    checkException(mcreg(data.x, data.y1, method.reg="WDeming", threshold=as.numeric(NA)))
    checkException(mcreg(data.x, data.y1, method.reg="WDeming", threshold=numeric(0)))
    checkException(mcreg(data.x, data.y1, method.reg="WDeming", threshold="L"))
    checkException(mcreg(data.x, data.y1, method.reg="WDeming", threshold=-1))
    checkException(mcreg(data.x, data.y1, method.reg="LinReg", alpha=0))
    checkException(mcreg(data.x, data.y1, method.reg="LinReg", alpha=as.numeric(NA)))
    checkException(mcreg(data.x, data.y1, method.reg="LinReg", alpha=numeric(0)))
    checkException(mcreg(data.x, data.y1, method.reg="LinReg", alpha="L"))
    checkException(mcreg(data.x, data.y1, method.reg="LinReg", method.ci="bootstrap", nsamples=0))
    checkException(mcreg(data.x, data.y1, method.reg="LinReg", method.ci="bootstrap", nsamples=as.numeric(NA)))
    checkException(mcreg(data.x, data.y1, method.reg="LinReg", method.ci="bootstrap", nsamples="L"))
    checkException(mcreg(data.x7, data.y7, method.reg="LinReg"))                                                        # NAs in data
    checkException(mcreg(data.x7, data.y7, method.reg="LinReg", na.rm=TRUE))                                            # only single pair of points
	    
	## Calls with incorrectly specified data
    
    checkException(mcreg(as.numeric(NA), as.numeric(NA)))
    checkException(mcreg(data.x, NA))
	checkException(mcreg(data.x[1:2], data.y[1:2]))
    checkException(mcreg(data.x, data.y2))
    
    ## Calls with incorrect combination of parameters
    
    checkException(mcreg(data.x, data.y3, method.reg="WDeming"))                                                        # negative values
    checkException(mcreg(data.x, data.y3, method.reg="PaBa"))
    checkException(mcreg(data.x, data.y3, method.reg="WLinReg"))
    
    checkException(mcreg(data.x, data.y6, method.reg="WDeming"))                                                        # values equal to zero
    checkException(mcreg(data.x, data.y6, method.reg="WLinReg")) 
    
    checkException(mcreg(data.x, data.y1, sample=paste("Sample", 1:9, sep="_")))                                        # not enough samples names
    
    checkException(mcreg(data.x, data.y1, method.reg="WDeming", method.ci="bootstrap", method.bootstrap.ci="tBoot"))
    checkException(mcreg(data.x, data.y1, method.reg="PaBa", method.ci="bootstrap", method.bootstrap.ci="tBoot"))
    checkException(mcreg(data.x, data.y1, method.reg="WDeming", method.ci="analytical")) 
}