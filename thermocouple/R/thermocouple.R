ThermocoupleEquationTypeB <- function(vT)
{
# vT vector with temperatures
sapply(vT,function(x){
if (x<630.615) w<-1:7 else w<-8:16
n <- 0:(length(w)-1)
round(sum(thermocouple::thermocoupleCoefficientsTypeB[w,1]*x^n), 3)
})
}

ThermocoupleEquationTypeE <- function(vT)
{
# vT vector with temperatures
sapply(vT,function(x){
if (x<0) w<-1:13 else w<-14:24
n <- 0:(length(w)-1)
round(sum(thermocouple::thermocoupleCoefficientsTypeE[w,1]*x^n), 3)
})
}

ThermocoupleEquationTypeJ <- function(vT)
{
# vT vector with temperatures
sapply(vT,function(x){
if (x<760) w<-1:8 else w<-9:14
n <- 0:(length(w)-1)
round(sum(thermocouple::thermocoupleCoefficientsTypeJ[w,1]*x^n), 3)
})
}

ThermocoupleEquationTypeK <- function(vT)
{
# vT vector with temperatures
a0 <-  0.118597600000E+00
a1 <- -0.118343200000E-03
a2 <-  0.126968600000E+03
sapply(vT,function(x){
if (x<0) w<-1:10 else  w<-11:20
n <- 0:(length(w)-1)
if (x<0) return(round(sum(thermocouple::thermocoupleCoefficientsTypeK[w,1]*x^n), 3)) else 
return(round(sum(thermocouple::thermocoupleCoefficientsTypeK[w,1]*x^n  + a0*exp(a1*(x - a2)^2) ), 3))
})
}

ThermocoupleEquationTypeN <- function(vT)
{
# vT vector with temperatures
sapply(vT,function(x){
if (x<0) w<-1:8 else  w<-9:19
n <- 0:(length(w)-1)
round(sum(thermocouple::thermocoupleCoefficientsTypeN[w,1]*x^n), 3)
})
}

ThermocoupleEquationTypeR <- function(vT)
{
# vT vector with temperatures
sapply(vT,function(x){
if (x<1064.180) w<-1:10 else {
if (x<1768.100) w<-11:16 else w<-17:21
}
n <- 0:(length(w)-1)
round(sum(thermocouple::thermocoupleCoefficientsTypeR[w,1]*x^n), 3)
})
}

ThermocoupleEquationTypeS <- function(vT)
{
# vT vector with temperatures
sapply(vT,function(x){
if (x<1064.180) w<-1:9 else {
if (x<1768.100) w<-10:14 else w<-15:19
}
n <- 0:(length(w)-1)
round(sum(thermocouple::thermocoupleCoefficientsTypeS[w,1]*x^n), 3)
})
}

ThermocoupleEquationTypeT <- function(vT)
{
# vT vector with temperatures
sapply(vT,function(x){
if (x<0) w<-1:16 else w<-17:24
n <- 0:(length(w)-1)
round(sum(thermocouple::thermocoupleCoefficientsTypeT[w,1]*x^n), 3)
})
}

ThermocoupleInverseEquationTypeB <- function(vV)
{
# vV vector with voltages
sapply(vV,function(x){
n <- 0:(dim(thermocouple::thermocoupleInverseCoefficientsTypeB)[1]-1)
k <- 1
if (x>=2.431) k <- 2
round(sum(thermocouple::thermocoupleInverseCoefficientsTypeB[,k]*x^n), 3)
})
}

ThermocoupleInverseEquationTypeE <- function(vV)
{
# vV vector with voltages
sapply(vV,function(x){
n <- 0:(dim(thermocouple::thermocoupleInverseCoefficientsTypeE)[1]-1)
k <- 1
if (x>=0) k <- 2
round(sum(thermocouple::thermocoupleInverseCoefficientsTypeE[,k]*x^n), 3)
})
}

ThermocoupleInverseEquationTypeJ <- function(vV)
{
# vV vector with voltages
sapply(vV,function(x){
n <- 0:(dim(thermocouple::thermocoupleInverseCoefficientsTypeJ)[1]-1)
k <- 1
if (x>=42.919) k <- 3 else {
if (x>=0) k <- 2 
}
round(sum(thermocouple::thermocoupleInverseCoefficientsTypeJ[,k]*x^n), 3)
})
}

ThermocoupleInverseEquationTypeK <- function(vV)
{
# vV vector with voltages
sapply(vV,function(x){
if (is.na(x)) return(NA)
n <- 0:(dim(thermocouple::thermocoupleInverseCoefficientsTypeK)[1]-1)
k <- 1
if (x>=20.644) k <- 3 else if (x>=0) k <- 2 
round(sum(thermocouple::thermocoupleInverseCoefficientsTypeK[,k]*x^n), 3)
})
}

ThermocoupleInverseEquationTypeN <- function(vV)
{
# vV vector with voltages
sapply(vV,function(x){
if (is.na(x)) return(NA)
n <- 0:(dim(thermocouple::thermocoupleInverseCoefficientsTypeN)[1]-1)
k <- 1
if (x>=20.613) k <- 3 else {
if (x>=0) k <- 2 
}
round(sum(thermocouple::thermocoupleInverseCoefficientsTypeN[,k]*x^n), 3)
})
}

ThermocoupleInverseEquationTypeR <- function(vV)
{
# vV vector with voltages
sapply(vV,function(x){
n <- 0:(dim(thermocouple::thermocoupleInverseCoefficientsTypeR)[1]-1)
k <- 1
if (x>=19.739) k <- 4 else {
if (x>=11.361) k <- 3 else if (x>=1.923) k <- 2 
}
round(sum(thermocouple::thermocoupleInverseCoefficientsTypeR[,k]*x^n), 3)
})
}

ThermocoupleInverseEquationTypeS <- function(vV)
{
# vV vector with voltages
sapply(vV,function(x){
n <- 0:(dim(thermocouple::thermocoupleInverseCoefficientsTypeS)[1]-1)
k <- 1
if (x>=17.536) k <- 4 else {
if (x>=10.332) k <- 3 else if (x>=1.874) k <- 2 
}
round(sum(thermocouple::thermocoupleInverseCoefficientsTypeS[,k]*x^n), 3)
})
}

ThermocoupleInverseEquationTypeT <- function(vV)
{
# vV vector with voltages
sapply(vV,function(x){
n <- 0:(dim(thermocouple::thermocoupleInverseCoefficientsTypeT)[1]-1)
k <- 1
if (x>=0) k <- 2 
round(sum(thermocouple::thermocoupleInverseCoefficientsTypeT[,k]*x^n), 3)
})
}

ThermistorResistance <- function(Tx, R0, betaTH, T0) R0*exp(betaTH*(1/Tx)-(1/T0))
# Estimate thermistor resistance from temperature
# Tx variable temperature
# Ro resistance at temperature To (25C, expressed in Kelvin)
# Beta parameter of the thermistor (calculated or from the data sheet)
# http://hydraraptor.blogspot.co.uk/2007/10/measuring-temperature-easy-way.html

ThermistorSensitivity <- function(T, beta) -beta/T^2
# Thermistor Sensitivity(relative change in resistance for a change in temperature)
# John G. Webster and Halit Eren, 2014
# Measurement, Instrumentation, and Sensors Handbook, Second Edition
# Spatial, Mechanical, Thermal, and Radiation Measurement
# CRC Press

ThermistorCalibrationEquation <- function(R, R0, thCoeffs)
{
# Thermistor calibration equation
# John G. Webster and Halit Eren, 2014
# Measurement, Instrumentation, and Sensors Handbook, Second Edition
# Spatial, Mechanical, Thermal, and Radiation Measurement
# CRC Press
n <- length(thCoeffs)
1/sum(thCoeffs * (log(R/R0)^(1:n-1)))
}

ThermistorCalibrationEquationHoge1 <- function(Rt, A0, A1, A2)
{
# Chiachung Chen, 2009
# Evaluation of resistance–temperature calibration equations for NTC thermistors
# Measurement 42, Elsevier
1/(A0+A1*log(Rt)+A2*log(Rt)^2)
}

ThermistorCalibrationEquationHoge2 <- function(Rt, A0, A1, A2, A3)
{
# Chiachung Chen, 2009
# Evaluation of resistance–temperature calibration equations for NTC thermistors
# Measurement 42, Elsevier
1/(A0+A1*log(Rt)+A2*log(Rt)^2+A3*log(Rt)^3)
}

ThermistorCalibrationEquationHoge3 <- function(Rt, A0, A1, A2, A3, A4)
{
# Chiachung Chen, 2009
# Evaluation of resistance–temperature calibration equations for NTC thermistors
# Measurement 42, Elsevier
1/(A0+A1*log(Rt)+A2*log(Rt)^2+A3*log(Rt)^3+A4*log(Rt)^4)
}

ThermistorCalibrationEquationHoge4 <- function(Rt, A0, A1, A2, A5)
{
# Chiachung Chen, 2009
# Evaluation of resistance–temperature calibration equations for NTC thermistors
# Measurement 42, Elsevier
1/(A0+A1*log(Rt)+A2*log(Rt)^2+A5/log(Rt))
}

ThermistorCalibrationEquationHoge5 <- function(Rt, C1, C2, C3)
{
# Chiachung Chen, 2009
# Evaluation of resistance–temperature calibration equations for NTC thermistors
# Measurement 42, Elsevier
1/((C1+C2*log(Rt))/(1+C3*log(Rt)))
}

ThermistorTemperature <- function(R, R0, betaTH, T0) 1/((1/betaTH)*log(R/R0)+1/(T0+273.15)) - 273.15
# Estimate thermistor temperature from resistance
# R variable resistance
# Ro resistance at temperature To (25C, expressed in Kelvin)
# Beta parameter of the thermistor (calculated or from the data sheet)
# http://www.mosaic-industries.com/embedded-systems/microcontroller-projects/temperature-measurement/ntc-thermistors/resistance-equation

ThermistorTemperatureFitPolynomial <- function(R, R0, A, B, C, D) 1/(A+B*log(R/R0)+C*(log(R/R0))^3+D*(log(R/R0))^5) - 273.15
# Estimate thermistor temperature from resistance
# R variable resistance
# Ro resistance at temperature To (25C, expressed in Kelvin)
# Beta parameter of the thermistor (calculated or from the data sheet)
# http://www.mosaic-industries.com/embedded-systems/microcontroller-projects/temperature-measurement/ntc-thermistors/resistance-equation

ThermistorConvertADCreadingToTemperatureC <- function(adc, R0, T0, betaTH, R1, R2, vadc = 5.0, vcc = 5.0, ADCbits=10)
{
# Convert ADC reading into a temperature in Celcius by using two resistors
# vadc ADC reference
# vcc supply voltage to potential divider
# http://hydraraptor.blogspot.co.uk/2007/10/measuring-temperature-easy-way.html
T0 = T0 + 273.15 # temperature at stated resistance, e.g. 25C
vs = R1 * vcc / (R1 + R2) # effective bias voltage
rs = R1 * R2 / (R1 + R2)# effective bias impedance
k = R0 * exp(-betaTH / T0)  # constant part of calculation
v = adc * vadc / 2^ADCbits   # convert the 10 bit ADC value to a voltage
r = rs * v / (vs - v)     # resistance of thermistor
(betaTH / log(r / k)) - 273.15 # temperature
}

ThermistorConvertTemperatureCtoADCreading <- function(T, R0, T0, R1, R2, betaTH, vadc = 5.0, vcc = 5.0, ADCbits=10)
{
# Convert a temperature into a ADC value by using two resistors
# vadc ADC reference
# http://hydraraptor.blogspot.co.uk/2007/10/measuring-temperature-easy-way.html
T0 = T0 + 273.15 # temperature at stated resistance, e.g. 25C
vs = R1 * vcc / (R1 + R2) # effective bias voltage
rs = R1 * R2 / (R1 + R2)# effective bias impedance
r = R0 * exp(betaTH * (1 / (T + 273.15) - 1 / T0)) # resistance of the thermistor
v = vs * r / (rs + r)     # the voltage at the potential divider
round(v / vadc * 2^ADCbits)  # the ADC reading
}

ThermistorCalculateBeta <- function(R0, T0, R1, T1) log(R1/R0) / ((1/(T1+273.15))-(1/(T0+273.15)))
# Estimate thermistor beta coefficient from two known resistance/temperature values
# RepRap wiki
# Measuring Thermistor Beta
# http://reprap.org/wiki/MeasuringThermistorBeta

ThermistorAlphaApproximatedFromBeta <- function(T, betaTH) -betaTH/T^2 *100
# http://www.daycounter.com/Calculators/Steinhart-Hart-Thermistor-Calculator.phtml

ThermistorResistanceDeviation <- function(deltaBetaTH, deltaR25) (((1+deltaR25)/100)*(1+deltaBetaTH/100) - 1)*100 
ThermistorTemperatureDeviation <- function(deltaBetaTH, deltaR25, alpha) ThermistorResistanceDeviation(deltaBetaTH, deltaR25) / alpha
# http://www.daycounter.com/Calculators/Steinhart-Hart-Thermistor-Calculator.phtml

ThermistorTemperatureSteinhartHart <- function(R, R0, A, B, C=0, D) 1/ (A + B*log(R/R0) + C*log(R/R0)^2 + D*log(R/R0)^3)
# Steinhart-Hart equation for thermistor temperature
# Steinhart-Hart Coefficient A (K^0)
# Steinhart-Hart Coefficient B (K^-1)
# Steinhart-Hart Coefficient C (K^-2)
# Steinhart-Hart Coefficient D (K^-3)
# Ro resistance at temperature To (25°C, expressed in ohms)
# R resistance at temperature T
# Daycounter, Inc. Engineering Services
# Steinhart-Hart Thermistor Calculator
# http://www.daycounter.com/Calculators/Steinhart-Hart-Thermistor-Calculator.phtml

ThermistorResistanceSteinhartHartUsing3T <- function(T, T2, T3, R0, A1, B1, C1=0, D1) R0*exp(A1+B1/T+C1/T2+D1/T3)
# Steinhart-Hart equation for thermistor resistance
# Steinhart-Hart Coefficient A1 (K^0)
# Steinhart-Hart Coefficient B1 (K^1)
# Steinhart-Hart Coefficient C1 (K^2)
# Steinhart-Hart Coefficient D1 (K^3)
# Ro resistance at temperature To (25°C, expressed in ohms)
# T measured temperature for resistance R
# Daycounter, Inc. Engineering Services
# Steinhart-Hart Thermistor Calculator
# http://www.daycounter.com/Calculators/Steinhart-Hart-Thermistor-Calculator.phtml

ThermistorResistanceSteinhartHart <- function(T, A, B, C) ((-27/2*((A-1/T)/C)+3/2*sqrt(3)*(sqrt(27*((A-1/T)/C)^2+4*(B/C)^3)))^(1/3) - 
(27/2*((A-1/T)/C)+3/2*sqrt(3)*(sqrt(27*((A-1/T)/C)^2+4*(B/C)^3)))^(1/3))/3
# Steinhart-Hart equation for thermistor resistance
# AVX Corporation, 2014
# AVX NTC Thermistors v11.4
# http://www.avx.com

ThermistorResistanceSteinhartHart2 <- function(T, A, B, C) exp((sqrt((4*B^3*T^2)/C+27*A^2*T^2-54*A*T+27)/(2*3^(3/2)*C*T)+1/(2*C*T)-A/(2*C))^(1/3)-
B/(3*C*(sqrt((4*B^3*T^2)/C+27*A^2*T^2-54*A*T+27)/(2*3^(3/2)*C*T)+1/(2*C*T)-A/(2*C))^(1/3)))
# Steinhart-Hart equation for thermistor resistance, calculated with Maxima

ThermistorSteinhartHartCoeffFromMeasurements <- function(resAndTemp)
{
# Calculate Steinhart-Hart coefficients A, B, C from measurements
# resAndTemp matrix with temperatures (C) in column 1 and resistance (ohm) in column 2
# NTC Thermistor theory
# BetaTHERM sensors
# www.betatherm.com
if (dim(resAndTemp)[1] != 3) stop('resAndTemp must be 3 x 2')
if (dim(resAndTemp)[2] != 2) stop('resAndTemp must be 3 x 2')
b <- 1/(resAndTemp[,2] + 273.15 )
b <- cbind(b)
a <- matrix(c(1, log(resAndTemp[1,1]), log(resAndTemp[1,1])^3, 1, log(resAndTemp[2,1]), log(resAndTemp[2,1])^3, 1, log(resAndTemp[3,1]), log(resAndTemp[3,1])^3),3,3,byrow=T)
x <- solve(a, b)
list(A=x[1], B=x[2], C=x[3])
}

ThermistorHoge1CoeffFromMeasurements <- function(resAndTemp)
{
# Calculate Hoge1 coefficients A0, A1, A2 from measurements
# resAndTemp matrix with temperatures (C) in column 1 and resistance (ohm) in column 2
if (dim(resAndTemp)[1] != 3) stop('resAndTemp must be 3 x 2')
if (dim(resAndTemp)[2] != 2) stop('resAndTemp must be 3 x 2')
b <- 1/(resAndTemp[,2] + 273.15 )
b <- cbind(b)
a <- matrix(c(1, log(resAndTemp[1,1]), log(resAndTemp[1,1])^2, 1, log(resAndTemp[2,1]), log(resAndTemp[2,1])^2, 1, log(resAndTemp[3,1]), log(resAndTemp[3,1])^2),3,3,byrow=T)
x <- solve(a, b)
list(A0=x[1], A1=x[2], A2=x[3])
}

ThermistorHoge2CoeffFromMeasurements <- function(resAndTemp)
{
# Calculate Hoge2 coefficients A0, A1, A2, A3 from measurements
# resAndTemp matrix with temperatures (C) in column 1 and resistance (ohm) in column 2
if (dim(resAndTemp)[1] != 4) stop('resAndTemp must be 4 x 2')
if (dim(resAndTemp)[2] != 2) stop('resAndTemp must be 4 x 2')
b <- 1/(resAndTemp[,2] + 273.15 )
b <- cbind(b)
a <- matrix(c(1, log(resAndTemp[1,1]), log(resAndTemp[1,1])^2, log(resAndTemp[1,1])^3,
1, log(resAndTemp[2,1]), log(resAndTemp[2,1])^2, log(resAndTemp[2,1])^3,
1, log(resAndTemp[3,1]), log(resAndTemp[3,1])^2, log(resAndTemp[3,1])^3,
1, log(resAndTemp[4,1]), log(resAndTemp[4,1])^2, log(resAndTemp[4,1])^3),3,3,byrow=T)
x <- solve(a, b)
list(A0=x[1], A1=x[2], A2=x[3], A3=x[4])
}

ThermistorHoge3CoeffFromMeasurements <- function(resAndTemp)
{
# Calculate Hoge3 coefficients A0, A1, A2, A3, A4 from measurements
# resAndTemp matrix with temperatures (C) in column 1 and resistance (ohm) in column 2
if (dim(resAndTemp)[1] != 5) stop('resAndTemp must be 5 x 2')
if (dim(resAndTemp)[2] != 2) stop('resAndTemp must be 5 x 2')
b <- 1/(resAndTemp[,2] + 273.15 )
b <- cbind(b)
a <- matrix(c(1, log(resAndTemp[1,1]), log(resAndTemp[1,1])^2, log(resAndTemp[1,1])^3, log(resAndTemp[1,1])^4,
1, log(resAndTemp[2,1]), log(resAndTemp[2,1])^2, log(resAndTemp[2,1])^3, log(resAndTemp[2,1])^4,
1, log(resAndTemp[3,1]), log(resAndTemp[3,1])^2, log(resAndTemp[3,1])^3, log(resAndTemp[3,1])^4,
1, log(resAndTemp[4,1]), log(resAndTemp[4,1])^2, log(resAndTemp[4,1])^3, log(resAndTemp[4,1])^4,
1, log(resAndTemp[5,1]), log(resAndTemp[5,1])^2, log(resAndTemp[5,1])^3, log(resAndTemp[5,1])^4),3,3,byrow=T)
x <- solve(a, b)
list(A0=x[1], A1=x[2], A2=x[3], A3=x[4], A4=x[5])
}

ThermistorHoge4CoeffFromMeasurements <- function(resAndTemp)
{
# Calculate Hoge4 coefficients A0, A1, A2, A5 from measurements
# resAndTemp matrix with temperatures (C) in column 1 and resistance (ohm) in column 2
if (dim(resAndTemp)[1] != 4) stop('resAndTemp must be 4 x 2')
if (dim(resAndTemp)[2] != 2) stop('resAndTemp must be 4 x 2')
b <- 1/(resAndTemp[,2] + 273.15 )
b <- cbind(b)
a <- matrix(c(1, log(resAndTemp[1,1]), log(resAndTemp[1,1])^2, log(resAndTemp[1,1]),
1, log(resAndTemp[2,1]), log(resAndTemp[2,1])^2, log(resAndTemp[2,1]),
1, log(resAndTemp[3,1]), log(resAndTemp[3,1])^2, log(resAndTemp[3,1]),
1, log(resAndTemp[4,1]), log(resAndTemp[4,1])^2, log(resAndTemp[4,1])),3,3,byrow=T)
x <- solve(a, b)
list(A0=x[1], A1=x[2], A2=x[3], A5=x[4])
}

ThermistorAnyNumberOfCoeffFromMeasurements <- function(resAndTemp)
{
# Calculate coefficients A0, A1, A2, ..., An from n+1 measurements
# resAndTemp matrix with temperatures (C) in column 1 and resistance (ohm) in column 2
if (dim(resAndTemp)[1]<3) stop('resAndTemp must be have 3 or more rows')
if (dim(resAndTemp)[2] != 2) stop('resAndTemp must be n x 2')
b <- 1/(resAndTemp[,2] + 273.15 )
b <- cbind(b)
a <-matrix(1,n,n)
for (n in 1:dim(resAndTemp)[1])
for (m in 2:dim(resAndTemp)[1]){
a[n, m] <- log(resAndTemp[n,1])^(m-1)
}
x <- solve(a, b)
list(A0=x[1], A1=x[2], A2=x[3], A3=x[4], A4=x[5])
}

DS1820CalcCRCbit <- function(shiftReg, dataBit)
{
#  Calculate 8-bit CRC for DS1820
# Peter H. Anderson, 98
# DS1820 Digital Thermometer - Calculating an 8-bit CRC Value
# http://www.phanderson.com/PIC/16C84/crc.html
fb <- bitwXor(bitwAnd(shiftReg, 0x01), dataBit) #exclusive or least sig bit of current shift reg with the data bit
shiftReg = shiftReg %/% 2 # shift one place to the right
        if (fb==1)
        {
           shiftReg <- bitwXor(shiftReg, 0x10001100 ) # CRC ^ binary 1000 1100 
           }
shiftReg 
}

ThermistorTemperatureAccuracy <- function(ResTol, alpha) ResTol / alpha
ThermistorResistanceTolerance <- function(TempAccy, alpha) TempAccy * alpha
# Thermistor relationship resistance tolerance vs temperature accuracy
# Spectrum Sensors & Controls Inc., 2014
# NTC Thermistors Engineering Notes
# www.SpecSensors.com

ThermistorVolumeResistivityFromRho <- function(Rho, Thck, L, W) Rho * Thck / (L * W)
# Rho material resistivity in ohm/cm
# Thck thickness of the conductor (chip) (cm)
# L length of the conductor (chip) (cm)
# W width of the conductor (chip) (cm)
# Equation #1
# NTC Thermistor theory
# BetaTHERM sensors
# www.betatherm.com

ThermistorVolumeResistivityFromR25 <- function(R25, Thck, L, W) L * W / Thck * R25
# R25 measured resistance 25C (ohms)
# Thck thickness of the conductor (chip) (cm)
# L length of the conductor (chip) (cm)
# W width of the conductor (chip) (cm)
# Equation #1
# NTC Thermistor theory
# BetaTHERM sensors
# www.betatherm.com

ThermistorSlope <- function(R0, R70) R0 / R70
# Thermistor Slope (Resistance Ratio)
# R0 resistance at 0C
# R70 resistance at 70C
# NTC Thermistor theory
# BetaTHERM sensors
# www.betatherm.com

AWGTOmm <- function(n) 0.127 * 92 ^((36-n)/39)
# convert American wire gauge (SWG) to mm
# n gauge number
# http://www.rapidtables.com/calc/wire/awg-to-mm.htm

ThermocoupleFundamentalRelation<-function(S, T0, T1) S * (T1 - T0)
# T0, T1 temperatures at both ends
# S Seebeck coefficient (uV/C) or Sab Seebeck coefficient between material a and b
# V voltage difference
# pp. 13 eq. 2.1

ThermocoupleFundamentalRelation2<-function(Sa, Sb, T0, T1) (Sa - Sb) * (T1 - T0)
# T0, T1 temperatures at both ends
# Sa Seebeck coefficient of material a
# Sb Seebeck coefficient of material b
# V voltage difference
# pp. 13 eq. 2.4

ThermocoupleVoltageContributionTwoHomogeneousWires<-function(Sab, T0, T1, T2) Sab * (T2 - T0) + Sab * (T1 - T2)
# Voltage Contribution of Two Homogeneous Wires
# T0, T1 temperatures at both ends
# T2 temperature at a point !=T0, T1
# Sab Seebeck coefficient between material a and b
# V voltage difference
# pp. 15 eq. 2.9

ThermocoupleWithReference<-function(Sa, Sb, T0, T1, T2) Sa * (T2 - T0) + Sb * (T1 - T2) + Sa * (T0 - T1)
# Thermocouple with Reference
# T0, T2 temperatures at both ends
# T1 temperature at a reference point
# Sa Seebeck coefficient of material a
# Sb Seebeck coefficient of material b
# V voltage difference
# pp. 17 eq. 2.13

ThermocoupleWithReference2<-function(Sab, T1, T2) Sab * (T2 - T1)
# Thermocouple with Reference
# T0, T2 temperatures at both ends
# T1 temperature at a reference point
# Sa Seebeck coefficient of material a
# Sb Seebeck coefficient of material b
# V voltage difference
# pp. 17 eq. 2.15

ThermocoupleStemLossErrorEstimate <- function(L, h, k, r0, ri){
# E = error (percent of difference between tip temperature and back-end temperature)
# L = sensor insertion depth (cm)
# h = surface heat transfer coefficient (watts.cm2 C)
# k = thermal conductivity of sheath material (watts.cm C)
# r0 = sheath outer radius
# ri = sheath inner radius
# pp. 47 eq. 3.8
alpha <- sqrt(2*r0 * h/(k*(r0^2 - ri^2)))
F <- k * alpha/h
2 * F / ((1 + F)*exp(alpha * L)-(1 - F)*exp(-alpha * L))
}

RTDmetalResistance <- function(R0, T, A, B, C, metal=NA)
{
# Metal RTD resistance
# R0 resistance at 0C
# T temperature in C
# A, B, C specific constants
# http://www.capgo.com/Resources/Temperature/RTDs/RTD.html
if (T==0) return(R0)
if (!is.na(metal)){
A <- B <- C <- 0
if (tolower(metal)=='copper') A <- 0.00427
if (tolower(metal)=='molybdenum') {
A <- 0.00427
B <- 0.00385
}
if (tolower(metal)=='nickel') A <- 0.00672
if (tolower(metal)=='nickel-iron') A <- 0.00518
if (tolower(metal)=='Platinum') {
A <- 0.00385
B <- 0.00392
C <- 0.00377
}
if (A==0) stop('Wrong metal!')
}
R0*( 1 + A*T + B*T^2 + C*T^3 )
}

RTDalpha <- function(R0, R100) (R100-R0)/(100*R0)
# RTD alpha coefficient
# R0 resistance at 0C
# R100 resistance at 100C
# http://www.capgo.com/Resources/Temperature/RTDs/RTD.html

RTDdelta <- function(R0, Rth, Th,alpha) (Th-(Rth-R0)/(R0*alpha))/((Th/100-1)*(Th/100))
# RTD delta coefficient
# R0 resistance at 0C
# Th highest temperature in the calibration range
# Rth resistance of the sensor at the highest temperature
# John G. Webster and Halit Eren, 2014
# Measurement, Instrumentation, and Sensors Handbook, Second Edition
# Spatial, Mechanical, Thermal, and Radiation Measurement
# CRC Press

RTDbeta <- function(R0, Rtl, Tl,alpha, delta)  (Tl-((Rtl-R0)/(R0*alpha)+delta*(Tl/100-1)*(Tl/100)))/((Tl/100-1)*(Tl/100)^3)
# RTD beta coefficient
# R0 resistance at 0C
# Tl lowest temperature in the calibration range
# Rtl resistance of the sensor at the lowest temperature
# John G. Webster and Halit Eren, 2014
# Measurement, Instrumentation, and Sensors Handbook, Second Edition
# Spatial, Mechanical, Thermal, and Radiation Measurement
# CRC Press

RTDequation <- function(R0, T, A, B, C=NA)
{
# RTD equation
# R0 resistance at 0C
# T temperature in C
# A, B, C RTD constants
# John G. Webster and Halit Eren, 2014
# Measurement, Instrumentation, and Sensors Handbook, Second Edition
# Spatial, Mechanical, Thermal, and Radiation Measurement
# CRC Press
if(T==0) return(R0)
if ((T>0)&(T<850)) return(R0*(1+A*T+B*T^2))
if ((T<0)&(T> -200)) return(R0*(1+A*T+B*T^2+C*(T-100)^3))
NA
}

RTDcoefficientA <- function(alpha, delta) alpha+alpha*delta/100
RTDcoefficientB <- function(alpha, delta) alpha*delta/100^2
RTDcoefficientC <- function(alpha, beta) alpha*beta/100^4
# John G. Webster and Halit Eren, 2014
# Measurement, Instrumentation, and Sensors Handbook, Second Edition
# Spatial, Mechanical, Thermal, and Radiation Measurement
# CRC Press

RTDplatinumResistance <- function(R0, T, A=NA, B=NA, C=NA, stdRTD='DIN43760')
{
# Callendar-Van Dusen equation for platinum RTD
# R0 resistance at 0C
# T temperature in C
# A, B, C specific constants (optional)
# Above 0 Celsius the constant C is equal to zero
# http://www.capgo.com/Resources/Temperature/RTDs/RTD.html
# 
# National Instruments, 2014
# Taking Temperature Measurements with RTDs: How-To Guide
if (T==0) return(R0)
if(!is.numeric(A) | !is.numeric(B) | !is.numeric(C))
{
if (stdRTD=='IEC751'){
if ((T<0)&(T>-200)){
A <- 3.90830E-3
B <- -5.77500E-7
C <- -4.18301E-12
}
if ((T>0)&(T<850)){
A <- 3.90830E-3
B <- -5.77500E-7
C <- 0
}
}
if (stdRTD=='SAMA'){
A <- 3.97869E-3
B <- -5.86863E-7
C <- -4.16696E-12
}
if (stdRTD=='DIN43760'){
A <- 3.9080E-3
B <- -5.8019E-7
C <- -4.2735E-12
}
if (stdRTD=='American'){
A <- 3.9692E-3
B <- -5.8495E-7
C <- -4.2325E-12
}
if (stdRTD=='ITS-90'){
A <- 3.9848E-3
B <- -5.870E-7
C <- -4.0000E-12
}
}
if (T<0) Rt <- R0 * ( 1 + A*T + B*T^2 + C*(T - 100)*T^3 )
if (T>0) Rt <- R0 *( 1 + A*T + B*T^2 )
Rt
}

RTDplatinumTemperature <- function(R0, R, alpha, beta, delta)
{
# Callendar-Van Dusen equation for platinum RTD temperature from resistance
# R0 resistance at 0C
# R Measured resistance
# alpha, beta, delta specific constants
# John G. Webster, 1999
# The Measurement, Instrumentation and Sensors Handbook
# CRC Press LLC
# eq (32.24), (32.25)
if(R==R0) return(0)
if(R>R0) return(((R-R0)/(alpha*R0))+delta*((T/100-1)*(T/100)))
((R-R0)/(alpha*R0))+delta*((T/100-1)*(T/100))+beta*((T/100-1)*(T/100)^3)
}

RTDplatinumResistanceFromAlpha <- function(R0, T, alpha=NA, stdRTD='DIN43760')
{
# Callendar-Van Dusen equation
# Omega Inc.,2014
# The RTD
# http://www.omega.com/
if (T==0) return(R0)
if (!is.numeric(alpha))
{
if (stdRTD=='IEC751') alpha <- 0.00385055
if (stdRTD=='SAMA') alpha <- 0.0039200
if (stdRTD=='DIN43760') alpha <- 0.003850
if (stdRTD=='American') alpha <- 0.003911
if (stdRTD=='ITS-90') alpha <- 0.003926
}
d <- 1.49
b <- ifelse(T>0,0,0.11)
R0+R0*alpha*(T - d*(T/100-1)*(T/100)-b*(T/100-1)(T/100)^3)
}

RTDnickelResistance <- function(R0, T, A=NA, B=NA, D=NA, F=NA)
{
# Callendar-Van Dusen equation for Nickel RTD
# R0 resistance at 0C
# T temperature in C
# A, B, C specific constants (optional)
# http://www.capgo.com/Resources/Temperature/RTDs/RTD.html
if (T==0) return(R0)
if(!is.numeric(A) | !is.numeric(B) | !is.numeric(D) | !is.numeric(F))
{
A <- 5.485E-3
B <- 6.650E-6
D <- 2.805E-11
F <- -2.000E-17
}
R0 *(1 + A*T + B*T^2 + D*T^4 + F*T^6 )
}

RTDmetalResistanceFromAlpha <- function(R0, T, alpha=NA, metal='nickel')
{
# 
# alpha (optional) resistance's temperature coefficient resistance's temperature coefficient
# R0 Resistance at 0C
# T temperature C
# Resistive temperature detectors PTxx
# www.madur.com
if (T==0) return(R0)
if (!is.numeric(alpha)){
if (tolower(metal)=='nickel') alpha <- 0.5866
if (tolower(metal)=='iron') alpha <- 0.5671
if (tolower(metal)=='molybdenum') alpha <- 0.4579
if (tolower(metal)=='tungsten') alpha <- 0.4403
if (tolower(metal)=='aluminium') alpha <- 0.4308
if (tolower(metal)=='copper') alpha <- 0.4041
if (tolower(metal)=='silver') alpha <- 0.3819
if (tolower(metal)=='platinum') alpha <- 0.3729
if (tolower(metal)=='gold') alpha <- 0.3715
}
R0 *(1 + alpha*T )
}

RTDnickelResistanceFromAlpha <- function(R0, T, alpha=NA)
{
# simplified equation for Nickel RTD resistance
# R0 resistance at 0C
# T temperature in C
# alpha (optional) resistance's temperature coefficient
# http://www.capgo.com/Resources/Temperature/RTDs/RTD.html
if (T==0) return(R0)
if (!is.numeric(alpha)) alpha <- 0.00672
R0 *(1 + alpha*T )
}

RTDnickelTemperatureFromAlpha <- function(R0, Rt, alpha=NA)
{
# simplified equation for Nickel RTD temperature
# R0 resistance at 0C
# Rt resistance at  temperature T
# alpha (optional) resistance's temperature coefficient
# http://www.capgo.com/Resources/Temperature/RTDs/RTD.html
if (T==0) return(R0)
if (!is.numeric(alpha)) alpha <- 0.00672
(Rt / R0 - 1) / alpha
}

RTDnickelIronResistanceFromAlpha <- function(R0, T, alpha=NA)
{
# simplified equation for Nickel-Iron RTD resistance
# R0 resistance at 0C
# T temperature in C
# alpha (optional) resistance's temperature coefficient
# http://www.capgo.com/Resources/Temperature/RTDs/RTD.html
if (T==0) return(R0)
if (!is.numeric(alpha)) alpha <- 0.00518
R0 *(1 + alpha*T )
}

RTDnickelIronTemperatureFromAlpha <- function(R0, Rt, alpha=NA)
{
# simplified equation for Nickel-Iron RTD temperature
# R0 resistance at 0C
# Rt resistance at  temperature T
# alpha (optional) resistance's temperature coefficient
# http://www.capgo.com/Resources/Temperature/RTDs/RTD.html
if (T==0) return(R0)
if (!is.numeric(alpha)) alpha <- 0.00518
(Rt / R0 - 1) / alpha
}

RTDmolybdenumResistanceFromAlpha <- function(R0, T, alpha=NA)
{
# simplified equation for Molybdenum RTD resistance
# R0 resistance at 0C
# T temperature in C
# alpha (optional) resistance's temperature coefficient
# http://www.capgo.com/Resources/Temperature/RTDs/RTD.html
if (T==0) return(R0)
if (!is.numeric(alpha)) alpha <- 0.00518
R0 *(1 + alpha*T )
}

RTDmolybdenumTemperatureFromAlpha <- function(R0, Rt, alpha=NA)
{
# simplified equation for Molybdenum RTD temperature
# R0 resistance at 0C
# Rt resistance at  temperature T
# alpha (optional) resistance's temperature coefficient
# http://www.capgo.com/Resources/Temperature/RTDs/RTD.html
if (T==0) return(R0)
if (!is.numeric(alpha)) alpha <- 0.00518
(Rt / R0 - 1) / alpha
}

RTDtemperatureFit <- function(R, R0, fitRTD='linear', alpha=0.00385)
{
# Mosaic Industries, Inc., 2014
# Relating resistance to temperature
# http://www.mosaic-industries.com/embedded-systems/microcontroller-projects/temperature-measurement/platinum-rtd-sensors/resistance-calibration-table
if (fitRTD=='linear'){
return((R/R0 - 1) / alpha)
}
if (fitRTD=='quadratic'){
return(-244.83 + 2.3419 * R + 0.0010664 * R^2)
}
if (fitRTD=='cubic'){
return(-247.29 + 2.3992 * R + .00063962 * R^2 + 1.0241E-6 * R^3)
}
if (fitRTD=='polynomial'){
c0 <- -245.19
c1 <- 2.5293
c2 <- -0.066046
c3 <- 4.0422E-3
c4 <- -2.0697E-6
c5 <- -0.025422
c6 <- 1.6883E-3
c7 <- -1.3601E-6
return(c0+R*(c1+R*(c2+R*(c3+c4*R)))/(1+R*(c5+R*(c6+c7*R))))
}
}

DiameterAWG <- function(AWG) 0.005*92^((36 - AWG)/39)
#American Wire Gauge (AWG) diameter from AWG number

ThermocoupleEquationTypeKrationalPolynomial<-function(vV, thermocoupleType='k')
{
# Mosaic Industries, Inc., 2014
# rational polynomial function approximation for Type K thermocouples
# http://www.mosaic-industries.com/embedded-systems/microcontroller-projects/temperature-measurement/thermocouple/calibration-table#computing-cold-junction-voltages
sapply(vV,function(x){
if (thermocoupleType=='b') if (x>=2.431) w<-2 else w<-1
if (thermocoupleType=='e') if (x>=53.112) w<-5 else if (x>=24.964) w<-4 else if (x>=0.591) w<-3 else if (x>=-5.237) w<-2 else w<-1
if (thermocoupleType=='j') if (x>=57.953) w<-5 else if (x>=45.494) w<-4 else if (x>=21.840) w<-3 else if (x>=0) w<-2 else w<-1
if (thermocoupleType=='k') if (x>=3.327500e+01) w<-5 else if (x>=1.639700e+01) w<-4 else if (x>=4.096000e+00) w<-3 else if (x>=-3.554000e+00 ) w<-2 else w<-1
if (thermocoupleType=='n') if (x>=20.613) w<-3 else if (x>=0) w<-2 else w<-1
if (thermocoupleType=='r') if (x>=14.277) w<-4 else if (x>=7.461) w<-3 else if (x>=1.469) w<-2 else w<-1
if (thermocoupleType=='s') if (x>=12.856 ) w<-4 else if (x>=6.913) w<-3 else if (x>=1.441) w<-2 else w<-1
if (thermocoupleType=='t') if (x>= 9.288) w<-4 else if (x>=0) w<-3 else if (x>=-4.648) w<-2 else w<-1
T0 <- thermocouple::thermocoupleCoefficientsTypeKrationalPolynomial[5,w]
v0 <- thermocouple::thermocoupleCoefficientsTypeKrationalPolynomial[6,w]
p1 <- thermocouple::thermocoupleCoefficientsTypeKrationalPolynomial[7,w]
p2 <- thermocouple::thermocoupleCoefficientsTypeKrationalPolynomial[8,w]
p3 <- thermocouple::thermocoupleCoefficientsTypeKrationalPolynomial[9,w]
p4 <- thermocouple::thermocoupleCoefficientsTypeKrationalPolynomial[10,w]
q1 <- thermocouple::thermocoupleCoefficientsTypeKrationalPolynomial[11,w]
q2 <- thermocouple::thermocoupleCoefficientsTypeKrationalPolynomial[12,w]
q3 <- thermocouple::thermocoupleCoefficientsTypeKrationalPolynomial[13,w]
v<-x
T0+(v-v0)*(p1+(v-v0)*(p2+(v-v0)*(p3+p4*(v-v0)))) / (1 + (v-v0)*(q1+(v-v0)*(q2+q3*(v-v0))))
})
}

ThermocoupleEquationTemperatureToVoltage<-function(vT, thermocoupleType='k')
{
# Computing cold junction voltages
# http://www.mosaic-industries.com/embedded-systems/microcontroller-projects/temperature-measurement/thermocouple/calibration-table#computing-cold-junction-voltages
sapply(vT,function(x){
if (thermocoupleType=='b') w<-2
if (thermocoupleType=='e') w<-3
if (thermocoupleType=='j') w<-4
if (thermocoupleType=='k') w<-5
if (thermocoupleType=='n') w<-6
if (thermocoupleType=='r') w<-7
if (thermocoupleType=='s') w<-8
if (thermocoupleType=='t') w<-9
T0 <- thermocouple::thermocoupleColdJunctionVoltageCoeff[6,w]
v0 <- thermocouple::thermocoupleColdJunctionVoltageCoeff[7,w]
p1 <- thermocouple::thermocoupleColdJunctionVoltageCoeff[8,w]
p2 <- thermocouple::thermocoupleColdJunctionVoltageCoeff[9,w]
p3 <- thermocouple::thermocoupleColdJunctionVoltageCoeff[10,w]
p4 <- thermocouple::thermocoupleColdJunctionVoltageCoeff[11,w]
q1 <- thermocouple::thermocoupleColdJunctionVoltageCoeff[12,w]
q2 <- thermocouple::thermocoupleColdJunctionVoltageCoeff[13,w]
T<-x
v0+((T-T0)*(p1+(T-T0)*(p2+(T-T0)*(p3+p4*(T-T0)))))/(1+(t-T0)*(q1+q2*(T-T0)))
})
}

RTDtemperatureFromResistance <- function(R, R0)
{
# R Sensor's resistance (ohm)
# R0 Resistance at temperature 0C
# R Sensor resistance
# Resistive temperature detectors PTxx
# www.madur.com
if (R>=R0){
A <- 3.9083
B <- -5.775E-7 
C <- -4.183E-12
(-R0*A+sqrt(R0^2*A^2-4*R0*B*(R0-R))) / (2*R0*B)
} else {
R100 <- R/R0*100
p <- -5.67E-6
q <- 2.4984E-2
r <- 2.22764 
s <- -242.078
p*R^3 + q*R^2 + r*R100 + s
}
}

ThermocoupleLeadWireExternalResistanceUS <- function(thermocoupleType, thermocoupleLength, thermocoupleGauge, leadWireType, leadWireLength, leadWireGauge) 
{
t1 <- 0
if (tolower(thermocoupleType)=='k') t1 <- 3 else if (tolower(thermocoupleType)=='j') t1 <- 4 else if (tolower(thermocoupleType)=='t') t1 <- 5 else if (tolower(thermocoupleType)=='e') t1 <- 6 else if (tolower(thermocoupleType)=='n') t1 <- 7 else if (tolower(thermocoupleType)=='s') t1 <- 8 else if (tolower(thermocoupleType)=='r') t1 <- 9
if (t1==0) stop('Wrong thermocouple Type')
t2 <- 0
if (tolower(leadWireType)=='k') t2 <- 3 else if (tolower(leadWireType)=='j') t2 <- 4 else if (tolower(leadWireType)=='t') t2 <- 5 else if (tolower(leadWireType)=='e') t2 <- 6 else if (tolower(leadWireType)=='n') t2 <- 7 else if (tolower(leadWireType)=='s') t2 <- 8 else if (tolower(leadWireType)=='r') t2 <- 9
if (t2==0) stop('Wrong lead wire Type')
r1 <- which(thermocouple::thermocoupleWireSizeResistanceImperial[["AWGNo"]]==thermocoupleGauge)
r2 <- which(thermocouple::thermocoupleWireSizeResistanceImperial[["AWGNo"]]==leadWireGauge)
thermocoupleLength * thermocouple::thermocoupleWireSizeResistanceImperial[r1,t1] + leadWireLength * thermocouple::thermocoupleWireSizeResistanceImperial[r2,t2] 
}

ThermocoupleTable10colsTo2 <- function(thermocoupleTable)
{
# converts the thermocouple table from n X 12 to m X 2
n <- dim(thermocoupleTable)[1]
T1 <- thermocoupleTable[1, "TempC",]
Tn <- thermocoupleTable[n, "TempC",]
thrK <- cbind(T1:Tn,0)
m <- which( seq(T1, Tn, by=10) %in% thermocoupleTable[["TempC"]])
mL <- length(m)
thrK[seq(1, n*10-19, by=10),2] <- thermocoupleTable[m, "T0C"]
thrK[seq(2, n*10-19, by=10),2] <- thermocoupleTable[m[-mL], "T1C"]
thrK[seq(3, n*10-19, by=10),2] <- thermocoupleTable[m[-mL], "T2C"]
thrK[seq(4, n*10-19, by=10),2] <- thermocoupleTable[m[-mL], "T3C"]
thrK[seq(5, n*10-19, by=10),2] <- thermocoupleTable[m[-mL], "T4C"]
thrK[seq(6, n*10-19, by=10),2] <- thermocoupleTable[m[-mL], "T5C"]
thrK[seq(7, n*10-19, by=10),2] <- thermocoupleTable[m[-mL], "T6C"]
thrK[seq(8, n*10-19, by=10),2] <- thermocoupleTable[m[-mL], "T7C"]
thrK[seq(9, n*10-19, by=10),2] <- thermocoupleTable[m[-mL], "T8C"]
thrK[seq(10, n*10-19, by=10),2] <- thermocoupleTable[m[-mL], "T9C"]
thrK
}

SensorSensitivity <- function(T1, E1, T2, E2) (T1-T2)/(E1-E2)
# Sensitivity of the sensor
# T1 measured temperature
# E1 resistance (platinum sensor) or the thermoelectric emf (thermocouple) for T1
# T2 measured temperature
# E2 resistance (platinum sensor) or the thermoelectric emf (thermocouple) for T2
# Gerd Scheller, 2003
# Error Analysis of a Temperature Measurement System 
# with worked examples
# JUMO, FAS 625, Edition 06.03

SelfHeatingError <- function(I, R, Ek) I^2*R/Ek
# self-heating error
# I intensity (A)
# R resistance (ohm)
# EK self-heating coefficient(mW/C)
# Gerd Scheller, 2003
# Error Analysis of a Temperature Measurement System 
# with worked examples
# JUMO, FAS 625, Edition 06.03

TminusT90CCT2008 <- function(T90K){
# T - T90 computed by a polynomial (CCT WG4 2008)
# Franco Pavese and Gianfranco Molinar Min Beciet, 2013
# Modern Gas-Based Temperature and Pressure Measurements
# Springer Science + Business Media
# pp. 42
sapply(T90K,function(x){
d <- c(-14.0651, 40.9970, -44.1079, 16.5315)
if ((x>=0.65) & (x<1)) return(sum(d * x^(0:3))/1e3)
a <- c(8.7999, -54.8216, 101.4590, -83.5816, 32.2307, -4.7513)
if ((x>=1) & (x<2)) return(sum(a * x^(0:5))/1e3)
if (x>=2 & x<8) return(0)
b <- c(4.42457E1, -1.76311E2, -1.53985E3, -3.63685E3, -4.19898E3, -2.61319E3, -8.41922E2, -1.10322E2)
if (x>=8 & x<273.16) return(sum(b * (log10(x/273.16))^(1:8))/1e3)
C <- c( 0.0497, -0.3032, 1.0254, -1.2895, 0.5176)
if (x>=273.16 & x<1358) return(x*sum(C * (273.16/x)^(2*0:4))/1e3)
})
}

TminusT90Pavese4CubicPolynomials<- function(T90K){
# T - T90 computed by a polynomial (CCT WG4 2008)
# Franco Pavese and Gianfranco Molinar Min Beciet, 2013
# Modern Gas-Based Temperature and Pressure Measurements
# Springer Science + Business Media
# pp. 42
sapply(T90K,function(x){
if ((x>=2.3) & (x<90)) {
a0 <- -0.000593965
a1 <- 0.000084027
a2 <- -0.000003974
a3 <- 0.000000052
h <- x - 2.3
return(a0 + h*(a1 + h/2 * (a2 * h + 1/3 * a3 * h^2 ))/1e3)
}
if ((x>=90) & (x<300)) {
a0 <- -0.002642255
a1 <- -0.000063857
a2 <- 0.000000602
a3 <- 0.000000004
h <- x - 90
return(a0 + h*(a1 + h/2 * (a2 * h + 1/3 * a3 * h^2 ))/1e3)
}
if ((x>=300) & (x<450)) {
a0 <- 0.003884825
a1 <- 0.000157792
a2 <- -0.000001213
a3 <- -0.000000002
h <- x - 300
return(a0 + h*(a1 + h/2 * (a2 * h + 1/3 * a3 * h^2 ))/1e3)
}
if ((x>=450) & (x<=1238)) {
a0 <- 0.012749953
a1 <- -0.000047281
a2 <- 0.000000582
a3 <- -0.000000001
h <- x - 450
return(a0 + h*(a1 + h/2 * (a2 * h + 1/3 * a3 * h^2 ))/1e3)
}
})
}

TminusT90Pavese6CubicPolynomials<- function(T90K){
# T - T90 computed by a polynomial (CCT WG4 2008)
# Franco Pavese and Gianfranco Molinar Min Beciet, 2013
# Modern Gas-Based Temperature and Pressure Measurements
# Springer Science + Business Media
# pp. 42
sapply(T90K,function(x){
if ((x>=2.3) & (x<30)) {
a0 <- 0.000302995
a1 <- -0.000113512
a2 <- 0.00001565
a3 <- -0.000000743
h <- x - 2.3
return(a0 + h*(a1 + h/2 * (a2 * h + 1/3 * a3 * h^2 ))/1e3)
}
if ((x>=30) & (x<90)) {
a0 <- 0.000530791
a1 <- 0.000034942
a2 <- -0.000004931
a3 <- 0.000000099
h <- x - 30
return(a0 + h*(a1 + h/2 * (a2 * h + 1/3 * a3 * h^2 ))/1e3)
}
if ((x>=90) & (x<273.16)) {
a0 <- -0.002667189
a1 <- -0.000081847
a2 <- 0.000001038
a3 <- 0.000000001
h <- x - 90
return(a0 + h*(a1 + h/2 * (a2 * h + 1/3 * a3 * h^2 ))/1e3)
}
if ((x>=273.16) & (x<380)) {
a0 <- 0.000458927
a1 <- 0.000119811
a2 <- 0.000000924
a3 <- -0.000000038
h <- x - 273.16
return(a0 + h*(a1 + h/2 * (a2 * h + 1/3 * a3 * h^2 ))/1e3)
}
if ((x>=380) & (x<800)) {
a0 <- 0.010793868
a1 <- 0.000001234
a2 <- 0.000000046
a3 <- 0.000000001
h <- x - 380
return(a0 + h*(a1 + h/2 * (a2 * h + 1/3 * a3 * h^2 ))/1e3)
}
if ((x>=800) & (x<=1238)) {
a0 <- 0.022102959
a1 <- 0.000068722
a2 <- 0.000000276
a3 <- -0.000000002
h <- x - 800
return(a0 + h*(a1 + h/2 * (a2 * h + 1/3 * a3 * h^2 ))/1e3)
}
})
}

ThermistorApproxDriftResistance <- function(Ri, T, a, b) Ri + Ri * (a+b*log(T))/100
# Approximation of Drift Resistance of NTC Thermistors
# Rt resistance at time T
# Ri initial resistance
# T aging time
# a intercept at T=1
# b slope (%deltaR per decade of time T)
# Quality Thermistor, Inc. 2108
# http://www.cornerstonesensors.com/About.asp?PageCode=Stability&Print=Page

ThermistorApproxDriftTime <- function(Ri, Rt, a, b) exp(((Rt - Ri)/Ri * 100 - a)/b)
# Approximation of Drift Time of NTC Thermistors
# T aging time
# Rt resistance at time T
# Ri initial resistance
# a intercept at T=1
# b slope (%deltaR per decade of time T)
# Quality Thermistor, Inc. 2108
# http://www.cornerstonesensors.com/About.asp?PageCode=Stability&Print=Page

BimaterialStripCurvatureRadiusFromTemperature <- function(T0, R0, T, m, n, alpha1, alpha2, thickn) 1/((6*(1+m)^2*(alpha2-alpha1)*(T-T0))/(thickn*(3*(1+m)^2*(1+m*n)*(m^2+1/m*n)))+1/R0)
BimaterialStripTemperatureFromCurvatureRadius <- function(R0, T0, R, m, n, alpha1, alpha2, thickn) (thickn*(3*(1+m)^2*(1+m*n)*(m^2+1/m*n)))*(1/R-1/R0)/(6*(1+m)^2*(alpha2-alpha1))+T0
# curvature radius of a bimetallic strip uniformly heated from T0 to T in the absence of external forces 
# 1/R0 Initial curvature of the strip at temperature T0
# alpha2, alpha1 Coefficients of expansion of the two elements
# n E1/E2, with E1 and E2 their respective Young’s moduli
# m t1/t2, with t1 and t2 their respective thicknesses
# thickn t1 + t2 thickness of the strip
# John G. Webster, 1999
# The Measurement, Instrumentation and Sensors Handbook
# CRC Press LLC, eq (32.1)



SplineEval<-function(x, knotsK, coeffsC){
# Spline algorithm used in The Observed Properties of Liquid Helium at the Saturated Vapor Pressure
# Donnelly, Donnelly and Hills [J. Low Temp. Phys. 44, 471 (1981)]
K=knotsK
C=coeffsC
sapply(x,function(X){
Ncap7=length(K)
N<-Ncap7 - 8
if (X<K[4] | (X>K[N+5])) stop('The Value of Temperature you entered is ouside the fitted range')
if ((X>=K[4]) & (X<K[Ncap7-3])){
J1 =0
J=Ncap7-7
while((J-J1) > 1){
L=(J1+J)/2
if (X>=K[L+4]) J1 =L else J=L
}
K1=K[J+1]
K2=K[J+2]
K3=K[J+3]
K4=K[J+4]
K5=K[J+5]
K6=K[J+6]
E2=X-K2
E3=X-K3
E4=K4-X
E5=K5-X
C11=((X-K1)*C[J+1]+E4*C[J])/(K4-K1)
Cd11=(C[J+1]-C[J])/(K4-K1)
C21=(E2*C[J+2]+E5*C[J+1])/(K5-K2)
Cd21=(C[J+2]-C[J+1])/(K5-K2)
C31=(E3*C[J+3]+(K6-X)*C[J+2])/(K6-K3)
Cd31=(C[J+3]-C[J+2])/(K6-K3)
C12=(E2*C21+E4*C11)/(K4-K2)
Cd12=(C21+E2*Cd21-C11+E4*Cd11)/(K4-K2)
Cdd12=2*(Cd21-Cd11)/(K4-K2)
C22=(E3*C31+E5*C21)/(K5-K3)
Cd22=(C31+E3*Cd31-C21+E5*Cd21)/(K5-K3)
Cdd22=2*(Cd31-Cd21)/(K5-K3)
S=(E3*C22+E4*C12)/(K4-K3)
Sd=(E3*Cd22+C22+E4*Cd12-C12)/(K4-K3)
Sdd=(E3*Cdd22+2*Cd22+E4*Cdd12-2*Cd12)/(K4-K3)
} else stop('Error!')
S
})
}


