#  "KANNAN CURRENT VOLTAGE PROPERTY OF GaAs MESFET -  870818"
# 3030 LET M=40: REM 40 data points for this particular data set
??? Try to find explanation of this -- too messy as is.
L9<-40
# m<-readline("How many points to use in fit ")
m<-40
if (m>L9) stop("Not enough data points")
#  LET N=6: REM 6 parameters
# 3041 PRINT "b1 = Idss, b2 = Rs, b3 = Rd, b4 = Vp0"
# 3043 PRINT "   b5 = gamma, b6 = alpha"
# 3070 FOR I=1 TO L9
# 3080 READ Y(I,1), Y(I,2), Y(I,3): REM Vds, Ids, Vgs
# 3090 NEXT I
#  0.0001<=b1<=1 # Idss
#  0.001<=b2<=100 #  Rs
#  0.0001<=b3<=100 # REM positive Rd
#  -100<=b4<=-0.001 #  REM NEGATIVE Vp0
#  0.0001<=b5<=1 : REM ALPHA
#  0.0001<=b6<=1 : REM GAMMA
#  3240 FOR J=1 TO 5: REM CASE SPECIFIC
#  3245 LET Z1=0: REM MAXIMUM SO FAR
#  3250 FOR I=1 TO M/5
#  3255 LET L=(J-1)*M/5+I
#  3260 IF Y(L,2)>Z1 THEN LET Z1=Y(L,2)
#  3265 NEXT I
#  3267 LET Z1=1/SQR(M*Z1)
#  3270 FOR I=1 TO M/5
#  3275 LET L=(J-1)*M/5+I
#  3280 LET Y(L,4)=Z1
#  3295 NEXT I: REM Y(I,4) HOLDS 1/SQRT(M*MAX CURRENT IN SET)
#  3300 NEXT J


lo<-c(0.0001, 0.001, 0.0001, -100, 0.0001, 0.0001)
up<-c(1, 100, 100, -0.001, 1, 1)



start<-c(b1=0.100, b2=2, b3=2, b4=-3, b5=-0.2, b6=2)
y1<-c(0.2, 0.4, 0.8, 1, 1.6, 2, 3, 4, 0.2, 0.4, 0.8, 1, 
     1.6, 2, 3, 4, 0.2, 0.4, 0.8, 1, 1.6, 2, 3, 4, 0.2, 
     0.4, 0.8, 1, 1.6, 2, 3, 4, 0.2, 0.4, 0.8, 1, 1.6, 2, 3, 4)
y2<-c(0.027, 0.052, 0.09, 0.098, 0.1, 0.099, 0.096, 0.092, 
     0.024, 0.045, 0.076, 0.081, 0.083, 0.082, 0.08, 0.078, 
     0.02, 0.038, 0.061, 0.063, 0.065, 0.065, 0.064, 0.064, 
     0.016, 0.029, 0.043, 0.045, 0.047, 0.047, 0.048, 0.049, 
     0.011, 0.019, 0.027, 0.028, 0.03, 0.031, 0.033, 0.035)
y3<-c(0, 0, 0, 0, 0, 0, 0, 0, -0.5, -0.5, -0.5, -0.5, -0.5,
     -0.5, -0.5, -0.5, -1, -1, -1, -1, -1, -1, -1, -1, 
     -1.5, -1.5, -1.5, -1.5, -1.5, -1.5, -1.5, -1.5, 
     -2, -2, -2, -2, -2, -2, -2, -2)
kandat<-data.frame(y1=y1, y2=y2, y3=y3)

formula <- "y 

3500 LET I3=0: REM start with computable function
3510 LET Z1=Y(I,1)-Y(I,2)*(b2+b3): REM VDS-ADJUSTED
3515 LET Z2=b4+b5*Z1
3520 IF ABS(Z2)<.00001 THEN 3800: REM AVOID ZERO DIVIDE
3530 LET Z3=Y(I,3)-Y(I,2)*b2-Z2
3540 IF ABS(Z3)<0.00001 THEN 3830
3550 LET Z4=b6*(Y(I,1)-Y(I,2)*(b2+b3))
3560 LET Z5=Z4/Z3
3565 IF ABS(Z5)>30 THEN 3860: REM AVOID LARGE EXPONENT
3570 LET Z6=EXP(Z5)
3575 LET Z7=1/Z6
3580 LET Z8=(Z6-Z7)/(Z6+Z7)
3590 LET Z9=1-(Y(I,3)-Y(I,2)*b2)/Z2
3600 LET Z0=b1*Z9*Z9*Z8
3610 LET Y(I,5)=Z0: REM SAVE MODEL VALUE
3620 LET R1=(Z0-Y(I,2))*Y(I,4)
3625 LET YEXP=Z0*Y(I,4): REM note scaling
3630 LET Y(I,6)=R1: REM SAVE RESIDUAL
3640 RETURN
3800 PRINT "Vp0 + gamma*(Vds-Id*(Rs+Rd)) near zero =";Z2;
3802 PRINT " data  point ";I
3810 PRINT #3,"Vp0 + gamma*(Vds-Id*(Rs+Rd)) near zero =";Z2;
3812 PRINT #3,"  data point ";I
3820 GOTO 3950
3830 PRINT "Denominator of tanh argument near zero = ";Z3;
3832 PRINT "  data point ";I
3840 PRINT #3,"Denominator of tanh argument near zero = ";Z3;
3842 PRINT #3,"  data point ";I
3850 GOTO 3950
3860 PRINT "Argument of exp() too large = ";Z5;" data point ";I
3870 PRINT #3,"Argument of exp() too large = ";Z5;" data point  ";I
3880 GOTO 3950
3950 LET I3=1: REM cannot compute the residual
3960 RETURN
4000 PRINT "DERIVATIVES NOT DEFINED"
4010 PRINT #3,"DERIVATIVES NOT DEFINED"
4030 STOP
4500 PRINT "Model for small Vds"
4520 PRINT #3,"Model for small Vds"
4530 PRINT "Data","Model","Residual"
4540 PRINT #3,"Data","Model","Residual"
4550 FOR I=1 TO 40 STEP 8
4560 PRINT Y(I,2),Y(I,5),Y(I,6)
4570 PRINT #3,Y(I,2),Y(I,5),Y(I,6)
4580 NEXT I
4590 PRINT
4600 PRINT #3,
4633 PRINT "number of data points available =";L9
4634 PRINT #3,"number of data points available =";L9
4635 INPUT "How many points to use in fit ";M
4636 IF M>L9 THEN 3035
4637 PRINT #3,"How many points to use in fit ";M
4640 LET N=6: REM 3 parameters
4650 RETURN
