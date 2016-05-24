`PreSet.Instr` <-
function()
  {
    ## set up a list of instruments from passcal Guralp
###  usage    Kal = PreSet.Instr()

    types = c("40T", "3T", "L28", "LE3D20s", "GEOSP1", "CMG3ESPC" ,  "60T", "TRIL120" )
    
    K = as.list(types)


    CMG40TPASS=c(
      "ZEROS 2",
      "0.0000E+00 0.0000E+00",
      "0.0000E+00 0.0000E+00",
      "POLES 3",
      "-0.1480E+00 0.1480E+00",
      "-0.1480E+00 -0.1480E+00",
      "-50.0 0.0",
      "CONSTANT 1.0",
      "SENSE 800")

    CMG40T=c(
      "ZEROS 3",
      "0.0000E+00 0.0000E+00",
      "0.0000E+00 0.0000E+00",
      "159.0 0.0",
      "POLES 3",
      "-0.1480E+00 0.1480E+00",
      "-0.1480E+00 -0.1480E+00",
      "-50.0 0.0",
      "CONSTANT -1.97292",
      "SENSE 800")

    CMG60T.greg=c(
      "ZEROS 2",
      "0 0",
      "0 0",
      "POLES 5",
      "-0.0740159 0.0740159",
      "-0.0740159 -0.0740159",
      "-1005.31 0",
      "-502.655 0",
      "-1130.97 0",
      "CONSTANT 1.115583e+12",
      "SENSE 2000")

    
    CMG3ESPC=c(
      "ZEROS 2",
      "0 0",
      "0 0",
      "POLES 5",
      "-11.78E-3  11.78E-3",
      "-11.78E-3 -11.78E-3",
      "-160 0",
      "-80  0",
      "-180 0",
      "CONSTANT 14476458.95",
      "SENSE 1962")
    
    CMG60T=c(
      "ZEROS 2",
      "0 0",
      "0 0",
      "POLES 5",
      "-11.78E-3  11.78E-3",
      "-11.78E-3 -11.78E-3",
      "-160 0",
      "-80  0",
      "-180 0",
      "CONSTANT 14476458.95",
      "SENSE 1962")
    

    L28=c(   
      "ZEROS 2",
      "0.0000E+00 0.0000E+00",
      "0.0000E+00 0.0000E+00",
      "POLES 2",
      "-19.99E+00 19.99E+00",
      "-19.99E+00 -19.99E+00",
      "CONSTANT 1",
      "SENSE 30.4")
    
    CMG3T=c(   
      "ZEROS 3",
      "0.0000E+00 0.0000E+00",
      "0.0000E+00 0.0000E+00",
      "9.2050E+02 0.0000E+00",
      "POLES 4",
      "-3.701E-02 +3.701E-02",
      "-3.701E-02 -3.701E-02",
      "-4.599E+02 +2.362E+02",
      "-4.599E+02 -2.362E+02",
      "CONSTANT -2.904E+02",
      "SENSE 1500")

    LE3D20s =c(
         "ZEROS 3",
         "0.0000E+00 0.0000E+00",
         "0.0000E+00 0.0000E+00",
         "0.0000E+00 0.0000E+00",     
         "POLES 3",    
         "-0.22E+00 +0.235E+00",
         "-0.22E+00 -0.235E+00",
         "-0.23E+00 0.000E+00",
         "CONSTANT 1",
         "SENSE 1000")

 GEOSP1=c(
   "ZEROS 3",
   "0.0000E+00 0.0000E+00",
   "0.0000E+00 0.0000E+00",
   "0.0000E+00 0.0000E+00",     
   "POLES 3",    
   "-0.5 0.8660254",
   "-1.0 0.0000000",
   "-0.5 -0.8660254",
   "CONSTANT 1", 
   "SENSE 700.7874"
   )
    
TRIL120 = c("ZEROS 5",
  "0 0",
  "0 0",
  "-90.0 0",
  "-160.7 0",
  "-3108 0",
  "POLES 7",
  "-0.03852 0.03658",
  "-0.03852 -0.03658",
  "-178 0",
  "-135 160",
  "-135 -160",
  "-671 1154",
  "-671 -1154",
  "CONSTANT 308000",
  "SENSE 1201")

    ###########   this is a filter I derived using
    ########    specs and matlab (yuk!)


#### library(signal)
    
####  -0.5+0.8660254i -1.0+0.0000000i -0.5-0.8660254i
####  bf = butter(3, 1 , type = "high", plane = "s")
####    as.Zpg(bf)

    GEOSP12=c(
      "ZEROS 1",
      "1.0000E+00 0.0000E+00",
      "POLES 1",
      "-0.939E+00 0.000E+00",
      "CONSTANT 1",
      "SENSE 700.7874"
    )
    
    GEOSP0=c(
    "ZEROS 0",
    "POLES 0", 
    "CONSTANT 1",
    "SENSE 700.7874"
    )
    

    
    K[[1]] = ReadSet.Instr(CMG40T)
    K[[2]] = ReadSet.Instr(CMG3T)
    K[[3]] = ReadSet.Instr(L28)
    K[[4]] = ReadSet.Instr(LE3D20s)
    K[[5]] = ReadSet.Instr(GEOSP1)
    K[[6]] = ReadSet.Instr(CMG3ESPC)
    K[[7]] = ReadSet.Instr(CMG60T)
    K[[8]] = ReadSet.Instr(TRIL120)

    names(K) = types
    ##  Others?
    return(K)
  }

