`SENSORsensitivity` <-
function(K=1)
  {
    
    ###  instruments
    ### AI = c("40T", "3T",  "LD",  "MC",  "EL")
    ###  1 2 101 102 103 104

    ###  MC = 80 mV/Pa  = 0.08 V/Pa
    ###  for microphones see:
    ### Johnson, J.B., R.C. Aster, M.C. Ruiz, S.D. Malone, P.J. McChesney,
    ### J.M. Lees, and P.R. Kyle, Interpretation and utility of infrasonic
    ### records from erupting volcanoes, J. Volc. Geoth. Res., 121 (1-2), 15-63, 2003. 

    #### 40T = CMG40T Guralp
    #### 3T =  CMG3T  Guralp
    #### L28 =  L28 Mark Products
    #### LD  = Larson-Davis Infrasonic Microphone
    ####  EL   = Electret Microphone
    ####  MC = McChesney Microphone
    ####  EL(SANGAY) = Electret Microphone used at Sangay 1998


    
    
    codes = c(1,    2,      3,  101,    102,    103,    104,   200 )
    
    AI = c("40T", "3T",  "L28", "LD",    "EL", "MC", "EL(SANGAY)", "TRIL120")
    
    II = c(0.8/1000, 1.5/1000, .0304/1000,   0.04841724,  .1,  0.08, 0.03, 1.201/1000  )
    
    unit = c("mu m/s", "mu m/s", "mu m/s", "Pa", "Pa","Pa","Pa", "mu m/s" )
    ###   updir 1 = positive is up motion, -1 positive is down motion
    
    updir  = c(1, 1, -1, 1, 1, 1, 1, 1 )

    if(missing(K))
      {
        print(data.frame(cbind(codes, AI, II, unit)))
        return(NULL)
      }

    i = match(K, codes)
   return(list(code=codes[i], inst=AI[i], sense=II[i], units=unit[i], updir=updir[i]))
    
  }

