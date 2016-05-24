# Constants object
getConstants = function() {
  a = list()
  
  # Board measurements
  a$R1 = 6.35  # center to DB wire
  a$R2 = 15.9  # center to SB wire
  a$R3 = 99    # center to inner triple ring
  a$R4 = 107   # center to outer triple ring
  a$R5 = 162   # center to inner double ring
  a$R = 170    # center to outer double ring
  
  # Dartboard scores arrangement, clockwise starting at top center
  # Standard arrangement
  a$standard = c(20,1,18,4,13,6,10,15,2,17,3,19,7,16,8,11,14,9,12,5) 

  # Curtis arrangement
  a$curtis = c(20,1,19,3,17,5,15,7,13,9,11,10,12,8,14,6,16,4,18,2)

  # Linear arrangement
  a$linear = c(20,19,18,17,16,15,14,13,12,11,10,9,8,7,6,5,4,3,2,1)

  # What is the arrangement of scores being used?
  a$S = a$standard

  return(a)
}

                                       



