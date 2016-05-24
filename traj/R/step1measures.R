step1measures <-
function(Data, Time, ID = FALSE, verbose = TRUE)
{
  data = Data
  time = Time
  
  input.data = data
  input.time = time

  if(dim(data)[1] != dim(time)[1] || dim(data)[2] != dim(time)[2])
    stop("data and time must be the same size.")
  
  sample.size = dim(data)[1]
  
  # Deal with IDs
  if(ID){
   IDvector = data[,1] 
   data = data[,-1]
   time = time[,-1]
  }

  max.num.obs  = dim(data)[2]
  

  # Clean input tables  
  clean.data = matrix(ncol= max.num.obs, nrow=sample.size)
  clean.time = matrix(ncol= max.num.obs, nrow=sample.size)
  
  num.obs = rep(999,sample.size)
  less.than.4.obs = NULL  

  # Check for appropriate amount of usable data
  for(i_sample in 1:sample.size)
  {
      real.obs.pos = which(!is.na(data[i_sample,]))
      num.obs[i_sample] = length(real.obs.pos)
      clean.data[i_sample,] = as.vector(c(unlist(data[i_sample,real.obs.pos]), rep(NA, max.num.obs - num.obs[i_sample])))
      clean.time[i_sample,] = as.vector(c(unlist(time[i_sample,real.obs.pos]), rep(NA, max.num.obs - num.obs[i_sample])))
      
      if(length(real.obs.pos) < 4)
        less.than.4.obs = c(less.than.4.obs, i_sample)
      
      clean.data.pos = which(!is.na(clean.data[i_sample,]))
      if(any(is.na(clean.time[i_sample,clean.data.pos])))
         stop(paste("There must be a time associated to every observation. Line: ", i_sample, sep = ""))
  }
  if(!is.null(less.than.4.obs)){
    clean.data = clean.data[-less.than.4.obs,]
    clean.time = clean.time[-less.than.4.obs,]
  }

  sample.size = nrow(clean.data)  

  # Generate new Id vector for usable data 
  if(ID){
    if(!is.null(less.than.4.obs))
      IDvector = IDvector[-less.than.4.obs]
  }else
    IDvector = seq(1:sample.size)

  # Setup output table
  data = clean.data
  time = clean.time

  output = data.frame(matrix(ncol = 25, nrow = sample.size))
  names = c("ID", "m1", "m2", "m3", "m4", "m5", "m6", "m7", "m8", "m9", "m10", "m11", "m12", "m13",
            "m14", "m15", "m16", "m17", "m18", "m19", "m20", "m21", "m22", "m23", "m24")
  colnames(output) = names
  
  output$ID = IDvector

  
  #################################
  # Analysis
  #################################
  
  ################################################################################################
  #1 RANGE (m1)
  ################################################################################################
  for(i in 1:sample.size) {
    output$m1[i]=max(data[i,], na.rm=TRUE)-min(data[i,], na.rm=TRUE)
  }
  
  ################################################################################################
  #2 MEAN-OVER-TIME (m2)
  ################################################################################################
  for(i in 1:sample.size) {
    output$m2[i]=mean(data[i,], na.rm=TRUE)
  }

  
  ################################################################################################
  #3 STANDARD DEVIATION (m3)
  ################################################################################################
  for(i in 1:sample.size) {
    output$m3[i]=sqrt(var(data[i,], na.rm=TRUE))
  }
  
  ################################################################################################
  #4 COEFFICIENT OF VARIATION (m4)
  ################################################################################################
  for(i in 1:sample.size) {
 
    output$m4[i]=100*output$m3[i]/output$m2[i]

  } 
  
  ################################################################################################
  #5 CHANGE (m5)
  ################################################################################################
  for(i in 1:sample.size) {
    output$m5[i]=last(data[i,], na.rm=TRUE)-first(data[i,], na.rm=TRUE)
  }
  
  ################################################################################################
  #6 MEAN CHANGE PER TIME UNIT (m6)
  ################################################################################################
  for(i in 1:sample.size) {
    output$m6[i]=(last(data[i,], na.rm=TRUE)-first(data[i,], na.rm=TRUE))/(last(time[i,], na.rm=TRUE)-first(time[i,], na.rm=TRUE)+1)
  }
  
  ################################################################################################
  #7 CHANGE RELATIVE TO THE FIRST SCORE (m7)
  ################################################################################################
  for(i in 1:sample.size) {

    output$m7[i]=(last(data[i,], na.rm=TRUE)-first(data[i,], na.rm=TRUE))/first(data[i,], na.rm=TRUE)
    #if(i == 5)print(first(data[i,], na.rm=TRUE))
  }
  
  ################################################################################################
  #8 CHANGE RELATIVE TO THE MEAN-OVER-TIME (m8)
  ################################################################################################
  for(i in 1:sample.size) {

    output$m8[i]=(last(data[i,], na.rm=TRUE)-first(data[i,], na.rm=TRUE))/output$m2[i]
  }
  
  ################################################################################################
  #9 SLOPE OF THE LINEAR MODEL (m9)
  ################################################################################################
  for(i in 1:sample.size) {
    b=coefficients(lm(data[i,]~time[i,]))  
    output$m9[i]=b[2]
  }

  
  ################################################################################################
  #10 PROPORTION OF VARIANCE EXPLAINED BY THE LINEAR MODEL y_i=a+bt_i+epsilon_i (m10)
  ################################################################################################
  for(i in 1:sample.size) {
    model=lm(data[i,]~time[i,])  
    r=resid(model)
    RSS=r%*%r
    Y=subset(data[i,], is.na(data[i,])==FALSE)
    m=length(Y)
    SYY=Y%*%Y-(sum(Y)^2)/m
    SSREG=SYY[1]-RSS[1]

    output$m10[i]=SSREG/SYY

  }
  
  ################################################################################################
  #11 MAXIMUM OF THE FIRST DIFFERENCE (m11)
  ################################################################################################
  # matrix of first difference FD
  FD=matrix(nrow=sample.size, ncol=max.num.obs-1)
  for(i in 1:sample.size) {
    for (j in 1:(max.num.obs-1)) {
        FD[i,j] = data[i,(j+1)]-data[i,j]
    }
  }
  
  for(i in 1:sample.size) {
    output$m11[i]=max(FD[i,], na.rm=TRUE)  		
  }
  
  ################################################################################################
  #12 SD OF THE FIRST DIFFERENCE (m12)
  ################################################################################################
  for(i in 1:sample.size) {
    output$m12[i]=sqrt(var(FD[i,], na.rm=TRUE))
  }  
  
  ################################################################################################
  #13 SD OF THE FIRST DIFFERENCE PER TIME UNIT (m13)
  ################################################################################################
  # Create a matrix of first difference per time unit
  FDunit=matrix(nrow=sample.size, ncol=max.num.obs-1)
  for(i in 1:sample.size) {
    for (j in 1:(max.num.obs-1)) {
      FDunit[i,j]=FD[i,j]/(time[i,j+1]-time[i,j])
    }
  }		
  
  for(i in 1:sample.size) {	
    output$m13[i]=sqrt(var(FDunit[i,], na.rm=TRUE))
  }
  
  ################################################################################################
  #14 MEAN OF THE ABSOLUTE FIRST DIFFERENCE (m14)
  ################################################################################################
  for(i in 1:sample.size) {
    output$m14[i]=mean(abs(FD[i,]), na.rm=TRUE) 
  }

  ################################################################################################
  #15 MAXIMUM OF THE ABSOLUTE FIRST DIFFERENCE (m15)
  ################################################################################################
  for(i in 1:sample.size) {
    output$m15[i]=max(abs(FD[i,]), na.rm=TRUE)
  }
  
  ################################################################################################
  #16 RATIO OF THE MAXIMUM ABSOLUTE FIRST DIFFERENCE TO THE MEAN-OVER-TIME (m16)
  ################################################################################################
  for(i in 1:sample.size) {

    output$m16[i]=output$m15[i]/output$m2[i]
  }
  
  ################################################################################################
  #17 RATIO OF THE MAXIMUM ABSOLUTE FIRST DIFFERENCE TO THE SLOPE (m17)
  ################################################################################################
  for(i in 1:sample.size) {

    output$m17[i]=output$m15[i]/output$m9[i]
  }
  
  ################################################################################################
  #18 RATIO OF THE SD OF THE FIRST DIFFERENCE TO THE SLOPE (m18)
  ################################################################################################
  for(i in 1:sample.size) {

    output$m18[i]=output$m12[i]/output$m9[i]
  }
  
  ################################################################################################
  #19 MEAN OF THE SECOND DIFFERENCE (m19)
  ################################################################################################
  # matrix of second difference SD
  SD=matrix(nrow=sample.size, ncol=max.num.obs-2)
  for(i in 1:sample.size) {
    for (j in 1:(max.num.obs-2)) {
      SD[i,j]=FD[i,(j+1)]-FD[i,j]
    }
  }

  for(i in 1:sample.size) {
    output$m19[i]=mean((SD[i,]), na.rm=TRUE)
  }
  
  ################################################################################################
  #20 MEAN OF THE ABSOLUTE SECOND DIFFERENCE (m20)
  ################################################################################################
  for(i in 1:sample.size) {
    output$m20[i]=mean(abs(SD[i,]), na.rm=TRUE)			
  }
  
  ################################################################################################
  #21 MAXIMUM OF THE ABSOLUTE SECOND DIFFERENCE (m21)
  ################################################################################################
  for(i in 1:sample.size) {
    output$m21[i]=max(abs(SD[i,]), na.rm=TRUE)
  }
  
  ################################################################################################
  #22 RATIO OF THE MAXIMUM ABSOLUTE SECOND DIFFERENCE TO THE MEAN-OVER-TIME (m22)
  ################################################################################################
  for(i in 1:sample.size) {
    
    output$m22[i]=output$m21[i]/output$m2[i]  		
  }
  
  ################################################################################################
  #23 RATIO OF THE MAXIMUM ABSOLUTE SECOND DIFFERENCE TO THE MEAN ABSOLUTE FIRST DIFFERENCE (m23)
  ################################################################################################
  # This measure creates problems is someone has a completely stable pattern solved with +0.0001
  for(i in 1:sample.size) {
    
    output$m23[i]=output$m21[i]/output$m14[i]
  }
  
  ################################################################################################
  #24 RATIO OF THE MEAN ABSOLUTE SECOND DIFFERENCE TO THE MEAN ABSOLUTE FIRST DIFFERENCE (m24)
  ################################################################################################
  # This measure creates problems is someone has a completely stable pattern solved with + 0.0001
  for(i in 1:sample.size) {
    
    output$m24[i]=mean(abs(SD[i,]), na.rm=TRUE)/output$m14[i]	
  }
  
  #########################
  # Deal with devide by 0  
  #########################

  # mean 0 
  temp.data = output$m2[(output$m2 != 0)]
  abs.temp.data = abs(temp.data)
  if(length(temp.data) == 0) # if this is the case, all your data has a mean of 0.
    mean.0 = 0.0001
  else
    mean.0 = temp.data[which.min(abs.temp.data)]/100
  # apply to calculations
  m4.na.pos = which(is.na(output$m4)| is.infinite(output$m4)) 
  if(length(m4.na.pos) != 0)
    output$m4[m4.na.pos]=100*output$m3[m4.na.pos]/mean.0

  m8.na.pos = which(is.na(output$m8)| is.infinite(output$m8)) 
  if(length(m8.na.pos) != 0){
    if(length(m8.na.pos) > 1)
      output$m8[m8.na.pos]=(apply(data[m8.na.pos,], 1, last, na.rm = TRUE) - apply(data[m8.na.pos,], 1, first,na.rm = TRUE))/mean.0
    else
      output$m8[m8.na.pos]=(last(data[m8.na.pos,], na.rm = TRUE) - first(data[m8.na.pos,],na.rm = TRUE))/mean.0
  }
  m16.na.pos = which(is.na(output$m16)| is.infinite(output$m16)) 
  if(length(m16.na.pos) != 0)
    output$m16[m16.na.pos]=output$m15[m16.na.pos]/mean.0
  
  m22.na.pos = which(is.na(output$m22)| is.infinite(output$m22)) 
  if(length(m22.na.pos) != 0)
    output$m22[m22.na.pos]=output$m21[m22.na.pos]/mean.0


  # y1 0
  temp.data = data[(data[,1] != 0),]
  abs.temp.data = abs(temp.data)
  if(nrow(temp.data) == 0) # if this is the case, the fist value of all your data is 0.
    y1.0 = 0.0001
  else
    y1.0 = temp.data[which.min(abs.temp.data)]/100

  m7.na.pos = which(is.na(output$m7) | is.infinite(output$m7))
  if(length(m7.na.pos) != 0){
    if(length(m7.na.pos) > 1)
      output$m7[m7.na.pos]= (apply(data[m7.na.pos,], 1, last, na.rm = TRUE) - apply(data[m7.na.pos,], 1, first, na.rm = TRUE)) / y1.0
    else
      output$m7[m7.na.pos]= (last(data[m7.na.pos,], na.rm = TRUE) - first(data[m7.na.pos,], na.rm = TRUE)) / y1.0
  }
  
  
  # SYY 0
  form10=vector(length=sample.size)
  for(i_test in 1:sample.size) {
    model=lm(data[i_test,]~time[i_test,])  
    r=resid(model)
    RSS=r%*%r
    Y=subset(data[i_test,], is.na(data[i_test,])==FALSE)
    m=length(Y)
    form10[i_test]=Y%*%Y-(sum(Y)^2)/m
    if (form10[i_test]==0) form10[i_test]=NA
  }

  if(!is.na(min(form10)))
    syy.0 = min(form10)
  else
    syy.0 = 0.0001

  m10.na.pos = which(is.na(output$m10) | is.infinite(output$m10)) 
  if(length(m10.na.pos) != 0)
  {
    for(i_na in m10.na.pos) {
      model=lm(data[i_na,]~time[i_na,])  
      r=resid(model)
      RSS=r%*%r
      Y=subset(data[i_na,], is.na(data[i_na,])==FALSE)
      m=length(Y)
      SSREG=SYY[1]-RSS[1]
    
      output$m10[i_na]=SSREG/syy.0
      
    }  
  }


  # slope 0 
  temp.data = output$m9[(output$m9 != 0)]
  abs.temp.data = abs(temp.data)
  if(length(temp.data) == 0) # if this is the case, all your data has a slope of 0.
    slope.0 = 0.0001
  else
    slope.0 = temp.data[which.min(abs.temp.data)]/100
  
  m17.na.pos = which(is.na(output$m17) | is.infinite(output$m17)) 
  if(length(m17.na.pos) != 0)
    output$m17[m17.na.pos]= output$m15[m17.na.pos]/slope.0

  m18.na.pos = which(is.na(output$m18) | is.infinite(output$m18)) 
  if(length(m18.na.pos) != 0)
    output$m18[m18.na.pos]= output$m12[m18.na.pos]/slope.0


  # mean absolute differences 0
  temp.data = output$m14[(output$m14 != 0)]
  abs.temp.data = abs(temp.data)
  if(length(temp.data) == 0) # if this is the case, all your data has a mean abs. first difference of 0.
    mean.abs.0 = 0.0001
  else
    mean.abs.0 = temp.data[which.min(abs.temp.data)]/100

  m23.na.pos = which(is.na(output$m23) | is.infinite(output$m23)) 
  if(length(m23.na.pos) != 0)
    output$m23[m23.na.pos]= output$m21[m23.na.pos]/mean.abs.0

  m24.na.pos = which(is.na(output$m24) | is.infinite(output$m24)) 
  if(length(m24.na.pos) != 0){
    abs.na.data = abs(SD[m24.na.pos,])
    if(length(m24.na.pos) > 1)
      output$m24[m24.na.pos] = apply(abs.na.data, 1, mean, na.rm = TRUE) / mean.abs.0
    else
      output$m24[m24.na.pos] = mean(abs.na.data, na.rm = TRUE) / mean.abs.0
  }

  # Check correlations
  check.correlation(output[,-1], verbose)

  trajMeasures = structure(list(measurments = output, data = cbind(IDvector, clean.data), 
                                time = cbind(IDvector, clean.time)), class = "trajMeasures")
  
  return(trajMeasures)

}
