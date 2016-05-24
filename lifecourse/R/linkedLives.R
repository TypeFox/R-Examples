

linkedLives = function(currentData,sequenceUnits,lifecourseSequences){
  courseLen = length(lifecourseSequences[1,])
  bbcseq = lifecourseSequences
  hen = length(lifecourseSequences[,1])
  slen = length(sequenceUnits)
  testStats = array(0,c(1000,1))
  totalScore = score = 0    			# calculating the mobility index for the observed sequences
  for(i in 1:length(bbcseq[,1])){
    myseq = bbcseq[i,1:courseLen]
    score = (mobility_index(myseq,sequenceUnits,slen)) 
    totalScore = totalScore + score
  }
  testStatistic = totalScore/(courseLen*hen) # 
  #--------------------------------------- 
  for(j in 1:1000){
    totalScore = score = 0
    Population = lifecourseSequences
        
      thePopulation = currentData[order(currentData[,1]),]
      selectYear = thePopulation[,1]
      U = length(unique(selectYear))
      UY = unique(selectYear)
      for(g in 1:U){
        candidates = which(selectYear%in%UY[g])
        Population[candidates,] = sample(bbcseq[candidates[1],])
      }
          
    for(i in 1:hen){
      myseq = Population[i,] 
      score = (mobility_index(myseq,sequenceUnits,courseLen))
      totalScore = totalScore + score
    }
    testStats[j] = totalScore/(hen*courseLen) # as the lifecourse becomes more fluid, we observe more changepoints and the test statistic is larger.
  }
  LeftTailed_pval = length(which(testStats<=testStatistic))/1000
  RightTailed_pval = length(which(testStats>=testStatistic))/1000
  list(testStats,LeftTailed_pval,RightTailed_pval,testStatistic) 
}
