# Calculates biotic index scores for individual samples.
# @description Calculates chosen biotic index on individual samples,
# based on vectors of abundances and taxon names passed as arguments.
#
# @param abundances An integer vector representing the abundances of taxa
# in a sample. Absences should be represented as \code{NA} rather than
# zero values but it will process data with zeros.
# @param taxonlist A string vector containing taxon names (family level)
# present in all the samples being analysed.
# @param index A string representing the choice of index value. Defaults to
# BMWP.
# @return A list containing the different components of each score. List
# dimensions depend on the index calculated.

calcscore<-function(abundances, taxonlist, index){

  # recombine taxon list and abundances
  sampledata<-as.data.frame(cbind.data.frame(taxonlist, abundances))

  # find rows where taxa are present (>0 abundance)
  taxapresent<-sampledata[sampledata[,2] > 0,]

  # check for BWMP composite families and remove, except for PSI and WHPT
    if (index!="WHPT" && index!="WHPT_AB" && index!="PSI"){

    # extract composite families from indextable
    composites<-indextable[indextable$Composite!="",]
    composites<-as.data.frame(cbind.data.frame(composites$Taxon, composites$Composite))

    # locate any rows to delete
    rowstodelete<-apply(composites, 1, findcomposites, df=taxapresent)

    # remove NAs from vector
    rowstodelete<-na.omit(rowstodelete)

    # remove double counting rows from taxapresent if there are any
    if (length(rowstodelete)>0){
      taxapresent<-taxapresent[-rowstodelete, ]
    }
  }

  # check that there are any taxa present in the sample, then extract scores
  if (nrow(taxapresent)!=0){
    if (index=="BMWP"){
      samplescores<-extractrows(taxapresent, indextable)
      scorelist<-c(sum(samplescores$BMWP, na.rm=TRUE), sum(!is.na(samplescores$BMWP)), round(sum(samplescores$BMWP, na.rm=TRUE)/sum(!is.na(samplescores$BMWP)),2))
    }
    if (index=="Whalley"){
      samplescores<-extractrows(taxapresent, indextable)
      scorelist<-c(sum(samplescores$Whalley, na.rm=TRUE), sum(!is.na(samplescores$Whalley)), round(sum(samplescores$Whalley, na.rm=TRUE)/sum(!is.na(samplescores$Whalley)),2))
    }
    if (index=="Riffle"){
      samplescores<-extractrows(taxapresent, indextable)
      scorelist<-c(sum(samplescores$Riffle, na.rm=TRUE), sum(!is.na(samplescores$Riffle)), round(sum(samplescores$Riffle, na.rm=TRUE)/sum(!is.na(samplescores$Riffle)),2))
    }
    if (index=="Pool"){
      samplescores<-extractrows(taxapresent, indextable)
      scorelist<-c(sum(samplescores$Pool, na.rm=TRUE), sum(!is.na(samplescores$Pool)), round(sum(samplescores$Pool, na.rm=TRUE)/sum(!is.na(samplescores$Pool)),2))
    }
    if (index=="RiffPool"){
      samplescores<-extractrows(taxapresent, indextable)
      scorelist<-c(sum(samplescores$RiffPool, na.rm=TRUE), sum(!is.na(samplescores$RiffPool)), round(sum(samplescores$RiffPool, na.rm=TRUE)/sum(!is.na(samplescores$RiffPool)),2))
   }
    if (index=="WHPT"){
      samplescores<-extractrows(taxapresent, indextable)
      scorelist<-c(round(sum(samplescores$WHPT, na.rm=TRUE)/sum(!is.na(samplescores$WHPT)),2), sum(!is.na(samplescores$WHPT)))
    }
    if (index=="WHPT_AB"){

      # slice abundance array and extract scores
      WHPTABscores<-logslice(taxapresent)

      # divide array based on abundance class
      WHPTclass1<-WHPTABscores[WHPTABscores$class ==1,]
      WHPTclass2<-WHPTABscores[WHPTABscores$class ==2,]
      WHPTclass3<-WHPTABscores[WHPTABscores$class ==3,]
      WHPTclass4<-WHPTABscores[WHPTABscores$class ==4,]

      # extract scores from relevant columns
      WHPTscore1<-sum(WHPTclass1$WHPT_AB1, na.rm=TRUE)
      WHPTscore2<-sum(WHPTclass2$WHPT_AB2, na.rm=TRUE)
      WHPTscore3<-sum(WHPTclass3$WHPT_AB3, na.rm=TRUE)
      WHPTscore4<-sum(WHPTclass4$WHPT_AB4, na.rm=TRUE)

      # calculate index (ASPT and N-taxa)
      scorelist<-c(round(sum(WHPTscore1, WHPTscore2, WHPTscore3, WHPTscore4)/sum(!is.na(WHPTABscores$WHPT_AB1)),2), sum(!is.na(WHPTABscores$WHPT_AB1)))
    }
    if (index=="AWIC"){
      samplescores<-extractrows(taxapresent, indextable)
      scorelist<-c( round(sum(samplescores$AWIC, na.rm=TRUE) / sum(!is.na(samplescores$AWIC)),2))
      }

    if (index=="PSI"){

      # slice abundance array and extract scores
      PSIscores<-logslice(taxapresent)

      # increment PSI scores based on log abundance class
      PSIscores$PSI<-ifelse(PSIscores$class==2, PSIscores$PSI+1, ifelse(PSIscores$class==3, PSIscores$PSI+2, ifelse(PSIscores$class==4, PSIscores$PSI+3, PSIscores$PSI)))
      PSI_ABscores<-PSIscores[PSIscores$PSIcat=="A" | PSIscores$PSIcat=="B",]
      PSI_AB<-sum(PSI_ABscores$PSI, na.rm=TRUE)
      PSI_All<-sum(PSIscores$PSI, na.rm=TRUE)

      scorelist<-c(round((PSI_AB/PSI_All)*100,2))
    }
    if (index=="LIFE"){
      # LIFE calculation

      # slice abundance array and extract scores
      LIFEscores<-logslice(taxapresent)

      # divide array based on LIFE score
      LIFEabove7<-LIFEscores[LIFEscores$LIFE > 7,]
      LIFE7<-LIFEscores[LIFEscores$LIFE==7,]
      LIFEbelow7<-LIFEscores[LIFEscores$LIFE < 7,]

      # modify LIFE scores based on log abundance class
      LIFEabove7$LIFE<-ifelse(LIFEabove7$class==2, LIFEabove7$LIFE+1, ifelse(LIFEabove7$class==3, LIFEabove7$LIFE+2, ifelse(LIFEabove7$class==4, LIFEabove7$LIFE+3, LIFEabove7$LIFE)))
      LIFEbelow7$LIFE<-ifelse(LIFEbelow7$class==2, LIFEbelow7$LIFE-1, ifelse(LIFEbelow7$class==3, LIFEbelow7$LIFE-2, ifelse(LIFEbelow7$class==4, LIFEbelow7$LIFE-3, LIFEbelow7$LIFE)))

      # recombine sections
      LIFEscores<-rbind(LIFEabove7, LIFE7, LIFEbelow7)

      # calculate LIFE index (sum/n)
      scorelist<-c(round(sum(LIFEscores$LIFE, na.rm=TRUE)/sum(!is.na(LIFEscores$LIFE)),2))
    }
  return(scorelist)
  # this closes the nrows loop
  }
  # this is the function closure
}
