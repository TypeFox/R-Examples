calculatedABCanalysis=function(Data){
# Vlist= calculatedABCanalysis(Data)
# Vlist$Aind
# Vlist$Bind
# Vlist$Cind
# Vlist$ABlimit
# Vlist$BClimit 
# Berechnung der ABC analyse ohne plots und sonstige extras
# berechnet ueber ABCanalysis
#
# INPUT
# Data[1:n]         data set,  it is cleaned using  CleanedData = ABCcleanData(Data,RemoveSmallYields)
#                    before results are calculated
# 
# OUTPUT
# Aind,Bind,Cind indices such that:
# Data[Aind]  ===  set A  the "critical few
# Data[Bind]  ===  set B
# Data[Cind]  ===  set C  the "trivial many
# ABlimit            the limit between sets A and B:  [SetA,Aind] = find(Data <=ABlimit );
# BClimit            the limit between sets B and C:  [SetC,Cind] = find(Data > BClimit );
# author: MT 07/2015

abcres = ABCanalysis(Data=ABCcleanData(Data)$CleanedData)

return(list(Aind=abcres$Aind,Bind=abcres$Bind,Cind=abcres$Cind,ABlimit=abcres$smallestAData,BClimit=abcres$smallestBData))
}