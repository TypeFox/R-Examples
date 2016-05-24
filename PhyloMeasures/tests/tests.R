################################################################################################
##    Copyright (C) 2015,  Constantinos Tsirogiannis and Brody Sandel.  
##
##    Email: analekta@gmail.com and brody.sandel@bios.au.dk
##
##    This file is part of PhyloMeasures.
##
##    PhyloMeasures is free software: you can redistribute it and/or modify
##    it under the terms of the GNU General Public License as published by
##    the Free Software Foundation, either version 3 of the License, or
##    (at your option) any later version.
##
##    PhyloMeasures is distributed in the hope that it will be useful,
##    but WITHOUT ANY WARRANTY; without even the implied warranty of
##    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##    GNU General Public License for more details.
##
##    You should have received a copy of the GNU General Public License
##    along with PhyloMeasures.  If not, see <http://www.gnu.org/licenses/>
################################################################################################

require(ape)
library(PhyloMeasures)

tree.filename = "test.tre"
input.tree = read.tree(tree.filename)
names = input.tree$tip.label

for( k in 0:length(names))
{
  all.samples = t(combn(names,k))

  input.data = matrix(0,nrow = nrow(all.samples),ncol = length(names))

  for( i in 1: nrow(all.samples) )
    for( j in 1: length(names) )
    {
      if(is.element(names[j], all.samples[i,]))
        input.data[i,j] = 1
    }

  colnames(input.data) = names

  ##########################################
  ########## Check PD functions ############ 
  ##########################################

  results.pd = pd.query(input.tree, input.data, is.standardised = FALSE)
  moments.pd = pd.moments(input.tree, c(k))

  expectation.check = 0
  deviation.check = 0

  for(l in 1:length(results.pd))
    expectation.check = expectation.check + results.pd[l]

   expectation.check = expectation.check/length(results.pd)

  for(l in 1:length(results.pd))
  {
    deviation = results.pd[l]-expectation.check
    deviation.check = deviation.check + (deviation*deviation)
  }

  deviation.check = sqrt(deviation.check/length(results.pd))

  if( abs(moments.pd[1] - expectation.check) > 0.01 )
    stop("There is an unexpected discrepancy in the value of the PD expectation.")

  if( abs(moments.pd[2] - deviation.check) > 0.01 )
    stop("There is an unexpected discrepancy in the value of the PD deviation.")


  ##########################################
  ########## Check MPD functions ########### 
  ##########################################

  results.mpd = mpd.query(input.tree, input.data, is.standardised = FALSE)
  moments.mpd = mpd.moments(input.tree, c(k))

  expectation.check = 0
  deviation.check = 0

  for(l in 1:length(results.mpd))
    expectation.check = expectation.check + results.mpd[l]

   expectation.check = expectation.check/length(results.mpd)

  for(l in 1:length(results.mpd))
  {
    deviation = results.mpd[l]-expectation.check
    deviation.check = deviation.check + (deviation*deviation)
  }

  deviation.check = sqrt(deviation.check/length(results.mpd))

  if( abs(moments.mpd[1] - expectation.check) > 0.01 )
    stop("There is an unexpected discrepancy in the value of the MPD expectation.")

  if( abs(moments.mpd[2] - deviation.check) > 0.01 )
    stop("There is an unexpected discrepancy in the value of the MPD deviation.")


  ##########################################
  ######### Check MNTD functions ########### 
  ##########################################

  results.mntd = mntd.query(input.tree, input.data, is.standardised = FALSE)
  moments.mntd = mntd.moments(input.tree, c(k))

  expectation.check = 0
  deviation.check = 0

  for(l in 1:length(results.mntd))
    expectation.check = expectation.check + results.mntd[l]

   expectation.check = expectation.check/length(results.mntd)

  for(l in 1:length(results.mntd))
  {
    deviation = results.mntd[l]-expectation.check
    deviation.check = deviation.check + (deviation*deviation)
  }

  deviation.check = sqrt(deviation.check/length(results.mntd))

  if( abs(moments.mntd[1] - expectation.check) > 0.01 )
    stop("There is an unexpected discrepancy in the value of the MNTD expectation.")

  if( abs(moments.mntd[2] - deviation.check) > 0.01 )
    stop("There is an unexpected discrepancy in the value of the MNTD deviation.")


  ##########################################
  ########## Check CAC functions ########### 
  ##########################################
  my.chi = 0.76
  results.cac = cac.query(input.tree, input.data, my.chi, is.standardised = FALSE)
  moments.cac = cac.moments(input.tree, my.chi, c(k))

  expectation.check = 0
  deviation.check = 0

  for(l in 1:length(results.cac))
    expectation.check = expectation.check + results.cac[l]

   expectation.check = expectation.check/length(results.cac)

  for(l in 1:length(results.cac))
  {
    deviation = results.cac[l]-expectation.check
    deviation.check = deviation.check + (deviation*deviation)
  }

  deviation.check = sqrt(deviation.check/length(results.cac))

  if( abs(moments.cac[1] - expectation.check) > 0.01 )
    stop("There is an unexpected discrepancy in the value of the CAC expectation.")

  if( abs(sqrt(moments.cac[2]) - deviation.check) > 0.01 )
    stop("There is an unexpected discrepancy in the value of the CAC deviation.")


  for( h in 0: length(names) )
  {
    all.samples.b = t(combn(names,h))

    input.data.b = matrix(0,nrow = nrow(all.samples.b),ncol = length(names))

    for( i in 1:nrow(all.samples.b) )
      for( j in 1: length(names) )
      {
        if(is.element(names[j], all.samples.b[i,]))
          input.data.b[i,j] = 1
      }

    colnames(input.data.b) = names

    input.sizes = matrix(0,nrow = 1, ncol = 2) 
    input.sizes[1,1] = k  
    input.sizes[1,2] = h

    ##########################################
    ########## Check CBL functions ########### 
    ##########################################
    results.cbl = cbl.query(input.tree, input.data, matrix.b = input.data.b, is.standardised = FALSE)
    moments.cbl= cbl.moments(input.tree, input.sizes)

    expectation.check = 0
    deviation.check = 0

    for(l in 1:nrow(results.cbl))
      for(r in 1:ncol(results.cbl))
      expectation.check = expectation.check + results.cbl[l,r]

    expectation.check = expectation.check/(nrow(results.cbl)*ncol(results.cbl))

    for(l in 1:nrow(results.cbl))
      for(r in 1:ncol(results.cbl))
      {
        deviation = results.cbl[l,r]-expectation.check
        deviation.check = deviation.check + (deviation*deviation)
      }

    deviation.check = sqrt(deviation.check/(nrow(results.cbl)*ncol(results.cbl)))

    if( abs(moments.cbl[1] - expectation.check) > 0.01 )
      stop("There is an unexpected discrepancy in the value of the CBL expectation.")

    if( abs(moments.cbl[2] - deviation.check) > 0.01 )
      stop("There is an unexpected discrepancy in the value of the CBL deviation.")


    ##########################################
    ########## Check CD functions ########### 
    ##########################################

    results.cd = cd.query(input.tree, input.data, matrix.b = input.data.b, is.standardised = FALSE)
    moments.cd= cd.moments(input.tree, input.sizes)

    expectation.check = 0
    deviation.check = 0

    for(l in 1:nrow(results.cd))
      for(r in 1:ncol(results.cd))
      expectation.check = expectation.check + results.cd[l,r]

    expectation.check = expectation.check/(nrow(results.cd)*ncol(results.cd))

    for(l in 1:nrow(results.cd))
      for(r in 1:ncol(results.cd))
      {
        deviation = results.cd[l,r]-expectation.check
        deviation.check = deviation.check + (deviation*deviation)
      }

    deviation.check = sqrt(deviation.check/(nrow(results.cd)*ncol(results.cd)))

    if( abs(moments.cd[1] - expectation.check) > 0.01 )
      stop("There is an unexpected discrepancy in the value of the CD expectation.")

    if( abs(moments.cd[2] - deviation.check) > 0.01 )
      stop("There is an unexpected discrepancy in the value of the CD deviation.")

  } # for( h in 0: length(names) )

} # for( k in 1:length(names))

cat("\n")
cat("---------- All tests were completed successfully ----------")
cat("\n")
cat("\n")





