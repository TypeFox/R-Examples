# demo.R 
# Part of the SEHmodel package.
#
# Copyright (C) 2015        Melen Leclerc <melen.leclerc@rennes.inra.fr>
#                           Jean-Francois Rey <jean-francois.rey@paca.inra.fr>
#                           Samuel Soubeyrand <Samuel.Soubeyrand@avignon.inra.fr>
#                           Emily Walker <emily.walker@avignon.inra.fr>
#                           INRA - BioSP Site Agroparc - 84914 Avignon Cedex 9
#
# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
#

#' @description Demonstration of SEHmodel package on Bt maize pollen (Genetically Modified crop) on non target Lepidoptera larvae.
#' In this example, 40\% of fields are sources (80 fields) and 100 individuals are exposed to toxic pollen.
#' @title SEHmodel package pollen demonstration
#' @name demo.pollen.run
#' @param nb_fields number of fields (sources and neutrals) in the landscape (default 200)
#' @param max_size landscape size (in meter) (default 5000)
#' @param raster_size raster size (default 2^10)
#' @param nb_ind number of individuals to simulate (default 100)
#' @aliases demo.pollen.run SEHmodel-demo
#' @usage demo.pollen.run(nb_fields=200,max_size=5000,raster_size=2^10,nb_ind=100)
#' @keywords demo.pollen.run demo
#' @export
demo.pollen.run<-function(nb_fields=200,max_size=5000,raster_size=2^10,nb_ind=100) {
  
  # Number of day of simulation
  time.min = 1  #first day
  time.max = 60  #last day
  
  # simulate a landscape of size 5000x5000 meters with 200 fields and 0.4% of toxic fields
  l<-simulateInitialPartition(nb_fields,0.4,10,0,max_size,0,max_size)

  # simulate Margins fields on the prefious landscape
  l2<-simulateThickMargins(l)
  
  # individuals parameters
  dob=sample(time.min:time.max,size=nb_ind,replace=T)  # sample date fo birth for each indivuduals (12 is the pollen issued period)
  life_expectancies=rep(20,nb_ind)  # life duration for each individuals
  toxic_gap=rep(800000,nb_ind)  # intern toxic threshold
  
  # create class individuals with parameters above
  ind<-simulateIndividuals(l2,n=nb_ind,mintime = time.min,maxtime = time.max, dob=dob,life_duration=life_expectancies,toxic_threshold=toxic_gap)
  plot(l2,ind)  # plot the landscape and individuals distribution
  
  # create toxic pollen issued period for each toxic fields
  list_sources<-getSPSources(l2) # get landscape sources fields (toxic)
  pollen.emis<-create.pollen.sources(0.4*nb_fields,time.max)  # 0.4% of 200 fields with time.max numbers of days
  pollen.emis<-data.frame(t=pollen.emis,row.names=row.names(list_sources))
  
  #precipitation on 60 days
  precip=c(0.0,0.0,0.0,0.0,6.5,1.5,1.5,1.5,0.0,0.0,0.0,0.5,4.0,0.0,0.0,4.0,5.0,3.0,1.5, 0.0,0.0,0.0,3.0,5.5,0.0,38.5,7.5,2.0,0.0,0.0,0.0,0.0,0.0,0.0,0.0,3.0,0.5,12.0,0.0,0.0,0.0,0.0,2.5,0.0,0.0,0.0,0.0,0.0,0.0,16.5,0.5,4.5,5.0,0.0,0.0,0.0,0.0,0.0,0.0,3.5)
  
  # Simulate toxic dispertion over the landscape
  td<-toxicIntensity(objectL=l2,toxic_emission=pollen.emis,mintime=time.min,maxtime=time.max,size_raster=raster_size,alpha=list(minalpha=0.1,maxalpha=0.95,covariate_threshold=30,simulate=F,covariate=precip))
  
  # Uncomment to active precipitation simulation
  # td<-toxicIntensity(objectL=l2,toxic_emission=pollen.emis,mintime=time.min,maxtime=time.max,size_raster=raster_size,alpha=list(minalpha=0.1,maxalpha=0.95,covariate_threshold=30,simulate=T,covariate=NULL))
  
  # Simulate toxic exposition on individuals
  ind2<-ecoToxic(l2,ind,td,time.min,time.max)
  
  # plot toxic exposition and assimilation over time for the first individu
  plot(x=l2,y=ind2,objectT=td,numind=1)
  
  # plot individuals info at the end of the simulation
  plot(l2,ind2,time=time.max,objectT=td)
  
  return(list(landscape_demo=l2,individuals_demo=ind2,toxicIntensity_demo=td))
}