# SWRC.R
#
# This file contains a script which plots the measurement sites on a GoogleMap support
#
#
# author: Emanuele Cordano on 08-11-2012

#This program is free software: you can redistribute it and/or modify
#it under the terms of the GNU General Public License as published by
#the Free Software Foundation, either version 3 of the License, or
#(at your option) any later version.
#
#This program is distributed in the hope that it will be useful,
#but WITHOUT ANY WARRANTY; without even the implied warranty of
#MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#GNU General Public License for more details.
#
#You should have received a copy of the GNU General Public License
#along with this program.  If not, see <http://www.gnu.org/licenses/>.

###############################################################################





rm(list=ls())
library(soilwater)

theta_sat <- 0.4
theta_res <- 0.0
alpha <- 0.004# 1/millimeters
n <- 1.8
m <- 1-1/n

psi <- -(0:10000) # millimiters

theta <- swc(psi=psi,alpha=alpha,n=n,m=m,theta_sat=theta_sat,theta_res=theta_res) 
	
plot(theta,psi)