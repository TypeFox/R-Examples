# Copyright (C) 2014 Mohammad H. Ferdosi
#
# HSPhase is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# HSPhase program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http:#www.gnu.org/licenses/>.

.rpohhsphase <- function(genotypeMatrix, oh, forwardVectorSize = 30, excludeFP = TRUE, nsap = 3, maxRec = 15)
{
	
	cat("id group \n",file="temp.txt")
	rhsr_rc <- function(genotypeMatrix ,oh)
	{
		d <- as.dist(.fastdist(oh)) 
		if(length(d)>2)
		{
			fit <- hclust(d, method = "ward")
			groups <- cutree(fit, k = 2)
			a <- which(groups==1)
			b <- which(groups==2)
			# print(length(a))
			if(length(a)>3)
			{
				rec_1 <- recombinations(bmh(genotypeMatrix[a,],excludeFP = excludeFP,nsap = nsap,forwardVectorSize = forwardVectorSize))
				rec_1 <- max(rec_1)
				#print(paste("rec_1",rec_1))
			}
			else
			{
				rec_1 <- 9
			}
			print(length(b))
			if(length(b)>3)
			{
				rec_2 <- recombinations(bmh(genotypeMatrix[b,],excludeFP = excludeFP,nsap = nsap,forwardVectorSize = forwardVectorSize))
				rec_2 <- max(rec_2)
				#print(paste("rec_1",rec_2))
			}
			else
			{
				rec_2 <- 9
			}
			if(rec_1>maxRec)
			{
				rhsr_rc(genotypeMatrix[a,],oh[a,a])
			}
			else
			{		
				write.table(data.frame(names(a),round(abs(rnorm(1)*10^5))),"temp.txt",append = TRUE,col.names = FALSE,row.names = FALSE)
			}
			if(rec_2>maxRec)
			{
				rhsr_rc(genotypeMatrix[b,],oh[b,b])
			}
			else
			{				
				write.table(data.frame(names(b),round(abs(rnorm(1)*10^6))),"temp.txt",append = TRUE,col.names = FALSE,row.names = FALSE)
			}
		}
		else
		{
			if(!is.integer(oh))
			{
				write.table(data.frame(rownames(oh),round(abs(rnorm(1)*10^6))),"temp.txt",append = TRUE,col.names = FALSE,row.names = FALSE)			
			}
			print(paste(class(oh),nrow(oh)))
		
		}
	}
	
	
	result <- rhsr_rc(genotypeMatrix, oh)
	
	result <-read.table("temp.txt",header = TRUE)
	file.remove("temp.txt")
	result
}
