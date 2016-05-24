# Copyright (C) 2013 Mohammad H. Ferdosi
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

.prSimple <- function (oh, snpNooh, intercept = 26.3415, coefficient = 77.3171 ) 
{
	maxsnpnooh <-  (intercept + snpNooh * coefficient) - (15  * snpNooh)
	cat("id group \n", file = "temp.txt")
	rhsr_rc <- function(oh, maxsnpnooh = maxsnpnooh) 
	{
		print("----")		
		## d <- dist(oh, method = "manhattan")
		d <- as.dist(.fastdist(oh))
		if (length(d) > 2) {
			fit <- hclust(d, method = "ward")
			groups <- cutree(fit, k = 2)
			a <- which(groups == 1)
			b <- which(groups == 2)
			
			
			if(length(a)>2)
			{
				subohA <- oh[a,a]
				maxSubohA <- max(subohA[lower.tri(subohA)])
			}
			else
			{
				maxSubohA = 0 
			}
			
			if(length(b)>2)
			{		
				subohB <- oh[b,b]
				maxSubohB <- max(subohB[lower.tri(subohB)])
			}
			else
			{
				maxSubohB = 0 				
			}
			
			if (maxSubohA > maxsnpnooh && length(a)>2) 
			{
				
				rhsr_rc(oh[a, a],maxsnpnooh)
			}
			else {
				
				write.table(data.frame(names(a), round(abs(rnorm(1) * 
														10^5))), "temp.txt", append = TRUE, col.names = FALSE, 
						row.names = FALSE)
			}
			if (maxSubohB > maxsnpnooh  && length(b)>2)
			{
				rhsr_rc(oh[b, b],maxsnpnooh)
			}
			else {
				write.table(data.frame(names(b), round(abs(rnorm(1) * 
														10^6))), "temp.txt", append = TRUE, col.names = FALSE, 
						row.names = FALSE)
			}
		}
		else {
			if (!is.integer(oh)) 
				write.table(data.frame(rownames(oh), round(abs(rnorm(1) * 
														10^6))), "temp.txt", append = TRUE, col.names = FALSE, 
						row.names = FALSE)
		}
	}
	result <- rhsr_rc(oh, maxsnpnooh)
	result <- read.table("temp.txt", header = TRUE)
	file.remove("temp.txt")
	result
}

