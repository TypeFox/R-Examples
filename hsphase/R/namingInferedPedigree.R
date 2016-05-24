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

pedigreeNaming <- function(inferredPedigree, realPedigree)
{
	inferredPedigree <- inferredPedigree[inferredPedigree[,1]%in%realPedigree[,1],]
	realPedigree <- realPedigree[realPedigree[,1]%in%inferredPedigree[,1],]
	realPedigree <- realPedigree[, 1:2]
	
	colnames(realPedigree) <- colnames(inferredPedigree) <- c("id","group")
	realPedigree$group <- as.character(realPedigree$group)
	inferredPedigree$group <- as.character(inferredPedigree$group)
	
	mergedPartition <- merge(realPedigree,inferredPedigree,by="id")
	mergedPartition <- mergedPartition[order(mergedPartition[,2]),]
	
	
	simfun<-function(x) {x}
	opedigree <- aggregate(mergedPartition,by = list(mergedPartition$group.y),function(x)x)
	
	
	index <- tapply(mergedPartition[,2],mergedPartition[,3],function(x)sort(table(as.character(x)),decreasing = T)[1],simplify = F)
	frqReal <- table(unlist(lapply(index,names)))
	morethan1 <- names(frqReal[frqReal>1])
	for(i in 1:length(index))
	{
		if(all(morethan1!=names(index[[i]]) ))
			inferredPedigree[inferredPedigree[,2] == names(index[i]),2] <- names(index[[i]]) 
		else
		{
			ind <- which(unlist(lapply(index,names)) == names(index[[i]]))
			sireName <- unlist(strsplit(names(sort(unlist(index[ind]),decreasing = T)[1]),split = "[.].*"))
			inferredPedigree[inferredPedigree[,2] == sireName,2] <- names(index[[i]]) 
		}
		
	}
	inferredPedigree
} 