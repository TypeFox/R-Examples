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

cs <- function(halfsib, mapPath, separator = " ")
{
  debut <- proc.time()
  if (!is.list(halfsib)) 
    stop("Halfsib should be a list of matrices")
  if(is.data.frame(mapPath))
    map <- mapPath
  else
    map <- as.matrix(read.table(mapPath, sep = separator, header = T,check.names=F))
  if (ncol(map) != 3)
  {
    stop("The map file must contain 3 columns and separated with space")
  }
  print(paste("There are ", length(unique(map[, 2])), " chromosome/s"))
  if (all(colnames(map) != c("Name", "Chr", "Position")))
  {
    stop("The map file first line must contain a header with Name Chr Position, check the header and separator")
  }
  if(length(which(duplicated(map[,1])))>0)
  {
    stop("There are duplicated SNP names in the input file")
  }
  ## map <- map[mixedorder(map[, 3]), ]
  ## map <- map[mixedorder(map[, 2]), ]
  map[, 3] <- as.integer(map[, 3])
  map[, 2] <- as.integer(map[, 2])
  map <- map[order(map[, 3]), ]
  map <- map[order(map[, 2]), ]
  map <- map[!is.na(map[,2]),]
  header <- as.matrix(colnames(halfsib[[1]]))
  selectedHeader <- map[map[,1]%in%header,]
  selectedHeader <- selectedHeader[order(selectedHeader$Position), ]
  selectedHeader <- selectedHeader[order(selectedHeader$Chr), ]
  chrNames <- aggregate(selectedHeader[,1], by = list(selectedHeader[,2]) , function(x) x, simplify = FALSE)
  chrNames <- chrNames[order(as.data.frame(chrNames$Group.1)[, 1]), ]
  chr <- numeric(length(chrNames$Group.1))
  print(paste("Chromosome", length(chr)))
  for (i in 1:length(chr))
  {
    chr[i] <- nrow(as.data.frame(chrNames$x[i]))
  }
  interval <- c(0, cumsum(chr))
  chrNO <- length(interval) - 1
  allhs <- list(chrNO * length(halfsib))
  n <- 1
  i=1
  j=1
  for (i in 1:length(halfsib))
  {
    for (j in 1:chrNO)
    {
      print(paste(n, chrNO * length(halfsib)))
      allhs[[n]] <- halfsib[[i]][, as.character(chrNames$x[[j]])]
      names(allhs)[n] <- paste(names(halfsib[i]), gsub(" ","",chrNames$Group.1[j], fixed=TRUE), sep = "_")
      n <- n + 1
    }
  }
  gc()
  print(proc.time()-debut)
  allhs
} 