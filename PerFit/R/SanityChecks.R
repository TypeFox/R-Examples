# Sanity.dma(): Sanity check - Data matrix adequacy ----
Sanity.dma <- function(matrix, N, I)
{
  if (!is.numeric(matrix) | (sum(matrix == 0 | matrix == 1, na.rm=TRUE) != (N*I - sum(is.na(matrix)))))
  {
    stop('The data matrix is not numeric with 0/1 entries only. Aborted.')
  }
}

# Sanity.dma.poly(): Sanity check - Data matrix adequacy (polytomous) ----
Sanity.dma.poly <- function(matrix, N, I, M)
{
  if (!is.numeric(matrix) | 
        (sum(matrix %in% (0:M), na.rm=TRUE) != (N*I - sum(is.na(matrix)))) | 
        sum(c(min(matrix, na.rm=TRUE), max(matrix, na.rm=TRUE)) != c(0, M), na.rm=TRUE))
  {
    stop('The data matrix is not numeric with entries {0, 1, ..., Ncat-1} only. Aborted.')
  }
}

# sanity.prv(): Sanity check - Perfect response vectors ----
Sanity.prv <- function(matrix, N, I)
{
  NC         <- rowSums(matrix, na.rm = TRUE)
  NC.withNAs <- NC + rowSums(is.na(matrix))
  all.0s     <- vector(mode = "numeric", length = 0)
  all.1s     <- vector(mode = "numeric", length = 0)
  perfect.vectors <- NULL
  if (min(NC)==0 | max(NC.withNAs)==I)
  {
    perfect.vectors <- noquote("Not all item response vectors were included in the analysis (all-0s and/or all-1s patterns removed).")
    all.0s  <- which(NC == 0)
    all.1s  <- which((NC > 0) & (NC.withNAs == I))
    NC      <- NC[-c(all.0s, all.1s)]
  } else
  {
    perfect.vectors <- noquote("All item response vectors were included in the analysis.")
  }
  # Data matrix without perfect vectors:
  if (length(c(all.0s, all.1s)) == 0)
  {
    matrix.red <- matrix
  } else
  {
    matrix.red <- matrix[(1:N)[-c(all.0s, all.1s)],]
  }
  return(list(perfect.vectors=perfect.vectors, all.0s=all.0s, all.1s=all.1s, NC=NC, matrix.red=matrix.red))
}

# Sanity.IPa(): Sanity check - IP matrix adequacy ----
Sanity.IPa <- function(ip, I)
{
  if (!is.numeric(ip) | dim(ip)[1] != I | dim(ip)[2] != 3)
  {
    stop('The item parameters matrix "IP" is not numeric with dimension I x 3 
           (I = number of items; columns = discrimination, difficulty, pseudo-guessing). 
           Aborted.')
  }
}

# Sanity.IPa.poly(): Sanity check - IP matrix adequacy (polytomous) ----
Sanity.IPa.poly <- function(IP, I, Ncat)
{
  if (!is.numeric(IP) | dim(IP)[1] != I | dim(IP)[2] != Ncat)
  {
    stop('The item parameters matrix "IP" is not numeric with dimension I x Ncat 
           (I = number of items; 
            first (Ncat-1) columns = thresholds [GRM] or difficulties [PCM, GPCM]; 
            column Ncat-th = slopes). 
           Aborted.')
  }
}

# Sanity.ABa(): Sanity check - Ability matrix adequacy ----
Sanity.ABa <- function(Ability, N)
{
  if (!is.numeric(Ability) | length(Ability) != N)
  {
    stop('The person parameters vector "Ability" is not numeric with length N 
           (N = number of respondents). 
           Aborted.')
  }
}

# Sanity.IPm(): Sanity check - IP model ----
Sanity.IPm <- function(model)
{
  if (!(model %in% c("1PL", "2PL", "3PL")))
  {
    stop('Parameter "model" can only be "1PL", "2PL", or "3PL". Aborted.')
  }
}

# Sanity.IPm.poly(): Sanity check - IP model (polytomous) ----
Sanity.IPm.poly <- function(model)
{
  if (!(model %in% c("PCM", "GPCM", "GRM")))
  {
    stop('Parameter "model" can only be "PCM", "GPCM", or "GRM". Aborted.')
  }
}

# Sanity.Abm(): Sanity check - Ability method ----
Sanity.Abm <- function(method)
{
  if (!(method %in% c("ML", "BM", "WL")))
  {
    stop('Parameter "method" can only be "ML", "BM", or "WL". Aborted.')
  }
}

# Sanity.Abm.poly(): Sanity check - Ability method (polytomous) ----
Sanity.Abm.poly <- function(method)
{
  if (!(method %in% c("EB", "EAP", "MI")))
  {
    stop('Parameter "method" can only be "EB", "EAP", or "MI". Aborted.')
  }
}

# Sanity.cls(): Sanity check - Class PerFit ----
Sanity.cls <- function(x)
{
  if (class(x) != "PerFit")
  {
    stop('Object "x" is not of class PerFit. Aborted.')
  }
}

# Sanity.clsPO(): Sanity check - Class PerFit.object ----
Sanity.clsPO <- function(x)
{
  if (class(x) != "PerFit.cutoff")
  {
    stop('Object "cutoff.obj" is not of class PerFit.cutoff. Aborted.')
  }
}