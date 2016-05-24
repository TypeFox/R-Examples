#' Generic Function to Export Output to Files
#' @param x object to export to a file
#' @param filename name of the file to which the output should be exported
#' @param ... further arguments for the method
#' @author Tobias Verbeke
#' @export
export <- function(x, filename, ...){
  UseMethod("export")
}

#' Export the summary output for an mpm object to a text file
#' Output the mpm summary to a tab-demimited file for processing by other
#' programs (Excel, Spotfire...)  If the filename is empty, return the data
#' instead of writing to file (useful for web services).
#' 
#' Polar (spherical) coordinates are added if the \code{summary.spm} object
#' contains 2 (3) dimensions.
#' 
#' @param x object of class \code{summary.mpm} as produced by the function of
#'   the same name
#' @param filename prefix used to name the output file following <filename>_xyz.txt
#' @param ... further arguments; currently none are used
#' @return the output is returned invisibly
#' @author Rudi Verbeeck, Tobias Verbeke
#' @seealso \code{\link{summary.mpm}}
#' @method export summary.mpm
#' @S3method export summary.mpm
#' @export
#' @keywords manip
export.summary.mpm <- function(
    x, # summary.mpm object
    filename = "",
    ...)
# Output the mpm summary to a tab-demimited file for processing by other programs (Excel, Spotfire...)
# If the filename is empty, return the data instead of writing to file (useful for web services).
# In this case call the function as (X=dataset to analyse, N=number of required dimensions):
#     mpmResult <- dump.summary.mpm(summary.mpm(mpm(X, row.weight="mean",col.weight="mean"), maxdim=N))
{
  if (!inherits(x, "summary.mpm"))
    stop("Use only with 'summary.mpm' objects.")
  
  # position of factor columns. They are labelled as Prf1, Prf2, Prf3 ... or Pcf1, Pcf2, Pcf3 ... for columns
  FdimsR <- grep("Prf*", names(x$Rows))
  FdimsC <- grep("Pcf*", names(x$Columns))
  if (length(FdimsR) != length(FdimsC))
    stop("Number of principal row factors not equal to number of principal column factors.\nCannot combine into a single dataset.")
  
  # Principal factor columns in x$Rows and x$Columns should have the same name
  # to be combined into a single dataset
  names(x$Rows)[FdimsR] <- paste("PF", 1:length(FdimsR), sep="")
  names(x$Columns)[FdimsC] <- paste("PF", 1:length(FdimsC), sep="")
  # Same is true for the column called "RowWeight" in x$Rows and "ColWeight" in x$Columns
  names(x$Rows)[grep("RowWeight", names(x$Rows))] <- "Weight"
  names(x$Columns)[grep("ColWeight", names(x$Columns))] <- "Weight"
  
  #
  # Add polar/spherical coordinates if summary.mpm structure contains 2 or 3 dimensions
  #
  if (length(FdimsR) == 2)
  {
    # Add columns in 2D: Radius, Angle, Type
    # Type can be R (row in the input matrix) or C (column in the input matrix)
    Rpolar <- complex(real = x$Rows[,FdimsR[1]], imag = x$Rows[,FdimsR[2]])
    Cpolar <- complex(real = x$Columns[,FdimsC[1]], imag = x$Columns[,FdimsC[2]])
    
    r <- rbind(
        cbind(x$Rows,
            Area = sqrt(x$Rows$Weight / max(x$Rows$Weight)),
            Radius = Mod(Rpolar), Angle = Arg(Rpolar)*180/pi,
            Type="R"),
        cbind(x$Columns,
            Area = sqrt(x$Columns$Weight / max(x$Columns$Weight)),
            Radius = Mod(Cpolar), Angle = Arg(Cpolar)*180/pi,
            Type="C"))
  }
  else if (length(FdimsR) == 3)
    # Add columns in 3D: Radius, Azimuth, Elevation, Type
    r <- rbind(
        cbind(x$Rows,
            Volume = (x$Rows$Weight / max(x$Rows$Weight))^(1/3),
            Radius = sqrt(x$Rows[,FdimsR[1]]^2 + x$Rows[,FdimsR[2]]^2 + x$Rows[,FdimsR[3]]^2),
            Azimuth = atan2(x$Rows[,FdimsR[2]], x$Rows[,FdimsR[1]]) * 180/pi,
            Elevation = atan2(x$Rows[,FdimsR[3]], sqrt(x$Rows[,FdimsR[1]]^2 + x$Rows[,FdimsR[2]]^2)) * 180/pi,
            Type="R"),
        cbind(x$Columns,
            Volume = (x$Columns$Weight / max(x$Columns$Weight))^(1/3),
            Radius = sqrt(x$Columns[,FdimsC[1]]^2 + x$Columns[,FdimsC[2]]^2 + x$Columns[,FdimsC[3]]^2),
            Azimuth = atan2(x$Columns[,FdimsC[2]], x$Columns[,FdimsC[1]]) * 180/pi,
            Elevation = atan2(x$Columns[,FdimsC[3]], sqrt(x$Columns[,FdimsC[1]]^2 + x$Columns[,FdimsC[2]]^2)) * 180/pi,
            Type="C"))
  else # Add columns as additional rows to the dataset
  {
    rRadius = x$Rows[,FdimsR[1]]^2
    cRadius = x$Columns[,FdimsC[1]]^2
    for (i in 2:length(FdimsR)) {
      rRadius = rRadius + x$Rows[,FdimsR[i]]^2
      cRadius = cRadius + x$Columns[,FdimsC[i]]^2
    }
    rRadius = sqrt(rRadius)
    cRadius = sqrt(cRadius)
    r <- rbind(
        cbind(x$Rows,
            Volume = (x$Rows$Weight / max(x$Rows$Weight))^(1/length(FdimsR)),
            Radius = rRadius,
            Type="R"),
        cbind(x$Columns,
            Volume = (x$Columns$Weight / max(x$Columns$Weight))^(1/length(FdimsC)),
            Radius = cRadius,
            Type="C"))
  }
  
  #
  # Add row names explicitely. Otherwise they will not be returned to the calling web service.
  #
  r <- cbind(ID = row.names(r), r)
  
  #
  # Write to output file
  #
  if (filename != "") {
    # SO, 20071129 : changed extension to "_xyz" : compatibility with MAP
    # write.table(r, file=paste(filename, "_spm.txt", sep=""), sep="\t", row.names=F, col.names=T)
    write.table(r, file=paste(filename, "_xyz.txt", sep=""), sep="\t", row.names=FALSE, col.names=TRUE)
    # this is how we wrote to file before we explicitely added the row names:
    # write.table(r, file=paste(filename, "_spm.txt", sep=""), sep="\t", col.names=NA)
  }
  
  invisible(r)
}

