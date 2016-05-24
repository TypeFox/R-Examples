##############################################################################
#' Performanced of a Sun SPARCcenter 2000 in the SPEC SDM91 benchmark
#'
#' A dataset containing performance data for a Sun SPARCcenter 2000 (16 CPUs)
#'
#' A Sun SPARCcenter 2000 with 16 CPUs was used for the SPEC SDM91 benchmark
#' in October 1994. The benchmark simulates a number of users working on the
#' UNIX server and measures the number of script executions per hour.
#'
#' The data frame contains the following variables:
#' \itemize{
#'   \item \code{load} The number of simulated users (1--216).
#'   \item \code{throughput} The achieved throughput in scripts per hour.
#' }
#'
#' @name specsdm91
#' @docType data
#' @keywords datasets
#' @format A data frame with 7 rows on 2 variables
#' @source Neil J. Gunther. Guerrilla Capacity Planning: A Tactical
#'   Approach to Planning for Highly Scalable Applications and Services.
#'   Springer, Heidelberg, Germany, 1st edition, 2007.
#'   Original dataset from
#'   \url{http://www.spec.org/osg/sdm91/results/results.html}
NULL


##############################################################################
#' Performance of a ray-tracing software on different hardware configurations
#'
#' A dataset containing performance data for a ray-tracing benchmark.
#'
#' The benchmark measured the number of ray-geometry intersections per second.
#' The data was gathered on an SGI Origin 2000 with 64 R12000 processors
#' running at 300 MHz.
#'
#' The data frame contains the following variables:
#' \itemize{
#'   \item \code{processors} The number of CPUs used for the benchmark (1--64).
#'   \item \code{throughput} The number of operations per second.
#' }
#'
#' @name raytracer
#' @docType data
#' @keywords datasets
#' @format A data frame with 11 rows on 2 variables
#' @source Neil J. Gunther. Guerrilla Capacity Planning: A Tactical
#'   Approach to Planning for Highly Scalable Applications and Services.
#'   Springer, Heidelberg, Germany, 1st edition, 2007.
#'   Original dataset from \url{http://sourceforge.net/projects/brlcad/}
NULL


##############################################################################
#' Performance of an Oracle database used for online transaction processing
#'
#' A dataset containing performance data for an Oracle OLTP database measured
#' between 8:00am and 8:00pm on January, 19th 2012. The measurements were
#' recorded for two minute intervals during this time and a timestamp indicates
#' the end of the measurement interval. The performance metrics were taken from
#' the \code{v$sysmetric} family of system performance views.
#'
#' The Oracle database was running on a 4-way server.
#'
#' The data frame contains different types of measurements:
#' \itemize{
#'   \item Variables of the "time" type are expressed in seconds per second.
#'   \item Variables of the "rate" type are expressed in events per second.
#'   \item Variables of the "util" type are expressed as a percentage.
#' }
#'
#' The data frame contains the following variables:
#' \itemize{
#'   \item \code{timestamp} The end of the two minute interval for which the
#'     remaining variables contain the measurements.
#'   \item \code{db_time} The time spent inside the database either working on
#'     a CPU or waiting  (I/O, locks, buffer waits ...). This time is expressed
#'     as seconds per second, so two sessions working for exactly one second
#'     each will contribute a total of two seconds per second of \code{db_time}.
#'     In Oracle this value is also known as \emph{Average Active Sessions}
#'     (AAS).
#'   \item \code{cpu_time} The CPU time used during the interval. This is also
#'     expressed as seconds per second. A 4-way machine has a theoretical
#'     capacity of four CPU seconds per second.
#'   \item \code{call_rate} The number of user calls (logins, parses, or
#'     execute calls) per second.
#'   \item \code{exec_rate} The number of statement executions per second.
#'   \item \code{lio_rate} The number of logical I/Os per second. A logical
#'     I/O is the Oracle term for a cache hit in the database buffer cache.
#'     This metric does not indicate if an additional physical I/O was
#'     necessary to load the buffer from disk.
#'   \item \code{txn_rate} The number of database transactions per second.
#'   \item \code{cpu_util} The CPU utilization of the database server in
#'     percent. This was also measured from within the database.
#' }
#'
#' @name oracledb
#' @docType data
#' @keywords datasets
#' @format A data frame with 360 rows on 8 variables
#'
NULL