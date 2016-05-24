'logicopt' <- function(in_tt=NULL,
                       n_in=0,
                       n_out=0,
                       find_dc=FALSE,
                       input_sizes=NULL,
                       exact_cover=TRUE,
                       esp_file="",
                       mode="espresso")
########################################################################
########################################################################
#' @title  Truth Table Logic Optimization
#
#' @description This function provides various options to optimize and
#' analyze an input truth table that represents a sum of Boolean or 
#' multi-valued input product terms. 
#' This function leverages the powerful logic minimization
#' algorithms from Espresso.  These algorithms are a standard for optimizing
#' large functions in digital logic synthesis and have been modified to
#' handle general logic minimization problems to serve the R community.
#' The input truth table is an R data frame or alteratively an 
#' Espresso format file.  See examples section for more details.  
# 
#' @return The logicopt function returns a list of two items: an 
#' truth table and a vector.  The vector represents the size and 
#' number of solutions in the output truth table.  For example, a vector [10] 
#' means there is a single solution in the truth table which has 10 rows.  
#' A vector [5 8 2] means there are three solutions in the truth table 
#' with 5, 8, and 2 rows respectively.
#
#' @param in_tt 
#' An R data frame table representing a sum of product terms (PTs) truth table.  
#' The PTs have one or more inputs with a positive integer 
#' value or a "-" which means the input is not specified. The outputs are 
#' Boolean and have possible values 1, 0, "-", or "~" and specify that the
#' corresponding PT is part of the ON set, the OFF set, the DC (don't care) set, 
#' or is unspecified (not a part of any set) respectively.  PTs should not
#' be in both the ON set and OFF set.  When logicopt optimizes, it attempts to
#' find the fewest number of PTs in the ON set and uses PTs in the
#' DC set to further reduce the solution. 
#'
#' @param n_in
#' Integer number of input columns in the truth table.  Inputs 
#' must come first in the truth table. 
#
#' @param  n_out
#' Integer number of output columns in the truth table.  
#' Outputs must come after the n_in inputs of the truth table.
#' The number of columns in the in_tt must be n_in + n_out.
#
#' @param  find_dc 
#' FALSE (default) means any unspecified input conditions are 
#' added to the OFF set for that output.  TRUE means that any unspecified input 
#' conditions are added to DC set for that output. The DC set 
#' is used to further optimize the ON set. Don't cares can also be explicitly
#' defined by using "-" for the output in the input truth table. 
#
#' @param input_sizes 
#' Integer vector which represents the number of possible values for each 
#' input.  Default is NULL which means the size for each input will be 
#' determined automatically by the software by analyzing the input truth table
#' and counting the number of values used.
#' Specifying input_sizes 
#' is important when the input table has unspecified ("-") or unused input
#' values and all possible values are not used.  
#'
#' @param exact_cover
#' Do an exact covering of prime implicant table.  Option applies to QM based
#' algorithms (mode = "qm", "multi-min", and "multi-full").  If FALSE, the covering
#' algoirthm is heuristic and runs faster but may not find an exact solution.
#' If TRUE, algorithm is exact and finds an exact minimum solution.  
#' Default is TRUE.
#'
#' @param esp_file 
#' File name for espresso format file to read and process.  If esp_file is 
#' specified, the input truth table options (in_tt, n_in, n_out, input_sizes,
#' and find_dc) are ignored.  The mode and exact_cover options still apply. 
#' 
#' @param mode
#' A string that specifies the mode to use for optimization:
#' \itemize{
#' \item{"espresso"}{ -- Use the classic espresso algorithm to optimize 
#' the input truth table in_tt or the espresso format table in esp_file. 
#' Returns a single solution of the optimized ON set.  This option should
#' be used for very large truth tables.}
#' \item{"qm"}{ -- Use Quine-McCluskey (QM) algorithm to optimize in_tt 
#' or esp_file and return a single solution of the optimized ON set.  Use
#' with caution for large truth tables.}
#' \item{"primes"}{ -- Return set of prime implicants for in_tt or esp_file. 
#' Three solutions are returned in a single truth table.  These represent the 
#' ESSENTIAL PRIMES, the PARTIALLY REDUNDANT PRIMES, and the TOTALLY REDUNDANT 
#' PRIMES.  Use the \code{\link{print_primes_tt}} function to print the results.}
#' \item{"multi-min"}{ -- Use QM to find the minimum set of non-redundant
#' solutions that cover all the ESSENTIAL PRIMES and the minimal set of
#' PARTIALLLY REDUNDANT PRIMES.  Solutions are ordered 
#' by size.  The number of solutions found is capped at 50.}
#' \item{"multi-full"}{ -- Find additional coverings of prime implicants
#' beyond what is found in multi-min.  A exhaustive covering of all PARTIALLY 
#' REDUNDANT PRIMES is found.  The number of solutions is capped at 50.}
#' \item{"echo"}{ -- Echo the ON, OFF, and DC sets for the in_tt or esp_file 
#' truth table without any optimization.  Note the resulting truth table
#' can be extremely large because it will cover the complete Boolean (or MV)
#' space for all inputs and outputs.  Use with caution.}
#' }
#
#' @examples 
#'
#' ######################### EXAMPLE #1 ###############################
#' # create a truth table with 4 inputs A, B, C, D and 2 outputs X and Y
#' e.ex1 <- data.frame(
#'  A = c(0,0,0,0,0,0,0,0,1,1,1,1,1,1,1,1), 
#'  B = c(0,0,0,0,1,1,1,1,1,1,1,1,0,0,0,0), 
#'  C = c(0,0,1,1,0,0,1,1,0,0,1,1,0,0,1,1), 
#'  D = c(0,1,1,0,0,1,1,0,0,1,1,0,0,1,1,0), 
#'  X = c(1,1,1,1,0,1,1,1,"-",1,1,0,1,1,0,0), 
#'  Y = c("-",1,1,1,1,1,1,1,1,1,1,0,0,0,0,"-"))
#' 
#' # show the unoptimized equations
#' tt2eqn(e.ex1,4,2)
#' 
#' # optimize the truth table
#' tte <- logicopt(e.ex1,4,2)
#' 
#' # show the optimized equations
#' tt2eqn(tte[[1]],4,2)
#' 
#' # generate and print the prime implicants from optimized tte
#' ttp <- logicopt(tte[[1]],4,2,mode="primes")
#' print_primes_tt(ttp,TRUE,4,2)
#' 
#' ######################### EXAMPLE #2 ###############################
#' # get path to an espresso format file 
#' file <- system.file("extdata/espresso/small2.esp", package="LogicOpt")
#'
#' # get the espresso truth table without optimization
#' small2 <- logicopt(esp_file=file,mode="echo")
#'
#' # print the unoptimized equations
#' print_multi_tt(small2,TRUE,4,3)
#' 
#' # optimize with espresso algorithm 
#' small2_opt <- logicopt(small2[[1]],4,3,mode="espresso")
#'
#' # print the optimized equations
#' print_multi_tt(small2_opt,TRUE,4,3)
#'
#' ######################### EXAMPLE #3 ###############################
#' # load up truth table created from a QCA dataset 
#' data(l.represent.1)
#' 
#' # read documentation on how truth table was created 
#' ?l.represent.1
#'
#' # find the set of minimum solutions that cover QCA dataset (mode="multi-min")
#' # treat unspecified input conditions as don't cares (find_dc=TRUE)
#' # find a exact covering (exact_cover=TRUE)
#' lomm <- logicopt(l.represent.1,n_in=5,n_out=1,find_dc=TRUE,
#'   exact_cover=1,mode="multi-min")
#' 
#' # print the solutions in equation format
#' print_multi_tt(lomm,TRUE,5,1,QCA=TRUE)
#' 
#' ######################### EXAMPLE #4 ###############################
#' # optimize a truth table from Genetic Programming  
#' inpath <- system.file("extdata/espresso/robot1_in.esp", package="LogicOpt")
#' robot1 <- logicopt(esp_file=inpath,mode="echo")
#'
#' # unoptimized truth table has 273 rows (256 in ON set and 18 in OFF set)
#' robot1[2]
#'
#' # optimize l.robot1
#' robot1_opt <- logicopt(robot1[[1]],8,3)
#'
#' # optimized results have 13 rows that cover outputs zero, one, and minus
#' robot1_opt[2]
#'
#' # print optimized equations (where each output is 1)
#' print_multi_tt(robot1_opt,TRUE,8,3)
#'
#' ######################### EXAMPLE #5 ###############################
#' # show how to use input_sizes 
#'
#' # get vector of number of unique values for each input 
#' data(l.partybans.1)
#' pb_in_vals <- num_input_values(l.partybans.1,5) 
#' pb_in_vals
#'
#' # optimize with mode=espresso
#' epb <- logicopt(l.partybans.1,5,1,find_dc=TRUE,mode="espresso")
#' epb
#' pb_in_opt_vals <- num_input_values(epb[[1]],5) 
#'
#' # note how some input values have been optimized away and are no longer used! 
#' pb_in_opt_vals
#'
#' # we need original input sizes to process the optimized truth table 
#' qmpb <- logicopt(epb[[1]],5,1,find_dc=FALSE, input_sizes=pb_in_vals,mode="qm")
#' print_multi_tt(epb,TRUE,5,1)
#' print_multi_tt(qmpb,TRUE,5,1)
#'
#' @keywords Espresso Logic Minimization Quine-McKluskey QCA 
#' @export
#' @importFrom utils write.table
#
##############################################################
{

##############################################################
# sanity checks
##############################################################
    if (esp_file == "") {
       if (is.data.frame(in_tt) == FALSE)
          stop('Bad input truth table in_tt.\n')

       if ((n_in + n_out) != ncol(in_tt))
          stop('Number of inputs (',n_in,') + number of outputs (',
           n_out,') is not equal to number of in_tt columns (',
           ncol(in_tt),').\n')

       if (nrow(in_tt) == 0)
          stop('in_tt has no rows.\n')

       file_name <- "esptemp"
       write_esp_file(file_name, input_sizes, in_tt, n_in, n_out, find_dc)
    }
    else {
       if ((find_dc == TRUE) || (n_in > 0) || (n_out > 0) || length(input_sizes) >0)
          warning("When esp_file is used, options find_dc, n_in, n_out, and input_sizes are ignored.")

       if (! file.exists(esp_file))
          stop('File "', esp_file, '" does not exist.\n')
       file_name <- esp_file 
    }
  
##############################################################
# call espresso
##############################################################
    if (exact_cover)
       Sys.setenv(ESP_EXACT_COVER=1)
    else
       Sys.unsetenv("ESP_EXACT_COVER")

    ret <- esp_system(mode, file_name, use_system=FALSE)

    if ((ret != 0) || (! file.exists("esptemp.out"))) {
       unlink("esptemp")
       stop('Call to espresso failed.\n')
    }

##############################################################
# create truth table from espresso results 
##############################################################
    esptt <- readLines("esptemp.out")
 
    # get solution sizes from .p meta commands 
    dot_p <- grep("^.p ", esptt)
    p_str <- esptt[dot_p]
    p_size <- as.numeric(gsub("\\D", "", p_str)) 

    # get input and output names from esp_file .ilb and .ob
    if (esp_file != "") {
       dot_ilb <- as.numeric(grep(".ilb ", esptt)[1])
       ilb_str <- esptt[dot_ilb]
       in_names <- unlist((strsplit(ilb_str,split=" ")))[-1]
       dot_ob <- grep(".ob ", esptt)[1]
       ob_str <- esptt[dot_ob]
       out_names <- unlist((strsplit(ob_str,split=" ")))[-1]
       col_names <- c(in_names,out_names)
    }
    # use names from input truth table
    else {
       col_names <- colnames(in_tt)[seq(ncol(in_tt))]
    }
    
    # ignore espresso meta commands and comments
    esptt <- esptt[!grepl("[.#]", esptt)]
    
    esptt <- t(sapply(esptt, function(x) {
       unlist(strsplit(x, split=" "))
    }))
    
    if ((ncol(esptt) == 0) || (nrow(esptt) == 0)) {
       unlink("esptemp")
       unlink("esptemp.out")
       stop('All terms of input truth table optimized away.\n')
    }

    rownames(esptt) <- seq(nrow(esptt))
    if (length(col_names) == ncol(esptt))
       colnames(esptt) <- col_names 
    else
       warning("No input and output names.  Will use system generated names."
          ,call.=FALSE)

##############################################################
# clean up and return results 
##############################################################
    if (Sys.getenv("LOGOPT_SAVE_TEMP_FILES") == "") {
       unlink("esptemp")
       unlink("esptemp.out")
    }
    
    # return list of truth table and solution sizes 
    return(list(as.data.frame(esptt),p_size))

}

'write_esp_file' <- function(file_name, input_sizes, in_tt, n_in, n_out, find_dc)
##############################################################
# write an espresso format file from a R truth table 
##############################################################
{
    espfile <- file(file_name, "w")
    cat(".mv", n_in+1,"0", file=espfile)

    n_sizes = length(input_sizes)
    if (n_sizes > 0) {
       if (n_sizes != n_in)
          stop('Vector paramater input_sizes has ', n_sizes,  
               ' elements -- it neeeds to have n_in=', 
               n_in, ' elements.\n')
    } 

    for (j in 1:n_in) {
        # if input_sizes specified, use it, otherwise set size must be at least 2 (Boolean) 
        if (n_sizes > 0) size = input_sizes[j] else size = 2 
        nval <- length(unique_int(in_tt[j],n_vals=size))
        if (nval < 2)
           stop('Input ', j, ' must have at least two values.  Specify input_sizes.\n')
        cat(" -", nval, sep="", file=espfile)
     }
    cat(" ", n_out, "\n", file=espfile)

    for (j in 1:n_in) {
       if (n_sizes > 0) size = input_sizes[j] else size = 2
       valj <- unique_int(in_tt[j], n_vals=size)
       cat('.label var=', j-1, valj,'\n', file=espfile)
    }

    if (find_dc)
       cat(".type fr\n", file=espfile) 
    else
       cat(".type fd\n", file=espfile) 
    close(espfile)

    write.table(in_tt, file=file_name,sep=" ",append=TRUE, row.names=FALSE, 
                col.names=FALSE, quote = FALSE)

    espfile <- file(file_name, "a")
    cat(".e\n", file=espfile) 
    close(espfile)
}

'num_input_values' <- function(tt, n_in)
############################################################################
#' @title Find size of input values 
#' @description Find number of unique input values for each input in
#' tt.  
#' @param tt a truth table where first n_in columns are inputs
#' @param n_in number of inputs
#' @return returns a vector of number of unique input values for each input
#' @export
############################################################################
{
    valv <- NULL
    for (j in 1:n_in) {
       valj <- length(unique_int(tt[j], strip_dc=TRUE))
       valv <- c(valv, valj)
    }
    return (valv)
}

'unique_int' <- function(col, strip_dc=TRUE, n_vals=0)
########################################################################
########################################################################
# Return a vector of the unique values in the input column 'col'.
# The col contains possible values from an Espresso optimized 
# truth table.  Don't care values ("-" or "dc") are converted to -1.
# If strip_dc=TRUE, we remove the don't care value from the return
# vector. If n_vals is > 0, the return vector is filled with values 
# (v1 .. vn) to make the vector size n_vals.  
########################################################################
{
    uniq <- unique(col) 
    uint <- as.data.frame(lapply(uniq, function(x) {
           x <- as.character(x)
           x[x %in% c("-", "dc")] <- -1
           return(as.numeric(x))
    }))

    vect <- uint[,1]
    if (strip_dc)
       vect <- vect[vect != -1]

    if (n_vals > 0) {
       vlen <- length(vect) 
       if ((n_vals-vlen) > 0)
          for (i in 1:(n_vals-vlen))
              vect <- c(vect,paste('v',i,sep=''))
    }

    return (sort(vect))
}

'print_multi_tt' <- function(esp_multi, eqn=FALSE, n_in, n_out, max_sol=50, QCA=FALSE)
############################################################################
#' @title Print logicopt() results 
#'
#' @description This function prints the the results from logicopt()
#' in truth table or equation format.
# 
#' @return None
#
#' @param esp_multi 
#' An R data frame table representing a truth table with 1 or more 
#' solutions.
#' 
#' @param eqn 
#' Print in equation format.  Default is FALSE. 
#'
#' @param n_in 
#' Number of inputs in the esp_multi truth table. 
#'
#' @param n_out 
#' Number of outputs in the esp_multi truth table. 
#'
#' @param max_sol 
#' Maximum number of solutions to print.  Default is 50.
#'
#' @param QCA 
#' Attempt to print out in a QCA like format 
#'
#' @examples
#' data(l.partybans.0)
#' tt <- logicopt(l.partybans.0,n_in=5,n_out=1,find_dc=TRUE,mode="multi-full")
#' print_multi_tt(tt,5,1,eqn=TRUE,max_sol=5) 
#'
#' @export
############################################################################
{
   ett <- esp_multi[[1]]
   esize <- esp_multi[[2]]

   n_solutions = length(esize)
   if (max_sol < n_solutions) {
      n_solutions = max_sol
      cat("\nPrinting", n_solutions,"of", length(esize), "solutions.\n")
   }
 
   first <- 1 
   for (i in 1:n_solutions) {
      if (QCA)
         cat("M", i, ": ", sep="")
      else
         cat("SOLUTION", i,"\n")
      last <- esize[i]+first-1
      tt <- ett[first:last,]
      if (eqn)
         cat(tt2eqn(tt,n_in,n_out,QCA),sep="\n")
      else
         print(tt)
      first <- last+1
      if (i == max_sol) return;
   }
}

'tt2eqn' <- function(tt,n_in,n_out,QCA=FALSE)
############################################################################
#' @title Equations from a Truth Table
#'
#' @description This function generates the ON set equations for a truth table.
#' Inputs are uppercase if they are positive and lowercase for negative.
# 
#' @return Vector of equations strings. One for each output. 
#
#' @param tt 
#' R data frame truth table. 
#' 
#' @param n_in 
#' Number of inputs in the tt. 
#'
#' @param n_out 
#' Number of outputs in the tt. 
#'
#' @param QCA 
#' Print in QCA format. 
#'
#' @examples
#' data(l.small)
#' tt <- logicopt(l.small,n_in=4,n_out=3)
#' eqn <- tt2eqn(tt[[1]],4,3)
#'
#' @export
############################################################################
{
   have_mv <- any(apply(tt,2,function(x) x > 1))
   eqn <- ""
   tot_eqn <- NULL
   for (j in 1:n_out) {
      out_nm <- colnames(tt)[n_in+j]
      out_on <- tt[tt[out_nm]=="1",]
      n_rows <- nrow(out_on)
      if (n_rows > 0) {
         if (! QCA) eqn <- paste(out_nm,"= ")
         for (i in 1:n_rows) {
            row <- as.character(t(out_on[i,]))
            first_and <- TRUE
            for (k in 1:n_in) {
                ch <- row[k]
                if (ch != "-") {
                   if (first_and == FALSE)
                    # eqn <- paste(eqn,sep="")
                      eqn <- paste(eqn,"*",sep="")
                   eqn <- paste(eqn,literal(colnames(tt)[k],ch,have_mv),sep="")
                   first_and <- FALSE
               }
            }
            if ((first_and == FALSE) && (i != n_rows))
               eqn <- paste(eqn," + ",sep="")
            #  eqn <- paste(eqn,"+",sep="")
         }
      }
      if (QCA) eqn <- paste(eqn,"<=>",out_nm)
      tot_eqn <- c(tot_eqn,eqn)
      eqn <- ""
   }

   return(tot_eqn)
}

'literal' <- function(name,value,have_mv)
{
   if (have_mv)
      lit <- paste(name,"{",value,"}", sep="")
   else {
      if (value == 0) 
         lit <- tolower(name)
      else if (value == 1)
         lit <- toupper(name)
   }

   return (lit)
}

'print_primes_tt' <- function(primes,eqn=FALSE,n_in,n_out)
############################################################################
#' @title Print the Primes from logicopt(mode="primes") 
#'
#' @description This function prints the results from logicopt(..,,mode="primes")
#' option.
# 
#' @return None
#
#' @param primes 
#' An R data frame table generated by logicopt(...,mode="primes").
#'
#' @param eqn 
#' Print in equation format.  Default is FALSE. 
#'
#' @param n_in 
#' Number of inputs in the primes truth table. 
#'
#' @param n_out 
#' Number of outputs in the primes truth table. 
#'
#'
#' @examples
#' data(l.small)
#' ptt <- logicopt(l.small,n_in=4,n_out=3,find_dc=TRUE,mode="primes")
#' print_primes_tt(ptt,eqn=TRUE,4,3) 
#'
#' @export
############################################################################
{
   ptt <- primes[[1]]
   psize <- primes[[2]]

   if (length(psize) != 3)
      stop("Don't have a primes truth table.\n")

   first <- 1 

   cat("\nESSENTIAL PRIMES\n")
   if (psize[1] > 0) {
      last <- psize[1]+first-1
      tt <- ptt[first:last,]
      if (eqn)
         cat(tt2eqn(tt,n_in,n_out),sep="\n")
      else
         print(tt)
      first <- last+1
   }
   else
      cat("   <empty>\n")

   cat("\nPARTIALLY REDUNDANT PRIMES\n")
   if (psize[2] > 0) {
      last <- psize[2]+first-1
      tt <-ptt[first:last,]
      if (eqn)
         cat(tt2eqn(tt,n_in,n_out),sep="\n")
      else
         print(tt)
      first <- last+1
   }
   else
      cat("   <empty>\n")

   cat("\nTOTALLY REDUNDANT PRIMES\n")
   if (psize[3] > 0) {
      last <- psize[3]+first-1
      tt <- ptt[first:last,]
      if (eqn)
         cat(tt2eqn(tt,n_in,n_out),sep="\n")
      else
         print(tt)
   }
   else
      cat("   <empty>\n")
}


'esp_system' <- function(mode, file_name, use_system=TRUE)
{
    if (use_system) {
      if (mode == "qm") 
         esp_command <- paste('yespresso -srmv -Dqm', file_name,'>esptemp.out')
      else if (mode == "primes") 
         esp_command <- paste('yespresso -srmv -Dprimes', file_name,'>esptemp.out')
      else if (mode == "espresso") 
         esp_command <- paste('yespresso -srmv ', file_name,'>esptemp.out')
      else if (mode == "multi-min")
         esp_command <- paste('yespresso -srmv -Dmulti-min', file_name,'>esptemp.out')
      else if (mode == "multi-full")
         esp_command <- paste('yespresso -srmv -Dmulti-full', file_name,'>esptemp.out')
      else if (mode == "echo")
         esp_command <- paste('yespresso -srmv -Decho', file_name,'>esptemp.out')
      else 
         stop('Mode "', mode, '" not supported.\n', call. = FALSE)

      ret = system(esp_command)
    }
    else {
       ret <- esp_R(mode, file_name)
       # need to fix return 
       ret <- 0
    }
    return(ret)
}


#' @useDynLib LogicOpt esp_main
esp_R <- function(mode, file_name) {
    #libpath <- system.file("libs", package="LogicOpt")
    #lib <- paste0(libpath,"/LogicOpt",.Platform$dynlib.ex)
    #dyn.load(lib)

    if (!is.character(mode)) {
       stop("mode must be of type character.\n")
    }
    if ((mode == "qm") ||
        (mode == "primes") ||
        (mode == "espresso") ||
        (mode == "multi-full") ||
        (mode == "multi-min") ||
        (mode == "echo")) {

        result <- .C(esp_main, esp_mode=mode, esp_file=file_name, PACKAGE="LogicOpt")
        return(result)

    }
    else {
       stop("Unsupported mode\n")
    }

}

