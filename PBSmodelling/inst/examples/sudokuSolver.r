#***************************************************************
# Sudok Solver
# Alex Couture-Beil <alex@mofo.ca>
#
# Requires:
#   PBSmodelling package
# 
# Usage:
#   source("sudokuSolver.r")
#   type "solveSudoku()" to solve the puzzle with animated steps (cooler)
#   or "solveSudoku(FALSE)" to solve it without displaying steps (faster)
#
# Details:
#     The aim of the puzzle is to enter the digits 1 through 9 in each cell of a 
#   9x9 grid made up of 3x3 subgrids (called "regions") so that each row, column, 
#   and region contains exactly one instance of each digit. A set of clues, or 
#   "givens", constrain the puzzle such that there is only one way to correctly 
#   fill in the remainder.
#   -- Extract from Wikipedia, Sudoku, http://en.wikipedia.org/wiki/Sudoku 
#
#---------------------------------------------------------------

# Checks if a given sudoku puzzle (does not have to be complete)
# is consistant with the constraints of sudoku
# returns TRUE if consistant, FALSE otherwise

local(envir=.PBSmodEnv,expr={
locale = sys.frame(sys.nframe() - 1) # local environment

is.consistant.master <- function(x)
{
	for(i in 1:9) {
		if (test.dup(x[i,]))
			return(FALSE)
		if (test.dup(x[,i]))
			return(FALSE)
	}
	
	for(i in 0:2)
		for(j in 0:2)
			if (test.dup(c(x[1+3*i,1+3*j], x[2+3*i,1+3*j], x[3+3*i,1+3*j], 
			               x[1+3*i,2+3*j], x[2+3*i,2+3*j], x[3+3*i,2+3*j], 
			               x[1+3*i,3+3*j], x[2+3*i,3+3*j], x[3+3*i,3+3*j])))
				return(FALSE)
	
	if (any(!is.na(x)))
		if (any(x<1, na.rm=TRUE) || any(x>9, na.rm=TRUE))
			stop("warning one value is out of the range of 0-9")
	
	return(TRUE)
}

#returns TRUE if duplicates are found (fail)
#FALSE otherwise
test.dup <- function(x)
{
	last <- 0
	x <- sort(x)
	for(i in x)
		if (i==last)
			return(TRUE)
		else
			last <- i
	return(FALSE)
}

#checks if the state is consistant
#with the assumption that the previous state (where x[i,j]=NA) was consistant
is.consistant <- function(x,i,j)
{
	if (any(x[i,][-j]==x[i,j], na.rm = TRUE))
		return(FALSE)
	if (any(x[,j][-i]==x[i,j], na.rm = TRUE))
		return(FALSE)
	
	ii <- ceiling(i/3)-1
	jj <- ceiling(j/3)-1
	
	
	for(a in (1:3)+3*ii)
		for(b in (1:3)+3*jj)
			if ((a!=i || b!=j) && !is.na(x[a,b]))
				if (x[i,j]==x[a,b])
					return(FALSE)
	
	return(TRUE)
}

#tests for completeness
#returns TRUE if complete, FALSE otherwise
is.complete <- function(x)
{
	return(!any(is.na(x)))
}

#returns the next unassigned variable location in the form c(i,j)
#if the state is complete (ie no unassigned vars)
#return NULL
getNextVar <- function(x)
{
	d <- dim(x)
	for(i in 1:d[1])
		for(j in 1:d[2])
			if (is.na(x[i,j]))
				return(c(i,j))
	return(NULL)
}

#clear out the matrix
clearPuz <- function()
{
	setWinVal(list(s=matrix(NA,9,9)))
}

#restore saved sudoku puzzle (ie before running solveSudoku())
resetPuz <- function()
{
	if (exists("sudokuMat",envir=locale))
		setWinVal(list(s=tget(sudokuMat,tenv=locale)))
}

#function to display message instructing user to run appropriate cmd
savePuz <- function()
{
	sudokuMat <- getWinVal()$s; tput(sudokuMat,tenv=locale)
	
	if (!is.consistant.master(sudokuMat))
		stop("initial sudoku puzzle contains errors")
	
	cat("------------------------------------------------------------\n")
	cat("Due to the limitations of R and tcl/tk in terms of threading\n")
	cat("you must execute the function solveSudoku() from the console.\n")
	cat("If the window invoked the function, the window would not\n")
	cat("be updated until the function has finished.\n\n")
	cat("Type \"tcall(solveSudoku)()\" in the R console prompt below.\n")
	cat("------------------------------------------------------------\n")
}

#Start a depth-first search for the solution
solveSudoku <- function(refreshWin=TRUE, solHalt=TRUE)
{
	sudokuMat <- getWinVal()$s; tput(sudokuMat,tenv=locale)

	if (!is.consistant.master(sudokuMat))
		stop("initial sudoku puzzle contains errors")

	nodes <- list()
	nodes[[1]] <- sudokuMat
	while(1) {
		if (length(nodes)==0) {
			print("Exhausted List.");
			break;
		}
		
		#pop off last added node
		node <- nodes[[length(nodes)]]
		nodes[[length(nodes)]] <- NULL
		
		#get next unassigned variable position
		ind <- getNextVar(node)
		
		#update window with new state
		if (refreshWin)
			setWinVal(list(s=node))
		
		#check for completeness
		if (is.null(ind)) {
			print("hooray - a solution was found.");
			print(node)
			cat(" -------------------------------------------------\n");
			if (solHalt)
				break;
			#pause("Solution found. press <enter> to continue. <esc> to halt. ")
			flush.console()
			
		}
		#expand state otherwise
		else {
			for(i in 9:1) {
				node[ind[1], ind[2]] <- i
				if (is.consistant(node, ind[1], ind[2])) {
					nodes[[length(nodes)+1]] <- node
				}
			}
		}
	}
}

#load one of two sample puzzles
loadPuz <- function(puzzle)
{
	if (missing(puzzle))
		puzzle <- getWinAct()[1]

	if (puzzle=="easy") {
		setWinVal(list(s=structure(c(5, 8, NA, NA, NA, 6, 9, 3, NA, 9, 7, NA, NA, 
		3, NA, 6, NA, NA, NA, NA, NA, NA, 9, 8, 5, NA, 2, NA, 4, NA, 
		NA, NA, 9, 3, NA, 5, NA, 6, NA, NA, NA, NA, NA, 2, NA, 8, NA, 
		3, 4, NA, NA, NA, 6, NA, 6, NA, 4, 7, 5, NA, NA, NA, NA, NA, 
		NA, 8, NA, 6, NA, NA, 4, 1, NA, 2, 9, 8, NA, NA, NA, 5, 6), .Dim = as.integer(c(9, 
		9)))))
	}
	else if (puzzle=="evil") {
		setWinVal(list(s=structure(c(3, 2, NA, 8, NA, 9, NA, NA, NA, NA, NA, NA, 
		2, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, NA, 9, 7, 2, NA, 
		1, NA, 8, NA, NA, 3, NA, NA, 5, NA, NA, NA, NA, NA, 4, NA, NA, 
		8, NA, NA, 6, NA, 5, NA, 9, 7, 4, NA, NA, NA, NA, NA, NA, NA, 
		NA, NA, NA, NA, NA, 8, NA, NA, NA, NA, NA, NA, 7, NA, 5, NA, 
		8, 1), .Dim = as.integer(c(9, 9)), .Dimnames = list(NULL, NULL))))
	}
}

#display help messages
sudokuHelp <- function(section=NULL)
{
	if (is.null(section)) section <- getWinAct()[1]
	
	if (section=="about") {
cat("-------------------------------------------------------------------------------
Sudoku, also known as Number Place or Nanpure, is a logic-based placement 
puzzle. The aim of the puzzle is to enter the digits 1 through 9 in each cell 
of a 9x9 grid made up of 3x3 subgrids (called \"regions\") so that each row, 
column, and region contains exactly one instance of each digit. A set of clues, 
or \"givens\", constrain the puzzle such that there is only one way to 
correctly fill in the remainder.

Completed sudoku puzzles are a type of Latin square, with the additional 
constraint on the contents of individual regions. Leonhard Euler is sometimes 
cited as the source of the puzzle based on his work with Latin squares.

The modern puzzle Sudoku was invented in Indianapolis in 1979 by Howard Garns. 
Garns' puzzles appeared in Dell Magazines, which published them under the title 
\"Number Place\". Sudoku became popular in Japan in 1986, when puzzle publisher 
Nikoli discovered the game in older Dell publications. The puzzles became an 
international hit in 2005.

-- Extract from Wikipedia, Sudoku, http://en.wikipedia.org/wiki/Sudoku 
(Retrieved on October 3, 2006)
-------------------------------------------------------------------------------
")
	}
	if (section=="usage") {
		cat("-------------------------------------------------------------------------------
The sudoku puzzle is represented as a 9x9 matrix, which contains 81 variables. 
The solver attempts to solve this constraint satisfaction problem by using a 
depth-first search technique. Since this problem can be expressed as a graph 
colouring problem, it should be no surprise that this algorithm requires 
exponential time. (given the matrix was an n x n matrix)

The solver saves all the defined constraints, as displayed in the GUI, as the 
initial state of the sudoku puzzle. This state is saved as the first element of 
a list of states.

The solver removes the last saved state from the list and then scans the matrix 
left to right, and top to bottom looking for unassigned variables. Upon finding 
an unassigned variable, the solver calculates all possible values that do not 
conflict with any existing assigned variables. These numbers are then inserted 
into N new states which are appended to the end of the states list.

The solver then reiterates this process of taking the last saved state, and 
replacing it with N new states until the matrix is complete, that is, no 
unassigned variables were found.

Updating the GUI requires a lot of processor power, and significantly slows 
down the solving speed. It is therefore possible to disable the animated 
solving by setting the \"refreshWin\" argument to FALSE.

If you would like to display all possible solutions by exhaustively searching 
the list of possible states, set the \"solHalt\" to FALSE.


Usage: solveSudoku(refreshWin=TRUE, solHalt=TRUE)
-------------------------------------------------------------------------------
")
}
}

#create window and initialize with a sample problem.
if (!require(PBSmodelling, quietly=TRUE)) stop("The PBSmodelling package is required for this example")
createWin("sudokuSolverWin.txt")
loadPuz("easy")

}) # end local scope
