##  RUnit : A unit test framework for the R programming language
##  Copyright (C) 2003-2009  Thomas Koenig, Matthias Burger, Klaus Juenemann
##
##  This program is free software; you can redistribute it and/or modify
##  it under the terms of the GNU General Public License as published by
##  the Free Software Foundation; version 2 of the License.
##
##  This program is distributed in the hope that it will be useful,
##  but WITHOUT ANY WARRANTY; without even the implied warranty of
##  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##  GNU General Public License for more details.
##
##  You should have received a copy of the GNU General Public License
##  along with this program; if not, write to the Free Software
##  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA


cat("\n\nRUnit test cases for 'printJUnitProtocol' function\n\n")


testRUnit.printJUnitProtocol <- function()
{
  ##  copy baseenv() logger
  tmp <- get(".testLogger", envir = RUnitEnv)
  testCaseDir <- file.path(system.file(package="RUnit"), "examples")
  testSuiteInternal <- defineTestSuite("RUnit Self Test", testCaseDir, "correctTestCase.r")
  testData2 <- runTestSuite(testSuiteInternal, useOwnErrorHandler=FALSE)

  timeStamp <- format(Sys.time(), "%y%m%d-%H%M")
  testProtocolFile <- file.path(tempdir(), paste(timeStamp, "test_printJUnitProtocol.xml", sep="_"))
  ret <- printJUnitProtocol(testData2, fileName=testProtocolFile)

  assign(".testLogger", tmp, envir = RUnitEnv)
  
  checkTrue( file.exists(testProtocolFile))


  ##  input argument error handling
  ##  missing 'testData' object
  checkException(printJUnitProtocol())

  ##  wrong class
  checkException(printJUnitProtocol("dummy"))

  
  ##  fileName arg errors
  testData <- list()
  class(testData) <- "RUnitTestData"
  ##  wrong type
  checkException(printJUnitProtocol(testData, fileName=numeric(1)))
  ##  wrong length
  checkException(printJUnitProtocol(testData, fileName=character(0)))
  checkException(printJUnitProtocol(testData, fileName=character(2)))
}

