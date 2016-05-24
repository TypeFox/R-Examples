## ------------------------------------------------------------------------
library("aoos")
# Standard definition for a generic without default:
numeric : strLength(x, ...) %g% standardGeneric("strLength")

# A method for x:character
strLength(x ~ character, ...) %m% { nchar(x) }

# Kind of the default method for x:ANY
strLength(x ~ ANY, ...) %m% { strLength(as.character(x)) }

# Check that it works:
strLength(123)
strLength("ab")

## ------------------------------------------------------------------------
Test(x ~ numeric, y = list()) %type% {
  stopifnot(length(.Object@x) == 1)
  stopifnot(.Object@x > 0)
  .Object
}

Test : character : Child() %type% .Object

Test(2)
Child("Hej", x = 2)

## ------------------------------------------------------------------------
'numeric | character' : Either() %type% .Object
Either(1)
Either("Hello World!")

## ------------------------------------------------------------------------
Either(x ~ numeric | character) %type% .Object
Either(1)
Either("Hello World!")

## ------------------------------------------------------------------------
'numeric | character' : complicatedFunction(x = 1) %g% as.character(x)
complicatedFunction(x ~ character | integer) %m% as.numeric(x)
complicatedFunction()
complicatedFunction("1")
complicatedFunction(1L)

## ------------------------------------------------------------------------
library("aoos")
Employee <- function(.name, .salary) {
  "Common base class for all employees"
  
  print <- function(x, ...) {
    cat("Name  : ", .self$.name, "\nSalary: ", .self$.salary)
  }
  
  getName <- function() .name
  getSalary <- function() .self$.salary
  
  retList(c("Employee", "Print"))
  
}

peter <- Employee("Peter", 5)
peter
peter$getName()
peter$getSalary()

## ------------------------------------------------------------------------
peter
peter$print()

## ------------------------------------------------------------------------
Manager <- function(.name, .salary, .bonus) {
  "Extending the Employee class"
  
  bonus <- function(x) {
    if (!missing(x)) .self$.bonus <- x
    .self$.bonus
  }
  
  print <- function(x, ...) {
    cat("Name  : ", .self$.name, "\nSalary: ", .self$.salary, 
        "\nBonus:", .self$.bonus)
  }
  
  retList("Manager", super = Employee(.name, .salary))
  
}

julia <- Manager("Julia", 5, 5 * 1e6)
julia
julia$getSalary()
julia$bonus(10)
julia

## ------------------------------------------------------------------------
Class <- function() {
  
  overloaded(x) %g% { 
    cat("This is the default ... \n")
    x 
  } 
  
  overloaded(x ~ numeric) %m% {
    cat("This is the method for 'numeric' values ... \n")
    x
  }
  
  retList("Class")
}

instance <- Class()
instance$overloaded(1)
instance$overloaded("a")

## ------------------------------------------------------------------------
Child <- function() {
  
  # Normally you would make the call to the parents constructor in the call
  # to retList. But here we need to access the elements directly during init...
  .super <- Class()
  
  # This points %m% to the generic (in .super) which should be extended:
  .super$overloaded(x ~ integer) %m% {
    cat("This is the method for 'integer' values ... \n")
    x
  }
  
  retList("Child", super = .super)
  
}

instance <- Child()
instance$overloaded(1)
instance$overloaded("a")
instance$overloaded(1L)

