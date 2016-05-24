message("TESTING: isBeingCreated()...")

library("R.oo")

setConstructorS3("Car", function(brand=NULL, nbrOfWheels=0) {
  if(!isBeingCreated(Car)) {
    if (is.null(brand))
      throw("A car must have a brand")

    if (nbrOfWheels <= 0)
      throw("A car must have one or more wheels: ", nbrOfWheels)
  }

  extend(Object(), "Car",
    .brand = brand,
    .nbrOfWheels = nbrOfWheels
  )
})

setMethodS3("as.character", "Car", function(this, ...) {
  cat(class(this)[1], ":", this$.brand, " with ",
                     this$.nbrOfWheels, " wheels.", sep="");
})

print(Car("Volvo", 4))
print(Car("BMW", 4))
print(Car("Tyrrell P34", 6))
print(Car("T-Rex", 3))

message("TESTING: isBeingCreated()...DONE")
