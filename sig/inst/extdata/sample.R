#Create some well named functions
calc_third_angle_of_isoceles_triangle_in_radians <- function(paired_angle)
{
  pi - 2 * paired_angle
}
append_csv_extension_to_filename <- function(file) paste0(file, ".CSV")


#and some badly named functions
f <- function() 1
g <- function(x, y, z) x + y + z

#and some other variables that shouldn't appear in output
a <- 1
b <- TRUE
ex <- expression(h <- function() 2)

