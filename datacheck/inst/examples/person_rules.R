is.integer(id)  # id
id > 0 & id < 11  # id

is.character(lastname)  #lastname
is_proper_name(lastname)  #lastname

is.character(firstname)  #firstname
is_proper_name(firstname)  #firstname

is.character(gender)  #gender
is_one_of(gender, c("m", "f"))  #gender

age > 10 & age < 70 
