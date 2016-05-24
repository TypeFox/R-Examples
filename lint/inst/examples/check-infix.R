# spacing around infix operators +, -, <, >, etc.
{ # should catch
# ------------
{ # Arithmatic operators
1-1   # minus
1 -1  # space before
1- 1  # space after
1+1   # plus 
1*1   # multiply
1/1   # divide
1^1   # exponent
1**1  # exponent 2
}
{ # Assignment
a<-1  # left assign
1->a  # right assign
a<<-1 # double left assign
1->>a # double right assign
}
{ # logical
1<1   # less
1<=1  # less or equal
1==1  # equal
1!=1  # not equal
1>=1  # greater or equal
1>1   # greater 
T|F   # binary or
T||F  # logical or
T&F   # binary and
T&&F  # binary and
}
{ # Specials
1%%1          # modulo
1%/%1         # integer divide
1%*%1         # matrix product
1%o%1         # outer product
1%x%1         # Kronecker
1%in%1        # matching operator
1%.%1         # special compose defined in package dostats
1%test%1      # included to check general tests.
}
}
{ # Should not catch
# ----------------
!T    # unary not
1:2   # sequence
a$b   # list element extraction
a@b   # slot extraction

{ # Arithmatic operators
1 - 1   # minus
1 + 1   # plus 
1 * 1   # multiply
1 / 1   # divide
1 ^ 1   # exponent
1 ** 1  # exponent
a = 1   # assign =
1 < 1   # less
1 <= 1  # less or equal
1 == 1  # equal
1 != 1  # not equal
1 >= 1  # Greater or equal
1 > 1   # greater
}
{ # Assignment
a <- 1  # left assign
1 -> a  # right assign
a <<- 1 # double left assign
1 ->> a # double right assign
}
{ # logical
T | F   # binary or
T || F  # logical or
T & F   # binary and
T && F  # binary and
}
{ # Specials
1 %% 1          # modulo
1 %/% 1         # integer divide
1 %*% 1         # matrix product
1 %o% 1         # outer product
1 %x% 1         # Kronecker
1 %in% 1        # matching operator
1 %.% 1         # special compose defined in package dostats
1 %test% 1      # included to check general tests.
}
}

{ # for equals check
# should catch
a=1
a =1
a= 1
# should not
a = 1
a(b=1) # function call
}
