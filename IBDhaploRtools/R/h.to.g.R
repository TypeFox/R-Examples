## function to change from 15 state (haplotypic/phased) to 9 state
## (genotypic/unphased).

## We are assuming that we have the traditional ordering of states,
## not the lexicographical ordering.

## This is the pattern that I am using
## 1 -- 1
## 2 -- 2
## 3 -- 2
## 4 -- 3
## 5 -- 3
## 6 -- 4
## 7 -- 5 
## 8 -- 6
## 9 -- 7
## 10 -- 7
## 11 -- 8
## 12 -- 8
## 13 -- 8
## 14 -- 8
## 15 -- 9


h.to.g <- function( hap.states ){

  gen.states <- hap.states

  gen.states <- replace(gen.states, hap.states==1, 1)
  gen.states <- replace(gen.states, hap.states==2, 2)
  gen.states <- replace(gen.states, hap.states==3, 2)
  gen.states <- replace(gen.states, hap.states==4, 3)
  gen.states <- replace(gen.states, hap.states==5, 3)
  gen.states <- replace(gen.states, hap.states==6, 4)
  gen.states <- replace(gen.states, hap.states==7, 5)
  gen.states <- replace(gen.states, hap.states==8, 6)
  gen.states <- replace(gen.states, hap.states==9, 7)
  gen.states <- replace(gen.states, hap.states==10, 7)
  gen.states <- replace(gen.states, hap.states==11, 8)
  gen.states <- replace(gen.states, hap.states==12, 8)
  gen.states <- replace(gen.states, hap.states==13, 8)
  gen.states <- replace(gen.states, hap.states==14, 8)
  gen.states <- replace(gen.states, hap.states==15, 9)

  return(gen.states)
  
}
