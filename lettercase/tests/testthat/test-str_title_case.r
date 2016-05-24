library(testthat)
library(lettercase)
library(magrittr)

context( 'str_title_case')

str_title_case( 'mission of burma' ) %>% expect_equal( 'Mission Of Burma' ) 
str_title_case( 'mission_of_burma' ) %>% expect_equal( 'Mission Of Burma' )
str_title_case( 'Mission of Burma' ) %>% expect_equal( 'Mission Of Burma' )

# str_title_case( 'mission-of-burma' ) %>% expect_equal( 'Mission Of Burma' )
# str_title_case( 'MISSION OF BURMA' ) %>% expect_equal( 'Mission Of Burma' )
# str_title_case( 'MissionOfBurma' )   %>% expect_equal( 'Mission Of Burma' )

