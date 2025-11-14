## check_strand(c(1,-1,0))
## check_strand(2)
## check_strand(c(1,NA,-1), na=0)

## check_strand(as.factor(c("+", "-", NA, ".")), na=".")
## check_strand("forward")

## check_strand(c(TRUE, FALSE, NA))
## check_strand(c(TRUE, FALSE, NA), na=TRUE)
## check_strand(c(TRUE, FALSE, NA), na=1) # na converted to lgl (TRUE)

## strand_int(c("+","-",".", NA)) %T>% print %>% is
## strand_int(c(TRUE, FALSE, NA)) %T>% print %>% is
## strand_int(c(1, -1, 0, NA)) %T>% print %>% is

## strand_chr(c("+","-",".", NA)) %T>% print %>% is
## strand_chr(c(TRUE, FALSE, NA)) %T>% print %>% is
## strand_chr(c(1, -1, 0, NA)) %T>% print %>% is

## strand_lgl(c("+","-",".", NA)) %T>% print %>% is
## strand_lgl(c(TRUE, FALSE, NA)) %T>% print %>% is
## strand_lgl(c(1, -1, 0, NA)) %T>% print %>% is
