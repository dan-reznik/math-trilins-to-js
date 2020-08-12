library(tidyverse)
library(fs)
source("p5js_support.R")

convert_formulas <- function(fname_math_in,
                             vars_db,
                             fname_js_support) {
  fname_math_csv <- fs::path_ext_set(fname_math_in,"csv")
  # output goes to js directory
  fname_js_out <- fname_math_in %>%
    path_file() %>%
    path_ext_set("js") %>%
    str_c("js/",.)
  prepare_formulas(fname_math_in,
                   fname_math_csv)
  df_formulas <- read_csv2(fname_math_csv,col_types = "ccc")
  # read dependencies into first dictionary
  vars_db <- read_delim(vars_db, delim = fixed("="), col_names = F)
  vars_dict <- new.env()
  vars_db %>% pwalk(~(vars_dict[[..1]] <- ..2))
  vars_db_intersect <- vars_db %>% 
    mutate(dependence = intersect_vars(X1,X2))
  # iverton: documentar
  #vars_dict_eq <- new.env()
  vars_dict_dependence <- new.env()
  #vars_db_intersect %>% pwalk(~(vars_dict_eq[[..1]] <- ..2))
  vars_db_intersect %>% pwalk(~(vars_dict_dependence[[..1]] <- ..3))
  
  df_formulas_vars <- df_formulas %>%
    mutate(vars=trilins %>% process_trilins) # so reporta unicos
  codigo_js <- create_js_code(df_formulas_vars,vars_dict,vars_dict_dependence) 
  
  # concatena com support
  read_file(fname_js_support) %>%
    str_c("\n\n", codigo_js) %>%
    write_file(fname_js_out)
}

# input .txt: {"X(3)", "cosA|cosB|cosC", "CIRCUMCENTER"} 

# Rscript p5js.R "data/x0001_0200 cform v2b.txt"
main <- function() {
  args <- commandArgs(trailingOnly = TRUE)
  if (length(args)==1) {
    args <- commandArgs(trailingOnly = TRUE)
    convert_formulas(args[1],
                     "data/vars_db.txt",
                     "js/support_functions.js")
  } else
    print("Error: usage: Rscript p5js.R 'fname_math_in.txt'")
}

main()

