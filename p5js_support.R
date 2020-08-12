clean_lines <- function(lines) {
  lines %>%
    str_replace(fixed('", '),'";') %>%
    str_replace(fixed(', "'),';"') %>%
    str_remove("^\\{") %>%
    str_remove("\\},?$") %>%
    str_remove_all(fixed('\"')) %>%
    str_replace_all(" ([+-]) ","\\1") %>%
    str_replace_all(fixed(", "),",") %>%
    str_replace_all("Power\\(([[:alnum:] +-/*]+), ?2\\)","(\\1)*(\\1)") %>%
    str_replace_all("\\(([:word:]+)\\)", "\\1") %>%
    str_replace_all(fixed("Power"),"Math.pow")
}

clean_lines_old <- function(lines) {
  lines %>%
    str_replace(fixed('", '),'";') %>%
    str_replace(fixed(', "'),';"') %>%
    str_remove("^\\{") %>%
    str_remove("\\},?$") %>%
    str_replace(fixed("List("),"[") %>%
    str_replace(fixed(');"'),'];"')%>%
    str_remove_all(fixed('\"')) %>%
    str_replace_all(" ([+-]) ","\\1") %>%
    str_replace_all(fixed(", "),",") %>%
    str_replace_all(fixed("Power"),"Math.pow") %>%
    #Substitui separador entre trilineares por |
    str_replace_all("Math.pow\\(([^,]+),([^,]+)\\)", "Math.pow\\(\\1&\\2)") %>%
    str_replace_all(fixed(","), "|") %>%
    str_replace_all(fixed("&"), ",")
}

prepare_formulas <- function(raw_formulas_file, formulas_csv_file){
  lines <- read_lines(raw_formulas_file) %>%
    head(-1)
  lines_clean <- clean_lines(lines)
  c("kimberling;trilins;name",
    lines_clean) %>% write_lines(formulas_csv_file)
}

intersect_vars <- function(X1, X2) {
  variaveis_equacao <- X2 %>%
    str_split(boundary("word"))
  
  inter <- variaveis_equacao %>% map(intersect, X1)
  inter
}

delete_me <- function(x) {
  (x%in%c("","Math","pow","a","b","c",
          "cPi3", "cPi6","sPi3","sPi6" ))|
    (str_detect(x,"^[:digit:]+$"))
}

process_trilins <- function(trilins){
  trilins %>%
    #strip_first_last %>% # remove colchetes [ ... ]
    str_split("[^[:alnum:]]") %>% # split por nao-alfanums (+,-,/)
    map(~discard(.x,delete_me)) %>% # remove indesejaveis
    map(unique)
}
strip_first_last <- function(s) s%>%str_sub(start=2)%>%str_sub(end=-2)

add_dependence_vars <- function(vars, vars_all,vars_dict_dependence){
  if(length(vars)>0){
    for(var in vars){
      dependencies <- vars_dict_dependence[[var]]#get(var, env=vars_dict_dependence)
      if(length(dependencies)>0){
        vars_all <- vars_all %>% append(dependencies)
        vars_all <- add_dependence_vars(dependencies, vars_all,vars_dict_dependence)
      }
    }
  }
  return(vars_all)
}


# n: 131
# triple: [
#    secA*(-(S*sec2A)+sumT2)*(-((sec2B+sec2C)*sumS2)+2*sumT2),
#    secB*(-(S*sec2B)+sumT2)*(-((sec2A+sec2C)*sumS2)+2*sumT2),
#    secC*(-(S*sec2C)+sumT2)*(-((sec2A+sec2B)*sumS2)+2*sumT2)
#         ]
# vars:  c("secA","S","sec2A","sumT2","sec2B","sec2C","sumS2","secB","secC") 
create_function_js <- function(n,trilins,vars,vars_dict,vars_dict_dependence) {
  #n <- 1### TESTES
  #vars <- df_formulas_vars$vars[[n]]### TESTES
  #trilins <- df_formulas_vars$trilins[n]### TESTES
  s3 <- str_split(trilins,fixed("|"))
  vars_with_dependence_vars <- add_dependence_vars(vars, vars,vars_dict_dependence)
  if(length(vars_with_dependence_vars)>0){
    vars_with_dependence_vars_unique_invert <- vars_with_dependence_vars %>%
      rev() %>%
      unique()
  } else{
    vars_with_dependence_vars_unique_invert <- 
      vars_with_dependence_vars 
  }
  vars_block <- vars_with_dependence_vars_unique_invert %>% 
    map_chr(~str_c("   let ",.x,"=",get(.x,env=vars_dict),";")) %>%
    str_c(collapse="\n")
  
  s3 %>% map_chr(~str_glue(
    "function trilin_X{n}(orbit, [a, b, c]) {{",
    # inserir bloco de variaveis usando vars_dict ou inner_join
    "   /* begin vars */",
    "{vars_block}",
    "   /* end vars */",
    "   let v1 = {.x[1]};",
    "   let v2 = {.x[2]};",
    "   let v3 = {.x[3]};",
    "   let tris = [v1,v2,v3];",
    "   return trilin_to_cartesian(orbit, [a, b, c], tris);",
    "}}",
    .sep="\n"))
}

create_js_code <- function(df_formulas_vars,vars_dict,vars_dict_dependence){
  df_formulas_vars %>%
    mutate(index = row_number(),
           js = pmap_chr(list(index, trilins, vars), 
                     ~create_function_js(..1,..2,..3,vars_dict,vars_dict_dependence))) %>%
    pull(js) %>%
    str_c(collapse = "\n\n")
}