### Imputación de datos usando Regression of order statistics (ROS)
### Autor: Tamara Ricardo
### Fecha de modificación:
# Tue May 28 14:27:28 2024 ------------------------------


# Carga paquetes ----------------------------------------------------------
# remotes::install_github("mikmart/censlm")
pacman::p_load(
  censlm,
  janitor,
  rio,
  tidyverse
)

# Carga datos -------------------------------------------------------------
data <- import("COVID_PAMI_clean.xlsx")

### Modifica datos
data <- data %>% 
  
  # tiempo a factor
  mutate(tiempo_cat = factor(tiempo), 
         .after = tiempo) %>% 
  
  # límite de detección en cada tiempo
  
  mutate(cut_off = min(ua[ua>0], na.rm = T), 
         .by = tiempo_cat, 
         .after = ua) %>% 
  
  # ua censurada
  mutate(ua_cen = if_else(ua == 0, cut_off, ua), 
         .after = cut_off) %>% 
  
  # quita NAs
  filter(!is.na(ua_cen))


## Imputa por censML
set.seed(123)


fit <- clm(log(data$ua_cen) ~ 1, left = log(data$cut_off))

### Añade datos imputados
data <- data %>% 
  
  ## reemplaza valores censurados por imputados
  mutate(ua_cen = if_else(ua_cen == cut_off, exp(imputed(fit)), ua) %>% 
           round(3))


### Guarda base
export(data, file = "COVID_PAMI_imputed.xlsx")
rm(list = ls())

