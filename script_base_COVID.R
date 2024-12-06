### Limpieza de datos 
### Autor: Tamara Ricardo
### Fecha de modificación:
# Tue May 28 14:22:30 2024 ------------------------------


# Carga paquetes ----------------------------------------------------------
pacman::p_load(
  rio,
  gtools,
  janitor,
  tidyverse
  )


# Carga datos crudos ------------------------------------------------------
data <- import("COVID_Hogares.xlsx")

### Explora datos
tabyl(data$Tiempo)


# Limpia datos ------------------------------------------------------------
data_clean <- data %>% 
  
  # filtra datos día 365 y edad < 60
  filter(Tiempo != 365 & Edad >= 60) %>% 
  
  # une con base hogares
  left_join(import("listado con hogares.xlsx", sheet = 2)) %>% 
  
  # une con base de datos fechas
  left_join(import("Fechas_obito_brote.xlsx")) %>% 
  
  # estandariza nombres de columnas
  clean_names() %>% 
  rename(covid_estudio = covid_durante_estudio,
         grupo = estado_al_inicio) %>% 
  
  # región sanitaria
  mutate(region = case_when(
    str_detect(id, "A") ~ "UGLXI",
    str_detect(id, "X_") ~ "UGLX",
    str_detect(id, "VII_") ~ "UGLVII"), 
    .after = id) %>% 
  
  # categoriza edad
  mutate(edad_cat = quantcut(edad, q = 2), 
         .after = edad) %>% 
  
  # remueve datos cambio de hogar
  separate_wider_delim(hogar, names = c("hogar", NA), 
                       delim = "/", too_few = "align_start") %>% 
  
  # modifica variables
  mutate(
    # grupo de exposición
    grupo = if_else(grupo == "naive", "naïve", "no-naïve"),
    
    # título de anticuerpos
    titulo = if_else(ua == 0, "DLD", as.character(titulo)),
    
    # comorbilidades
    across(.cols = c(diabetes:fragilidad), .fns = ~ na_if(.x, "Nr")),
    
    # fechas
    across(.cols = starts_with("fecha"), .fns = ~ ymd(.x))
  ) %>% 
  
  # anonimiza hogar 
  mutate(hogar = paste(region, 
                       factor(hogar) %>% as.numeric() %>% sprintf("%02d", .), sep = "_"), 
         .by = region) %>% 
  
  # crea variables para el análisis
  mutate(
    # fallecido durante el estudio
    obito_estudio = if_else(fecha_obito <= max(fecha_muestra, na.rm = T),
                            "Si", "No", missing = "No"),
    
    # brote de COVID durante el estudio
    brote_estudio = if_else(fecha_brote <= max(fecha_muestra, na.rm = T),
                            "Si", "No", missing = "No"),
    
    .by = id) %>% 
  
  # ordena columnas
  select(id, region, hogar, contains("brote"), 
         sexo, contains("edad"), contains("obito"), everything())

# # Guarda base
# export(data_clean, file = "COVID_PAMI-60.xlsx")


## Explora datos
tabyl(data_clean$muestra_tomada)

tabyl(data_clean$censurado_t1)

tabyl(data_clean$titulo)

# Filtra fallecidos t1 ----------------------------------------------------
data_clean <- data_clean %>% 
  
  # filtra fallecidos t1
  filter(muestra_tomada == "Si" & censurado_t1 == "No") %>% 
  
  # crea variable para datos por debajo del límite de detección
  mutate(censored = if_else(titulo == "DLD", TRUE, FALSE)) %>% 
  
  # descarta columnas innecesarias
  select(!c(fragilidad:censurado_t1))


# ## Guarda base limpia
# export(data_clean, file = "COVID_PAMI_clean.xlsx")

# Crea base wide ----------------------------------------------------------
data_wide <- data_clean %>% 
  
  # descarta columnas innecesarias
  select(-fecha_muestra, -titulo) %>% 
  
  # variables caracter a factor
  mutate(across(.cols = where(is.character), .fns = ~ as.factor(.x))) %>% 
  
  # formato wide
  pivot_wider(values_from = "ua", names_from = "tiempo", names_prefix = "T")


# ## Guarda base wide
# export(data_wide, file = "COVID_PAMI_wide.xlsx")
