### Análisis de datos para el manuscrito: Efecto de la exposición previa a COVID-19,
### ocurrencia de brotes y tipo de vacuna en la respuesta inmune humoral de 
### adultos mayores institucionalizados
### Autor: Tamara Ricardo
### Última modificación
# Tue May 28 14:27:34 2024 ------------------------------


# Carga paquetes ----------------------------------------------------------
pacman::p_load(
  # GLS
  nlme,
  # GAMM
  mgcv, 
  gratia,
  itsadug,
  # Residuales
  DHARMa,
  performance,
  # Selección variables
  MuMIn,
  # Análisis exploratorio
  gtsummary,
  # Formato tablas
  flextable,
  # Manejo de datos
  gridExtra,
  rio,
  gtools,
  janitor,
  tidyverse
)


# Carga dataset -----------------------------------------------------------
data <- import("COVID_PAMI_imputed.xlsx") %>% 
  
  # caracteres a factor
  mutate(across(.cols = where(is.character),
                .fns = ~ factor(.x))) %>% 
  
  # ua a logaritmo
  mutate(ua_log = log(ua_cen)) %>% 
  
  # tratamiento
  mutate(tratamiento = fct_cross(grupo, vacuna) %>% 
           fct_relevel("naïve:Sinopharm","naïve:Sputnik-V", after = 1),
         .after = grupo) %>% 
  
  # filtra datos faltantes comorbilidades
  filter(!is.na(diabetes)) %>% 
  
  # descarta niveles vacíos
  droplevels()


# Análisis exploratorio ---------------------------------------------------
### nro pacientes
fct_count(data$id) %>% nrow()

### nro de hogares
fct_count(data$hogar) %>% nrow()

### hogares x región
data %>% 
  filter(tiempo == 0) %>% 
  select(hogar, region) %>% 
  distinct() %>% 
  tabyl(region) %>% 
  adorn_pct_formatting()
  

### Descripción grupos exposición
data %>% 
  filter(tiempo == 0) %>% 
  tabyl(grupo, vacuna) %>% 
  adorn_totals() %>% 
  adorn_percentages() %>% 
  adorn_pct_formatting()

### Descripción edades
data %>% 
  filter(tiempo == 0) %>% 
  tbl_summary(by = sexo,
              include = edad_cat) %>% 
  add_p()


### Comorbilidades
data %>% 
  filter(tiempo == 0) %>% 
  
  tbl_summary(
    by = covid_estudio,
    include = c(diabetes:inmunodeficiencia)) %>% 
  add_overall() %>% 
  add_p()

### Spaghetti plots
ggpubr::ggarrange(
  # variable respuesta: ua_imp
  ggplot(data, 
         aes(x = tiempo_cat, 
             y = ua_cen, 
             group = id, 
             color = tratamiento)) +
    
    geom_line() +
    
    labs(y = "ua", x = "tiempo (días)") +
    
    scale_color_brewer(palette = "Set2") +
    
    theme_light()
  ,
  # Variable respuesta: log(ua)
  ggplot(data, 
         aes(x = tiempo_cat, 
             y = ua_log, 
             group = id, 
             color = tratamiento)) +
    
    geom_line() +
    
    labs(y = "log(ua)", x = "tiempo (días)") +

    scale_color_brewer(palette = "Set2") +
    
    theme_light(), 
  
  common.legend = T, legend = "bottom", labels = "AUTO", font.label = list(size=12))


# Figura 1: gráfico de perfiles -------------------------------------------
g1 <- ggplot(data = data, 
             aes(x = tiempo, 
                 y = ua_log, 
                 color = tratamiento, 
                 group = tratamiento)) +
  
  geom_smooth(method = "loess", se = F) + 
  
  # Colorblind friendly
  scale_color_brewer(palette = "Set2") +

  labs(y = "log(ua)") +
  
  scale_x_continuous(breaks = c(0,21,42,120,180), minor_breaks = F) +
  
  theme_light() + 
  theme(legend.title = element_blank())


### Gráfico con leyenda
g1 + theme(legend.position = "bottom")

## Guarda gráfico
ggsave("FIGS/Fig1cl.svg", width = 16, height = 10, units = "cm", dpi = 300)


# Modelos marginales (GLS) ------------------------------------------------
### Autocorrelación temporal
# Plot time series
plot(data$tiempo, data$ua_log, 
     type = "b", 
     xlab = "Time", ylab = "Antibody Levels")

# Calculate autocorrelation
acf_result <- acf(data$ua_log, 
                  main = "Autocorrelation Function")


## test for independance 
Box.test(data$ua_log, type = "Ljung-Box")


### Autorregresiva continua de orden 1 c/varianza homogénea
gls1 <- gls(ua_log ~ tratamiento * tiempo_cat +
              sexo + edad_cat + covid_estudio + brote_estudio + 
              diabetes + hta + insuf_cardiaca + inmunodeficiencia,
            correlation = corCAR1(form = ~1|hogar/id), 
            data = data)


### Autorregresiva continua de orden 1 c/varianza heterogénea
gls2 <-  gls(ua_log ~ tratamiento * tiempo_cat +
               sexo + edad_cat + covid_estudio + brote_estudio + 
               diabetes + hta + insuf_cardiaca + inmunodeficiencia,
             correlation = corCAR1(form = ~1|hogar/id), 
             weights = varIdent(form= ~ 1|tiempo_cat),
             data = data)


### Compara modelos marginales----
compare_performance(gls1, gls2, rank = T)

AIC(gls1) - AIC(gls2)

### Selección variables
anova(gls2)

## (-) covid_estudio
gls2a <- update(gls2, ~.-covid_estudio)

anova(gls2a)

## (-) sexo
gls2b <- update(gls2a, ~.-sexo)

anova(gls2b)

## (-) edad_cat
gls2c <- update(gls2b, ~.-edad_cat)

anova(gls2c)

## (-) inmunodeficiencia
gls2d <- update(gls2c, ~.-inmunodeficiencia)

anova(gls2d)

## (-) insuf_cardiaca
gls2e <- update(gls2d, ~.-insuf_cardiaca)

anova(gls2e)

## (-) hta
gls2f <- update(gls2e, ~.-hta)

anova(gls2f)

## (-) diabetes
gls2g <- update(gls2f, ~.-diabetes)

anova(gls2g)

## (Compara modelos)
compare_performance(gls2, gls2a, gls2b, gls2c, gls2d, gls2e, gls2f, gls2g, 
                    rank = F)


# Tabla 1 -----------------------------------------------------------------
# Selección de variables explicativas en el modelo GLS con estructura de correlación temporal autorregresiva continua de primer orden (corCAR1) y heterogeneidad de varianzas (varIdent). Todos los modelos incluyen la interacción entre tiempo y tratamiento Los signos + y - indican inclusión o exclusión de la variable explicativa.

## Crea datos para la tabla
fit <- model.sel(gls2, gls2a, gls2b, gls2c, gls2d, gls2e, gls2f, 
                 gls2g, rank = "AIC") %>% 
  
  as_tibble(rownames = "Modelo") %>%
  
  mutate(across(logLik:delta, .fns = ~ round(.x, 2)),
         
         across(where(is.factor), .fns = ~ as.character(.x) %>% replace_na("-")))


## Exporta tabla
fit %>% select(Modelo, brote_estudio, covid_estudio, sexo, edad_cat,
               diabetes, hta, insuf_cardiaca, inmunodeficiencia, AIC, delta) %>% 
  
  flextable() %>% 
  
  autofit() %>% 
  
  font(fontname = "Times New Roman", part = "all") %>% 
  fontsize(size = 12, part = "all") %>% 
  bold(part = "header") %>% 
  
  save_as_docx(path = "./FIGS/tab1.docx")


# Modela respuesta no lineal (GAMM) ---------------------------------------
### GAMM corCar1 y varianzas heterogéneas
gamm <- gamm(ua_log ~ brote_estudio + tratamiento + # Términos paramétricos
               
               s(tiempo, by = tratamiento, k = 5), # Curva suavizado interacción
             
             random = list(id = ~ 1, hogar = ~ 1), # Efecto aleatorio paciente y EEP
             
             correlation = corCAR1(form = ~1|id), # Autocorrelación temporal
             
             
             weights = varIdent(form= ~ 1|tiempo_cat), # Varianzas heterogéneas
             
             data = data, 
             
             method = "REML")


# Resumen términos no lineales
summary(gamm$gam)

# Resumen términos paramétricos
summary(gamm$lme)

## Análisis de residuales
acf(gamm$lme %>% residuals())

pacf(gamm$lme %>% residuals())

appraise(gamm$gam)

## Grados de libertad efectivos
edf(gamm$gam)


# Figura 2: curvas suavizado ----------------------------------------------
svg(filename = "FIGS/Fig2cl.svg", width = 6.3, height = 7.9)
par(mfrow = c(2,1), cex = .75)
# Brote: Si
plot_smooth(gamm$gam, view = "tiempo", plot_all = "tratamiento", 
            cond = list(brote_estudio="Si"), v0 = c(0,21,42,120,180),
            rug = F, rm.ranef = F, legend_plot_all = T,
            ylim = c(2,11), xlim = c(0,200), xlab = "",
            ylab = "", hide.label = F, main = "(A)", adj = 0,
            # Colorblind friendly
            col = c("#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854","#FFD92F"))


# Brote: No
plot_smooth(gamm$gam, view = "tiempo", plot_all = "tratamiento", 
            cond = list(brote_estudio="No"), v0 = c(0,21,42,120,180),
            rug = F, rm.ranef = F,legend_plot_all = F,
            ylim = c(2,11), xlim = c(0,200), xlab = "",
            ylab = "", hide.label = F, main = "(B)", adj = 0,
            # Colorblind friendly
            col = c("#66C2A5","#FC8D62","#8DA0CB","#E78AC3","#A6D854","#FFD92F"))
dev.off()



# Figura 3: términos paramétricos -----------------------------------------
svg(filename = "FIGS/Fig3.svg", width = 6.3, height = 7.9)
par(mfrow = c(2,1), cex = .75)
plot_parametric(gamm$gam, pred = list(tratamiento = levels(data$tratamiento)),
                cond = list(brote_estudio="Si"), xlim = c(4, 12),
                xlab = "log(ua)", main = "")
title(main = "(A)", adj = 0)

plot_parametric(gamm$gam, pred = list(tratamiento = levels(data$tratamiento)),
                cond = list(brote_estudio="No"),  xlim = c(4, 12),
                xlab = "log(ua)", main = "")
title(main = "(B)", adj = 0)


dev.off()

rm(list = ls())
