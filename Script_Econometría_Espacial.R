#-------------------------------------------------------------------------------
#                     TRABAJO ECONOMETRÍA ESPACIAL
#-------------------------------------------------------------------------------
#NOMBRE: Fernando Yunga
#FECHA: 12 de mayo de 2025

#Cargar la librerias
library(spatialreg) # Econometria espacial
library(spdep) # Econometria espacial
library(sf) # leer archivos shapfiles y elaborar mapas 
library(RColorBrewer) # Eleccion de colores
library(classInt) # metodos para clasificar
library(readxl)
library(relimp)
library(stats)
library(lattice)

#Cambiar el directorio de trabajo
setwd("D:/Ferchito/Desktop/UNAM/CICLO 4/ECONOMETRIA ESPACIAL/PROYECTO/nxcantones")

#Leer archivos shapes y transformarlo en objeto Shape y DataFrame 
cant_base1 <- read_sf("nxcantones.shp")
View(cant_base1)
cant_base2 <- cant_base1[-c(219, 220, 221), ]# quitamos cantones sin información

#III.	Análisis Exploratorio de Datos--------------------------------------------

#a.	Estadísticas básicas de indicadores------------- 
attach(cant_base2)
summary(cant_base2[, c("POBREZA", "GR_ESC", "NEMPRESAS")])

#b.	Función kernel para cada variable---------------
lpob <- log(POBREZA+1)
lesc <- log(GR_ESC+1)
lnempr <- log(NEMPRESAS+1)
densityplot( ~ lpob, data=cant_base2, bw="SJ", adjust=1, kernel="gaussian")
densityplot( ~ lesc, data=cant_base2, bw="SJ", adjust=1, kernel="gaussian")
densityplot( ~ lnempr, data=cant_base2, bw="SJ", adjust=1, kernel="gaussian")

#c.	Correlación de variable endógena con cada una de las exógenas-------
cor.test(lpob, lesc,method = c("pearson"))
cor.test(lpob, lnempr,method = c("pearson"))

#d.	Diagrama de relación (dispersión)------------
plot(lesc,lpob)
plot(lnempr, lpob)

#e.	Construcción y análisis de matriz W_ Queen o K-vecinos-------------
# Coordenadas
coords <- st_centroid(cant_base2)
# Vecinos cercanos
kv <- knn2nb(knearneigh(coords, k=5))
# Matriz de ponderacion 5 Vecinos estandarizada
wkv <- nb2listw(kv,style="W", zero.policy=TRUE)
View(wkv)
print(wkv)
#  Geometría de los cantones
coords <- st_coordinates(st_centroid(cant_base2))
plot(st_geometry(cant_base2), border = "gray", main = "Red de 5 Vecinos más Cercanos")
# Enlaces de vecindad
plot(kv, coords, add = TRUE, col = "blue", lwd = 1.2)
# Puntos con centroides
points(coords, pch = 20, col = "red")

#f.	Análisis con estadísticos de Moran----------------------------
##i.	Diagramas de dispersión (variables en logaritmo)
moran.plot(lpob, wkv, pch=20, main="Diagrama de dispersion de Moran para pobreza")
moran.plot(lesc, wkv, pch=20, main="Diagrama de dispersion de Moran para escolaridad")
moran.plot(lnempr, wkv, pch=20,main="Diagrama de dispersion de Moran para número de empresas")

##ii.	Test de correlación espacial (variables en logaritmo)
moran.test(lpob, wkv,randomisation=TRUE, alternative="two.sided", na.action=na.exclude)
moran.test(lesc, wkv,randomisation=TRUE, alternative="two.sided", na.action=na.exclude)
moran.test(lnempr, wkv,randomisation=TRUE, alternative="two.sided", na.action=na.exclude)

#g.	Análisis de dependencia local (LISA)--------------------------
##i.	Mapas de significancia
##ii.	Mapas con clúster espacial

### Para el log de pobreza----

lmoran <- localmoran(lpob, wkv)
View(lmoran)
Prm <- lmoran[,"Pr(z != E(Ii))"]
summary(Prm)
Prb <- Prm
Prb <- replace(Prb,Prm<1,0)
summary(Prb)
# Signficativos al 90% (valor 2)
Prb <- replace(Prb,Prm<=0.1,2)
Prb <- replace(Prb,Prm<=0.075,2.5)
# Signficativos al 95% (valor 1)
Prb <- replace(Prb,Prm<=0.05,1)
Prb <- replace(Prb,Prm<=0.025,1.5)
table(Prb)
pal_map <- c("white", "green", "darkgreen")
pal <- c("white", "green", "darkgreen")
z3.Ii <- classIntervals(Prb, n=3, style="fixed", fixedBreaks=c(0,1,2,3))
cols.Ii <- findColours(z3.Ii, pal_map)
plot(cant_base2$geometry,border= "gray")
plot(cant_base2$geometry, col=cols.Ii, border= FALSE,cex =.1,add=T,)
invisible(title(main=paste("Probabilidades con Analisis LISA para el logaritmo de pobreza", sep="\n")))
leg <- c("no.sig (169)","Sig. 95% (29)","Sig. 90% (20)")
legend("topleft", fill=pal, legend=leg, bty="n",cex =1)    
box()

#### Mapa de los grupos de cluster
Prm <- lmoran[,"Pr(z != E(Ii))"]
summary(Prm)
lm_pob <- Prm
summary(lm_pob)
lm_pob <- replace(lm_pob,Prm<1,0)
# Rezagos espaciales de las variables en logaritmo para construir Analisis LISA
wlpob <- lag.listw(wkv,lpob)
#wlch <- lag.listw(wkv,log(1+ch))
# High-High  (Vecinos (High)-localidad (High))
lm_pob <- replace(lm_pob,Prm<=0.1 & lpob>mean(lpob) & wlpob>mean(wlpob),1)
lm_pob <- replace(lm_pob,Prm<=0.05 & lpob>mean(lpob) & wlpob>mean(wlpob),1.5)
# Low-Low  (Vecinos (Low)-localidad (Low))
lm_pob <- replace(lm_pob,Prm<=0.1 & lpob<=mean(lpob) & wlpob<=mean(wlpob),2)
lm_pob <- replace(lm_pob,Prm<=0.05 & lpob<=mean(lpob) & wlpob<=mean(wlpob),2.5)
# Low-High  (Vecinos (Low)-localidad (High))
lm_pob <- replace(lm_pob,Prm<=0.1 & lpob>mean(lpob) & wlpob<=mean(wlpob),3)
lm_pob <- replace(lm_pob,Prm<=0.05 & lpob>mean(lpob) & wlpob<=mean(wlpob),3.5)
# High-Low  (Vecinos (High)-Localidad (Low))
lm_pob <- replace(lm_pob,Prm<=0.1 & lpob<=mean(lpob) & wlpob>mean(wlpob),4)
lm_pob <- replace(lm_pob,Prm<=0.05 & lpob<=mean(lpob) & wlpob>mean(wlpob),4.5)
table(lm_pob)

# Mapa de Clusteres
pal_map <- c("white", "firebrick1", "darkblue","slategray3","rosybrown2")
pal <- c("white", "firebrick1", "darkblue","slategray3","rosybrown2")
#pal_map <- c("white", "firebrick1", "darkblue")
#pal <- c("white", "firebrick1", "darkblue","slategray3","rosybrown2")
z3.Ii <- classIntervals(lm_pob, n=5, style="fixed", fixedBreaks=c(0,1,2,3,4,5))
cols.Ii <- findColours(z3.Ii, pal_map)
plot(cant_base2$geometry,border= "gray")
plot(cant_base2$geometry, col=cols.Ii, border= FALSE,cex =.1, add=T)
invisible(title(main=paste("Clusteres Espaciales (LISA) para el logaritmo de pobreza", sep="\n")))
leg <- c("no.sig (169)","High-High (17)","Low-Low (23)","Low-High (3)","High-Low (6)")
legend("topleft", fill=pal, legend=leg, bty="n",cex =1)    
box()

### Para el log de escolaridad--------

lmoran <- localmoran(lesc, wkv)
View(lmoran)
Prm <- lmoran[,"Pr(z != E(Ii))"]
summary(Prm)
Prb <- Prm
Prb <- replace(Prb,Prm<1,0)
summary(Prb)
# Signficativos al 90% (valor 2)
Prb <- replace(Prb,Prm<=0.1,2)
Prb <- replace(Prb,Prm<=0.075,2.5)
# Signficativos al 95% (valor 1)
Prb <- replace(Prb,Prm<=0.05,1)
Prb <- replace(Prb,Prm<=0.025,1.5)
table(Prb)

pal_map <- c("white", "green", "darkgreen")
pal <- c("white", "green", "darkgreen")
z3.Ii <- classIntervals(Prb, n=3, style="fixed", fixedBreaks=c(0,1,2,3))
cols.Ii <- findColours(z3.Ii, pal_map)
plot(cant_base2$geometry,border= "gray")
plot(cant_base2$geometry, col=cols.Ii, border= FALSE,cex =.1,add=T,)
invisible(title(main=paste("Probabilidades con Analisis LISA para el logaritmo de escolaridad", sep="\n")))
leg <- c("no.sig (180)","Sig. 95% (20)","Sig. 90% (18)")
legend("topleft", fill=pal, legend=leg, bty="n",cex =1)    
box()

#### Mapa de los grupos de cluster
Prm <- lmoran[,"Pr(z != E(Ii))"]
summary(Prm)
lm_esc <- Prm
summary(lm_esc)
lm_esc <- replace(lm_esc,Prm<1,0)
# Rezagos espaciales de las variables en logaritmo para construir Analisis LISA
wlesc <- lag.listw(wkv,lesc)
# High-High  (Vecinos (High)-localidad (High))
lm_esc <- replace(lm_esc,Prm<=0.1 & lesc>mean(lesc) & wlesc>mean(wlesc),1)
lm_esc <- replace(lm_esc,Prm<=0.05 & lesc>mean(lesc) & wlesc>mean(wlesc),1.5)
# Low-Low  (Vecinos (Low)-localidad (Low))
lm_esc <- replace(lm_esc,Prm<=0.1 & lesc<=mean(lesc) & wlesc<=mean(wlesc),2)
lm_esc <- replace(lm_esc,Prm<=0.05 & lesc<=mean(lesc) & wlesc<=mean(wlesc),2.5)
# Low-High  (Vecinos (Low)-localidad (High))
lm_esc <- replace(lm_esc,Prm<=0.1 & lesc>mean(lesc) & wlesc<=mean(wlesc),3)
lm_esc <- replace(lm_esc,Prm<=0.05 & lesc>mean(lesc) & wlesc<=mean(wlesc),3.5)
# High-Low  (Vecinos (High)-Localidad (Low))
lm_esc <- replace(lm_esc,Prm<=0.1 & lesc<=mean(lesc) & wlesc>mean(wlesc),4)
lm_esc <- replace(lm_esc,Prm<=0.05 & lesc<=mean(lesc) & wlesc>mean(wlesc),4.5)
table(lm_esc)

# Mapa de Clusteres
pal_map <- c("white", "firebrick1", "darkblue","slategray3","rosybrown2")
pal <- c("white", "firebrick1", "darkblue","slategray3","rosybrown2")
#pal_map <- c("white", "firebrick1", "darkblue")
#pal <- c("white", "firebrick1", "darkblue","slategray3","rosybrown2")
z3.Ii <- classIntervals(lm_esc, n=5, style="fixed", fixedBreaks=c(0,1,2,3,4,5))
cols.Ii <- findColours(z3.Ii, pal_map)
plot(cant_base2$geometry,border= "gray")
plot(cant_base2$geometry, col=cols.Ii, border= FALSE,cex =.1, add=T)
invisible(title(main=paste("Clusteres Espaciales (LISA) para el logaritmo de escolaridad", sep="\n")))
leg <- c("no.sig (180)","High-High (13)","Low-Low (17)","Low-High (3)","High-Low (5)")
legend("topleft", fill=pal, legend=leg, bty="n",cex =1)    
box()

### Para el log de número de empresas------

lmoran <- localmoran(lnempr, wkv)
View(lmoran)
Prm <- lmoran[,"Pr(z != E(Ii))"]
summary(Prm)
Prb <- Prm
Prb <- replace(Prb,Prm<1,0)
summary(Prb)
# Signficativos al 90% (valor 2)
Prb <- replace(Prb,Prm<=0.1,2)
Prb <- replace(Prb,Prm<=0.075,2.5)
# Signficativos al 95% (valor 1)
Prb <- replace(Prb,Prm<=0.05,1)
Prb <- replace(Prb,Prm<=0.025,1.5)
table(Prb)

pal_map <- c("white", "green", "darkgreen")
pal <- c("white", "green", "darkgreen")
z3.Ii <- classIntervals(Prb, n=3, style="fixed", fixedBreaks=c(0,1,2,3))
cols.Ii <- findColours(z3.Ii, pal_map)
plot(cant_base2$geometry,border= "gray")
plot(cant_base2$geometry, col=cols.Ii, border= FALSE,cex =.1,add=T,)
invisible(title(main=paste("Probabilidades con Analisis LISA para el número de empresas", sep="\n")))
leg <- c("no.sig (190)","Sig. 95% (15)","Sig. 90% (13)")
legend("topleft", fill=pal, legend=leg, bty="n",cex =1)    
box()

#### Mapa de los grupos de cluster
Prm <- lmoran[,"Pr(z != E(Ii))"]
summary(Prm)
lm_nempr <- Prm
summary(lm_nempr)
lm_nempr <- replace(lm_nempr,Prm<1,0)
# Rezagos espaciales de las variables en logaritmo para construir Analisis LISA
wlnempr <- lag.listw(wkv,lnempr)
# High-High  (Vecinos (High)-localidad (High))
lm_nempr <- replace(lm_nempr,Prm<=0.1 & lnempr>mean(lnempr) & wlnempr>mean(wlnempr),1)
lm_nempr <- replace(lm_nempr,Prm<=0.05 & lnempr>mean(lnempr) & wlnempr>mean(wlnempr),1.5)
# Low-Low  (Vecinos (Low)-localidad (Low))
lm_nempr <- replace(lm_nempr,Prm<=0.1 & lnempr<=mean(lnempr) & wlnempr<=mean(wlnempr),2)
lm_nempr <- replace(lm_nempr,Prm<=0.05 & lnempr<=mean(lnempr) & wlnempr<=mean(wlnempr),2.5)
# Low-High  (Vecinos (Low)-localidad (High))
lm_nempr <- replace(lm_nempr,Prm<=0.1 & lnempr>mean(lnempr) & wlnempr<=mean(wlnempr),3)
lm_nempr <- replace(lm_nempr,Prm<=0.05 & lnempr>mean(lnempr) & wlnempr<=mean(wlnempr),3.5)
# High-Low  (Vecinos (High)-Localidad (Low))
lm_nempr <- replace(lm_nempr,Prm<=0.1 & lnempr<=mean(lnempr) & wlnempr>mean(wlnempr),4)
lm_nempr <- replace(lm_nempr,Prm<=0.05 & lnempr<=mean(lnempr) & wlnempr>mean(wlnempr),4.5)
table(lm_nempr)

# Mapa de Clusteres
pal_map <- c("white", "firebrick1", "darkblue","slategray3","rosybrown2")
pal <- c("white", "firebrick1", "darkblue","slategray3","rosybrown2")
#pal_map <- c("white", "firebrick1", "darkblue")
#pal <- c("white", "firebrick1", "darkblue","slategray3","rosybrown2")
z3.Ii <- classIntervals(lm_nempr, n=5, style="fixed", fixedBreaks=c(0,1,2,3,4,5))
cols.Ii <- findColours(z3.Ii, pal_map)
plot(cant_base2$geometry,border= "gray")
plot(cant_base2$geometry, col=cols.Ii, border= FALSE,cex =.1, add=T)
invisible(title(main=paste("Clusteres Espaciales (LISA) para el logaritmo de número de empresas", sep="\n")))
leg <- c("no.sig (190)","High-High (14)","Low-Low (5)","Low-High (7)","High-Low (2)")
legend("topleft", fill=pal, legend=leg, bty="n",cex =1)    
box()

#IV. Análisis Confirmatorio modelo de corte transversal-------------------------

##a.	Modelo Estático (lineal o log-log)---------------

###i.	Análisis econométricos y económicos 
mod_ols <- lm(lpob ~ lesc+lnempr, data=cant_base2)
summary(mod_ols)

###ii.	Mapas de terciles y quintiles de los errores
cant_base2$residuos <- residuals(mod_ols)
# Quintiles
brks_quint <- round(quantile(cant_base2$residuos, probs = seq(0.00, 1, 0.201), na.rm = TRUE), digits = 2)
colores_quint <- brewer.pal(5, "Reds")
plot(cant_base2$geometry, col = colores_quint[findInterval(cant_base2$residuos, brks_quint, all.inside = TRUE)],
     axes = FALSE, border = FALSE)
legend("topleft", legend = brks_quint, fill = colores_quint, bty = "n", cex = 1)
title(main = "Quintiles de residuos del modelo OLS")
box()
# Terciles
brks_terc <- round(quantile(cant_base2$residuos, probs=seq(0,1,0.34), na.rm = TRUE), digits=2)
colores_terc <- brewer.pal(3, "Blues")
plot(cant_base2$geometry, col = colores_terc[findInterval(cant_base2$residuos, brks_terc, all.inside = TRUE)],
     axes = FALSE, border = FALSE)
legend("topleft", legend = brks_terc, fill = colores_terc, bty = "n", cex = 1)
title(main = "Terciles de residuos del modelo OLS")
box()

###iii.	Test para dependencia espacial: MORAN, LM-Lag
# Prueba de Moran a residuales del modelo OLS
I_Moran <- lm.morantest(mod_ols,wkv)
print(I_Moran)
#Pruebas de Multiplicadores de Lagranges lm.LMtests(columbus.lm,col.listw,test=c("LMerr","RLMerr","LMlag","RLMlag","SARMA"))
lm.LMtests(mod_ols,wkv,test=c("LMlag", "LMerr", "RLMlag", "RLMerr","SARMA"))

#b.	Modelo Rezago Espacial-------------------

modelo_lag <- lagsarlm(lpob ~ lesc+lnempr, data=cant_base2,wkv)
summary(modelo_lag)

#c.	Modelo Error Espacial--------------------

modelo_err <- errorsarlm(lpob ~ lesc+lnempr, data=cant_base2,wkv)
summary(modelo_err)

#d.	Modelo Durbin Espacial------------------

modelo_durbin <- lagsarlm(lpob ~ lesc+lnempr, data=cant_base2,wkv, type="mixed")
summary(modelo_durbin)

#e.	Modelo de Rezago y Error Espacial (SARMA)---------------

modelo_sarar <- sacsarlm(lpob ~ lesc+lnempr, data=cant_base2,wkv, type="sac")
summary(modelo_sarar)

#f.	Análisis de efectos espaciales (directo, indirecto y total)------------

# Evaluar impactos en modelos espaciales
#Matrices para el calculo de la traza
W <- as(as_dgRMatrix_listw(wkv), "CsparseMatrix") #Matrices dispersas (sparse)
trMatc <- trW(W, type="mult") # traza con potencia de matrices dispersas
trMC <- trW(W, type="MC") # Simulacion de Monte Carlo de las traza

# Simulacion de impactos modelo rezago espacial
# Impactos 
summary(modelo_lag)
impacts(modelo_lag, listw=wkv)
impacts(modelo_lag, tr=trMatc)
impacts(modelo_lag, tr=trMC)

# Simulacion de impactos modelo Durbin Rezago Espacial
# Impactos
summary(modelo_durbin)
impacts(modelo_durbin, listw=wkv)
impacts(modelo_durbin, tr=trMatc)
impacts(modelo_durbin, tr=trMC)

# Simulacion de impactos modelo SARAR
# Impactos 
summary(modelo_sarar)
impacts(modelo_sarar, listw=wkv)
impacts(modelo_sarar, tr=trMatc)
impacts(modelo_sarar, tr=trMC)
#--------------------------------------------------------------------------------

