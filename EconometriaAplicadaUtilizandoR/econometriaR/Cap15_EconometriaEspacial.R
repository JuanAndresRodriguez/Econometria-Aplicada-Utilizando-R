#Econometria Espacial, 2015
#Profesor: Miguel Angel Mendoza, FE-UNAM 
#Estimaciones de Modelos Espaciales 

#Cargar la libreria spdep

library(spdep) # Econometria espacial 
library(maptools) # leer archivos shapfiles y elaborar mapas 
library(RColorBrewer) # Eleccion de colores
library(classInt) # m??todos para clasificar


#Cambiar el directorio de trabajo
setwd("/Volumes/LACIE SHARE/Academico/LibroEconometria_R/Capitulo_14/BaseDatos_Capitulo14")
#Leer archivos shapes y transformarlo en objeto Shape y DataFrame 
empleo <- readShapePoly("Zona_Centro.shp")
summary(empleo)

# Logaritmo del Empleo 
lempleo <- log(empleo$POCUPADA)
# Capital Humano y logaritmo del capital humano, a??os de escolaridad
ch <- empleo$ESCOLA_15
lch <- log(empleo$ESCOLA_15)

# Caracteristicas estadisticas del empleo y capital humano
summary(empleo$POCUPADA)
summary(empleo$ESCOLA_15)
summary(lempleo)
summary(ch)
summary(lch)

# Construir lista de vecinos tipo Queen de poligonos 
pr.nb <- poly2nb(empleo, queen=TRUE)
# Matriz de ponderacion W estandarizada
wqueen <- nb2listw(pr.nb, style="W")

# Caracteristicas de la Matriz W tipo Queen
summary(wqueen)

# Grafica con la conexicion espacial  
cent <- coordinates(empleo)
plot(empleo, border="grey", lwd=1.5)
plot(pr.nb,cent, add=T, col="darkred")


# Estadistico de Moran 
moran_lempleo <- moran.test(lempleo, wqueen,randomisation=TRUE, alternative="two.sided", na.action=na.exclude)
moran_ch <- moran.test(ch, wqueen,randomisation=TRUE, alternative="two.sided", na.action=na.exclude)
moran_lch <- moran.test(lch, wqueen,randomisation=TRUE, alternative="two.sided", na.action=na.exclude)

#Ver resultados
print(moran_lempleo)
print(moran_ch)
print(moran_lch)

# Grafica de diagrama de dispersion de Moran 
moran.plot(lempleo, wqueen, pch=20)
moran.plot(ch, wqueen, pch=20)
moran.plot(lch, wqueen, pch=20)

#Clasificador con cuartiles y colores especificos
# Mapa de cuartile del logaritmo del empleo
brks <- round(quantile(lempleo, probs=seq(0,1,0.25)), digits=2)
colours <- brewer.pal(4,"Reds")

plot(empleo, col=colours[findInterval(lempleo, brks, all.inside=TRUE)],
     axes=F)
legend(x=-87.9, y=25.2, legend=leglabs(brks), fill=colours, bty="n")
invisible(title(main=paste("EMPLEO", sep="\n")))
box()

# Mapa de cuartiles del logaritmo del capital humano
brks <- round(quantile(lch, probs=seq(0,1,0.25)), digits=2)
colours <- brewer.pal(4,"Blues")

plot(empleo, col=colours[findInterval(lch, brks, all.inside=TRUE)],
     axes=F)
legend(x=-87.9, y=25.2, legend=leglabs(brks), fill=colours, bty="n")
invisible(title(main=paste("CAPITAL HUMANO", sep="\n")))
box()

# An??lisis LISA

# Valores de referencia z de la distribuci??n t
z <- c(1.65, 1.96)
zc <- c(2.8284, 3.0471)

# Estimaci??n de indice de Moran local (Ii)
f.Ii <- localmoran(lempleo, wqueen)
zIi <- f.Ii[,"Z.Ii"] # Asignaci??n de la distribuci??n Z del Ii
mx <- max(zIi) 
mn <- min(zIi)

# Mapa de signficancia para los z-scores
pal <- c("white", "green", "darkgreen")
z3.Ii <- classIntervals(zIi, n=3, style="fixed", fixedBreaks=c(mn, z, mx))
cols.Ii <- findColours(z3.Ii, pal)  
plot(empleo, col=cols.Ii)
brks <- round(z3.Ii$brks,4)
leg <- paste(brks[-4], brks[-1], sep=" - ")
legend(x=-87.9, y=25.2, fill=pal, legend=leg, bty="n")         

# Mapa de los grupos de cluster
pal.rb <- c("skyblue1","white","blue","red")
z4.Ii <- classIntervals(zIi, n=4, style="fixed", fixedBreaks=c(min(f.Ii), -z[1], z[1], z[2], max(f.Ii)))
cols.Ii <- findColours(z4.Ii, pal.rb)               
plot(empleo, col=cols.Ii)
brks <- round(z4.Ii$brks,4)
leg <- paste(brks[-5], brks[-1], sep=" - ")
legend(x=-87.9, y=25.2, fill= pal.rb, legend=leg, bty="n")


# Modelo OLS 
ModeloEmpleo_OLS <- lm(lempleo ~ lch , data=empleo)
summary(ModeloEmpleo_OLS)

# Prueba de Moran a residuales del modelo OLS
I_Moran <- lm.morantest(ModeloEmpleo_OLS,wqueen)
print(I_Moran)

#Pruebas de Multiplicadores de Lagranges lm.LMtests(columbus.lm,col.listw,test=c("LMerr","RLMerr","LMlag","RLMlag","SARMA"))
lm.LMtests(ModeloEmpleo_OLS,wqueen,test=c("LMerr","RLMerr","LMlag","RLMlag","SARMA"))

# Modelos Espaciales 
# Estimar el  Modelo Rezago Espacial 
ModeloEmpleo_lag <- lagsarlm(lempleo ~ lch , data=empleo,wqueen)
summary(ModeloEmpleo_lag)

# Estimar el modelo de Error Espacial 
ModeloEmpleo_err <- errorsarlm(lempleo ~ lch , data=empleo,wqueen)
summary(ModeloEmpleo_err)

# Estimar modelo SARAR
ModeloEmpleo_sarar <- sacsarlm(lempleo ~ lch , data=empleo, wqueen, type="sac")
summary(ModeloEmpleo_sarar)

#Estimar el modelo de Durbin Rezago Espacial
ModeloEmpleo_lag_durbin <- lagsarlm(lempleo ~ lch , data=empleo,wqueen, type="mixed")
summary(ModeloEmpleo_lag_durbin)
