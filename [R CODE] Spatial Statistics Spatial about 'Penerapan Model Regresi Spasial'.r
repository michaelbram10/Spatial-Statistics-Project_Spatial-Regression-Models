##### TUGAS 4: Analisis Regresi Spasial #####

### Library ###
library(sf)
library(spdep)
library(ggplot2)
library(lmtest)
library(rgdal)
library(spdep)
library(raster)
library(rgeos)
library(readxl)
library(GWmodel)
library(car)
library(lmtest)
library(maptools)
library(RColorBrewer)
library(nortest)
library(DescTools)
library(spatialreg)

### Import Data ###
library(readxl)
ppm_jabar <- read_excel("D:/Document Bram/Michael Bram UI/SEMESTER 6/Spasial/Tugas/Tugas 4 (Regresi Spasial)/ppm_jabar.xlsx", 
                        sheet = "Data Regresi Spasial")
View(ppm_jabar)
head(ppm_jabar)

### Describe Data ###
str(ppm_jabar)
dim(ppm_jabar)
names(ppm_jabar)

### Define Variabel ###
Y = ppm_jabar$`Persentase Penduduk Miskin`
X1 = ppm_jabar$`Pengeluaran Per Kapita (Ribu Rupiah/Orang/Tahun)`
X2 = ppm_jabar$`Angka Harapan Hidup`

### Import Shapefile ###
Jabar <- readOGR("D:/Document Bram/Michael Bram UI/SEMESTER 6/Spasial/Tugas/Tugas 3 (GWR)/petajawabaratminuswaduk/petajawabaratminuswaduk.shp")
plot(Jabar, main="Peta Jawa Barat")
head(Jabar)

### Uji Moran Y ###
neighbors <- poly2nb(Jabar)
neighbors
neighbors_listw <- nb2listw(neighbors)
neighbors_listw
moran.test(Y, neighbors_listw, randomisation=T, alternative="greater")

plot(neighbors, coordinates(Jabar), col='red', lwd=2, add=TRUE)

### Model Regresi Global ###
model = Y ~ X1+X2

model_lm = lm(model, data=ppm_jabar)
summary(model_lm)
model_lm_res = residuals(model_lm)   #residual model
model_lm_res

### Uji Asumsi Regresi Global ###
##UJI ASUMSI RESIDUAL
#Menampilkan plot
par(mfrow=c(2,2)) 
plot(model_lm)

#1. UJI ASUMSI HETEROSKESDATISITAS (Breusch Pagan Test)
library("quantmod")
library("lmtest")
bp_lm = bptest(model_lm)
bp_lm
spreadLevelPlot(model_lm)

#2. UJI ASUMSI NORMALITAS (Shapiro Wilk Test)
#Cara 1: Uji Shapiro Wilk
norm_lm = shapiro.test(model_lm_res)
norm_lm

#Cara 2: Uji Kolmogorov-Smirnov
ks.test(model_lm$residuals, "pnorm",
        mean=mean(model_lm$residuals), sd=sd(model_lm$residuals))

#Cara 3: Anderson Darling
ad.test(model_lm_res)

#Menampilkan Normal Q-Q Plot
qqnorm(model_lm_res, ylab="Residuals", xlab="Normal Scores")
qqline(model_lm_res)

#3. UJI AUTOKORELASI (Durbin Watson Test)
#Catatan : plot ambil dari Residuals vs Fitted

#Cara 1: Durbin Watson Test
dw_lm = dwtest(model, data = ppm_jabar)
dw_lm

#Cara 2: RunsTest
library(DescTools)
RunsTest(model_lm_res)

### Uji Moran terhadap Residual OLS ###
moran_y <- moran.test(Y, neighbors_listw)
moran_y

moran_X1 <- moran.test(X1, neighbors_listw)
moran_X1

moran_X2 <- moran.test(X2, neighbors_listw)
moran_X2

moran_res <- moran.test(model_lm_res, neighbors_listw)
moran_res

moran.stat <- c(moran_y$statistic, moran_X1$statistic, moran_X2$statistic,
                moran_res$statistic)
moran.stat

moran.pval <- c(moran_y$p.value, moran_X1$p.value, moran_X2$p.value,
                moran_res$p.value)
moran.pval

#Menampilkan statistic dan p-value dari moran terhadap residual ols
moran.table <- data.frame(rbind(unname(moran.stat, force = FALSE), moran.pval))
colnames(moran.table) <- c('Y', 'X1', 'X2', 'Residual')
rownames(moran.table) <- c('statistic', 'p.value')
moran.table

### Langrange Multiplier Test ###
ppm_lagrange <- lm.LMtests(model_lm, neighbors_listw, test=c("LMerr","RLMerr","LMlag","RLMlag","SARMA"))
summary(ppm_lagrange)
#Catatan: yang signifikan hanya LMlag dan SARMA

### SLX ###
model_SLX = lmSLX(model, data=ppm_jabar, neighbors_listw, Durbin = TRUE)
summary(model_SLX)
model_SLX_res = residuals(model_SLX)
model_SLX_res

#1.Asumsi Heteroskedatisitas
bp_SLX = bptest(model_SLX)
bp_SLX
#2. Asumsi Normalitas
norm_SLX = shapiro.test(model_SLX_res)
norm_SLX
#3. Asumsi Autokorelasi
dw_SLX = dwtest(model_SLX, data = ppm_jabar)
dw_SLX

### SAR/SLM (LMlag) (SIGNIFIKAN) ###
model_SLM = lagsarlm(model ,data=ppm_jabar, neighbors_listw)
summary(model_SLM)
model_SLM_res = residuals(model_SLM)
model_SLM_res

#1.Asumsi Heteroskedatisitas
bp_SLM = bptest.Sarlm(model_SLM)
bp_SLM
#2. Asumsi Normalitas
norm_SLM = shapiro.test(model_SLM_res)
norm_SLM
#3. Asumsi Autokorelasi
library(DescTools)
RunsTest(model_SLM_res)

### SEM ### 
model_SEM = errorsarlm(model,data=ppm_jabar,neighbors_listw)
summary(model_SEM)
model_SEM_res = residuals(model_SEM)
model_SEM_res

#1.Asumsi Heteroskedatisitas
bp_SEM = bptest.Sarlm(model_SEM)
bp_SEM
#2. Asumsi Normalitas
norm_SEM = shapiro.test(model_SEM_res)
norm_SEM
#3. Asumsi Autokorelasi
library(DescTools)
RunsTest(model_SEM_res)

### GSM/SARMA (SIGNIFIKAN) ###
model_SARMA = sacsarlm(model, data=ppm_jabar,neighbors_listw)
summary(model_SARMA)
model_SARMA_res = residuals(model_SARMA)
model_SARMA_res

#1.Asumsi Heteroskedatisitas
bp_SARMA = bptest.Sarlm(model_SARMA)
bp_SARMA
#2. Asumsi Normalitas
norm_SARMA = shapiro.test(model_SARMA_res)
norm_SARMA
#3. Asumsi Autokorelasi
library(DescTools)
RunsTest(model_SARMA_res)

### SDM ###
model_SDM = lagsarlm(model,data=ppm_jabar,neighbors_listw)
summary(model_SDM)
model_SDM_res = residuals(model_SDM)
model_SDM_res

#1.Asumsi Heteroskedatisitas
bp_SDM = bptest.Sarlm(model_SDM)
bp_SDM
#2. Asumsi Normalitas
norm_SDM = shapiro.test(model_SDM_res)
norm_SDM
#3. Asumsi Autokorelasi
library(DescTools)
RunsTest(model_SDM_res)

### SDEM ###
model_SDEM = errorsarlm(model, data=ppm_jabar, neighbors_listw, etype="mixed")
summary(model_SDEM)
model_SDEM_res = residuals(model_SDEM)
model_SDEM_res

#1.Asumsi Heteroskedatisitas
bp_SDEM = bptest.Sarlm(model_SDEM)
bp_SDEM
#2. Asumsi Normalitas
norm_SDEM = shapiro.test(model_SDEM_res)
norm_SDEM
#3. Asumsi Autokorelasi
library(DescTools)
RunsTest(model_SARMA_res)

### GNSM ###
model_GNSM = sacsarlm(model,data=ppm_jabar,neighbors_listw,type="sacmixed")
summary(model_GNSM)
model_GNSM_res = residuals(model_GNSM)
model_GNSM_res

#1.Asumsi Heteroskedatisitas
bp_GNSM = bptest.Sarlm(model_GNSM)
bp_GNSM
#2. Asumsi Normalitas
norm_GNSM = shapiro.test(model_GNSM_res)
norm_GNSM
#3. Asumsi Autokorelasi
library(DescTools)
RunsTest(model_GNSM_res)

### AIC ###
AICs = c(AIC(model_lm), AIC(model_SLX), AIC(model_SLM), AIC(model_SEM), 
        AIC(model_SARMA), AIC(model_SDM), AIC(model_SDEM), AIC(model_GNSM))
plot(AICs, type="l", lwd=1.5, xaxt="n", xlab="",col="red")
axis(1, at=1:7,labels=F) #6= number of models
labels<-c("global", "SLX", "SLM","SEM","SARMA", "SDM","SDEM","GNSM")
text(1:8, par("usr")[3]-.25, srt=50, adj=1, labels=labels, xpd=T)
mtext(side=1, text="Model Specification", line=3)

### Tabel Residual ###
model_lm_res = residuals(model_lm)
model_SLM_res = residuals(model_SLM)
model_SARMA_res = residuals(model_SARMA)
#er.ols
df.res.global <- data.frame(ID = c(0:26), model_lm_res)
df.res.SLM <- data.frame(ID = c(0:26), model_SLM_res)
df.res.SARMA <- data.frame(ID = c(0:26), model_SARMA_res)

### Menggabungkan Data Peta2an ###
library(tidyverse)
df.ress <- list(df.res.global, df.res.SLM, df.res.SARMA)
datres <- data.frame(df.ress %>% reduce(full_join, by = 'ID'))
head(datres)

resid <- data.frame(ID = c(0:26), model_lm$residuals)
datres2 <- merge(ppm_jabar[,1:4], datres, by = 'ID')
View(datres2)

#Jadiin Data txt
write.table(datres2, file = "dataresidual.txt")

### Import Data Residual ###
library(readxl)
dataresidual <- read_excel("D:/Document Bram/Michael Bram UI/SEMESTER 6/Spasial/Tugas/Tugas 4 (Regresi Spasial)/dataresidual.xlsx", 
                           sheet = "dataresidual")
View(dataresidual) 

### Merubah variabel pada data Shapefile menjadi variabel numerik ###
Jabar@data$ID <- as.numeric(row.names(Jabar@data))
Jabar@data$ID
Jabar@data$row <- as.numeric(row.names(Jabar@data))
Jabar@data$row 
str(Jabar@data)
length(Jabar@data$ADM2_EN)
length(Jabar@data$ID)

### Menggabungkan Shapefile dengan Data di Excel ###
temp <- merge(Jabar@data, dataresidual, by="ID", all.Jabar=T, sort=F)
temp
Jabar@data
dim(temp)
Jabar@data <- temp[order(temp$row),]
Jabar@data$ADM2_EN
str(Jabar@data)

### Layout ###
Data <- read_excel("D:/Document Bram/Michael Bram UI/SEMESTER 6/Spasial/Tugas/Tugas 4 (Regresi Spasial)/ppm_jabar.xlsx", 
                   sheet = "Data Regresi Spasial")
Data.spdf<-SpatialPointsDataFrame(temp[,19:18], temp) #3 dan 4 adalah kolom Lintang Bujur
head(Data.spdf)
kabkot <- list('sp.pointLabel', Data.spdf, label=Jabar@data$ADM2_EN, cex=0.7)
titik <- list('sp.points', Data.spdf, pch=16, cex=.8, col = "blue")

#Plot Peta Residual Regresi Global
spplot(Jabar, 'model_lm_res', col.regions = brewer.pal(7,"RdBu"), 
       main="Plot Peta Residual Regresi Global", 
       scales = list(draw = TRUE), 
       sp.layout=list(kabkot, titik), 
       xlab="Longitude", 
       ylab="Latitude", 
       at=seq(-max(abs(Jabar@data$model_lm_res))-0.001,
              max(abs(Jabar@data$model_lm_res))+0.001,
              2*max(abs(Jabar@data$model_lm_res))/7))  

#Plot Peta Residual Model SAR/LMlag
spplot(Jabar, 'model_SLM_res', col.regions = brewer.pal(7,"PuOr"), 
       main="Plot Peta Residual Model SAR", 
       scales = list(draw = TRUE), 
       sp.layout=list(kabkot, titik), 
       xlab="Longitude", 
       ylab="Latitude",
       at=seq(-max(abs(Jabar@data$model_SLM_res))-0.001,
              max(abs(Jabar@data$model_SLM_res))+0.001,
              2*max(abs(Jabar@data$model_SLM_res))/7))

#Plot Peta Residual Model SARMA
spplot(Jabar, 'model_SARMA_res', col.regions = brewer.pal(7,"PiYG"), 
       main="Plot Peta Residual Model SARMA", 
       scales = list(draw = TRUE), 
       sp.layout=list(kabkot, titik), 
       xlab="Longitude", 
       ylab="Latitude",
       at=seq(-max(abs(Jabar@data$model_SARMA_res))-0.001,
              max(abs(Jabar@data$model_SARMA_res))+0.001,
              2*max(abs(Jabar@data$model_SARMA_res))/7))