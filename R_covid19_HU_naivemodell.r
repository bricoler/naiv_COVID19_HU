# ver 1.0 -- marc, 11, 2020
# naiv, egyerszeru modell a magyar COVID-19 fertozott betegek szamanak rovidtavu alakulasara

# a modell azért naiv, mert NEM veszi figyelembe 
#   kezdo gocpontok es azok tavolsaganak hatasak
#   a mar valamennyire tisztazott parametereket: onset delay, incubation period (5-8nap), reprodukcios szam (2,5), overdispersion stb.
#   orszagok egymasra hatasat, infrastrukturalis valtozokat es a beavatkozasok hatasat
#   SIR modell korlatait. stb.
#   valamint determinisztikus.

#memoria uritese 
rm(list = ls()) 

#adatok betoltese Johns Hopkins egyetem folyamatos frissitesu online adataibol
rawcovid <- read.csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_19-covid-Confirmed.csv", header=TRUE, as.is=T, sep=",", fileEncoding = "UTF-8") 

#aggregalas orszagokra 
l2 <- dim(rawcovid)[2]
covid <- aggregate(rawcovid[5:l2], by = list(rawcovid$Country.Region), FUN = sum)
rownames(covid) <- covid[,1]
covid <- covid[,-1]

#Jan 30 ota fertozott EU orszagok kivalasztasa
orszagok <- c("Spain", "Italy", "Germany", "France", "Finland", "Sweden", "UK", "Belgium")
covid_eu <- covid[orszagok,]

#nullak valojaban nem ismertnek tekinthetoek
covid_eu[covid_eu==0]<-NA

# az elso megjelenesre relativizalt vektorok
ES1 <- na.omit(as.numeric(covid_eu["Spain",]))
IT1 <- na.omit(as.numeric(covid_eu["Italy",]))
DE1 <- na.omit(as.numeric(covid_eu["Germany",]))
FR1 <- na.omit(as.numeric(covid_eu["France",]))
FI1 <- na.omit(as.numeric(covid_eu["Finland",]))
SE1 <- na.omit(as.numeric(covid_eu["Sweden",]))
UK1 <- na.omit(as.numeric(covid_eu["UK",]))
BE1 <- na.omit(as.numeric(covid_eu["Belgium",]))

# kb. 17 felett lobban be a virusterjedes
# vektorok, amelyek csak a 17 feletti szakasz mutatjak

ES2 <- ES1[which(ES1>17)]
IT2 <- IT1[which(IT1>17)]
DE2 <- DE1[which(DE1>17)]
FR2 <- FR1[which(FR1>17)]
FI2 <- FI1[which(FI1>17)]
SE2 <- SE1[which(SE1>17)]
UK2 <- UK1[which(UK1>17)]
BE2 <- BE1[which(BE1>17)]

#nem kezelheto data.frame-kent, mert eltero hosszauk

#lakossagszammal leosztas, igy a fertozottek lakosagaranya jon ki

ES <- ES2 / 46549045
IT <- IT2 / 60391000
DE <- DE2 / 82576900
FR <- FR2 / 67187000
FI <- FI2 / 5516224
SE <- SE2 / 10142686
UK <- UK2 / 65648100
BE <- BE2 / 11574255

HU_pop <- 9771000

# Haztartasi mutato
# Average number of persons per household by household composition
# https://ec.europa.eu/eurostat/web/products-datasets/product?code=lfst_hhanwhtc
# Nem bizonyult 95% CI mellett sziginfikiansnak, de azert bennehagytam

ES_d <-  2.5
IT_d <- 2.3 
DE_d <- 2 
FR_d <- 2.2
FI_d <- 2.1
SE_d <- 1.8
UK_d <- 2.3
BE_d <- 2.3

HU_d <- 2.3

# Number of big International airports 
# Nagy (eves forgalom 1 millio felett)
# Forras: https://en.wikipedia.org/wiki/List_of_international_airports_by_country
# Terjedes elso szakaszaban kiemelt szerepe van a reptereknek
# https://www.pnas.org/content/103/7/2015.short
# Ez bizonyult szignifikansnak 95% CI mellett

ES_a <-  17
IT_a <- 22
DE_a <- 17
FR_a <- 15
FI_a <-  1 
SE_a <-  7 
UK_a <- 13 
BE_a <- 2 

HU_a <- 1

c_names <- c("ES", "IT", "DE", "FR", "FI", "SW", "UK", "BE") # orszagnevek
c_popd <-  c(ES_a, IT_a, DE_a, FR_a, FI_a, SE_a, UK_a, BE_a) # nagyrepterek
c_np <- c(ES_d, IT_d, DE_d, FR_d, FI_d, SE_d, UK_d, BE_d) # haztartasi zsufoltsag


#egyszeru linlog modell, repter es felfutas a 17. nap utan

ES_C <- coef(summary(lm(log(ES)~c(1:(length(ES))))))[1,1]
ES_beta <- coef(summary(lm(log(ES)~c(1:(length(ES))))))[2,1]
ES_R2 <- summary(lm(log(ES)~c(1:(length(ES)))))$adj.r.squared

IT_C <- coef(summary(lm(log(IT)~c(1:(length(IT))))))[1,1]
IT_beta <- coef(summary(lm(log(IT)~c(1:(length(IT))))))[2,1]
IT_R2 <- summary(lm(log(IT)~c(1:(length(IT)))))$adj.r.squared

DE_C <- coef(summary(lm(log(DE)~c(1:(length(DE))))))[1,1]
DE_beta <- coef(summary(lm(log(DE)~c(1:(length(DE))))))[2,1]
DE_R2 <- summary(lm(log(DE)~c(1:(length(DE)))))$adj.r.squared

FR_C <- coef(summary(lm(log(FR)~c(1:(length(FR))))))[1,1]
FR_beta <- coef(summary(lm(log(FR)~c(1:(length(FR))))))[2,1]
FR_R2 <- summary(lm(log(FR)~c(1:(length(FR)))))$adj.r.squared

FI_C <- coef(summary(lm(log(FI)~c(1:(length(FI))))))[1,1]
FI_beta <- coef(summary(lm(log(FI)~c(1:(length(FI))))))[2,1]
FI_R2 <- summary(lm(log(FI)~c(1:(length(FI)))))$adj.r.squared

SE_C <- coef(summary(lm(log(SE)~c(1:(length(SE))))))[1,1]
SE_beta <- coef(summary(lm(log(SE)~c(1:(length(SE))))))[2,1]
SE_R2 <- summary(lm(log(SE)~c(1:(length(SE)))))$adj.r.squared

UK_C <- coef(summary(lm(log(UK)~c(1:(length(UK))))))[1,1]
UK_beta <- coef(summary(lm(log(UK)~c(1:(length(UK))))))[2,1]
UK_R2 <- summary(lm(log(UK)~c(1:(length(UK)))))$adj.r.squared

BE_C <- coef(summary(lm(log(BE)~c(1:(length(BE))))))[1,1]
BE_beta <- coef(summary(lm(log(UK)~c(1:(length(UK))))))[2,1]
BE_R2 <- summary(lm(log(UK)~c(1:(length(UK)))))$adj.r.squared

c_Cs <- c(ES_C, IT_C, DE_C, FR_C, FI_C, SE_C, UK_C, BE_C)
c_betas <- c(ES_beta, IT_beta, DE_beta, FR_beta, FI_beta, SE_beta, UK_beta, BE_beta)
c_R2s <- c(ES_R2, IT_R2, DE_R2, FR_R2, FI_R2, SE_R2, UK_R2, BE_R2)

# most mar dataframezheto

covid_res <- data.frame(c_names, c_Cs, c_popd, c_np, c_betas, c_R2s)

View(covid_res)


#repter es beta 

modellm_beta <- lm(covid_res$c_betas ~ covid_res$c_popd) 
coef_lm_beta <- coef(summary(modellm_beta))
res_lm_beta <- predict(modellm_beta) - covid_res$c_betas
pred_beta <- predict(modellm_beta, HU_a, interval= "confidence" , level=0.9)

HU_beta <- coef_lm_beta[1] + HU_d *  coef_lm_beta[2]
HU_betamax <- HU_beta + max(res_lm_beta) # ez kb. 90% CI 
HU_betamin <- HU_beta - max(res_lm_beta) # ez kb. 90% CI

#Joslas a HU_C-re
modellm_c <- lm(covid_res$c_Cs ~ covid_res$c_popd)
coef_lm_c <- coef(summary(modellm_c))
res_lm_c <- predict(modellm_c) - covid_res$c_Cs

HU_C <- coef_lm_c[1] + HU_d *  coef_lm_c[2]
HU_Cmax <- HU_C + max(res_lm_c) # ez kb. 90% CI 
HU_Cmin <- HU_C - max(res_lm_c) # ez kb. 90% CI 

summary(modellm_beta) 
summary(modellm_c) # 


#mivel a p ertek gyenge, ezert becsles elvetesebecsles
# becsles a HU_beta-ra sima atlag es szoras hasznalataval

shapiro.test(covid_res$c_Cs)
shapiro.test(covid_res$c_betas)

HU_beta <- mean(covid_res$c_betas)
error_beta <- qnorm(0.975)*sd(covid_res$c_betas)/sqrt(8)

HU_betamax <- HU_beta + error_beta
HU_betamin <- HU_beta - error_beta

HU_C <- mean(covid_res$c_Cs)
error_c <- qnorm(0.975)*sd(covid_res$c_Cs)/sqrt(8)

HU_Cmax <- HU_C + error_c
HU_Cmin <- HU_C - error_c


# hazai adatok szimulalasa shossz napra

shossz <- 40

HU_p <- NULL
HU_pmin <- NULL
HU_pmax <- NULL

# e kitevojeben a josolt beta 

for (i in 1:shossz) {
HU_p[i] <- exp(1)**(HU_C+(HU_beta*i))
HU_pmin[i] <- exp(1)**(HU_Cmin+(HU_betamin*i))
HU_pmax[i] <- exp(1)**(HU_Cmax+(HU_betamax*i))
}

# A lakossanyaranyos becslecst fel kell szorozni a lakossaggal
# elivleg itt meg az R2 is josolni lehett volna es annak 
# a bizonytalansaga is be lehetett volna epiteni
# ez Finnorszag eseteben 0,6, amire szinten szamolhato CI 90% 0.4-0.84

HU_p <- HU_p*HU_pop
HU_pmin <- (HU_pmin*HU_pop) * mean(covid_res$c_R2s)
HU_pmax <- (HU_pmax*HU_pop) * (1 + (1-mean(covid_res$c_R2s)))



#Kerekites

HU_pi <- ceiling(HU_p)
HU_pmini <- ceiling(HU_pmin) 
HU_pmaxi <- ceiling(HU_pmax) 

pnap <- length(HU_pmini)

#vizualizacio 
par(mfrow=c(1,1))


#kirajzolas1: nyers adatok
plot(c(1:length(FR1)), FR1, xlab = "Eltelt napok száma (az országban az első igazolt fertőzés megjelenésétől)", ylab = "Igazolt ferőtőzött betegek száma (log)", col=1, lty = 1, log="y", pch = 16, cex=1)
points (ES1, col =2, pch = 16, cex=1)
points (UK1, col =3, pch = 16, cex=1)
points (SE1, col =5, pch = 16, cex=1)
points (DE1, col =4, pch = 16, cex=1)
points (IT1, col =6, pch = 16, cex=1)
points (BE1, col =7, pch = 16, cex=1)
points (FI1, col =8, pch = 16, cex=1)
legend("topleft", legend=c("Olaszország" , "Spanyolorszag" , "Egyesült Királyság", "Svédország", "Németország", "Franciaország", "Belgium", "Finnország"), pch=rep(16,8), col = c(6,2,3,5,4,1,8), y.intersp = 0.7, inset = 0.01, cex=1, box.lwd = 0,box.col = "white",bg = "white")


#kirajzolas2 : szimulalalt adatok
#hany napot rajzoljon ki : kir
kir <- 17
pnap <- c(1:kir)

plot(pnap, HU_pmaxi[1:kir], xlab = "Napok, 17. beteg megjelenése utáni értékek", ylab = "Igazolt fertőzött betegek száma, becslés (CI~0.92)", col=8, lty = 2, pch = 16, cex=0, ylim=c(1,max(HU_pmaxi[1:kir])))
points (HU_pmini[1:kir], col =8, pch = 16, cex=0)
polygon(c(pnap[1:kir], rev(pnap[1:kir])), c(HU_pmaxi[1:kir] ,rev(HU_pmini[1:kir])), col = rgb(1, 0, 0,0.5) )
points(pnap, HU_pi[1:kir], col=1, lty = 2, pch = 16, cex=1, type="b")

#SIR kod egy resze nem az en szellemi termekem, ezert nem oszthatom meg nyilvanosan, de ha vki ir emailt, kuldom
