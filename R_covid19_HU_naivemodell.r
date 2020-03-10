# ver 1.0 -- marc, 10, 2020
# naiv, egyerszeru modell a magyar COVID-19 fertozott betegek szamanak rovidtavu alakulasara

# a modell azért naiv, mert NEM veszi figyelembe 
#   kezdo gocpontok es azok tavolsaganak hatasak
#   a mar valamennyire tisztazott parametereket: onset delay, incubation period (5-8nap), reprodukcios szam (2,5), overdispersion stb.
#   orszagok egymasra hatasat, infrastrukturalis valtozokat es a beavatkozasok hatasat
#   SIR modell korlatait. stb.

# a virusterjedes korai szakaszait exponencialis emelkedes jellemzi
# (de ez az emelkedes is erzekeny a kezdoparameterekre, pl. a gocpontok szamara)

# ez modell egyszeruen linearis regresszioval becslest ad a nepsuruseg alapjan a magyar exponencialis fuggvenyre
# tovabbi valtozok es fejlettebb modszertannal jelentosen javithato

# ez egy elozetes, alacsony megbizhatosagu amator modell, varhatoan professzionalis epidemiologusoktol es matematikusok jobb modelleket epitenek
 

#memoria uritese 
rm(list = ls()) 

#adatok betoltese Johns Hopkins egyetem folyamatos frissitesu adatai
rawcovid <- read.csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_19-covid-Confirmed.csv", header=TRUE, as.is=T, sep=",", fileEncoding = "UTF-8") 

#aggregalas orszagokra 
l2 <- dim(rawcovid)[2]
covid <- aggregate(rawcovid[5:l2], by = list(rawcovid$Country.Region), FUN = sum)
rownames(covid) <- covid[,1]
covid <- covid[,-1]

#Feb 1 ota fertozott EU orszagok kivalasztasa
orszagok <- c("Spain", "Italy", "Germany", "France", "Finland", "Sweden", "UK")
covid_eu <- covid[orszagok,]

#nullak valojaban nem ismertnek tekinthetoek
covid_eu[covid_eu==0]<-NA

# az elso megjelenesre relativizalt vektorok
ES <- na.omit(as.numeric(covid_eu["Spain",]))
IT <- na.omit(as.numeric(covid_eu["Italy",]))
DE <- na.omit(as.numeric(covid_eu["Germany",]))
FR <- na.omit(as.numeric(covid_eu["France",]))
FI <- na.omit(as.numeric(covid_eu["Finland",]))
SE <- na.omit(as.numeric(covid_eu["Sweden",]))
UK <- na.omit(as.numeric(covid_eu["UK",]))

#megjegyzes nem kezelheto data.frame-kent

#lakossagszammal leosztas

ES <- ES / 46549045
IT <- IT / 60391000
DE <- DE / 82576900
FR <- FR / 67187000
FI <- FI / 5516224
SE <- SE / 10142686
UK <- UK / 65648100

HU_pop <- 9771000

# nepsurusegek

ES_d <-  92
IT_d <- 200 
DE_d <- 232  
FR_d <- 118 
FI_d <-  16
SE_d <-  23
UK_d <- 272 

HU_d <- 105

c_names <- c("ES", "IT", "DE", "FR", "FI", "SE", "UK")
c_popd <-  c(92, 200, 232, 118, 16, 23, 272)

#egyszeru linlog modell, ide lehetne valami okosabb nls parameterkereso csodat is -- simple lmlreg, could be optimized
#(bar csodaparameterkeresovel kiszamoltam es nem volt lenyeges kulonbseg -- atmenetileg elfogadhato, de tovabbfejlesztheto)

ES_beta <- coef(summary(lm(log(ES)~c(1:(length(ES))))))[2,1]
ES_R2 <- summary(lm(log(ES)~c(1:(length(ES)))))$adj.r.squared

IT_beta <- coef(summary(lm(log(IT)~c(1:(length(IT))))))[2,1]
IT_R2 <- summary(lm(log(IT)~c(1:(length(IT)))))$adj.r.squared

DE_beta <- coef(summary(lm(log(DE)~c(1:(length(DE))))))[2,1]
DE_R2 <- summary(lm(log(DE)~c(1:(length(DE)))))$adj.r.squared

FR_beta <- coef(summary(lm(log(FR)~c(1:(length(FR))))))[2,1]
FR_R2 <- summary(lm(log(FR)~c(1:(length(FR)))))$adj.r.squared

FI_beta <- coef(summary(lm(log(FI)~c(1:(length(FI))))))[2,1]
FI_R2 <- summary(lm(log(FI)~c(1:(length(FI)))))$adj.r.squared

SE_beta <- coef(summary(lm(log(SE)~c(1:(length(SE))))))[2,1]
SE_R2 <- summary(lm(log(SE)~c(1:(length(SE)))))$adj.r.squared

UK_beta <- coef(summary(lm(log(UK)~c(1:(length(UK))))))[2,1]
UK_R2 <- summary(lm(log(UK)~c(1:(length(UK)))))$adj.r.squared

c_betas <- c(ES_beta, IT_beta, DE_beta, FR_beta, FI_beta, SE_beta, UK_beta)
c_R2s <- c(ES_R2, IT_R2, DE_R2, FR_R2, FI_R2, SE_R2, UK_R2)

# most mar dataframezheto

covid_res <- data.frame(c_names, c_popd, c_betas, c_R2s)

View(covid_res)


# magyar parameterek becslese a c_popd hasznalataval, a becsles kiegeszitheto tovabbi adatokkal, ide otleteket varnank! 
# igy a pre es a poszteszteken is megbukik a modell, de Turkey-szemleletben dolgozva elfogadhato

modellm <- lm(covid_res$c_betas ~ covid_res$c_popd)
coef_lm <- coef(summary(modellm))
res_lm <- predict(modellm) - covid_res$c_betas

#Joslas a HU_beta-ra
HU_beta <- coef_lm[1] + HU_d *  coef_lm[2]
HU_betamax <- HU_beta + max(res_lm)
HU_betamin <- HU_beta - max(res_lm)

# Magyar adatok szimulalasa shossz napra

shossz <- 40

HU_p <- NULL
HU_pmin <- NULL
HU_pmax <- NULL

for (i in 1:shossz) {
HU_p[i] <- exp(1)**(HU_beta*i)
HU_pmin[i] <- exp(1)**(HU_betamin*i)
HU_pmax[i] <- exp(1)**(HU_betamax*i)
}

pnap <- c(1:length(HU_p))

#vizualizacio 
par(mfrow=c(1,2))


#kirajzolas1: nyers adatok

plot(c(1:length(IT)), IT, xlab = "Relatív napok - az országban az első igazolt fertőzés megjelenésétől", ylab = "Igazolt ferőtőzött betegek száma", col=1, lty = 1, log="y", pch = 16, cex=1)
points (ES, col =2, pch = 16, cex=1)
points (UK, col =3, pch = 16, cex=1)
points (SE, col =5, pch = 16, cex=1)
points (DE, col =4, pch = 16, cex=1)
points (FR, col =6, pch = 16, cex=1)
points (FI, col =7, pch = 16, cex=1)
legend("topleft", legend=c("Olaszország" , "Spanyolorszag" , "Egyesült Királyság", "Svédország", "Németország", "Franciaország", "Finnország"), pch=rep(16,7), col = 1:8, inset = 0.05)


#kirajzolas2 : szimulalalt adatok
plot(pnap, HU_pmax, xlab = "Napok, első 40 nap", ylab = "Igazolt fertőzött betegek száma", col=8, lty = 2, pch = 16, cex=0, log="y")
points (HU_pmin, col =8, pch = 16, cex=0)
polygon(c(pnap, rev(pnap)), c(HU_pmax ,rev(HU_pmin)), col = rgb(1, 0, 0,0.5) )
lines(pnap, HU_p, type="b")
legend("topleft", legend="Naiv szimuláció")