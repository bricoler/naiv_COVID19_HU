# ver 1.0 -- marc, 10, 2020
# naiv, egyerszeru modell a magyar COVID-19 fertozott betegek szamanak rovidtavu alakulasara
# a nagy repterek szamanak bevonasaval

# a modell azért naiv, mert NEM veszi figyelembe 
#   kezdo gocpontok es azok tavolsaganak hatasak
#   a mar valamennyire tisztazott parametereket: onset delay, incubation period (5-8nap), reprodukcios szam (2,5), overdispersion stb.
#   orszagok egymasra hatasat, infrastrukturalis valtozokat es a beavatkozasok hatasat
#   SIR modell korlatait. stb.

# a virusterjedes korai szakaszait exponencialis emelkedes jellemzi
# (de ez az emelkedes is erzekeny a kezdoparameterekre, pl. a gocpontok szamara)

# ez modell egyszeruen linearis regresszioval becslest ad 
# a nagy repterek szama alapjan a magyar exponencialis fuggvenyre
# tovabbi valtozok es fejlettebb modszertannal javithato
# gyakorlatilag a Finorszagi terjedest, de a magyar lakossagszamra

# ez egy elozetes, alacsony megbizhatosagu amator, naiv modell
 

#memoria uritese 
rm(list = ls()) 

#adatok betoltese Johns Hopkins egyetem folyamatos frissitesu online adataibol
rawcovid <- read.csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_19-covid-Confirmed.csv", header=TRUE, as.is=T, sep=",", fileEncoding = "UTF-8") 

#aggregalas orszagokra 
l2 <- dim(rawcovid)[2]
covid <- aggregate(rawcovid[5:l2], by = list(rawcovid$Country.Region), FUN = sum)
rownames(covid) <- covid[,1]
covid <- covid[,-1]

#Feb 1 ota fertozott EU orszagok kivalasztasa
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

#nem kezelheto data.frame-kent, mert eltero hosszauk

#lakossagszammal leosztas, igy a fertozottek lakosagaranya jon ki

ES <- ES1 / 46549045
IT <- IT1 / 60391000
DE <- DE1 / 82576900
FR <- FR1 / 67187000
FI <- FI1 / 5516224
SE <- SE1 / 10142686
UK <- UK1 / 65648100
BE <- BE1 / 11574255

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

c_names <- c("ES", "IT", "DE", "FR", "FI", "SE", "UK", "BE") # orszagnevek
c_popd <-  c(ES_a, IT_a, DE_a, FR_a, FI_a, SE_a, UK_a, BE_a) # nagyrepterek
c_np <- c(ES_d, IT_d, DE_d, FR_d, FI_d, SE_d, UK_d, BE_d) # haztartasi zsufoltsag


#egyszeru linlog modell, ide lehetne valami okosabb nls parameterkereso csodat is  
#(bar csodaparameterkeresovel kiszamoltam es nem volt lenyeges kulonbseg )

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

# magyar parameterek becslese a nagy repterek hasznalataval
# egyes pre es a poszteszteken is megbukik a modell, de Turkey-szemleletben dolgozva elfogadhato


#Joslas a HU_beta-ra

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

summary(modellm_beta) # p ertek elfogadhato
summary(modellm_c) # p ertek elfogadhato

# Magyar adatok szimulalasa shossz napra, 
# e hatvanyara kell emelni, a kitevore tettunk joslast a linreggel

shossz <- 365

HU_p <- NULL
HU_pmin <- NULL
HU_pmax <- NULL

for (i in 1:shossz) {
HU_p[i] <- exp(1)**(HU_C+(HU_beta*i))
HU_pmin[i] <- exp(1)**(HU_Cmin+(HU_betamin*i))
HU_pmax[i] <- exp(1)**(HU_Cmax+(HU_betamax*i))
}

#A lakossanyaranyos becslecst fel kell szorozni a lakossaggal
# elivleg itt meg az R2 is josolni lehett volna es annak 
# a bizonytalansaga is be lehetett volna epiteni
# ez Finnorszag eseteben 0,6, amire szinten szamolhato CI 90% 0.4-0.84

HU_p <- HU_p*HU_pop
HU_pmin <- (HU_pmin*HU_pop) * mean(covid_res$c_R2s)
HU_pmax <- (HU_pmax*HU_pop) * (1 + (1-mean(covid_res$c_R2s)))



# A harom joslast egymashoz illeszteni
# hany fertozottszamtol kezdje az illesztett vektort, mf

mf <- 11

HU_pi <- ceiling(HU_p[which(HU_p>mf)])
HU_pmini <- ceiling(HU_pmin[which(HU_pmin>mf)]) 
HU_pmaxi <- ceiling(HU_pmax[which(HU_pmax>mf)]) 

pnap <- length(HU_pmini)

#vizualizacio 
par(mfrow=c(1,1))


#kirajzolas1: nyers adatok

plot(c(1:length(IT)), IT, xlab = "Eltelt napok száma (az országban az első igazolt fertőzés megjelenésétől)", ylab = "Igazolt ferőtőzött betegek száma (log)", col=1, lty = 1, log="y", pch = 16, cex=1)
points (ES, col =2, pch = 16, cex=1)
points (UK, col =3, pch = 16, cex=1)
points (SE, col =5, pch = 16, cex=1)
points (DE, col =4, pch = 16, cex=1)
points (FR, col =6, pch = 16, cex=1)
points (BE, col =7, pch = 16, cex=1)
points (FI, col =8, pch = 16, cex=1)
legend("topleft", legend=c("Olaszország" , "Spanyolorszag" , "Egyesült Királyság", "Svédország", "Németország", "Franciaország", "Belgium", "Finnország"), pch=rep(16,8), col = 1:8, y.intersp = 0.7, inset = 0.01, cex=1, box.lwd = 0,box.col = "white",bg = "white")


#kirajzolas2 : szimulalalt adatok
#hany napot rajzoljon ki : kir
kir <- 30
pnap <- c(1:kir)

plot(pnap, HU_pmaxi[1:kir], xlab = "Napok, 2020.03.10 - 2020.04.10", ylab = "Igazolt fertőzött betegek száma, becslés (qCI~0.72)", col=8, lty = 2, pch = 16, cex=0, log="y")
points (HU_pmini[1:kir], col =8, pch = 16, cex=0)
polygon(c(pnap[1:kir], rev(pnap[1:kir])), c(HU_pmaxi[1:kir] ,rev(HU_pmini[1:kir])), col = rgb(1, 0, 0,0.5) )
points(pnap, HU_pi[1:kir], col=1, lty = 2, pch = 16, cex=1, type="b")

##SIR
## Nem az en szellemi termekem volt az alapkod, ezert nem tudtam felrakni. Emailos megkereses eseten szivesen elkuldom.
