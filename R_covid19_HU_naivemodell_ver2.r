# ver 2.0 -- marc, 11, 2020
# naiv, egyerszeru modell a magyar COVID-19 fertozott betegek szamanak rovidtavu alakulasara

# a modell azért naiv, mert NEM veszi figyelembe 
#   kezdo gocpontok es azok tavolsaganak hatasak
#   a mar valamennyire tisztazott parametereket: onset delay, incubation period (5-8nap), reprodukcios szam (2,5), overdispersion stb.
#   orszagok egymasra hatasat, infrastrukturalis valtozokat es a beavatkozasok hatasat
#   SIR modell korlatait. stb.
#   valamint determinisztikus.

#memoria uritese 
rm(list = ls()) 

# Felepites:
# Betolti az adatokat online
# Aggregal orszagokra
# Minden orszagbol kidobja min_fert erteknel kisebb szakaszt
# Kiszamol minden orszagra exp es log fuggvenyertekeket
# Exponencialis EU-sokbol josol HU-ra 
# Kirajzolja

min_fert <- 17
min_hossz <- 14 

# ha min_hossznal nagyobb a vektor akkor kiszamol exponencialis es logaritmikus illesztes
# minden orszagra visszadobja az adatokat


#BETOLTES
#adatok betoltese Johns Hopkins egyetem folyamatos frissitesu online adataibol
rawcovid <- read.csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_19-covid-Confirmed.csv", header=TRUE, as.is=T, sep=",", fileEncoding = "UTF-8") 

#AGGREGALAS
#aggregalas orszagokra 
# itt eldobjuk a GEO-t, de a jovoben lehetne kezdeni valamit a tavolsagi adatokkal, ketdimenzioban (ido, ter)
l2 <- dim(rawcovid)[2] # napok hossza

covid <- aggregate(rawcovid[5:l2], by = list(rawcovid$Country.Region), FUN = sum) # orszagok szerinti aggregalas
rownames(covid) <- covid[,1] #sornevek javitasa
covid <- covid[,-1] # elso oszlop torles


#dif matrix eloallitasa (delta)
d_covid <- covid[1:dim(covid)[1],2:dim(covid)[2]] - covid[1:dim(covid)[1],(1:dim(covid)[2]-1)]


#SZAMITAS
# kulonbozo orszagok metaadatok szamitas, exp es log fuggvenyillesztes
# hianyzik: orszaglakossag es lakossagaranyos emelkedesszamitas

covid[covid==0]<-NA #0-kbol NA 

cname <- NULL
chossz <-NULL
cmin <- NULL
cmax <- NULL
cmean <- NULL
cmedian <-NULL
csd <- NULL

tempv <-NULL
tempv2 <- NULL

c_expC <- NULL
c_expbeta <- NULL
c_expR2 <- NULL

c_logC <- NULL
c_logbeta <- NULL
c_logR2 <- NULL

c_tend <- NULL

for (i in 1:dim(covid)[1]) 
{
cname[i] <- rownames(covid)[i] # orszag neve
chossz[i] <- sum(as.numeric(covid[i,])>0, na.rm=T) # hany nap van tobb mint 0 fertozott
cmin[i] <- min(as.numeric(covid[i,]), na.rm=T) # minimum fertozottszam
cmax[i] <- max(as.numeric(covid[i,]), na.rm=T) # maximum fertozottszam
cmean[i] <- mean(as.numeric(covid[i,]), na.rm=T) # atlagos fertozottszam
cmedian[i] <- median(as.numeric(covid[i,]), na.rm=T) #median fertozottszam
csd[i] <- sd(as.numeric(covid[i,]), na.rm=T) #sd fertozottszam

tempv <- NULL 
tempv2 <- NULL

tempv <- na.omit(as.numeric(covid[i,])) # atmeneti vektor a nullak kidobalasa
tempv2 <- tempv[which(tempv>min_fert)] # 17-nel kisebbek kidobalasa

#ha 10-nel kevesebb a megfigyeles akkor nem szamolunk
#utana exp illesztes
#utana log illesztes

if (length(tempv2)>min_hossz) {
#exp
c_expC[i] <- coef(summary(lm(log(tempv2)~c(1:(length(tempv2))))))[1,1]
c_expbeta[i] <- coef(summary(lm(log(tempv2)~c(1:(length(tempv2))))))[2,1]
c_expR2[i] <- summary(lm(log(tempv2)~c(1:(length(tempv2)))))$adj.r.squared 
#log
c_logC[i] <- coef(summary(lm((tempv2)~log(c(1:(length(tempv2)))))))[1,1]
c_logbeta[i] <- coef(summary(lm((tempv2)~log(c(1:(length(tempv2)))))))[2,1]
c_logR2[i] <- summary(lm((tempv2)~log(c(1:(length(tempv2))))))$adj.r.squared 


if (c_logR2[i]<c_expR2[i]) {c_tend[i] <- "exp"} else {c_tend[i] <- "log"}

} 
else  
{
c_expC[i] <- 0
c_expbeta[i] <- 0
c_expR2[i] <- 0 
c_logC[i] <- 0
c_logbeta[i] <- 0
c_logR2[i] <- 0 
c_tend[i]<- 0
}

}

cmeta <- data.frame(cname,chossz,cmin,cmax,cmean,cmedian,csd, c_expC, c_expbeta, c_expR2, c_logC, c_logbeta, c_logR2, c_tend) 

#cmeta, orszagokrol fertozesmetaadatok

View(cmeta)

# exp es a log orszagok szetvalogatasa ket kulon data.frame-be.

cmeta_exp <- cmeta[as.numeric(which(cmeta$c_tend=="exp")),]
cmeta_log <- cmeta[as.numeric(which(cmeta$c_tend=="log")),]

# Magyarorszag vonzaskorezetebe talalhato orszagok

europa <- c("Italy", "Germany", "France", "Finland", "Sweden", "United Kingdom", "Denmark", "Austria", "Greece", "Iceland", "Ireland", "Norway", "Portugal", "Switzerland")

cmeta_exp_eu <- cmeta_exp[(cmeta_exp$cname %in% europa),]

#tesztek exponencialis parameterek eloszlasara

shapiro.test(cmeta_exp_eu$c_expC)[2]
shapiro.test(cmeta_exp_eu$c_expbeta)[2]
shapiro.test(cmeta_exp_eu$c_expR2)[2]


#exp_beta
HU_beta <- mean(cmeta_exp_eu$c_expbeta)
error_beta <- qnorm(0.975)*sd(cmeta_exp_eu$c_expbeta)/sqrt(15)

HU_betamax <- HU_beta + error_beta
HU_betamin <- HU_beta - error_beta

#exp_C
HU_C <- mean(cmeta_exp_eu$c_expC)
error_c <- qnorm(0.975)*sd(cmeta_exp_eu$c_expC)/sqrt(15)

HU_Cmax <- HU_C + error_c
HU_Cmin <- HU_C - error_c

#expR2
HU_R2 <- mean(cmeta_exp_eu$c_expR2)
error_R2 <- qnorm(0.975)*sd(cmeta_exp_eu$c_expR2)/sqrt(15)

HU_R2max <- HU_R2 + error_R2
HU_R2min <- HU_R2 - error_R2


# hazai adatok szimulalasa shossz napra

shossz <- 14

HU_p <- NULL
HU_pmin <- NULL
HU_pmax <- NULL

# e kitevojeben a josolt beta 

for (i in 1:shossz) {
HU_p[i] <- exp(1)**(HU_C+(HU_beta*i))
HU_pmin[i] <- exp(1)**(HU_Cmin+(HU_betamin*i))
HU_pmax[i] <- exp(1)**(HU_Cmax+(HU_betamax*i))
}

#R2 bizonytalansag fapados hozzaszamitasa
HU_pmin <- (HU_pmin) * HU_R2
HU_pmax <- (HU_pmax) * (1+ (1-HU_R2))


#Kerekites egesz emberre
HU_pi <- ceiling(HU_p)
HU_pmini <- ceiling(HU_pmin) 
HU_pmaxi <- ceiling(HU_pmax) 

pnap <- length(HU_pmini)

#vizualizacio 
par(mfrow=c(1,1))


#kirajzolas2 : szimulalalt adatok
#hany napot rajzoljon ki : kir
kir <- 14
pnap <- c(1:kir)

plot(pnap, HU_pmaxi[1:kir], xlab = "Napok, 17. beteg megjelenése utáni 14 nap", ylab = "Igazolt fertőzött betegek száma, becslés (CI~0.92)", col=8, lty = 2, pch = 16, cex=0, ylim=c(1,max(HU_pmaxi[1:kir])))
points (HU_pmini[1:kir], col =8, pch = 16, cex=0)
polygon(c(pnap[1:kir], rev(pnap[1:kir])), c(HU_pmaxi[1:kir] ,rev(HU_pmini[1:kir])), col = rgb(1, 0, 0,0.5) )
points(pnap, HU_pi[1:kir], col=1, lty = 2, pch = 16, cex=1, type="b")

#SIR kod egy resze nem az en szellemi termekem, ezert nem oszthatom meg nyilvanosan, de ha vki ir emailt, kuldom
