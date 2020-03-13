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


# parameterek
min_fert <- 19 #honnantol veszi figyelembe az adatsort
min_hossz <- 4 #milyen hosszo adatsort tekint elemezhetonek min
shossz <- 21  # hany napra vetitse elore az eredmenyeket a modell


#BETOLTES
#adatok betoltese Johns Hopkins egyetem folyamatos frissitesu online adataibol
rawcovid <- read.csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_19-covid-Confirmed.csv", header=TRUE, as.is=T, sep=",", fileEncoding = "UTF-8") 

#OECD orszagok alapadatai
oecdraw <- read.csv("https://raw.githubusercontent.com/bricoler/naiv_COVID19_HU/R_code/OECD_2019.csv", header=TRUE, as.is=T, sep=",", fileEncoding = "UTF-8") 


#AGGREGALAS
#aggregalas orszagokra 
#E: itt eldobjuk a GEO-t, de a jovoben lehetne kezdeni valamit a tavolsagi adatokkal, ketdimenzioban (ido, ter)
l2 <- dim(rawcovid)[2] # napok hossza

rcovid <- aggregate(rawcovid[5:l2], by = list(rawcovid$Country.Region), FUN = sum) # orszagok szerinti aggregalas rcovid
rownames(rcovid) <- rcovid[,1] #sornevek javitasa rcovidban
colnames(rcovid)[1] <- "cname"

rownames(oecdraw) <- oecdraw[,1] #sornevek oecdrawban is


# OECD orszagok vizsgalata csak

rrcovid <- rcovid[(rcovid$cname %in% oecdraw$cname),] ##rrcovid: csak OECD orszagok
rrcovid <- rrcovid[,-1] # elso oszlop a cname torlese rrcovidban


#dif matrix eloallitasa (delta elozo nap)
d_covid <- rrcovid[1:dim(rrcovid)[1],2:dim(rrcovid)[2]] - rrcovid[1:dim(rrcovid)[1],(1:(dim(rrcovid)[2]-1))]

#az OECD lakossagszam leosztas
#covid : lakossagszamaranyos fertozottek aranya OECD orszagokban

covid <- rrcovid

for (i in 1:dim(rrcovid)[1]) 
{
covid[i,] <- rrcovid[i,] / oecdraw[rownames(rrcovid)[i],"pop"]
} 

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
cmin[i] <- min(as.numeric(covid[i,]), na.rm=T) * oecdraw[rownames(rrcovid)[i],"pop"] # minimum fertozottszam emberben
cmax[i] <- max(as.numeric(covid[i,]), na.rm=T) * oecdraw[rownames(rrcovid)[i],"pop"] # maximum fertozottszam
cmean[i] <- mean(as.numeric(covid[i,]), na.rm=T) * oecdraw[rownames(rrcovid)[i],"pop"] # atlagos fertozottszam
cmedian[i] <- median(as.numeric(covid[i,]), na.rm=T) * oecdraw[rownames(rrcovid)[i],"pop"] #median fertozottszam
csd[i] <- sd(as.numeric(covid[i,]), na.rm=T) * oecdraw[rownames(rrcovid)[i],"pop"] #sd fertozottszam emberben

tempv <- NULL 
tempv2 <- NULL

tempv <- na.omit(as.numeric(covid[i,])) # atmeneti vektor a nullak kidobalasa

tempv2 <- tempv[which( tempv * oecdraw[rownames(rrcovid)[i],"pop"] > min_fert) ] # 17-nel kisebbek kidobalasa - (tempv * vissza kell szorozni lakossagszamra)
chossz[i] <- length(tempv2) # hany nap van tobb mint min_fert fertozott

#ha minhossznak kevesebb a megfigyeles akkor nem szamolunk
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
c_tend[i]<- "ned2c" # not enough data to calculate
}

}

cmeta <- data.frame(cname,chossz,cmin,cmax,cmean,cmedian,csd, c_expC, c_expbeta, c_expR2, c_logC, c_logbeta, c_logR2, c_tend) 

#cmeta, orszagokrol fertozesmetaadatok

View(cmeta)

#Ahol meg keves adat volt fuggvenyszamitashoz, annak kidobasa

cmeta_exp <- cmeta[as.numeric(which(cmeta$c_tend!="ned2c")),]


#oecd orszagadatok hozzafuzzese a szamitott adatokhoz

cmeta_exp_oecd <- merge(cmeta_exp, oecdraw) # orszagok metadatai exp es log fuggveny es OECD adatok


#tesztek fuggvenyparameterek eloszlasara

shapiro.test(cmeta_exp_oecd$c_expC)[2]
shapiro.test(cmeta_exp_oecd$c_expbeta)[2]
shapiro.test(cmeta_exp_oecd$c_expR2)[2]

shapiro.test(cmeta_exp_oecd$c_logC)[2]
shapiro.test(cmeta_exp_oecd$c_logbeta)[2]
shapiro.test(cmeta_exp_oecd$c_logR2)[2]


#exp_beta
HU_beta <- mean(cmeta_exp_oecd$c_expbeta)
error_beta <- qnorm(0.95)*sd(cmeta_exp_oecd$c_expbeta)/sqrt(dim(cmeta_exp_oecd)[2]-1) # 90% CI

HU_betamax <- HU_beta + error_beta
HU_betamin <- HU_beta - error_beta

#exp_C
HU_C <- mean(cmeta_exp_oecd$c_expC)
error_c <- qnorm(0.95)*sd(cmeta_exp_oecd$c_expC)/sqrt(dim(cmeta_exp_oecd)[2]-1) # 90% CI

HU_Cmax <- HU_C + error_c
HU_Cmin <- HU_C - error_c

#expR2
HU_R2 <- mean(cmeta_exp_oecd$c_expR2)
error_R2 <- qnorm(0.975)*sd(cmeta_exp_oecd$c_expR2)/sqrt(dim(cmeta_exp_oecd)[2]-1)

HU_R2max <- HU_R2 + error_R2
HU_R2min <- HU_R2 - error_R2


# hazai adatok szimulalasa shossz n
pnap <- c(1:shossz)

HU_p <- NULL
HU_pmin <- NULL
HU_pmax <- NULL

# e kitevojeben a josolt beta kiszamitasa

for (i in 1:shossz) {
HU_p[i] <- exp(1)**(HU_C+(HU_beta*i))
HU_pmin[i] <- exp(1)**(HU_Cmin+(HU_betamin*i))
HU_pmax[i] <- exp(1)**(HU_Cmax+(HU_betamax*i))
}

#R2 bizonytalansag fapados hozzaszamitasa
HU_pmin <- (HU_pmin) * HU_R2
HU_pmax <- (HU_pmax) * (1+ (1-HU_R2))



#Kerekites egesz emberre, hazai populacioval felszorzas
HU_pi <- ceiling((HU_p) * oecdraw["Hungary","pop"])
HU_pmini <- ceiling((HU_pmin) * oecdraw["Hungary","pop"])
HU_pmaxi <- ceiling((HU_pmax) * oecdraw["Hungary","pop"])

# eredmenyek dataframe

results <- data.frame(pnap, HU_p, HU_pmin, HU_pmax, HU_pi, HU_pmini, HU_pmaxi)

#vizualizacio 
par(mfrow=c(1,1))


#kirajzolas2 : szimulalalt adatok
#hany napot rajzoljon ki : kir


plot(pnap, HU_pmax[1:shossz], xlab = "Napok, X. beteg megjelenése után", ylab = "Igazolt fertőzött betegek aránya, becslés (CI 0.82 - 0.86)", col=8, lty = 2, pch = 16, cex=0)
points (HU_pmin[1:shossz], col =8, pch = 16, cex=0)
polygon(c(pnap[1:shossz], rev(pnap[1:shossz])), c(HU_pmax[1:shossz] ,rev(HU_pmin[1:shossz])), col = rgb(1, 0, 0,0.5) )
points(pnap, HU_p[1:shossz], col=1, lty = 2, pch = 16, cex=1, type="b")

#eredmenyek kiiratasa
View(results)
#technikai reszletek kiirasa
#write.csv(cmeta_exp_oecd, file = "D:/R/Data/covid19_hu_tec_res1.csv", fileEncoding = "UTF-8")
#joslas kiirasa
#write.csv(results, file = "D:/R/Data/covid19_hu_res1.csv",fileEncoding = "UTF-8")


#SIR kod egy resze nem az en szellemi termekem, ezert nem oszthatom meg nyilvanosan, de ha vki ir emailt, kuldom
