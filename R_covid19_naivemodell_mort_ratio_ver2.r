#31-03-2020
#Elozetes halalozasi aranyok szamitasa idosoros adatokon
#Explorativ, vizualizacio

#memoria uritese 
rm(list = ls()) 

#BETOLTES
#adatok betoltese Johns Hopkins egyetem folyamatos frissitesu online adataibol
rawcovid <- read.csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_confirmed_global.csv", header=TRUE, as.is=T, sep=",", fileEncoding = "UTF-8")  #rawcovid az igazolt betegek

rawcovidd <- read.csv("https://raw.githubusercontent.com/CSSEGISandData/COVID-19/master/csse_covid_19_data/csse_covid_19_time_series/time_series_covid19_deaths_global.csv", header=TRUE, as.is=T, sep=",", fileEncoding = "UTF-8") #rawcovidd az igazolt halalesetek

#OECD orszagok alapadatai -- Ebben az elemzesben nem hasznalom, de kesobb bovitheto, ezert bennehagyom
# oecdraw <- read.csv("https://raw.githubusercontent.com/bricoler/naiv_COVID19_HU/R_code/OECD_2019.csv", header=TRUE, as.is=T, sep=",", fileEncoding = "UTF-8") 


#AGGREGALAS
#aggregalas orszagokra igazolt betegek
#E: itt eldobjuk a GEO-t, de a jovoben lehetne kezdeni valamit a tavolsagi adatokkal, ketdimenzioban (ido, ter)
l2 <- dim(rawcovid)[2] # napok hossza
rcovid <- aggregate(rawcovid[5:l2], by = list(rawcovid$Country.Region), FUN = sum) # orszagok szerinti aggregalas rcovid
rownames(rcovid) <- rcovid[,1] #sornevek javitasa rcovidban
colnames(rcovid)[1] <- "cname"

# Dél-Korea nevének javítása
rownames(rcovid)[which(rownames(rcovid) %in% "Korea, South")] <- "Korea South"
rcovid$cname[rcovid$cname=="Korea, South"] <- "Korea South"
rrcovid <- rcovid[,-1] # elso oszlop a cname torlese rrcovidban

#aggregalas a halalesetekre
#E: itt eldobjuk a GEO-t, de a jovoben lehetne kezdeni valamit a tavolsagi adatokkal, ketdimenzioban (ido, ter)
l2 <- dim(rawcovidd)[2] # napok hossza
rcovidd <- aggregate(rawcovidd[5:l2], by = list(rawcovidd$Country.Region), FUN = sum) # orszagok szerinti aggregalas rcovidd
rownames(rcovidd) <- rcovidd[,1] #sornevek javitasa rcovidban
colnames(rcovidd)[1] <- "cname"
# Dél-Korea nevének javítása
rownames(rcovidd)[which(rownames(rcovidd) %in% "Korea, South")] <- "Korea South"
rcovidd$cname[rcovidd$cname=="Korea, South"] <- "Korea South"
# OECD orszagok vizsgalata csak
##rrcovidd <- rcovidd[(rcovidd$cname %in% oecdraw$cname),] ##rrcovidd: csak OECD orszagok
colnames(rcovidd)[1] <- "cname"
rrcovidd <- rcovidd[,-1] # elso oszlop a cname torlese rrcovidban
#az OECD lakossagszam leosztas
#covid : lakossagszamaranyos fertozottek aranya OECD orszagokban


rcovid[rcovid==0]<-NA #0-kbol NA 
rcovidd[rcovidd==0]<-NA #0-kbol NA 


#covid az igazolt halalesetek es az igazol beteszam aranya
covid <- rrcovidd / rrcovid

covid[covid==0]<-NA #0-kbol NA 

#Populacioval valo osztas, nem ertelmes
#for (i in 1:dim(rrcovid)[1]) 
#{
#covid[i,] <- rrcovid[i,] / oecdraw[rownames(rrcovid)[i],"pop"]
#} 


#HALALOZASI ARANY ORSZAGONKENT
cname <- NULL # orszag neve
rmin <- NULL # legalacsonyabb arany (igazolt halalozas es igazolt betegszam hanyadoasa)
rmax <- NULL # legmagasabb arany (igazolt halalozas es igazolt betegszam hanyadoasa)
rmean <-NULL # atlag
rmedian <-NULL #median
rsd <-NULL # szorasa

rni2 <- NULL # shapiroval normal p ertek

rlen <- NULL # hossza

rsu <- NULL # igazolt halalozas es igazolt betegszam hanyadoasanak 95% UP
rsd <- NULL # igazolt halalozas es igazolt betegszam hanyadoasanak 95% DOWN


for (i in 1:dim(covid)[1]) 
{
cname[i] <- rownames(covid)[i] # orszag neve
rmin[i] <- min(as.numeric(covid[i,]), na.rm=T) 
rmax[i] <- max(as.numeric(covid[i,]), na.rm=T) 
rmean[i] <- mean(as.numeric(covid[i,]), na.rm=T) 
rmedian[i] <- median(as.numeric(covid[i,]), na.rm=T) 
rsd[i] <- sd(as.numeric(covid[i,]), na.rm=T) 
rlen[i]  <- sum(as.numeric(as.numeric(covid[i,])>0), na.rm=T)

rni2[i] <- ifelse (rlen[i]>5,as.numeric(shapiro.test(as.numeric(covid[i,]))[2]),0)

rsu[i] <- rmean[i] + qnorm(0.95)*rsd[i]/sqrt(rlen[i]) # 90% CI
rsd[i] <- rmean[i] + qnorm(0.95)*rsd[i]/sqrt(rlen[i]) # 90% CI



}

rresults <- data.frame(cname,rlen,rmean,rmax,rmin,rmedian,rsd,rsu,rsd,rni2, stringsAsFactors = TRUE)
colnames(rresults)[1] <- "cname"
View(rresults) 

# write.csv(rresults, file = "D:/R/Data/covid19_IEF.csv",fileEncoding = "UTF-8")


#Kornyezo orszagok tukreben a ggplot vizualizacio
#Ukrajnat azert hagyom ki, mert extrem magas

raw_kornyezo_oecd <- covid[rownames(covid) %in% c("Hungary", "Romania", "Serbia", "Slovakia", "Croatia", "Slovenia","Austria"),]

raw_kornyezo_oecd2 <- t(raw_kornyezo_oecd)

#install.packages("ggplot2")
library("ggplot2")
#install.packages("reshape")
library(reshape)

g_raw <- melt(raw_kornyezo_oecd2)
colnames(g_raw) <- c("Datum" , "Orszag", "COVID_19_halálozasi_arany_atlaga_es_eloszlasa")

# Box plot
boxp <- ggplot(g_raw, aes(Orszag, COVID_19_halálozasi_arany_atlaga_es_eloszlasa)) + 
  geom_boxplot(aes(fill = Orszag)) +
  theme_minimal() +
  theme(legend.position = "none"
  )
boxp

#VEGE
