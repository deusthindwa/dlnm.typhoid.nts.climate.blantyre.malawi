#02/08/2018-07/11/2019
#by Deus Thindwa

#----------Install required packages.
dlnm.analysis.packages <- c("tidyverse", "lubridate","xts","ggthemes","PerformanceAnalytics","reshape2","rugarch","timetk","parallel","timeSeries","tseries","data.table","ggplot2","dlnm","broom","caret","gridExtra","splines","splines2","pspline","cowplot","mgcv","spi","chron","gridGraphics","grid","pscl","MASS", "AER", "Hmisc", "MuMIn", "VGAM", "forecast", "seasonal", "plotly", "ggmap", "rgeos", "tmap", "maptools", "maps", "ggfortify", "htmltools","webshot","knitr","flexdashboard", "imager", "httr", "gmodels", "curl", "here")

#----------load required packages.
lapply(dlnm.analysis.packages, library, character.only=TRUE)

#----------load shape file of malawi map.
dlnmtmp <- tempfile()
download.file("https://raw.githubusercontent.com/deusthindwa/dlnm.typhoid.nts.climate.blantyre.malawi/master/data/malawi_map.zip", destfile=dlnmtmp)
unzip(dlnmtmp, exdir = ".")
malawi.map <- rgdal::readOGR(".","malawi_map")


#----------subsetting to get blantyre map only.
blantyre1.map <- malawi.map@data$OBJECTID >289 & malawi.map@data$OBJECTID <297 #id from 290 to 296 
blantyre2.map <- malawi.map@data$OBJECTID >308 & malawi.map@data$OBJECTID <311 #id from 309 to 310
blantyre3.map <- malawi.map@data$OBJECTID >342  #id fom 243


#----------convert shape file map into dataframe for ggplot.
blantyre.map <- rbind(fortify(malawi.map[blantyre1.map,]), fortify(malawi.map[blantyre2.map,]), fortify(malawi.map[blantyre3.map,]))
blantyre.map$id <- as.integer(blantyre.map$id)

#----------merge blantyre map dataset with location attributes dataset.
blantyre.demog <- read.csv(curl("https://raw.githubusercontent.com/deusthindwa/dlnm.typhoid.nts.climate.blantyre.malawi/master/data/blantyre_demog.csv"))
map.features <- read.csv(curl("https://raw.githubusercontent.com/deusthindwa/dlnm.typhoid.nts.climate.blantyre.malawi/master/data/blantyre_features.csv"))
blantyre.demog$id <- as.integer(blantyre.demog$id)
map.all <- merge(x=blantyre.map, y=blantyre.demog, by="id", x.all=TRUE)
rm(list = ls()[grep("^blantyre", ls())])

#----------plot blantyre map with 1998-2008 population census.
ggplot() + 
  geom_polygon(data=map.all, aes(x=long, y=lat, group=group, fill=popc), colour="gray50") + 
  theme_classic() + 
  theme(axis.text.x = element_text(face="bold", size=10, color="black"), axis.text.y = element_text(face="bold", size=10, color="black")) + 
  labs(fill="(1998 - 2008) Population censuses") + xlab("Longitude") + ylab("Latitude") + 
  geom_point(data=map.features, aes(x =long, y =lat, shape=Geolocation, size=Geolocation), color="black") +
  scale_shape_manual(values=c(17, 16, 3)) +
  scale_size_manual(values=c(2,4,3)) + 
  theme(legend.key.height=unit(0.8,"line")) + 
  theme(legend.key.width=unit(0.8,"line"))

#----------load typhoid and NTS cases dataset.
case <-read.csv(curl("https://raw.githubusercontent.com/deusthindwa/dlnm.typhoid.nts.climate.blantyre.malawi/master/data/case.csv"))
case$case_date <- dmy(case$case_date)
case$case_count <- c(1)

#----------create separate datasets for typhi and NTS cases.
case.typhi <- subset(case, organism == "typhi")
case.iNTS <- subset(case, organism == "iNTS")

#----------assign 0 to case_count when a date has no typhi case.
case.typhi <-aggregate(case.typhi$case_count, by=list(case.typhi$case_date), FUN=sum, na.rm=TRUE)
names(case.typhi) <- c("date", "case_count")
case.typhi <- merge(case.typhi, data.table(date=seq.Date(min(case$case_date), max(case$case_date), by="day")), by="date", all=TRUE)
case.typhi[is.na(case.typhi)] <- 0

#----------assign 0 to case_count when dates have no iNTS case.
case.iNTS <-aggregate(case.iNTS$case_count, by=list(case.iNTS$case_date), FUN=sum, na.rm=TRUE)
names(case.iNTS) <- c("date", "case_count")
case.iNTS <- merge(case.iNTS, data.table(date=seq.Date(min(case$case_date), max(case$case_date), by="day")), by="date", all=TRUE)
case.iNTS[is.na(case.iNTS)] <- 0

#----------convert dataframes to xts objects for time series plot.
case.typhi = as.xts(case.typhi[,-1,drop = FALSE], order.by = as.Date(case.typhi[,1]))
case.iNTS = as.xts(case.iNTS[,-1,drop = FALSE], order.by = as.Date(case.iNTS[,1]))

#----------weekly or monthly sum for typhi and iNTS cases.
case.typhiW <- apply.weekly(case.typhi, FUN = sum)
case.iNTSW <- apply.weekly(case.iNTS, FUN = sum)
case.typhi <- apply.monthly(case.typhi, FUN = sum)
case.iNTS <- apply.monthly(case.iNTS, FUN = sum)

#----------load climate dataset.
climate <- read.csv(curl("https://raw.githubusercontent.com/deusthindwa/dlnm.typhoid.nts.climate.blantyre.malawi/master/data/climate.csv"))

#----------average daily values for temperature and rainfall from 2 stations.
climate$climate_date <- dmy(climate$date)
climate$rainfall <- (climate$chil_r + climate$chic_r)/2 
climate$temperature <- ((climate$chil_mint + climate$chil_maxt)/2 + (climate$chic_mint + climate$chic_maxt)/2)/2
climate <- subset(climate, select = c(climate_date, rainfall, temperature))

#----------create separate datasets for daily temperature and rainfall.
climate.rain <- subset(climate, select = c(climate_date, rainfall))
climate.temp <- subset(climate, select = c(climate_date, temperature))

#----------convert temperature and rainfall data frames to xts objects for use in time series plotting.
climate.rain = as.xts(climate.rain[,-1,drop = FALSE], order.by = as.Date(climate.rain[,1]))
climate.temp = as.xts(climate.temp[,-1,drop = FALSE], order.by = as.Date(climate.temp[,1]))

#----------monthly mean for rainfall and temperature
climate.rain <- apply.monthly(climate.rain, FUN = mean)
climate.temp <- apply.monthly(climate.temp, FUN = mean)

#----------create tibbles for case and climate seasonal-adjustment
case.typhi <-tk_tbl(case.typhi, preserve_index = TRUE, rename_index = "date") 
case.iNTS <-tk_tbl(case.iNTS, preserve_index = TRUE, rename_index = "date") 
climate.rain <-tk_tbl(climate.rain, preserve_index = TRUE, rename_index = "date") 
climate.temp <-tk_tbl(climate.temp, preserve_index = TRUE, rename_index = "date") 

#----------seasonally-adjusted cases and climate
case.iNTS.ts <- ts(na.omit(case.iNTS$case_count), frequency = 12)
trend_n <-tk_tbl(exp(seasadj(mstl(log(case.iNTS.ts)))), preserve_index = FALSE) #multiplicative series (log-transform); already has at least 1 NTS case.
case.iNTS$case_count_sea <- trend_n$value  #seasonally-adjusted cases: trend+remainder
case.iNTS$case_count_sea[case.iNTS$case_count_sea < 0] <- 0
setnames(case.iNTS, old="case_count", new="case_count_obs")

case.typhi.ts <- ts(na.omit(case.typhi$case_count), frequency = 12)
trend_n <-tk_tbl(exp(seasadj(mstl(log(case.typhi.ts+1)))), preserve_index = FALSE) #multiplicative series (log-transform); add 1 since log(0) is not defined.
case.typhi$case_count_sea <- trend_n$value #seasonally-adjusted cases: trend+remainder
case.typhi$case_count_sea[case.typhi$case_count_sea < 0] <- 0
setnames(case.typhi, old="case_count", new="case_count_obs")

climate.rain.ts <- ts(na.omit(climate.rain$rainfall), frequency = 12)
trend_n <-tk_tbl(seasadj(mstl(climate.rain.ts)), preserve_index = FALSE) #additive series
climate.rain$rainfall_sea <- trend_n$value #seasonally-adjusted rainfall: trend+remainder
climate.rain$rainfall_sea[climate.rain$rainfall_sea < 0] <- 0
setnames(climate.rain, old="rainfall", new="rainfall_obs")

climate.temp.ts <- ts(na.omit(climate.temp$temperature), frequency = 12)
trend_n <-tk_tbl(seasadj(mstl(climate.temp.ts)), preserve_index = FALSE) #additive series
climate.temp$temperature_sea <- trend_n$value #seasonally-adjusted temperature: trend+remainder
climate.temp$temperature_sea[climate.temp$temperature_sea < 0] <- 0
setnames(climate.temp, old="temperature", new="temperature_obs")

#----------plot decomposed all series
x<-mstl(case.iNTS.ts,s.window="period") %>% ggfortify:::autoplot.ts(main="A",xlab="Years (2000-2015)",ylab="Number of iNTS cases",size=1,colour="orange2",is.date=FALSE) + theme_bw()
y<-mstl(case.typhi.ts,s.window="periodic") %>% ggfortify:::autoplot.ts(main="B",xlab="Years (2000-2015)",ylab="Number of typhoid cases",size=1,colour="red2",is.date=FALSE) + theme_bw()
z<-mstl(climate.rain.ts,s.window="periodic") %>% ggfortify:::autoplot.ts(main="C",xlab="Years (2000-2015)",ylab="Rainfall (mm)",size=1,colour="blue2",is.date=FALSE) + theme_bw()
v<-mstl(climate.temp.ts,s.window="periodic") %>% ggfortify:::autoplot.ts(main="D",xlab="Years (2000-2015)",ylab="Temperature (°C)",size=1,colour="green2",is.date=FALSE) + theme_bw()
grid.arrange(grobs=list(x,y,z,v), ncol=4, nrow=1)

#----------plot observed vs seasonal-adjusted series
S1<-ggplot(as.data.frame(case.iNTS)) + 
  geom_line(aes(date, case_count_obs, color="Observed data"), size=0.8) + 
  geom_line(aes(date, case_count_sea, color="Seasonal-adjusted"), size=0.8) + 
  scale_color_manual(values = c("Observed data"="black","Seasonal-adjusted"="orange2")) +
  labs(title="A", x ="", y = "iNTS cases") + 
  theme(plot.title = element_text(hjust = 0)) + 
  theme(axis.title.x = element_text(size = 10)) + 
  theme(axis.title.y = element_text(size = 10)) +
  theme(axis.text.x = element_text(face="bold", size=10), axis.text.y = element_text(face="bold", size=10)) + 
  theme(legend.justification=c(0.5,0), legend.position = c(0.7,0.7), legend.text = element_text(size = 10), legend.title = element_text(face="bold", size=0)) + 
  theme(legend.key.height=unit(1,"line")) + 
  theme(legend.key.width=unit(1,"line"))

S2<-ggplot(as.data.frame(case.typhi)) + 
  geom_line(aes(date, case_count_obs, color="Observed data"), size=0.8) + 
  geom_line(aes(date, case_count_sea, color="Seasonal-adjusted"), size=0.8) + 
  scale_color_manual(values = c("Observed data"="black","Seasonal-adjusted"="red2")) +
  labs(title="B", x ="", y = "Typhoid cases") + 
  theme(plot.title = element_text(hjust = 0)) + 
  theme(axis.title.x = element_text(size = 10)) + 
  theme(axis.title.y = element_text(size = 10)) +
  theme(axis.text.x = element_text(face="bold", size=10), axis.text.y = element_text(face="bold", size=10)) + 
  theme(legend.justification=c(0.5,0), legend.position = c(0.3,0.7), legend.text = element_text(size = 10), legend.title = element_text(face="bold", size=0)) + 
  theme(legend.key.height=unit(1,"line")) + 
  theme(legend.key.width=unit(1,"line"))

S3<-ggplot(as.data.frame(climate.rain)) + 
  geom_line(aes(date, rainfall_obs, color="Observed data"), size=0.8) + 
  geom_line(aes(date, rainfall_sea, color="Seasonal-adjusted"), size=0.8) + 
  scale_color_manual(values = c("Observed data"="black","Seasonal-adjusted"="blue2")) +
  labs(title="C", x ="Year", y = "Rainfall (mm)") + 
  theme(plot.title = element_text(hjust = 0)) + 
  theme(axis.title.x = element_text(size = 10)) + 
  theme(axis.title.y = element_text(size = 10)) +
  theme(axis.text.x = element_text(face="bold", size=10), axis.text.y = element_text(face="bold", size=10)) + 
  theme(legend.justification=c(0.5,0), legend.position = c(0.3,0.7), legend.text = element_text(size = 10), legend.title = element_text(face="bold", size=0)) + 
  theme(legend.key.height=unit(1,"line")) + 
  theme(legend.key.width=unit(1,"line"))

S4<-ggplot(as.data.frame(climate.temp)) + 
  geom_line(aes(date, temperature_obs, color="Observed data"), size=0.8) + 
  geom_line(aes(date, temperature_sea, color="Seasonal-adjusted"), size=0.8) + 
  scale_color_manual(values = c("Observed data"="black","Seasonal-adjusted"="green2")) + 
  labs(title="D", x ="Year", y = "Temperature (°C)") + 
  theme(plot.title = element_text(hjust = 0)) + 
  theme(axis.title.x = element_text(size = 10)) + 
  theme(axis.title.y = element_text(size = 10)) + 
  theme(axis.text.x = element_text(face="bold", size=10), axis.text.y = element_text(face="bold", size=10)) + 
  theme(legend.justification=c(0.5,0), legend.position = c(0.3,0.7), legend.text = element_text(size = 10), legend.title = element_text(face="bold", size=0)) + 
  theme(legend.key.height=unit(1,"line")) + 
  theme(legend.key.width=unit(1,"line"))

grid.arrange(grobs=list(S1, S2, S3, S4), ncol=2, nrow=2)
rm(list = ls()[grep("^trend_n", ls())])

#----------linear interpolation and extrapolation
census.year <-c(1998, 2008)
census.popn <-c(809397, 1022680)
census.count.1998.2008 <- approx(census.year, census.popn, n=11) 
census.count.2009.2015 <- approxExtrap(census.year, census.popn, xout=c(2009, 2010, 2011, 2012, 2013, 2014, 2015))

#----------calculated incidence of monthly iNTS
case.iNTS$census[year(case.iNTS$date)==2000]<- 852054; case.iNTS$census[year(case.iNTS$date)==2001]<-873382
case.iNTS$census[year(case.iNTS$date)==2002]<- 894710; case.iNTS$census[year(case.iNTS$date)==2003]<-916039
case.iNTS$census[year(case.iNTS$date)==2004]<- 937367; case.iNTS$census[year(case.iNTS$date)==2005]<-958695
case.iNTS$census[year(case.iNTS$date)==2006]<- 980023; case.iNTS$census[year(case.iNTS$date)==2007]<-1001352
case.iNTS$census[year(case.iNTS$date)==2008]<-1022680; case.iNTS$census[year(case.iNTS$date)==2009]<-1044008
case.iNTS$census[year(case.iNTS$date)==2010]<-1065337; case.iNTS$census[year(case.iNTS$date)==2011]<-1086665
case.iNTS$census[year(case.iNTS$date)==2012]<-1107993; case.iNTS$census[year(case.iNTS$date)==2013]<-1129322
case.iNTS$census[year(case.iNTS$date)==2014]<-1150650; case.iNTS$census[year(case.iNTS$date)==2015]<-1171978
case.iNTS$incid_sea <-case.iNTS$case_count_sea*100000/case.iNTS$census 
case.iNTS$incid_obs <-case.iNTS$case_count_obs*100000/case.iNTS$census 

#----------calculated incidence of monthly typhoid
case.typhi$census[year(case.typhi$date)==2000]<- 852054; case.typhi$census[year(case.typhi$date)==2001]<-873382
case.typhi$census[year(case.typhi$date)==2002]<- 894710; case.typhi$census[year(case.typhi$date)==2003]<-916039
case.typhi$census[year(case.typhi$date)==2004]<- 937367; case.typhi$census[year(case.typhi$date)==2005]<-958695
case.typhi$census[year(case.typhi$date)==2006]<- 980023; case.typhi$census[year(case.typhi$date)==2007]<-1001352
case.typhi$census[year(case.typhi$date)==2008]<-1022680; case.typhi$census[year(case.typhi$date)==2009]<-1044008
case.typhi$census[year(case.typhi$date)==2010]<-1065337; case.typhi$census[year(case.typhi$date)==2011]<-1086665
case.typhi$census[year(case.typhi$date)==2012]<-1107993; case.typhi$census[year(case.typhi$date)==2013]<-1129322
case.typhi$census[year(case.typhi$date)==2014]<-1150650; case.typhi$census[year(case.typhi$date)==2015]<-1171978
case.typhi$incid_sea <-case.typhi$case_count_sea*100000/case.typhi$census 
case.typhi$incid_obs <-case.typhi$case_count_obs*100000/case.typhi$census 

#----------boxplots of weekly and month seasonal dynamics of iNTS and typhoid
j <- seq(from=1, to=53, by=4)
i <- seq(from=1, to=12, by=1)

pox2 <- ggplot(subset(case.iNTS, year(date)<2011), aes(x=month(date), y=incid_obs))  + 
  geom_boxplot(aes( group=month(date)), color="orange2", fill="orange2", alpha=0.2) + 
  labs(title="B",x="Month (Jan-Dec)", y = "Monthly iNTS incidence") + 
  scale_x_discrete(limits = i) + 
  ylim(c(0,25)) +
  theme(axis.title.y = element_text(size = 11)) + 
  theme(axis.title.x = element_text(size = 11)) + 
  theme(axis.text.x = element_text(face="bold", size=11), axis.text.y = element_text(face="bold", size=11)) 

pox4 <- ggplot(subset(case.typhi, year(date)>2010), aes(x=month(date), y=incid_obs))  + 
  geom_boxplot(aes( group=month(date)), color="red2", fill="red2", alpha=0.2) + 
  labs(title="D", x="Month (Jan-Dec)", y = "Monthly typhoid incidence") + 
  scale_x_discrete(limits =i) + 
  ylim(c(0,25)) +
  theme(axis.title.y = element_text(size = 11)) + 
  theme(axis.title.x = element_text(size = 11)) + 
  theme(axis.text.x = element_text(face="bold", size=11), axis.text.y = element_text(face="bold", size=11)) 

grid.arrange(grobs=list(pox2, pox4), ncol=1, nrow=2)

#----------distributions of typhi and iNTS cases by sex and age.
case$sex[case$sex == ""] <- NA
case$age[case$age == ""] <- NA
case$date <- ymd(case$case_date)
case$year <- year(case$date)
dat<-case
dat$ageRounded<-floor(dat$age)
datNoUnknowns<-dat[dat$sex!="Unknown",]
datNoUnknowns$sex<-factor(datNoUnknowns$sex)

dev.off()
agesex.p1<-datNoUnknowns %>%
  filter(organism=="iNTS") %>%
  count(sex,ageRounded) %>%
  complete(sex,ageRounded,fill=list(n=0)) %>%
  ggplot(mapping=aes(fill=sex,y=n,x=ageRounded)) +
  geom_bar(position="dodge",stat="identity") +
  theme_bw() + 
  scale_fill_manual(values=c("steelblue","orange")) +
  scale_x_continuous(breaks = seq(0, 91, 5)) +
  labs(title="A", x ="Age (years)", y = "Number of iNTS cases") + 
  theme(axis.text.x = element_text(face="bold", size=10, color="black"), axis.text.y = element_text(face="bold", size=10, color="black")) + 
  theme(legend.justification=c(0.5,0), legend.position = c(0.5, 0.6), legend.text = element_text(size = 10), legend.title = element_text(size = 10)) + 
  labs(fill="Missing age and sex: 2,596 (32.3%)") + 
  theme(legend.key.height=unit(0.8,"line")) + 
  theme(legend.key.width=unit(0.8,"line"))

agesex.p2<-datNoUnknowns %>%
  filter(organism=="typhi") %>%
  count(sex,ageRounded) %>%
  complete(sex,ageRounded,fill=list(n=0)) %>%
  ggplot(mapping=aes(fill=sex,y=n,x=ageRounded)) +
  geom_bar(position="dodge",stat="identity") +
  theme_bw() + 
  scale_fill_manual(values=c("steelblue","orange")) +
  scale_x_continuous(breaks = seq(0, 91, 5)) +
  labs(title="B", x ="Age (years)", y = "Number of typhoid cases") +
  theme(axis.text.x = element_text(face="bold", size=10, color="black"), axis.text.y = element_text(face="bold", size=10, color="black")) + 
  theme(legend.justification=c(0.5,0), legend.position = c(0.5, 0.6), legend.text = element_text(size = 10), legend.title = element_text(size = 10)) + 
  labs(fill="Missing age and sex: 61 (2.4%)") + 
  theme(legend.key.height=unit(0.8,"line")) + 
  theme(legend.key.width=unit(0.8,"line"))

grid.arrange(agesex.p1,agesex.p2,nrow=2)

#----------contour plots of (de)seasonal dynamics of iNTS and typhoid
case.iNTS.spi <- subset(case.iNTS, year(case.iNTS$date)<2011, select=c(date,incid_sea)) 
case.iNTS.spi$month <- month(case.iNTS.spi$date)
case.iNTS.spi$year <- year(case.iNTS.spi$date)
case.iNTS.spi$date <- NULL
case.iNTS.spi <- spread(case.iNTS.spi, year, incid_sea)
case.iNTS.spi <- as.matrix(case.iNTS.spi)[,-1]
plot_ly(x = c(2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010), y = c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"), z = ~case.iNTS.spi, type = "contour", colorscale = 'heatmap', contours = list(showlabels = TRUE)) %>% 
colorbar(title = "Seasonal-adjusted \n iNTS incidence per \n 100,000 population") %>%
layout(title="<b>A</b>", xaxis=list(title ="Year",color="black"), yaxis=list(title="Month",color="black"), font=list(size = 13))

climate.rain.spin <- subset(climate.rain, year(climate.rain$date)<2011, select=c(date,rainfall_sea)) 
climate.rain.spin$month <- month(climate.rain.spin$date)
climate.rain.spin$year <- year(climate.rain.spin$date)
climate.rain.spin$date <- NULL
climate.rain.spin <- spread(climate.rain.spin, year, rainfall_sea)
climate.rain.spin <- as.matrix(climate.rain.spin)[,-1]
plot_ly(x = c(2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010), y = c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"), z = ~climate.rain.spin, type = "contour", colorscale = 'Earth', contours = list(showlabels = TRUE)) %>% 
colorbar(title = "Seasonal-adjusted \n Rainfall (mm)") %>%
layout(title="<b>B</b>", xaxis=list(title ="Year",color="black"), yaxis=list(title="Month",color="black"), font=list(size = 13))

climate.temp.spin <- subset(climate.temp, year(climate.temp$date)<2011, select=c(date,temperature_sea)) 
climate.temp.spin$month <- month(climate.temp.spin$date)
climate.temp.spin$year <- year(climate.temp.spin$date)
climate.temp.spin$date <- NULL
climate.temp.spin <- spread(climate.temp.spin, year, temperature_sea)
climate.temp.spin <- as.matrix(climate.temp.spin)[,-1]
plot_ly(x = c(2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010), y = c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"), z = ~climate.temp.spin, type = "contour", colorscale = 'Viridis', contours = list(showlabels = TRUE)) %>% 
colorbar(title = "Seasonal-adjusted \n Temperature (°C)") %>%
layout(title="<b>C</b>", xaxis=list(title ="Year",color="black"), yaxis=list(title="Month",color="black"), font=list(size = 13))

case.typhi.spi <- subset(case.typhi, year(case.typhi$date)>2010, select=c(date,incid_sea)) 
case.typhi.spi$month <- month(case.typhi.spi$date)
case.typhi.spi$year <- year(case.typhi.spi$date)
case.typhi.spi$date <- NULL
case.typhi.spi <- spread(case.typhi.spi, year, incid_sea)
case.typhi.spi <- as.matrix(case.typhi.spi)[,-1]
plot_ly(x = c(2011,2012,2013,2014,2015), y = c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"), z = ~case.typhi.spi, type = "contour", colorscale = 'heatmap', contours = list(showlabels = TRUE)) %>% 
colorbar(title = "Seasonal-adjusted \n typhoid incidence per \n 100,000 population") %>%
layout(title="<b>D</b>", xaxis=list(title ="Year",color="black"), yaxis=list(title="Month",color="black"), font=list(size = 13))

climate.rain.spit <- subset(climate.rain, year(climate.rain$date)>2010, select=c(date,rainfall_sea)) 
climate.rain.spit$month <- month(climate.rain.spit$date)
climate.rain.spit$year <- year(climate.rain.spit$date)
climate.rain.spit$date <- NULL
climate.rain.spit <- spread(climate.rain.spit, year, rainfall_sea)
climate.rain.spit <- as.matrix(climate.rain.spit)[,-1]
plot_ly(x = c(2011,2012,2013,2014,2015), y = c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"), z = ~climate.rain.spit, type = "contour", colorscale = 'Earth', contours = list(showlabels = TRUE)) %>% 
colorbar(title = "Seasonal-adjusted \n Rainfall (mm)") %>%
layout(title="<b>E</b>", xaxis=list(title ="Year",color="black"), yaxis=list(title="Month",color="black"), font=list(size = 13))

climate.temp.spit <- subset(climate.temp, year(climate.temp$date)>2010, select=c(date,temperature_sea)) 
climate.temp.spit$month <- month(climate.temp.spit$date)
climate.temp.spit$year <- year(climate.temp.spit$date)
climate.temp.spit$date <- NULL
climate.temp.spit <- spread(climate.temp.spit, year, temperature_sea)
climate.temp.spit <- as.matrix(climate.temp.spit)[,-1]
plot_ly(x = c(2011,2012,2013,2014,2015), y = c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec"), z = ~climate.temp.spit, type = "contour", colorscale = 'Viridis', contours = list(showlabels = TRUE)) %>% 
colorbar(title = "Seasonal-adjusted \n Temperature (°C)") %>%
layout(title="<b>F</b>", xaxis=list(title ="Year",color="black"), yaxis=list(title="Month",color="black"), font=list(size = 13))

#----------prepare final monthly NTS dataset for use in DLNM
mo.dlnmN <- bind_cols(case.iNTS, climate.rain, climate.temp, id=NULL)
mo.dlnmN$date1 <- mo.dlnmN$date2 <- NULL
mo.dlnmN <- subset(mo.dlnmN, year(date) < 2011)
mo.dlnmN$time <- seq.int(from = 1, to=132, by=1)
mo.dlnmN$year <- year(mo.dlnmN$date)
mo.dlnmN$month <- month(mo.dlnmN$date)
mo.dlnmN$incid_seaX<-round(mo.dlnmN$incid_sea, digits = 0)
mo.dlnmN$incid_obsX<-round(mo.dlnmN$incid_obs, digits = 0)

#----------prepare final monthly typhoid dataset for use in DLNM
mo.dlnmT <- bind_cols(case.typhi, climate.rain, climate.temp, id=NULL)
mo.dlnmT$date1 <- mo.dlnmT$date2 <- NULL
mo.dlnmT <- subset(mo.dlnmT, year(date) > 2010)
mo.dlnmT$time <- seq.int(from = 1, to=60, by=1)
mo.dlnmT$year <- year(mo.dlnmT$date)
mo.dlnmT$month <- month(mo.dlnmT$date)
mo.dlnmT$incid_seaX<-round(mo.dlnmT$incid_sea, digits = 0)
mo.dlnmT$incid_obsX<-round(mo.dlnmT$incid_obs, digits = 0)

#----------change AIC to QAIC for model comparisons
nts.quasipoisson <- function(...) { 
  res <- quasipoisson(...)
  res$aic <- poisson(...)$aic 
  res
}
typ.quasipoisson <- function(...) { 
  res <- quasipoisson(...) 
  res$aic <- poisson(...)$aic 
  res
}

#----------test all possible dfs combinations for rainfall, temperature and lag
QAICtable <- data.frame(model.no=rep(NA,27), lag.df=rep(NA,27), fx1.df=rep(NA,27), fx2.df=rep(NA,27), 
                      QAIC.ntsR=rep(NA,27), QAIC.ntsT=rep(NA,27), QAIC.nts=rep(NA,27), QAIC.typR=rep(NA,27), 
                      QAIC.typT=rep(NA,27), QAIC.typ=rep(NA,27))
l=1
for(k in 3:5){
  for(j in 3:5){
    for(i in 3:5){
  nts.lagknots <- logknots(8, fun="ns", df=i)
  nts.varknots.r=equalknots(mo.dlnmN$rainfall_obs, fun="ns", df=j)
  nts.varknots.t=equalknots(mo.dlnmN$temperature_obs, fun="ns", df=k)
  nts.mo.cb.rain <- crossbasis(mo.dlnmN$rainfall_obs, lag=8, argvar=list(fun="ns", knots=nts.varknots.r), arglag=list(knots=nts.lagknots))
  nts.mo.cb.temp <- crossbasis(mo.dlnmN$temperature_obs, lag=8, argvar=list(fun="ns", knots=nts.varknots.t), arglag=list(knots=nts.lagknots))
  nts.modelR <- glm(mo.dlnmN$incid_obsX ~ nts.mo.cb.rain + month + year, family=nts.quasipoisson(), na.action=na.delete, mo.dlnmN)
  nts.modelT <- glm(mo.dlnmN$incid_obsX ~ nts.mo.cb.temp + month + year, family=nts.quasipoisson(), na.action=na.delete, mo.dlnmN)
  nts.model <-  glm(mo.dlnmN$incid_obsX ~ nts.mo.cb.rain + month + nts.mo.cb.temp + year, family = nts.quasipoisson(), na.action=na.delete, mo.dlnmN)
  
  typ.lagknots <- logknots(8, fun="ns", df=i)
  typ.varknots.r=equalknots(mo.dlnmT$rainfall_obs, fun="ns", df=j)
  typ.varknots.t=equalknots(mo.dlnmT$temperature_obs, fun="ns", df=k)
  typ.mo.cb.rain <- crossbasis(mo.dlnmT$rainfall_obs, lag=8, argvar=list(fun="ns", knots=typ.varknots.r), arglag=list(knots=typ.lagknots))
  typ.mo.cb.temp <- crossbasis(mo.dlnmT$temperature_obs, lag =8, argvar = list(fun="ns", knots=typ.varknots.t), arglag=list(knots=typ.lagknots))
  typ.modelR <- glm(mo.dlnmT$incid_obsX ~ typ.mo.cb.rain + month + year, family = typ.quasipoisson(), na.action=na.delete, mo.dlnmT)
  typ.modelT <- glm(mo.dlnmT$incid_obsX ~ typ.mo.cb.temp + month + year, family = typ.quasipoisson(), na.action=na.delete, mo.dlnmT)
  typ.model <- glm(mo.dlnmT$incid_obsX ~ typ.mo.cb.rain + month + typ.mo.cb.temp + year, family = typ.quasipoisson(), na.action=na.delete, mo.dlnmT)
  
  QAICtable[l,] <- c(l,i,j,k, QAIC(nts.modelR, chat=summary(nts.modelR)$dispersion), QAIC(nts.modelT, chat=summary(nts.modelT)$dispersion), 
                     QAIC(nts.model, chat=summary(nts.model)$dispersion), QAIC(typ.modelR, chat=summary(typ.modelR)$dispersion), QAIC(typ.modelT, chat=summary(typ.modelT)$dispersion), 
                     QAIC(typ.model, chat=1))
  l=l+1  
  }
  }
}

#----------construct cross-basis for iNTS using optimal dfs to predict rainfall effect
sort(mo.dlnmN$rainfall_obs, decreasing=FALSE)
varknots=equalknots(mo.dlnmN$rainfall_obs, fun="ns", df=3)
lagknots <- logknots(8, fun="ns", df=3)
mo.cb.rain.iNTS <- crossbasis(mo.dlnmN$rainfall_obs, lag =8, argvar = list(knots=varknots), arglag = list(knots=lagknots))
summary(mo.cb.rain.iNTS)

#----------construct cross-basis for iNTS using optimal dfs to predict temperature effect
sort(mo.dlnmN$temperature_obs, decreasing=FALSE)
varknots=equalknots(mo.dlnmN$temperature_obs, fun="ns", df=3)
lagknots <- logknots(8, fun="ns", df=3)
mo.cb.temp.iNTS <- crossbasis(mo.dlnmN$temperature_obs, lag=8, argvar=list(knots=varknots), arglag=list(knots=lagknots))
summary(mo.cb.temp.iNTS)

#----------model fitting for iNTS
mo.model.iNTS <- glm(incid_obsX ~  mo.cb.rain.iNTS + mo.cb.temp.iNTS + month + year, family = quasipoisson(), na.action=na.exclude, mo.dlnmN)

#----------model validation check for iNTS
dev.off()
pacf(residuals(mo.model.iNTS,type="deviance"),na.action=na.omit,main="Partial autocorrelation (original model)",xlim=c(0,8))
mo.model.iNTS <- update(mo.model.iNTS,.~.+Lag(residuals(mo.model.iNTS,type="deviance"),1)) #add residuals at lag 1 to significantly reduce partial autocorrelation
pacf(residuals(mo.model.iNTS,type="deviance"),na.action=na.omit,main="Autocorrelation (adjusted model)",xlim=c(0,8))

#----------validated model predictions for iNTS
mo.pred.rain.iNTS <- crosspred(mo.cb.rain.iNTS, mo.model.iNTS, cen = 0, by=0.2)
mo.pred.temp.iNTS <- crosspred(mo.cb.temp.iNTS, mo.model.iNTS, cen = 23, by=0.2)

#----------plotting countour and curves for rainfall on iNTS
dev.off()
par(mar=c(5,5,2,2)+0.1)
plot(mo.pred.rain.iNTS, "contour", key.title=title("iNTS"), plot.title=title("", xlab ="Mean monthly rainfall (mm)", ylab = "Monthly lag", cex.lab=1.3, cex.axis=1.5,main="A"))
plot(mo.pred.rain.iNTS, "slices", xlab="Monthly lag (given mean 9 mm/month)", var=c(9), col="orange2", ci.arg=list(col=topo.colors(70, alpha = 1)), ci.level=0.95, ci='b', lwd=4.5, ylab="RR of iNTS", cex.lab=1.3, cex.axis=1.5,main="B")
plot(mo.pred.rain.iNTS, "slices", xlab="Monthly lag (given mean 13 mm/month)", var=c(13), col="orange2", ci.arg=list(col=topo.colors(70, alpha = 1)), ci.level=0.95, ci='b', lwd=4.5, ylab="RR of iNTS", cex.lab=1.3, cex.axis=1.5,main="C")

#----------plotting countour and curves for temperature on iNTS
dev.off()
par(mar=c(5,5,2,2)+1)
plot(mo.pred.temp.iNTS, "contour", key.title=title("iNTS"), plot.title=title("", xlab ="Mean monthly temperature (°C)", ylab = "Monthly lag", cex.lab=1.3, cex.axis=1.5,main="A"))
plot(mo.pred.temp.iNTS, xlab="Monthly lag (given mean 19 °C/month)", "slices", var=c(19), col="orange2", ci.arg=list(col=terrain.colors(70, alpha = 1)), ci.level=0.95, ci='b',lwd=4.5, ylab="RR of iNTS", cex.lab=1.3, cex.axis=1.5,main="B")
plot(mo.pred.temp.iNTS, xlab="Monthly lag (given mean 29 °C/month)", "slices", var=c(29), col="orange2", ci.arg=list(col=terrain.colors(70, alpha = 1)), ci.level=0.95, ci='b',lwd=4.5, ylab="RR of iNTS", cex.lab=1.3, cex.axis=1.5,main="C")

#----------construct cross-basis for typhoid using optimal dfs to predict rainfall effect
sort(mo.dlnmT$rainfall_obs, decreasing = FALSE)
varknots=equalknots(mo.dlnmT$rainfall_obs, fun = "ns", df=4)
lagknots <- logknots(8, df=3)
mo.cb.rain.typhoid1 <- crossbasis(mo.dlnmT$rainfall_obs, lag =8, argvar = list(knots=varknots), arglag = list(knots=lagknots))
summary(mo.cb.rain.typhoid1)

#----------model fitting for typhoid
mo.model.typhoid1 <- glm(incid_obsX ~ mo.cb.rain.typhoid1 + month + year, family = quasipoisson(), na.action=na.exclude, mo.dlnmT)

#----------model validation check for typhoid
dev.off()
pacf(residuals(mo.model.typhoid1,type="deviance"),na.action=na.omit,main="Autocorrelation from original model",xlim=c(0,8))
mo.model.typhoid1 <- update(mo.model.typhoid1,.~.+Lag(residuals(mo.model.typhoid1,type="deviance"),4)) #add residuals at lag 1 to significantly reduce partial autocorrelation
pacf(residuals(mo.model.typhoid1,type="deviance"),na.action=na.omit,main="Autocorrelation (adjusted model)",xlim=c(0,8))

#----------validated model predictions for typhoid
mo.pred.rain.typhoid1 <- crosspred(mo.cb.rain.typhoid1, mo.model.typhoid1, cen = 0, by=0.2)

#----------plotting countour and curves for rainfall on typhoid
dev.off()
par(mar=c(5,5,2,2)+1)
plot(mo.pred.rain.typhoid1, "contour", key.title=title("typhoid"), plot.title=title("", xlab ="Mean monthly rainfall (mm)", ylab = "Monthly lag", cex.lab=1.3, cex.axis=1.5,main="D"))
plot(mo.pred.rain.typhoid1, xlab="Monthly lag (given mean 9 mm/month)", "slices",var=c(9), col="red2", ci.arg=list(col=terrain.colors(70, alpha = 1)), ci.level=0.95, ci='b',lwd=4.5, ylab="RR of typhoid", cex.lab=1.3, cex.axis=1.5,main="E")
plot(mo.pred.rain.typhoid1, xlab="Monthly lag (given mean 13 mm/month)", "slices",var=c(13), col="red2", ci.arg=list(col=terrain.colors(70, alpha = 1)), ci.level=0.95, ci='b',lwd=4.5, ylab="RR of typhoid", cex.lab=1.3, cex.axis=1.5,main="F")

#----------construct cross-basis for typhoid using optimal dfs to predict temperature effect
sort(mo.dlnmT$rainfall_obs, decreasing = FALSE)
varknots=equalknots(mo.dlnmT$rainfall_obs, fun = "ns", df=3)
lagknots <- logknots(8, df=3)
mo.cb.rain.typhoid2 <- crossbasis(mo.dlnmT$rainfall_obs, lag =8, argvar = list(knots=varknots), arglag = list(knots=lagknots))
summary(mo.cb.rain.typhoid2)

sort(mo.dlnmT$temperature_obs, decreasing = FALSE)
varknots=equalknots(mo.dlnmT$temperature_obs, fun = "ns", df=3)
lagknots <- logknots(8, df=3)
mo.cb.temp.typhoid <- crossbasis(mo.dlnmT$temperature_obs, lag =8, argvar = list(knots=varknots), arglag = list(knots=lagknots))
summary(mo.cb.temp.typhoid)

#----------model fitting for typhoid
mo.model.typhoid2 <- glm(incid_obsX ~ mo.cb.rain.typhoid2 + mo.cb.temp.typhoid + month + year, family = quasipoisson(), na.action=na.exclude, mo.dlnmT)

#----------model validation check for typhoid
dev.off()
pacf(residuals(mo.model.typhoid2,type="deviance"),na.action=na.omit,main="Autocorrelation from original model",xlim=c(0,8))

#----------validated model predictions for typhoid
mo.pred.temp.typhoid2 <- crosspred(mo.cb.temp.typhoid, mo.model.typhoid2, cen = 23, by=0.2)

#----------plotting countour and curves for temperature on typhoid
dev.off()
par(mar=c(5,5,2,2)+0.1)
plot(mo.pred.temp.typhoid2, "contour", key.title=title("typhoid"), plot.title=title("", xlab ="Mean monthly temperature (°C)", ylab = "Monthly lag", cex.lab=1.3, cex.axis=1.5,main="D"))
plot(mo.pred.temp.typhoid2, "slices", xlab="Monthly lag (given mean 19 °C/month)", var=c(19), col="red2", ci.arg=list(col=topo.colors(70, alpha = 1)), ci.level=0.95, ci='b', lwd=4.5, ylab="RR of typhoid", cex.lab=1.3, cex.axis=1.5,main="E")
plot(mo.pred.temp.typhoid2, "slices", xlab="Monthly lag (given mean 25 °C/month)", var=c(25), col="red2", ci.arg=list(col=topo.colors(70, alpha = 1)), ci.level=0.95, ci='b', lwd=4.5, ylab="RR of typhoid", cex.lab=1.3, cex.axis=1.5,main="F")
plot(mo.pred.temp.typhoid2, "slices", xlab="Monthly lag (given mean 29 °C/month)", var=c(29), col="red2", ci.arg=list(col=topo.colors(70, alpha = 1)), ci.level=0.95, ci='b', lwd=4.5, ylab="RR of typhoid", cex.lab=1.3, cex.axis=1.5,main="G")

#----------plot iNTS model diagnostic plots (rainfall+temperature)
dev.off()
par(mfrow=c(3,4))
plot(mo.dlnmN$date, residuals(mo.model.iNTS,type="deviance"), pch=19, cex=0.8, col=grey(0.6),main="A", ylab="Residuals",xlab="",cex.lab=1.2,cex.axis=1.2) 
abline(h=0,lty=2,lwd=3)
pacf(residuals(mo.model.iNTS,type="deviance"),na.action=na.omit,main="B",xlim=c(0,8),xlab="", cex.lab=1.2,cex.axis=1.2)
acf(residuals(mo.model.iNTS,type="deviance"),na.action=na.omit,main="C",xlim=c(0,8),xlab="", cex.lab=1.2,cex.axis=1.2)
plot(matrix(mo.dlnmN$incid_obsX),matrix(predict(mo.model.iNTS, type="response")), col=c("orange2","gray30"),main="D",xlab="Predicted iNTS incidence", ylab="Observed iNTS incidence",pch=19,cex=0.8,cex.lab=1.2,cex.axis=1.2)
legend("topleft", legend=c("observed", "predicted"), col=c("grey30", "orange2"), cex=0.8, pch=19)

plot(mo.dlnmT$date, residuals(mo.model.typhoid1,type="deviance"), pch=19, cex=0.8, col=grey(0.6),main="E", ylab="Residuals",xlab="",cex.lab=1.2,cex.axis=1.2)
abline(h=0,lty=2,lwd=3)
pacf(residuals(mo.model.typhoid1,type="deviance"),na.action=na.omit,main="F",xlim=c(0,8),xlab="",cex.lab=1.2, cex.axis=1.2)
acf(residuals(mo.model.typhoid1,type="deviance"),na.action=na.omit,main="G",xlim=c(0,8),xlab="",cex.lab=1.2, cex.axis=1.2)
plot(matrix(mo.dlnmT$incid_obsX),matrix(predict(mo.model.typhoid1, type="response")), col=c("red2","gray30"),main="H",xlab="Predicted typhoid incidence", ylab="Observed typhoid incidence",pch=19,cex=0.8,cex.lab=1.2,cex.axis=1.2)
legend("topleft", legend=c("observed", "predicted"), col=c("grey30", "red2"), cex=0.8, pch=19)

plot(mo.dlnmT$date, residuals(mo.model.typhoid2,type="deviance"), pch=19, cex=0.8, col=grey(0.6),main="I", ylab="Residuals", xlab="Year",cex.lab=1.2,cex.axis=1.2) 
abline(h=0,lty=2,lwd=3)
pacf(residuals(mo.model.typhoid2,type="deviance"),na.action=na.omit,main="J",xlim=c(0,8),xlab="Lag (month)",cex.lab=1.2, cex.axis=1.2)
acf(residuals(mo.model.typhoid2,type="deviance"),na.action=na.omit,main="K",xlim=c(0,8),xlab="Lag (month)",cex.lab=1.2, cex.axis=1.2)
plot(matrix(mo.dlnmT$incid_obsX),matrix(predict(mo.model.typhoid2, type="response")), col=c("red2","gray30"),main="L",xlab="Predicted typhoid incidence", ylab="Observed typhoid incidence",pch=19,cex=0.8,cex.lab=1.2,cex.axis=1.2)
legend("topleft", legend=c("observed", "predicted"), col=c("grey30", "red2"), cex=0.8, pch=19)

#----------Plots of predictions due to rainfall/temperature
RRTable <- data.frame(rbind(
mo.pred.rain.iNTS$matRRfit["9",], mo.pred.rain.iNTS$matRRlow["9",], mo.pred.rain.iNTS$matRRhigh["9",],
mo.pred.rain.iNTS$matRRfit["13",],mo.pred.rain.iNTS$matRRlow["13",],mo.pred.rain.iNTS$matRRhigh["13",],
mo.pred.temp.iNTS$matRRfit["19",],mo.pred.temp.iNTS$matRRlow["19",],mo.pred.temp.iNTS$matRRhigh["19",],
mo.pred.temp.iNTS$matRRfit["29",],mo.pred.temp.iNTS$matRRlow["29",],mo.pred.temp.iNTS$matRRhigh["29",],
mo.pred.rain.typhoid1$matRRfit["9",],mo.pred.rain.typhoid1$matRRlow["9",],mo.pred.rain.typhoid1$matRRhigh["9",],
mo.pred.rain.typhoid1$matRRfit["13",],mo.pred.rain.typhoid1$matRRlow["13",],mo.pred.rain.typhoid1$matRRhigh["13",],
mo.pred.temp.typhoid2$matRRfit["19",],mo.pred.temp.typhoid2$matRRlow["19",],mo.pred.temp.typhoid2$matRRhigh["19",],
mo.pred.temp.typhoid2$matRRfit["25",],mo.pred.temp.typhoid2$matRRlow["25",],mo.pred.temp.typhoid2$matRRhigh["25",]
))
kable(RRTable)

#----------vary degrees of freedom for sensitivity analyses of NTS/typhoid predictions
nts.cb.rain <- list()
nts.cb.temp <- list()
typ.cb.rain <- list()
typ.cb.temp <- list()
sort(mo.dlnmN$rainfall_obs, decreasing=FALSE)
sort(mo.dlnmN$temperature_obs, decreasing=FALSE)
sort(mo.dlnmT$rainfall_obs, decreasing=FALSE)
sort(mo.dlnmT$temperature_obs, decreasing=FALSE)
l=1
for(k in 3:5){
  for(j in 3:5){
    for(i in 3:5){
      nts.lagknots <- logknots(8, fun="ns", df=i)
      nts.varknots.r=equalknots(mo.dlnmN$rainfall_obs, fun="ns", df=j)
      nts.varknots.t=equalknots(mo.dlnmN$temperature_obs, fun="ns", df=k)
      nts.cb.rain[[l]] <- crossbasis(mo.dlnmN$rainfall_obs, lag=8, argvar=list(fun="ns", knots=nts.varknots.r), arglag=list(knots=nts.lagknots))
      nts.cb.temp[[l]] <- crossbasis(mo.dlnmN$temperature_obs, lag=8, argvar=list(fun="ns", knots=nts.varknots.t), arglag=list(knots=nts.lagknots))
      
      typ.lagknots <- logknots(8, fun="ns", df=i)
      typ.varknots.r=equalknots(mo.dlnmT$rainfall_obs, fun="ns", df=j)
      typ.varknots.t=equalknots(mo.dlnmT$temperature_obs, fun="ns", df=k)
      typ.cb.rain[[l]] <- crossbasis(mo.dlnmT$rainfall_obs, lag=8, argvar=list(fun="ns", knots=typ.varknots.r), arglag=list(knots=typ.lagknots))
      typ.cb.temp[[l]] <- crossbasis(mo.dlnmT$temperature_obs, lag =8, argvar = list(fun="ns", knots=typ.varknots.t), arglag=list(knots=typ.lagknots))
      
      l=l+1  
    }
  }
}

#----------formulate alternative sensitivity models and plot predictions of NTS (Supplementay Figure S5)
nts.cb.rain.s1 <- nts.cb.rain[[2]]; nts.cb.temp.s1 <- nts.cb.temp[[2]]
nts.cb.rain.s2 <- nts.cb.rain[[3]]; nts.cb.temp.s2 <- nts.cb.temp[[3]]
nts.cb.rain.s3 <- nts.cb.rain[[4]]; nts.cb.temp.s3 <- nts.cb.temp[[4]]
nts.cb.rain.s4 <- nts.cb.rain[[7]]; nts.cb.temp.s4 <- nts.cb.temp[[7]]
nts.cb.rain.s5 <- nts.cb.rain[[10]]; nts.cb.temp.s5 <- nts.cb.temp[[10]]
nts.cb.rain.s6 <- nts.cb.rain[[19]]; nts.cb.temp.s6 <- nts.cb.temp[[19]]

nts.model.s1 <- glm(incid_obsX~nts.cb.rain.s1 + nts.cb.temp.s1 + month + year, family=quasipoisson(), na.action=na.exclude, mo.dlnmN)
nts.model.s2 <- glm(incid_obsX~nts.cb.rain.s2 + nts.cb.temp.s2 + month + year, family=quasipoisson(), na.action=na.exclude, mo.dlnmN)
nts.model.s3 <- glm(incid_obsX~nts.cb.rain.s3 + nts.cb.temp.s3 + month + year, family=quasipoisson(), na.action=na.exclude, mo.dlnmN)
nts.model.s4 <- glm(incid_obsX~nts.cb.rain.s4 + nts.cb.temp.s4 + month + year, family=quasipoisson(), na.action=na.exclude, mo.dlnmN)
nts.model.s5 <- glm(incid_obsX~nts.cb.rain.s5 + nts.cb.temp.s5 + month + year, family=quasipoisson(), na.action=na.exclude, mo.dlnmN)
nts.model.s6 <- glm(incid_obsX~nts.cb.rain.s6 + nts.cb.temp.s6 + month + year, family=quasipoisson(), na.action=na.exclude, mo.dlnmN)

nts.model.s1 <- update(nts.model.s1,.~.+Lag(residuals(nts.model.s1,type="deviance"),1)) #add residuals at lag 1
nts.model.s2 <- update(nts.model.s2,.~.+Lag(residuals(nts.model.s2,type="deviance"),1)) #to correct for partial
nts.model.s3 <- update(nts.model.s3,.~.+Lag(residuals(nts.model.s3,type="deviance"),1)) #autocorrelation
nts.model.s4 <- update(nts.model.s4,.~.+Lag(residuals(nts.model.s4,type="deviance"),1)) 
nts.model.s5 <- update(nts.model.s5,.~.+Lag(residuals(nts.model.s5,type="deviance"),1)) 
nts.model.s6 <- update(nts.model.s6,.~.+Lag(residuals(nts.model.s6,type="deviance"),1)) 

nts.pred.rain.s1 <- crosspred(nts.cb.rain.s1, nts.model.s1, cen=0, by=0.2)
nts.pred.rain.s2 <- crosspred(nts.cb.rain.s2, nts.model.s2, cen=0, by=0.2)
nts.pred.rain.s3 <- crosspred(nts.cb.rain.s3, nts.model.s3, cen=0, by=0.2)
nts.pred.rain.s4 <- crosspred(nts.cb.rain.s4, nts.model.s4, cen=0, by=0.2)
nts.pred.rain.s5 <- crosspred(nts.cb.rain.s5, nts.model.s5, cen=0, by=0.2)
nts.pred.rain.s6 <- crosspred(nts.cb.rain.s6, nts.model.s6, cen=0, by=0.2)

nts.pred.temp.s1 <- crosspred(nts.cb.temp.s1, nts.model.s1, cen=23, by=0.2)
nts.pred.temp.s2 <- crosspred(nts.cb.temp.s2, nts.model.s2, cen=23, by=0.2)
nts.pred.temp.s3 <- crosspred(nts.cb.temp.s3, nts.model.s3, cen=23, by=0.2)
nts.pred.temp.s4 <- crosspred(nts.cb.temp.s4, nts.model.s4, cen=23, by=0.2)
nts.pred.temp.s5 <- crosspred(nts.cb.temp.s5, nts.model.s5, cen=23, by=0.2)
nts.pred.temp.s6 <- crosspred(nts.cb.temp.s6, nts.model.s6, cen=23, by=0.2)

par(mfrow=c(2,6))
plot(nts.pred.rain.s1,"slices",var=c(9),col="orange2",ci.arg=list(col=topo.colors(70,alpha=1)),ci.level=0.95,ci='b',lwd=4.5,xlab="Monthly lag (mean 9 mm/M)",ylab="RR of iNTS",cex.lab=1.1,cex.axis=1.1,main="A")
plot(nts.pred.rain.s2,"slices",var=c(9),col="orange2",ci.arg=list(col=topo.colors(70,alpha=1)),ci.level=0.95,ci='b',lwd=4.5,xlab="Monthly lag (mean 9 mm/M)",ylab="RR of iNTS",cex.lab=1.1,cex.axis=1.1,main="C")
plot(nts.pred.rain.s3,"slices",var=c(9),col="orange2",ci.arg=list(col=topo.colors(70,alpha=1)),ci.level=0.95,ci='b',lwd=4.5,xlab="Monthly lag (mean 9 mm/M)",ylab="RR of iNTS",cex.lab=1.1,cex.axis=1.1,main="E")
plot(nts.pred.rain.s4,"slices",var=c(9),col="orange2",ci.arg=list(col=topo.colors(70,alpha=1)),ci.level=0.95,ci='b',lwd=4.5,xlab="Monthly lag (mean 9 mm/M)",ylab="RR of iNTS",cex.lab=1.1,cex.axis=1.1,main="G")
plot(nts.pred.rain.s5,"slices",var=c(9),col="orange2",ci.arg=list(col=topo.colors(70,alpha=1)),ci.level=0.95,ci='b',lwd=4.5,xlab="Monthly lag (mean 9 mm/M)",ylab="RR of iNTS",cex.lab=1.1,cex.axis=1.1,main="I")
plot(nts.pred.rain.s6,"slices",var=c(9),col="orange2",ci.arg=list(col=topo.colors(70,alpha=1)),ci.level=0.95,ci='b',lwd=4.5,xlab="Monthly lag (mean 9 mm/M)",ylab="RR of iNTS",cex.lab=1.1,cex.axis=1.1,main="K")

plot(nts.pred.temp.s1,"slices",var=c(19),col="orange2",ci.arg=list(col=terrain.colors(70,alpha=1)),ci.level=0.95,ci='b',lwd=4.5,xlab="Monthly lag (mean 19 °C/M)",ylab="RR of iNTS",cex.lab=1.1,cex.axis=1.1,main="B")
plot(nts.pred.temp.s2,"slices",var=c(19),col="orange2",ci.arg=list(col=terrain.colors(70,alpha=1)),ci.level=0.95,ci='b',lwd=4.5,xlab="Monthly lag (mean 19 °C/M)",ylab="RR of iNTS",cex.lab=1.1,cex.axis=1.1,main="D")
plot(nts.pred.temp.s3,"slices",var=c(19),col="orange2",ci.arg=list(col=terrain.colors(70,alpha=1)),ci.level=0.95,ci='b',lwd=4.5,xlab="Monthly lag (mean 19 °C/M)",ylab="RR of iNTS",cex.lab=1.1,cex.axis=1.1,main="F")
plot(nts.pred.temp.s4,"slices",var=c(19),col="orange2",ci.arg=list(col=terrain.colors(70,alpha=1)),ci.level=0.95,ci='b',lwd=4.5,xlab="Monthly lag (mean 19 °C/M)",ylab="RR of iNTS",cex.lab=1.1,cex.axis=1.1,main="H")
plot(nts.pred.temp.s5,"slices",var=c(19),col="orange2",ci.arg=list(col=terrain.colors(70,alpha=1)),ci.level=0.95,ci='b',lwd=4.5,xlab="Monthly lag (mean 19 °C/M)",ylab="RR of iNTS",cex.lab=1.1,cex.axis=1.1,main="J")
plot(nts.pred.temp.s6,"slices",var=c(19),col="orange2",ci.arg=list(col=terrain.colors(70,alpha=1)),ci.level=0.95,ci='b',lwd=4.5,xlab="Monthly lag (mean 19 °C/M)",ylab="RR of iNTS",cex.lab=1.1,cex.axis=1.1,main="L")

#----------formulate alternative sensitivity models and plot predictions of Typhoid (Supplementay Figure S6)
typ.cb.rain.s1 <- typ.cb.rain[[2]]; typ.cb.temp.s1 <- typ.cb.temp[[2]]
typ.cb.rain.s2 <- typ.cb.rain[[3]]; typ.cb.temp.s2 <- typ.cb.temp[[3]]
typ.cb.rain.s3 <- typ.cb.rain[[4]]; typ.cb.temp.s3 <- typ.cb.temp[[4]]
typ.cb.rain.s4 <- typ.cb.rain[[7]]; typ.cb.temp.s4 <- typ.cb.temp[[7]]
typ.cb.rain.s5 <- typ.cb.rain[[10]]; typ.cb.temp.s5 <- typ.cb.temp[[10]]
typ.cb.rain.s6 <- typ.cb.rain[[19]]; typ.cb.temp.s6 <- typ.cb.temp[[19]]

typ.modelR.s1 <- glm(incid_obsX~typ.cb.rain.s1 + month + year, family=quasipoisson(), na.action=na.exclude, mo.dlnmT)
typ.modelR.s2 <- glm(incid_obsX~typ.cb.rain.s2 + month + year, family=quasipoisson(), na.action=na.exclude, mo.dlnmT)
typ.modelR.s3 <- glm(incid_obsX~typ.cb.rain.s3 + month + year, family=quasipoisson(), na.action=na.exclude, mo.dlnmT)
typ.modelR.s4 <- glm(incid_obsX~typ.cb.rain.s4 + month + year, family=quasipoisson(), na.action=na.exclude, mo.dlnmT)
typ.modelR.s5 <- glm(incid_obsX~typ.cb.rain.s5 + month + year, family=quasipoisson(), na.action=na.exclude, mo.dlnmT)
typ.modelR.s6 <- glm(incid_obsX~typ.cb.rain.s6 + month + year, family=quasipoisson(), na.action=na.exclude, mo.dlnmT)
typ.modelT.s1 <- glm(incid_obsX~typ.cb.rain.s1 + typ.cb.temp.s1 + month + year, family=quasipoisson(), na.action=na.exclude, mo.dlnmT)
typ.modelT.s2 <- glm(incid_obsX~typ.cb.rain.s2 + typ.cb.temp.s2 + month + year, family=quasipoisson(), na.action=na.exclude, mo.dlnmT)
typ.modelT.s3 <- glm(incid_obsX~typ.cb.rain.s3 + typ.cb.temp.s3 + month + year, family=quasipoisson(), na.action=na.exclude, mo.dlnmT)
typ.modelT.s4 <- glm(incid_obsX~typ.cb.rain.s4 + typ.cb.temp.s4 + month + year, family=quasipoisson(), na.action=na.exclude, mo.dlnmT)
typ.modelT.s5 <- glm(incid_obsX~typ.cb.rain.s5 + typ.cb.temp.s5 + month + year, family=quasipoisson(), na.action=na.exclude, mo.dlnmT)
typ.modelT.s6 <- glm(incid_obsX~typ.cb.rain.s6 + typ.cb.temp.s6 + month + year, family=quasipoisson(), na.action=na.exclude, mo.dlnmT)

typ.modelR.s1 <- update(typ.modelR.s1,.~.+Lag(residuals(typ.modelR.s1,type="deviance"),1))
typ.modelR.s1 <- update(typ.modelR.s1,.~.+Lag(residuals(typ.modelR.s1,type="deviance"),4))
typ.modelR.s2 <- update(typ.modelR.s2,.~.+Lag(residuals(typ.modelR.s2,type="deviance"),1)) 
typ.modelR.s5 <- update(typ.modelR.s5,.~.+Lag(residuals(typ.modelR.s5,type="deviance"),1)) 
typ.modelR.s5 <- update(typ.modelR.s5,.~.+Lag(residuals(typ.modelR.s5,type="deviance"),4))
typ.modelR.s6 <- update(typ.modelR.s6,.~.+Lag(residuals(typ.modelR.s6,type="deviance"),1)) 
typ.modelR.s6 <- update(typ.modelR.s6,.~.+Lag(residuals(typ.modelR.s6,type="deviance"),4)) 

typ.modelT.s1 <- update(typ.modelT.s1,.~.+Lag(residuals(typ.modelT.s1,type="deviance"),3))
typ.modelT.s1 <- update(typ.modelT.s1,.~.+Lag(residuals(typ.modelT.s1,type="deviance"),1))
typ.modelT.s1 <- update(typ.modelT.s1,.~.+Lag(residuals(typ.modelT.s1,type="deviance"),2))
typ.modelT.s3 <- update(typ.modelT.s3,.~.+Lag(residuals(typ.modelT.s3,type="deviance"),2))
typ.modelT.s3 <- update(typ.modelT.s3,.~.+Lag(residuals(typ.modelT.s3,type="deviance"),1)) 
typ.modelT.s4 <- update(typ.modelT.s4,.~.+Lag(residuals(typ.modelT.s4,type="deviance"),2)) 
typ.modelT.s5 <- update(typ.modelT.s5,.~.+Lag(residuals(typ.modelT.s5,type="deviance"),1)) 
typ.modelT.s6 <- update(typ.modelT.s6,.~.+Lag(residuals(typ.modelT.s6,type="deviance"),3)) 
typ.modelT.s6 <- update(typ.modelT.s6,.~.+Lag(residuals(typ.modelT.s6,type="deviance"),1)) 
typ.modelT.s6 <- update(typ.modelT.s6,.~.+Lag(residuals(typ.modelT.s6,type="deviance"),2))

typ.pred.rain.s1 <- crosspred(typ.cb.rain.s1, typ.modelR.s1, cen=0, by=0.2)
typ.pred.rain.s2 <- crosspred(typ.cb.rain.s2, typ.modelR.s2, cen=0, by=0.2)
typ.pred.rain.s3 <- crosspred(typ.cb.rain.s3, typ.modelR.s3, cen=0, by=0.2)
typ.pred.rain.s4 <- crosspred(typ.cb.rain.s4, typ.modelR.s4, cen=0, by=0.2)
typ.pred.rain.s5 <- crosspred(typ.cb.rain.s5, typ.modelR.s5, cen=0, by=0.2)
typ.pred.rain.s6 <- crosspred(typ.cb.rain.s6, typ.modelR.s6, cen=0, by=0.2)

typ.pred.temp.s1 <- crosspred(typ.cb.temp.s1, typ.modelT.s1, cen=23, by=0.2)
typ.pred.temp.s2 <- crosspred(typ.cb.temp.s2, typ.modelT.s2, cen=23, by=0.2)
typ.pred.temp.s3 <- crosspred(typ.cb.temp.s3, typ.modelT.s3, cen=23, by=0.2)
typ.pred.temp.s4 <- crosspred(typ.cb.temp.s4, typ.modelT.s4, cen=23, by=0.2)
typ.pred.temp.s5 <- crosspred(typ.cb.temp.s5, typ.modelT.s5, cen=23, by=0.2)
typ.pred.temp.s6 <- crosspred(typ.cb.temp.s6, typ.modelT.s6, cen=23, by=0.2)

par(mfrow=c(2,6))
plot(typ.pred.rain.s1,"slices",var=c(9),col="red2",ci.arg=list(col=topo.colors(70,alpha=1)),ci.level=0.95,ci='b',lwd=4.5,xlab="Monthly lag (mean 9 mm/M)",ylab="RR of typhoid",cex.lab=1.1,cex.axis=1.1,main="A")
plot(typ.pred.rain.s2,"slices",var=c(9),col="red2",ci.arg=list(col=topo.colors(70,alpha=1)),ci.level=0.95,ci='b',lwd=4.5,xlab="Monthly lag (mean 9 mm/M)",ylab="RR of typhoid",cex.lab=1.1,cex.axis=1.1,main="C")
plot(typ.pred.rain.s3,"slices",var=c(9),col="red2",ci.arg=list(col=topo.colors(70,alpha=1)),ci.level=0.95,ci='b',lwd=4.5,xlab="Monthly lag (mean 9 mm/M)",ylab="RR of typhoid",cex.lab=1.1,cex.axis=1.1,main="E")
plot(typ.pred.rain.s4,"slices",var=c(9),col="red2",ci.arg=list(col=topo.colors(70,alpha=1)),ci.level=0.95,ci='b',lwd=4.5,xlab="Monthly lag (mean 9 mm/M)",ylab="RR of typhoid",cex.lab=1.1,cex.axis=1.1,main="G")
plot(typ.pred.rain.s5,"slices",var=c(9),col="red2",ci.arg=list(col=topo.colors(70,alpha=1)),ci.level=0.95,ci='b',lwd=4.5,xlab="Monthly lag (mean 9 mm/M)",ylab="RR of typhoid",cex.lab=1.1,cex.axis=1.1,main="I")
plot(typ.pred.rain.s6,"slices",var=c(9),col="red2",ci.arg=list(col=topo.colors(70,alpha=1)),ci.level=0.95,ci='b',lwd=4.5,xlab="Monthly lag (mean 9 mm/M)",ylab="RR of typhoid",cex.lab=1.1,cex.axis=1.1,main="K")

plot(typ.pred.temp.s1,"slices",var=c(19),col="red2",ci.arg=list(col=terrain.colors(70,alpha=1)),ci.level=0.95,ci='b',lwd=4.5,xlab="Month-lag (mean 19 °C/M)",ylab="RR oftyphoid",cex.lab=1.1,cex.axis=1.1,main="B")
plot(typ.pred.temp.s2,"slices",var=c(19),col="red2",ci.arg=list(col=terrain.colors(70,alpha=1)),ci.level=0.95,ci='b',lwd=4.5,xlab="Month-lag (mean 19 °C/M)",ylab="RR of typhoid",cex.lab=1.1,cex.axis=1.1,main="D")
plot(typ.pred.temp.s3,"slices",var=c(19),col="red2",ci.arg=list(col=terrain.colors(70,alpha=1)),ci.level=0.95,ci='b',lwd=4.5,xlab="Month-lag (mean 19 °C/M)",ylab="RR of typhoid",cex.lab=1.1,cex.axis=1.1,main="F")
plot(typ.pred.temp.s4,"slices",var=c(19),col="red2",ci.arg=list(col=terrain.colors(70,alpha=1)),ci.level=0.95,ci='b',lwd=4.5,xlab="Month-lag (mean 19 °C/M)",ylab="RR of typhoid",cex.lab=1.1,cex.axis=1.1,main="H")
plot(typ.pred.temp.s5,"slices",var=c(19),col="red2",ci.arg=list(col=terrain.colors(70,alpha=1)),ci.level=0.95,ci='b',lwd=4.5,xlab="Month-lag (mean 19 °C/M)",ylab="RR of typhoid",cex.lab=1.1,cex.axis=1.1,main="J")
plot(typ.pred.temp.s6,"slices",var=c(19),col="red2",ci.arg=list(col=terrain.colors(70,alpha=1)),ci.level=0.95,ci='b',lwd=4.5,xlab="Month-lag (mean 19 °C/M)",ylab="RR of typhoid",cex.lab=1.1,cex.axis=1.1,main="L")

#----------descriptive stats of the study population seasonal-unadjusted
ci(mo.dlnmN$incid_obsX)
ci(mo.dlnmT$incid_obsX)

#END SCRIPT.
