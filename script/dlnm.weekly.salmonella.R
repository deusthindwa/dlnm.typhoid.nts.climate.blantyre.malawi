#02/08/2018-07/11/2019
#by Deus Thindwa

#----------Install required packages.
dlnm.analysis.packages <- c("tidyverse", "lubridate","xts","ggthemes","PerformanceAnalytics","reshape2","rugarch","timetk","parallel","timeSeries","tseries","data.table","ggplot2","dlnm","broom","caret","gridExtra","splines","splines2","pspline","cowplot","mgcv","spi","chron","gridGraphics","grid","pscl","MASS", "AER", "Hmisc", "MuMIn", "VGAM", "forecast", "seasonal", "plotly", "ggmap", "rgeos", "tmap", "maptools", "maps", "ggfortify", "htmltools","webshot","knitr","flexdashboard", "imager", "httr", "gmodels", "curl", "here")

#----------load required packages.
lapply(dlnm.analysis.packages, library, character.only=TRUE)

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
case.typhi <- apply.weekly(case.typhi, FUN = sum)
case.iNTS <- apply.weekly(case.iNTS, FUN = sum)

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

climate.rain <- apply.weekly(climate.rain, FUN = mean)
climate.temp <- apply.weekly(climate.temp, FUN = mean)

#----------create tibbles for case and climate seasonal-adjustment
case.typhi <-tk_tbl(case.typhi, preserve_index = TRUE, rename_index = "date") 
case.iNTS <-tk_tbl(case.iNTS, preserve_index = TRUE, rename_index = "date") 
climate.rain <-tk_tbl(climate.rain, preserve_index = TRUE, rename_index = "date") 
climate.temp <-tk_tbl(climate.temp, preserve_index = TRUE, rename_index = "date") 

#----------seasonally-adjusted cases and climate
case.iNTS.ts <- ts(na.omit(case.iNTS$case_count), frequency = 53)
case.typhi.ts <- ts(na.omit(case.typhi$case_count), frequency = 53)
climate.rain.ts <- ts(na.omit(climate.rain$rainfall), frequency = 53)
climate.temp.ts <- ts(na.omit(climate.temp$temperature), frequency = 53)

#----------plot decomposed all series
dev.off()
x<-mstl(case.iNTS.ts,s.window="period") %>% ggfortify:::autoplot.ts(main="A",xlab="Years (2000-2015)",ylab="Number of iNTS cases",size=1,colour="orange2",is.date=FALSE) + theme_bw()
y<-mstl(case.typhi.ts,s.window="periodic") %>% ggfortify:::autoplot.ts(main="B",xlab="Years (2000-2015)",ylab="Number of typhoid cases",size=1,colour="red2",is.date=FALSE) + theme_bw()
z<-mstl(climate.rain.ts,s.window="periodic") %>% ggfortify:::autoplot.ts(main="C",xlab="Years (2000-2015)",ylab="Rainfall (mm)",size=1,colour="blue2",is.date=FALSE) + theme_bw()
v<-mstl(climate.temp.ts,s.window="periodic") %>% ggfortify:::autoplot.ts(main="D",xlab="Years (2000-2015)",ylab="Temperature (°C)",size=1,colour="green2",is.date=FALSE) + theme_bw()
grid.arrange(grobs=list(x,y,z,v), ncol=4, nrow=1)

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
case.iNTS$incid_obs <-case.iNTS$case_count*100000/case.iNTS$census 

#----------calculated incidence of monthly typhoid
case.typhi$census[year(case.typhi$date)==2000]<- 852054; case.typhi$census[year(case.typhi$date)==2001]<-873382
case.typhi$census[year(case.typhi$date)==2002]<- 894710; case.typhi$census[year(case.typhi$date)==2003]<-916039
case.typhi$census[year(case.typhi$date)==2004]<- 937367; case.typhi$census[year(case.typhi$date)==2005]<-958695
case.typhi$census[year(case.typhi$date)==2006]<- 980023; case.typhi$census[year(case.typhi$date)==2007]<-1001352
case.typhi$census[year(case.typhi$date)==2008]<-1022680; case.typhi$census[year(case.typhi$date)==2009]<-1044008
case.typhi$census[year(case.typhi$date)==2010]<-1065337; case.typhi$census[year(case.typhi$date)==2011]<-1086665
case.typhi$census[year(case.typhi$date)==2012]<-1107993; case.typhi$census[year(case.typhi$date)==2013]<-1129322
case.typhi$census[year(case.typhi$date)==2014]<-1150650; case.typhi$census[year(case.typhi$date)==2015]<-1171978
case.typhi$incid_obs <-case.typhi$case_count*100000/case.typhi$census 

#----------boxplots of weekly and month seasonal dynamics of iNTS and typhoid
j <- seq(from=1, to=53, by=4)
i <- seq(from=1, to=12, by=1)

pox1 <- ggplot(subset(case.iNTS, year(date)<2011), aes(x=week(date), y=incid_obs))  + 
  geom_boxplot(aes( group=week(date)),color="orange2", fill="orange2", alpha=0.2) + 
  labs(title="A",x="Week (Jan-Dec)", y = "Weekly iNTS incidence") + 
  scale_x_discrete(limits = j) + 
  ylim(c(0,25)) +
  theme(axis.title.y = element_text(size = 11)) + 
  theme(axis.title.x = element_text(size = 11)) + 
  theme(axis.text.x = element_text(face="bold", size=11), axis.text.y = element_text(face="bold", size=11)) 

pox2 <- ggplot(subset(case.typhi, year(date)>2010), aes(x=week(date), y=incid_obs))  + 
  geom_boxplot(aes( group=week(date)), color="red2", fill="red2", alpha=0.2) + 
  labs(title="C", x="Week (Jan-Dec)", y = "Weekly typhoid incidence") + 
  scale_x_discrete(limits = j) + 
  ylim(c(0,25)) +
  theme(axis.title.y = element_text(size = 11)) + 
  theme(axis.title.x = element_text(size = 11)) + 
  theme(axis.text.x = element_text(face="bold", size=11), axis.text.y = element_text(face="bold", size=11)) 

grid.arrange(grobs=list(pox1, pox2), ncol=2, nrow=1)

#----------contour plots of (un)seasonal dynamics of iNTS and typhoid
case.iNTS.spi <- subset(case.iNTS, year(case.iNTS$date)<2011, select=c(date,incid_obs)) 
case.iNTS.spi$week <- week(case.iNTS.spi$date)
case.iNTS.spi$year <- year(case.iNTS.spi$date)
case.iNTS.spi$date <- NULL
case.iNTS.spi <- spread(case.iNTS.spi, year, incid_obs)
case.iNTS.spi <- as.matrix(case.iNTS.spi)[,-1]
plot_ly(x = c(2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010), z = ~case.iNTS.spi, type = "contour", colorscale = 'heatmap', contours = list(showlabels = TRUE)) %>% 
colorbar(title = "Seasonal-unadjusted \n iNTS incidence per \n 100,000 population") %>%
layout(title="<b>A</b>", xaxis=list(title ="Year",color="black"), yaxis=list(title="Week (Jan-Dec)",color="black"), font=list(size = 13))

climate.rain.spin <- subset(climate.rain, year(climate.rain$date)<2011, select=c(date,rainfall)) 
climate.rain.spin$week <- week(climate.rain.spin$date)
climate.rain.spin$year <- year(climate.rain.spin$date)
climate.rain.spin$date <- NULL
climate.rain.spin <- spread(climate.rain.spin, year, rainfall)
climate.rain.spin <- as.matrix(climate.rain.spin)[,-1]
plot_ly(x = c(2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010), z = ~climate.rain.spin, type = "contour", colorscale = 'Earth', contours = list(showlabels = TRUE)) %>% 
colorbar(title = "Seasonal-unadjusted \n Rainfall (mm)") %>%
layout(title="<b>B</b>", xaxis=list(title ="Year",color="black"), yaxis=list(title="Week (Jan-Dec)",color="black"), font=list(size = 13))

climate.temp.spin <- subset(climate.temp, year(climate.temp$date)<2011, select=c(date,temperature)) 
climate.temp.spin$week <- week(climate.temp.spin$date)
climate.temp.spin$year <- year(climate.temp.spin$date)
climate.temp.spin$date <- NULL
climate.temp.spin <- spread(climate.temp.spin, year, temperature)
climate.temp.spin <- as.matrix(climate.temp.spin)[,-1]
plot_ly(x = c(2000,2001,2002,2003,2004,2005,2006,2007,2008,2009,2010), z = ~climate.temp.spin, type = "contour", colorscale = 'Viridis', contours = list(showlabels = TRUE)) %>% 
colorbar(title = "Seasonal-unadjusted \n Temperature (°C)") %>%
layout(title="<b>C</b>", xaxis=list(title ="Year",color="black"), yaxis=list(title="Week (Jan-Dec)",color="black"), font=list(size = 13))

case.typhi.spi <- subset(case.typhi, year(case.typhi$date)>2010, select=c(date,incid_obs)) 
case.typhi.spi$week <- week(case.typhi.spi$date)
case.typhi.spi$year <- year(case.typhi.spi$date)
case.typhi.spi$date <- NULL
case.typhi.spi <- spread(case.typhi.spi, year, incid_obs)
case.typhi.spi <- as.matrix(case.typhi.spi)[,-1]
plot_ly(x = c(2011,2012,2013,2014,2015), y = c("1","2","3","4","5","6","7","8","9","10","11","12"), z = ~case.typhi.spi, type = "contour", colorscale = 'heatmap', contours = list(showlabels = TRUE)) %>% 
colorbar(title = "Seasonal-unadjusted \n typhoid incidence per \n 100,000 population") %>%
layout(title="<b>D</b>", xaxis=list(title ="Year",color="black"), yaxis=list(title="Week (Jan-Dec)",color="black"), font=list(size = 13))

climate.rain.spit <- subset(climate.rain, year(climate.rain$date)>2010, select=c(date,rainfall)) 
climate.rain.spit$week <- week(climate.rain.spit$date)
climate.rain.spit$year <- year(climate.rain.spit$date)
climate.rain.spit$date <- NULL
climate.rain.spit <- spread(climate.rain.spit, year, rainfall)
climate.rain.spit <- as.matrix(climate.rain.spit)[,-1]
plot_ly(x = c(2011,2012,2013,2014,2015), z = ~climate.rain.spit, type = "contour", colorscale = 'Earth', contours = list(showlabels = TRUE)) %>% 
colorbar(title = "Seasonal-unadjusted \n Rainfall (mm)") %>%
layout(title="<b>E</b>", xaxis=list(title ="Year",color="black"), yaxis=list(title="Week (Jan-Dec)",color="black"), font=list(size = 13))

climate.temp.spit <- subset(climate.temp, year(climate.temp$date)>2010, select=c(date,temperature)) 
climate.temp.spit$week <- week(climate.temp.spit$date)
climate.temp.spit$year <- year(climate.temp.spit$date)
climate.temp.spit$date <- NULL
climate.temp.spit <- spread(climate.temp.spit, year, temperature)
climate.temp.spit <- as.matrix(climate.temp.spit)[,-1]
plot_ly(x = c(2011,2012,2013,2014,2015), z = ~climate.temp.spit, type = "contour", colorscale = 'Viridis', contours = list(showlabels = TRUE)) %>% 
colorbar(title = "Seasonal-unadjusted \n Temperature (°C)") %>%
layout(title="<b>F</b>", xaxis=list(title ="Year",color="black"), yaxis=list(title="Week (Jan-Dec)",color="black"), font=list(size = 13))

#----------prepare final monthly NTS dataset for use in DLNM
mo.dlnmN <- bind_cols(case.iNTS, climate.rain, climate.temp, id=NULL)
mo.dlnmN$date1 <- mo.dlnmN$date2 <- NULL
mo.dlnmN <- subset(mo.dlnmN, year(date) < 2011)
mo.dlnmN$time <- seq.int(from = 1, to=574, by=1)
mo.dlnmN$year <- year(mo.dlnmN$date)
mo.dlnmN$month <- month(mo.dlnmN$date)
mo.dlnmN$week <- week(mo.dlnmN$date)
mo.dlnmN$incid_obsX<-round(mo.dlnmN$incid_obs, digits = 0)

#----------prepare final monthly typhoid dataset for use in DLNM
mo.dlnmT <- bind_cols(case.typhi, climate.rain, climate.temp, id=NULL)
mo.dlnmT$date1 <- mo.dlnmT$date2 <- NULL
mo.dlnmT <- subset(mo.dlnmT, year(date) > 2010)
mo.dlnmT$time <- seq.int(from = 1, to=262, by=1)
mo.dlnmT$year <- year(mo.dlnmT$date)
mo.dlnmT$month <- month(mo.dlnmT$date)
mo.dlnmT$week <- week(mo.dlnmT$date)
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
  nts.lagknots <- logknots(5, fun="ns", df=i)
  nts.varknots.r=equalknots(mo.dlnmN$rainfall, fun="ns", df=j)
  nts.varknots.t=equalknots(mo.dlnmN$temperature, fun="ns", df=k)
  nts.mo.cb.rain <- crossbasis(mo.dlnmN$rainfall, lag=5, argvar=list(fun="ns", knots=nts.varknots.r), arglag=list(knots=nts.lagknots))
  nts.mo.cb.temp <- crossbasis(mo.dlnmN$temperature, lag=5, argvar=list(fun="ns", knots=nts.varknots.t), arglag=list(knots=nts.lagknots))
  nts.modelR <- glm(mo.dlnmN$incid_obsX ~ nts.mo.cb.rain + week + year, family=nts.quasipoisson(), na.action=na.delete, mo.dlnmN)
  nts.modelT <- glm(mo.dlnmN$incid_obsX ~ nts.mo.cb.temp + week + year, family=nts.quasipoisson(), na.action=na.delete, mo.dlnmN)
  nts.model <-  glm(mo.dlnmN$incid_obsX ~ nts.mo.cb.rain + nts.mo.cb.temp + week + year, family = nts.quasipoisson(), na.action=na.delete, mo.dlnmN)
  
  typ.lagknots <- logknots(5, fun="ns", df=i)
  typ.varknots.r=equalknots(mo.dlnmT$rainfall, fun="ns", df=j)
  typ.varknots.t=equalknots(mo.dlnmT$temperature, fun="ns", df=k)
  typ.mo.cb.rain <- crossbasis(mo.dlnmT$rainfall, lag=5, argvar=list(fun="ns", knots=typ.varknots.r), arglag=list(knots=typ.lagknots))
  typ.mo.cb.temp <- crossbasis(mo.dlnmT$temperature, lag =5, argvar = list(fun="ns", knots=typ.varknots.t), arglag=list(knots=typ.lagknots))
  typ.modelR <- glm(mo.dlnmT$incid_obsX ~ typ.mo.cb.rain + week + year, family = typ.quasipoisson(), na.action=na.delete, mo.dlnmT)
  typ.modelT <- glm(mo.dlnmT$incid_obsX ~ typ.mo.cb.temp + week + year, family = typ.quasipoisson(), na.action=na.delete, mo.dlnmT)
  typ.model <- glm(mo.dlnmT$incid_obsX ~ typ.mo.cb.rain + typ.mo.cb.temp + week + year, family = typ.quasipoisson(), na.action=na.delete, mo.dlnmT)
  
  QAICtable[l,] <- c(l,i,j,k, QAIC(nts.modelR, chat=summary(nts.modelR)$dispersion), QAIC(nts.modelT, chat=summary(nts.modelT)$dispersion), 
                     QAIC(nts.model, chat=summary(nts.model)$dispersion), QAIC(typ.modelR, chat=summary(typ.modelR)$dispersion), QAIC(typ.modelT, chat=summary(typ.modelT)$dispersion), 
                     QAIC(typ.model, chat=1))
  l=l+1 
  }
  }
}

#----------construct cross-basis for iNTS using optimal dfs to predict rainfall effect
sort(mo.dlnmN$rainfall, decreasing=FALSE)
varknots=equalknots(mo.dlnmN$rainfall, fun="ns", df=3)
lagknots <- logknots(5, fun="ns", df=3)
mo.cb.rain.iNTS <- crossbasis(mo.dlnmN$rainfall, lag=5, argvar = list(knots=varknots), arglag = list(knots=lagknots))
summary(mo.cb.rain.iNTS)

#----------construct cross-basis for iNTS using optimal dfs to predict temperature effect
sort(mo.dlnmN$temperature, decreasing=FALSE)
varknots=equalknots(mo.dlnmN$temperature, fun="ns", df=3)
lagknots <- logknots(5, fun="ns", df=3)
mo.cb.temp.iNTS <- crossbasis(mo.dlnmN$temperature, lag=5, argvar=list(knots=varknots), arglag=list(knots=lagknots))
summary(mo.cb.temp.iNTS)

#----------model fitting for iNTS
mo.model.iNTSR <- glm(incid_obsX ~  mo.cb.rain.iNTS + week + year, family = quasipoisson(), na.action=na.exclude, mo.dlnmN)
mo.model.iNTST <- glm(incid_obsX ~  mo.cb.temp.iNTS + week + year, family = quasipoisson(), na.action=na.exclude, mo.dlnmN)

#----------model validation check for iNTS
dev.off()
pacf(residuals(mo.model.iNTSR,type="deviance"),na.action=na.omit,main="Partial autocorrelation (original model)",xlim=c(0,5))
mo.model.iNTSR <- update(mo.model.iNTSR,.~.+Lag(residuals(mo.model.iNTSR,type="deviance"),1))
mo.model.iNTSR <- update(mo.model.iNTSR,.~.+Lag(residuals(mo.model.iNTSR,type="deviance"),5))
mo.model.iNTSR <- update(mo.model.iNTSR,.~.+Lag(residuals(mo.model.iNTSR,type="deviance"),4))
mo.model.iNTSR <- update(mo.model.iNTSR,.~.+Lag(residuals(mo.model.iNTSR,type="deviance"),3))
mo.model.iNTSR <- update(mo.model.iNTSR,.~.+Lag(residuals(mo.model.iNTSR,type="deviance"),2))

dev.off()
pacf(residuals(mo.model.iNTST,type="deviance"),na.action=na.omit,main="Partial autocorrelation (original model)",xlim=c(0,5))
mo.model.iNTST <- update(mo.model.iNTST,.~.+Lag(residuals(mo.model.iNTST,type="deviance"),1))
mo.model.iNTST <- update(mo.model.iNTST,.~.+Lag(residuals(mo.model.iNTST,type="deviance"),5))
mo.model.iNTST <- update(mo.model.iNTST,.~.+Lag(residuals(mo.model.iNTST,type="deviance"),4))
mo.model.iNTST <- update(mo.model.iNTST,.~.+Lag(residuals(mo.model.iNTST,type="deviance"),3))
mo.model.iNTST<- update(mo.model.iNTST,.~.+Lag(residuals(mo.model.iNTST,type="deviance"),2))

#----------validated model predictions for iNTS
mo.pred.rain.iNTSR <- crosspred(mo.cb.rain.iNTS, mo.model.iNTSR, cen = 0, by=0.2)
mo.pred.temp.iNTST <- crosspred(mo.cb.temp.iNTS, mo.model.iNTST, cen = 23, by=0.2)

#----------plotting countour and curves for rainfall on iNTS
dev.off()
par(mar=c(5,5,2,2)+0.1)
plot(mo.pred.rain.iNTSR, "contour", key.title=title("iNTS"), plot.title=title("", xlab ="Mean weekly rainfall (mm)", ylab = "Weekly lag", cex.lab=1.3, cex.axis=1.5,main="A"))
plot(mo.pred.rain.iNTSR, "slices", xlab="Weekly lag (given mean 9 mm/wk)", var=c(9), col="orange2", ci.arg=list(col=topo.colors(70, alpha = 1)), ci.level=0.95, ci='b', lwd=4.5, ylab="RR of iNTS", cex.lab=1.3, cex.axis=1.5,main="B")
plot(mo.pred.rain.iNTSR, "slices", xlab="Weekly lag (given mean 25 mm/wk)", var=c(25), col="orange2", ci.arg=list(col=topo.colors(70, alpha = 1)), ci.level=0.95, ci='b', lwd=4.5, ylab="RR of iNTS", cex.lab=1.3, cex.axis=1.5,main="C")

#----------plotting countour and curves for temperature on iNTS
dev.off()
par(mar=c(5,5,2,2)+1)
plot(mo.pred.temp.iNTST, "contour", key.title=title("iNTS"), plot.title=title("", xlab ="Mean weekly temperature (°C)", ylab = "Weekly lag", cex.lab=1.3, cex.axis=1.5,main="B"))
plot(mo.pred.temp.iNTST, xlab="Weekly lag (given mean 19 °C/wk)", "slices", var=c(19), col="orange2", ci.arg=list(col=terrain.colors(70, alpha = 1)), ci.level=0.95, ci='b',lwd=4.5, ylab="RR of iNTS", cex.lab=1.3, cex.axis=1.5,main="B")
plot(mo.pred.temp.iNTST, xlab="Weekly lag (given mean 29 °C/wk)", "slices", var=c(29), col="orange2", ci.arg=list(col=terrain.colors(70, alpha = 1)), ci.level=0.95, ci='b',lwd=4.5, ylab="RR of iNTS", cex.lab=1.3, cex.axis=1.5,main="C")

#----------construct cross-basis for typhoid using optimal dfs to predict temperature effect
sort(mo.dlnmT$rainfall, decreasing = FALSE)
varknots=equalknots(mo.dlnmT$rainfall, fun = "ns", df=4)
lagknots <- logknots(5, df=3)
mo.cb.rain.typhoid <- crossbasis(mo.dlnmT$rainfall, lag =5, argvar = list(knots=varknots), arglag = list(knots=lagknots))
summary(mo.cb.rain.typhoid)

sort(mo.dlnmT$temperature, decreasing = FALSE)
varknots=equalknots(mo.dlnmT$temperature, fun = "ns", df=4)
lagknots <- logknots(5, df=3)
mo.cb.temp.typhoid <- crossbasis(mo.dlnmT$temperature, lag=5, argvar = list(knots=varknots), arglag = list(knots=lagknots))
summary(mo.cb.temp.typhoid)

#----------model fitting for typhoid
mo.model.typhoidR <- glm(incid_obsX ~ mo.cb.rain.typhoid + month + year, family = quasipoisson(), na.action=na.exclude, mo.dlnmT)
mo.model.typhoidT <- glm(incid_obsX ~ mo.cb.temp.typhoid + month + year, family = quasipoisson(), na.action=na.exclude, mo.dlnmT)

#----------model validation check for typhoid
dev.off()
pacf(residuals(mo.model.typhoidR,type="deviance"),na.action=na.omit,main="Autocorrelation from original model",xlim=c(0,5))
mo.model.typhoidR <- update(mo.model.typhoidR,.~.+Lag(residuals(mo.model.typhoidR,type="deviance"),1))
mo.model.typhoidR <- update(mo.model.typhoidR,.~.+Lag(residuals(mo.model.typhoidR,type="deviance"),3))
mo.model.typhoidR <- update(mo.model.typhoidR,.~.+Lag(residuals(mo.model.typhoidR,type="deviance"),2))

dev.off()
pacf(residuals(mo.model.typhoidT,type="deviance"),na.action=na.omit,main="Autocorrelation from original model",xlim=c(0,5))
mo.model.typhoidT <- update(mo.model.typhoidT,.~.+Lag(residuals(mo.model.typhoidT,type="deviance"),1))
mo.model.typhoidT <- update(mo.model.typhoidT,.~.+Lag(residuals(mo.model.typhoidT,type="deviance"),3))
mo.model.typhoidT <- update(mo.model.typhoidT,.~.+Lag(residuals(mo.model.typhoidT,type="deviance"),2))

#----------validated model predictions for typhoid
mo.pred.rain.typhoidR <- crosspred(mo.cb.rain.typhoid, mo.model.typhoidR, cen = 0, by=0.2)
mo.pred.temp.typhoidT <- crosspred(mo.cb.temp.typhoid, mo.model.typhoidT, cen = 23, by=0.2)

#----------plotting countour and curves for rainfall on iNTS
dev.off()
par(mar=c(5,5,2,2)+0.1)
plot(mo.pred.rain.typhoidR, "contour", key.title=title("typhoid"), plot.title=title("", xlab ="Mean weekly rainfall (mm)", ylab = "Weekly lag", cex.lab=1.3, cex.axis=1.5,main="C"))
plot(mo.pred.rain.typhoidR, "slices", xlab="Weekly lag (given mean 9 mm/wk)", var=c(9), col="red2", ci.arg=list(col=topo.colors(70, alpha = 1)), ci.level=0.95, ci='b', lwd=4.5, ylab="RR of typhoid", cex.lab=1.3, cex.axis=1.5,main="E")
plot(mo.pred.rain.typhoidR, "slices", xlab="Weekly lag (given mean 35 mm/wk)", var=c(35), col="red2", ci.arg=list(col=topo.colors(70, alpha = 1)), ci.level=0.95, ci='b', lwd=4.5, ylab="RR of typhoid", cex.lab=1.3, cex.axis=1.5,main="F")

#----------plotting countour and curves for temperature on typhoid
dev.off()
par(mar=c(5,5,2,2)+0.1)
plot(mo.pred.temp.typhoidT, "contour", key.title=title("typhoid"), plot.title=title("", xlab ="Mean weekly temperature (°C)", ylab = "Weekly lag", cex.lab=1.3, cex.axis=1.5,main="D"))
plot(mo.pred.temp.typhoidT, "slices", xlab="Weekly lag (given mean 19 °C/wk)", var=c(19), col="red2", ci.arg=list(col=topo.colors(70, alpha = 1)), ci.level=0.95, ci='b', lwd=4.5, ylab="RR of typhoid", cex.lab=1.3, cex.axis=1.5,main="E")
plot(mo.pred.temp.typhoidT, "slices", xlab="Weekly lag (given mean 28 °C/day)", var=c(28), col="red2", ci.arg=list(col=topo.colors(70, alpha = 1)), ci.level=0.95, ci='b', lwd=4.5, ylab="RR of typhoid", cex.lab=1.3, cex.axis=1.5,main="F")
plot(mo.pred.temp.typhoidT, "slices", xlab="Weekly lag (given mean 29 °C/day)", var=c(34), col="red2", ci.arg=list(col=topo.colors(70, alpha = 1)), ci.level=0.95, ci='b', lwd=4.5, ylab="RR of typhoid", cex.lab=1.3, cex.axis=1.5,main="G")

#END SCRIPT.
