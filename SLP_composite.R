rm(list=ls())
library(raster)

#Netcdf
#slp <- brick("../ERA_Download/era_interim_moda_SLP_all.nc")
slp <- brick("../../Climate_analysis/ERA_Download/era_interim_moda_SLP_all.nc")

## Tree ring data
trDat <- read.table("../../KBP_South/KBPS_cull_gap.rwl_tabs.txt", header = TRUE)

## Decide start year and end year based on target and tree ring data
F_yr <- min(as.numeric(substr(names(slp), 2, 5)))
L_yr <- as.numeric(max(trDat$year))

#Crop data to it
trDat <-trDat[which(trDat$year >= F_yr-1 & trDat$year<= L_yr),]

## Smallest x% years, Largest x% years
quants <- quantile(trDat$ars, probs = c(0.25, 0.75))
lq_yrs <- trDat$year[which(trDat$ars<quants[1])]
uq_yrs <- trDat$year[which(trDat$ars>quants[2])]

#rWi <- trDat[order(trDat$ars, decreasing=T)[1:5],]
#rNa <- trDat[order(trDat$ars)[1:5],]

## set spatial extent for area interested in. Longitude min - max then latitude min - max
## Creating this way allows for other created rasters to recognize as an extent with 
## the proper spatial extent, otherwise have to set xmin, xmax, ymin, ymax individually.

ext <- extent(60, 200.25, -80.25, -4.50)

#Spatial crop using extent
datC <- crop(slp, ext)
datC <- datC[[which(as.numeric(substr(names(datC), 2, 5)) >= F_yr & 
                      as.numeric(substr(names(datC), 2, 5)) <= L_yr)]]
datC <- datC[[-c(1:2, (nlayers(datC)-3):nlayers(datC))]] #removes first incomplete season JF and last SON from year

##### Seasonal Indexing ######
library(chron)

yr_mo_dy <- substr(names(datC), 2, 11)
d <- as.Date(gsub(".", '/', yr_mo_dy, fixed = T)) #fix the format by replacing "." with "/"
## changed $year==12 to $mon<8 (POSIXlt indexes months from 0
## and growing season starts in SEP) 


yr_season <- paste( 1900 + # this is the base year for POSIXlt year numbering 
                      as.POSIXlt( d )$year - 
                      1*(as.POSIXlt( d )$mon<8) ,   # offset needed for grwoing season in SH
                    c('DJF', 'MAM', 'JJA', 'SON')[          # indexing from 0-based-mon
                      1+((as.POSIXlt(d)$mon+1) %/% 3)%%4] 
                    , sep="-")

## Get yearly seasonal means
datM <- stackApply(datC, yr_season, mean) #raster with mean for each season
names(datM) <- unique(yr_season)

## create seasonal means.
ssn <- substring(names(datM), 7) #season names for 30 yr seasonal mean
datMs <- stackApply(datM, ssn, mean) #create the mean
names(datMs) <- unique(ssn) #give them meaningful names
datMs <- subset(datMs, c(3,4,1,2)) #reorder because they are setup as MAM, JJA, SON, DJF


## Extract wide years from TR data; clunky right now.
# datW <- datM[[which(as.numeric(substr(names(datM), 2, 5)) ==  rWi$year[1] )]]
# datW <- addLayer(datW, datM[[which(as.numeric(substr(names(datM), 2, 5)) ==  rWi$year[2] )]],
#                  datM[[which(as.numeric(substr(names(datM), 2, 5)) ==  rWi$year[3])]],
#                  datM[[which(as.numeric(substr(names(datM), 2, 5)) ==  rWi$year[4])]],
#                  datM[[which(as.numeric(substr(names(datM), 2, 5)) ==  rWi$year[5])]])
datW <- datM[[which(as.numeric(substr(names(datM), 2, 5)) %in% lq_yrs) ]]

ssn <- substring(names(datW), 7) #season names
datWm <- stackApply(datW, ssn, mean) #seasonal mean for wide years
names(datWm) <- unique(ssn) #meaningful names

SeasonsW <- datWm - datMs #composite difference

## Second verse same as the first.
# datN <- datM[[which(as.numeric(substr(names(datM), 2, 5)) ==  rNa$year[1] )]]
# datN <- addLayer(datN, datM[[which(as.numeric(substr(names(datM), 2, 5)) ==  rNa$year[2] )]],
#                  datM[[which(as.numeric(substr(names(datM), 2, 5)) ==  rNa$year[3])]],
#                  datM[[which(as.numeric(substr(names(datM), 2, 5)) ==  rNa$year[4])]],
#                  datM[[which(as.numeric(substr(names(datM), 2, 5)) ==  rNa$year[5])]])             
datN <- datM[[which(as.numeric(substr(names(datM), 2, 5)) %in% uq_yrs) ]]

ssn <- substring(names(datN), 7)
datNm <- stackApply(datN, ssn, mean)
names(datNm) <- unique(ssn)

SeasonsN <- datNm - datMs 

#### Plotting the composites

library(rgdal)
library(rgeos)
coast_shapefile <- crop(readOGR("../../Climate_analysis/GISData/ne_10m_coastline.shp"), ext)

#Create color ramps for mapping and number of colors to use
library(colorRamps)
col5 <- colorRampPalette(c('#08519c', 'gray96', "#fee0d2", "firebrick3"))

#Load the lattice packages to display the maps
library(rasterVis)
library(gridExtra)


#### 30yr Mean - Narrow years
levelplot(SeasonsN, layout=c(2,2), col.regions = col5, pretty=TRUE, main="Narrow Composite Mean SLP Diff 1979 - 2011",
          colorkey=list(space="bottom"),
          par.settings = list(layout.heights=list(xlab.key.padding=1),
                              strip.background=list(col="lightgrey")
          ), par.strip.text = list(font="bold")) + 
  layer(sp.lines(coast_shapefile))


#### 30yr Mean - Wide years
levelplot(SeasonsW, layout=c(2,2), col.regions = col5, pretty=TRUE, main="Wide Composite Mean SLP Diff 1979 - 2011",
          colorkey=list(space="bottom"),
          par.settings = list(layout.heights=list(xlab.key.padding=1),
                              strip.background=list(col="lightgrey")
          ), par.strip.text = list(font="bold")) + 
  layer(sp.lines(coast_shapefile))

#### Narrow - Wide
SeasonsNW <- datNm - datWm

levelplot(SeasonsNW, layout=c(2,2), col.regions = col5, pretty=TRUE, main="Narrow-Wide Composite Mean SLP Diff 1979 - 2011",
          colorkey=list(space="bottom"),
          par.settings = list(layout.heights=list(xlab.key.padding=1),
                              strip.background=list(col="lightgrey")
          ), par.strip.text = list(font="bold")) + 
  layer(sp.lines(coast_shapefile))

#######Can we delete stuff below here? Or move to separate script?#############
#Load the lattice packages to display the maps
library(rasterVis)
library(gridExtra)


#### 30yr Mean - Narrow years
levelplot(SeasonsN, layout=c(2,2), col.regions = col5, pretty=TRUE, main="Narrow Composite Mean SLP Diff 1979 - 2011",
          colorkey=list(space="bottom"),
          par.settings = list(layout.heights=list(xlab.key.padding=1),
                              strip.background=list(col="lightgrey")
          ), par.strip.text = list(font="bold")) + 
  layer(sp.lines(coast_shapefile))


#### 30yr Mean - Wide years
levelplot(SeasonsW, layout=c(2,2), col.regions = col5, pretty=TRUE, main="Wide Composite Mean SLP Diff 1979 - 2011",
          colorkey=list(space="bottom"),
          par.settings = list(layout.heights=list(xlab.key.padding=1),
                              strip.background=list(col="lightgrey")
          ), par.strip.text = list(font="bold")) + 
  layer(sp.lines(coast_shapefile))

#### Narrow - Wide
SeasonsNW <- datNm - datWm

levelplot(SeasonsNW, layout=c(2,2), col.regions = col5, pretty=TRUE, main="Narrow-Wide Composite Mean SLP Diff 1979 - 2011",
          colorkey=list(space="bottom"),
          par.settings = list(layout.heights=list(xlab.key.padding=1),
                              strip.background=list(col="lightgrey")
          ), par.strip.text = list(font="bold")) + 
  layer(sp.lines(coast_shapefile))


#### Looking at climate correlation - not relevant to composite ####

# If data is not continuous you must
# replace all NA with a number outside of the bounds of data in order to do correlations.
# For temperature, 0 will not work since values do go to 0 or below. 
# Done at this stage because data can be too large.
# Also subsets seasons at same time
#Creates SONw etc. for mean of wide years
for (i in unique(substring(yr_season, 6))){
  d <- values(subset(datW, grep(i, names(datW), value=T)))
  d[is.na(d[])] <- -9999
  assign(paste0(i, "w"), d)
  #Replace with first differencing using diff()
  lm_x <- seq(1:dim(get(i))[2])
  r <- t(resid(lm(t(get(i)) ~ lm_x))) 
  assign(paste0(i, "w"), r)
  rm(r, d, lm_x)
}

#Plot climate field composite for wide years

#create rasters to ingest the seasonal composites
CorT <- setExtent(raster(nrow = nrow(datW), ncol = ncol(datW)),ext)
Cor <- setExtent(raster(nrow = nrow(datW), ncol = ncol(datW)),ext)

##### DJF
for(i in 1:dim(DJF.w)[1]){
  Cor[i] <- cor(x=DJF.w[i,], y = rWi$mr_kbp, method = 'pearson') 
  CorT[i] <- cor.test(x=DJF.w[i,], y = rWi$mr_kbp, method = 'pearson')$p.value 
}
CorT[CorT > 0.05] <- NA
## Using the altered data frame, create the cropped field
DJF.wC <- mask(Cor, CorT)


#### MAM
for(i in 1:dim(MAM.w)[1]){
  Cor[i] <- cor(x=MAM.w[i,], y = rWi$mr_kbp, method = 'pearson') 
  CorT[i] <- cor.test(x=MAM.w[i,], y = rWi$mr_kbp, method = 'pearson')$p.value 
}
CorT[CorT > 0.05] <- NA
## Using the altered data frame, create the cropped field
MAM.wC <- mask(Cor, CorT)

#### JJA
for(i in 1:dim(JJA.w)[1]){
  Cor[i] <- cor(x=JJA.w[i,], y = rWi$mr_kbp, method = 'pearson') 
  CorT[i] <- cor.test(x=JJA.w[i,], y = rWi$mr_kbp, method = 'pearson')$p.value 
}
CorT[CorT > 0.05] <- NA
## Using the altered data frame, create the cropped field
JJA.wC <- mask(Cor, CorT)

#### SON
for(i in 1:dim(SON.w)[1]){
  Cor[i] <- cor(x=SON.w[i,], y = rWi$mr_kbp, method = 'pearson') 
  CorT[i] <- cor.test(x=SON.w[i,], y = rWi$mr_kbp, method = 'pearson')$p.value 
}
CorT[CorT > 0.05] <- NA
## Using the altered data frame, create the cropped field
SON.wC <- mask(Cor, CorT)
