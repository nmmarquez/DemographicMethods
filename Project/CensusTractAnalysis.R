rm(list=ls())

library(sp)
library(rgr)
library(spdep)
library(rgdal)
library(raster)
library(dplyr)
library(ar.matrix)
library(leaflet)
library(DemographicSimulation)
library(surveillance)
setwd("~/Documents/Classes/DemographicMethods/Project/Tract2010")

file_zip <- tempfile()
extr_dir <- tempdir()
url_ <- "http://www2.census.gov/geo/tiger/GENZ2015/shp/cb_2015_us_county_5m.zip"
download.file(url_, file_zip)
unzip(file_zip, exdir=extr_dir)

# pull the shape file into memory
USCounty <- readOGR(extr_dir)
USCounty <- USCounty[as.numeric(as.character(USCounty$STATEFP))<=56,]
unlink(file_zip)
unlink(extr_dir)
USCounty@data$GEOID <- paste0(USCounty$STATEFP, USCounty$COUNTYFP) 


DF <- readOGR(".")
DF@data$StateID <- substring(DF@data$GEOID10, 1, 2) # State code
DF@data$GEOID <- substring(DF@data$GEOID10, 1, 5) # County Code
DFGeoStates <- DF[DF$StateID != "72",] # Remove Puerto Rico
DFGeoStates <- DFGeoStates[DFGeoStates$DP0010001 != 0,] # Remove 0 Pop Locations
length(table(DF$GEOID))# check number of counties to see if we match expected

ruralCutoff <- 80000

# Whats the population split living urban rural in the US? about 20%
DFGeoStates@data %>% group_by(GEOID) %>% 
    summarize(pop=sum(DP0010001)) %>%
    mutate(rural=pop<ruralCutoff) %>% 
    group_by(rural) %>%
    summarize(tpop=sum(pop)) %>%
    mutate(tpop=tpop/sum(tpop))

# Cool now what about the county distribution whats that look like
DFGeoStates@data %>% group_by(GEOID) %>% 
    summarize(pop=sum(DP0010001)) %>%
    mutate(rural=pop<ruralCutoff) %>% 
    group_by(rural) %>%
    summarize(prural=n()) %>%
    mutate(prural=prural/sum(prural))

# Almost the exact opposite 80 percent of counties are rural

# lets combine the data
USCounty@data <- DFGeoStates@data %>% group_by(GEOID) %>% 
    summarize(pop=sum(DP0010001)) %>%
    mutate(Rural=pop<ruralCutoff) %>%
    select(GEOID, Rural) %>%
    right_join(USCounty@data, by="GEOID") %>%
    as.data.frame

# simulate some values
set.seed(123)
USCounty.graph <- poly2adjmat(USCounty)
USCounty$zeta <- Q.lCAR(USCounty.graph, .6, .93) %>% sim.AR(1, .) %>% c
USCounty$prob <- expit(USCounty$zeta*2) * .2 + ((!USCounty$Rural) * .6)
USCounty@data[is.na(USCounty$prob),"prob"] <- .5
USCounty$prob %>% density %>% plot

# Join this to the census tract dataset
DFGeoStates@data <- USCounty@data %>% select(GEOID, prob) %>% 
    right_join(DFGeoStates@data, by="GEOID")
DFGeoStates@data[is.na(DFGeoStates$prob),"prob"] <- .5
DFGeoStates@data$Population <- DFGeoStates@data$DP0010001
DFGeoStates@data <- DFGeoStates@data %>% 
    select(GEOID, GEOID10, prob, Population, StateID)
DFGeoStates@data$PopIndexEnd <- cumsum(DFGeoStates@data$Population)
DFGeoStates@data$PopIndexStart <- 
    c(1, DFGeoStates@data$PopIndexEnd[1:(nrow(DFGeoStates) -1)] + 1)
DFGeoStates@data$LxType <- rbinom(nrow(DFGeoStates@data), 1, DFGeoStates$prob)
DFGeoStates$prob %>% density %>% plot

Npop <- sum(DFGeoStates$Population)
JPFitFunc <- paramFitFunc(start_params = c(0, 0, 1.1, 0.57, 80),
             location_ = "Japan", year_ = 2016, sex_ = "Both",
             max_age = 140, returnParams = FALSE)

MXFitFunc <- paramFitFunc(start_params = c(0, 0, 1.1, 0.57, 80),
                          location_ = "Mexico", year_ = 2000, sex_ = "Both",
                          max_age = 140, returnParams = FALSE)

system.time(MxSims <- unlist(lapply(1:nrow(DFGeoStates@data), function(i){
    if(DFGeoStates$LxType[i] == 1){
        z <- JPFitFunc(DFGeoStates$Population[i])
        if(length(z) != DFGeoStates$Population[i]){
            print("Jp problme")
            print(i)
            print("Diff length")
            print(length(z) - DFGeoStates$Population[i])
        }
        return(z)
        #return(rnorm(DFGeoStates$Population[i]))
    }
    else{
        z <- MXFitFunc(DFGeoStates$Population[i])
        if(length(z) != DFGeoStates$Population[i]){
            print("Mx problme")
            print(i)
            print("Diff length")
            print(length(z) - DFGeoStates$Population[i])
        }
        return(z)
        #return(rnorm(DFGeoStates$Population[i]))
    }
        })))

DFGeoStates$CensusID <- substring(DFGeoStates$GEOID10, 6, 12)

# sample randomly from a polygon

MxSims[subset(DFGeoStates@data, CensusID == "074111")[,"PopIndexStart"]]
MxSims[subset(DFGeoStates@data, CensusID == "000502")[,"PopIndexStart"]][1]
MxSims[subset(DFGeoStates@data, CensusID == "011201")[,"PopIndexStart"]][1]

USCounty$CountyAvgLx <- sapply(USCounty$GEOID, function(g){
    DFsub <- subset(DFGeoStates@data, GEOID == g)
    if(nrow(DFsub) == 0){
        return(NA)
    }
    mean(unlist(lapply(1:nrow(DFsub), function(i){
        six <- DFsub$PopIndexStart[i]
        eix <- DFsub$PopIndexEnd[i]
        MxSims[six:eix]
    })))
})

USCounty@data[is.na(USCounty$CountyAvgLx),"CountyAvgLx"] <-
    mean(USCounty$CountyAvgLx, na.rm=TRUE)

pal <- colorNumeric(palette="YlGnBu", domain=USCounty@data$CountyAvgLx)

popup <- paste0("Loc Name: ", USCounty@data$NAME,
                "<br>Lx: ", USCounty@data$CountyAvgLx,
                "<br>Rural: ", USCounty@data$Rural)
# see map
map1 <-leaflet() %>%
    addProviderTiles("CartoDB.Positron") %>%
    addPolygons(data=USCounty, fillColor=~pal(CountyAvgLx), color="#b2aeae",
                fillOpacity=0.7, weight=0.3, smoothFactor=0.2,
                popup=popup) %>%
    addLegend("bottomright", pal=pal, values=USCounty$CountyAvgLx, 
              title="Avg Life<br>Expectancy", opacity=1)
map1

DFGeoStates$CensusAvgLx <- sapply(1:nrow(DFGeoStates), function(i){
        six <- DFGeoStates$PopIndexStart[i]
        eix <- DFGeoStates$PopIndexEnd[i]
        mean(MxSims[six:eix])})

CAGeoStates <- DFGeoStates[DFGeoStates$StateID == "06",]

popup2 <- paste0("Loc Name: ", CAGeoStates@data$GEOID10,
                "<br>Lx: ", CAGeoStates@data$CensusAvgLx,
                "<br>Rural: ", CAGeoStates@data$prob)
pal2 <- colorNumeric(palette="YlGnBu", domain=CAGeoStates@data$CensusAvgLx)

# see map
map2 <-leaflet() %>%
    addProviderTiles("CartoDB.Positron") %>%
    addPolygons(data=CAGeoStates, fillColor=~pal2(CensusAvgLx), color="#b2aeae",
                fillOpacity=0.7, weight=0.3, smoothFactor=0.2,
                popup=popup2) %>%
    addLegend("bottomright", pal=pal2, values=CAGeoStates$CensusAvgLx, 
              title="Avg Life<br>Expectancy", opacity=1)
map2


# map southern california points
CountyZoom <- c("06059")#, "06037")

SocalGeoStates <- CAGeoStates[CAGeoStates$GEOID %in% CountyZoom,]

spsample(DFGeoStates[1,], 10, "random")

DFpoints <- lapply(1:nrow(SocalGeoStates@data), function(i){
    DFp <- spsample(SocalGeoStates[i,], SocalGeoStates$Population[i], "random")
    six <- SocalGeoStates$PopIndexStart[i]
    eix <- SocalGeoStates$PopIndexEnd[i]
    SpatialPointsDataFrame(
        DFp,
        data.frame(Type=rep(SocalGeoStates$LxType[i],
                            SocalGeoStates$Population[i])) %>%
            mutate(death=MxSims[six:eix]))
})

DFpointsSpatial <- SpatialPoints(
    do.call(rbind, lapply(DFpoints, function(x) x@coords)),
    proj4string=DFpoints[[1]]@proj4string)

DFpointsSpatial <- SpatialPointsDataFrame(
    DFpointsSpatial, 
    do.call(rbind, lapply(DFpoints, function(x) x@data))
)

DFsub <- DFpointsSpatial[sample(1:nrow(DFpointsSpatial), 10^4),]

map3 <- leaflet() %>%
    addProviderTiles("CartoDB.Positron") %>%
    addCircleMarkers(
        data=DFsub,
        radius=.5, color=c("red","blue")[DFsub$Type+1],
        popup=as.character(round(DFsub$death, 4)))

MapList <- list(map1, map2, map3)

saveRDS(MapList, "~/Documents/Classes/DemographicMethods/Project/Maps.Rds")

rm(list=ls())
MapList2 <- readRDS("~/Documents/Classes/DemographicMethods/Project/Maps.Rds")
MapList2[[1]]
MapList2[[2]]
MapList2[[3]]
