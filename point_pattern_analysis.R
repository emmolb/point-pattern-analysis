# Point pattern analysis
# Example 1 : Drumlins data
# source data : spatialpoints shapefile
# converting data
setwd(xxxxx)
library(maptools)
library(sp)
library(spatstat)
drumlins <- readShapePoints("drumlins.shp")
class(drumlins)
plot(drumlins)
drumlins_SP <- as(drumlins,"SpatialPoints")
drumlins_ppp <- as(drumlins_SP,"ppp")
class(drumlins_ppp)
plot(drumlins_ppp)
plot(density(drumlins_ppp))
plot(drumlins_ppp,add=TRUE)
# specification of observation window
# bounding box
bb <- bounding.box(drumlins_ppp)
bb
W <- as(bb,"owin")
class(W)
plot(W)
# alternatives for bounding box : convexhull, ripras
ch <- convexhull.xy(drumlins_ppp)
rr <- ripras(drumlins_ppp)
# add observation window to point pattern
drumlins_rr <- ppp(drumlins_ppp$x,drumlins_ppp$y,window=rr)
# density estimates
k025 <- (density(drumlins_rr,0.25))
plot(k025)
plot(drumlins_rr,add=TRUE,col="darkgray")
SG <- as(k025,"SpatialGridDataFrame")
k050 <- (density(drumlins_rr,0.50))
plot(k050)
plot(drumlins_rr,add=TRUE,col="darkgray")

# Example 2
# source data : spatialpolygon shapefile for window and csv-file for points
setwd(xxxxx)
library(maptools)
library(sp)
library(spatstat)
S <- readShapePoly("city_limits_km.shp")
plot(S)
class(S)
SP <- as(S,"SpatialPolygons")
W <- as(SP,"owin")
plot(W)
class(W)
class(SP)
xy <- read.table("StLouisCrime1982.txt",header=T,sep="\t")
str(xy)
head(xy)
tail(xy)
attach(xy)
pp <- ppp(X,Y,window=W,marks=CRIME)
summary(pp)
plot(pp)
plot(split(pp))
plot(density(pp))         
gun <- pp[CRIME=="gun_homicide"]
hit <- pp[CRIME=="hit_and_run"]
rob <- pp[CRIME=="street_robbery"]
plot(density(gun))
contour(density(gun),col="red",add=T)
plot(gun,col="gray",add=T)
plot(density(hit))
contour(density(hit),col="red",add=T)
plot(hit,col="gray",add=T)
plot(density(rob))
contour(density(rob),col="red",add=T)
plot(rob,col="gray",add=T)
# R provides a function that will suggest on optimal bandwidth to use
r <- bw.diggle(gun)
r
plot(density(gun,0.533))
contour(density(gun),col="red",add=T)
plot(gun,col="gray",add=T)

# Example 3 (London transport points)
# source data : R workspace (2 shapefiles : SpatialPolygonsDataFrame & SpatialPointsDataFrame)
# the aim of this exercise is to demonstrate how clusters can be identified in spatial points 
# these clusters are then used as the basis of aggregation to reduce the total number of points
library(maptools)
library(sp)
library(spatstat)
setwd(xxxxx)
load("lnd.Rdata")
class(lnd)
plot(lnd)
load("stations.Rdata")
class(stations)
plot(stations)
stations_SP <- as(stations,"SpatialPoints")
stations_ppp <- as(stations_SP,"ppp")
class(stations_ppp)
plot(stations_ppp)
plot(density(stations_ppp))
plot(stations_ppp,add=T)
contour(density(stations_ppp))
dens <- density(stations_ppp)
class(dens)
plot(dens)
r <- bw.diggle(stations_ppp)
r
dens <- density(stations_ppp,1493.611)
class(dens)
plot(dens)
plot(density(stations_ppp,1493.611))
contour(density(stations_ppp,1493.611),add=T)
# the density plot illustrates that there are some areas of relatively high density 
# this concentration is also illustrated by the contour lines
# having identified the "islands" of hot dot density, the next stage is to save them into polygons
# convert to spatial grid
Dsg <- as(dens,"SpatialGridDataFrame")
class(Dsg)
# convert again to an image
Dim <- as.image.SpatialGridDataFrame(Dsg)
Dcl <- contourLines(Dim,nlevels=9)
# convert to SpatialLinesDataFrame
SLDF <- ContourLines2SLDF(Dcl,CRS(proj4string(lnd)))
plot(SLDF,col = terrain.colors(8))
library(rgeos)
Polyclust <- gPolygonize(SLDF[5, ])
gas <- gArea(Polyclust,byid=T)/10000
Polyclust <- SpatialPolygonsDataFrame(Polyclust,data=data.frame(gas),match.ID=F)
plot(Polyclust)
cag <- aggregate(stations,by=Polyclust,FUN=length)
plot(dens,main=" ")
plot(lnd,border="grey",lwd=2,add=T)
plot(SLDF,col=terrain.colors(8),add=T)
plot(cag,col="red",border="white",add=T)
graphics::text(coordinates(cag) + 1000, labels = cag$CODE)
# points inside and outside the clusters
sIn <- stations [cag, ]
plot(sIn)
sOut <- stations[!row.names(stations) %in% row.names(sIn), ]
plot(sOut,col="blue",add=T)
plot(cag,border="red",lwd=3,add=T)
nrow(sIn)/nrow(stations)

# Example 4
# source data : spatialpolygon shapefile for window and csv-file for points
setwd(xxxxx)
# read shapefile Belgium
library(maptools)
Belgie <- readShapePoly("Bel_adm0.shp")
plot(Belgie)
# pointdevue <- readShapePoints("pointdevue.shp")
# class(pointdevue)
# plot(pointdevue)
# read pointdata
pos <- read.csv("poslatlon.csv",header=TRUE,sep=";")
str(pos)
head(pos)
points(pos$lon,pos$lat,pch=".",col="blue")
# define observation window based on shapefile
library(spatstat)
library(sp)
SP <- as(Belgie,"SpatialPolygons")
W <- as(SP,"owin")
plot(W)
# define point pattern with marks as a continuous variable
PP <- ppp(pos$lon,pos$lat,window=W,marks=pos$sales)
plot(PP)
class(PP)
summary(PP)
r <- bw.diggle(PP)
r
plot(density(PP,sigma=0.005878732))
plot(density(PP))
plot(Smooth(PP))
PP2 <- ppp(pos$lon,pos$lat,window=W,marks=pos$unsolds)
plot(PP2)
plot(Smooth(PP2))

# Example 5
# source data : csv-file (with north/east coordinates)
setwd(xxxxx)
data <- scan("2005.csv",skip=1,what=list(tear=0,north=0,east=0),sep=",",quiet=T)
str(data)
plot(data$north,data$east)
library(MASS)
# create a density plot using kde2d-function
# kernel is a list containing 3 items : x and y representing the grid and z a 500 * 500 matrix
# containing the kernel density for each cell of the grid
kernel <- kde2d(data$north,data$east,n=500)
image(kernel,col=c(0,rainbow(50)))
# to add a colour legend to the z-plot, use package "fields"
library(fields)
image.plot(kernel,col=c(0,rainbow(50)))
# applying the same logic to crime data Chicago area, 2013
setwd("c:/R/Rdata/chicago")
chicago <- read.csv("chicagocrimes_2013.csv",header=T,sep=",")
data <- subset(chicago,Primary.Type=="ARSON" | Primary.Type=="BURGLARY" | Primary.Type=="ASSAULT" | Primary.Type=="ROBBERY" | 
                 Primary.Type=="NARCOTICS" | Primary.Type=="HOMICIDE" | Primary.Type=="SEX OFFENSE",
               select=c(Longitude,Latitude,Primary.Type))
data <- na.omit(data)
library(MASS)
kernel <- kde2d(data$Longitude,data$Latitude,n=500)
image(kernel,col=c(0,rainbow(50)))
str(kernel)
class(kernel)
head(kernel)
data <- subset(chicago,Primary.Type=="ARSON" | Primary.Type=="BURGLARY" | Primary.Type=="ASSAULT" | Primary.Type=="ROBBERY" | 
                 Primary.Type=="NARCOTICS" | Primary.Type=="HOMICIDE" | Primary.Type=="SEX OFFENSE",
               select=c(X.Coordinate,Y.Coordinate,Primary.Type))
data <- na.omit(data)
str(data)
plot(data$X.Coordinate,data$Y.Coordinate)
library(MASS)
kernel <- kde2d(data$X.Coordinate,data$Y.Coordinate,n=500)
image(kernel,col=c(0,rainbow(50)))
# to add a colour legend to the z-plot, use package "fields"
library(fields)
image <- image.plot(kernel,col=c(0,rainbow(50)))
# add shapefile boundary to the plots
# read shapefile, convert kernel density estimates to raster and finally clip raster
# read shapefile Chicago
library(maptools)
ch <- readShapePoly("City_Boundary.shp")
plot(ch)
library(raster)
r <- raster(kernel)
plot(r)
projection(r) <- CRS("+init=epsg:4326")
writeRaster(r,"chicagocrime.tiff","GTiff",overwrite=TRUE)
r2 <- crop(r,extent(ch))
r3 <- mask(r2,ch)
plot(r3)
plot(ch, add=TRUE, lwd=2)
# clipping the raster file with a shapefile containing multiple polygons
ch2 <- readShapePoly("wards.shp")
plot(ch2)
plot(r)
rr <- mask(r,ch2)
plot(rr)
plot(ch2,add=TRUE)

# Example 6 : clipping a raster (1)
library(maptools)  ## For wrld_simpl
library(raster)
## Example SpatialPolygonsDataFrame
data(wrld_simpl)
SPDF <- subset(wrld_simpl, NAME=="Brazil")
plot(SPDF)
## Example RasterLayer
r <- raster(nrow=1e3, ncol=1e3)
r[] <- 1:length(r)
plot(r)
## crop and mask
r2 <- crop(r, extent(SPDF))
r3 <- mask(r2, SPDF)
## Check that it worked
plot(r2)
plot(r3)
plot(SPDF, add=TRUE, lwd=2)

# Example 7 : clipping a raster (2)
setwd(xxxxx)
library(rgdal)
library(raster)
# shapefile US = "USA_adm1.shp"
# get this shapefile from www.gadm.org
us <- getData("GADM",country="USA",level=1)
summary(us)
plot(us)
str(us@data)
# subset US shapefile to desired states
selectstates <- c("Maine","Vermont","Massachusetts","New Hapshire","Connecticut","Rhode Island",
                  "New York","Pennsylvania","New Jersey","Maryland","Delaware","Virginia","West Virginia")
us.sub <- us[as.character(us@data$NAME_1) %in% selectstates, ]
plot(us.sub)
# create a random raster over the space
r <- raster(xmn=-85,xmx=-65,ymn=36,ymx=48,nrow=100,ncol=100)
r[] <- runif(100*100)
plot(r)
plot(us.sub,lwd=2,add=T)
# use the mask function
rr <- mask(r,us.sub)
plot(rr)
plot(us.sub,lwd=2,add=T)
# use precipitation rasterfile
precipitation <- raster("prec1.bil")
plot(precipitation)    
# crop precipitation data by extent of states subset
precipitation.sub <- crop(precipitation,extent(us.sub))
plot(precipitation.sub)
# as a final step we need to identify those pixels of the precipitation rater file
# that lie within the borders of the states subset file
precipitation.sub <- mask(precipitation.sub,us.sub)
plot(precipitation.sub)
plot(us.sub,lwd=2,add=T)

# Example 8 : smooth.ppp
df = data.frame(depth = c(-5, -10, -15, -20, -25, -30, -35, 
                          -5, -10, -15, -20, -25, -30, -35, 
                          -5, -10, -15, -20, -25, -30, -35),
                date = as.POSIXlt(c("2014-02-02", "2014-02-02", "2014-02-02", "2014-02-02", "2014-02-02",
                                    "2014-02-02", "2014-02-02", "2014-04-14", "2014-04-14", "2014-04-14", 
                                    "2014-04-14", "2014-04-14", "2014-04-14", "2014-04-14", "2014-07-21",
                                    "2014-07-21","2014-07-21","2014-07-21","2014-07-21","2014-07-21","2014-07-21")),
                temperature = c(11.9071,11.9657,11.9751,11.9813,
                                11.9972,12.0255,12.1483,17.0442,
                                14.3784,14.2104,14.2206,14.1834,
                                14.1979,14.2189,18.4762,16.3302,
                                15.1438,14.0497,13.346,13.0996,13.0504))
df$month <- as.numeric(format(df$date, "%m"))
str(df)
head(df)
library(spatstat)
C <- ppp(df$month, df$depth, c(2,7), c(-35, -5))
marks(C) <- df[,"temperature"]
C.smooth <- Smooth.ppp(C, sigma = 2)
plot(Smooth(C), col = topo.colors(128), 
     main="Temperature", ylab="Depth (m)", 
     ylim=TRUE,xlab="Month")
contour(C.smooth, add = TRUE)

# Example 9
# kernel density estimation map of point pattern with package spatstat
# data : polygon shapefile & csv-pointfile
setwd(xxxxx)
library(maptools)
library(sp)
library(spatstat)
city_shp <- readShapePoly("lond_city.shp")
city_win <- as(as(city_shp,"SpatialPolygons"),"owin")
plot(city_win)                
bikedata <- read.csv("London_cycle_hire_locs.csv",header=TRUE,sep=",")
str(bikedata)
head(bikedata)
summary(bikedata$X)
summary(bikedata$Y)
bikes <- ppp(bikedata$X,bikedata$Y,marks=bikedata$Capacity,window=city_win)
bikes <- ppp(bikedata$X,bikedata$Y,c(524400,534700),c(177900,183600),marks=bikedata$Capacity)
points(bikedata$X,bikedata$Y)
plot(bikes)
plot(density(bikes))
r <- bw.diggle(bikes)
r
plot(density(bikes,102.9335))
k125 <- density.ppp(bikes,sigma=125,weights=bikes$Capacity)
plot(k125)
plot(k125,col=heat.colors(256))
plot(Smooth.ppp(bikes))
plot(bikes,add=T)
plot(city_win,lwd=2,add=T)

# Example 10
# the following example will outline how to create a surface using kernel density estimation (KDE)
# and then clip the surface so that it is constrained within the limits of a polygon
# data : polygon shapefile & csv-pointfile
setwd(xxxxx)
cycle <- read.csv("London_cycle_hire_locs.csv",header=TRUE,sep=",")
str(cycle)
head(cycle)
plot(cycle$X,cycle$Y)
library(maptools)
lon_city <- readShapePoly("lond_city.shp")
summary(lon_city)
str(lon_city@data)
plot(cycle$X,cycle$Y)
plot(lon_city,add=TRUE,lwd=2)
# create a kernel density surface based on the locations of the points
# we use the sm.density function in the sm package
library(sm)
cycle_dens <- sm.density(data.frame(cycle$X,cycle$Y),weights=cycle$Capacity,display="image",ngrid=100)
# add the points and the polygon
points(cycle$X,cycle$Y)
plot(lon_city,add=T,lwd=2)
str(cycle_dens)
## We can convert the cycle_dens output into a spatial grid for further spatial analysis. 
temp=SpatialPoints(expand.grid(x=cycle_dens$eval.points[,1], y=cycle_dens$eval.points[,2]))
temp = SpatialPixelsDataFrame(temp, data.frame(kde = array(cycle_dens$estimate, 
                                                           length(cycle_dens$estimate))))
sel=over(temp,lon_city)
str(sel)
temp2 <- temp[!is.na(sel[, 1]), ]
temp2 <- temp2[lon_city, ]
image(temp2)
plot(lon_city, add=T, lwd=2)
points(cycle$X, cycle$Y)
title("Density of London Cycle Hire Bikes in the City of London")

# Example 11
# same excercise as example 10
# data : polygon shapefile & csv-pointfile
setwd(xxxxx)
outlets <- read.csv("outletlatlon.csv",header=T,sep=";")
str(outlets)
plot(outlets$longitude,outlets$latitude)
library(maptools)
belgie <- readShapePoly("bel_adm0.shp")
summary(belgie)
plot(outlets$longitude,outlets$latitude)
plot(belgie,add=TRUE,lwd=2)
# create a kernel density surface based on the locations of the points
# we use the sm.density function in the sm package
library(sm)
outlets_dens <- sm.density(data.frame(outlets$longitude,outlets$latitude),display="image",ngrid=100)
# add the points and the polygon
points(outlets$longitude,outlets$latitude)
plot(belgie,add=T,lwd=2)
## We can convert the cycle_dens output into a spatial grid for further spatial analysis. 
temp=SpatialPoints(expand.grid(x=outlets_dens$eval.points[,1], y=outlets_dens$eval.points[,2]))
temp = SpatialPixelsDataFrame(temp, data.frame(kde = array(outlets_dens$estimate, 
                                                           length(outlets_dens$estimate))))
sel=over(temp,belgie)
str(sel)
temp2 <- temp[!is.na(sel[, 1]), ]
temp2 <- temp2[belgie, ]
image(temp2)
plot(belgie, add=T, lwd=2)
points(outlets$longitude, outlets$latitude)
title("Density of retail outlets in Belgium")

# Example 13 : point analysis Wikileaks Irad Data

setwd(xxxxx)
deaths <- read.csv("Deaths only.csv",header=TRUE,sep=",")
str(deaths)
levels(deaths$Type)
deaths$Type <- sub(deaths$Type,pattern="CRIMINAL EVENT",replacement="Criminal Event")
deaths$Type <- sub(deaths$Type,pattern="criminal event",replacement="Criminal Event")
deaths$Type <- sub(deaths$Type,pattern="EXPLOSIVE HAZARD",replacement="Explosive Hazard")
deaths$Type <- factor(deaths$Type)
levels(deaths$Type)
# create a map with GGMAP
library(ggmap)
library(sp)
library(rgdal)
map <- get_map(location="Iraq",maptype="hybrid",zoom=6)
ggmap(map,extent="panel")
# put the data on the map
ggmap(map,extent="panel") + geom_point(data=deaths,aes(x=Longitude,y=Latitude),alpha=0.1,color="red")
ggmap(map,extent="panel") + geom_point(data=deaths,aes(x=Longitude,y=Latitude,size=Total.deaths),alpha=0.1,color="red")
ggmap(map,extent="panel") + geom_point(data=deaths[deaths$Enemy.detained > 0, ],aes(x=Longitude,y=Latitude,size=Enemy.detained),alpha=0.1,color="green")
# create a spatialpointsdataframe from the coordinates in deaths
class(deaths)
deaths.spp <- SpatialPointsDataFrame(coords=deaths[,18:19],proj4string=CRS("+proj=longlat +ellps=WGS84 +datum=WGS84 +no_defs"),data=deaths)
class(deaths.spp)
str(deaths.spp)
proj4string(deaths.spp)
# reproject data into meters
deaths.spp.rp <- spTransform(deaths.spp,CRS("+init=epsg:3839"))
proj4string(deaths.spp.rp)
str(deaths.spp.rp)
deaths.spp.rp$Easting <- coordinates(deaths.spp.rp)[,1]
deaths.spp.rp$Northing <- coordinates(deaths.spp.rp)[,2]
str(deaths.spp.rp)
plot(deaths.spp.rp@data$Longitude,deaths.spp.rp@data$Latitude)
plot(deaths.spp.rp@data$Easting,deaths.spp.rp@data$Northing)
# focus on the region of Baghdad
deaths.df.rp.bag <- deaths.spp.rp@data[deaths.spp.rp@data$Region=="MND-BAGHDAD",]
# note that this is NOT a spatial data frame (the selection is done on ..@data !)
class(deaths.df.rp.bag)
str(deaths.df.rp.bag)
plot(deaths.df.rp.bag$Easting,deaths.df.rp.bag$Northing)
# select data within a given bounds
deaths.df.rp.bag <- deaths.df.rp.bag[(deaths.df.rp.bag$Northing < 4950000 & deaths.df.rp.bag$Northing > 4940000) &
                                       (deaths.df.rp.bag$Easting > 10000000 & deaths.df.rp.bag$Easting < 10010000),]
str(deaths.df.rp.bag)
plot(deaths.df.rp.bag$Easting,deaths.df.rp.bag$Northing)
deaths.df.rp.bag <- na.omit(deaths.df.rp.bag)
library(spatstat)
deaths.chull <- convexhull.xy(x=deaths.df.rp.bag$Easting,y=deaths.df.rp.bag$Northing)
plot(deaths.chull)
deaths.ppp <- ppp(x=deaths.df.rp.bag$Easting,y=deaths.df.rp.bag$Northing,window=deaths.chull,
                  marks=deaths.df.rp.bag[,c(3,8:17)])
summary(deaths.ppp)
plot(deaths.ppp)
plot(split(deaths.ppp))
dens <- density(deaths.ppp)
plot(dens)
contour(density(deaths.ppp))
plot(density(split(deaths.ppp)))
contour(density(split(deaths.ppp)))
tdeaths.ppp <- ppp(x=deaths.df.rp.bag$Easting,y=deaths.df.rp.bag$Northing,window=deaths.chull,
                   marks=deaths.df.rp.bag[,c(17)])
plot(tdeaths.ppp,maxsize=3000)
plot(density(tdeaths.ppp))
plot(Smooth(tdeaths.ppp))
points <- ppp(x=deaths.df.rp.bag$Easting,y=deaths.df.rp.bag$Northing,window=deaths.chull,
              marks=deaths.df.rp.bag[,c(17)])
plot(points)
plot(density(points))
a <- density(points,weights=deaths.df.rp.bag$Enemy.detained)
plot(a)

# Example 14 : stat density plots of Baltimore crime data
# http://www.obscureanalytics.com/2012/12/07/visualizing-baltimore-with-r-and-ggplot2-crime-data
# http://www.obscureanalytics.com/2012/12/10/visualizing-baltimore-2-vacant-property-and-some-more-crime
library(ggplot2)
library(foreign)
library(stringr)
library(lubridate)
library(plyr)
library(xtable)
library(scales)
library(RColorBrewer)
library(ggmap)
library(maptools)
library(sp)
library(rgdal)
library(spatstat)
setwd(xxxxx)
city_shp <- readOGR(".","baltcity_line")
class(city_shp)
summary(city_shp)
proj4string(city_shp)
# if we want to convert this shapefile to lat/lon :
# city_shp <- spTransform(city_shp,CRS("+proj=longlat +datum=WGS84"))
# we store the orginal projection of the shapefile :
origProj <- city_shp@proj4string
origProj
# ggplot only takes data frames, so we convert the shapefile to a data frame
str(city_shp)
city_shp_df <- fortify(city_shp,region="LABEL")
str(city_shp_df)
# we use this data frame as the first layer
bound_plot <- ggplot(data=city_shp_df,aes(x=long,y=lat,group=group)) +
  geom_polygon(color="gray",fill="lightblue") +
  coord_equal() + theme_nothing()
bound_plot
# neighborhoods
nbhds_shp <- readOGR(".","Neighborhoods_2010")
class(nbhds_shp)
summary(nbhds_shp)
str(nbhds_shp@data)
# we store the orginal projection of the shapefile :
origProj <- nbhds_shp@proj4string
origProj
# convert to data frame with fortify
nbhds_pl_df <- fortify(nbhds_shp,region="Name")
plot_nbhds_pl_df <- ggplot(data=nbhds_pl_df,aes(x=long,y=lat,group=group)) + geom_polygon(color="black",fill="white") +
  coord_equal() + theme_nothing()
plot_nbhds_pl_df
plot <- ggplot(data=city_shp_df,aes(x=long,y=lat,group=group)) + geom_polygon(fill="lightblue")
plot
neighborhoods <- plot + geom_polygon(data=nbhds_pl_df,aes(x=long,y=lat,grou=group),color="gray",fill="white") +
  coord_equal() + theme_nothing()
neighborhoods
# crime data
crimeData <- read.csv("BPD_Part_1_Victim_Based_Crime_Data.csv",header=T,sep=",")
str(crimeData)
table(crimeData$description)
str(crimeData)
crimeData <- na.omit(crimeData)
# the coordinates of the crime points are given as text, so we need to convert
# the data in crimeData$Location.1
latlng <- str_replace(str_replace(crimeData$Location.1,'\\(',''),')','')
latlng <- str_split(latlng,', ')
latlng_df <- ldply(latlng)
str(latlng_df)
head(latlng_df)
crimeData$lat <- as.numeric(latlng_df[,1])
crimeData$long <- as.numeric(latlng_df[,2])
str(crimeData)
head(crimeData)
# convert lat/long coordinates to the original projection
latlng_df2 <- crimeData[,c("long","lat")]
latlng_spdf <- SpatialPoints(latlng_df2,proj4string=CRS("+proj=longlat +datum=WGS84"))
latlng_spdf <- spTransform(latlng_spdf,origProj)
proj4string(latlng_spdf)
str(latlng_spdf)
latlng_spdf_coords <- coordinates(latlng_spdf)
str(latlng_spdf_coords)
crimeData$long <- latlng_spdf_coords[,1]
crimeData$lat <- latlng_spdf_coords[,2]
str(crimeData)
# create a plot of 2d kernel density estimates for crime types
# for example : burglary
table(crimeData$description)
p <- neighborhoods + 
  geom_point(data=crimeData[crimeData$description=="BURGLARY",],aes(x=long,y=lat,group=1),shape="x",color="red",alpha=0.8)
p
p <- neighborhoods + 
  stat_density2d(data=crimeData[crimeData$description=="BURGLARY",],aes(x=long,y=lat,group=1,fill= ..level..,alpha=..level..),
                 size=2,bins=4,geom="polygon") 
p
p <- neighborhoods + 
  stat_density2d(data=crimeData[crimeData$description=="BURGLARY",],aes(x=long,y=lat,group=1,fill= ..level..,alpha=..level..),
                 size=2,bins=4,geom="polygon") + scale_fill_gradient(low="gray",high="orange") +
  annotate("text",x=1405000,y=565000,label="Burglary",size=6)
p
p <- neighborhoods + 
  stat_density2d(data=crimeData[crimeData$description=="AUTO THEFT",],aes(x=long,y=lat,group=1,fill= ..level..,alpha=..level..),
                 size=2,bins=4,geom="polygon") + scale_fill_gradient(low="gray",high="orange") +
  annotate("text",x=1405000,y=565000,label="Auto Teft",size=6)
p
p <- neighborhoods + 
  stat_density2d(data=crimeData[crimeData$description=="RAPE",],aes(x=long,y=lat,group=1,fill= ..level..,alpha=..level..),
                 size=2,bins=4,geom="polygon") + scale_fill_gradient(low="gray",high="orange") +
  annotate("text",x=1405000,y=565000,label="Rape",size=6)
p
p <- neighborhoods + 
  stat_density2d(data=crimeData,aes(x=long,y=lat,group=1,fill= ..level..,alpha=..level..),
                 size=2,bins=4,geom="polygon") + scale_fill_gradient(low="lightblue",high="darkred") +
  facet_wrap(~ description)
p
# visualizations of 2-dimensional kernel density estimates
# burglay incidents
burg_df <- subset(crimeData,description=="BURGLARY")
str(burg_df)
ggplot(data=burg_df,aes(x=long,y=lat)) + geom_point(color="pink",size=0.5) + coord_equal()
# aggravated assault incidents
aggass_df <- subset(crimeData,description=="AGG. ASSAULT")
str(aggass_df)
ggplot(data=aggass_df,aes(x=long,y=lat)) + geom_point(color="red",size=0.5) + coord_equal()
kde2dRange <- c(apply(burg_df[,c("long","lat")], 2,range))
kde2dRange
library(MASS)
getKde <- function(in_df, N=400, Lims=kde2dRange){
  pts <- as.matrix(in_df[,c('long','lat')])
  dens <- kde2d(pts[,1],pts[,2], n=N, lims=Lims)
  dens_df <- data.frame(expand.grid(dens$x, dens$y), z = c(dens$z))
  colnames(dens_df) <- c('x','y','z')
  return(dens_df)
}
plotKde2d <- function(in_df){
  fillCols <- rev(brewer.pal(11,'Spectral'))
  return(
    ggplot() + 
      geom_tile(data = in_df, aes(x=x, y=y, fill=z, group=1)) + 
      scale_fill_gradientn(colours=fillCols) + 
      theme_bw() +
      coord_equal()
  )
}
# get 2d kernel density estimate
burgDens <- getKde(burg_df)
summary(burgDens)
# plot 2d kernel density estimate
plotKde2d(burgDens)
plotKde2d(burgDens) + 
  geom_path(data=nbhds_pl_df,aes(x=long,y=lat,group=group),color='white',alpha=0.4) +
  ggtitle("Density map of burglary incidents")
kde2dRange <- c(apply(aggass_df[,c("long","lat")], 2,range))
kde2dRange
# get 2d kernel density estimate
aggassDens <- getKde(aggass_df)
summary(aggassDens)
# plot 2d kernel density estimate
library(lattice)
plotKde2d(aggassDens)
plotKde2d(aggassDens) + 
  geom_path(data=nbhds_pl_df,aes(x=long,y=lat,group=group),color='white',alpha=0.4) +
  ggtitle("Density map of aggravated assault incidents")

# Example 15
# heatmap of traffic signals Toronto
setwd(xxxxx)
library(ggmap)
library(ggplot2)
map <- get_map(location="Toronto",zoom=11)
Toronto <- ggmap(map,color="bw")
Toronto
rawdata <- read.csv("traffic_signals.csv",header=T,sep=",")
data <- data.frame(as.numeric(rawdata$Longitude),as.numeric(rawdata$Latitude))
names(data) <- c("lon","lat")
data <- na.omit(data)
str(data)
head(data)
p <- Toronto + geom_point(data=data,aes(x=lon,y=lat,col="red"))
p
p <- Toronto + geom_point(data=data,aes(x=lon,y=lat),col="red")
p
p <- p + stat_density2d(aes(x=lon,y=lat,fill=..level..,alpha=..level..),size=2,bins=4,data=data,geom="polygon") +
  scale_fill_gradient(low="yellow",high="darkred") + theme(legend.position="none")
p
geocode(location="Toronto")
Toronto <- c(lon=-79.38,lat=43.65)
map <- get_map(Toronto,zoom=11,color="bw")
Toronto_map <- ggmap(map,extent="panel",maprange=FALSE)
Toronto_map
p <- Toronto_map + geom_point(data=data,aes(x=lon,y=lat),col="blue",size=.4)
p
p <- p + geom_density2d(data=data,aes(x=lon,y=lat)) +
  stat_density2d(aes(x=lon,y=lat,fill=..level..,alpha=..level..),size=2,bins=4,data=data,geom="polygon") +
  scale_fill_gradient(low="yellow",high="darkred") + scale_alpha(range=c(0.00,0.40),guide=FALSE) +
  theme(legend.position="none")
p
p <- Toronto_map + geom_point(data=data,aes(x=lon,y=lat),col="blue",size=.4)
p
p <- p + geom_density2d(data=data,aes(x=lon,y=lat)) +
  stat_density2d(aes(x=lon,y=lat,fill=..level..,alpha=..level..),size=2,bins=4,data=data,geom="polygon") +
  scale_fill_gradient(low="yellow",high="darkred") + scale_alpha(range=c(0.00,0.40),guide=FALSE) +
  theme(legend.position="none")
p

# Example 16
# stat_summary 2D plot of points with lat/lon-coordinates
setwd(xxxxx)
NYdata <- read.csv("NYsample.csv",header=T,sep=";")
str(NYdata)
library(ggplot2)
# basic plot
pred.plot <- ggplot(data=NYdata,aes(x=lon,y=lat,colour=prediction)) + geom_point()
pred.plot
# binned data
# plot the mean in a 2D-region using stat_summary2D
pred.plot <- ggplot(data=NYdata,aes(x=lon,y=lat,z=prediction)) + stat_summary2d(fun=mean)
pred.plot
colormap <- c("Violet","Blue","Green","Yellow","Red","White")
pred.plot <- ggplot(data=NYdata,aes(x=lon,y=lat,z=prediction)) + 
  stat_summary2d(fun=median) +
  scale_fill_gradientn(name="Median",colours=colormap,space="Lab")
pred.plot
pred.plot <- ggplot(data=NYdata,aes(x=lon,y=lat,z=prediction)) + 
  stat_summary2d(fun=median,binwidth=c(0.010,0.010)) +
  scale_fill_gradientn(name="Median",colours=colormap,space="Lab")
pred.plot
library(ggmap)
loc <- get_map(location=c(min(NYdata$lon),min(NYdata$lat),max(NYdata$lon),max(NYdata$lat)),source="google")
theme_set(theme_bw(base_size=6))
map <- ggmap(loc) %+% NYdata +
  aes(x=lon,y=lat,z=prediction) +
  stat_summary2d(fun=median,binwidth=c(0.005,0.005),alpha=0.6) +
  scale_fill_gradientn(name="Median",colours=colormap,space="Lab") +
  labs(x="Longitude",y="Latitude") +
  coord_map()
map                    

# Example 17
# convert lat-lon coordinates / UTM coordinates (easting-northing)
setwd(xxxxx)
municipalities <- read.csv("coord_gemeenten_Belgie.csv",header=T,sep=";")
# conversion from lat/lon to UTM with xls-file latlon to UTM and vice versa
# see c:\R\Rdata\UTMconversions.xls
munB <- subset(municipalities,select=c(XCOORD,YCOORD))
write.csv(munB,file="c:/temp/munb.csv")
# import of c:/R/Rdata/munb_latlon_utm (latlon coordinates converted to utm)
munb <- read.csv("munb_latlon_utm.csv",header=T,sep=";")
str(munb)
head(munb)
plot(munb$xcoord,munb$ycoord)
plot(munb$easting,munb$northing,ylim=c(5450000,5720000),col="red")

# Example 18
# point pattern analysis of 311 requests NYC
setwd(xxxxx)
data <- read.csv("311_Requests_jan_2013.csv",header=T,sep=",")
str(data)
NYC311 <- subset(data,select=c(Unique.Key,Created.Date,Complaint.Type,Location.Type,Borough,X.Coordinate..State.Plane.,
                               Y.Coordinate..State.Plane.,Latitude,Longitude,Location))
# rename columns
names(NYC311) [names(NYC311)=="Created.Date"] <- c("cdate")
names(NYC311) [names(NYC311)=="X.Coordinate..State.Plane."] <- c("xcoord")
names(NYC311) [names(NYC311)=="Y.Coordinate..State.Plane."] <- c("ycoord")
str(NYC311)
# delete records with missing value(s)
NYC311 <- na.omit(NYC311)
dim(NYC311)
# date and time handling
library(lubridate)
library(plyr)
NYC311$cdate1 <- strptime(NYC311$cdate,format="%m/%d/%Y %I:%M:%S %p")
head(NYC311$cdate1)
NYC311$date <- as.Date(NYC311$cdate1,format="%d/%m/%Y")
NYC311$Year <- year(NYC311$date)
NYC311$Month <- month(NYC311$date)
NYC311$Weekday <- wday(NYC311$date)
NYC311$Monthday <- mday(NYC311$date)
NYC311$Monthweek <- ceiling(NYC311$Monthday/7)
NYC311$month <- month(NYC311$date,label=TRUE,abbr=TRUE)
NYC311$weekday <- wday(NYC311$date,label=TRUE,abbr=TRUE)
NYC311$wday <- factor(NYC311$weekday,order=TRUE,levels=c("Mon","Tues","Wed","Thurs","Fri","Sat","Sun"))
NYC311$wday <- revalue(NYC311$wday,c("Mon"="Mon","Tues"="Tue","Wed"="Wed","Thurs"="Thu","Fri"="Fri","Sat"="Sat","Sun"="Sun"))
NYC311$cdate2 <- as.character(NYC311$cdate1)
NYC311$hour <- substr(NYC311$cdate2,12,13)
NYC311$hour <- as.numeric(NYC311$hour)
str(NYC311)
# limiting the data to january 2013
str(NYC311)
plot(NYC311$Longitude,NYC311$Latitude)
plot(NYC311$xcoord,NYC311$ycoord)
library(hexbin)
bin <- hexbin(NYC311$Longitude,NYC311$Latitude,xbin=100)
plot(bin)
library(IDPmisc)
iplot(NYC311$Longitude,NYC311$Latitude)
iplot(NYC311$Longitude,NYC311$Latitude,pixs=2)
mean(NYC311$Longitude)
mean(NYC311$Latitude)
library(ggmap)
NYC <- c(lon=-73.92507,lat=40.73083)
map <- get_map(NYC,zoom=11,color="bw")
map <- ggmap(map,extent="panel")
map
library(ggplot2)
map <- map + geom_point(data=NYC311,aes(x=Longitude,y=Latitude),alpha=0.2,color="pink")
map
map <- map + 
  stat_density2d(data=NYC311,aes(x=Longitude,y=Latitude,group=1,fill=..level..,alpha=..level..),
                 size=2,bins=4,geom="polygon") +
  scale_fill_gradient(low="gray",high="orange")
map
map <- map + 
  stat_density2d(data=nyc311,aes(x=Longitude,y=Latitude,group=1,fill=..level..,alpha=..level..),
                 size=2,bins=4,geom="polygon") +
  scale_fill_gradient(low="gray",high="orange")
map
library(rgdal)
NYC <- readOGR(".","nymcwi")
plot(NYC)
proj4string(NYC)
origProj <- NYC@proj4string
# to work with ggplot we fortify the shapefile
NYC_df <- fortify(NYC)
str(NYC_df)
baseplot <- ggplot(data=NYC_df,aes(x=long,y=lat,group=group)) + geom_polygon(color="gray",fill="lightblue") +
  coord_equal() + theme_nothing()
baseplot
# convert the lat/long coordinates of 311 datafile to the original projection
str(nyc311)
nyc311_spdf <- SpatialPoints(nyc311,proj4string=CRS("+proj=longlat +datum=WGS84"))
class(nyc311_spdf)
nyc311_spdf <- spTransform(nyc311_spdf,origProj)
proj4string(nyc311_spdf)
nyc311_spdf_coords <- coordinates(nyc311_spdf)
str(nyc311_spdf_coords)
nyc311$long <- nyc311_spdf_coords[,1]
nyc311$lat <- nyc311_spdf_coords[,2]
str(nyc311_)
p <- baseplot + geom_point(data=nyc311,aes(x=long,y=lat,group=1),shape="x",color="red",alpha=0.8)
p
# visualization of 2-dimensional kernel density estimates
kde2dRange <- c(apply(nyc311[,c("long","lat")],2,range))
kde2dRange
library(MASS)
getKde <- function(in_df, N=400, Lims=kde2dRange){
  pts <- as.matrix(in_df[,c('long','lat')])
  dens <- kde2d(pts[,1],pts[,2], n=N, lims=Lims)
  dens_df <- data.frame(expand.grid(dens$x, dens$y), z = c(dens$z))
  colnames(dens_df) <- c('x','y','z')
  return(dens_df)
}
plotKde2d <- function(in_df){
  fillCols <- rev(brewer.pal(11,'Spectral'))
  return(
    ggplot() + 
      geom_tile(data = in_df, aes(x=x, y=y, fill=z, group=1)) + 
      scale_fill_gradientn(colours=fillCols) + 
      theme_bw() +
      coord_equal()
  )
}
# get 2d kernel density estimate
nyDens <- getKde(nyc311)
summary(nyDens)
library(RColorBrewer)
# plot 2d kernel density estimate
plotKde2d(nyDens)
plotKde2d(nyDens) + 
  geom_path(data=NYC_df,aes(x=long,y=lat,group=group),color='white',alpha=0.4) +
  ggtitle("Density map of 311 requests New York City, January 2013") +
  annotate("text",x=1000000,y=130000,label="NYC Open Data",size=4.5)

# Example 19
# point patterns can also be visualized with package "squash"
# this package provides functions for color-based visualization of multivariate data
# base question : given a large number of 3-dimensional points (x,y,z),
# how does z vary as a function of x and y ?
setwd(xxxxx)
catalog <- read.csv("earthquake-catalog.csv",header=T,sep=",",stringsAsFactors=FALSE)
# date and time handling
library(lubridate)
library(plyr)
catalog$time <- sub("T"," ",catalog$time)
catalog$time <- sub("Z","",catalog$time)
catalog$time <- strptime(catalog$time, format="%Y-%m-%d %H:%M:%S")
catalog$date <- as.Date(catalog$time)
catalog$year <- year(catalog$time)
str(catalog)
library(squash)
attach(catalog)
# the color in the squashgram indicates a summary (in this case the mean) of all z values
# of the points falling into the bin
squashgram(depth ~ longitude + latitude, FUN=mean,main = "Depth of earthquakes (2004-2014)")
# a larger quare indicates more points faling into the rectangular interval
squashgram(depth ~ longitude + latitude, FUN=mean, shrink=10,main = "Depth of earthquakes (2004-2014)")

# Example 20
# Dot density maps in R
# www.flowingdata.com
# write introduction in guide based on document flowingdata.com
# create a dot density map for population in California
# source : American Community Survey at the census tract level
setwd(xxxxx)
# data Nebraska
nepop <- read.csv("ACS_12_5YR_B01003_with_ann.csv",stringsAsFactors=FALSE,colClasses=acsClasses,sep=",")
str(nepop)
# we only need 2 colums : the census tract geographic id to link with shapefile data and
# population estimate
nepop.sub <- nepop[,c(2,4)]
head(nepop.sub)
str(nepop.sub)
# Nebraska shapefile for census tracts
netract.shp <- readShapePoly("tl_2013_31_tract.shp",proj4string=CRS("+proj=longlat"))
class(netract.shp)
nepolys <- SpatialPolygonsDataFrame(netract.shp,data=as(netract.shp,"data.frame"))
class(nepolys)
summary(nepolys)
str(nepolys@data)
head(nepolys@data)
plot(nepolys)
# merge datasets
nedata <- merge(nepolys@data,nepop.sub,by.x="GEOID",by.y="geoid2",sort=FALSE)
str(nedata)
head(nedata)
# we will make a dot density map of "totalpop" and define one dot per 1000 people
plotvar <- nedata$totalpop/100
# generate randon dots in polygons (census tracts)
nedots.rand <- dotsInPolys(nepolys,as.integer(plotvar),f="random")
str(nedots.rand)
nedots.rand@coords
# now we can plot the counties and draw the dots on top of the existing county map
par(mar=c(0,0,0,0))
plot(nepolys,lwd=0.1)
plot(nedots.rand,add=TRUE,pch=19,cex=0.1,col="#00880030")
# in many parts of the state of Nebraska, the population is sparse
# therefore we will focus on a single area to see more than a blob of dots
# plot of a specific region (zoom on Omaha region)
plot(nepolys,lwd=0.1,xlim=c(-96.187303,-95.847993),ylim=c(41.073844,41.388218))
plot(nedots.rand,add=TRUE,pch=19,cex=0.1,col="#008800")
# grid of dots instead of random
nedots.reg <- dotsInPolys(nepolys, as.integer(plotvar), f="regular")
plot(nepolys, lwd=0.1, xlim=c(-96.187303,-95.847993), ylim=c(41.073844,41.388218))
plot(nedots.reg, add=TRUE, pch=19, cex=0.3, col="#008800")
# color options
randomCols <- colors()[sample(1:502, 12)]
par(mfrow=c(2,3), mar=c(0,0,1,0))
for (i in 1:6) {
  plot(nepolys, lwd=0.1, xlim=c(-96.187303,-95.847993), ylim=c(41.073844,41.388218), col=randomCols[i])
  plot(nedots.reg, add=TRUE, pch=19, cex=0.1, col=randomCols[6+i])
}
# using color to show another dimension
# we can color dots by category
# for example, for the state of Nebraska we have data of the population by race at the tract level
nerace <- read.csv("race-nebraska-truncated.csv", 
                   stringsAsFactors=FALSE, 
                   colClasses=c("character", "numeric", "numeric", "numeric", "numeric"))
nedata <- merge(nepolys@data, nerace, by.x="GEOID", by.y="geoid2", sort=FALSE)
races <- c("white", "black", "amind", "asian")
dotCols <- c("#da6678", "#647eee", "#c29219", "#09900d")
par(mar=c(0,0,0,0))
plot(nepolys, lwd=0.2, xlim=c(-96.187303,-95.847993), ylim=c(41.073844,41.388218))
for (i in 2:length(races)) {
  nevar <- nedata[,races[i]] / 10
  nedots.race <- dotsInPolys(nepolys, as.integer(nevar), f="random")
  plot(nedots.race, add=TRUE, pch=19, cex=0.01, col=dotCols[i])
}
# for whites only :
plot(nepolys, lwd=0.2, xlim=c(-96.187303,-95.847993), ylim=c(41.073844,41.388218))
for (i in 1:1) {
  nevar <- nedata[,races[i]] / 10
  nedots.race <- dotsInPolys(nepolys, as.integer(nevar), f="random")
  plot(nedots.race, add=TRUE, pch=19, cex=0.01, col=dotCols[i])
}
# dot density plot of race with ggplot 2
nemap <- fortify(nepolys)
map <- ggplot(nemap,aes(x=long,y=lat)) + geom_polygon(aes(group=group),color=I("grey65"),fill="lightblue") + coord_equal()
map
str(nedata)
dots.bl <- dotsInPolys(nepolys,as.integer(nedata$black/10))
dots.bl$ethnicity <- "black"
dots.as <- dotsInPolys(nepolys,as.integer(nedata$asian/10))
dots.as$ethnicity <- "asian"
dots.am <- dotsInPolys(nepolys,as.integer(nedata$amind/10))
dots.am$ethnicity <- "amind"
# spRbind takes only two parameters at a time
dots.all <- spRbind(dots.bl,dots.as)
dots.all <- spRbind(dots.all,dots.am)
str(dots.all)
ethno.df <- data.frame(coordinates(dots.all)[,1:2],ethnicity=dots.all$ethnicity)
str(ethno.df)
df <- data.frame(coordinates(cadots.rand)[,1:2])
str(df)
map <- map + geom_point(data=ethno.df,aes(x=x,y=y,colour=factor(ethnicity)),size=0.1) +
  scale_colour_manual(values=c("grey","cyan","magenta"))
map

# Example 21
setwd(xxxxx)
library(raster)
library(maptools)
NYC <- readShapePoly("nybb.shp")
plot(NYC)
proj4string(NYC)
str(NYC@data)
class(NYC)
import <- raster("SF2013_raster.tif")
plot(import)
a <- flip(import,direction="y")
plot(a)
# crop and mask raster file
r2 <- crop(a,extent(NYC))
plot(r2)
r3 <- mask(r2,NYC)
plot(r3)
plot(NYC,add=TRUE,lwd=1)
# the best way to extract coordinates from a SPDF is to use ggplot2
library(ggplot2)
dum = fortify(NYC)
str(dum)
map <- ggplot(dum, aes(x = long, y = lat,group=group)) + geom_path()
map
# convert raster to dataframe
b <- rasterToPoints(a)
df <- data.frame(b)
str(df)
colnames(df) <- c("X","Y","layer")
# convert raster to dataframe
b <- rasterToPoints(r3)
df <- data.frame(b)
str(df)
colnames(df) <- c("X","Y","layer")
map <-   ggplot(df)+
  geom_raster(data=df,aes(X,Y,fill=layer))+
  scale_fill_gradientn(name="Stop and Frisk",colours = rainbow(20))
map
map + geom_path(data=dum,aes(x=long,y=lat,group=group)) +
  coord_equal() 

# Example 23
# concentric circles as points in ggplot2 / ggmap
setwd(xxxxx)
concentric <- read.csv("concentric.csv",header=T,sep=";")
str(concentric)
head(concentric)
library(ggplot)
library(ggmap)
library(maps)
# with ggplot2 :
all_states <- map_data("state")
str(all_states)
p <- ggplot()
p <- ggplot() + geom_polygon(data=all_states,aes(x=long,y=lat,group=group),colour="grey",fill="white") +
  coord_equal()
p
# the basic idea is to use separate geoms for two variables, making sure that the smaller one is plotted
# after the larger one
# add total population
p <- p + geom_point(data=concentric,aes(x=longitude,y=latitude,size=totalPop),colour="#b5e521")
p
# add subpopulation
p <- p + geom_point(data=concentric,aes(x=longitude,y=latitude,size=subPop),colour="#00a3e8")
p
# change name of legend
p <- p + guides(size=guide_legend(title="Population"))
p
# with ggmap :
map <- get_map(location="united states",zoom=3,maptype="terrain",source="google")
p <- ggmap(map)
p
# add total population
p <- p + geom_point(data=concentric,aes(x=longitude,y=latitude,size=totalPop),colour="#D55E00")
p
# add subpopulation
p <- p + geom_point(data=concentric,aes(x=longitude,y=latitude,size=subPop),colour="#00a3e8")
p
# change name of legend
p <- p + guides(size=guide_legend(title="Population"))
p

# Example 25
# Combine data and spatial locations
# analysis of air quality data from several measurement stations which measure pollution in
# different areas of the urban environment
# source : https://github.com/oscarperpinan/spacetime-vis
setwd(xxxxx)
library(sp)
## spatial location of stations
airStations <- read.csv("airStations.csv",header=T,sep=";")
str(airStations)
head(airStations)
# the dataframe is converted to a SpatialPointsDataFrame object
coordinates(airStations) <- ~ long + lat
# geographical projection
proj4string(airStations) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
str(airStations)
# measurements data
airQuality <- read.csv("airQuality.csv",header=T,sep=";")
str(airQuality)
# only interested in NO2 
NO2 <- airQuality[airQuality$codParam==8, ]
head(NO2)
# aggregate data
NO2agg <- aggregate(dat ~ codEst, data=NO2,
                    FUN = function(x) {
                      c(mean=signif(mean(x), 3),
                        median=median(x),
                        sd=signif(sd(x), 3))
                    })
NO2agg
NO2agg <- do.call(cbind, NO2agg)
NO2agg <- as.data.frame(NO2agg)
str(NO2agg)
library(maptools)
library(RColorBrewer)
# link aggregated data with stations to obtain a SpatialPointsDataFrame.
# Codigo and codEst are the stations codes
idxNO2 <- match(airStations$Codigo, NO2agg$codEst)
idxNO2
# the aggregated data (a data frame) and the spatial information (a SpatialPointsDataFrame)
# are combined with the spCbind method from the maptools package
NO2sp <- spCbind(airStations[, c('Nombre', 'alt')], NO2agg[idxNO2, ])
str(NO2sp)
save(NO2sp, file="NO2sp.RData")
load("NO2sp.RData")
# proportional symbol mapping
library(classInt)
plotvar <- NO2sp@data$mean
nclr <- 5
plotclr <- brewer.pal(nclr,"Greens")
class <- classIntervals(plotvar,nclr,style="quantile")
class
spplot(NO2sp["mean"],col.regions=plotclr,at=round(class$brks,digits=2),cex=sqrt(1:5))
class <- classIntervals(plotvar,nclr,style="equal")
class
spplot(NO2sp["mean"],col.regions=plotclr,at=round(class$brks,digits=2),cex=sqrt(1:5))
spplot(NO2sp,c("mean"),col.regions=bpy.colors(10),edge.col="black",scales=list(draw=TRUE),
       key.space="right",cex=sqrt(1:5))
# to plot the data with ggplot2 we need to transform the SpatialPointsDataFrame
# into a conventional data frame (which will contain two columns with
# latitude and longitude values)
airPal <- colorRampPalette(c("springgreen1","sienna3","gray5"))(5)
NO2df <- data.frame(NO2sp)
str(NO2df)
NO2df$Mean <- cut(NO2df$mean,5)
library(ggplot2)
p <- ggplot(data=NO2df,aes(x=long,y=lat,size=Mean,fill=Mean)) +
  geom_point(pch=21,col="black") +
  theme_bw() +
  scale_fill_manual(values=airPal)
p
# two main improvements can be added :
# define classes dependent on the data structure (instead of a uniform distribution assumed with cut)
# encode each group with a symbol size (circle area) such that visual discrimination is enhanced
library(classInt)
nClasses <- 5
intervals <- classIntervals(NO2sp@data$mean,n=nClasses,style="fisher")
tab <- print(intervals)
dent <- c(0.64,1.14,1.65,2.79,4.32,6.22,9.65,12.95,15.11)
dentAQ <- dent[seq_len(nClasses)]
idx <- findCols(intervals)
cexNO2 <- dentAQ[idx]
NO2sp$classNO2 <- factor(names(tab)[idx])
# ggplot2 version :
NO2df <- data.frame(NO2sp)
ggplot(data=NO2df,aes(x=long,y=lat,size=classNO2,fill=classNO2)) +
  geom_point(pch=21,col="black") +
  theme_bw() +
  scale_fill_manual(values=airPal) +
  scale_size_manual(values=dentAQ*2)
# spplot version : 
spplot(NO2sp["classNO2"],col.regions=airPal,cex=dentAQ,edge.col="black",scales=list(draw=TRUE))
# improved legend :
NO2key <- list(x=0.98, y=0.02,corner=c(1,0),title=expression(NO[2]~~(paste(mu,plain(g))/m^3)),
               cex.title=.75,cex=0.7,background="gray92")
pNO2 <- spplot(NO2sp["classNO2"],col.regions=airPal,cex=dentAQ,edge.col="black",scales=list(draw=TRUE),
               key.space=NO2key)
pNO2
# the spatial distribution of points is better understood if we add underlying layers with
# information about the spatial context
# a suitable method is to make use of the ggmap package
str(NO2sp)
madridBox <- bbox(NO2sp)
madridBox
library(ggmap)
madridGG <- get_map(c(madridBox),source="google")
madridMap <- ggmap(madridGG)
madridMap
madridMap + 
  geom_point(data=NO2df,aes(x=long,y=lat,size=classNO2,fill=classNO2),pch=21,col="black") +
  scale_fill_manual(values=airPal) +
  scale_size_manual(values=dentAQ*2)
# the result of ggmap is a raster object with attributes
# it can be displayed with grid.raster as an underlying layer of the previous spplot result
library(grid)
library(latticeExtra)
bbMap <- attr(madridGG,"bb")
height <- with(bbMap,ur.lat - ll.lat)
width <- with(bbMap,ur.lon - ll.lon)
pNO2 + layer(grid.raster(madridGG,width=width,height=height,default.units="native"),under=TRUE)
# spatial interpolation
# the measurement at discrete points give limited information about the underlying process
# it is quite common to approximate the spatial distribution of the measured variable with
# the interpolation between measurement locations
str(NO2df)
d <- NO2df
# use of interp() function of akima package
library(akima)
akima_li <- interp(x=d$long,
                   y=d$lat,
                   z=d$mean,
                   yo=seq(min(d$lat), max(d$lat), length=250),
                   xo=seq(min(d$long), max(d$long), length=250),
                   linear=FALSE,
                   extrap=TRUE)
dInterp <- data.frame(expand.grid(x=akima_li$x, y=akima_li$y), z=c(akima_li$z))
str(dInterp)
plot <- ggmap(madridGG) + geom_tile(data=dInterp, aes(x, y, fill=z), alpha=0.5, color=NA) +
  scale_fill_gradientn(colours=brewer.pal(9,"Reds")) +
  labs(fill="") +
  theme_nothing(legend=TRUE)
plot
# using inverse distance weighted interpolation (IDW) with the gstat package
library(gstat)
class(NO2sp)
proj4string(NO2sp)
airGrid <- spsample(NO2sp,type="regular",n=1e5)
gridded(airGrid) <- TRUE
airKrige <- krige(mean ~ 1,NO2sp,airGrid)
NO2sp_t <- spTransform(NO2sp, CRS=CRS("+proj=merc +ellps=WGS84"))
str(NO2sp_t)
airGrid <- spsample(NO2sp_t,type="regular",n=1e5)
gridded(airGrid) <- TRUE
airKrige <- krige(mean ~ 1,NO2sp_t,airGrid)
str(airKrige)
class(airKrige)
spplot(airKrige["var1.pred"])    
# shapefile districts of Madrid
library(rgdal)
madridshapefile <- readOGR(".","Distritos13")
plot(madridshapefile)
proj4string(madridshapefile)
proj4string(madridshapefile) <- CRS("+proj=utm +zone30")
distritosMadrid <- spTransform(madridshapefile,CRS=CRS("+proj=longlat +ellps=WGS84"))
library(grid)
library(latticeExtra)
spplot(airKrige["var1.pred"],col.regions=colorRampPalette(airPal)) +
  layer({
    sp.polygons(distritosMadrid,fill="transparent",lwd=.3)
    sp.points(NO2sp,pch=21,alpha=0.8,fill="gray50",col="white")
  })

# Example 26
# spatiotemporal point observations
setwd(xxxxx)
library(sp)
## spatial location of stations
airStations <- read.csv("airStations.csv",header=T,sep=";")
# the dataframe is converted to a SpatialPointsDataFrame object
coordinates(airStations) <- ~ long + lat
# geographical projection
proj4string(airStations) <- CRS("+proj=longlat +ellps=WGS84 +datum=WGS84")
str(airStations)
# measurements data
airQuality <- read.csv("airQuality.csv",header=T,sep=";")
str(airQuality)
# only interested in NO2 
NO2 <- airQuality[airQuality$codParam==8, ]
head(NO2)
library(zoo)
library(spacetime)
# in NO2 each line corresponds with a measurement at one station during a day of the year
# we need a data frame with each station as a different variable during a day of the year
# this means we have to convert the data from long to wide
NO2$time <- with(NO2,ISOdate(year,month,day))
NO2 <- subset(NO2,select=c(codEst,dat,time))
str(NO2)
head(NO2)
library(reshape2)
NO2wide <- dcast(NO2,time ~ codEst,value.var="dat")
str(NO2wide)
head(NO2wide)
NO2zoo <- zoo(NO2wide[,-1],NO2wide$time)
str(NO2zoo)
head(NO2zoo)
class(NO2zoo)
dats <- data.frame(vals=as.vector(t(NO2zoo)))
str(dats)
head(dats)
class(dats)
NO2st <- STFDF(airStations,index(NO2zoo),dats)
str(NO2st)
class(NO2st)
airPal <- colorRampPalette(c("springgreen","sienna3","gray5"))(5)
stplot(NO2st[,1:12],cuts=5,col.regions=airPal,edge.col="black")
stplot(NO2st,mode="xt",col.regions=colorRampPalette(airPal)(15),
       scales=list(x=list(rot=45)),xlab="",ylab="")

