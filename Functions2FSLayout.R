## Functions used to analyze the data regarding 
## the paper Bispo et al. (2023)  Using spatial point process models, 
## clustering and space partitioning to reconfigure fire stations layout


## Function that identifies spatial points falling outside the district boundaries
outlimitdist<-function(name){
   # name="Aveiro"
Dist<-dist[dist$NAME_1==name,]
spdf_error <- SpatialPointsDataFrame(coords = Data[Data$Distrito==name, c("Longitude", "Latitude")], data = Data[Data$Distrito==name,], proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
out<-sum(is.na(over(spdf_error, as(Dist, "SpatialPolygons"))))
result<-list() 
if(out==0) {
 result<-print(paste0(name,": all points inside bounderies"))
 } else {
 points_out <- spdf_error[is.na(over(spdf_error, as(Dist, "SpatialPolygons"))), ]
 dfout<-data.frame(dist_error=rep(name,length(points_out$Longitude)),Longitude=points_out$Longitude,Latitude=points_out$Latitude)
 spdf_true<-SpatialPointsDataFrame(coords = dfout[,c("Longitude", "Latitude")], data = dfout, proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
 dfout$dist_true<-as.character(dist$NAME_1[over(spdf_true, as(dist, "SpatialPolygons"))])
 result<-dfout
 }
return(result)
}

## Function to find the grid neighborhood for a number of FS given by the actual existent FS
ncell<-function(n,d=dens){#n is the true FS in the district
r <- raster(dens)
nmax<-numeric()
s=seq(3,31,2)
for (i in 1:length(s)) {
  r2 <- r==focal(raster(dens), matrix(1,s[i],s[i]), fun = function(X) max(X, na.rm=TRUE))
  maxXY <- xyFromCell(r2, Which(r2==1, cells=TRUE))
  nmax[i]<-nrow(maxXY)
}
result<-s[which.min(abs(nmax-n))]
return(result)
}


## Function that define the border for districts
dborder<-function(df){
boundD<-dist[dist$NAME_1==paste(unique(df$Distrito)),]
spdf <- SpatialPointsDataFrame(coords = df[, c("LongitudeCB", "LatitudeCB")], data = df, proj4string = CRS("+proj=longlat +datum=WGS84 +no_defs +ellps=WGS84 +towgs84=0,0,0"))
df <- df %>%  filter(!is.na(over(spdf, as(boundD, "SpatialPolygons"))))
}


## Function to set the geometries for each district
setregion<-function(df){
r<-lapply(list(slot(df, "polygons")),SpatialPolygons)
window<-lapply(r, as.owin)
return(window)
}

## Function to create a point pattern
pppfun<-function(df,w){ppp(df$Longitude,df$Latitude,window=w)}

## Function from https://carsonfarmer.com/2009/09/voronoi-polygons-with-r/
## To create a nice bounded Voronoi polygons tessellation of a point layer in R, 
## we need two libraries: sp and deldir. The following function takes a 
## SpatialPointsDataFrame as input, and returns a SpatialPolygonsDataFrame 
## that represents the Voronoi tessellation of the input point layer.

voronoipolygons <- function(x, poly) {
  require(deldir)
  if (.hasSlot(x, 'coords')) {
    crds <- x@coords  
  } else crds <- x
  bb = bbox(poly)
  rw = as.numeric(t(bbox(poly)))
  z <- deldir(crds[,1], crds[,2],rw=rw)
  w <- tile.list(z)
  polys <- vector(mode='list', length=length(w))
  require(sp)
  for (i in seq(along=polys)) {
    pcrds <- cbind(w[[i]]$x, w[[i]]$y)
    pcrds <- rbind(pcrds, pcrds[1,])
    polys[[i]] <- Polygons(list(Polygon(pcrds)), ID=as.character(i))
  }
  SP <- SpatialPolygons(polys)
  SpatialPolygonsDataFrame(
    SP, data.frame(x=crds[,1], y=crds[,2], 
                   row.names=sapply(slot(SP, 'polygons'), 
                                    function(x) slot(x, 'ID'))))  
}


## Functions used to calculate new and true average dist.
## 1)
newdist<-function(df.reg,city,maxima){
coordOPT<-data.frame(lon=maxima[,1],lat=maxima[,2])
events<-data.frame(lon=unique(city)$Longitude,lat=unique(city)$Latitude)
d<-pointDistance(events[,c(1,2)],coordOPT,lonlat=FALSE) 
dm <- as.matrix(d)
dmin <- apply(dm, 1, min, na.rm=TRUE)
wdmin <- apply(dm, 1, which.min)
tesselation <- deldir(coordOPT$lon,coordOPT$lat)
tiles <- tile.list(tesselation)
distances<-data.frame(lat=city$Latitude,lon=city$Longitude,tile=wdmin)
centers<-data.frame(coordOPT,tile=1:nrow(coordOPT))
res<-merge(distances,centers,by="tile")
for (i in 1:nrow(res)){res$dnew[i]=distm(c(res$lon.x[i],res$lat.x[i]),c(res$lon.y[i],res$lat.y[i]),fun=distHaversine)/1000}#distancias dos eventos aos novos clusters
return(res$dnew)
}

## 2)
distTP<-function(reg,city,max){
df<-data.frame(conc=city$Concelho,lat=city$Latitude,lon=city$Longitude, actual=city$euc_dist)#distancias reais
df$new<-newdist(df.reg = reg,city = city,maxima = max)
result<-round(colMeans(df[,c("actual","new")]),2)
return(result)
}




