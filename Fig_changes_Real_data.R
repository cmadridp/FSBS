
library(ncdf4) # package for netcdf manipulation
library(raster) # package for raster manipulation
library(ggplot2) # package for plotting
library(lubridate)
library(fields)

this.dir = dirname(rstudioapi::getActiveDocumentContext()$path)
setwd(this.dir)

data=stack("sst.mon.mean.nc")

m1=6 #June
i.year=1940
f.year=2019
t.years=f.year-i.year+1
iyear.pos=i.year-1850+1

m1.ipos=(iyear.pos-1)*12+m1 #m1 from i.year
m1.indices=sequence(t.years,from=m1.ipos,by=12)

m1.data=data[[m1.indices]]

#extracting the data
data=stack("sst.mon.mean.nc")
data1=nc_open("sst.mon.mean.nc")
lat=ncvar_get(data1,'lat')
lon=ncvar_get(data1,'lon')
time=ncvar_get(data1,'time')
time2=as.Date('1891-01-01')+time
sst=ncvar_get(data1,"sst")
image.plot(lon,rev(lat),sst[,length(lat):1,m1.indices[80]],xlab='Longitud',ylab='Latitude')
title("Sea Surface Temperature April 2019")
rect(xleft=251,xright = 300,ybottom=1,ytop=50, density=0, col = "green") 
rect(xleft=153.5,xright = 162.5,ybottom=-28.5,ytop=-24.5, density=0, col = "black",lwd = 2.) #Regions of interest 
rect(xleft=110.5,xright = 157.5,ybottom=-39.5,ytop=-10.5, density=0, col = "black",lwd=2) #Regions of interest 


#plot of australia
lon_australia=matrix(NaN,48,1)
aux=-1
for(i in 1:48){
  aux=aux+1
  lon_australia[i]=110.5+aux
}
lat_australia=matrix(NaN,30,1)
aux=-1
for(i in 1:30){
  aux=aux+1
  lat_australia[i]=-10.5-aux
}
image.plot(lon_australia,rev(lat_australia),sst[111:158,130:101,m1.indices[80]],zlim=c(12,28.6),xlab='Longitude',ylab='Latitude',xlim=c(110.5,157.5))
title("SST Australia June 2019")


########### Plotting the estimated periods #####################################
aust.sst.jun=sst[111:158,130:101,m1.indices]

means.map1=matrix(NA, nrow=dim(aust.sst.jun)[1], ncol=dim(aust.sst.jun)[2])
for(i in 1:(dim(means.map1)[1])){
  for(j in 1: (dim(means.map1)[2])){
    means.map1[i,j]=mean(aust.sst.jun[i,j,1:43])
  }
}

image.plot(lon_australia,rev(lat_australia),means.map1,zlim=c(12,28.6),xlab='Longitude',ylab='Latitude',xlim=c(110.5,157.5))
title("Average SST Australia June 1940-1982")


means.map2=matrix(NA, nrow=dim(aust.sst.jun)[1], ncol=dim(aust.sst.jun)[2])
for(i in 1:(dim(means.map2)[1])){
  for(j in 1: (dim(means.map2)[2])){
    means.map2[i,j]=mean(aust.sst.jun[i,j,44:58])
  }
}

image.plot(lon_australia,rev(lat_australia),means.map2,zlim=c(12,28.6),xlab='Longitude',ylab='Latitude',xlim=c(110.5,157.5))
title("Average SST Australia June 1983-1997")

means.map3=matrix(NA, nrow=dim(aust.sst.jun)[1], ncol=dim(aust.sst.jun)[2])
for(i in 1:(dim(means.map3)[1])){
  for(j in 1: (dim(means.map3)[2])){
    means.map3[i,j]=mean(aust.sst.jun[i,j,59:80])
  }
}

image.plot(lon_australia,rev(lat_australia),means.map3,zlim=c(12,28.6),xlab='Longitude',ylab='Latitude',xlim=c(110.5,157.5))
title("Average SST Australia June 1998-2019")


###############################################################################################################
