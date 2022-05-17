rm(list = ls())

library(fda)
library(ks)
library(mvtnorm)
library(ncdf4) # package for netcdf manipulation
library(raster) # package for raster manipulation
library(ggplot2) # package for plotting
library(lubridate)
library(fields)

source("FSBS.R")
source("Auxiliar_Functions_FSBS.R")


#The one I want
data=stack("sst.mon.mean.nc")
data1=nc_open("sst.mon.mean.nc")
lat=ncvar_get(data1,'lat')
lon=ncvar_get(data1,'lon')
time=ncvar_get(data1,'time')
time2=as.Date('1891-01-01')+time
sst=ncvar_get(data1,"sst")


#Sub-data by time
m1=6 #June
i.year=1940
f.year=2019
t.years=f.year-i.year+1
iyear.pos=i.year-1850+1

m1.ipos=(iyear.pos-1)*12+m1 #m1 from i.year
m1.indices=sequence(t.years,from=m1.ipos,by=12)

m1.data=data[[m1.indices]]





image.plot(lon,rev(lat),sst[,length(lat):1,m1.indices[80]],xlab='Longitud',ylab='Latitude')
title("Sea Surface Temperature April 2019")
rect(xleft=251,xright = 300,ybottom=1,ytop=50, density=0, col = "green") #Red indicates the reference region to show in paper
rect(xleft=153.5,xright = 162.5,ybottom=-28.5,ytop=-24.5, density=0, col = "black",lwd = 2.) #Regions of interest corresponding to the Caribbean Sea 
rect(xleft=110.5,xright = 157.5,ybottom=-39.5,ytop=-10.5, density=0, col = "black",lwd=2) #Regions of interest corresponding to the Caribbean Sea 



#plot only close to australia
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
image.plot(lon_australia,rev(lat_australia),sst[111:158,130:101,m1.indices[80]],xlab='Longitude',ylab='Latitude',xlim=c(110.5,157.5))
title("SST Australia June 2019")

#The region of interest is a 10x5 grid, i.e. n=50
lat_grid=seq(-39.5,-10.5,1)
lon_grid=seq(110.5,157.5,1)
#lat_grid=c(-28.5,-27.5,-26.5,-25.5,-24.5)
#lon_grid=c(153.5,154.5,155.5,156.5,157.5,158.5,159.5,160.5,161.5,162.5)
matrix_lon_lat=matrix=matrix(NaN,1440,2)
for(i in 0:29){
  matrix_lon_lat[(i*48+1):((i+1)*48),]=lon_grid
  matrix_lon_lat[(i*48+1):((i+1)*48),2]=lat_grid[i+1]
}

#as.vector(sst[104:113,119:115,m1.indices[50]])
#sst[154:163,119:110,m1.indices[50]]
aux_matrix_sst=matrix(NaN,1440,80)
for(i in 1:80){
  for(j in 1:30){
    for(k in 1:48){
      if(is.na(sst[111:158,130:101,m1.indices[i]][k,j])){
        sst[111:158,130:101,m1.indices[i]][k,j]=0
      }
    }
  }
  
}

#-------Putting the data together as in my format-------#

Tb=80
r=1
C=3
#the responses lies in dimension 3
n=1440
#the observations lies in dimension 1
d=2
#we define the n_ts
nts<-matrix(NaN,nrow=Tb, ncol=1)
for(i in 1:Tb){
  nts[i]=n
}
N_t=sum(nts)
N_t

#creating X_s
X_s=matrix(,0,2)
for(aux in 1:80){
  X_s=rbind(X_s,matrix_lon_lat)
}
#creating y_t_is
y_t_is=matrix(NaN,N_t,1)
for(i in 0:79){
  aux=as.vector(sst[111:158,130:101,m1.indices[i+1]])
  #print(aux)
  y_t_is[((i)*n+1):((i+1)*n)]=aux
}

#creating matrix data
labels_Ts=matrix(NaN,nrow=N_t,ncol=1)
aux1=0
for (i in 1:Tb) {
  aux1=nts[i]+aux1
  for (j in 1:nts[i]) {
    labels_Ts[aux1-nts[i]+j]=i
  }
}
#labels_Ts
dat=matrix(NaN,nrow=N_t)
dat=cbind(labels_Ts,X_s)
dat=cbind(dat,y_t_is)
dat



#using the data---------------------------------------------------
#we create the seeded intervals
s.inter<-seeded.intervals(Tb,6)
s.inter


#the choice of bandwidth h_bar
H_bar_stimated=Hpi.diag(X_s)
H_bar_stimated
h_bar_aux=H_bar_stimated
h_bar=sqrt(mean(diag(h_bar_aux)))



h=60
tau=0.00145#june
start_time <- Sys.time()
S_1=seedBS(dat, s.inter, h, h_bar, tau, d,n )
end_time <- Sys.time()
end_time - start_time
S_1+1939


##############################
########## July ##############
##############################
m1=7
h=60
tau=0.00144



