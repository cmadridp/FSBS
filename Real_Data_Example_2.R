rm(list = ls())
#library(zoo)
#library(splines)
library(fda)
#library(foreach)
#library(doParallel)
library(ks)
#library(tidyverse)
#library(GGally)
#library(fields)
#library(ggplot2)
#library(dplyr)
#library(tidyr)
#library(faux)
library(mvtnorm)
#library(ggplot2)
#library(ggpubr)
#library(ggthemes)
#library(extrafont)
#library(datasets)
#library(ggpubr)
#library(gridExtra)  


library(ncdf4) # package for netcdf manipulation
library(raster) # package for raster manipulation
#library(rgdal) # package for geospatial analysis
library(ggplot2) # package for plotting
library(lubridate)
library(fields)



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

image.plot(lon[200:360],rev(lat),sst[200:360,length(lat):1,m1.indices[80]],xlab='Longitud',ylab='Latitude',cex.axis=.98)
title("Sea Surface Temperature April 2019")
rect(xleft=251,xright = 300,ybottom=1,ytop=50, density=0, col = "black",lwd=2)
rect(xleft=277.5,xright = 286.5,ybottom=13.5,ytop=18.5, density=0, col = "blue",lwd=2)


#plot only close to honduras
lon_Carbibbean=matrix(NaN,10,1)
aux=-1
for(i in 1:10){
  aux=aux+1
  lon_Carbibbean[i]=277.5+aux
}
lat_Caribbean=matrix(NaN,6,1)
aux=-1
for(i in 1:6){
  aux=aux+1
  lat_Caribbean[i]=18.5-aux
}
image.plot(lon_Carbibbean,rev(lat_Caribbean),sst[278:287,77:72,m1.indices[80]],xlab='Longitude',ylab='Latitude',xlim=c(277.5,286.5))
title("SST Australia April 2019")

#The region of interest is a 10x5 grid, i.e. n=50
lat_grid=seq(13.5,18.5,1)
lon_grid=seq(277.5,286.5,1)
#lat_grid=c(-28.5,-27.5,-26.5,-25.5,-24.5)
#lon_grid=c(153.5,154.5,155.5,156.5,157.5,158.5,159.5,160.5,161.5,162.5)
matrix_lon_lat=matrix=matrix(NaN,60,2)
for(i in 0:6){
  matrix_lon_lat[(i*10+1):((i+1)*10),]=lon_grid
  matrix_lon_lat[(i*10+1):((i+1)*10),2]=lat_grid[i+1]
}

#as.vector(sst[104:113,119:115,m1.indices[50]])
#sst[154:163,119:110,m1.indices[50]]
aux_matrix_sst=matrix(NaN,60,80)
for(i in 1:80){
  for(j in 1:6){
    for(k in 1:10){
      if(is.na(sst[278:287,77:72,m1.indices[i]][k,j])){
        sst[278:287,77:72,m1.indices[i]][k,j]=0
      }
    }
  }
  
}
#
#punto=c(103.5,-27.5)
#lonind=which.min(abs(lon-punto[1]))
#latind=which.min(abs(lat-punto[2]))
#sst[lonind,latind,m1.indices[50]]
#-------Putting the data together as in my format-------#

Tb=80
r=1
C=3
#the responses lies in dimension 3
n=60
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
  aux=as.vector(sst[278:287,77:72,m1.indices[i+1]])
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



seeded.intervals <- function(Tb,C, decay = 2, unique.int = T){
  aux <- log(Tb)
  depth <- C*log(aux, base = decay)
  depth <- ceiling(depth)
  M <- sum(2^(1:depth)-1)
  
  boundary_mtx <- matrix(NA, ncol = 3)
  colnames(boundary_mtx) <- c("layer","st", "end")
  boundary_mtx[1, ] <- c(1,0, Tb) #first interval when k=1 in the definition
  
  if(depth>=2)
  {
    for(i in 2:depth){
      int_length <- Tb * (1/decay)^(i-1) #my lk in def seeded intervals
      n_int <- ceiling(Tb/int_length)*2-1 #my nk bar in def seeded intervals sometimes very slight numerical inaccuracies
      aux=floor(seq(0, Tb-int_length, length.out = (n_int)))
      lab=matrix(i,nrow=length(aux),ncol=1)
      aux=cbind(lab,aux)
      boundary_mtx <- rbind(boundary_mtx,cbind(aux, ceiling(seq(int_length, Tb, length.out = (n_int)))))
    }
  }
  
  if(unique.int){return(unique(boundary_mtx))}
  boundary_mtx
}


######################################################################################################################
#CUSUM statistic aux functions (used in algorithm)

#K<-function(x,d){
# aux=(1/sqrt(2*pi)^d)*exp((-1/2)*(norm(x, type = "2"))^2)
#return(aux)
#}

#K_h<-function(x,h,d){
# res<-(1/(h)^d)*K(x/h,d)
#return(res)
#}

K_h<-function(x,h,d){
  res<-(1/(h)^d)*(1/sqrt(2*pi)^d)*exp((-1/2)*(norm(x/h, type = "2"))^2)
  return(res)
}

#X is the data matrix. sum(n_t) observations, 1+d+n columns, first column= time label, 2:d+1 col= x_t,i vector, d+2:1+d+n col= y_t,i vec
#x is a d-dim function
#acortar tiempo (los for deben ser eliminados)
p_hat<-function(x, dat,h_bar,d){
  N_t<-nrow(dat)
  aux<-matrix(NaN,nrow=N_t, ncol=1)
  for(i in 1:N_t){
    aux[i]=K_h(x-dat[i,2:(d+1)],h_bar,d)
  }
  #aux=K_h(x-dat[,2:(d+1)],h_bar,d)
  aux<-as.vector(aux)
  res<-mean(aux)
  return(res)
}


#acortar tiempo (los for deben ser eliminados)
statistic<-function(x,dat,t,h,h_bar,d,phat){
  indices=which(dat[,1]==t)
  n_t=length(indices)
  X_t=dat[indices, 2:(d+1)]
  Y_t=dat[indices,(d+2):(1+d+1) ]
  aux<-matrix(NaN,nrow=n_t, ncol=1)
  for (i in 1:n_t) {
    aux[i]=Y_t[i]*K_h(x-X_t[i],h,d)
  }
  res<-mean(aux)/phat
  return(res)
}


#acortar tiempo (los for deben ser eliminados)
CUSUM<-function(s,e,x,t,h,h_bar,dat,d,phat){
  aux_1=0
  aux=s+1
  if(aux<=t){
    for(i in aux:t){
      aux_1=aux_1+statistic(x,dat,i,h,h_bar,d,phat)
    }
  }
  aux_1=aux_1 * sqrt((e-t)/((e-s)*(t-s)))
  aux_2=0
  aux1=t+1
  if(aux1<=e)
  {
    for(i in aux1:e){
      aux_2=aux_2+statistic(x,dat,i,h,h_bar,d,phat)
    }
  }
  aux_2=aux_2 *sqrt((t-s)/((e-s)*(e-t)))
  return(aux_1-aux_2)
}


######################################################################################################################
#Algorithm
#
#necesito ver la manera que me regrese S
##############(s,e,h,h_bar,tau,rho,dat,s.inter,samp,d,CUSUM_mat)
FSBS<-function(s,e,h,h_bar,tau,rho,dat,s.inter,samp,d,CUSUM_mat){
  ###----------original------------------
  #print("Here")
  n_inter=length(s.inter[,1])
  samp_size=length(samp)
  AD=matrix(ncol=2)
  AD[1,1]=-1
  AD[1,2]=0
  for(i in 1:n_inter){
    alpha=s.inter[i,2]
    beta=s.inter[i,3]
    
    
    #for(j in 1:samp_size){
    alp_bet=beta-alpha
    if((s<= alpha) && (beta<=e) && (alp_bet > 2*rho) && (ceiling(alpha+rho) <= floor(beta-rho) )  ){
      #indices_1=which((CUSUM_mat[,1]==alpha))
      #indices_1
      #indices_2=which((CUSUM_mat[,2]==beta))
      #indices_2
      #indices=which(indices_1==indi.ces_2)
      #indices
      indices1=which(CUSUM_mat[,1]==alpha)
      #indices1
      mat_aux=CUSUM_mat[indices1,]
      indices2=which(mat_aux[,2]==beta)
      #indices2
      mat_aux=mat_aux[indices2,]
      #mat_aux
      #indices
      #indices=which((CUSUM_mat[indices,2]==beta))
      #indices
      #beta
      A_aux=max(mat_aux[,5])
      #A_aux
      #A_aux=max(CUSUM_mat[indices_1,5])
      #CUSUM_mat[indices_1,5] #Aqui esta el error!!!!!
      D_aux=mat_aux[which(mat_aux[,5]==A_aux)[1],4]
      #D_aux
      #D_aux=CUSUM_mat[which(CUSUM_mat[indices,5]==A_aux)[1],4]
    }
    else{
      A_aux=-1
      D_aux=0
    }
    AD=rbind(AD,c(A_aux,D_aux))
    #AD
    #}
  }
  # AD
  
  A_star=max(AD[,1])
  ind=which(AD[,1]==A_star)[1] #First element that reaches the maximum
  if (is.na(A_star)) {
    print('Missing')
    return(NULL)
  }
  else{
    if(A_star > tau){
      S=AD[ind,2]
      aux=AD[ind,2]
      S=append(S,FSBS(s,aux,h,h_bar,tau, rho, dat,s.inter,samp,d,CUSUM_mat))
      #print("Here")
      S=append(S,FSBS(aux,e,h,h_bar,tau, rho, dat,s.inter,samp,d,CUSUM_mat))
      
      return(S)
    }
  }
}


#s.inter= seeded intervals. Has three columns. The first one is a label of the layer the interval belongs to
#d= dim of x_t,i's
#n_bar=dim of y_t,i's
#S empty vector to store changepoints

seedBS <- function (dat, s.inter, h, h_bar, tau, d,n ){
  #Obtain params from dat
  S=0
  S=as.vector(S)
  Tb=max(dat[,1]) #number of times
  N_t=length(dat[,1]) #total number of observations
  
  #Obtain params from s.inter
  n_inter=length(s.inter[,1]) #total number of seeded intervals
  
  #Sample log(Tb) points
  samp_size=floor(log(Tb))
  set.seed(1)
  samp=sample(N_t, samp_size, replace=F)
  
  #initialization of interval
  s=0
  e=Tb
  
  #set rho
  rho=log(Tb)/(n*h^{d})
  
  #Creating statistic matrix
  statistic_matrix=matrix(NaN,nrow=Tb,ncol=samp_size)
  for(i in 1:Tb){
    for(j in 1:samp_size)
    {
      x=dat[samp[j],2:(d+1)]
      phat=p_hat(x,dat,h_bar,d)
      statistic_matrix[i,j]=statistic(x,dat,i,h,h_bar,d,phat)
    }
  }
  
  
  
  #Compute CUSUM statistics 
  n_inter=length(s.inter[,1])
  CUSUM_mat=matrix(NaN, ncol=6)
  aux_count=0
  for(i in 1:n_inter){
    alpha=s.inter[i,2]
    beta=s.inter[i,3]
    if((beta-alpha) > (2*rho)){
      # print('here1')
      t.lo=ceiling(alpha+rho)
      t.up=floor(beta-rho)
      if(t.lo <= t.up){
        # print('here2')
        for(j in 1:samp_size)
        {
          
          #x=dat[samp[j],2:(d+1)]
          #phat=p_hat(x,dat,h_bar,d)
          for(t in t.lo:t.up)
          {
            aux_count=1+aux_count
            
            val_from_statis_matrix_1=statistic_matrix[(alpha+1):t,j]
            #val_from_statis_matrix_1
            sum_val_from_statis_matrix_1=sum(val_from_statis_matrix_1)
            sum_val_from_statis_matrix_1=sum_val_from_statis_matrix_1*sqrt((beta-t)/((beta-alpha)*(t-alpha)))
            #sum_val_from_statis_matrix_1
            
            val_from_statis_matrix_2=statistic_matrix[(t+1):beta,j]
            #val_from_statis_matrix_2
            #beta
            sum_val_from_statis_matrix_2=sum(val_from_statis_matrix_2)
            sum_val_from_statis_matrix_2=sum_val_from_statis_matrix_2*sqrt((t-alpha)/((beta-alpha)*(beta-t)))
            
            CUSUM_x_alpha_beta_t=sum_val_from_statis_matrix_1-sum_val_from_statis_matrix_2
            CUSUM_mat=rbind(CUSUM_mat, c( alpha, beta, j, t , abs(CUSUM_x_alpha_beta_t),aux_count ) )
          }
          
        }
        
      }
      
    }
  }
  CUSUM_mat=CUSUM_mat[-1,]
  
  S=FSBS(s,e,h,h_bar,tau,rho,dat,s.inter,samp,d,CUSUM_mat)
  return(S)
}


#using the data---------------------------------------------------
#we create the seeded intervals
s.inter<-seeded.intervals(Tb,6)
s.inter

#s.inter_train<-seeded.intervals(Tb/2,C)
#s.inter_train


#the choice of bandwidth h_bar
H_bar_stimated=Hpi.diag(X_s)
H_bar_stimated
h_bar_aux=H_bar_stimated
h_bar=sqrt(mean(diag(h_bar_aux)))

####June
tau=0.00001
h=1000
start_time <- Sys.time()
S_1=seedBS(dat, s.inter, h, h_bar, tau, d,n )
end_time <- Sys.time()
end_time - start_time
S_1+1939


####July
tau=0.000007
h=1000
start_time <- Sys.time()
S_1=seedBS(dat, s.inter, h, h_bar, tau, d,n )
end_time <- Sys.time()
end_time - start_time
S_1+1939





