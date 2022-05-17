#################
#################
### Functions ###
#################
#################



#Seeded intervals are created
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

#Gaussian Kernel 
K_h<-function(x,h,d){
  res<-(1/(h)^d)*(1/sqrt(2*pi)^d)*exp((-1/2)*(norm(x/h, type = "2"))^2)
  return(res)
}


#dat is the data matrix. sum(n_t) observations, 1+d+n columns, first column= time label, 2:d+1 col= x_t,i vector, d+2:1+d+n col= y_t,i vec
#x is a d-dim vector
p_hat<-function(x, dat,h_bar,d){
  N_t<-nrow(dat)
  aux<-matrix(NaN,nrow=N_t, ncol=1)
  for(i in 1:N_t){
    aux[i]=K_h(x-dat[i,2:(d+1)],h_bar,d)
  }
  aux<-as.vector(aux)
  res<-mean(aux)
  return(res)
}

#Statistic F_t,h
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


#CUSUM statistic
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

