
###############################################################
# Extract change-points
post_cp_extract <- function(result){
  if(is.null(result)){
    return(list(cp_result90=NULL, cp_result95=NULL, cp_result99=NULL))
  }
  tmp <- matrix(unlist(result), ncol=6, byrow=T)
  cp_result90 <- tmp[,4]
  cp_result95 <- tmp[,5]
  cp_result95[is.na(tmp[,2])] <- NA
  cp_result95 <- cp_result95[!is.na(cp_result95)]
  cp_result99 <- tmp[,6]
  cp_result99[is.na(tmp[,3])] <- NA
  cp_result99 <- cp_result99[!is.na(cp_result99)]
  if(length(cp_result90)==0){cp_result90 <- NULL}
  if(length(cp_result95)==0){cp_result95 <- NULL}
  if(length(cp_result99)==0){cp_result99 <- NULL}
  return(list(cp_result90=cp_result90, cp_result95=cp_result95, cp_result99=cp_result99))
}
vari <- function(Y){
  te <- dim(Y)[1]
  p <- dim(Y)[2]
  Yd <- t(t(Y)-apply(Y,2,mean))
  eb <- matrix(Yd[1:(te-1),],te-1,p)
  ef <- matrix(Yd[2:te,],te-1,p)
  ae <- apply(eb*ef,2,sum)/apply(eb^2,2,sum)
  ee <- ef-(matrix(rep(1,te-1),te-1,1)%*%ae)*eb
  se <- apply(ee^2,2,mean)
  temp <- (1-ae)^2
  ad <- sum((se/temp)^2)
  aese <- ae*se
  den <- ((1-ae)^3*(1+ae))*length(t(t(aese))[1,])
  if(ad!=0 & prod(den)!=0){ 
    a1 <- 4*sum((aese/den)^2)/ad;
    bw <- ceiling(1.1447*((a1*te)^(1/3)))
  }
  if(ad==0 | prod(den)==0) bw <- te-2
  if(bw>(te-2)) bw <- te-2
  if(bw>0) {kern <- 1-((1:bw)/bw)} else kern <- 0
  Yb <- t(Yd)
  lam <- matrix(0,p,p)
  for(h in 1:length(kern)){
    lam <- lam+{matrix(Yb[,1:(te-h)],p,te-h)%*%t(matrix(Yb[,(1+h):te],p,te-h))}/te*kern[h]
  }
  lonrun <- lam+t(lam)+Yb%*%t(Yb)/te
  return(lonrun)
}


BGHK_BS_functional <- function(data, d=2, num_basis=20, rangeval=c(0,1), no_grids=1000, start, end, critical_values,
                               previous_est_cps=c(0,0,0)){
  N <- end-start+1
  if(N<2){
    return(c(NULL, NULL, NULL))
  }
  current_data <- data[,start:end]
  basis <- create.bspline.basis(rangeval=c(rangeval[1], rangeval[2]), nbasis=num_basis, norder=4)
  data.fd <- smooth.basis(0:no_grids/no_grids, current_data, basis)$fd
  data.pca <- pca.fd(data.fd, nharm=d, fdPar(data.fd), centerfns=TRUE)
  scores <- data.pca$scores
  values <- data.pca$values[1:d]
  
  if(d>1){
    forward <- apply(scores,2,cumsum)
    backward <- apply(scores[N:1,],2,cumsum)
  }
  if(d==1){
    forward <- matrix(cumsum(scores),N,1)
    backward <- matrix(cumsum(scores[N:1]),N,1)
  }
  
  # ZhangShao(2011)EJS
  x <- matrix(rep(1:N,d),N,d)/N
  numer <- t(t(forward)-t(x)*forward[N,])/sqrt(N)
  demean <- numer*sqrt(N)
  total <- t(t(demean^2)*values^{-1})
  T_stat <- apply(total,1,sum)/N
  # plot(T_stat)
  if(mean(T_stat)<critical_values[1]){
    return(c(NULL, NULL, NULL))
  }else{
    est_cp <- which.max(T_stat)+start-1
    est_cps <- c(NA, NA, NA)
    est_cps[mean(T_stat)>critical_values] <- est_cp
    return(c(BGHK_BS_functional(data,d,num_basis,rangeval,no_grids,start,est_cp,critical_values,est_cps),
             list(c(previous_est_cps, est_cps)),
             BGHK_BS_functional(data,d,num_basis,rangeval,no_grids,est_cp+1,end,critical_values,est_cps)))
  }
}

HK_BS_functional <- function(data, d=2, num_basis=20, rangeval=c(0,1), no_grids=1000, start, end, critical_values,
                             previous_est_cps=c(0,0,0)){
  N <- end-start+1
  if(N<2*d){
    return(c(NULL, NULL, NULL))
  }
  current_data <- data[,start:end]
  basis <- create.bspline.basis(rangeval=c(rangeval[1], rangeval[2]), nbasis=num_basis, norder=4)
  data.fd <- smooth.basis(0:no_grids/no_grids, current_data, basis)$fd
  data.pca <- pca.fd(data.fd, nharm=d, fdPar(data.fd), centerfns=TRUE)
  scores <- data.pca$scores
  values <- data.pca$values[1:d]
  Sigma <- vari(scores)
  
  if(d>1){
    forward <- apply(scores,2,cumsum)
    backward <- apply(scores[N:1,],2,cumsum)
  }
  if(d==1){
    forward <- matrix(cumsum(scores),N,1)
    backward <- matrix(cumsum(scores[N:1]),N,1)
  }
  
  # ZhangShao(2011)EJS
  x <- matrix(rep(1:N,d),N,d)/N
  numer <- t(t(forward)-t(x)*forward[N,])/sqrt(N)
  total.k<-numer%*%solve(Sigma)%*%t(numer)
  T_stat <- diag(total.k)
  
  if(mean(T_stat)<critical_values[1]){
    return(c(NULL, NULL, NULL))
  }else{
    est_cp <- which.max(T_stat)+start-1
    est_cps <- c(NA, NA, NA)
    est_cps[mean(T_stat)>critical_values] <- est_cp
    return(c(HK_BS_functional(data,d,num_basis,rangeval,no_grids,start,est_cp,critical_values,est_cps),
             list(c(previous_est_cps, est_cps)),
             HK_BS_functional(data,d,num_basis,rangeval,no_grids,est_cp+1,end,critical_values,est_cps)))
  }
}


###############################################################
# SN CP detection for functional time series
SN_BS_functional <- function(data, d=2, num_basis=20, rangeval=c(0,1), no_grids=1000, start, end, critical_values,
                             previous_est_cps=c(0,0,0)){
  N <- end-start+1
  if(N<2*d){
    return(c(NULL, NULL, NULL))
  }
  current_data <- data[,start:end]
  basis <- create.bspline.basis(rangeval=c(rangeval[1], rangeval[2]), nbasis=num_basis, norder=4)
  data.fd <- smooth.basis(0:no_grids/no_grids, current_data, basis)$fd
  data.pca <- pca.fd(data.fd, nharm=d, fdPar(data.fd), centerfns=TRUE)
  scores <- data.pca$scores
  values <- data.pca$values[1:d]
  
  if(d>1){
    forward <- apply(scores,2,cumsum);
    backward <- apply(scores[N:1,],2,cumsum);
  }
  if(d==1){
    forward <- matrix(cumsum(scores),N,1)
    backward <- matrix(cumsum(scores[N:1]),N,1)
  }
  
  # ZhangShao(2011)EJS
  denomi <- array(0,c(d,d,N))
  x <- matrix(rep(1:N,d),N,d)/N
  numer <- t(t(forward)-t(x)*forward[N,])/sqrt(N)
  for(k in 1:(N-1)){
    x1 <- matrix(rep(1:k,d),k,d)/k
    x2 <- matrix(rep(1:(N-k),d),(N-k),d)/(N-k)
    
    if (k>1){
      denomi1 <- t(t(forward[1:k,])-t(x1)*forward[k,])
    }else{
      denomi1 <- t(t(t(forward[1:k,]))-t(x1)*forward[k,])
    }
    
    if (k<(N-1)){
      denomi2 <- t(t(backward[1:(N-k),])-t(x2)*backward[N-k,])
    }else{
      denomi2 <- t(t(t(backward[1:(N-k),]))-t(x2)*backward[N-k,])
    }
    
    denomi[,,k] <- {t(denomi1)%*%denomi1+t(denomi2)%*%denomi2}/N^2
  }
  # print(N)
  G <- rep(NA,N-1)
  for(k in 1:(N-1)){
    G[k] <- numer[k,]%*%solve(denomi[,,k])%*%t(t(numer[k,]))
  }
  if(max(G)<critical_values[1]){
    return(c(NULL, NULL, NULL))
  }else{
    est_cp <- which.max(G)+start-1
    est_cps <- c(NA, NA, NA)
    est_cps[max(G)>critical_values] <- est_cp
    return(c(SN_BS_functional(data,d,num_basis,rangeval,no_grids,start,est_cp,critical_values,est_cps),
             list(c(previous_est_cps, est_cps)),
             SN_BS_functional(data,d,num_basis,rangeval,no_grids,est_cp+1,end,critical_values,est_cps)))
  }
}

###############################################################
# Gaussian kernel 
psi1 <- function(t,s){
  val <- 1.462652
  C1 <- 1/val/2   # The norm of the kernel is 0.5, this controls the temporal dependence
  C1*exp(t^2/2+s^2/2)
}
# Wiener kernel
psi2 <- function(t,s){
  C2 <- 3/sqrt(6) # The norm of the kernel is 0.5, this controls the temporal dependence
  C2*min(t,s)
}


