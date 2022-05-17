#################
#################
### Functions ###
#################
#################

#|K-\hat{K}| distance
K_distance=function(vec1,vec2){
  if(length(vec1)==0)
  {
    dH=2
  }
  else{
    dH=abs(length(vec2)-length(vec1))
  }
  return(dH)
}


#Hausdorff distance
Hausdorff = function(vec1,vec2){
  if(length(vec1)==0)
  {
    dH=200
  }
  else{
    dist = matrix(0,nrow = length(vec1),ncol = length(vec2))
    for (i in 1:nrow(dist)){
      for (j in 1:ncol(dist)){
        dist[i,j] = abs(vec1[i] - vec2[j])
      }
    }
    dH = max(max(apply(dist, 2, function(x) min(x) )), max(apply(dist, 1, function(x) min(x) )))
  }
  return(dH)
}   