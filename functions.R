f.infmat=function(var.e,DD,N,m,T,omega,gamma)
{
  ### design matrices XX and ZZ and XX and covariancematrix VV
  intercept=rep(1,m)
  time=seq(0,m-1)
  time1=time
  time1[time>T]=T
  time2=rep(0,m)
  time2[time>T]=time[time>T]-T
  
  ZZ=XX=cbind(intercept,time1,time2)
  #print(XX)
  VV=ZZ%*%DD%*%t(ZZ)+var.e*diag(1,m)
  #  print(ZZ)
  
  ### survival function 
  time=seq(0,1,length=m)
  remain = (1 - omega)^time^gamma
  prob.vec <- rep(0, m)
  for(kk in 1:(m-1))
    prob.vec[kk] <- remain[kk] - remain[kk + 1]
  prob.vec[m]=1-sum(prob.vec[1:(m-1)])
  number.vec=prob.vec*N
  
  ### calculate variance under attrition
  infmat <- 0
  for(j in m:2) 
  {
    XX <- XX[1:j,  ]
    ZZ <- ZZ[1:j,  ]
    VV <- VV[1:j,1:j]
    VVinv <- solve(VV)
    infmat <- infmat +  number.vec[j] * t(XX) %*% VVinv %*% XX
  }
  infmat <- infmat +number.vec[1] * outer(XX[1,],XX[1,]) /VV[1,1]
}

