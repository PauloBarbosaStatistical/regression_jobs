library(mvtnorm)
par(mfrow=c(1,1))

Sigma=cbind(c(1,0.1),c(0.1,1)) # matriz de variancias e covariancias
mu=c(0,0) # Vetor de medias
dmvnorm(x=c(0,0), mean=mu, sigma=Sigma)

xy<-seq(-4,4,0.1)
mu=c(0,0)
Sigma=diag(2)
Z1<-matrix(NA,length(xy),length(xy))
for(i in 1:length(xy)){
  for(j in 1:length(xy)){
    Z1[i,j]<-dmvnorm(x=c(xy[i],xy[j]), mean=mu, sigma=Sigma)
  }
}

persp(Z1, phi = 30, theta = -45,col=heat.colors(ncol(Z1)*nrow(Z1),.8),
      ticktype = "detailed",xlab="X")
#.....................................................................#
Sigma=cbind(c(1,0.5),c(0.5,1)) # matriz de variancias e covariancias
mu=c(0,0) # Vetor de m?dias
dmvnorm(x=c(0,0), mean=mu, sigma=Sigma)

xy<-seq(-4,4,0.1)
mu=c(0,0)
Sigma=diag(2)
Z1<-matrix(NA,length(xy),length(xy))
for(i in 1:length(xy)){
  for(j in 1:length(xy)){
    Z1[i,j]<-dmvnorm(x=c(xy[i],xy[j]), mean=mu, sigma=Sigma)
  }
}

persp(Z1, phi = 30, theta = -45,col=heat.colors(ncol(Z1)*nrow(Z1),.8),
      ticktype = "detailed",xlab="X")

#.....................................................................#
Sigma=cbind(c(1,0.9),c(0.9,1)) # matriz de variancias e covariancias
mu=c(0,0) # Vetor de medias
dmvnorm(x=c(0,0), mean=mu, sigma=Sigma)

xy<-seq(-4,4,0.1)
mu=c(0,0)
Sigma=diag(2)
Z1<-matrix(NA,length(xy),length(xy))
for(i in 1:length(xy)){
  for(j in 1:length(xy)){
    Z1[i,j]<-dmvnorm(x=c(xy[i],xy[j]), mean=mu, sigma=Sigma)
  }
}

persp(Z1, phi = 30, theta = -45,col=heat.colors(ncol(Z1)*nrow(Z1),.8),
      ticktype = "detailed",xlab="X")
#......................................................................#

