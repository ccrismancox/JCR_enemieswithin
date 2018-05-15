
NBgrad <- function(b, X, Y){ ##DO NOT USE
     if(!all(X[,1]==1)){
          X<-cbind(1, X)
     }
     
     eta <- b[length(b)]
     alpha<-exp(eta)
     a2 <- exp(2*eta)
     beta<-b[-length(b)]
     XB <- as.numeric(X %*% beta)
     e.XB<-exp(XB)
     e.2XB<-exp(2*XB)
     eXBA <- exp(XB+eta)
     
     dB <- X* (Y-e.XB)/(eXBA+1)
     dA <-( exp(-eta)/(eXBA+1)) * (Y * alpha + 
                                        eXBA * log(eXBA+1) - 
                                        eXBA * digamma(Y + 1/alpha) +
                                        eXBA * digamma(1/alpha) -
                                        eXBA+
                                        log(eXBA + 1)-
                                        digamma(Y + 1/alpha) + 
                                        digamma(1/alpha))
     return(cbind(dB,dA))
     
}


NBreg<-function(b, X, Y, single=FALSE){
     if(!all(X[,1]==1)){
          X<-cbind(1, X)
     }
     
     alpha<-exp(b[length(b)])
     beta<-b[-length(b)]
     e.XB<-exp(X %*% beta)
     LL<-Y*log((alpha*e.XB)/(1+alpha*e.XB))-(1/alpha)*log(1+alpha*e.XB) + 
          lgamma(Y + 1/alpha)-lgamma(Y+1)-lgamma(1/alpha)
     if(single){
          return(sum(LL))
     }else{
          return(LL)
     }
}



FEnbreg <- function(b, X, Y, FE){
     if(!all(X[,1]==1)){X <- cbind(1, X)}
     llik <- matrix(ncol=1, nrow=length(unique(FE)))
     l <- 1
     XB <- X %*% b 
     for(i in unique(FE)){
          lam <- exp(XB[FE==i])
          llik[l] <- lgamma(sum(lam))  + 
               lgamma(sum(Y[FE==i])+1)- 
               # sum(lgamma(Y[FE==i]+1))- 
               lgamma(sum(Y[FE==i])+ sum(lam))+
               sum(lgamma(lam+Y[FE==i]) -
                        lgamma(lam)- lgamma(Y[FE==i]+1))
          l <- l+1
     }
     return(llik)
}

