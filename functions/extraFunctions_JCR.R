ClusterNegbin <- function(mod, cl,data, vcov){     
     b <- coef(mod)
     t <- mod$theta
     scoret <-  function(th, b, X, y){
          mu <- exp(X%*%b)
          dt <- (digamma(th + y) - digamma(th) + log(th) + 1 - 
                      log(th + mu) - (y + th)/(mu + th))
          return(dt)
     }
     scoreb <- function(th, b, X, y){
          return(y*X + (X*as.numeric((-y-th)*exp(X%*%b)))/as.numeric(t+exp(X%*%b)))
     }
     
     
     data <- data[row.names(mod$x),]
     u <- cbind(scoreb(t, b, mod$x, mod$y), scoret(t, b, mod$x, mod$y))
     k <- ncol(u)
     m <- dim(table(data[,cl]))
     clust <- matrix(0L, nrow=m, ncol=k)
     fc <- factor(data[,cl])
     for(i in 1:k){
          clust[,i] <- tapply(u[,i], fc, sum)
     }
     
     
     
     
     clustcov <- vcov %*% ((m/(m-1)) * t(clust)%*% clust)[1:length(b), 1:length(b)] %*% vcov
     ans = list(vcov = clustcov,
                se = sqrt(diag(clustcov)))
     return(ans)
}



ClusterLogit <- function(mod, cl,data){
     score <- function(X,Y, b){
          Xb <- b %*% t(X) ##write it this way to help make grad easier
          grad <- as.numeric((Y - plogis(Xb)))*X ##Grad function
          return(grad)
     }
     data <- data[row.names(mod$x),]
     u <- score(mod$x, mod$y, mod$coef)
     k <- ncol(u)
     m <- dim(table(data[,cl]))
     clust <- matrix(nrow=m, ncol=k)
     fc <- factor(data[,cl])
     for(i in 1:k){
          clust[,i] <- tapply(u[,i], fc, sum)
     }
     clustcov <- vcov(mod) %*% ((m/(m-1)) * t(clust)%*% clust) %*% vcov(mod)
     ans = list(vcov = clustcov,
                se = sqrt(diag(clustcov)))
     return(ans)
} 







av.data<-function(x, x1=NULL, length.out=25){ ##New function to set columns in
     ##a matrix or data.frame at their mean (for continuous variables)
     ##or median (dummies)
     ##ARGUMENTS:
     ##x=Data
     ##x1= variable that should vary (if any)
     ##length.out=nrow of output
     
     
     if(class(x)=="data.frame"){
          x<-as.matrix(x) 
     }##convert data.frames into matrices (I just like them better)
     
     newdata<-x ##want dim(newdata)==dim(x)
     for(i in 1:ncol(x)){
          if(all(x[,i] %in% c(0,1))){ ##find dummies
               newdata[,i]<-median(x[,i]) ##set to median
          }else{ ##not a dummy?
               newdata[,i]<-mean(x[,i]) ##Set to mean
          }     
     }
     
     newdata<-newdata[1:length.out,] ##Trim to desired length
     
     if(!is.null(x1)){ ##if user inputs an x1, then
          X1<-seq(min(x[,x1]),  ##Create a sequence from min
                  max(x[,x1]),  ##to max
                  length.out=length.out) ##with desired length
          newdata[,x1]<-X1 ##And put it in the right column
     }
     
     
     
     return(newdata)
}



NBgrad <- function(b, X, Y){
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
