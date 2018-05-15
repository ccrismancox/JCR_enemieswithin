 stepwise.ols <- function(Y, X, eigen, alpha=0.1) {
	if(length(colnames(X))!= ncol(X)){stop("provide colnames for X")}
	if(length(colnames(eigen))!= ncol(eigen)){stop("provide colnames for eigen")}
	if(!(all(X[,1]==1))){
		X <-cbind(1,X)
		colnames(X)[1] <- "(Intercept)"
	}
	n <- nrow(X)
	ols <- function(X, Y){
			b <- solve(t(X) %*% X) %*% t(X) %*% Y
			resid <- Y - X %*% b
			out <- list(B =b, resid=resid)
	}
	Ftest <- function(R, UR){
				SSR.R <- sum(R$resid^2)
				SSR.U <- sum(UR$resid^2)
				q <- length(UR$B) -length(R$B)
				n <- length(R$resid)
				k <- length(UR$B)
				F <- ((SSR.R-SSR.U)/q ) / (SSR.U/(n-k))
				p <- pf(F, df1=q, df2=(n-k), lower.tail=FALSE)
				return(list(F=F, p=p))
			}
			
    current <- ols(X, Y)
	first <- current
	rnames <- row.names(current$B)
	cand <- eigen[, !colnames(eigen) %in% rnames]
	sel <- cand[,cor(cand,Y)==max(cor(cand, Y))]
	Xnew <- cbind(X, sel)
	colnames(Xnew)[ncol(Xnew)]<-colnames(cand)[cor(cand,Y)==max(cor(cand, Y))]
	pr.F <- Ftest(current, ols(Xnew, Y))
	if(pr.F$p < alpha){
		current <- ols(Xnew, Y)
		rnames <- row.names(current$B)
		cand <- eigen[, !colnames(eigen) %in% rnames]
		iter<-0
		repeat{
		for(i in 1:ncol(cand)){
			trip <- colnames(cand)[i]
			Xtest <- cbind(Xnew, cand[,i])
			colnames(Xtest)[ncol(Xtest)]<-colnames(cand)[i]
			pr.F <- Ftest(ols(Xnew,Y), ols(Xtest, Y))
			if(pr.F$p < alpha){
			print(pr.F$p)
				new2 <- ols(Xtest, Y)
				drop<-c()
				for(j in (length(first$B)+1):(length(new2$B)-1)){
					Xtry <- Xtest[, -j]
					pr.F <- Ftest(ols(Xtry, Y), new2)
					if(pr.F$p > alpha){
						drop <- c(drop, j)
					}
				}
				if(length(drop)>0){
					Xnew <- Xtest[, -drop]
				}else{
					Xnew <- Xtest
				}
				cand <- eigen[, !colnames(eigen) %in% colnames(Xnew)]
				break
			}
		}
		iter <- iter+1
		cat(iter, "\n")
		if(trip==colnames(cand)[ncol(cand)]){break}}
		return(Xnew)
	}else{return(X)}
}

stval= function(x){return(runif(ncol(x))*0.1)}
stepwise.glm <- function(Y, X, lik, stval=stval , method="BFGS", eigen, alpha=0.1, iterlim=200) {
	require(maxLik)
	if(length(colnames(X))!= ncol(X)){stop("provide colnames for X")}
	if(!(all(X[,1]==1))){
		X <-cbind(1,X)
		colnames(X)[1] <- "(Intercept)"
	}
	n <- nrow(X)
	
	LRT <- function(R, UR){
				LRr <- logLik(R)
				LRu <- logLik(UR)
				LR <- -2 * LRr + 2*LRu
				df <- length(coef(UR))- length(coef(R))
				p <- pchisq(LR, df= df,lower.tail=FALSE)
				return(list(LR=LR, p=p))
			}
		
    current <- maxLik(lik, X=X, Y=Y, start=stval(X), method=method , iterlim=iterlim)
	first <- current
	rnames <- colnames(X)
	cand <- eigen[, !colnames(eigen) %in% rnames]
	sel <- cand[,cor(cand,Y)==max(cor(cand, Y))]
	Xnew <- cbind(X, sel)
	colnames(Xnew)[ncol(Xnew)]<-colnames(cand)[cor(cand,Y)==max(cor(cand, Y))]
	test <- maxLik(lik, X=Xnew, Y=Y, start=stval(Xnew), method=method, iterlim=iterlim)
	pr.L <- LRT(current, test)
	if(pr.L$p < alpha){
		current <-  test
		rnames <- colnames(Xnew)
		cand <- eigen[, !colnames(eigen) %in% rnames]
		iter<-0
		repeat{
		for(i in 1:ncol(cand)){
			trip <- colnames(cand)[i]
			Xtest <- cbind(Xnew, cand[,i])
			colnames(Xtest)[ncol(Xtest)]<-colnames(cand)[i]
			restricted <- maxLik(lik, X=Xnew, Y=Y, start=stval(Xnew), method=method, iterlim=iterlim)
			unrestricted <- maxLik(lik, X=Xtest, Y=Y, start=stval(Xtest), method=method, iterlim=iterlim)
			pr.L <- LRT(restricted, unrestricted)
			if(pr.L$p < alpha){
			print(pr.L$p)
				drop<-c()
				for(j in (length(first$estimate)+1):(length(unrestricted$estimate)-1)){
					Xtry <- Xtest[, -j]
					pr.L <- LRT(maxLik(lik, X=Xtry, Y=Y, start=stval(Xtry), method=method, iterlim=iterlim), unrestricted)
					if(pr.L$p > alpha){
						drop <- c(drop, j)
					}
				}
				if(length(drop)>0){
					Xnew <- Xtest[, -drop]
				}else{
					Xnew <- Xtest
				}
				cand <- eigen[, !colnames(eigen) %in% colnames(Xnew)]
				break
			}
		}
		iter <- iter+1
		cat(iter, "\n")
		if(trip==colnames(cand)[ncol(cand)]){break}}
		return(Xnew)
	}else{return(X)}
}

	



