rm(list = ls())

library("matrixcalc")
library("MBESS")

# set.seed(1)
dimFixedEff = 2
dimU1 = 20
dimU2 = 5
dimU3 = 4
dimU4 = 10 

numInd = 200 ## individual in one study (N) 
numSamp = 2000 ## number of repeat (number of studies)  
u1 = matrix( runif(numInd*dimU1,-2,2) ,ncol=dimU1,nrow=numInd) 
for ( i in 1: round(dimU1/5) ) {
	u1[,i] = u1[,dimU1-i] ## cor col. 
}
for ( i in 1: round(numInd/5) ) {
	u1[i,1:round(dimU1/2)] = u1[numInd-i,1:round(dimU1/2)] ## row col. 
}


u2 = matrix( round (runif(numInd*dimU2,0,1)),ncol=dimU2,nrow=numInd) + 1 
# for ( i in 1: round(dimU2/5) ) {
	# u2[,i] = u2[,dimU2-i] ## cor col. 
# }
# for ( i in 1: round(numInd/10) ) {
	# u2[i,1:round(dimU2/2)] = u2[numInd-i,1:round(dimU2/2)] ## row col. 
# }

u3 = NULL 
for ( i in dimU1:(dimU1-1) ){
for ( j in dimU2:(dimU2-1) ){
	u3 = cbind (u3, u1[,i]*u2[,j] )
}
}

u4 = NULL 
for ( i in 1:4 ){
for ( j in (i+1):5 ){
	u4 = cbind (u4, u1[,i]*u1[,j])
}
}

U1 = u1%*%t(u1) # u1 u1^t 
U2 = u2%*%t(u2) 
U3 = u3%*%t(u3)
U4 = u4%*%t(u4)

I = diag(1,numInd)

designX =  matrix( round( runif(numInd*dimFixedEff,0,2) ) ,ncol=dimFixedEff,nrow=numInd) 
designX = cbind ( rep(1,numInd), designX ) 

control_L1penalise = 10 ## high means strong penalise 
control_L1penalise2 = 0 

A = designX ## u1 and u2 append 
dim = numInd

beta1 = rep(0.25,dimU1)
beta2 = rep(0.25,dimU2)
beta3 = rep(1,dimU3)
beta4 = rep(.5,dimU4) 

mu = c(1,.1,.1)
BETA = c(beta1,beta2,beta3,beta4)

sigma1 = 1 # square root version 
sigma2 = 2
sigma3 = 3 # square root version 
sigma4 = 4
sigmaE = 2.5 # square root version

dataY = NULL # generate some data 
for (i in 1:numSamp) {
	b1 = rnorm(dimU1,beta1[1],sd=sigma1)
	b2 = rnorm(dimU2,beta2[1],sd=sigma2) 
	b3 = rnorm(dimU3,beta3[1],sd=sigma3)
	b4 = rnorm(dimU4,beta4[1],sd=sigma4) 
	
	ep = rnorm(numInd,mean=0,sd=sigmaE)
	## make data. 
	y = A%*%mu + u1%*%b1 + u2%*%b2 + u3%*%b3 + u4%*%b4 + ep 
	dataY = cbind(dataY,y)
}

meanY = rowMeans(dataY) ## mean of all vectors 
dataY = dataY - meanY # (y-mu) 
C = outer(dataY[,1],dataY[,1]) ## (y-mu)^t (y-mu) of sample #1
if (numSamp>1) {
	for (i in 2:numSamp) { ## next sample. 
		C = C + outer(dataY[,i],dataY[,i]) ## compute covariance 
	}
}
C = C / numSamp
# if ( is.positive.definite(C) == F ){
	# print("Sample cov is not positive definite.")
# } 

objfunc = function (A,mu,sigma,sigma1,sigma2,sigma3,sigma4) {
	# maximize this original problem 

	mat = (sigma1*U1+sigma2*U2+ sigma3*U3+sigma4*U4 + sigma*I) # covariance divide sig_eps 
	
	a = A%*%mu 
	au = as.vector(a-meanY)
	auR = t(au) %*% solve(mat) %*% au # (a-u)^t Cov^-1 (a-u)
	
	mat2 = solve( mat + diag(.0001,ncol(mat)) ) ## inverse covariance 
	
	ret =  - log(det( mat ) ) - 
			matrix.trace( mat2 %*% C ) - 
			auR -
			control_L1penalise*sum(abs(mu)) 

	return ( ret )
	 
}


findRoots = function(x,tZ){
## x is the eigen 
## tz d^2 + e d - 1
	return ( ( x + sqrt(x^2+4*tZ) ) /2/tZ ) ## first root ? 
}

updateR = function(Gammainv,tY,Y,C,auouter,sigma){
	h = eigen( Gammainv*tY - Y - C/sigma - auouter/sigma  )  
	# svd(X) then X = h$vectors %*% diag(h$values) %*% t(h$vectors)
	# find diagonal entries of M; what if these are not real numbers ?	
	e = sapply(h$values,findRoots,tY)
	rhat = h$vectors %*% diag(e) %*% t(h$vectors)
	rhat = (rhat + t(rhat))/2 ## make symmetric due to rounding errors 
	return (rhat) 
}


updateGamma1 = function(tY,Y,Rinv,gamma2,gamma3,gamma4) {
	x1 = matrix.trace(Y%*%U1) 
	x2 = sum(U1*U1)*tY 
	x3 = tY*sum(U1*(gamma2*U2+gamma3*U3+gamma4*U4+I-Rinv))
	x4 = (x1 - x3) / (x2)
	if (x4<0) {return (0.01)}
	# if (x4>10) {return (10)}
	return( x4 )
}

updateGamma2 = function(tY,Y,Rinv,gamma1,gamma3,gamma4) {
	x1 = matrix.trace(Y%*%U2) 
	x2 = sum(U2*U2)*tY
	x3 = tY*sum(U2*(gamma1*U1+gamma3*U3+gamma4*U4+I-Rinv))
	x4 = (x1 - x3) / (x2)
	if (x4<0) {return (0.01)}
	# if (x4>10) {return (10)}
	return( x4 )
}

updateGamma3 = function(tY,Y,Rinv,gamma1,gamma2,gamma4) {
	x1 = matrix.trace(Y%*%U3) 
	x2 = sum(U3*U3)*tY
	x3 = tY*sum(U3*(gamma1*U1+gamma2*U2+gamma4*U4+I-Rinv))
	x4 = (x1 - x3) / (x2)
	if (x4<0) {return (0.01)}
	# if (x4>10) {return (10)}
	return( x4 )
}

updateGamma4 = function(tY,Y,Rinv,gamma1,gamma2,gamma3) {
	x1 = matrix.trace(Y%*%U4) 
	x2 = sum(U4*U4)*tY
	x3 = tY*sum(U4*(gamma1*U1+gamma2*U2+gamma3*U3+I-Rinv))
	x4 = (x1 - x3) / (x2)
	if (x4<0) {return (0.01)}
	# if (x4>10) {return (10)}
	return( x4 )
}


updateSigmaEp = function( R,C,au,dim ){
	return ( as.numeric ( ( matrix.trace(R%*%C) + t(au)%*%R%*%(au))/dim ) )
}

updateBeta = function( A,R,sigma,dataY ) { ## solve a generalized least square 
	R = R / sigma 
	x1 = solve ( t(A)%*%R%*%A )
	x2 = t(A)%*%R%*%dataY
	return (x1%*%x2)
}

updateY = function(Y,Rinv,gamma1,gamma2,gamma3,gamma4,tY) {
	return ( Y + tY * (Rinv - gamma1*U1-gamma2*U2-gamma3*U3-gamma4*U4-I) )
}

testbound1 = function (){
	x1 = matrix.trace(Y%*%U1)
	x2 = sum(U1*(gamma2_guess*U2+I-Rinv))
	if (x1 < 0 & x2 < 0){
		if (tY > x1/x2) { return (1) } # tY > x1/x2 
	}
	if (x1 > 0 & x2 > 0) {
		if (tY < x1/x2) { return (1) } ## tY < x1/x2 
	}
	if (x1 < 0 & x2 > 0) {
		return (0) 
	}
	if (x1 > 0 & x2 < 0) {
		return (1) # any tY is okay 
	}
	return (0)
}

testbound2 = function (){
	x1 = matrix.trace(Y%*%U2)
	x2 = sum(U2*(gamma1_guess*U1+I-Rinv))
	if (x1 < 0 & x2 < 0){
		if (tY > x1/x2) { return (1) } # tY > x1/x2 
	}
	if (x1 > 0 & x2 > 0) {
		if (tY < x1/x2) { return (1) } ## tY < x1/x2 
	}
	if (x1 < 0 & x2 > 0) {
		return (0) 
	}
	if (x1 > 0 & x2 < 0) {
		return (1) # any tY is okay 
	}
	return (0)
}

##!! add in, gamma sum |x| . this is L1 norm on the effects (betas)

##!! do a FISTA line search 
##!! solve ||Ax - b||^2_2 + gamma ||x||_1 

func_g = function(A,R,x,b) {
	x1 = A%*%x-b
	return ( t(x1) %*% R %*% x1 )
}

gradient_g = function(A,R,x,b) {
	x1 = t(A)%*%R%*%A %*% x
	x2 = t(A) %*% R %*% b 
	return ( x1-x2 )  
}

prox_th = function (xk,tk) {
## prox of L1 norm. 
	tk = tk*control_L1penalise
	xnew = xk 
	gt_t = which ( xk >= tk ) ## where x >= t
	xnew[gt_t] = xk[gt_t] - tk 
	lt_t = which ( xk <= -tk )
	xnew[lt_t] = xk[lt_t] + tk 
	between = which ( xk <= tk && xk >= -tk )
	xnew[between] = 0 
	return (xnew)
}

func_GtX = function(A,R,x,b,t){
## do the line search on "t", set t=1 backtrack until cond#3 is met
	u = x - t* as.vector(gradient_g(A,R,x,b)) ## input the prox 
	return ( ( x - prox_th( u ,t ) ) /t )
}

cond3 = function( x,tk,A,R,b ) {
	left = func_g ( A,R, x - tk*func_GtX(A,R,x,b,tk), b )
	right = func_g ( A,R,x,b ) - tk*t(gradient_g(A,R,x,b))%*%func_GtX(A,R,x,b,tk) + tk/2*norm(func_GtX(A,R,x,b,tk),"2")^2
	return ( left <= right )
}

prox_gradient = function (A,R,b,guessX,tk,numberIter) {
	xk = guessX
	for (iter in 1:numberIter) {
		while (cond3(xk,tk,A,R,b) == FALSE) {
			tk = tk/2
			if (tk < 10^-4){
				return (xk) 
			}
		}
	## !! can modify to have fista here. 
		xk = prox_th ( xk - tk*gradient_g(A,R,xk,b), tk )
	}
	return (xk)
}

## !! ## !! ## !! ## !! ## !!

obval = NULL 
dual_obval = NULL
numberIter = 100
control_residual = 10 # control residual 
tau_incr = tau_decr = 2
movement = betaval = NULL 
err_rg = err_rr = bound = tYv = NULL 

## !! optim ...  
countit = 1 
for ( random_init in 1:1 ) {

	gamma1_guess = runif(1,1,2) ## guess the squared valued 
	gamma2_guess = runif(1,1,2)
	gamma3_guess = runif(1,1,2) ## guess the squared valued 
	gamma4_guess = runif(1,1,2)
	sigma = runif(1,1,2) 
	R = I #
	Y = I ## init values
	mu_guess = runif(dimFixedEff+1,-1,1) 
	beta_guess = runif(dimU1+dimU2+dimU3+dimU4,-1,1) 
	
	movement = rbind (movement, 
		c(gamma1_guess*sigma, gamma2_guess*sigma,
		gamma3_guess*sigma, gamma4_guess*sigma,
		sigma,mu_guess) ) # bind starting points 
	betaval = rbind(betaval,beta_guess)
	tY = .5
	counter = 0 
	step = .01

	initob = objfunc(cbind(A,u1,u2,u3,u4),c(mu_guess,beta_guess),
	sigma,gamma1_guess*sigma,gamma2_guess*sigma,gamma3_guess*sigma,gamma4_guess*sigma )
		
	obval = c(obval,initob) 
	obvali = NULL ; obvali = c(obvali,initob)

	for ( iter in 1:numberIter ) {
		
		countit = countit + 1 
		a = A%*%mu_guess + cbind(u1,u2,u3,u4)%*%beta_guess
		au = as.vector(a-meanY)
		auouter = outer(au,au)  

		# update R
		Gammainv = solve(
			gamma1_guess*U1+gamma2_guess*U2
			+ gamma3_guess*U3+gamma4_guess*U4 +
			+ I + diag(.0001,ncol(R)) )
		R0 = R 
		R = updateR(Gammainv,tY,Y,C,auouter,sigma) # real parts ? 
		# auR = au %*% R %*% au / sigma /numInd
		# if (auR > 1.25 | auR < .75 ) { 
			# R = R / as.numeric(auR)
		# } 
		
		Rinv = solve ( R + diag(.0001,ncol(R)) )
		
		# update gammas 
		gamma1_0 = gamma1_guess
		gamma1_guess = updateGamma1(tY,Y,Rinv,gamma2_guess,gamma3_guess,gamma4_guess)
		
		gamma2_0 = gamma2_guess
		gamma2_guess = updateGamma2(tY,Y,Rinv,gamma1_guess,gamma3_guess,gamma4_guess)
		
		gamma3_0 = gamma3_guess
		gamma3_guess = updateGamma3(tY,Y,Rinv,gamma1_guess,gamma2_guess,gamma4_guess)
		
		gamma4_0 = gamma4_guess
		gamma4_guess = updateGamma4(tY,Y,Rinv,gamma1_guess,gamma2_guess,gamma3_guess)
		
		
		sigma0 = sigma 
		sigma = updateSigmaEp ( R,C,au,dim )
		sig1 = sigma # unchanged 
		
		# update beta 
		tempY = meanY-cbind(u1,u2,u3,u4)%*%beta_guess
		if (control_L1penalise != 0){
			mu_guess = prox_gradient(A,R,tempY,mu_guess,1,80) # penalise 
		} else {
			mu_guess = updateBeta( A,R,sig1,tempY ) 
		}
		
		tempY = meanY-A%*%mu_guess
		tempMat = cbind(u1,u2,u3,u4)
		# if (control_L1penalise != 0){
			# beta_guess = prox_gradient(tempMat,R,tempY,beta_guess,1,80) # penalise 
		# } else {
			beta_guess = updateBeta( tempMat,R,sig1,tempY ) 
		# }
		
		
		res_rr = norm(c(
			gamma1_guess-gamma1_0,
			gamma2_guess-gamma2_0,
			sigma-sigma0,
			gamma3_guess - gamma3_0, 
			gamma4_guess - gamma4_0  ),"2")/sqrt(5) 
		err_rr = c(err_rr,res_rr)

		estCov = gamma1_guess*U1+gamma2_guess*U2
			+ gamma3_guess*U3+gamma4_guess*U4 +
			+ I
		res_rg = norm( Rinv - estCov,"F" )/ numInd ## primal residual 
		err_rg = c(err_rg,res_rg)

		if (iter > 2 ) {
		## keep penalty term constant until end. 
			if ( res_rg >= control_residual * res_rr ) {
				tY = tau_incr * tY
			} 
			if ( res_rr >= control_residual * res_rg ) {
				tY = tY/ tau_decr
			}
			# if ( 0 %in% c(testbound1(),testbound2()) ) {
				# tY = 1 # reset 
			# } 
		}
		
		Y0 = Y 
		if (iter == 1) { 
			Y = updateY(Y0,Rinv,gamma1_guess,gamma2_guess,gamma3_guess,gamma4_guess,10^-6 ) # very small step 
		}
		else if (iter > 2 & iter < 25) { 
			step = step / 2 
			Y = updateY(Y0,Rinv,gamma1_guess,gamma2_guess,gamma3_guess,gamma4_guess,step )
		} else {
			Y = updateY(Y0,Rinv,gamma1_guess,gamma2_guess,gamma3_guess,gamma4_guess,10^-5 )
		}
		
		var1 = gamma1_guess*sigma 
		var2 = gamma2_guess*sigma 
		var3 = gamma3_guess*sigma 
		var4 = gamma4_guess*sigma 
		movement = rbind(movement,c(var1,var2,var3,var4,sigma,mu_guess)) # var
		betaval = rbind(betaval,c(beta_guess))
		
		primal_i = objfunc (cbind(A,u1,u2,u3,u4),c(mu_guess,beta_guess),sig1,var1,var2,var3,var4) 
		best_obval = max(obvali) ## best uptil now. in this one path 
		obvali = c( obvali, primal_i ) ## append 
		obval = c( obval, primal_i )
		if ( abs ((primal_i - best_obval)/best_obval ) < 10^-5 ){
			break 
		}
	}
}

objfunc(cbind(A,u1,u2,u3,u4),c(mu,BETA),sigmaE^2,sigma1^2,sigma2^2,sigma3^2,sigma4^2)
max(obval)
movement = as.data.frame(movement)
movement = cbind (movement,obvali)
pos = which( obval == max(obval) ) 

movement = movement[ 1:pos, ]
tail (movement)
sigma1^2;sigma2^2;sigma3^2;sigma4^2;sigmaE^2;mu

library('xtable')
xtable(tail (movement),digits=3)

j = betaval[pos,]
ind = c ( rep (1,20),rep(2,5),rep(3,4),rep(4,10) )
j = data.frame(j,ind)
j = cbind(j,seq(0,.2,length=nrow(j)) )
names(j) = c("j","ind","y")

library("ggplot2")
theme_set(theme_bw(24))
p = ggplot() +
	geom_boxplot(data=j,aes(y=j,x=factor(ind)),alpha=1,size=2) +
	xlab(expression(mu[b] ) ) +
	ylab("")+ scale_y_continuous(limits=c(.1,1),breaks=seq(0,1,.25)) +
	theme(text = element_text(size=22) , 
	axis.text.x = element_text(colour="grey20",size=20,face="plain") , 
	axis.text.y = element_text(colour="grey20",size=20,face="plain") , 
	plot.title = element_text(lineheight=.8, face="plain", size=24), legend.position="none",
	axis.ticks.x=element_blank()
	) +  
	ggtitle( expression ( paste(mu[b] ) ) ) 
	

png(width=600,height=600,file="muB1penal.png")
p
dev.off()


	
